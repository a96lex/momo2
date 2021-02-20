program final3a
    implicit none
    double precision, dimension(:,:), allocatable :: r,f,v,vini
    double precision :: L,Ulj,dt,T,kin,kinetic,cutoff,mom,momentum,pre,Tbany
    double precision :: kintot,pottot,pretot,kbT,eps,sig,m
    double precision, dimension(5) :: ro
    integer :: N,i,seed,iter1,iter2,Nparts,j
    character(len=100)::radial
    common/props/ kbT
    
    seed = 136465
    call srand(seed)

    ! Parametres
    N=5
    ro=(/0.1d0, 0.2d0, 0.4d0, 0.6d0, 0.8d0/)
    Ulj=0d0
    dt=0.0001
    T=1.5d0
    Nparts=N**3
    eps=1.837d0 !kJ/mol
    sig= 4.1d0 !˚A
    m = 131.29d0 !g/mol
    kbT=1.5*eps

    allocate(r(Nparts,3),f(Nparts,3),v(Nparts,3),vini(Nparts,3))
    open (unit=1, file = "traj.xyz", action="write")
    open (unit=2, file = "thermodynamics.dat", action="write")

    ! main loop
    do j=1,5
        write(radial,"(A,F2.1,A)") "g(r)",ro(j),".dat"
        open (unit=3, file = radial, action="write")
        v=0d0
        f=0d0
        Tbany=1000d0
        call sccLatice(r,N,ro(j),L)
        cutoff=L/2d0
        v=0d0

        ! Making sure all particles are inside the simulation box
        do i=1,Nparts
            call pbc(r(i,:),  L)
        enddo

        ! Heat bath
        kintot=0d0
        pottot=0d0
        pretot=0d0
        call ljpot(r,L,Nparts,F,Ulj,cutoff,pre,ro(j))

        iter1=2000
        iter2=10000
        do i=1,iter2+iter1
            print*
            print*
            print*
            print*
            print*,"Densitat nº",j,i*100/(iter1+iter2),"%"
            call verlet(r,L,Nparts,F,Ulj,v,dt,cutoff,pre,ro(j))
            call thermAndersen(v,0.1d0,Tbany,Nparts)
            kin = kinetic(v,Nparts)
            mom = momentum(v,Nparts)
            if(i.gt.iter1) then 
                kintot=kintot+kin
                pottot=pottot+Ulj
                pretot=pretot+pre
            endif
            T = 2d0*kin/(3d0*dble(Nparts)-3d0)
            if(i.eq.iter1) then
                call writeXyz(Nparts,r)
                Tbany=1.5d0
            endif
        enddo
        call d_distribution(Nparts,r,L,3,100)
        write(2,fmt=*)ro(j)*m/sig**3,pottot*eps/iter2,kintot*eps/iter2,(pottot+kintot)*eps/iter2,pretot/iter2
    enddo

    close(1)
    close(2)
end program final3a 

subroutine inizalize_velocities(v,particles,T)
    implicit none
    integer :: particles,p,c
    double precision :: v(particles,3),T,kin
    kin = 0.d0
    do p = 1 , particles
        do c = 1, 3
            v(p,c) = rand() - 0.5d0
            kin = kin + (1.d0/2.d0) * v(p,c) ** 2.d0
        end do
    end do
    do p = 1 , particles
        do c = 1, 3
            v(p,c) = v(p,c) * sqrt(dble(particles)*T/(2.d0*kin)) 
        end do
    end do
end subroutine

subroutine sccLatice(r,N,ro,L)
    implicit none
    integer :: N,x,y,z,i
    double precision :: ro, r(N**3,3),a,L
    a=1.d0/(ro**(1.d0/3.d0))
    L=a*N
    i=0
    do x=1,N
        do y=1,N
            do z=1,N
                r(i,1)=dble(x)*a
                r(i,2)=dble(y)*a
                r(i,3)=dble(z)*a
                i=i+1
            enddo
        enddo
    enddo
end subroutine

subroutine verlet(r,L,Nparts,F,Ulj,v,dt,cutoff,pre,ro)
    implicit none
    ! i/o variables
    double precision :: r(Nparts,3),F(Nparts,3),L,Ulj,v(Nparts,3),dt,cutoff,pre,ro
    integer :: Nparts,i
    r = r + v * dt + 0.5d0 * F * dt ** 2.0d0
    do i=1,Nparts
        call pbc(r(i,:),L)
    enddo
    v = v + F * 0.5d0 * dt
    call ljpot(r,L,Nparts,F,Ulj,cutoff,pre,ro)
    v = v + F * 0.5d0 * dt
end subroutine verlet

subroutine ljpot(r,L,Nparts,F,Ulj,cutoff,pre,ro)
    implicit none
    double precision :: r(Nparts,3),f(Nparts,3),L,Ulj,ro,kbT
    double precision :: dr(3),d,cutoff,pre,preint
    integer :: i,j,Nparts,k,count
    common/props/ kbT
    F=0.d0
    Ulj=0.d0
    count=0
    pre=0d0
    do i=1,Nparts
        do j=i+1,Nparts
            dr=r(i,:)-r(j,:)
            call pbc(dr,L)
            d=(dr(1)**2.d0+dr(2)**2.d0+dr(3)**2.d0)**0.5d0
            if(d**2.lt.cutoff) then
                preint=0
                count=count+1
                do k=1,3
                    preint=preint+((48.d0/(d**14.d0)-24d0/(d**8.d0))*dr(k))**2
                    f(i,k)=f(i,k)+(48.d0/(d**14.d0)-24d0/(d**8.d0))*dr(k)
                    f(j,k)=f(j,k)-(48.d0/(d**14.d0)-24d0/(d**8.d0))*dr(k)
                enddo
                Ulj = Ulj+ 4.0d0 * (1d0 / d ** 12d0 - 1d0 / d ** 6d0) - 4.0d0 * (1d0 / cutoff ** 12d0 - 1d0 / cutoff ** 6d0)
                pre=pre+d*preint**0.5d0
            endif
        enddo
    enddo
    pre=pre/dble(count)+ro*kbT/(3d0*L**3d0)
end subroutine ljpot

subroutine pbc(dr,  L)
    implicit none
    double precision :: dr(3),L
    integer :: i
    do i=1,3
        if(dr(i).gt.L/2.d0) then
            dr(i)=dr(i)-L
        else if(dr(i).lt.-L/2.d0) then
            dr(i)=dr(i)+L
        endif
    enddo
end subroutine pbc

subroutine thermAndersen(v,nu,T,Nparts)
    implicit none
    double precision :: v(Nparts,3),nu,T
    double precision :: x2,x3,pi,sigma
    integer :: Nparts,i,j
    pi=4d0*datan(1d0)
    sigma=T**0.5d0
    do i=1,Nparts
        if (rand().lt.nu) then
            do j=1,3
                x2=rand()
                x3=rand()
                v(i,j)=sigma*dsqrt(-2d0*(dlog(1d0-x2)))*dcos(2d0*pi*x3)
            enddo
        endif
    enddo
end subroutine thermAndersen

function kinetic(v,Nparts) result(kin)
    implicit none
    integer :: Nparts,i,j
    double precision :: v(Nparts,3),kin
    kin=0
    do i=1,Nparts
        do j=1,3
            kin=kin+0.5d0*v(i,j)**2.d0
        enddo
    enddo
end function kinetic

function momentum(v,Nparts) result(mom)
    implicit none
    integer :: Nparts,i,j
    double precision :: v(Nparts,3),mom
    mom=0
    do i=1,Nparts
        do j=1,3
            mom=mom+v(i,j)
        enddo
    enddo
end function momentum

subroutine writeXyz(Nparts,r)
    implicit none
    integer :: Nparts,i
    double precision :: r(Nparts,3)
    write(unit=1,fmt=*) Nparts
    write(unit=1,fmt=*)
    do i=1,Nparts
        write(unit=1,fmt=*) "Xe",r(i,1),r(i,2),r(i,3)
    enddo    
end subroutine writeXyz

subroutine writeVel(Nparts,vini,v)
    implicit none
    integer :: Nparts,i
    double precision :: vini(Nparts,3),v(Nparts,3)
    do i=1,Nparts
        write(unit=3,fmt=*) (vini(i,1)**2d0+vini(i,2)**2d0+vini(i,3)**2d0)**0.5d0,(v(i,1)**2d0+v(i,2)**2d0+v(i,3)**2d0)**0.5d0
    enddo    
end subroutine writeVel

subroutine d_distribution(Nparts,r,L,unit,Nbins)
    implicit none
    integer :: Nparts, i,j,Nbins,bin,unit
    double precision :: r(Nparts,3),g(Nbins),dr(3),L,maxdist,d,pi
    maxdist=0d0
    g=0d0
    pi=4d0*datan(1d0)

    ! Find the highest distance
    do i=1,Nparts
        do j=1,Nparts
            if (i.ne.j) then
                dr=r(i,:)-r(j,:)
                call pbc(dr,L)
                d=(dr(1)**2.d0+dr(2)**2.d0+dr(3)**2.d0)**0.5d0
                if (d.gt.maxdist) maxdist=d
            endif  
        enddo
    enddo

    ! assign each distance to a bin
    do i=1,Nparts
        do j=1,Nparts
            if (i.ne.j) then
                dr=r(i,:)-r(j,:)
                call pbc(dr,L)
                d=(dr(1)**2.d0+dr(2)**2.d0+dr(3)**2.d0)**0.5d0
                if(d.gt.0d0) then
                    bin=int(d*dble(Nbins)/maxdist)
                    if(bin.gt.Nbins) bin=Nbins
                    if(bin.lt.1) bin=1
                    g(bin)=g(bin)+0.75d0/(pi*(d/2)**3.d0)
                endif
            endif  
        enddo
    enddo

    ! save to file
    do i=1,Nbins
        write(unit=unit,fmt=*) dble(i)*maxdist/Nbins, g(i)/Nparts
    enddo
    

end subroutine d_distribution