program final
    implicit none
    double precision, dimension(:,:), allocatable :: r,f,v,vini
    double precision :: L, ro,Ulj,dt,T,kin,kinetic,cutoff
    integer :: N,i,seed,iter,Nparts
    
    seed = 136465
    call srand(seed)

    ! Files
    open (unit=1, file = "traj.xyz", action="write")
    open (unit=2, file = "thermodynamics.dat", action="write")
    open (unit=3, file = "vel.dat", action="write")
    open (unit=4, file = "g(r)_ini.dat", action="write")
    open (unit=5, file = "g(r)_melted.dat", action="write")
    open (unit=7, file = "g(r)_isolated.dat", action="write")

    ! Parameters
    N=6
    ro=0.8442d0
    Ulj=0d0
    dt=0.0001
    T=1000d0
    Nparts=N**3

    ! initialization
    allocate(r(Nparts,3),f(Nparts,3),v(Nparts,3),vini(Nparts,3))
    v=0d0
    f=0d0
    call sccLatice(r,N,ro,L)
    cutoff=L/2d0
    v=0d0
    call d_distribution(Nparts,r,L,4,100)

    ! Making sure all particles are inside the simulation box
    do i=1,Nparts
        call pbc(r(i,:),  L)
    enddo

    ! Melting loop!
    call ljpot(r,L,Nparts,F,Ulj,cutoff)
    call writeXyz(Nparts,r)
    iter=2000
    do i=1,iter
        print*
        print*
        print*
        print*
        print*,"melting: ",i*100/iter,"%"
        call verlet(r,L,Nparts,F,Ulj,v,dt,cutoff)
        call thermAndersen(v,0.1d0,T,Nparts)
    enddo
    call writeXyz(Nparts,r)
    call d_distribution(Nparts,r,L,5,100)
    
    ! cooling it down
    call inizalize_velocities(v,Nparts,0.728d0)
    vini=v

    ! isolated system simulation
    iter=10000
    do i=1,iter
        print*
        print*
        print*
        print*
        print*,"Isolated system: ",i*100/iter,"%"
        call verlet(r,L,Nparts,F,Ulj,v,dt,cutoff)
        kin = kinetic(v,Nparts)
        T = 2d0*kin/(3d0*dble(Nparts)-3d0)
        if(mod(i,1).eq.0) write(2,fmt=*)i*dt,Ulj,kin,Ulj+kin,T
        if(mod(i,50).eq.0) call writeXyz(Nparts,r)
    enddo
    call writeVel(Nparts,vini,v)
    call d_distribution(Nparts,r,L,7,100)
    close(1)
    close(2)
    close(3)
    close(4)
    close(5)
    close(7)
end program final 

subroutine inizalize_velocities(v,particles,T)
    implicit none
    integer :: particles , p , c
    double precision :: v(particles,3) , T , kin
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

subroutine verlet(r,L,Nparts,F,Ulj,v,dt,cutoff)
    implicit none
    ! i/o variables
    double precision :: r(Nparts,3),F(Nparts,3),L,Ulj,v(Nparts,3),dt,cutoff
    integer :: Nparts,i
    r = r + v * dt + 0.5d0 * F * dt ** 2.0d0
    do i=1,Nparts
        call pbc(r(i,:),L)
    enddo
    v = v + F * 0.5d0 * dt
    call ljpot(r,L,Nparts,F,Ulj,cutoff)
    v = v + F * 0.5d0 * dt
end subroutine verlet

subroutine ljpot(r,L,Nparts,F,Ulj,cutoff)
    implicit none
    ! i/o variables
    double precision :: r(Nparts,3),f(Nparts,3),L,Ulj
    ! other variables
    double precision :: dr(3),d,cutoff
    integer :: i,j,Nparts,k
    F=0.d0
    Ulj=0.d0
    do i=1,Nparts
        do j=i+1,Nparts
            dr=r(i,:)-r(j,:)
            call pbc(dr,L)
            d=(dr(1)**2.d0+dr(2)**2.d0+dr(3)**2.d0)**0.5d0
            if(d**2.lt.cutoff) then
                do k=1,3
                    f(i,k)=f(i,k)+(48.d0/(d**14.d0)-24d0/(d**8.d0))*dr(k)
                    f(j,k)=f(j,k)-(48.d0/(d**14.d0)-24d0/(d**8.d0))*dr(k)
                enddo
                Ulj = Ulj+ 4.0d0 * (1d0 / d ** 12d0 - 1d0 / d ** 6d0) - 4.0d0 * (1d0 / cutoff ** 12d0 - 1d0 / cutoff ** 6d0)
            endif
        enddo
    enddo
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


REAL*8 FUNCTION KIN(NPARTS,V)
IMPLICIT NONE
  INTEGER I,NPARTS
  REAL*8 V(1:NPARTS,1:3)  
  KIN = 0.0
  DO I=1,NPARTS
   KIN = KIN+0.5D0*(V(I,1)*V(I,1)+V(I,2)*V(I,2)+V(I,3)*V(I,3)) 
  ENDDO
RETURN
END FUNCTION

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

subroutine writeXyz(Nparts,r)
    implicit none
    integer :: Nparts,i
    double precision :: r(Nparts,3)
    write(unit=1,fmt=*) Nparts
    write(unit=1,fmt=*)
    do i=1,Nparts
        write(unit=1,fmt=*) "C",r(i,1),r(i,2),r(i,3)
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