program final
    implicit none
    double precision, dimension(:,:), allocatable :: r,f,v
    double precision :: L, ro,Ulj,dt,T,kin,sigma,epsilon,m,kinetic
    integer :: N,i,seed,unit,iter
    character(len=100)  :: filename   
    seed = 136465
    unit=1
    filename="test.xyz"
    call srand(seed)
    sigma=4.1d0
    epsilon=1.837d0
    m=131.29d0
    N=4
    ro=0.8*epsilon/sigma**3d0
    Ulj=0
    dt=0.00001
    T=1000
    allocate(r(N**3,3),f(N**3,3),v(N**3,3))
    v=0d0
    f=0d0
    call sccLatice(r,N,ro,L)

    open (unit=unit, file = filename, action="write")
    open (unit=2, file = "temp.dat", action="write")
    open (unit=3, file = "ener.dat", action="write")
    write(unit,fmt=*) N**3
    write(unit,fmt=*)
    do i=1,N**3
        write(unit,fmt=*) "Xe",r(i,1),r(i,2),r(i,3)
    enddo    
    iter=1000
    do i=1,iter
        print*
        print*
        print*
        print*
        print*,"Melting: ",dble(i)*100/dble(iter),"%"
        call verlet(r,L,N,F,Ulj,v,dt*10,T,kin,sigma,epsilon)
        call thermAndersen(v,0.1d0,1000d0,N)
        kin = kinetic(v,N,m)
        T = kin * 2 / (3d0 * dble(N**3d0))
    enddo
    call writeXyz(N,r,unit)
    do i=1,10000
        ! print*
        ! print*
        ! print*
        ! print*
        print*,"Verlet: ",dble(i)*100/dble(10000),"%"
        call verlet(r,L,N,F,Ulj,v,dt,sigma,epsilon)
        call thermAndersen(v,0.1d0,1.5d0,N)
        kin = kinetic(v,N,m)
        T = kin * 2 / (3d0 * dble(N**3d0))
        if(T.lt.100) then
            print*,"heheE",T
            write(2,fmt=*)i*dt,T
        endif
        write(3,fmt=*)i*dt,Ulj,kin,Ulj+kin
    enddo
    call writeXyz(N,r,unit)
    close(unit)
    close(2)
    close(3)
end program final 

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
                i=i+1
                r(i,1)=dble(x-1)*a
                r(i,2)=dble(y-1)*a
                r(i,3)=dble(z-1)*a
            enddo
        enddo
    enddo
end subroutine

subroutine verlet(r,L,N,F,Ulj,v,dt,sigma,epsilon)
    implicit none
    ! i/o variables
    double precision :: r(N**3,3),F(N**3,3),L,Ulj,v(N**3,3),dt,sigma,epsilon,m
    integer :: N
    call ljpot(r,L,N,F,Ulj,sigma,epsilon)   
    r = r + v * dt + 0.5d0 * F * dt ** 2.0d0
    v = v + F * 0.5d0 * dt
    call ljpot(r,L,N,F,Ulj,sigma,epsilon)
    v = v + F * 0.5d0 * dt
end subroutine verlet

subroutine ljpot(r,L,N,F,Ulj,sigma,epsilon)
    implicit none
    ! i/o variables
    double precision :: r(N**3,3),f(N**3,3),L,Ulj,test,sigma,epsilon
    integer :: N
    ! other variables
    double precision :: dr(3),d
    integer :: i,j,npart,k
    npart=N**3
    F=0.d0
    Ulj=0.d0
    test=1
    do i=1,npart
        do j=i+1,npart
            dr=r(i,:)-r(j,:)
            call pbc(dr,L)
            d=(dr(1)**2.d0+dr(2)**2.d0+dr(3)**2.d0)**0.5d0
            if(d.lt.L/2.d0) then
                do k=1,3
                    f(i,k)=f(i,k)+(48.d0/(d**14.d0)-24d0/(d**8.d0))*dr(k)
                    f(j,k)=f(j,k)-(48.d0/(d**14.d0)-24d0/(d**8.d0))*dr(k)
                enddo
            endif
            Ulj = Ulj+ 4.0d0 * epsilon * (( d) ** 12d0 - ( d) ** 6d0)*1e-9
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

subroutine thermAndersen(v,nu,T,N)
    implicit none
    ! i/o variables
    double precision :: v(N**3,3),nu,T
    integer :: N
    ! other variables
    double precision :: x2,x3,pi,sigma
    integer :: npart,i,j
    npart=N**3
    pi=4d0*datan(1d0)
    sigma=T**0.5d0
    do i=1,npart
        if (rand().lt.nu) then
            do j=1,3
                x2=rand()
                x3=rand()
                v(i,j)=sigma*dsqrt(-2d0*(dlog(1d0-x2)))*dcos(2d0*pi*x3)
            enddo
        endif
    enddo
end subroutine thermAndersen

function kinetic(v,N,m) result(kin)
    implicit none
    integer :: N,i,j
    double precision :: v(N**3,3),kin,m
    kin=0
    do i=1,N**3
        do j=1,3
            kin=kin+0.5d0*m*v(i,j)**2.d0
        enddo
    enddo
end function kinetic

subroutine writeXyz(N,r,unit)
    implicit none
    integer :: N,unit,i
    double precision :: r(N**3,3)
    write(unit,fmt=*) N**3
    write(unit,fmt=*)
    do i=1,N**3
        write(unit,fmt=*) "Xe",r(i,1),r(i,2),r(i,3)
    enddo    
end subroutine writeXyz