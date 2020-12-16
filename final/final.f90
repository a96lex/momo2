program final
    implicit none
    double precision, dimension(:,:), allocatable :: r,f
    double precision :: L, ro
    integer :: N,i,seed
    seed = 165432156
    call srand(seed)
    N=4
    ro=0.8
    allocate(r(N**3,3),f(N**3,3))
    call sccLatice(r,N,ro,L)
    open (unit=1, file = "test.xyz", status = "old", action = "write")
    write(1,fmt=*) N**3
    write(1,fmt=*)
    do i=1,N**3
        write(1,fmt=*) "C",r(i,1),r(i,2),r(i,3)
    enddo
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

subroutine verletTherm(r,L,N,F,Ulj,v,dt,T)
    implicit none
    ! i/o variables
    double precision :: r(N**3,3),F(N**3,3),L,Ulj,v(N**3,3),dt,T,kinetic
    integer :: N
    call ljpot(r,L,N,F,Ulj)
    r = r + v * dt + 0.5d0 * F * dt ** 2.0d0
    v = v + F * 0.5d0 * dt
    call ljpot(r,L,N,F,Ulj)
    v = v + F * 0.5d0 * dt
    call thermAndersen(v,0.1d0,T,N)
    T = kinetic(v,N) * 2 / (3d0 * dble(N)**3d0)
end subroutine verletTherm

subroutine ljpot(r,L,N,F,Ulj)
    implicit none
    ! i/o variables
    double precision :: r(N**3,3),f(N**3,3),L,Ulj
    integer :: N
    ! other variables
    double precision :: dr(3),d
    integer :: i,j,npart,k
    npart=N**3
    F=0.d0
    Ulj=0.d0
    do i=1,npart
        do j=i+1,npart
            dr=r(i,:)-r(j,:)
            call pbc(dr,L)
            d=(dr(1)**2.d0+dr(2)**2.d0+dr(3)**2.d0)**0.5d0
            if(d.lt.L/2.d0) then
                do k=1,3
                    f(i,k)=f(i,k)+(48d0/(d**14d0)-24d0/(d**8d0))*dr(k)
                    f(j,k)=f(j,k)-(48d0/(d**14d0)-24d0/(d**8d0))*dr(k)
                enddo
            endif
            Ulj = Ulj+ 4.0d0 * (1d0 / d ** 12d0 - 1d0 / d ** 6d0)- 4.0d0 * (2d0 / L ** 12d0 - 1d0 / L ** 6d0)
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
    sigma=T**2.d0
    do i=1,npart
        if (rand().lt.nu) then
            x2=rand()
            x3=rand()
            do j=1,3
                v(i,j)=sigma*((-2d0*dlog(1d0-x2))*dcos(2d0*pi*x3))**0.5d0
            enddo
        endif
    enddo
end subroutine thermAndersen

function kinetic(v,N) result(kin)
    implicit none
    integer :: N,i,j
    double precision :: v(N**3,3),kin
    kin=0
    do i=1,N**3
        do j=1,3
            kin=kin+v(i,j)*0.5
        enddo
    enddo
end function kinetic