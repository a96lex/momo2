program final
    implicit none
    double precision, dimension(:,:), allocatable :: r,f,v
    double precision :: L, ro,Ulj,dt,T,kin,kinetic
    integer :: N,i,seed,unit,iter,j
    character(len=100)  :: filename   
    seed = 136465
    unit=1
    filename="test.xyz"
    call srand(seed)
    N=6
    ro=0.8
    Ulj=0
    dt=0.0001
    T=1000
    allocate(r(N**3,3),f(N**3,3),v(N**3,3))
    v=0d0
    f=0d0
    call sccLatice(r,N,ro,L)
    call inizalize_velocities(v,N**3,T)
    v=0d0
    call forceLJ(r,L,N,F,Ulj)
    open (unit=unit, file = filename, status = "old", action="write")
    open (unit=2, file = "temp.dat", status = "old", action="write")
    open (unit=3, file = "ener.dat", status = "old", action="write")
    ! write(unit,fmt=*) N**3
    ! write(unit,fmt=*)
    ! do i=1,N**3
    !     write(unit,fmt=*) "C",r(i,1),r(i,2),r(i,3)
    ! enddo    
    iter=500
    do i=1,iter
        print*
        print*
        print*
        print*
        print*,"melting: ",i*100/iter,"%"
        call verletTherm(r,L,N,F,Ulj,v,dt)
        call thermAndersen(v,0.1d0,1000d0,N)
        kin = kinetic(v,N)
        T = kin * 2 / (3d0 * dble(N**3d0))
        !write(2,fmt=*)i*dt,T
        ! write(3,fmt=*)i*dt,Ulj,kin,Ulj+kin
    enddo
    
    iter=2000
    do i=1,iter
        print*
        print*
        print*
        print*
        print*,"other thing: ",i*100/iter,"%"
        call verletTherm(r,L,N,F,Ulj,v,dt)
            do j=1,N**3
                call pbc(r(j,:),  L)
            enddo
        kin = kinetic(v,N)
        T = kin * 2 / (3d0 * dble(N**3d0))
        write(2,fmt=*)i*dt,T
        write(3,fmt=*)i*dt,Ulj,kin,Ulj+kin
        if(mod(i,50).eq.0) call writeXyz(N,r,unit)
    enddo
    close(unit)
    close(2)
    close(3)
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
                i=i+1
                r(i,1)=dble(x-1)*a
                r(i,2)=dble(y-1)*a
                r(i,3)=dble(z-1)*a
            enddo
        enddo
    enddo
end subroutine

subroutine verletTherm(r,L,N,F,Ulj,v,dt)
    implicit none
    ! i/o variables
    double precision :: r(N**3,3),F(N**3,3),L,Ulj,v(N**3,3),dt
    integer :: N
    r = r + v * dt + 0.5d0 * F * dt ** 2.0d0
    v = v + F * 0.5d0 * dt
    call ljpot(r,L,N,F,Ulj)
    v = v + F * 0.5d0 * dt
end subroutine verletTherm

subroutine ljpot(r,L,N,F,Ulj)
    implicit none
    ! i/o variables
    double precision :: r(N**3,3),f(N**3,3),L,Ulj,test
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
            Ulj = Ulj+ 4.0d0 * (1d0 / d ** 12d0 - 1d0 / d ** 6d0)- 4.0d0 * (2d0 / L ** 12d0 - 1d0 / L ** 6d0)
        enddo
    enddo
end subroutine ljpot

SUBROUTINE forceLJ(R,L,N, F, ULJ)
    IMPLICIT NONE
    INTEGER I,J, NPARTS,N
    REAL*8 R(1:N**3,1:3),F(1:N**3,1:3), L
    REAL*8 DX, DY, DZ, DR(1:3), CUTOFF,cf6, cf12, ULJ
    REAL*8 D2, D4, D6,D8, D12,D14, MODF
    cutoff=L/2d0
    NPARTS=N**3
     F = 0D0
     ULJ = 0D0
     DO I=1,NPARTS
      DO J=I+1,NPARTS
       DR(:) = R(I,:)-R(J,:)
       call PBC(DR(:), L)
!          D2 = DX**2 + DY**2 + DZ**2
       D2 = DR(1)**2 + DR(2)**2 + DR(3)**2
       IF(D2.LT.CUTOFF*CUTOFF) THEN
        D4 = D2*D2
        D6 = D4*D2
        D8 = D6*D2
        D12 = D6*D6
        D14 = D12*D2
        cf6 = cutoff**6d0
        cf12 = cf6*cf6
        MODF = (48D0/D14 - 24D0/D8)
        F(I,:) = F(I,:) + MODF*DR(:)
        F(J,:) = F(J,:) - MODF*DR(:)
!           F(I,2) = F(I,2) + MODF*DY
!           F(I,3) = F(I,3) + MODF*DZ
        ULJ = ULJ + 4.*(1./D12 - 1./D6) - 4.*(1./cf12 - 1./cf6)
      ENDIF
      ENDDO
     ENDDO
    RETURN
END SUBROUTINE forceLJ




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

function kinetic(v,N) result(kin)
    implicit none
    integer :: N,i,j
    double precision :: v(N**3,3),kin
    kin=0
    do i=1,N**3
        do j=1,3
            kin=kin+0.5d0*v(i,j)**2.d0
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
        write(unit,fmt=*) "C",r(i,1),r(i,2),r(i,3)
    enddo    
end subroutine writeXyz