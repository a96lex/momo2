program estadistica
    implicit none
    double precision :: time,pot,kin,sum,tem
    double precision :: avgpot,avgkin,avgsum,avgtem
    double precision :: allpot,allkin,allsum,alltem
    double precision :: errpot,errkin,errsum,errtem
    integer :: N,N2,io,i
    logical :: equilibrium

    open (unit=101, file = "thermodynamics.dat", action="read")
    ! count of the amount of data points after equilibrium
    equilibrium=.false.
    N=0
    N2=0
    do
        read (101,*, iostat=io) time,pot,kin,sum,tem
        if(tem.ge.60d0) equilibrium=.true.
        if(equilibrium) then 
            N = N + 1 
            allpot=allpot+pot
            allkin=allkin+kin
            allsum=allsum+sum
            alltem=alltem+tem
        endif
        if(.not.equilibrium) N2 = N2 + 1
        if (io.ne.0) exit
    enddo
    avgpot=allpot/N
    avgkin=allkin/N
    avgsum=allsum/N
    avgtem=alltem/N
    errpot=0
    errkin=0
    errsum=0
    errtem=0

    rewind(101)
    print*,N2
    do i=1,N2
        read(101,fmt=*,iostat = io) time,pot,kin,sum,tem
    enddo
    do
        read (101,*, iostat=io) time,pot,kin,sum,tem
        errpot=errpot+(avgpot-pot)**2
        errkin=errkin+(avgkin-kin)**2
        errsum=errsum+(avgsum-sum)**2
        errtem=errtem+(avgtem-tem)**2
        if (io.ne.0) exit
    enddo
    errpot = dsqrt(errpot/(N*(N-1)))
    errkin = dsqrt(errkin/(N*(N-1)))
    errsum = dsqrt(errsum/(N*(N-1)))
    errtem = dsqrt(errtem/(N*(N-1)))
    print*,avgpot,errpot
    print*,avgkin,errkin
    print*,avgsum,errsum
    print*,avgtem,errtem
end program estadistica