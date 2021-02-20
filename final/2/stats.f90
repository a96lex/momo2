program estadistica
    implicit none
    double precision :: time,pot,kin,sum,tem,mom
    double precision :: avgpot,avgkin,avgsum,avgtem,avgmom
    double precision :: allpot,allkin,allsum,alltem,allmom
    double precision :: errpot,errkin,errsum,errtem,errmom
    integer :: N,N2,io,i
    logical :: equilibrium

    open (unit=101, file = "thermodynamics.dat", action="read")
    ! count of the amount of data points after equilibrium
    equilibrium=.false.
    N=0
    N2=0
    i=0
    allpot=0d0
    allkin=0d0
    allsum=0d0
    alltem=0d0
    allmom=0d0
    do
        read (101,*, iostat=io) time,pot,kin,sum,tem,mom
        if(tem.ge.60d0) equilibrium=.true.
        if(equilibrium.and.mod(i,50).eq.0) then 
            N = N + 1 
            allpot=allpot+pot
            allkin=allkin+kin
            allsum=allsum+sum
            alltem=alltem+tem
            allmom=allmom+mom
        endif
        if(.not.equilibrium) N2 = N2 + 1
        if (io.ne.0) exit
    enddo
    avgpot=allpot/N
    avgkin=allkin/N
    avgsum=allsum/N
    avgtem=alltem/N
    avgmom=allmom/N
    errpot=0d0
    errkin=0d0
    errsum=0d0
    errtem=0d0
    errmom=0d0

    rewind(101)
    do i=1,N2
        read(101,fmt=*,iostat = io) time,pot,kin,sum,tem,mom
    enddo
    do
        read (101,*, iostat=io) time,pot,kin,sum,tem,mom
        errpot=errpot+(avgpot-pot)**2
        errkin=errkin+(avgkin-kin)**2
        errsum=errsum+(avgsum-sum)**2
        errtem=errtem+(avgtem-tem)**2
        errmom=errmom+(avgmom-mom)**2
        if (io.ne.0) exit
    enddo
    open(unit=102, file = "thermodynamic_stats.dat", action="write")
    errpot = (errpot/(N-1)/N)**0.5d0
    errkin = (errkin/(N-1)/N)**0.5d0
    errsum = (errsum/(N-1)/N)**0.5d0
    errtem = (errtem/(N-1)/N)**0.5d0
    errmom = (errmom/(N-1)/N)**0.5d0
    write(102,*) "Average potential energy:",avgpot,"Error:",errpot
    write(102,*) "Average kinetic energy:",avgkin,"Error:",errkin
    write(102,*) "Average total energy:",avgsum,"Error:",errsum
    write(102,*) "Average temperature:",avgtem,"Error:",errtem
    write(102,*) "Average momentum:",avgmom,"Error:",errmom
end program estadistica