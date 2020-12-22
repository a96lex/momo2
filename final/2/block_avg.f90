program estadistica
    implicit none
    double precision :: time,pot,kin,sum,tem
    double precision :: avgpot,avgkin
    double precision :: avgpotint,avgkinint
    double precision :: allpot,allkin
    double precision :: varpot,varkin
    integer :: N,N2,io,i,bLength,max_bLength,Nblocks
    logical :: equilibrium

    open (unit=101, file = "thermodynamics.dat", action="read")
    open(unit=102, file = "thermodynamic_variances.dat", action="write")
    write(102,*) "Block length","Pot","Kin"

    equilibrium=.false.
    N=0
    N2=0
    i=0
    allpot=0d0
    allkin=0d0
    do
        read (101,*, iostat=io) time,pot,kin,sum,tem
        if(tem.ge.70d0) equilibrium=.true.
        if(equilibrium) then 
            N = N + 1 
            allpot=allpot+pot
            allkin=allkin+kin
        endif
        if(.not.equilibrium) N2 = N2 + 1
        if (io.ne.0) exit
        i=i+1
    enddo
    avgpot=allpot/N
    avgkin=allkin/N

    max_bLength=300

    do bLength=2,max_bLength
        Nblocks=int(dble(N)/dble(bLength))
        rewind(101)
        print*,""
        print*,""
        print*,""
        print*,""
        print*,""
        print*,100*dble(bLength)/dble(max_bLength),"%"

        varpot=0
        varkin=0
        avgpotint=0
        avgkinint=0

        do i=1,N2
            read(101,fmt=*,iostat = io) time,pot,kin,sum,tem
        enddo

        i=0
        do
            i=i+1
            read (101,*, iostat=io) time,pot,kin,sum,tem

            avgpotint=avgpotint+pot/dble(bLength)
            avgkinint=avgkinint+kin/dble(bLength)

            if(mod(i,bLength).eq.0) then
                varpot=varpot+(avgpot-avgpotint)**2/dble(Nblocks)
                varkin=varkin+(avgkin-avgkinint)**2/dble(Nblocks)
                avgpotint=0
                avgkinint=0
            endif

            if (io.ne.0) exit
        enddo
        varpot = dsqrt(varpot)/Nblocks
        varkin = dsqrt(varkin)/Nblocks
        write(102,*) bLength,varpot,varkin
    enddo

end program estadistica