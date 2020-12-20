
       IMPLICIT NONE
       INTEGER NPARTS,M, I, Nsave_traj, Nsave_thermo, SEED
       INTEGER Ntimesteps1, Ntimesteps2
       PARAMETER(NPARTS=27)
       REAL*8 R(1:NPARTS,1:3),V(1:NPARTS,1:3), F(1:NPARTS,1:3)
       REAL*8 L, Lcell, rho, cutoff, dt, ELJ, ULJ, Ekin, Etot, KIN, T
       REAL*8 Tinst, Tinstant
       COMMON/DADES/ cutoff, L, dt

       SEED = 136465
       CALL SRAND(SEED)
c      Parameters:

       T = 1000D0
       dt = 3d-4
       
       Ntimesteps1 = 5000
       Ntimesteps2 = 5000
       Nsave_thermo = 10
       Nsave_traj = 10

   
       OPEN(UNIT=101, FILE='trajectoryF.xyz')
       OPEN(UNIT=102, FILE='thermodynamicsF.dat')
    1  FORMAT(A1,2X,3(F14.8,2X))
    2  FORMAT(I8,2X,4(F14.8,2X))
       
C -----	INITIALIZATION -------------------------------
       
C       R(1,:) = (/0d0, 0d0,0d0/ )
C       R(2,:) = (/1.2d0, 0d0,0d0/)
C       V(1,:) = (/0d0, 0d0,0d0/)
C       V(2,:) = (/0d0, 0d0,0d0/)

       rho = 0.8442d0
       cutoff = 2.d0
       
       L = (dble(NPARTS)/rho)**(1./3.)
       M = INT(NPARTS**(1./3.))
       Lcell = L/float(M)

       write(*,*) L, M, Lcell

       CALL INIT_POS_SC(NPARTS,Lcell, R)
       V = 0d0
       

C ---- SIMULATION LOOP -----------------------

c        DO I=1, Ntimesteps
c          CALL forceLJ(NPARTS,R, F, ULJ)
c          CALL timeEuler(NPARTS,F,R,V)
        CALL forceLJ(NPARTS,R, F, ULJ)
        DO I=1, Ntimesteps1
          print*,i,Ntimesteps1
          CALL timevV(NPARTS,F,ULJ,R,V)
          CALL ANDERSEN(NPARTS, T, V)
         if(mod(I, Nsave_traj)==0) then
	    WRITE(101, '(I4)') NPARTS
	    WRITE(101, *) ''
	    DO M=1,NPARTS
	      WRITE(101, 1) 'A', R(M,:)
	    ENDDO
	    endif
          
	ENDDO
	
       T = 0.728D0
       dt = 1d-3
       
       DO I=1, 1000
        print*,i,1000
          CALL timevV(NPARTS,F,ULJ,R,V)
          CALL ANDERSEN(NPARTS, T, V)
	ENDDO
	




        DO I=1, Ntimesteps2
          print*,i,Ntimesteps2
          CALL timevV(NPARTS,F,ULJ,R,V)
c          CALL ANDERSEN(NPARTS, T, V)
          if(mod(I, Nsave_traj)==0) then
	       WRITE(101, '(I4)') NPARTS
	       WRITE(101, *) ''
	       DO M=1,NPARTS
	        WRITE(101, 1) 'A', R(M,:)
	       ENDDO
	      endif
	    
          if(mod(I, Nsave_thermo)==0) then
	    Ekin = KIN(NPARTS,V)
            Tinst = Tinstant(NPARTS, Ekin)
	    Etot = Ekin + ULJ
        WRITE(102, 2) I, ULJ, Ekin, Etot, Tinst
        WRITE(*, 2) I, ULJ, Ekin, Etot, Tinst
	  endif
	ENDDO
 
       CLOSE(101)

       END

C ---- SUBROUTINES AND FUNCTIONS----------------------       
       SUBROUTINE INIT_POS_SC(NPARTS,Lcell, R)
       IMPLICIT NONE
         INTEGER NPARTS,I, NX, NY, NZ, M
         REAL*8 Lcell, R(1:NPARTS,1:3)
         
         M = INT(NPARTS**(1./3.))
	 I = 0
         DO NX=1,M
	      DO NY=1,M
	       DO NZ=1,M
	         I = I + 1
	         R(I, 1) = NX*Lcell
	         R(I, 2) = NY*Lcell
	         R(I, 3) = NZ*Lcell	         
          ENDDO
	     ENDDO
        ENDDO
       RETURN
       END SUBROUTINE INIT_POS_SC
       
       SUBROUTINE timeEuler(NPARTS,F, R,V)
       IMPLICIT NONE
       REAL*8 R(1:NPARTS,1:3),V(1:NPARTS,1:3), L
       REAL*8 F(1:NPARTS,1:3), ULJ, cutoff, dt
       INTEGER I, J, NPARTS
       COMMON/DADES/ cutoff, L, dt

       DO I=1, NPARTS
        R(I,:) = R(I,:)+V(I,:)*dt + 0.5d0*F(I,:)*dt*dt
c	Application of Periodic boundary conditions to return particle to
c 	simulation box
        call PBC(R(I,:), L) 
        V(I,:) = V(I,:)+F(I,:)*dt
       ENDDO

       RETURN
       END SUBROUTINE timeEuler

       SUBROUTINE timevV(NPARTS,F,ULJ, R,V)
       IMPLICIT NONE
       REAL*8 R(1:NPARTS,1:3),V(1:NPARTS,1:3), L
       REAL*8 F(1:NPARTS,1:3), ULJ, cutoff, dt
       INTEGER I, J, NPARTS
       COMMON/DADES/ cutoff, L, dt

       DO I=1, NPARTS
        R(I,:) = R(I,:)+V(I,:)*dt + 0.5d0*F(I,:)*dt*dt
c    Application of Periodic boundary conditions to return particle to
c     simulation box
        call PBC(R(I,:), L)
        V(I,:) = V(I,:)+F(I,:)*dt*0.5d0
        CALL forceLJ(NPARTS,R, F, ULJ)
        V(I,:) = V(I,:)+F(I,:)*dt*0.5d0
       ENDDO
       
       RETURN
       END SUBROUTINE timevV

       SUBROUTINE forceLJ(NPARTS, R, F, ULJ)
       IMPLICIT NONE
       REAL*8 R(1:NPARTS,1:3),F(1:NPARTS,1:3), L
       INTEGER I,J, NPARTS
       REAL*8 DX, DY, DZ, DR(1:3), CUTOFF,cf6, cf12, ULJ
       REAL*8 D2, D4, D6,D8, D12,D14, MODF, dt
       COMMON/DADES/ cutoff, L, dt

        F = 0D0
        ULJ = 0D0
        DO I=1,NPARTS
         DO J=I+1,NPARTS
          DR(:) = R(I,:)-R(J,:)
          call PBC(DR(:), L)
C          D2 = DX**2 + DY**2 + DZ**2
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
C           F(I,2) = F(I,2) + MODF*DY
C           F(I,3) = F(I,3) + MODF*DZ
           ULJ = ULJ + 4.*(1./D12 - 1./D6) - 4.*(1./cf12 - 1./cf6)
         ENDIF
         ENDDO
        ENDDO
       RETURN
       END SUBROUTINE forceLJ

      SUBROUTINE pbc(R, L)
       IMPLICIT NONE
       REAL*8 R(1:3), L
       INTEGER I
       DO I=1,3
         IF(R(i).GT.L/2D0) THEN
	  R(I) = R(I) - L
	 ELSE IF(R(i).LT.-L/2D0) THEN
	  R(I) = R(I) + L
	 ENDIF
       ENDDO
       RETURN
       END SUBROUTINE PBC
       
       
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

         REAL*8 FUNCTION Tinstant(NPARTS,Ekin)
         IMPLICIT NONE
           INTEGER NPARTS
           REAL*8 Ekin
           Tinstant = 2d0*Ekin/(3D0*dble(NPARTS)-3d0)
         RETURN
         END FUNCTION
         
       REAL*8 FUNCTION ULJ(NPARTS,R, CUTOFF)
       IMPLICIT NONE
         INTEGER M, I,J, NPARTS
         REAL*8 R(1:NPARTS,1:3), D2, D6, D12, CUTOFF  
         ULJ = 0.0
         DO I=1,NPARTS
	      DO J=I+1,NPARTS
	       D2 = (R(I,1)-R(J,1))**2 + (R(I,2)-R(J,2))**2 +
     &           (R(I,3)-R(J,3))**2 
	      IF(D2.LT.CUTOFF*CUTOFF) THEN
	        D6 = D2*D2*D2
	        D12 = D6*D6
            ULJ = ULJ + 4.*(1./D12 - 1./D6)
	      ENDIF
	      ENDDO
         ENDDO         
       RETURN 
       END FUNCTION


       SUBROUTINE BOXMULLER(SIGMA, X1,X2, XOUT1, XOUT2)
       IMPLICIT NONE
       double precision PI, sigma, x1, x2, xout1, xout2
       PI = 4d0*datan(1d0)
     
       XOUT1 = sigma*dsqrt(-2d0*dlog(1d0-x1))*dcos(2d0*PI*x2)
       XOUT2 = sigma*dsqrt(-2d0*dlog(1d0-x1))*dsin(2d0*PI*x2)
     
       RETURN
       END SUBROUTINE BOXMULLER

       SUBROUTINE ANDERSEN(NPARTS, T, V)
       IMPLICIT NONE
       INTEGER NPARTS, I
       REAL*8 T, sigma, NU, V(1:NPARTS, 1:3), v1, v2, v3
       SIGMA = DSQRT(T)
       NU = 0.1D0
       DO I=1, NPARTS
        IF(rand()<NU) THEN
            CALL BOXMULLER(SIGMA, dble(rand()), dble(rand()), v1, v3)
            CALL BOXMULLER(SIGMA, dble(rand()), dble(rand()), v2, v3)
            V(I, :) = (/v1, v2, v3/)
        ENDIF
       ENDDO
       RETURN
       END SUBROUTINE ANDERSEN