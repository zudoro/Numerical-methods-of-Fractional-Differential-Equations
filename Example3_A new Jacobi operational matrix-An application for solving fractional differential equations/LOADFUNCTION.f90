MODULE LOADFUNCTION

	USE GLBVAR
	USE GSQUAD
	
	IMPLICIT NONE
	
CONTAINS

! REAL*8 FUNCTION EX_DISP(PT2D)
! 
! 	TYPE(POINT2D), INTENT(IN) :: PT2D
! 	TYPE(POINT2D) :: PHY_PT
! 	REAL*8 :: R, THETA, X, Y
! 
! 	PHY_PT = GET_PHY_PT(PT2D)
! 	
! 	X = PHY_PT%X; Y = PHY_PT%Y
! 	
! 	CALL GET_RTHETA(R, THETA, PHY_PT)
! 
! 	IF (PROBLEM.EQ.1) THEN
! 		EX_DISP = X*Y
! 	ELSEIF (PROBLEM.EQ.2) THEN
! 		EX_DISP = (X**2-1.0D0)*(Y**2-1.0D0) + DSIN(X)*DCOS(Y)
! 	ELSEIF (PROBLEM.EQ.3) THEN
! 		EX_DISP = DSQRT(R)*(1.D0-R)*(DSIN(0.5D0*THETA) + DSIN(1.5D0*THETA))
! ! 	ELSEIF (PROBLEM==4) THEN
! ! 		EX_DISP = R**(2.D0/3.D0)*DSIN(2.D0*THETA/3.D0)
! 	ENDIF
! 	
! END FUNCTION EX_DISP
! 
! 
! REAL*8 FUNCTION LDFT2D(PT2D)
! 
! 	TYPE(POINT2D), INTENT(IN) :: PT2D
! 	TYPE(POINT2D) :: PHY_PT
! 	REAL*8 :: X, Y, R, THETA
! 
! 	PHY_PT = GET_PHY_PT(PT2D)
! 	
! 	X = PHY_PT%X; Y = PHY_PT%Y
! 	
! 	CALL GET_RTHETA(R, THETA, PHY_PT)
! 	
! 	IF (PROBLEM.EQ.1) THEN
! 		LDFT2D = 0.D0
! 	ELSEIF (PROBLEM.EQ.2) THEN
! 		LDFT2D = -2.0D0*(-2.0D0 + X**2 + Y**2 - DCOS(Y)*DSIN(X))
! 	ELSEIF (PROBLEM.EQ.3) THEN
! 		LDFT2D = 2.D0*(1.D0 + R + 2.D0*DCOS(THETA))*DSIN(0.5D0*THETA)/R**(1.5D0)
! ! 	ELSEIF (PROBLEM==4) THEN
! ! 		LDFT2D = 0.D0
! 	ENDIF
! 
! 	IF (DABS(R)<=EPS) THEN
! 		LDFT2D = 0.0D0
! 	ENDIF
! 	
! END FUNCTION LDFT2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TYPE(FUNCTION_1D) FUNCTION EXSOL(PHYX)

	REAL*8, INTENT(IN) :: PHYX
	INTEGER :: I, J, K, N, R
	
	
	IF (PROBLEM==0) THEN  ! Test problem of IVP
		EXSOL%VAL(0) = 2.0D0*PHYX**3 + 1.0D0
		
	ELSEIF (PROBLEM==1) THEN ! Example 2. in "On shifted Jacobi spectral approximation for solving fractional DE"
		EXSOL%VAL(0) = PHYX**8 - PHYX**7
		EXSOL%VAL(1) = 8.0D0*PHYX**7 - 7.0D0*PHYX**6
		
	ELSEIF (PROBLEM==2) THEN ! Example 1. in Master Thesis "Multistage Shifted Jacobi Spectral Method for solving linear and nonlinear fractional differential equations"
		EXSOL%VAL(0) = PHYX**4 - PHYX**1.50D0
		EXSOL%VAL(1) = 4.0D0*PHYX**3 - 1.50D0*PHYX**0.50D0
		
	ELSEIF (PROBLEM==3) THEN ! Example 3. in "A high order schema for the numerical solution of the fractional ordinary differential equations"
		EXSOL%VAL(0) = DABS(PHYX**2 - 0.50D0**2)
		IF (PHYX>=0.50D0) THEN 
			EXSOL%VAL(1) = 2.0D0*PHYX
		ELSE 
			EXSOL%VAL(1) = -2.0D0*PHYX
		ENDIF
	
	ELSEIF (PROBLEM==4) THEN
		EXSOL%VAL(0) = DSIN(LAMBDA*PHYX)
		EXSOL%VAL(1) = LAMBDA*DCOS(LAMBDA*PHYX)
	
	ELSEIF (PROBLEM==5) THEN ! Fourier Series
		EXSOL%VAL(:) = 0.0D0
		DO N = 1, FOURIER_TERM
			EXSOL%VAL(0) = EXSOL%VAL(0) + DSIN(2.0D0*PI*(2.0D0*N - 1.0D0)*PHYX)/((2.0D0*N - 1.0D0)*PI)
			EXSOL%VAL(1) = EXSOL%VAL(1) + 2.0D0*DCOS((2.0D0*N - 1.0D0)*PHYX)
		ENDDO
		EXSOL%VAL(0) = 4.0D0*EXSOL%VAL(0)
		EXSOL%VAL(1) = 4.0D0*EXSOL%VAL(1)
	
	ELSEIF (PROBLEM==6) THEN ! Example 3. in "A new Jacobi operational matrix-An application for solving fractional differential equations"
		IF (DABS(PHYX)<=EPS) THEN 
			EXSOL%VAL(0) = 1.0D0
			EXSOL%VAL(1) = 0.0D0
		ELSE 
			EXSOL%VAL(:) = 0.0D0
			DO R = 0, NUM_MF
				EXSOL%VAL(0) = EXSOL%VAL(0) + (-PHYX**NU)**R/DEXP(GAMMLN(NU*R + 1.0D0))
				EXSOL%VAL(1) = EXSOL%VAL(1) + 1.0d0*R*(-PHYX**NU)**(1.0D0*R-1.0D0)*(-NU*PHYX**(NU-1.0D0))/DEXP(GAMMLN(NU*R + 1.0D0))
			ENDDO
		ENDIF
	ENDIF
	
END FUNCTION EXSOL


REAL*8 FUNCTION LDF1D(PHYX)

	REAL*8, INTENT(IN) :: PHYX
	
	IF (PROBLEM==0) THEN ! Test problem of IVP
		LDF1D = 12.0D0*PHYX + 2.0D0*PHYX**3 + 1.0D0
		
	ELSEIF (PROBLEM==1) THEN ! Example 2. in "On shifted Jacobi spectral approximation for solving fractional DE"
		LDF1D = PHYX**9 - PHYX**8 + 56.0D0*PHYX**6 - 42.0D0*PHYX**5 + DSIN(PHYX)*(1.0D0/DSQRT(PI))*((32768.0D0/6435.0D0)*PHYX**(7.50D0) - (2048.0D0/429.0D0)*PHYX**(6.50D0))
		
	ELSEIF (PROBLEM==2) THEN ! Example 1. in Master Thesis "Multistage Shifted Jacobi Spectral Method for solving linear and nonlinear fractional differential equations"
		LDF1D = (DEXP(GAMMLN(5.0D0))/DEXP(GAMMLN(5.0D0 - NU)))*PHYX**(4.0D0 - NU) - (DEXP(GAMMLN(2.50D0))/DEXP(GAMMLN(2.50D0 - NU)))*PHYX**(1.50D0 - NU)
		
	ELSEIF (PROBLEM==3) THEN ! Example 3. in "A high order schema for the numerical solution of the fractional ordinary differential equations"
		IF (PHYX<=0.50D0) THEN 
			LDF1D = (-2.0D0/DEXP(GAMMLN(3.0D0 - NU)))*PHYX**(2.0D0 - NU) + DABS(PHYX**2 - 0.50D0**2)
		ELSE 
			LDF1D = (2.0D0/DEXP(GAMMLN(3.0D0 - NU)))*((2.0D0 - NU)*(PHYX - 0.50D0)**(1.0D0 - NU) + 2.0D0*(PHYX - 0.50D0)**(2.0D0 - NU) - PHYX**(2.0D0 - NU))  + DABS(PHYX**2 - 0.50D0**2)
		ENDIF
		
	ELSEIF (PROBLEM==4) THEN
		LDF1D = GET_FD_SINE(PHYX)
		
	ELSEIF (PROBLEM==5) THEN 
		LDF1D = GET_FD_FOURIER_SINE(PHYX, FOURIER_TERM)
	
	ELSEIF (PROBLEM==6) THEN ! Example 3. in "A new Jacobi operational matrix-An application for solving fractional differential equations"
		LDF1D = 0.0D0
	ENDIF

END FUNCTION LDF1D


REAL*8 FUNCTION GET_FD_SINE(X)

  REAL*8, INTENT(IN) :: X
  
  REAL :: S1, RAT
  INTEGER :: I, J, K, N
  
  S1 = 0.0D0
  N = 300
  
  DO K = 1, N
		RAT = 2.0D0*(-1.0D0)**K*(LAMBDA*X)**(2.0D0*K)/DEXP(GAMMLN(2.0D0*K + M - NU + 1.0D0))
		S1 = S1 + RAT
  ENDDO
  S1 = S1 + (2.0D0/DEXP(GAMMLN(M - NU + 1.0D0)))
  
  GET_FD_SINE = 0.50D0*LAMBDA*X**(M-NU)*S1
  
END FUNCTION GET_FD_SINE


REAL*8 FUNCTION GET_FD_FOURIER_SINE(X, MAX_N)

	REAL*8, INTENT(IN) :: X
	INTEGER, INTENT(IN) :: MAX_N
	
	REAL :: S1, RAT
  INTEGER :: I, J, K, N
  
  N = 150
  GET_FD_FOURIER_SINE = 0.0D0
  
  DO I = 1, MAX_N
		S1 = 0.0D0
		DO K = 1, N
			RAT = 2.0D0*(-1.0D0)**K*(2.0D0*PI*(2.0D0*I - 1.0D0)*X)**(2.0D0*K)/DEXP(GAMMLN(2.0D0*K + M - NU + 1.0D0))
			S1 = S1 + RAT
		ENDDO
		S1 = 0.50D0*2.0D0*PI*(2.0D0*I - 1.0D0)*X**(M-NU)*(S1 + (2.0D0/DEXP(GAMMLN(M - NU + 1.0D0))))
		GET_FD_FOURIER_SINE = GET_FD_FOURIER_SINE + (4.0D0/(2.0D0*I - 1.0D0)*PI)*S1
  ENDDO
  
END FUNCTION GET_FD_FOURIER_SINE

END MODULE LOADFUNCTION
