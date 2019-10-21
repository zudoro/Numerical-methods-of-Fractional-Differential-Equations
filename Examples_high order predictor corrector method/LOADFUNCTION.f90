MODULE LOADFUNCTION

	USE GLBVAR
	USE GSQUAD
	USE BASIS
	
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

!! NURBS PU W/ FLAT-TOP 1D

TYPE(FUNCTION_1D) FUNCTION EXSOL(PHYX)

	REAL*8, INTENT(IN) :: PHYX
	
	REAL*8 :: T
	
	T = PHYX
	
	IF (PROBLEM==1) THEN ! Example 1. in "High order predictor corrector method for solving DE of fractional order"
		EXSOL%VAL(0) = T**8.0D0 - 3.0D0*T**(4.0D0 + 0.50D0*NU) + (9.0D0/4.0D0)*T**(nu)
		EXSOL%VAL(1) = 8.0D0*T**7.0D0 - 3.0D0*(4.0D0 + 0.50D0*NU)*T**(3.0D0 + 0.50D0*NU) + (9.0D0/4.0D0)*NU*T**(NU - 1.0D0)
	ELSEIF (PROBLEM==2) THEN ! Example 2 in "High order predictor corrector method for solving DE of fractional order"
		EXSOL%VAL(0) = T**(3.0D0 + NU)
		EXSOL%VAL(1) = (3.0D0 + NU)*T**(2.0D0 + NU)
	ENDIF

END FUNCTION EXSOL


REAL*8 FUNCTION LDF1D(PHYX, LC_COEFF, PATCH, OPTION)

	REAL*8, INTENT(IN) :: PHYX, LC_COEFF(JORDER+1)
	INTEGER, INTENT(IN) :: PATCH
	CHARACTER(LEN=3), INTENT(IN) :: OPTION
	
	TYPE(FUNCTION_1D) :: LC_BS
	REAL*8 :: T, AP_SOL
	INTEGER :: I, J, K
	
	T = PHYX
	AP_SOL = 0.0D0
	
	DO J = 0, JORDER
		LC_BS = GET_PHY_JACOBI_POLY(T, J, PATCH)
		AP_SOL = AP_SOL + LC_COEFF(J+1)*LC_BS%VAL(0)
	ENDDO
	
	IF (OPTION=='FUL') THEN 
		IF (PROBLEM==1) THEN ! Example 1. in "High order predictor corrector method for solving DE of fractional order"
			LDF1D = (40320.0D0/DEXP(GAMMLN(9.0D0 - NU)))*T**(8.0D0 - NU) - 3.0D0*(DEXP(GAMMLN(5.0D0 + 0.50D0*NU))/DEXP(GAMMLN(5.0D0 - 0.50D0*NU)))*T**(4.0D0 - 0.50D0*NU) + (9.0D0/4.0D0)*DEXP(GAMMLN(NU + 1.0D0)) + (1.50D0*T**(0.50D0*NU) - T**(4.0D0))**(3.0D0) - AP_SOL**(1.50D0)
		ELSEIF (PROBLEM==2) THEN ! Example 2. "High order predictor corrector method for solving DE of fractional order"
			LDF1D = (DEXP(GAMMLN(4.0D0 + NU))/6.0D0)*T**3 + T**(3.0D0 + NU) - AP_SOL
! 			LDF1D = (DEXP(GAMMLN(4.0D0 + NU))/6.0D0)*T**3
		ENDIF
	ELSEIF (OPTION=='OFF') THEN 
		IF (PROBLEM==1) THEN ! Example 1. in "High order predictor corrector method for solving DE of fractional order"
			LDF1D = (40320.0D0/DEXP(GAMMLN(9.0D0 - NU)))*T**(8.0D0 - NU) - 3.0D0*(DEXP(GAMMLN(5.0D0 + 0.50D0*NU))/DEXP(GAMMLN(5.0D0 - 0.50D0*NU)))*T**(4.0D0 - 0.50D0*NU) + (9.0D0/4.0D0)*DEXP(GAMMLN(NU + 1.0D0)) + (1.50D0*T**(0.50D0*NU) - T**(4.0D0))**(3.0D0)
		ELSEIF (PROBLEM==2) THEN ! Example 2. "High order predictor corrector method for solving DE of fractional order"
			LDF1D = (DEXP(GAMMLN(4.0D0 + NU))/6.0D0)*T**3 + T**(3.0D0 + NU)
		ENDIF
	ENDIF

END FUNCTION LDF1D


TYPE(FUNCTION_1D) FUNCTION APSOL(PHYPT, COEFF_SOL, PATCH)
	
	REAL*8, INTENT(IN) :: PHYPT
	REAL*8, INTENT(IN) :: COEFF_SOL(DOF)
	INTEGER, INTENT(IN) :: PATCH
	
	TYPE(FUNCTION_1D) :: LOCBS, PU
	
	REAL*8 :: REFPT
	INTEGER :: I, II, J, JJ
	
! 	TYPE(POINT2D) :: PHY_PT
! 	TYPE(MATRIX_22) :: JACOB
	
	REFPT = PHYPT - PATCHBDPT(PATCH,1)
	APSOL%VAL(:) = 0.0D0
	
	DO I = 0, JORDER
! 		LOCBS = JACOBI_POLY(REFPT, I)
! 		LOCBS = GET_LOCBASIS(PHYPT, I, PATCH)
		LOCBS = GET_PHY_JACOBI_POLY(PHYPT, I, PATCH)
! 		PU = PHYPU1D(PHYPT, PATCH)
! 		APSOL = APSOL + COEFF_SOL(II)*LOCBS%VAL(0)*PU%VAL(0)
		APSOL%VAL(0) = APSOL%VAL(0) + COEFF_SOL(SUM(LC_NUMBS(1:PATCH-1))+1+I)*LOCBS%VAL(0)
		APSOL%VAL(1) = APSOL%VAL(1) + COEFF_SOL(SUM(LC_NUMBS(1:PATCH-1))+1+I)*LOCBS%VAL(1)
! 		print*, COEFF_SOL(SUM(LC_NUMBS(1:PATCH-1))+1+I), locbs%val(1)
	ENDDO
	
! 	print*, apsol%val(1)
END FUNCTION APSOL

END MODULE LOADFUNCTION
