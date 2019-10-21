MODULE ERRORESTIMATE

	USE BASIS
	USE LOADFUNCTION

    IMPLICIT NONE
    
CONTAINS

!----MAX. NORM ESTIMATE----
	SUBROUTINE MAXNORM(COEFF_SOL, MAXERR)	! MAXERR = (/ABS. NORM, REL. NORM(%)/)

	REAL*8, INTENT(IN) :: COEFF_SOL(DOF)
	REAL*8, INTENT(OUT) :: MAXERR(2)
	
	TYPE(FUNCTION_1D) EX_SOL, AP_SOL
	REAL*8 :: PHY_PT, TMP_ERR, ERR, MAXU
	INTEGER :: I, II, J, JJ, PATCH


	TMP_ERR = 0.0D0
	ERR = 0.0D0
	MAXU = 0.0D0
	
	MAXERR(:) = 0.0D0
	
	OPEN(11, FILE = './data/ext')
	OPEN(21, FILE = './data/app')
	OPEN(31, FILE = './data/err')
	
	DO PATCH = 1, NUMPATCH
		DO I = 1, 1001
			PHY_PT = PATCHBDPT(PATCH,1) + 1.0D0*(I-1)*((PATCHBDPT(PATCH,2) - PATCHBDPT(PATCH,1))/1000.0D0)
			AP_SOL = APSOL(PHY_PT, COEFF_SOL, PATCH)
			EX_SOL = EXSOL(PHY_PT)
			IF (DABS(EX_SOL%VAL(0))>MAXU) THEN
				MAXU = DABS(EX_SOL%VAL(0))
			ENDIF
			TMP_ERR = DABS(EX_SOL%VAL(0) - AP_SOL%VAL(0))
			WRITE(11,*) PHY_PT, EX_SOL%VAL(0)
			WRITE(21,*) PHY_PT, AP_SOL%VAL(0)
			WRITE(31,*) PHY_PT, TMP_ERR
			IF (TMP_ERR.GE.ERR) THEN
				ERR = TMP_ERR
			ENDIF
		ENDDO
	ENDDO
	CLOSE(11);CLOSE(21);CLOSE(31)
	
	MAXERR(1) = ERR
	MAXERR(2) = 100.0D0*ERR/MAXU
	
END SUBROUTINE MAXNORM


SUBROUTINE H_NORM_ERROR(COEFF_SOL, H_NORM)
	
	REAL*8, INTENT(OUT) :: H_NORM(2,2)	! H_NORM(1,:) = (/ABS. L2-NORM, REL. L2-NORM(%)/),  H_NORM(2,:) = (/ABS. H1-NORM, REL. H1-NORM(%)/)
	REAL*8, INTENT(IN) :: COEFF_SOL(DOF)

	TYPE(FUNCTION_1D) :: P, EX_SOL, AP_SOL
	REAL*8 :: H0(2), H1(2)
	INTEGER :: I, J, K, II, JJ, KK, PATCH
	
	H0(:) = 0.0D0
	H1(:) = 0.0D0
	
	DO PATCH = 1, NUMPATCH
		DO K = 1, NUMGSPT
			AP_SOL = APSOL(GSPT(1, K) + PATCHBDPT(PATCH,1), COEFF_SOL, PATCH)
			EX_SOL = EXSOL(GSPT(1, K) + PATCHBDPT(PATCH,1))
! 			print*, patch, k, ap_sol%val(0:1), ex_sol%val(0:1)
			H0(1) = H0(1) + (AP_SOL%VAL(0) - EX_SOL%VAL(0))**2*GSW(1, K)
			H0(2) = H0(2) + EX_SOL%VAL(0)**2*GSW(1, K)
			H1(1) = H1(1) + (AP_SOL%VAL(1) - EX_SOL%VAL(1))**2*GSW(1, K)
			H1(2) = H1(2) + EX_SOL%VAL(1)**2*GSW(1, K)
		ENDDO
	ENDDO
	
	H_NORM(1,1) = DSQRT(DABS(H0(1)))
	H_NORM(1,2) = DSQRT(DABS(H0(1)/H0(2)))*100.0D0
	
	H_NORM(2,1) = DSQRT(DABS(H0(1) + H1(1)))
	H_NORM(2,2) = DSQRT(DABS(H0(1) + H1(1))/(H0(2) + H1(2)))*100.0D0

END SUBROUTINE H_NORM_ERROR


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
		LOCBS = JACOBI_POLY(REFPT, I)
! 			LOCBS = GET_LOCBASIS(PHYPT, I, PATCH)
! 			LOCBS = GET_PHY_JACOBI_POLY(PHYPT, I, PATCH)
! 			PU = PHYPU1D(PHYPT, PATCH)
! 			APSOL = APSOL + COEFF_SOL(II)*LOCBS%VAL(0)*PU%VAL(0)
		APSOL%VAL(0) = APSOL%VAL(0) + COEFF_SOL(SUM(LC_NUMBS(1:PATCH-1))+1+I)*LOCBS%VAL(0)
		APSOL%VAL(1) = APSOL%VAL(1) + COEFF_SOL(SUM(LC_NUMBS(1:PATCH-1))+1+I)*LOCBS%VAL(1)
! 		print*, COEFF_SOL(SUM(LC_NUMBS(1:PATCH-1))+1+I), locbs%val(1)
	ENDDO
! 	print*, apsol%val(1)
END FUNCTION APSOL

END MODULE ERRORESTIMATE