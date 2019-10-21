MODULE GEOMETRY

	USE GSQUAD
	USE KNOT_HANDLING
	USE PATCH_MAPPING

	IMPLICIT NONE

CONTAINS

SUBROUTINE GET_GEO()
	
	INTEGER :: I, J, II, JJ, K, KK, PATCH
	
	! Set degree and knot vector of Bernstein polynomials 
	!--------------------------------------------------------------------------------------!
	BASIS_KVEC%POLY_ORDER = 2
	
	BASIS_KVEC%LENGTH = 5
	
! 	DO I = 1, NUMPATCH
		BASIS_KVEC%KNOTS(0:BASIS_KVEC%LENGTH) = (/0.D0, 0.D0, 0.D0, 1.D0, 1.D0, 1.D0/)
! 	ENDDO
	
! ---------------------- [[[[[ DEGREE ELEVATION OF B-SPLINE ]]]]] ---------------------	
! 	DO I = 1, NUMPATCH
		IF (JORDER<=BASIS_KVEC%POLY_ORDER) THEN
		ELSE
			DO K = 1, JORDER - BASIS_KVEC%POLY_ORDER
				BASIS_KVEC = DEGREE_ELEVATION(BASIS_KVEC, .TRUE.)
			ENDDO
		ENDIF
! 	ENDDO
	!--------------------------------------------------------------------------------------!
	
	
	! Set local number of basis, degree of freedom, and global index with index of local basis
	!--------------------------------------------------------------------------------------!
	DO PATCH = 1, NUMPATCH
		LC_NUMBS(PATCH) = JORDER + 1
	ENDDO
	DOF = SUM(LC_NUMBS(:))
	
! GLOBAL INDEX
	ALLOCATE(NDX(2,DOF))
	
	K = 0
	DO PATCH = 1, NUMPATCH
		DO I = 1, LC_NUMBS(PATCH)
			K = K + 1
			NDX(1,K) = PATCH
			NDX(2,K) = I-1
		ENDDO
	ENDDO
!--------------------------------------------------------------------------------------

	WRITE(*,*)
	WRITE(*,*) '<<< SET GLOBAL INDEX : DONE >>>'
	WRITE(*,*)
	
	!! SET INDEX OF BASIS FUNCTIONS IMPOSED ESSENTIAL BOUNDARY CONDITION
!--------------------------------------------------------------------------------------
! 	K = 0
! 	BD_DOF = 0
! 	DO I = 1, DOF
! ! 		IF (NDX(1,I)==0 .OR. NDX(1,I)==NUMBS(1)-1 .OR. NDX(2,I)==0 .OR. NDX(2,I)==NUMBS(2)-1) THEN
! 		IF (NDX(1,I)==0 .OR. NDX(1,I)==NUMBS(1)-1 .OR. NDX(2,I)==NUMBS(2)-1) THEN
! 			IF (NDX(2,I)==0) THEN
! 			ELSE
! 				K = K + 1
! 				BDNDX(K)%LC_NDX(:) = NDX(:,I)
! 				BDNDX(K)%GL_NDX = I
! 				BD_DOF = BD_DOF + 1
! 			ENDIF
! 		ENDIF
! 	ENDDO
! 	
! 	BDNDX(:)%LC_NUM = K
! 	
! 	ALLOCATE(LST_SOL(K), BD_COL_PT(K))

! 	BD_DOF = 2
! 	BDNDX(1)%LC_NDX(:) = NDX(:, 1)
! 	BDNDX(1)%GL_NDX = 1
! 	BDNDX(2)%LC_NDX(:) = NDX(:, DOF)
! 	BDNDX(2)%GL_NDX = DOF
!--------------------------------------------------------------------------------------

! 	WRITE(*,*)
! 	WRITE(*,*) '<<< SET LOCAL INDEX CORRESPONDING BASIS FUNCTIONS ON BOUNDARY IMPOSED ESSENTIAL BOUNDARY CONDITION : DONE >>>'
! 	WRITE(*,*)

	!! GENERATE GREEVILLE POINTS
!--------------------------------------------------------------------------------------
! 	ALLOCATE(GREV_PTX(NUMBS(1)), GREV_PTY(NUMBS(2)))
! 	
! 	CALL FIND_GREVILLE_PT(GREV_PTX, BASIS_KVEC(1))
! 	CALL FIND_GREVILLE_PT(GREV_PTY, BASIS_KVEC(2))
	
	!! FIND GREVILLE POINTS ON BOUNDARY
! 	K = 0
! 	DO I = 1, NUMBS(1)
! 		DO J = 1, NUMBS(2)
! ! 			IF (DABS(GREV_PTX(I))<=EPS .OR. DABS(GREV_PTX(I)-1.0D0)<=EPS .OR. DABS(GREV_PTY(J))<=EPS .OR. DABS(GREV_PTY(J)-1.0D0)<=EPS) THEN
! 			IF (DABS(GREV_PTX(I))<=EPS .OR. DABS(GREV_PTX(I)-1.0D0)<=EPS .OR. DABS(GREV_PTY(J)-1.0D0)<=EPS) THEN
! 				IF (DABS(GREV_PTY(J))<=EPS) THEN
! 				ELSE
! 					K = K + 1
! 					BD_COL_PT(K) = POINT2D(GREV_PTX(I), GREV_PTY(J))
! 				ENDIF
! 			ENDIF
! 		ENDDO
! 	ENDDO
!--------------------------------------------------------------------------------------


!! Set Collocation Points
!--------------------------------------------------------------------------------------
! 	LOC_NUMCOL = INT(DOF/(UBOUND(JETA, 1) - 1))
	
! 	ALLOCATE(COL_PT(DOF), COL_W(DOF))
! 	
! 	DO I = 1, UBOUND(JETA, 1) - 1
! 		CALL GAULEG(JETA(I), JETA(I+1), COL_PT(1 + (I-1)*LOC_NUMCOL:I*LOC_NUMCOL), COL_W(1 + (I-1)*LOC_NUMCOL:I*LOC_NUMCOL), LOC_NUMCOL)
! 	ENDDO
! 	COL_PT(1) = 0.0D0
! 	COL_PT(DOF) = 1.0D0
!--------------------------------------------------------------------------------------


	!! FIND END POINTS OF EACH PATCHES ON THE PHYSICAL DOMAIN
!--------------------------------------------------------------------------------------
	PATCHBDPT(1,1) = DOMAIN(1)
	PATCHBDPT(NUMPATCH,2) = DOMAIN(2)
	
	DELTAH = (PATCHBDPT(NUMPATCH,2) - PATCHBDPT(1,1))/(1.0D0*NUMPATCH)
	PATCHBDPT(1,2) = PATCHBDPT(1,1) + DELTAH
	
	IF (NUMPATCH==1) THEN 
		GOTO 76
	ENDIF
	
	IF (NUMPATCH>=3) THEN 
		DO I = 2, NUMPATCH - 1
			PATCHBDPT(I,1) = PATCHBDPT(I-1,2)
			PATCHBDPT(I,2) = PATCHBDPT(I,1) + DELTAH
		ENDDO
	ENDIF

	PATCHBDPT(NUMPATCH,1) = PATCHBDPT(NUMPATCH - 1, 2)
	
	76 CONTINUE
!--------------------------------------------------------------------------------------
! 	ALLOCATE(JETA(2*NUMPATCH - 1,2))
! 	
! 	IF (NUMPATCH==1) THEN 
! 		JETA(1,1) = PATCHBDPT(1,1)
! 		JETA(1,2) = PATCHBDPT(1,2)
! 		GOTO 77
! 	ENDIF
! 	
! 	JETA(1,1) = PATCHBDPT(1,1)
! 	JETA(1,2) = PATCHBDPT(1,2) - DELTA
! 	
! 	JETA(UBOUND(JETA,1),1) = PATCHBDPT(NUMPATCH,1) + DELTA
! 	JETA(UBOUND(JETA,1),2) = PATCHBDPT(NUMPATCH,2)
! 	
! 	DO I = 2, UBOUND(JETA,1) - 1
! 		JETA(I,1) = JETA(I-1,2)
! 		IF (MOD(I,2)==0) THEN
! 			JETA(I,2) = JETA(I,1) + 2.0D0*DELTA
! 		ELSE 
! 			JETA(I,2) = PATCHBDPT(INT(0.5*(I+1)),2) - DELTA
! 		ENDIF
! 	ENDDO
	
! 	77 CONTINUE 
	
!--------------------------------------------------------------------------------------
	
	
	!! GENERATE GAUSS POINTS ON THE PHYSICAL PATCHES
!--------------------------------------------------------------------------------------
	!! PU-FEM
	!--------------------------------------------------------------------!
! 	NUMGSPT = INT(0.5*(JORDER + PUORDER)) + 5
! 	ALLOCATE(GSPT(UBOUND(JETA,1),NUMGSPT), GSW(UBOUND(JETA,1),NUMGSPT))
	!--------------------------------------------------------------------!
	!! Discontinuous Galerkin
	!--------------------------------------------------------------------!
	NUMGSPT(1) = 11
	NUMGSPT(2) = 11
	NUMGSPT(3) = 11
	NUMGSPT(4) = 11
	
	INNER_NUMGSPT = 11
	
	ALLOCATE(GSPT(NUMPATCH, NUMGSPT(1)), GSW(NUMPATCH, NUMGSPT(1)))
	ALLOCATE(REFGSPT(NUMGSPT(2)), REFGSW(NUMGSPT(2)))
	ALLOCATE(SPGSPT(NUMPATCH, NUMGSPT(4)), SPGSW(NUMPATCH, NUMGSPT(4)))
	
	ALLOCATE(INNER_GSPT(INNER_NUMGSPT), INNER_GSW(INNER_NUMGSPT))
	
! 	ALLOCATE(GSPT(1,NUMGSPT), GSW(1,NUMGSPT))
	
	ALLOCATE(GJPT(NUMGSPT(3)), GJW(NUMGSPT(3)))
	!--------------------------------------------------------------------!
	
	DO I = 1, NUMPATCH
! 	DO I = 1, 1
		CALL GAULEG(PATCHBDPT(I,1), PATCHBDPT(I,2), GSPT(I,:), GSW(I,:), NUMGSPT(1))
		CALL GAULEG(PATCHBDPT(I,1), PATCHBDPT(I,2), SPGSPT(I,:), SPGSW(I,:), NUMGSPT(4))
! 		CALL GAUJAC(GJPT(I,:), GJW(I,:), NUMGSPT, ALF, BET)
	ENDDO
	
	CALL GAULEG(-1.0D0, 1.0D0, REFGSPT(:), REFGSW(:), NUMGSPT(2))
	CALL GAUJAC(GJPT(:), GJW(:), NUMGSPT(3), ALF, BET)
	
	!! FIND SUB-LIMITS OF INTEGRAL INDEX
!--------------------------------------------------------------------------------------
! 	ALLOCATE(SUBIRNDX(NUMPATCH,NUMPATCH,3))
! 	SUBIRNDX(:,:,:) = 0
! 	
! 	IF (NUMPATCH==1) THEN 
! 		SUBIRNDX(1,1,1) = 1
! 		GOTO 78
! 	ENDIF
! 	
! 	DO I = 1, NUMPATCH
! 		DO J = 1, NUMPATCH
! 			IF (I==J) THEN
! 				IF (I==1) THEN
! 					SUBIRNDX(I,J,1:2) = (/1, 2/)
! ! 					PRINT*, 'SUBIRNDX(1,1,:)', I,J, SUBIRNDX(I,J,:)
! 				ELSEIF (I==NUMPATCH) THEN
! 					SUBIRNDX(I,J,1:2) = (/UBOUND(JETA,1)-1, UBOUND(JETA,1)/)
! 				ELSE
! 					SUBIRNDX(I,J,1:3) = (/2*I-2, 2*I-1, 2*I/)
! 				ENDIF
! 			ELSEIF (I==(J-1)) THEN
! 				SUBIRNDX(I,J,1) = 2*I
! 			ELSEIF (I==(J+1)) THEN
! 				SUBIRNDX(I,J,1) = 2*J
! 			ENDIF
! 		ENDDO
! 	ENDDO
! 	
! 	78 CONTINUE
!--------------------------------------------------------------------------------------


!!  GENERATE BINOMIAL COEFFICIENTS
!--------------------------------------------------------------------------------------
! 	DO I = 1, MAX_DIFF_ORDER
! 		DO J = 0, MAX_DIFF_ORDER
! 			BINOM(I,J) = 1.0D0
! 		ENDDO
! 	ENDDO
! 	
! 	BINOM(1,0) = 1.0D0
! 	BINOM(1,1) = 1.0D0
! 	DO I = 2, MAX_DIFF_ORDER
! 		BINOM(I,0) = BINOM(I-1,0)
! 		DO J = 1, I-1
! 			BINOM(I,J) = BINOM(I-1,J-1) + BINOM(I-1,J)
! 		ENDDO
! 		BINOM(I,I) = BINOM(I-1,I-1)
! 	ENDDO
!--------------------------------------------------------------------------------------

END SUBROUTINE GET_GEO

! SUBROUTINE  FIND_GREVILLE_PT(GREV_PTS, KVEC)
! 	TYPE(KNOT_VECTOR), INTENT(IN) :: KVEC
! 	REAL*8, INTENT(OUT) :: GREV_PTS(KVEC%LENGTH - KVEC%POLY_ORDER)
! 	
! 	REAL*8 :: SUM1
! 	INTEGER :: I, J
! 	
! 	GREV_PTS(:) = 0.0D0
! 				
! 	DO I = 0, UBOUND(GREV_PTS,1) - 1
! 		SUM1 = 0.0D0
! 		DO J = 1, KVEC%POLY_ORDER
! 			SUM1 = SUM1 + KVEC%KNOTS(I+J)
! 		END DO
! 		SUM1 = SUM1/ KVEC%POLY_ORDER
! 		GREV_PTS(I+1) = SUM1
! 	ENDDO
! 	
! END SUBROUTINE FIND_GREVILLE_PT

END MODULE GEOMETRY
