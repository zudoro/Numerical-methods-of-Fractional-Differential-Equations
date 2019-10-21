	MODULE ASSEMBLE_GAL

	USE GSQUAD
	USE LUDECOMPOSITION
	USE BASIS
	USE LOADFUNCTION
	
	IMPLICIT NONE

CONTAINS



!------------------INTEGRAL CODE FOR THE STIFFNESS MATRIX ELEMENT--------------
SUBROUTINE GEN_KF(SOL)

	REAL*8, INTENT(OUT) :: SOL(DOF)
	
	REAL*8 :: LC_M(LC_NUMBS(1), LC_NUMBS(1)), LC_F(LC_NUMBS(1))
! 	REAL*8 :: DMY_M(LC_NUMBS(1), LC_NUMBS(1))
	
	TYPE(FUNCTION_1D) :: LOCBS(2), LC_UH(NUMIC, LC_NUMBS(1)), EX_SOL
! 	TYPE(FUNCTION_1D) :: PU(2)

	REAL*8 :: FDBS, D, INPROD_IC, INPROD_FD, INPROD_LS
	INTEGER :: I, J, II, JJ, K, KK, G, PATCHROW, PATCHCOLUMN, ROW, COLUMN
	INTEGER :: LC_INDX(LC_NUMBS(1))
	
	REAL*8 :: INNER_PRINT(1:JORDER+1), INNER_FD_PRINT(1:JORDER+1)
	
	LC_M(:,:) = 0.D0; 
	LC_F(:) = 0.D0; 
	SOL(:) = 0.0D0
! 	DMY_M(:,:) = 0.0D0
	
	!! Construct the stiffness matrix M_0 and load vector F_0 corresponding to the first patch Omega_0
	!-------------------------------------------------------------------------------------------------!
	ROW = 0
	DO I = 0, JORDER ! ROW
		ROW = ROW + 1
		COLUMN = 0
		DO J = 0, JORDER ! COLUMN
			COLUMN = COLUMN + 1
			DO K = 1, NUMGSPT(1)
				LOCBS(1) = GET_PHY_JACOBI_POLY(GSPT(1,K), I, 1)	! Row
				LOCBS(2) = GET_PHY_JACOBI_POLY(GSPT(1,K), J, 1)	! Column
! 					FDBS = GET_FD_JACOBI_POLY(GSPT(1,K), J)	!Column
! 									PU(1) = PHYPU1D(GSPT(SUBIRNDX(PATCHROW,PATCHCOLUMN,II),K), PATCHROW)
! 									PU(2) = PHYPU1D(GSPT(SUBIRNDX(PATCHROW,PATCHCOLUMN,II),K), PATCHCOLUMN)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				LC_M(ROW, COLUMN) = LC_M(ROW, COLUMN) + LOCBS(1)%VAL(0)*LOCBS(2)%VAL(0)*GSW(1,K)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			ENDDO
			INPROD_LS = GET_INNER_FD_LS(ROW, COLUMN, 1)
			
			INNER_PRINT(COLUMN) = LC_M(ROW, COLUMN)
			INNER_FD_PRINT(COLUMN) =  INPROD_LS
			
			LC_M(ROW,COLUMN) = LC_M(ROW, COLUMN) + INPROD_LS
		ENDDO
!-------------------------------------------------------------------------------------	!
! 		write(132,*) (inner_print(j), j=1, jorder + 1)
! 		write(135,*) (INNER_FD_PRINT(j), j=1, jorder + 1)
!-------------------------------------------------------------------------------------	!
		! Evaluate the inner product of the right hand side of the Voltera equation
		INPROD_IC = GET_INNER_IC(ROW, 1)
		INPROD_FD = GET_INNER_FD(SOL(:), ROW, 1)
		
		LC_F(ROW) = INPROD_IC + INPROD_FD
!-------------------------------------------------------------------------------------	!
! 		write(133,*) lc_f(row)
!-------------------------------------------------------------------------------------	!
	ENDDO
!-------------------------------------------------------------------------------------	!
! 	write(132,*) ("", j=1, jorder + 1)
! 	write(133,*) ""
! 	write(135,*) ("", j=1, jorder + 1)
!-------------------------------------------------------------------------------------	!
    
    OPEN(1, FILE = './data/f1')
	DO II=1, LC_NUMBS(1)
		WRITE(1,*) LC_F(II)
	ENDDO
	CLOSE(1)
! 
	OPEN(1, FILE = './data/k1')
	DO II=1, LC_NUMBS(1)
		WRITE(1,*) (LC_M(II,JJ), JJ=1,LC_NUMBS(1))
	ENDDO
	CLOSE(1)
	
	!! Copy M_0 & F_0 and then solve the linear system M_0[sol]=F_0
	!---------------------------------------------------------!
	CALL LUDCMP(LC_M, LC_NUMBS(1), LC_NUMBS(1), LC_INDX, D)
	CALL LUBKSB(LC_M, LC_NUMBS(1), LC_NUMBS(1), LC_INDX, LC_F)
	
	SOL(1:LC_NUMBS(1)) = LC_F(:)
	
	OPEN(1, FILE = './data/sol1')
	DO II = 1, LC_NUMBS(1)
		WRITE(1,*) SOL(II)
	ENDDO
	CLOSE(1)
	
!-------------------------------------------------------------------------------------	!	
! 	do i=1, lc_numbs(1)
!       write(134,*) sol(i)
!     enddo
!     write(134,*) ""
!-------------------------------------------------------------------------------------	!
    
    
	IF (NUMPATCH==1) THEN 
		GOTO 17
	ENDIF
	
	!! Construct linear system [M_k][sol]=[F_k] and solve it for each patches Omega_k
	!---------------------------------------------------------------------------------------------------------!
	DO K = 2, NUMPATCH
		!! Update the stiffness matrix M_k
		!---------------------------------------------------------------------------------------------------------!
		LC_M(:,:) = 0.0D0
		LC_F(:) = 0.0D0
		ROW = 0
		DO I = 0, JORDER ! ROW
			ROW = ROW + 1
			COLUMN = 0
			DO J = 0, JORDER ! COLUMN
				COLUMN = COLUMN + 1
				DO G = 1, NUMGSPT(1)
					LOCBS(1) = GET_PHY_JACOBI_POLY(GSPT(K, G), I, K)	! Row
					LOCBS(2) = GET_PHY_JACOBI_POLY(GSPT(K, G), J, K)	! Column
					
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					LC_M(ROW, COLUMN) = LC_M(ROW, COLUMN) + LOCBS(1)%VAL(0)*LOCBS(2)%VAL(0)*GSW(K, G)
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				ENDDO
				INPROD_LS = GET_INNER_FD_LS(ROW, COLUMN, K)
				
				INNER_PRINT(COLUMN) = LC_M(ROW, COLUMN)
                INNER_FD_PRINT(COLUMN) =  INPROD_LS
			
				LC_M(ROW, COLUMN) = LC_M(ROW, COLUMN) + INPROD_LS
			ENDDO
			
! 			write(132,*) (inner_print(j), j=1, jorder + 1)
!             write(135,*) (INNER_FD_PRINT(j), j=1, jorder + 1)
			
			! Evaluate the inner product of the right hand side of the Voltera equation
			INPROD_IC = GET_INNER_IC(ROW, K)
			INPROD_FD = GET_INNER_FD(SOL(:), ROW, K)
			
			LC_F(ROW) = INPROD_IC + INPROD_FD
! 			write(133,*) lc_f(row)
		ENDDO
		
! 		write(132,*) ("", j=1, jorder + 1)
!         write(133,*) ""
!         write(135,*) ("", j=1, jorder + 1)
		
        write(CHAR_PATCH, fmt='(i1)') K
        CHAR_MATRIX = trim('./data/k') // CHAR_PATCH
        CHAR_LOADF = trim('./data/f') // CHAR_PATCH
        CHAR_SOL = trim('./data/sol') // CHAR_PATCH
        OPEN(5, FILE=CHAR_MATRIX, STATUS='UNKNOWN')
        OPEN(6, FILE=CHAR_LOADF, STATUS='UNKNOWN')
        OPEN(7, FILE=CHAR_SOL, STATUS='UNKNOWN')
        
		DO II = 1, LC_NUMBS(K)
			WRITE(6,*) LC_F(II)
            WRITE(5,*) (LC_M(II, JJ), JJ = 1, LC_NUMBS(K))
		ENDDO
		CLOSE(5)
		CLOSE(6)
		
		!! Copy M_k & F_k and then solve the linear system M_k[sol]=F_k
		!---------------------------------------------------------------------------------------------------------!
		
		CALL LUDCMP(LC_M, LC_NUMBS(K), LC_NUMBS(K), LC_INDX, D)
		CALL LUBKSB(LC_M, LC_NUMBS(K), LC_NUMBS(K), LC_INDX, LC_F)
		
		SOL(SUM(LC_NUMBS(1:K-1)) + 1:SUM(LC_NUMBS(1:K))) = LC_F(:)
		
        DO II = 1, LC_NUMBS(K)
			WRITE(7,*) LC_F(II)
		ENDDO
		CLOSE(7)
		!---------------------------------------------------------------------------------------------------------!
	ENDDO
	!---------------------------------------------------------------------------------------------------------!
	
	17 CONTINUE
	
	WRITE(*,*)
	WRITE(*,*) '<<< ASSEMBLE STIFFNESS MATRIX AND LOAD VECTOR : DONE >>>'
	WRITE(*,*)
	
END SUBROUTINE GEN_KF

! Inner product of the function g(t) which is the initial condition, with a test function.
REAL*8 FUNCTION GET_INNER_IC(ROW, PATCH)
	
	INTEGER, INTENT(IN) :: ROW, PATCH
	
	TYPE(FUNCTION_1D) :: LC_BS
	REAL*8 :: INNER_PROD, IC_UPDATED(0:INT_M-1)
	INTEGER :: I, J, K, II, JJ, KK, FCTRY
	
	GET_INNER_IC = 0.0D0
	
	DO K = 0, INT_M-1
		INNER_PROD = 0.0d0
		IC_UPDATED(K) = IC(K+1)

		! Compute the inner product
		DO I = 1, NUMGSPT(1)
			LC_BS = GET_PHY_JACOBI_POLY(GSPT(PATCH, I), ROW-1, PATCH)
			INNER_PROD = INNER_PROD + (GSPT(PATCH, I) - PATCHBDPT(PATCH,1))**(1.0D0*K)*LC_BS%VAL(0)*GSW(PATCH, I)
		ENDDO
		FCTRY = FACTORIAL(K)
		
		GET_INNER_IC = GET_INNER_IC + IC_UPDATED(K)*INNER_PROD/(1.0D0*FCTRY)
	ENDDO
		
END FUNCTION GET_INNER_IC


! Inner product of the Voltera integral with a test function.
REAL*8 FUNCTION GET_INNER_FD(SOL, ROW, PATCH)
	
	REAL*8, INTENT(IN) :: SOL(DOF)
	INTEGER, INTENT(IN) :: PATCH, ROW
	
	TYPE(FUNCTION_1D) :: LC_BS, LM(2) ! Linear Mapping
	
	REAL*8 :: INNER_INT, MEMORY_INT, MEMORY_PATCH_INT(PATCH-1), CURRENT_INT
	INTEGER :: I, J, K, II, JJ, KK, Q, G
	
	GET_INNER_FD = 0.0D0
	MEMORY_INT = 0.0D0
	MEMORY_PATCH_INT(:) = 0.0D0
	CURRENT_INT = 0.0D0
	
	IF (PATCH==1) THEN 
		DO I = 1, NUMGSPT(2) ! Number of Gauss Quadrature points
			INNER_INT = 0.0D0
			LM(1) = REFLINE_TO_PHYLINE(REFGSPT(I), PATCHBDPT(PATCH,1:2))
			LC_BS = GET_PHY_JACOBI_POLY(LM(1)%VAL(0), ROW-1, PATCH)
			DO J = 1, NUMGSPT(3) ! Number of Gauss Jacobi points
				LM(2) = REFLINE_TO_PHYLINE(GJPT(J), (/PATCHBDPT(PATCH,1), LM(1)%VAL(0)/))
				INNER_INT = INNER_INT + LDF1D(LM(2)%VAL(0), SOL(1:LC_NUMBS(PATCH)), PATCH, 'OFF')*GJW(J)
			ENDDO
			GET_INNER_FD = GET_INNER_FD + (0.50D0*(LM(1)%VAL(0) - PATCHBDPT(PATCH,1)))**(NU)*INNER_INT*LC_BS%VAL(0)*REFGSW(I)
		ENDDO
		GET_INNER_FD = (1.0D0/DEXP(GAMMLN(NU)))*0.50D0*(PATCHBDPT(PATCH,2) - PATCHBDPT(PATCH,1))*GET_INNER_FD
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	!
!         DO I = 1, NUMGSPT(1) ! Number of Gauss Quadrature points of the inner product
!             INNER_INT = 0.0D0
!             LC_BS = GET_PHY_JACOBI_POLY(GSPT(PATCH, I), ROW-1, PATCH)
!             CALL GAULEG((GSPT(PATCH, I) - PATCHBDPT(PATCH, 1))**((NU - 1.0D0)/P), 0.0D0, INNER_GSPT(:), INNER_GSW(:), INNER_NUMGSPT)
!             DO J = 1, INNER_NUMGSPT ! Number of Gauss Quadrature points of the Voltera integral
!                 INNER_INT = INNER_INT + LDF1D(GSPT(PATCH, I) - INNER_GSPT(J)**(P/(NU - 1.0D0)), SOL(1:LC_NUMBS(PATCH)), PATCH, 'OFF')*INNER_GSPT(J)**(P + (P/(NU - 1.0D0)) - 1.0D0)*INNER_GSW(J)
!             ENDDO
!             INNER_INT = (P/(1.0D0 - NU))*INNER_INT
!             GET_INNER_FD = GET_INNER_FD + INNER_INT*LC_BS%VAL(0)*GSW(PATCH, I)
!         ENDDO
!         GET_INNER_FD = (1.0D0/DEXP(GAMMLN(NU)))*GET_INNER_FD
	ELSE 
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	!
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	!
! 		DO K = 1, PATCH-1
! 			DO I = 1, NUMGSPT(2) ! Number of Gauss Quadrature points
! 				INNER_INT = 0.0D0
! 				LM(1) = REFLINE_TO_PHYLINE(REFGSPT(I), PATCHBDPT(PATCH,1:2))
! 				LC_BS = GET_PHY_JACOBI_POLY(LM(1)%VAL(0), ROW-1, PATCH)
! 				DO J = 1, NUMGSPT(2) ! Number of Gauss Quadrature points
! 					LM(2) = REFLINE_TO_PHYLINE(REFGSPT(J), PATCHBDPT(K,1:2))
! 					INNER_INT = INNER_INT + (LM(1)%VAL(0) - LM(2)%VAL(0))**(NU - 1.0D0)*LDF1D(LM(2)%VAL(0), SOL(SUM(LC_NUMBS(1:K-1)) + 1 : SUM(LC_NUMBS(1:K))), K, 'FUL')*REFGSW(J)
! 				ENDDO
! 				MEMORY_PATCH_INT(K) = MEMORY_PATCH_INT(K) + INNER_INT*LC_BS%VAL(0)*REFGSW(I)
! 			ENDDO
! 			MEMORY_PATCH_INT(K) = 0.50D0*(PATCHBDPT(K,2) - PATCHBDPT(K,1))*MEMORY_PATCH_INT(K)
! 		ENDDO
! 		MEMORY_INT = (1.0D0/DEXP(GAMMLN(NU)))*0.50D0*(PATCHBDPT(PATCH,2) - PATCHBDPT(PATCH,1))*SUM(MEMORY_PATCH_INT(:))
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	!
        DO K = 1, PATCH - 1
            DO I = 1, NUMGSPT(1) ! Number of Gauss Quadrature points of the inner product
                INNER_INT = 0.0D0
                LC_BS = GET_PHY_JACOBI_POLY(GSPT(PATCH, I), ROW-1, PATCH)
                CALL GAULEG((GSPT(PATCH, I) - PATCHBDPT(K, 1))**((NU - 1.0D0)/P), (GSPT(PATCH, I) - PATCHBDPT(K, 2))**((NU - 1.0D0)/P), INNER_GSPT(:), INNER_GSW(:), INNER_NUMGSPT)
                DO J = 1, INNER_NUMGSPT ! Number of Gauss Quadrature points of the Voltera integral
                    INNER_INT = INNER_INT + LDF1D(GSPT(PATCH, I) - INNER_GSPT(J)**(P/(NU - 1.0D0)), SOL(SUM(LC_NUMBS(1:K-1)) + 1 : SUM(LC_NUMBS(1:K))), K, 'FUL')*INNER_GSPT(J)**(P + (P/(NU - 1.0D0)) - 1.0D0)*INNER_GSW(J)
                ENDDO
                INNER_INT = (P/(1.0D0 - NU))*INNER_INT
                MEMORY_PATCH_INT(K) = MEMORY_PATCH_INT(K) + INNER_INT*LC_BS%VAL(0)*GSW(PATCH, I)
            ENDDO
        ENDDO
        MEMORY_INT = (1.0D0/DEXP(GAMMLN(NU)))*SUM(MEMORY_PATCH_INT(:))
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	!
! 		DO K = 1, PATCH-2
! 			DO I = 1, NUMGSPT(1) ! Number of Gauss Quadrature points
! 				INNER_INT = 0.0D0
! 				LC_BS = GET_PHY_JACOBI_POLY(GSPT(PATCH, I), ROW-1, PATCH)
! 				DO J = 1, NUMGSPT(1) ! Number of Gauss Quadrature points
! 					INNER_INT = INNER_INT + (GSPT(PATCH, I) - GSPT(K, J))**(NU - 1.0D0)*LDF1D(GSPT(K, J), SOL(SUM(LC_NUMBS(1:K-1)) + 1 : SUM(LC_NUMBS(1:K))), K, 'FUL')*GSW(K, J)
! 				ENDDO
! 				MEMORY_PATCH_INT(K) = MEMORY_PATCH_INT(K) + INNER_INT*LC_BS%VAL(0)*GSW(PATCH, I)
! 			ENDDO
! 		ENDDO
! 		
! 	`	K = PATCH-1
! 		DO I = 1, NUMGSPT(4) ! Number of Gauss Quadrature points
! 			INNER_INT = 0.0D0
! 			LC_BS = GET_PHY_JACOBI_POLY(SPGSPT(PATCH, I), ROW-1, PATCH)
! 			DO J = 1, NUMGSPT(4) ! Number of Gauss Quadrature points
! 				INNER_INT = INNER_INT + (SPGSPT(PATCH, I) - SPGSPT(K, J))**(NU - 1.0D0)*LDF1D(SPGSPT(K, J), SOL(SUM(LC_NUMBS(1:K-1)) + 1 : SUM(LC_NUMBS(1:K))), K, 'FUL')*SPGSW(K, J)
! 			ENDDO
! 			MEMORY_PATCH_INT(K) = MEMORY_PATCH_INT(K) + INNER_INT*LC_BS%VAL(0)*SPGSW(PATCH, I)
! 		ENDDO
! 		
! 		MEMORY_INT = (1.0D0/DEXP(GAMMLN(NU)))*SUM(MEMORY_PATCH_INT(:))
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	!
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	!
		DO I = 1, NUMGSPT(2) ! Number of Gauss Quadrature points
			INNER_INT = 0.0D0
			LM(1) = REFLINE_TO_PHYLINE(REFGSPT(I), PATCHBDPT(PATCH,1:2))
			LC_BS = GET_PHY_JACOBI_POLY(LM(1)%VAL(0), ROW-1, PATCH)
			DO J = 1, NUMGSPT(3) ! Number of Gauss Jacobi points
				LM(2) = REFLINE_TO_PHYLINE(GJPT(J), (/PATCHBDPT(PATCH,1), LM(1)%VAL(0)/))
				INNER_INT = INNER_INT + LDF1D(LM(2)%VAL(0), SOL(SUM(LC_NUMBS(1:PATCH-1)) + 1 : SUM(LC_NUMBS(1:PATCH))), PATCH, 'OFF')*GJW(J)
			ENDDO
			CURRENT_INT = CURRENT_INT + (0.50D0*(LM(1)%VAL(0) - PATCHBDPT(PATCH,1)))**(NU)*INNER_INT*LC_BS%VAL(0)*REFGSW(I)
		ENDDO
		CURRENT_INT = (1.0D0/DEXP(GAMMLN(NU)))*0.50D0*(PATCHBDPT(PATCH,2) - PATCHBDPT(PATCH,1))*CURRENT_INT
		GET_INNER_FD = MEMORY_INT + CURRENT_INT
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	!
!         DO I = 1, NUMGSPT(1) ! Number of Gauss Quadrature points of the inner product
!             INNER_INT = 0.0D0
!             LC_BS = GET_PHY_JACOBI_POLY(GSPT(PATCH, I), ROW-1, PATCH)
!             CALL GAULEG((GSPT(PATCH, I) - PATCHBDPT(PATCH, 1))**((NU - 1.0D0)/P), 0.0D0, INNER_GSPT(:), INNER_GSW(:), INNER_NUMGSPT)
!             DO J = 1, INNER_NUMGSPT ! Number of Gauss Quadrature points of the Voltera integral
!                 INNER_INT = INNER_INT + LDF1D(GSPT(PATCH, I) - INNER_GSPT(J)**(P/(NU - 1.0D0)), SOL(1:LC_NUMBS(PATCH)), PATCH, 'OFF')*INNER_GSPT(J)**(P + (P/(NU - 1.0D0)) - 1.0D0)*INNER_GSW(J)
!             ENDDO
!             INNER_INT = (P/(1.0D0 - NU))*INNER_INT
!             CURRENT_INT = CURRENT_INT + INNER_INT*LC_BS%VAL(0)*GSW(PATCH, I)
!         ENDDO
!         CURRENT_INT = (1.0D0/DEXP(GAMMLN(NU)))*CURRENT_INT
        
        GET_INNER_FD = MEMORY_INT + CURRENT_INT
        
	ENDIF
	
END FUNCTION GET_INNER_FD


! Inner Product of the term included in the Volteral integral involving unknown function u(t), with test function
REAL*8 FUNCTION GET_INNER_FD_LS(ROW, COLUMN, PATCH) ! Inner product of the Volteral integral on the left side of the variational formulation
	
	INTEGER, INTENT(IN) :: PATCH, ROW, COLUMN
	
	TYPE(FUNCTION_1D) :: LC_BS(2), LM(2) ! Linear Mapping, LC_BS(1) = test function (row), LC_BS(2) = trial function (column)
	
	REAL*8 :: INNER_INT
	INTEGER :: I, J
	
	GET_INNER_FD_LS = 0.0D0
	
	DO I = 1, NUMGSPT(2) ! Number of Gauss Quadrature points
		INNER_INT = 0.0D0
		LM(1) = REFLINE_TO_PHYLINE(REFGSPT(I), PATCHBDPT(PATCH,1:2))
		LC_BS(1) = GET_PHY_JACOBI_POLY(LM(1)%VAL(0), ROW-1, PATCH)
		DO J = 1, NUMGSPT(3) ! Number of Gauss Jacobi points
			LM(2) = REFLINE_TO_PHYLINE(GJPT(J), (/PATCHBDPT(PATCH,1), LM(1)%VAL(0)/))
			LC_BS(2) = GET_PHY_JACOBI_POLY(LM(2)%VAL(0), COLUMN-1, PATCH)
			
			IF (PROBLEM==2) THEN 
				INNER_INT = INNER_INT + LC_BS(2)%VAL(0)*GJW(J)
			ENDIF
! 			print*, inner_int
		ENDDO
		GET_INNER_FD_LS = GET_INNER_FD_LS + (0.50D0*(LM(1)%VAL(0) - PATCHBDPT(PATCH,1)))**(NU)*INNER_INT*LC_BS(1)%VAL(0)*REFGSW(I)
	ENDDO
	GET_INNER_FD_LS = (1.0D0/DEXP(GAMMLN(NU)))*0.50D0*(PATCHBDPT(PATCH,2) - PATCHBDPT(PATCH,1))*GET_INNER_FD_LS
		
END FUNCTION GET_INNER_FD_LS

END MODULE ASSEMBLE_GAL
