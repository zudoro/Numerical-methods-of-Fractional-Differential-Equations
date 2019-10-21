PROGRAM MAIN

  USE GEOMETRY
  USE PLOT
  USE ERRORESTIMATE
	USE ASSEMBLE_GAL
	USE LUDECOMPOSITION
	
  IMPLICIT NONE

	REAL*8, ALLOCATABLE :: SOL(:)

  REAL*8 :: MAX_NORM(2), H_NORM(2,2), NU_ARRAY(8), AB_ARRAY(4,2)
  REAL :: CPTS, CPTE
 	INTEGER :: I, J, K, II, JJ, NP_ARRAY(1)
	
! 	NP_ARRAY = (/1, 2, 4, 8, 16, 32, 64, 128/)
! 	NP_ARRAY = (/2, 8, 16, 32, 64, 128/)
! 	NP_ARRAY = (/1, 3, 9, 27, 81/)
! 	NP_ARRAY(1) = 2
! 	NP_ARRAY = (/5, 10, 15, 20, 30/)
	NP_ARRAY = (/32/)
	
	NU_ARRAY = (/0.20D0, 0.40D0, 0.60D0, 0.80D0, 1.20D0, 1.40D0, 1.60D0, 1.80D0/)
	
	AB_ARRAY(1,:) = (/0.0D0, 0.0D0/)
	AB_ARRAY(2,:) = (/1.0D0, 1.0D0/)
	AB_ARRAY(3,:) = (/0.50D0, 0.50D0/)
	AB_ARRAY(4,:) = (/-0.50D0, -0.50D0/)
	
!!---------------------------------!!
	DO K = 1, 4
		ALPHA = AB_ARRAY(K,1)
		BETA = AB_ARRAY(K,2)
		
		print*, alpha, beta
		
		write(char_problem,fmt='(i1)') PROBLEM
	! 	write(char_puorder,fmt='(i1)') PUORDER
		write(char_jorder,fmt='(i1)') K
		write(char_numpatch,fmt='(i3)') NUMPATCH	
		
		write(char_j,fmt='(i2)') J

! 		FILENAME = trim('./data/test') // char_problem // trim('k') // char_numpatch
! 		FILENAME = trim('./data/test') // char_problem // trim('k') // char_numpatch // trim('j') // char_jorder
		FILENAME = trim('./data/test') // char_problem // trim('ab') // char_jorder
! 		FILENAME = trim('./data/test') //  char_problem // char_j
		
		OPEN(131, FILE=FILENAME, STATUS='UNKNOWN')
		
!!---------------------------------!!
! 	DO J = 1, UBOUND(NP_ARRAY,1)
! 		NUMPATCH = NP_ARRAY(J)
!!---------------------------------!!
! 	DO I = 2, 8, 2
! 		JORDER = I
		
!!---------------------------------!!
! 	DO J = 1, UBOUND(NU_ARRAY,1)
	DO J = 1, 8
		NU = NU_ARRAY(J)
		
		IF (NU<=1.0D0) THEN 
			M = 1.0D0
			INT_M = 1
			NUMIC = 1
		ELSEIF (NU>1.0D0 .AND. NU<=2.0D0) THEN 
			M = 2.0D0
			INT_M = 2
			NUMIC = 2
		ENDIF
		ALLOCATE(IC(NUMIC))
!!---------------------------------!!


! 		write(char_problem,fmt='(i1)') PROBLEM
! 	! 	write(char_puorder,fmt='(i1)') PUORDER
! 		write(char_jorder,fmt='(i1)') JORDER
! 		write(char_numpatch,fmt='(i3)') NUMPATCH	
! 		
! 		write(char_j,fmt='(i2)') J
! 
! ! 		FILENAME = trim('./data/test') // char_problem // trim('k') // char_numpatch
! ! 		FILENAME = trim('./data/test') // char_problem // trim('k') // char_numpatch // trim('j') // char_jorder
! ! 		FILENAME = trim('./data/test') // char_problem // trim('j') // char_jorder
! 		FILENAME = trim('./data/test') //  char_problem // char_j
! 		
! 		OPEN(131, FILE=FILENAME, STATUS='UNKNOWN')
! 		!--------------------------------------------------------------------------------------------------------------------------------------------------------------------!
! 		write(131,*) 'Example 1. in Master Thesis, Multistage Shifted Jacobi Spectral Method for solving linear & nonlinear Fractional Differential Equations'
! 		write(131,*)
! 		write(131,*) 'Alpha=', alpha
! 		write(131,*) 'Beta=', beta
! 		write(131,*) 'Nu=', nu
! 		WRITE(131,*) ""
! 		write(131,*) '				', 'degree', '	', 'number of patches', '	', 'dof', '			', 'Rel. Max. Norm', '						', 'Rel. L2 Norm', '									', 'Rel. H1 Norm'
! 		!--------------------------------------------------------------------------------------------------------------------------------------------------------------------!

!!---------------------!!
		DO I = 2, 2
			JORDER = I
!!---------------------!!
! 		DO J = 1, UBOUND(NP_ARRAY,1)
! 			NUMPATCH = NP_ARRAY(J)
!!---------------------!!
! 		DO J = 1, 1
			NUMPATCH = NP_ARRAY(1)
!!---------------------!!
			ALLOCATE(PATCHBDPT(NUMPATCH, 2), LC_NUMBS(NUMPATCH))
			
	! 		EXTRA_KNOTS = (/I,0/)
	! 		NURBS_REGULARITY = (/BS_ORDER(1)-1, BS_ORDER(2)-1/)

! 			CALL CPU_TIME(CPTS)

			CALL GET_GEO()

			ALLOCATE(SOL(DOF))

			CALL PLOTGM()

			CALL GEN_KF(SOL)
			
! 			do ii = 1, dof
! 				print*, sol(ii)
! 			enddo
			
			CALL MAXNORM(SOL, MAX_NORM)
			CALL H_NORM_ERROR(SOL, H_NORM)
			
! 			CALL CPU_TIME(CPTE)

			PRINT*, JORDER, NUMPATCH, DOF, MAX_NORM(1), H_NORM(1,1), H_NORM(2,1), MAX_NORM(2), H_NORM(1,2), H_NORM(2,2)
			WRITE(131, *) JORDER, NUMPATCH, DOF, MAX_NORM(1), H_NORM(1,1), H_NORM(2,1), MAX_NORM(2), H_NORM(1,2), H_NORM(2,2)

! 			PRINT*, ""
! 			PRINT*, "ELAPSED CPU TIME : ", CPTE - CPTS
! 			PRINT*, ""
			
			DEALLOCATE(SOL)
			DEALLOCATE(NDX)
			DEALLOCATE(GSPT, GSW)
			DEALLOCATE(PATCHBDPT, LC_NUMBS)
	! 		DEALLOCATE(SUBIRNDX)
	! 		DEALLOCATE(JETA)
		ENDDO
		
! 		write(131,*) ""
! 		write(131,*) 'Step size=', deltah
		
! 		CLOSE(131)
		DEALLOCATE(IC)
	ENDDO
	CLOSE(131)
	ENDDO
	
END PROGRAM MAIN
