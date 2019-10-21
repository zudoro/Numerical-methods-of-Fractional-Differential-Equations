PROGRAM MAIN

  USE GEOMETRY
  USE PLOT
  USE BOUNDARY
  USE ERRORESTIMATE
	USE ASSEMBLE_GAL
	USE LUDECOMPOSITION
	
  IMPLICIT NONE	

	REAL*8 :: MAX_NORM(2), H_NORM(2,2)
	REAL :: CPTS, CPTE
 	INTEGER :: I, J, II, JJ, NP_ARRAY(1)
	
	REAL*8 :: L2ERROR_H, L2ERROR_HALF_H
	
! 	NP_ARRAY = (/2, 4, 8, 16, 32, 64, 128, 256/)
! 	NP_ARRAY(1) = 16
	NP_ARRAY = (/3/)

!!---------------------------------!!
! 	DO J = 1, UBOUND(NP_ARRAY,1)
! 		NUMPATCH = NP_ARRAY(J)
!!---------------------------------!!
	DO I = 2, 2
		JORDER = I
!!---------------------------------!!

		write(char_problem,fmt='(i1)') PROBLEM
	! 	write(char_puorder,fmt='(i1)') PUORDER
		write(char_jorder,fmt='(i1)') JORDER

! 		FILENAME = trim('./data/test') // char_problem // trim('k') // char_numpatch
! 		FILENAME = trim('./data/test') // char_problem // trim('k') // char_numpatch // trim('j') // char_jorder
		FILENAME = trim('./data/test') // char_problem // trim('j') // char_jorder

		OPEN(131, FILE=FILENAME, STATUS='UNKNOWN')

!!---------------------!!
! 		DO I = 2, 8, 2
! 			JORDER = I
!!---------------------!!
		DO J = 1, UBOUND(NP_ARRAY,1)
			NUMPATCH = NP_ARRAY(J)
!!---------------------!!

		write(char_numpatch,fmt='(i3)') NUMPATCH
		write(char_patch,fmt='(i3)') J
		
		STIFNAME = trim('./data/stif1') // char_problem // trim('jorder') // char_jorder // trim("patch") // char_patch
		LOADNAME = trim('./data/load') // char_problem // trim('jorder') // char_jorder // trim("patch") // char_patch
		SOLNAME = trim('./data/sol') // char_problem // trim('jorder') // char_jorder // trim("patch") // char_patch
		EXTRANAME = trim('./data/stif2') // char_problem // trim('jorder') // char_jorder // trim("patch") // char_patch
		
		OPEN(132, FILE=STIFNAME, STATUS='UNKNOWN')
		OPEN(133, FILE=LOADNAME, STATUS='UNKNOWN')
		OPEN(134, FILE=SOLNAME, STATUS='UNKNOWN')
		OPEN(135, FILE=EXTRANAME, STATUS='UNKNOWN')
		
			ALLOCATE(PATCHBDPT(NUMPATCH, 2), LC_NUMBS(NUMPATCH))

! 			CALL CPU_TIME(CPTS)

			CALL GET_GEO()
			
			PRINT*, 'DOF=',DOF
			print*, 'NUMPATCH=', NUMPATCH
			
			ALLOCATE(SOL(DOF))

			CALL PLOTGM()

			CALL GEN_KF(SOL)
			
			CALL MAXNORM(SOL, MAX_NORM)
			CALL H_NORM_ERROR(SOL, H_NORM)
			
			IF (J==1) THEN 
                L2ERROR_H = 0.0D0
                L2ERROR_HALF_H = H_NORM(1,1)
            ELSE
                L2ERROR_H = L2ERROR_HALF_H
                L2ERROR_HALF_H = H_NORM(1,1)
            ENDIF
! 			CALL CPU_TIME(CPTE)

			PRINT*, JORDER, NUMPATCH, DOF, MAX_NORM(1), H_NORM(1,1), H_NORM(2,1), MAX_NORM(2), H_NORM(1,2), H_NORM(2,2), DLOG(L2ERROR_H/L2ERROR_HALF_H)/DLOG(2.0D0)
			WRITE(131, *) JORDER, NUMPATCH, DOF, MAX_NORM(1), H_NORM(1,1), H_NORM(2,1), MAX_NORM(2), H_NORM(1,2), H_NORM(2,2), DLOG(L2ERROR_H/L2ERROR_HALF_H)/DLOG(2.0D0)

! 			PRINT*, ""
! 			PRINT*, "ELAPSED CPU TIME : ", CPTE - CPTS
! 			PRINT*, ""
			
			DEALLOCATE(SOL)
			DEALLOCATE(NDX)
			DEALLOCATE(GSPT, GSW, GJPT, GJW, REFGSPT, REFGSW, SPGSPT, SPGSW, INNER_GSPT, INNER_GSW)
			DEALLOCATE(PATCHBDPT, LC_NUMBS)
	! 		DEALLOCATE(SUBIRNDX)
	! 		DEALLOCATE(JETA)
            
            CLOSE(132); CLOSE(133); CLOSE(134)
		ENDDO
		
! 		write(131,*) ""
! 		write(131,*) 'Step size=', deltah
		
		CLOSE(131)
	ENDDO
	
END PROGRAM MAIN
