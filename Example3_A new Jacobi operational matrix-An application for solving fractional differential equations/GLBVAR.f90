! NAME : HYUNJU KIM
! DATE : 11 / 8 / 2010

MODULE GLBVAR

	USE NEWTYPE

  IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   GLOBAL PARAMETER

	REAL(8), PARAMETER :: PI=DACOS(-1.0D0)
	REAL(8), PARAMETER :: DEGREE = PI/180.D0
	REAL(8), PARAMETER :: EPS=5.0D-15
	REAL(8), PARAMETER :: DELTA = 0.10D0
	INTEGER, PARAMETER :: PUORDER = 3
	
	! Problem number & domain
	INTEGER, PARAMETER :: PROBLEM = 6
	
	REAL(8), PARAMETER :: LAMBDA = PI ! PROBLEM 4
	INTEGER, PARAMETER :: FOURIER_TERM = 3 ! PROBLEM 5
	INTEGER, PARAMETER :: NUM_MF = 10**(3) ! PROBLEM 6
	
	REAL(8), DIMENSION(1:2), PARAMETER :: DOMAIN = (/0.0D0, 1.0D0/)

	! The order of fractional derivative and its ceiling number
	REAL*8 :: NU
	REAL*8 :: M
	INTEGER :: INT_M
	
	! PARAMETERS OF JACOBI POLYNOMIAL
	! Legendre polynomial basis functions
	REAL*8 :: ALPHA
	REAL*8 :: BETA
	
	! Chebyshev polynomial basis functions
! 	REAL*8, PARAMETER :: ALPHA = -0.50D0
! 	REAL*8, PARAMETER :: BETA  = -0.50D0
	
! 	REAL*8, PARAMETER :: ALPHA = 1.0D0
! 	REAL*8, PARAMETER :: BETA  = 1.0D0
	
! 	REAL*8, PARAMETER :: ALPHA = 0.50D0
! 	REAL*8, PARAMETER :: BETA  = 0.50D0
	
	! Number of initial conditions
	INTEGER :: NUMIC
	REAL*8, ALLOCATABLE :: IC(:)
	
	INTEGER :: NUMPATCH

	REAL*8, ALLOCATABLE :: PATCHBDPT(:, :)

	REAL*8 :: DELTAH
	
	! order of the shifted Jacobi polynomials
	INTEGER :: JORDER
	
	! order of the Lagrange interpolating function
	INTEGER :: IORDER
	
	!! variables used to construct Berstain polynomials
	TYPE(KNOT_VECTOR) :: BASIS_KVEC
	!! GLOBAL INDEX
	INTEGER :: DOF
	INTEGER, ALLOCATABLE :: LC_NUMBS(:)
! 	TYPE(BD_INFO) :: BDNDX(MAX_LENGTH)

	!! LOCAL BASIS INDEX
	INTEGER, ALLOCATABLE :: NDX(:,:)

	!! SUB-LIMITS OF INTEGRAL INDEX
! 	REAL*8, ALLOCATABLE :: JETA(:,:)
! 	INTEGER, ALLOCATABLE :: SUBIRNDX(:,:,:)
	
		! GAUSS POINTS
	INTEGER :: NUMGSPT
	REAL*8, ALLOCATABLE :: GSPT(:,:), GSW(:,:)
	
	!! SOLUTION ARRAY IN THE LEAST SQUARE
! 	REAL*8, ALLOCATABLE :: LST_SOL(:)

	!!  BINOMIAL COEFFICIENTS
	REAL(8) :: BINOM(MAX_DIFF_ORDER,0:MAX_DIFF_ORDER)
	
	!!  CONNECTIVITY ARRAY
! 	INTEGER, ALLOCATABLE :: CARRAY(:,:)
	
	CHARACTER(LEN=1) :: CHAR_PUORDER
	CHARACTER(LEN=1) :: CHAR_JORDER
	CHARACTER(LEN=3) :: CHAR_NUMPATCH
	CHARACTER(LEN=1) :: CHAR_PROBLEM
	CHARACTER(LEN=2) :: CHAR_J
	CHARACTER(LEN=100) :: FILENAME
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	MATERIAL COEFFICIENTS

	INTEGER, DIMENSION(1:2), PARAMETER :: ORDER1=(/ 1, 2 /)
	INTEGER, DIMENSION(1:2), PARAMETER :: ORDER2=(/ 2, 1 /)

	TYPE(POINT2D), PARAMETER :: PT0 = POINT2D(0.D0,0.D0)

END MODULE GLBVAR
