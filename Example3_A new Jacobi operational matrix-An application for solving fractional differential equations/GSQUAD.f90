MODULE GSQUAD

	implicit integer (i-n)
	implicit real(8) (a-h,o-z)

CONTAINS

SUBROUTINE GAULEG(X1,X2,X,W,N)
	INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: X1, X2
	REAL*8, INTENT(INOUT) :: X(N), W(N)
	REAL*8 :: PI, P1,P2,P3,PP,XL,XM,Z,Z1
  REAL*8, PARAMETER :: EPS=5.0D-16
  INTEGER :: I,J,M

	PI = DACOS(-1.0D0)
	M=(N+1)/2
	XM=0.5D0*(X2+X1)
	XL=0.5D0*(X2-X1)
	DO I = 1, M
		Z=DCOS(PI*(I-0.25D0)/(N+0.5D0))
1   CONTINUE
		P1=1.0D0
		P2=0.0D0
		DO J = 1, N
			P3=P2
			P2=P1
			P1=((2.0D0*J-1.0D0)*Z*P2-(J-1.0D0)*P3)/J
		ENDDO
		PP=N*(Z*P1-P2)/(Z*Z-1.0D0)
		Z1=Z
		Z=Z1-P1/PP
		IF(DABS(Z-Z1).GT.EPS)GOTO 1
		X(I)=XM-XL*Z
		X(N+1-I)=XM+XL*Z
		W(I)=2.0D0*XL/((1.0D0-Z*Z)*PP*PP)
		W(N+1-I)=W(I)
	ENDDO
END SUBROUTINE GAULEG


! Extended Trapezoidal Rule
SUBROUTINE TRAPZD(FUNC,A,B,S,N)

	IMPLICIT NONE

	REAL(8), INTENT(IN) :: A,B
	REAL(8), INTENT(INOUT) :: S
	INTEGER, INTENT(IN) :: N

	INTERFACE
		FUNCTION FUNC(X)
		
		REAL(8), INTENT(IN) :: X
		REAL(8) :: FUNC
		
		END FUNCTION FUNC
	END INTERFACE

	REAL(8) :: DEL,FSUM, XX
	INTEGER :: IT, I, J, II, JJ, K, KK
	
	S = 0.0D0
	FSUM = 0.0D0
	
	IF (N == 1) THEN
		S = 0.50D0*(B-A)*(FUNC(A) + FUNC(B))
	ELSE
		IT = 2.0D0**(N-2)
		DEL = (B-A)/IT	! This Is The Spacing Of The Points To Be Added.
		XX = A + 0.50D0*DEL
		
		DO J = 1, INT(IT)
			FSUM = FSUM + FUNC(XX)
			XX = XX + DEL
		ENDDO
		
		S = 0.50D0*(S+DEL*FSUM)	!This Replaces S By Its Refined Value.
	END IF

	END SUBROUTINE TRAPZD

! REAL*8 FUNCTION TRAPZD_FUNC(X)
! 	
! 	REAL*8, INTENT(IN) :: X
! 	
! 	TRAPZD_FUNC = 
! END FUNCTION TRAPZD_FUNC


SUBROUTINE POLINT(XA, YA, N, X, Y, DY)

!--------------------------------------------------------------------------------------------!
! Given arrays xa and ya , each of length n , and given a value x , this routine returns a
! value y , and an error estimate dy . If P (x) is the polynomial of degree N − 1 such that
! P ( xa i ) = ya i , i = 1, . . . , n , then the returned value y = P ( x ).
!--------------------------------------------------------------------------------------------!

	INTEGER, INTENT(IN) :: N
	REAL*8, INTENT(IN) :: X, XA(N), YA(N)
	
	REAL*8, INTENT(INOUT) :: DY
	REAL*8, INTENT(OUT) :: Y
	
	INTEGER, PARAMETER :: NMAX = 10
	
	INTEGER :: I, MM, NS
	REAL*8 :: DEN, DIF, DIFT, HO, HP, W, C(NMAX), D(NMAX)
	
	NS = 1
	DIF = DABS(X - XA(1))
	
	DO I = 1, N
		DIFT = DABS(X - XA(I))
		IF (DIFT<DIF) THEN
			NS = I
			DIF = DIFT
		ENDIF
		C(I) = YA(I)
		D(I) = YA(I)
	ENDDO
	
	Y = YA(NS)
	NS = NS - 1
	
	DO MM = 1, N-1
		DO I = 1, N - MM
			HO = XA(I) - X
			HP = XA(I+MM) - X
			W = C(I+1) - D(I)
			DEN = HO - HP
			IF (DABS(DEN)<=EPS) THEN
				PRINT*, 'failure in polint' ! This error can occur only if two input xa’s are (to within roundoff) identical.
				STOP 
			ENDIF
			DEN = W/DEN
			D(I) = HP*DEN
			C(I) = HO*DEN
		ENDDO
		IF (2*NS<(N-MM)) THEN
			DY = C(NS + 1)
		ELSE 
			DY = D(NS)
			NS = NS - 1
		ENDIF
		Y = Y + DY
	ENDDO
	
END SUBROUTINE POLINT


REAL*8 FUNCTION GAMMLN(XX)
	
	! Returns The Value Ln[Γ( Xx )] For XX > 0.
	
	REAL*8, INTENT(IN) :: XX
	
	INTEGER :: J
	REAL*8 :: SER,STP,TMP,X,Y,COF(6)
	
! 	Internal Arithmetic Will Be Done In Double Precision, A Nicety That You Can Omit If Five-Figure
! 	Accuracy Is Good Enough.
	
	COF = (/ 76.180091729471460D0,-86.505320329416770D0, 24.014098240830910D0,-1.2317395724501550D0,0.12086509738661790D-2, -0.53952393849530D-5 /)
	STP = 2.50662827463100050D0
	
	X=XX
	Y=X
	TMP=X+5.50D0
	TMP=(X+0.50D0)*DLOG(TMP)-TMP
	SER=1.0000000001900150D0
	
	DO J=1,6
		Y=Y+1.0D0
		SER=SER+COF(J)/Y
	ENDDO
	
	GAMMLN = TMP + DLOG(STP*SER/X)

END FUNCTION GAMMLN



!gamma(alpha)=dexp(LogGamma(alpha))

RECURSIVE FUNCTION LOGGAMMA(X) &
          RESULT (LOGGAMMA_RESULT)


! LOG OF GAMMA FUNCTION BY ASYMPTOTIC SERIES

        REAL*8 X,XTEMP,RECXS,SUM,TERM,X10, &
          HLFLN2PI,B(0:5),LOGGAMMA_RESULT
        INTEGER XINT,K,X11

        DATA HLFLN2PI / 0.918938533205D00/, &
          B / 0.0D00,0.0833333333333D00, &
        -2.77777777778D-03,7.93650793651D-04, &
        -5.95238095238D-04,8.41750841751D-04 /

        IF ( X > 10.0D00 ) THEN
!         USE ASYMPTOTIC SERIES

          RECXS = 1.0D00/(X*X)
          TERM = X
          SUM = (X-0.5D00)*DLOG(X)-X+HLFLN2PI
          DO  K = 1,5
            TERM = TERM*RECXS
            SUM = SUM+B(K)*TERM
          END DO
          LOGGAMMA_RESULT = SUM
          RETURN
        END IF

        IF ( X > 0.0D00 ) THEN
!         RECURRENCE TO  X > 10
          XTEMP = X - 1.0D00
          X11 = 11.0D00 - XTEMP
          X10 = XTEMP + X11
          XINT = XTEMP
          SUM = 0.0D00
          DO K = 1, X11-1
            SUM = SUM - DLOG(X10 - K)
          END DO
          LOGGAMMA_RESULT = SUM+LOGGAMMA(X10)
        RETURN
        END IF

        IF ( X == 0.0D00 ) THEN
          WRITE(*,*) &
           "!!X=0 IN LOGGAMMA; ZERO RETURNED"
          LOGGAMMA_RESULT = 0.0D00
          RETURN
        END IF

        WRITE(*,*) &
         "!!X < 0 IN LOGGAMMA; ZERO RETURNED"
        WRITE(*,*) "USE REFLECTION FORMULA"
        LOGGAMMA_RESULT = 0.0D00
        RETURN

END FUNCTION LOGGAMMA



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Evaluate an incomplete beta function 
!! betai(a,b,x)=\int_0^x t^{a-1} (1-t)^{b-1} dt

FUNCTION betai(a,b,x) 

	REAL*8, INTENT(IN) :: A,B,X 

	REAL*8 :: BETAI, BT

	IF (DABS(X)<=EPS .OR. DABS(X-1.0D0)<=EPS) THEN
		BT=0.D0 
	ELSE
		BT=(X**A)*(1.0D0-X)**B 
	END IF 

	IF (X < (A+1.0D0)/(A+B+2.D0)) THEN
		BETAI = BT*BETACF(A,B,X)/A 
	ELSE

		BETAI = DEXP(GAMMLN(A)+GAMMLN(B)-GAMMLN(A+B))-BT*BETACF(B,A,1.0D0-X)/B
! 		BETAI = GAMMA(A)*GAMMA(B)/GAMMA(A+B)-BT*BETACF(B,A,1.0D0-X)/B 
	END IF

END FUNCTION betai


FUNCTION BETACF(A,B,X) 

	REAL*8 , INTENT(IN) :: A,B,X 

	REAL*8 :: BETACF 

	INTEGER, PARAMETER :: MAXIT=100 
	REAL*8, PARAMETER :: EPSI=EPSILON(X), FPMIN=TINY(X)/EPSI

	REAL*8:: AA,C,D,DEL,H,QAB,QAM,QAP 
	INTEGER :: M,M2 

	QAB=A+B 
	QAP=A+1.0D0 
	QAM=A-1.0D0 
	C=1.0d0
	D=1.0D0-QAB*X/QAP 
	IF (DABS(D) < FPMIN) D=FPMIN 
	D=1.0D0/D 
	H=D 

	DO M=1,MAXIT
		M2=2*M 
		AA=M*1.0D0*(B-M*1.0D0)*X/((QAM+M2*1.0D0)*(A+M2*1.0D0)) 
		D=1.0D0+AA*D 
		IF (DABS(D) < FPMIN) D=FPMIN 
		C=1.0D0+AA/C 
		IF (DABS(C) < FPMIN) C=FPMIN 
		D=1.0D0/D 
		H=H*D*C 
		AA=-(A+M*1.0D0)*(QAB+M*1.0D0)*X/((A+M2*1.0D0)*(QAP+M2*1.0D0)) 
		D=1.0D0+AA*D 
		IF (DABS(D) < FPMIN) D=FPMIN 
		C=1.0D0+AA/C 
		IF (DABS(C) < FPMIN) C=FPMIN 
		D=1.0D0/D 
		DEL=D*C 
		H=H*DEL 
		IF (DABS(DEL-1.0D0) <= EPSI) EXIT
	ENDDO 

	BETACF=H

END FUNCTION BETACF


REAL*8 FUNCTION INCOB(A,B,X)

!       ========================================================
!       Purpose: Compute the incomplete beta function Ix(a,b)
!       Input :  a --- Parameter
!                b --- Parameter
!                x --- Argument ( 0 ף x ף 1 )
!       Output:  INCOB --- Ix(a,b)
!       Routine called: BETA for computing beta function B(p,q)
!       ========================================================

	REAL*8, INTENT(IN) :: A, B, X

	REAL*8 :: DK(51),FK(51), S0, BT
	REAL*8 :: T1, T2, TA, TB
	INTEGER :: K

	S0=(A+1.0D0)/(A+B+2.0D0)
	BT = DEXP(GAMMLN(A))*DEXP(GAMMLN(B))/DEXP(GAMMLN(A+B))
	
	IF (X.LE.S0) THEN
		DO K=1,20
			DK(2*K)=K*1.0D0*(B-K*1.0D0)*X/(A+2.0D0*K-1.0D0)/(A+2.0D0*K)
		ENDDO
		DO K=0,20
			DK(2*K+1)=-(A+K*1.0D0)*(A+B+K*1.0D0)*X/(A+2.D0*K)/(A+2.0D0*K+1.0D0)
		ENDDO
		T1=0.0D0
		DO K=20,1,-1
			T1=DK(K)/(1.0D0+T1)
		ENDDO
		TA=1.0D0/(1.0D0+T1)
		INCOB=X**A*(1.0D0-X)**B/(A)*TA
	ELSE
		DO K=1,20
			FK(2*K)=K*1.0D0*(A-K*1.0D0)*(1.0D0-X)/(B+2.0D0*K-1.00D0)/(B+2.0D0*K)
		ENDDO
		DO K=0,20
			FK(2*K+1)=-(B+K)*(A+B+K)*(1.0D0-X)/(B+2.D0*K)/(B+2.0D0*K+1.0D0)
		ENDDO
		T2=0.0D0
		DO K=20,1,-1
			T2=FK(K)/(1.0D0+T2)
		ENDDO
		TB=1.0D0/(1.0D0+T2)
		INCOB=1.0D0-X**A*(1.0D0-X)**B/(B)*TB
	ENDIF

END FUNCTION INCOB




REAL*8 FUNCTION BETA_FT(P,Q)

!       ==========================================
!       Purpose: Compute the beta function B(p,q)
!       Input :  p --- Parameter  ( p > 0 )
!                q --- Parameter  ( q > 0 )
!       Output:  BT --- B(p,q)
!       Routine called: GAMMA for computing ג(x)
!       ==========================================
	
	REAL*8, INTENT(IN) :: P, Q
	
	REAL*8 :: GP, GQ, GPQ
	
	GP = GAMMA(P)
	GQ = GAMMA(Q)
	GPQ = GAMMA(P + Q)
	BETA_FT = GP*GQ/GPQ

END FUNCTION BETA_FT




REAL*8 FUNCTION GAMMA(X)

!       ==================================================
!       Purpose: Compute gamma function ג(x)
!       Input :  x  --- Argument of ג(x)
!                       ( x is not equal to 0,-1,-2,תתת)
!       Output:  GAMMA --- ג(x)
!       ==================================================

	REAL*8, INTENT(IN) :: X
	
	REAL*8 :: G(26), GR
	REAL*8 :: Z, R
	INTEGER :: M, M1
	INTEGER :: I, J, K

	IF (X.EQ.INT(X)) THEN
		IF (X > EPS) THEN
			GAMMA = 1.0D0
			M1 = INT(X) - 1
			DO K = 2, M1
				GAMMA = GAMMA*K*1.0D0
			ENDDO
		ELSE
				GAMMA = 1.0D+300
		ENDIF
	ELSE
		IF (DABS(X).GT.1.0D0) THEN
			Z = DABS(X)
			M = INT(Z)
			R = 1.0D0
			DO K = 1, M
				R = R*(Z-K*1.0D0)
			ENDDO
			Z = Z - M*1.0D0
		ELSE
			Z = X
		ENDIF
		
		DATA G /1.0D0,0.5772156649015329D0, &
          -0.6558780715202538D0, -0.420026350340952D-1, &
          0.1665386113822915D0,-0.421977345555443D-1, &
          -0.96219715278770D-2, 0.72189432466630D-2, &
          -0.11651675918591D-2, -0.2152416741149D-3, &
          0.1280502823882D-3, -0.201348547807D-4, &
          -0.12504934821D-5, 0.11330272320D-5, &
          -0.2056338417D-6, 0.61160950D-8, &
          0.50020075D-8, -0.11812746D-8, &
          0.1043427D-9, 0.77823D-11, &
          -0.36968D-11, 0.51D-12, &
          -0.206D-13, -0.54D-14, 0.14D-14, 0.1D-15/
		GR = G(26)
		DO K = 25, 1, -1
			GR = GR*Z + G(K)
		ENDDO
		GAMMA = 1.0D0/(GR*Z)
		IF (DABS(X).GT.1.0D0) THEN
				GAMMA = GAMMA*R
				IF (X.LT.0.0D0) THEN 
					GAMMA = -PI/(X*GAMMA*DSIN(PI*X))
				ENDIF
		ENDIF
	ENDIF

END FUNCTION GAMMA




SUBROUTINE GAUJAC(X,W,N,ALF,BET)
	
!--------------------------------------------------------------------------------------------!
!	Given alf and bet , the parameters α and β of the Jacobi polynomials, this routine returns
!	arrays x(1:n) and w(1:n) containing the abscissas and weights of the n -point Gauss-Jacobi
!	quadrature formula. The largest abscissa is returned in x(1) , the smallest in x(n) .
!--------------------------------------------------------------------------------------------!

	INTEGER, INTENT(IN) :: N
	REAL*8, INTENT(IN) :: ALF, BET
	
	REAL*8, INTENT(OUT) :: X(N), W(N)
	
	INTEGER, PARAMETER :: MAXIT = 20
	
	INTEGER :: I, ITS, J
	REAL*8 :: ALFBET, AN, BN, R1, R2, R3
	REAL*8 :: A, B, C, P1, P2, P3, PP, TEMP, Z, Z1

	DO I = 1, N
		! LOOP OVER THE DESIRED ROOTS.
		IF(I.EQ.1) THEN
			! INITIAL GUESS FOR THE LARGEST ROOT.
			AN = ALF/N
			BN = BET/N
			R1 = (1.0D0+ALF)*(2.780D0/(4.0D0+N*N)+0.7680D0*AN/N)
			R2 = 1.0D0 + 1.480D0*AN + 0.960D0*BN + 0.4520D0*AN*AN + 0.830D0*AN*BN
			Z=1.0D0-R1/R2
		ELSE IF(I.EQ.2)THEN
			! INITIAL GUESS FOR THE SECOND LARGEST ROOT.
			R1 = (4.10D0+ALF)/((1.0D0+ALF)*(1.0D0+0.1560D0*ALF))
			R2 = 1.0D0 + 0.060D0*(N-8.0D0)*(1.+0.120D0*ALF)/N
			R3 = 1.0D0 + 0.0120D0*BET*(1.0D0+0.250D0*ABS(ALF))/N
			Z=Z-(1.0D0-Z)*R1*R2*R3
		ELSE IF(I.EQ.3)THEN
			! INITIAL GUESS FOR THE THIRD LARGEST ROOT.
			R1 = (1.670D0 + 0.280D0*ALF)/(1.0D0 + 0.370D0*ALF)
			R2 = 1.0D0 + 0.220D0*(N-8.0D0)/N
			R3 = 1.0D0 + 8.0D0*BET/((6.280D0+BET)*N*N)
			Z = Z-(X(1)-Z)*R1*R2*R3
		ELSE IF(I.EQ.N-1)THEN
			! INITIAL GUESS FOR THE SECOND SMALLEST ROOT.
			R1 = (1.0D0+0.2350D0*BET)/(0.7660D0 + 0.1190D0*BET)
			R2 = 1.0D0/(1.0D0 + 0.6390D0*(N-4.0D0)/(1.0D0 + 0.710D0*(N-4.0D0)))
			R3 = 1.0D0/(1.0D0 + 20.0D0*ALF/((7.50D0 + ALF)*N*N))
			Z = Z+(Z-X(N-3))*R1*R2*R3
		ELSE IF(I.EQ.N)THEN
			! INITIAL GUESS FOR THE SMALLEST ROOT.
			R1 = (1.0D0 + 0.370D0*BET)/(1.670D0 + 0.280D0*BET)
			R2 = 1.0D0/(1.0D0 + 0.220D0*(N-8.0D0)/N)
			R3 = 1.0D0/(1.0D0 + 8.0D0*ALF/((6.280D0+ALF)*N*N))
			Z = Z+(Z-X(N-2))*R1*R2*R3
		ELSE
			! INITIAL GUESS FOR THE OTHER ROOTS.
			Z = 3.0D0*X(I-1)-3.0D0*X(I-2)+X(I-3)
		ENDIF
		
		ALFBET = ALF + BET
		
		DO ITS = 1, MAXIT
			! REFINEMENT BY NEWTON’S METHOD.
			TEMP = 2.0D0 + ALFBET
			! START THE RECURRENCE WITH P 0 AND P 1 TO AVOID A DIVISION BY ZERO WHEN Α + Β = 0 OR −1.
			P1 = (ALF - BET + TEMP*Z)/2.0D0
			P2 = 1.0D0
			
			DO J = 2, N
				! LOOP UP THE RECURRENCE RELATION TO GET THE JACOBIPOLYNOMIAL EVALUATED AT Z.
				P3=P2
				P2=P1
				TEMP=2.0d0*J+ALFBET
				A=2*J*(J+ALFBET)*(TEMP-2.0D0)
				B=(TEMP-1.0D0)*(ALF*ALF-BET*BET+TEMP*(TEMP-2.0D0)*Z)
				C=2.0D0*(J-1+ALF)*(J-1+BET)*TEMP
				P1=(B*P2-C*P3)/A
			ENDDO
			
			PP=(N*(ALF-BET-TEMP*Z)*P1+2.0D0*(N+ALF)*(N+BET)*P2)/(TEMP*(1.0D0-Z*Z))
			
			!P1 IS NOW THE DESIRED JACOBI POLYNOMIAL. WE NEXT COMPUTE PP, ITS DERIVATIVE, BY A
			!STANDARD RELATION INVOLVING ALSO P2, THE POLYNOMIAL OF ONE LOWER ORDER.
			
			Z1=Z
			Z=Z1-P1/PP	! NEWTON’S FORMULA.
			IF(DABS(Z-Z1).LE.EPS) THEN 
				GOTO 1
			ENDIF
		ENDDO
		
		PRINT*, 'TOO MANY ITERATIONS IN GAUJAC'
		STOP 
		
		1 CONTINUE 
		
		X(I)=Z	! STORE THE ROOT AND THE WEIGHT.
		
		W(I) = DEXP(GAMMLN(ALF+N)+GAMMLN(BET+N)-GAMMLN(N+1.0D0)- GAMMLN(N+ALFBET+1.0D0))*TEMP*2.0D0**ALFBET/(PP*P2)
	ENDDO

END SUBROUTINE GAUJAC

END MODULE GSQUAD
