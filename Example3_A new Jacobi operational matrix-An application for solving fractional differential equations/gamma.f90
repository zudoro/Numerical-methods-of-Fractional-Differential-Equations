!gamma(alpha)=dexp(LogGamma(alpha))

RECURSIVE FUNCTION LOGGAMMA(X) &
          RESULT (LOGGAMMA_RESULT)


! LOG OF GAMMA FUNCTION BY ASYMPTOTIC SERIES

        DOUBLE PRECISION X,XTEMP,RECXS,SUM,TERM,X10, &
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
            SUM = SUM - LOG(X10 - K)
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