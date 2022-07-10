C     ****************************************************************** 
C     * ROUTINES TO CONVERT TO MACH NUMBER FROM OTHER NON-DIMENSIONAL  *
C     * COMPRESSIBLE FLOW QUANTITIES. USE F2PY TO COMPILE INTO A       *
C     * PYTHON MODULE.                                                 *
C     *                                                                *
C     * ALL ACCEPT THE PARAMETERS,                                     *
C     * M, DOUBLE 1-D ARRAY : OUTPUT MACH NUMBER;                      *
C     * X, DOUBLE 1-D ARRAY : INPUT QUANTITY;                          *
C     * G, DOUBLE SCALAR : INPUT SPECIFIC HEAT RATIO;                  *
C     * N, INTEGER SCALAR : LENGTH OF ARRAYS (HIDDEN BY F2PY).         *
C     *                                                                *
C     * MCPTO_APO, A_ACRIT ALSO OPTIONALLY ACCEPT,                     *
C     * SUP, LOGICAL : LOOK FOR SUPERSONIC SOLUTIONS (DEFAULT FALSE).  *
C     *                                                                *
C     * WHERE POSSIBLE, USE ANALYTIC INVERSION OF THE FUNCTION,        *
C     * OTHERWISE USE NEWTON/HALLEY ITERATION TO FIND ROOT.            *
C     *                                                                *
C     * JAMES BRIND, DECEMBER 2020                                     *
C     ******************************************************************
C
      SUBROUTINE TO_T(M,X,G,N)
C     ANALYTIC INVERSION SQRT((TO_T-1)*2./(G-1))
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1_2
      GM1_2 = (G-1.0D0)*0.5D0
C     MAIN LOOP
      DO I=1,N
            M(I) = SQRT((X(I)-1.0D0)/GM1_2)
      ENDDO
      END
C
      SUBROUTINE PO_P(M,X,G,N)
C     ANALYTIC INVERSION SQRT((PO_P^((G-1)/G)-1)*2./(G-1))
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1_2
      REAL*8 GM1_G
      GM1_2 = (G-1.0D0)*0.5D0
      GM1_G = (G-1.0D0)/G
C     MAIN LOOP
      DO I=1,N
            M(I) = SQRT((X(I)**GM1_G-1.0D0)/GM1_2)
      ENDDO
      END
C
      SUBROUTINE RHOO_RHO(M,X,G,N)
C     ANALYTIC INVERSION SQRT((RHOO_RHO^(G-1)-1)*2./(G-1))
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1_2
      REAL*8 GM1
      GM1_2 = (G-1.0D0)*0.5D0
      GM1 = G-1.0D0
C     MAIN LOOP
      DO I=1,N
            M(I) = SQRT((X(I)**GM1-1.0D0)/GM1_2)
      ENDDO
      END
C
      SUBROUTINE V_CPTO(M,X,G,N)
C     ANALYTIC INVERSION SQRT( V_cpTo^2/(G-1)/(1-0.5*V_cpTo^2))
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1
      REAL*8 XSQ
      GM1 = G-1.0D0
C     MAIN LOOP
      DO I=1,N
        XSQ = X(I)*X(I)
        M(I) = SQRT(XSQ/GM1/(1.0D0 - XSQ/2.0D0))
      ENDDO
      END
C
      SUBROUTINE F_MCPTO(M,X,G,SUP,N)
C     NUMERIC INVERSION USING NEWTON'S METHOD
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
      LOGICAL SUP
Cf2py optional :: SUP = .FALSE.
Cf2py intent(in) sup
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 SQ_GM1_G
      REAL*8 GM1_2
      REAL*8 C
C     ITERATION VARS
      INTEGER K
      REAL*8 MA
      REAL*8 SQ_TO_T
      REAL*8 MANEW
      REAL*8 ERR
      REAL*8 F
      REAL*8 DF
C     CONSTANTS
      REAL*8 NAN
      REAL*8 TOL
C     CALC GAMMA VARIANTS
      GM1_2 = (G-1.0D0)/2.0D0
      SQ_GM1_G = SQRT(G-1.0D0)/G
C     INITIALISE CONSTANTS
      NAN = 0.0D0
      NAN = 0.0D0/NAN
      TOL = 1.0D-6
C     MAIN LOOP
      DO I=1,N
C         INITIAL GUESS
          IF (SUP) THEN
              MA = 1.5D0
          ELSE
              MA = 0.5D0
          ENDIF
          ERR = HUGE(ERR)
          K = 0
C         UP TO 100 ITERATIONS UNTIL TOLERANCE MET 
          DO WHILE (ERR.gt.TOL.and.K.lt.100)
          K = K + 1
C         USE NEWTON'S METHOD FOR NEW GUESS OF MA
          SQ_TO_T = SQRT(1.0D0 + GM1_2 * MA * MA )
          F = SQ_GM1_G * (1.0D0/MA + G*MA) / SQ_TO_T
          C = SQ_TO_T*SQ_TO_T*SQ_TO_T
          DF = SQ_GM1_G * (MA*MA - 1.0D0) / MA / MA / C 
          MANEW = MA - (F-X(I)) / DF * 0.1D0
C         IF WE CROSS ZERO, THEN SHRINK MA GEOMETRICALLY
          IF (MANEW.lt.0.0D0) THEN
              MANEW = 0.5D0 * MA
          ENDIF
C         GET ERROR AND UPDATE MA
          ERR = ABS(MANEW - MA)
          MA = MANEW
          ENDDO
C         RETURN NAN IF NOT CONVERGED
          IF (ERR.lt.TOL) THEN
              M(I) = MA
          ELSE
              M(I) = NAN
          ENDIF
      ENDDO
      END
C
      SUBROUTINE MCPTO_AP(M,X,G,N)
C     NUMERIC INVERSION USING NEWTON'S METHOD
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1
      REAL*8 GP1
      REAL*8 M_GP1_GM1_2
      REAL*8 G_SQ_GM1
      REAL*8 GM1_2
      REAL*8 GP1_2
C     ITERATION VARS
      INTEGER K
      REAL*8 MA
      REAL*8 TO_T
      REAL*8 MANEW
      REAL*8 ERR
      REAL*8 F
      REAL*8 DF
C     CONSTANTS
      REAL*8 NAN
      REAL*8 TOL
C     CALC GAMMA VARIANTS
      GM1 = G-1.0D0
      GP1 = G+1.0D0
      GM1_2 = GM1/2.0D0
      GP1_2 = GP1/2.0D0
      M_GP1_GM1_2 = GP1/GM1/(-2.0D0)
      G_SQ_GM1 = G/SQRT(GM1)
C     INITIALISE CONSTANTS
      NAN = 0.0D0
      NAN = 0.0D0/NAN
      TOL = 1.0D-6
C     MAIN LOOP
      DO I=1,N
C         INITIAL GUESS
          MA = 0.5D0
C         UP TO 100 ITERATIONS UNTIL TOLERANCE MET 
          ERR = HUGE(ERR)
          K = 0
          DO WHILE (ERR.gt.TOL.and.K.lt.100)
          K = K + 1
C         USE NEWTON'S METHOD FOR NEW GUESS OF MA
          TO_T = (1.0D0 + GM1_2 * MA * MA )
          F = G_SQ_GM1 * MA * SQRT(TO_T)
          DF = G_SQ_GM1 * ( TO_T **2.0D0 - GM1_2 *MA*MA / SQRT(TO_T))
          MANEW = MA - (F-X(I)) / DF
C         GET ERROR AND UPDATE MA
          ERR = ABS(MANEW - MA)
          MA = MANEW
          ENDDO
C         RETURN NAN IF NOT CONVERGED
          IF (ERR.lt.TOL) THEN
              M(I) = MA
          ELSE
              M(I) = NAN
          ENDIF
      ENDDO
      END
C
      SUBROUTINE MCPTO_APO(M,X,G,SUP,N)
C     NUMERIC INVERSION USING NEWTON'S METHOD
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
      LOGICAL SUP
Cf2py optional :: SUP = .FALSE.
Cf2py intent(in) sup
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1
      REAL*8 GP1
      REAL*8 M_GP1_GM1_2
      REAL*8 G_SQ_GM1
      REAL*8 GM1_2
      REAL*8 GP1_2
C     ITERATION VARS
      INTEGER K
      REAL*8 MA
      REAL*8 MANEW
      REAL*8 ERR
      REAL*8 F
      REAL*8 DF
      REAL*8 TO_T
C     CONSTANTS
      REAL*8 FCRIT
      REAL*8 NAN
      REAL*8 TOL
C     CALC GAMMA VARIANTS
      GM1 = G-1.0D0
      GP1 = G+1.0D0
      GM1_2 = GM1/2.0D0
      GP1_2 = GP1/2.0D0
      M_GP1_GM1_2 = GP1/GM1/(-2.0D0)
      G_SQ_GM1 = G/SQRT(GM1)
C     CALC CONSTANTS
      NAN = 0.0D0
      NAN = 0.0D0/NAN
      TOL = 1.0D-6
C     CHOKING VALUE
      FCRIT = G_SQ_GM1 * (1.0D0 + GM1_2)**M_GP1_GM1_2
C     MAIN LOOP
      DO I=1,N
C       RETURN NAN IF CHOKED       
        IF (X(I).gt.FCRIT) THEN
            M(I) = NAN
        ELSE
            IF (SUP) THEN
                MA = 1.5D0
            ELSE
                MA = 0.5D0
            ENDIF
            K = 0
            ERR = HUGE(ERR)
            DO WHILE (ERR.gt.TOL.and.K.lt.100)

            K = K + 1
            TO_T = (1.0D0 + GM1_2 * MA*MA)
            F = G_SQ_GM1 * MA * TO_T**M_GP1_GM1_2
            DF = F *(1.0D0 - GP1_2*MA*MA/TO_T)/MA 
            MANEW = MA - (F-X(I)) / DF
            ERR = ABS(MANEW - MA)
            MA = MANEW
            ENDDO
            IF (ERR.lt.TOL) THEN
                M(I) = MA
            ELSE
                M(I) = NAN
            ENDIF
        ENDIF
      ENDDO
      END
C
      SUBROUTINE A_ACRIT(M,X,G,SUP,N)
C     NUMERIC INVERSION USING NEWTON'S METHOD
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
      LOGICAL SUP
Cf2py optional :: SUP = .FALSE.
Cf2py intent(in) sup
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1
      REAL*8 GP1
      REAL*8 GP1_GM1_2
      REAL*8 GM1_2
      REAL*8 GP1_2
C     ITERATION VARS
      INTEGER K
      REAL*8 MA
      REAL*8 MANEW
      REAL*8 ERR
      REAL*8 F
      REAL*8 DF
      REAL*8 TO_T
      REAL*8 REL_FAC
C     CONSTANTS
      REAL*8 NAN
      REAL*8 TOL
C     CALC GAMMA VARIANTS
      GM1 = G-1.0D0
      GP1 = G+1.0D0
      GM1_2 = GM1/2.0D0
      GP1_2 = GP1/2.0D0
      GP1_GM1_2 = GP1/GM1/(2.0D0)
C     SET CONSTANTS
      NAN = 0.0D0
      NAN = 0.0D0/NAN
      TOL = 1.0D-6
C     MAIN LOOP
      DO I=1,N
        IF (SUP) THEN
            MA = 1.5D0
        ELSE
            MA = 0.5D0
        ENDIF
        K = 0
        ERR = HUGE(ERR)
        DO WHILE (ERR.gt.TOL.and.K.lt.100000)

        K = K + 1

        TO_T = 1.0D0 + GM1_2 * MA*MA 
        F = ((TO_T/GP1_2)**GP1_GM1_2)/MA
        DF = F * (MA - 1.0D0/MA) / TO_T

        ! TO COPE WITH STEEP SLOPES AT LOW MA, REDUCE RELAXATION FACTOR
        IF (MA.lt.0.1D0) THEN
            REL_FAC = 0.01D0
        ELSE
            REL_FAC = 0.1D0
        ENDIF

        ! CATCH NEGATIVE MA AND RESET
        IF (MA.lt.0.0D0) THEN
            MA = 1.0D-3
        ENDIF

        MANEW = MA - (F-X(I)) / DF * REL_FAC
        ERR = ABS(MANEW - MA)
        MA = MANEW
        ENDDO
        IF (ERR.lt.TOL) THEN
            M(I) = MA
        ELSE
            M(I) = NAN
        ENDIF
      ENDDO
      END

      SUBROUTINE POSH_PO(M,X,G,N)
C     NUMERIC INVERSION USING NEWTON'S METHOD
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1
      REAL*8 GP1
      REAL*8 GM1_GP1
      REAL*8 GM1_2
      REAL*8 GP1_2
      REAL*8 G_GM1
      REAL*8 P_GM1
C     ITERATION VARS
      INTEGER K
      REAL*8 A
      REAL*8 B
      REAL*8 C
      REAL*8 MA
      REAL*8 MANEW
      REAL*8 ERR
      REAL*8 F
      REAL*8 DF
      REAL*8 TO_T
      REAL*8 MSQ
C     CONSTANTS
      REAL*8 NAN
      REAL*8 TOL
C     CALC GAMMA VARIANTS
      GM1 = G-1.0D0
      GP1 = G+1.0D0
      GM1_2 = GM1/2.0D0
      GP1_2 = GP1/2.0D0
      GM1_GP1 = GM1/GP1
      G_GM1 = G/GM1
      P_GM1 = 1.0D0/GM1
C     SET CONSTANTS
      NAN = 0.0D0
      NAN = 0.0D0/NAN
      TOL = 1.0D-6
C     MAIN LOOP
      DO I=1,N
        IF (X(I).gt.1.0D0) THEN
            M(I) = NAN
        ELSE
        MA = 1.5D0
        K = 0
        ERR = HUGE(ERR)
        DO WHILE (ERR.gt.TOL.and.K.lt.100)

        K = K + 1

        MSQ = MA*MA
        TO_T = (1.0D0 + GM1_2*MSQ)
        A = GP1_2*MSQ/TO_T
        B = 1.0D0/(G/GP1_2*MSQ - GM1_GP1)
        C = -G*MA*(MSQ-1.0D0)*(MSQ-1.0D0)/TO_T/TO_T
        F = A**G_GM1 * B**P_GM1
        DF = C *(A**P_GM1) *(B**G_GM1)
        MANEW = MA - (F-X(I)) / DF
        ERR = ABS(MANEW - MA)
        MA = MANEW
        ENDDO
        IF (ERR.lt.TOL) THEN
            M(I) = MA
        ELSE
            M(I) = NAN
        ENDIF
        ENDIF
      ENDDO
      END

      SUBROUTINE MASH(M,X,G,N)
C     NUMERIC INVERSION USING NEWTON'S METHOD
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1
      REAL*8 GP1
      REAL*8 GM1_2
      REAL*8 GP1_2
C     ITERATION VARS
      INTEGER K
      REAL*8 MA
      REAL*8 MANEW
      REAL*8 ERR
      REAL*8 F
      REAL*8 DF
      REAL*8 MSQ
      REAL*8 C
      REAL*8 SQRT_TO_T
C     CONSTANTS
      REAL*8 NAN
      REAL*8 TOL
C     CALC GAMMA VARIANTS
      GM1 = G-1.0D0
      GP1 = G+1.0D0
      GM1_2 = GM1/2.0D0
      GP1_2 = GP1/2.0D0
C     SET CONSTANTS
      NAN = 0.0D0
      NAN = 0.0D0/NAN
      TOL = 1.0D-6
C     MAIN LOOP
      DO I=1,N
C       RETURN NAN FOR PRE SHOCK MA LESS THAN UNITY
        IF (X(I).gt.1.0D0) THEN
            M(I) = NAN
        ELSE
        MA = 0.5D0
        K = 0
        ERR = HUGE(ERR)
        DO WHILE (ERR.gt.TOL.and.K.lt.100)

        K = K + 1

        MSQ = MA*MA
        SQRT_TO_T = SQRT(1.0D0 + GM1_2*MSQ)
        F = SQRT_TO_T/SQRT(G*MSQ - GM1_2)
        C = G * (2.0D0*MSQ - 1.0D0) + 1.0D0
        DF = -GP1**2.0D0 * MA / SQRT(2.0D0)/SQRT_TO_T  / C / SQRT(C)

        MANEW = MA - (F-X(I)) / DF

        ERR = ABS(MANEW - MA)
        MA = MANEW


        ENDDO

        IF (ERR.lt.TOL) THEN
            M(I) = MA
        ELSE
            M(I) = NAN
        ENDIF
        ENDIF
      ENDDO
      END
