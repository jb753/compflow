      SUBROUTINE TO_T(M,X,G,N)
c    return np.sqrt((To_T - 1.) * 2. / (ga - 1.))
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
      REAL*8 GM1_2
      GM1_2 = (G-1.0D0)*0.5D0
      DO I=1,N
            M(I) = SQRT((X(I)-1.0D0)/GM1_2)
      ENDDO
      END

      SUBROUTINE PO_P(M,X,G,N)
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
      REAL*8 GM1_2
      REAL*8 GM1_G
      GM1_2 = (G-1.0D0)*0.5D0
      GM1_G = (G-1.0D0)/G
      DO I=1,N
            M(I) = SQRT((X(I)**GM1_G-1.0D0)/GM1_2)
      ENDDO
      END

      SUBROUTINE RHOO_RHO(M,X,G,N)
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
      REAL*8 GM1_2
      REAL*8 GM1
      GM1_2 = (G-1.0D0)*0.5D0
      GM1 = G-1.0D0
      DO I=1,N
            M(I) = SQRT((X(I)**GM1-1.0D0)/GM1_2)
      ENDDO
      END

      SUBROUTINE V_CPTO(M,X,G,N)
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
      REAL*8 GM1
      REAL*8 XSQ
      GM1 = G-1.0D0
C    return np.sqrt( V_cpTo **2. / (ga - 1.) / (1. - 0.5 * V_cpTo **2.))
      DO I=1,N
        XSQ = X(I)*X(I)
        M(I) = SQRT(XSQ/GM1/(1.0D0 - XSQ/2.0D0))
      ENDDO
      END

      SUBROUTINE MCPTO_AP(M,X,G,N)
      INTEGER N
      INTEGER K
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
      REAL*8 GM1
      REAL*8 GP1
      REAL*8 M_GP1_GM1_2
      REAL*8 G_SQ_GM1
      REAL*8 GM1_2
      REAL*8 GP1_2
      REAL*8 FCRIT
      REAL*8 MA
      REAL*8 MANEW
      REAL*8 ERR
      REAL*8 F
      REAL*8 DF


      REAL*8 NAN
      REAL*8 TOL

      GM1 = G-1.0D0
      GM1_2 = (G-1.0D0)/2.0D0
      GP1_2 = (G+1.0D0)/2.0D0
      GP1 = G+1.0D0
      M_GP1_GM1_2 = GP1/GM1/(-2.0D0)
      G_SQ_GM1 = G/SQRT(GM1)

      NAN = 0.0D0
      NAN = 0.0D0/NAN
      TOL = 1.0D-6

      FCRIT = G_SQ_GM1 * (1.0D0 + GM1_2)**M_GP1_GM1_2

      DO I=1,N
        IF (X(I).gt.FCRIT) THEN
            M(I) = NAN
        ELSE
            MA = 0.5D0
            K = 0
            ERR = HUGE(ERR)
            DO WHILE (ERR.gt.TOL.and.K.lt.100)

            K = K + 1

            F = G_SQ_GM1 * MA * SQRT(1.0D0 + GM1_2 * MA*MA )
            DF = G_SQ_GM1 * ( (1.0D0 + GM1_2 * MA * MA )**2.0D0 
     &       - GM1_2 * MA * MA / SQRT(1.0D0 + GM1_2 * MA*MA))

            MANEW = MA - (F-X(I)) / DF

            ERR = ABS(MANEW - MA)
            MA = MANEW

            END DO

            IF (ERR.lt.TOL) THEN
                M(I) = MA
            ELSE
                M(I) = NAN
            ENDIF
        ENDIF
      ENDDO
      END


      SUBROUTINE MCPTO_APO(M,X,G,SUP,N)
      INTEGER N
      INTEGER K
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
      LOGICAL SUP
Cf2py optional,intent(in) :: SUP = .FALSE.
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
      REAL*8 GM1
      REAL*8 GP1
      REAL*8 M_GP1_GM1_2
      REAL*8 G_SQ_GM1
      REAL*8 GM1_2
      REAL*8 GP1_2
      REAL*8 MA
      REAL*8 MANEW
      REAL*8 ERR
      REAL*8 F
      REAL*8 DF


      REAL*8 NAN
      REAL*8 TOL

      GM1 = G-1.0D0
      GM1_2 = (G-1.0D0)/2.0D0
      GP1_2 = (G+1.0D0)/2.0D0
      GP1 = G+1.0D0
      M_GP1_GM1_2 = GP1/GM1/(-2.0D0)
      G_SQ_GM1 = G/SQRT(GM1)

      NAN = 0.0D0
      NAN = 0.0D0/NAN
      TOL = 1.0D-6


      DO I=1,N
        IF (SUP) THEN
            MA = 1.5D0
        ELSE
            MA = 0.5D0
        ENDIF
        K = 0
        ERR = HUGE(ERR)
        DO WHILE (ERR.gt.TOL.and.K.lt.100)

        K = K + 1

        F = G_SQ_GM1 * MA * (1.0D0 + GM1_2 * MA*MA )**M_GP1_GM1_2
        DF = G_SQ_GM1 * ( (1.0D0 + GM1_2 * MA * MA )**M_GP1_GM1_2 
     &       - GP1_2 * MA * MA * (1.0D0 + GM1_2 * MA*MA)**M_GP1_GM1_2 )

        MANEW = MA - (F-X(I)) / DF

        ERR = ABS(MANEW - MA)
        MA = MANEW

        END DO

        IF (ERR.lt.TOL) THEN
            M(I) = MA
        ELSE
            M(I) = NAN
        ENDIF
      ENDDO
      END




      SUBROUTINE A_ACRIT(M,X,G,SUP,N)
      INTEGER N
      INTEGER K
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
      LOGICAL SUP
Cf2py optional,intent(in) :: SUP = .FALSE.
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
      REAL*8 GM1
      REAL*8 GP1
      REAL*8 GP1_GM1_2
      REAL*8 G_SQ_GM1
      REAL*8 GM1_2
      REAL*8 GP1_2
      REAL*8 MA
      REAL*8 MANEW
      REAL*8 ERR
      REAL*8 F
      REAL*8 DF


      REAL*8 NAN
      REAL*8 TOL

      GM1 = G-1.0D0
      GM1_2 = (G-1.0D0)/2.0D0
      GP1_2 = (G+1.0D0)/2.0D0
      GP1 = G+1.0D0
      GP1_GM1_2 = GP1/GM1/(2.0D0)
      G_SQ_GM1 = G/SQRT(GM1)

      NAN = 0.0D0
      NAN = 0.0D0/NAN
      TOL = 1.0D-6


      DO I=1,N
        IF (SUP) THEN
            MA = 1.5D0
        ELSE
            MA = 0.5D0
        ENDIF
        K = 0
        ERR = HUGE(ERR)
        DO WHILE (ERR.gt.TOL.and.K.lt.100)

        K = K + 1

        F = (((1.0D0 + GM1_2 * MA*MA )/GP1_2)**GP1_GM1_2)/MA
        DF = ( (1.0D0 + GM1_2 * MA * MA)/GP1_2 )**GP1_GM1_2 
     &       * -1.0D0/MA/MA + GP1_2 /(1.0D0 + GM1_2 * MA*MA)

        MANEW = MA - (F-X(I)) / DF

        ERR = ABS(MANEW - MA)
        MA = MANEW

        END DO

        IF (ERR.lt.TOL) THEN
            M(I) = MA
        ELSE
            M(I) = NAN
        ENDIF
      ENDDO
      END

      SUBROUTINE POSH(M,X,G,N)
      INTEGER N
      INTEGER K
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
      REAL*8 GM1
      REAL*8 GP1
      REAL*8 GM1_GP1
      REAL*8 G_SQ_GM1
      REAL*8 GM1_2
      REAL*8 GP1_2
      REAL*8 G_GM1
      REAL*8 MG_GM1
      REAL*8 M_GM1
      REAL*8 P_GM1
      REAL*8 MA
      REAL*8 MANEW
      REAL*8 ERR
      REAL*8 F
      REAL*8 DF


      REAL*8 NAN
      REAL*8 TOL

      REAL*8 A
      REAL*8 B
      REAL*8 C


      GM1 = G-1.0D0
      GP1 = G+1.0D0
      GM1_2 = GM1/2.0D0
      GP1_2 = GP1/2.0D0
      GM1_GP1 = GM1/GP1
      G_SQ_GM1 = G/SQRT(GM1)
      G_GM1 = G/GM1
      MG_GM1 = -1.0D0*G/GM1
      M_GM1 = -1.0D0/GM1
      P_GM1 = 1.0D0/GM1

      NAN = 0.0D0
      NAN = 0.0D0/NAN
      TOL = 1.0D-6

      DO I=1,N
        IF (X(I).gt.1.0D0) THEN
            M(I) = NAN
        ELSE
        MA = 1.5D0
        K = 0
        ERR = HUGE(ERR)
        DO WHILE (ERR.gt.TOL.and.K.lt.100)

        K = K + 1

        A = GP1_2*MA*MA/(1.0D0 + GM1_2 * MA * MA)
        B = G/GP1_2*MA*MA - GM1_GP1

        F = A**G_GM1 * B**M_GM1

        A = (G*MA*((MA*MA-1.0D0)/(1.0D0+GM1_2*MA*MA))**2.0D0)
        B = (GP1_2*MA*MA/(1.0D0+GM1_2*MA*MA))
        C = (G/GP1_2*MA*MA-GM1_GP1)

        DF = -1.0D0*(G*MA*((MA*MA-1.0D0)/(1.0D0+GM1_2*MA*MA))**2.0D0) 
     &   *((GP1_2*MA*MA/(1.0D0+GM1_2*MA*MA))**P_GM1)
     &   *((G/GP1_2*MA*MA-GM1_GP1)**MG_GM1)

        MANEW = MA - (F-X(I)) / DF

        ERR = ABS(MANEW - MA)
        MA = MANEW


        END DO

        IF (ERR.lt.TOL) THEN
            M(I) = MA
        ELSE
            M(I) = NAN
        ENDIF
        ENDIF
      ENDDO
      END

      SUBROUTINE MASH(M,X,G,N)
      INTEGER N
      INTEGER K
      REAL*8 G
      REAL*8 M(N)
      REAL*8 X(N)
Cf2py intent(in) x
Cf2py intent(in) g
Cf2py intent(out) m
Cf2py intent(hide) n
      REAL*8 GM1
      REAL*8 GP1
      REAL*8 GP1_GM1_2
      REAL*8 G_SQ_GM1
      REAL*8 GM1_2
      REAL*8 GP1_2
      REAL*8 MA
      REAL*8 MANEW
      REAL*8 ERR
      REAL*8 F
      REAL*8 DF


      REAL*8 NAN
      REAL*8 TOL

      REAL*8 A
      REAL*8 B
      REAL*8 C


      GM1 = G-1.0D0
      GM1_2 = (G-1.0D0)/2.0D0
      GP1_2 = (G+1.0D0)/2.0D0
      GP1 = G+1.0D0
      GP1_GM1_2 = GP1/GM1/(2.0D0)
      G_SQ_GM1 = G/SQRT(GM1)

      NAN = 0.0D0
      NAN = 0.0D0/NAN
      TOL = 1.0D-6

      DO I=1,N
        IF (X(I).gt.1.0D0) THEN
            M(I) = NAN
        ELSE
        MA = 0.5D0
        K = 0
        ERR = HUGE(ERR)
        DO WHILE (ERR.gt.TOL.and.K.lt.100)

        K = K + 1

        F = SQRT((1.0D0 + GM1_2*MA*MA)/(G*MA*MA - GM1_2))

        A = GP1**2.0D0 * MA / SQRT(2.0D0)
        B = -1.0D0/SQRT(1.0D0 + GM1_2 * MA * MA)
        C = G * (2.0D0*MA*MA - 1.0D0) + 1.0D0

        DF = A * B / C / SQRT(C)

        MANEW = MA - (F-X(I)) / DF

        ERR = ABS(MANEW - MA)
        MA = MANEW


        END DO

        IF (ERR.lt.TOL) THEN
            M(I) = MA
        ELSE
            M(I) = NAN
        ENDIF
        ENDIF
      ENDDO
      END
