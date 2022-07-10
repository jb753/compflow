C     ****************************************************************** 
C     * ROUTINES TO CONVERT FROM MACH NUMBER TO OTHER NON-DIMENSIONAL  *
C     * COMPRESSIBLE FLOW QUANTITIES. USE F2PY TO COMPILE INTO A       *
C     * PYTHON MODULE.                                                 *
C     *                                                                *
C     * ALL ACCEPT THE SAME PARAMETERS,                                *
C     * Y, DOUBLE 1-D ARRAY : OUTPUT QUANTITY;                         *
C     * M, DOUBLE 1-D ARRAY : INPUT MACH NUMBER;                       *
C     * G, DOUBLE SCALAR : INPUT SPECIFIC HEAT RATIO;                  *
C     * N, INTEGER SCALAR : LENGTH OF ARRAYS (HIDDEN BY F2PY).         *
C     *                                                                *
C     * JAMES BRIND, DECEMBER 2020                                     *
C     ******************************************************************
C
      SUBROUTINE TO_T(Y,M,G,N)
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1
      GM1 = G-1.0D0
C     MAIN LOOP
      DO I=1,N
            Y(I) = GM1*M(I)
      ENDDO
      END
C
      SUBROUTINE PO_P(Y,M,G,N)
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1_2
      REAL*8 RECIP_GM1
      GM1_2 = (G-1.0D0)*0.5D0
      RECIP_GM1 = 1.0D0/(G-1.0D0)
C     MAIN LOOP
      DO I=1,N
            Y(I) = G * M(I) *(1.0D0 + GM1_2 *M(I)*M(I) )**RECIP_GM1
      ENDDO
      END
C
      SUBROUTINE RHOO_RHO(Y,M,G,N)
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1_2
      REAL*8 M_GM1
      GM1_2 = (G-1.0D0)*0.5D0
      M_GM1 = -1.0D0/(G-1.0D0)
C     MAIN LOOP
      DO I=1,N
            Y(I) = M(I)*(1.0D0 + GM1_2*M(I)*M(I))**M_GM1
      ENDDO
      END
C
      SUBROUTINE V_CPTO(Y,M,G,N)
C     sqrt( (G-1)*M^2 / (1-(G-1)/2*M^2) )
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1_2
      REAL*8 GM1
      REAL*8 SQRT_GM1
      REAL*8 RECIP_SQRT_TO_T
      GM1_2 = (G-1.0D0)*0.5D0
      GM1 = G-1.0D0
      SQRT_GM1 = SQRT(GM1)
C     MAIN LOOP
      DO I=1,N
          RECIP_SQRT_TO_T = 1.0D0/SQRT(1.0D0+GM1_2*M(I)*M(I))
          Y(I)=SQRT_GM1*RECIP_SQRT_TO_T*RECIP_SQRT_TO_T*RECIP_SQRT_TO_T
      ENDDO
      END
C
      SUBROUTINE MCPTO_APO(Y,M,G,N)
C     G/SQRT(G-1)*M*(1+(G-1)/2*M^2)^(-(G+1)/(G-1)/2)
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 M_GP1_GM1_2
      REAL*8 GM1_2
      REAL*8 GP1_2
      REAL*8 G_SQ_GM1
      REAL*8 TO_T
      REAL*8 F
      GM1_2 = (G-1.0D0)/2.0D0
      GP1_2 = (G+1.0D0)/2.0D0
      M_GP1_GM1_2 = (G+1.0D0)/(G-1.0D0)/(-2.0D0)
      G_SQ_GM1 = G / SQRT(G-1.0D0)
C     MAIN LOOP
      DO I=1,N
        TO_T = (1.0D0 + GM1_2 * M(I)*M(I))
        F = G_SQ_GM1 * TO_T**M_GP1_GM1_2
        Y(I) = F *(1.0D0 - GP1_2*M(I)*M(I)/TO_T) 
      ENDDO
      END
C
      SUBROUTINE MCPTO_AP(Y,M,G,N)
C     G/SQRT(G-1)*M*SQRT(1+(G-1)/2*M^2)
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1_2
      REAL*8 G_SQ_GM1
      REAL*8 SQRT_TO_T
      GM1_2 = (G-1.0D0)/2.0D0
      G_SQ_GM1 = G / SQRT(G-1.0D0)
C     MAIN LOOP
      DO I=1,N
        SQRT_TO_T = SQRT(1.0D0 + GM1_2 * M(I) * M(I) )
        Y(I) = G_SQ_GM1* ( SQRT_TO_T + GM1_2 *M(I)*M(I) / SQRT_TO_T)
      ENDDO
      END
C
      SUBROUTINE F_MCPTO(Y,M,G,N)
C     
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1_2
      REAL*8 MSQ
      REAL*8 SQ_GM1_G
      REAL*8 D
      GM1_2 = (G-1.0D0)/2.0D0
      SQ_GM1_G = SQRT(G-1.0D0)/G
C     MAIN LOOP
      DO I=1,N
        MSQ = M(I) * M(I)
        D = SQRT(1.0D0 + GM1_2 * MSQ )
        Y(I) = SQ_GM1_G * (MSQ-1.0D0) / D / D / D  / MSQ
      ENDDO
      END
C
      SUBROUTINE A_ACRIT(Y,M,G,N)
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GP1_GM1_2
      REAL*8 GM1_2
      REAL*8 GP1_2
      REAL*8 TO_T
      REAL*8 MSQ
      GM1_2 = (G-1.0D0)/2.0D0
      GP1_GM1_2 = (G+1.0D0)/(G-1.0D0)/2.0D0
      GP1_2 = (G+1.0D0)/2.0D0
C     MAIN LOOP
      DO I=1,N
        MSQ = M(I)*M(I)
        TO_T = 1.0D0 + GM1_2 * MSQ
        Y(I) = ((TO_T/GP1_2)**GP1_GM1_2)/TO_T*(1.0D0 - 1.0D0/MSQ)
      ENDDO
      END
C
      SUBROUTINE MASH(Y,M,G,N)
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1_2
      REAL*8 GP1
      REAL*8 MSQ
      REAL*8 SQRT_TO_T
      REAL*8 C
      GM1_2 = (G-1.0D0)/2.0D0
      GP1 = G+1.0D0
C     MAIN LOOP
      DO I=1,N
        MSQ = M(I)*M(I)
        SQRT_TO_T = SQRT(1.0D0 + GM1_2*MSQ)
        C = G * (2.0D0*MSQ - 1.0D0) + 1.0D0
        Y(I)= -GP1**2.0D0 * M(I) / SQRT(2.0D0)/SQRT_TO_T  / C / SQRT(C)
      ENDDO
      END
C
      SUBROUTINE POSH_PO(Y,M,G,N)
C     ARGUMENTS
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
C     INTERMEDIATE VARS
      REAL*8 GM1_2
      REAL*8 G_GM1
      REAL*8 M_GM1
      REAL*8 P_GM1
      REAL*8 GP1_2
      REAL*8 G_GP1
      REAL*8 GM1_GP1
      REAL*8 MSQ
      REAL*8 A
      REAL*8 B
      REAL*8 C
      REAL*8 TO_T
      GM1_2 = (G-1.0D0)/2.0D0
      G_GM1 = G/(G-1.0D0)
      G_GP1 = 2.0D0*G/(G+1.0D0)
      M_GM1 = (-1.0D0)/(G-1.0D0)
      p_gm1 = (1.0D0)/(G-1.0D0)
      GP1_2 = (G+1.0D0)/2.0D0
      GM1_GP1 = (G-1.0D0)/(G+1.0D0)
C     MAIN LOOP
      DO I=1,N
        MSQ = M(I)*M(I)
        TO_T = (1.0D0 + GM1_2*MSQ)
        A = GP1_2*MSQ/TO_T
        B = 1.0D0/(G/GP1_2*MSQ - GM1_GP1)
        C = -G*M(I)*(MSQ-1.0D0)*(MSQ-1.0D0)/TO_T/TO_T
        Y(I) = C *(A**P_GM1) *(B**G_GM1)
      ENDDO
      END
