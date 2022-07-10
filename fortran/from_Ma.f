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
C     1-(G-1)/2*M^2
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
      GM1_2 = (G-1.0D0)*0.5D0
C     MAIN LOOP
      DO I=1,N
            Y(I) = 1.0D0 + GM1_2 * M(I) *M(I)
      ENDDO
      END
C
      SUBROUTINE PO_P(Y,M,G,N)
C     (1-(G-1)/2*M^2)^(G/(G-1))
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
      GM1_2 = (G-1.0D0)*0.5D0
      G_GM1 = G/(G-1.0D0)
C     MAIN LOOP
      DO I=1,N
            Y(I) = (1.0D0 + GM1_2 *M(I)*M(I) )**G_GM1
      ENDDO
      END
C
      SUBROUTINE RHOO_RHO(Y,M,G,N)
C     (1-(G-1)/2*M^2)^(1/(G-1))
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
            Y(I) = (1.0D0 + GM1_2*M(I)*M(I))**RECIP_GM1
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
      GM1_2 = (G-1.0D0)*0.5D0
      GM1 = G-1.0D0
C     MAIN LOOP
      DO I=1,N
C           NUMERIC FORM FOR LOW MA
            IF (M(I).lt.1.0D-2) THEN
                Y(I)=SQRT(GM1*M(I)*M(I)/(1.0D0+GM1_2*M(I)*M(I)))
            ELSE
C           NUMERIC FORM FOR HIGH MA
                Y(I)=SQRT(GM1/(1.0D0/M(I)/M(I)+GM1_2))
            ENDIF
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
      REAL*8 G_SQ_GM1
      GM1_2 = (G-1.0D0)/2.0D0
      M_GP1_GM1_2 = (G+1.0D0)/(G-1.0D0)/(-2.0D0)
      G_SQ_GM1 = G / SQRT(G-1.0D0)
C     MAIN LOOP
      DO I=1,N
        Y(I)=G_SQ_GM1*M(I)*(1.0D0+GM1_2*M(I)*M(I))**M_GP1_GM1_2
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
      GM1_2 = (G-1.0D0)/2.0D0
      G_SQ_GM1 = G / SQRT(G-1.0D0)
C     MAIN LOOP
      DO I=1,N
        Y(I)=G_SQ_GM1*M(I)*SQRT(1.0D0+GM1_2*M(I)*M(I))
      ENDDO
      END
C
      SUBROUTINE F_MCPTO(Y,M,G,N)
C     SQRT(G-1)/G/M*(1+G*M^2)/SQRT(1+(G-1)/2*M^2)
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
      REAL*8 SQ_GM1_G
      REAL*8 C
      REAL*8 INF
      GM1_2 = (G-1.0D0)/2.0D0
      SQ_GM1_G = SQRT(G-1.0D0)/G
C     MAIN LOOP
      DO I=1,N
        IF (M(I).eq.0.0D0) THEN
            INF = 0.0D0
            INF = 1.0D0/INF
            Y(I) = INF
        ELSE
            C = SQRT(1.0D0/M(I)/M(I)+GM1_2)
            Y(I)=SQ_GM1_G*(1.0D0/M(I)/M(I)+G)/C
        ENDIF
      ENDDO
      END
C
      SUBROUTINE A_ACRIT(Y,M,G,N)
C     (2/(G+1)*(1+(G-1)/2*M^2))^((G+1)/(G-1)/2)/M
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
      GM1_2 = (G-1.0D0)/2.0D0
      GP1_GM1_2 = (G+1.0D0)/(G-1.0D0)/2.0D0
      GP1_2 = (G+1.0D0)/2.0D0
C     MAIN LOOP
      DO I=1,N
        Y(I)=(((1.0D0+GM1_2*M(I)*M(I))/GP1_2)**GP1_GM1_2)/M(I)
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
      REAL*8 MSQ
      GM1_2 = (G-1.0D0)/2.0D0
C     MAIN LOOP
      DO I=1,N
        MSQ = M(I)*M(I)
C       NUMERIC FORM FOR LOW MA
        IF (MSQ.lt.1.0D-2) THEN
            Y(I)=SQRT((1.0D0+GM1_2*MSQ)/(G*MSQ-GM1_2))
C       NUMERIC FORM FOR HIGH MA
        ELSE
            Y(I)=SQRT((1.0D0/MSQ+GM1_2)/(G-GM1_2/MSQ))
        ENDIF
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
      REAL*8 GP1_2
      REAL*8 G_GP1
      REAL*8 GM1_GP1
      REAL*8 MSQ
      GM1_2 = (G-1.0D0)/2.0D0
      G_GM1 = G/(G-1.0D0)
      G_GP1 = 2.0D0*G/(G+1.0D0)
      M_GM1 = (-1.0D0)/(G-1.0D0)
      GP1_2 = (G+1.0D0)/2.0D0
      GM1_GP1 = (G-1.0D0)/(G+1.0D0)
C     MAIN LOOP
      DO I=1,N
        MSQ = M(I)*M(I)
        Y(I)=(GP1_2/(1.0D0/MSQ+GM1_2))**G_GM1*(G_GP1*MSQ-GM1_GP1)**M_GM1
      ENDDO
      END
