      SUBROUTINE TO_T(Y,M,G,N)
C	  1-(g-1)/2 * M^2
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
      REAL*8 GM1_2
      GM1_2 = (G-1.0D0)*0.5D0
      DO I=1,N
            Y(I) = 1.0D0 + GM1_2 * M(I) *M(I)
      ENDDO
      END

      SUBROUTINE PO_P(Y,M,G,N)
C	  (1-(G-1)/2 * M^2)^(G/(G-1))
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
      REAL*8 GM1_2
      REAL*8 G_GM1
      GM1_2 = (G-1.0D0)*0.5D0
      G_GM1 = G/(G-1.0D0)
      DO I=1,N
            Y(I) = (1.0D0 + GM1_2 * M(I) *M(I))**G_GM1
      ENDDO
      END

      SUBROUTINE RHOO_RHO(Y,M,G,N)
C	  (1-(G-1)/2 * M^2)^(1/(G-1))
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
      REAL*8 GM1_2
      REAL*8 RECIP_GM1
      GM1_2 = (G-1.0D0)*0.5D0
      RECIP_GM1 = 1/(G-1.0D0)
      DO I=1,N
            Y(I) = (1.0D0 + GM1_2 * M(I) *M(I))**RECIP_GM1
      ENDDO
      END

      SUBROUTINE V_CPTO(Y,M,G,N)
C	  (1-(G-1)/2 * M^2)^(1/(G-1))
C    V_cpTo = np.sqrt(2.) * np.ones_like(Ma)  # Limit for Ma -> inf
C    V_cpTo[ii] = np.sqrt( (ga - 1.0) * Ma[ii] ** 2. / (1. + (ga - 1.) * Ma[ii] **2. /2.) )
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
      REAL*8 GM1_2
      REAL*8 GM1
      GM1_2 = (G-1.0D0)*0.5D0
      GM1 = G-1.0D0
      DO I=1,N
            IF (M(I).lt.1.0D-2) THEN
                Y(I)=0.0D0
            ELSE
                Y(I)=SQRT(GM1/(1.0D0/M(I)/M(I)+GM1_2))
            ENDIF
      ENDDO
      END


      SUBROUTINE MCPTO_APO(Y,M,G,N)
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
      REAL*8 M_GP1_GM1_2
      REAL*8 GM1_2
      REAL*8 G_SQ_GM1
      GM1_2 = (G-1.0D0)/2.0D0
      M_GP1_GM1_2 = (G+1.0D0)/(G-1.0D0)/(-2.0D0)
      G_SQ_GM1 = G / SQRT(G-1.0D0)
      DO I=1,N
        Y(I)=G_SQ_GM1*M(I)*(1.0D0+GM1_2*M(I)*M(I))**M_GP1_GM1_2
      ENDDO
      END


      SUBROUTINE MCPTO_AP(Y,M,G,N)
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
      REAL*8 GM1_2
      REAL*8 G_SQ_GM1
      GM1_2 = (G-1.0D0)/2.0D0
      G_SQ_GM1 = G / SQRT(G-1.0D0)
      DO I=1,N
        Y(I)=G_SQ_GM1*M(I)*SQRT(1.0D0+GM1_2*M(I)*M(I))
      ENDDO
      END


      SUBROUTINE A_ACRIT(Y,M,G,N)
C    A_Acrit[ii] = ( 1./Ma[ii] * (2. / (ga + 1.0) * To_T[ii])
C                    ** (0.5 * (ga + 1.0) / (ga - 1.0)))
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
      REAL*8 GP1_GM1_2
      REAL*8 GM1_2
      REAL*8 GP1_2
      GM1_2 = (G-1.0D0)/2.0D0
      GP1_GM1_2 = (G+1.0D0)/(G-1.0D0)/2.0D0
      GP1_2 = (G+1.0D0)/2.0D0
      DO I=1,N
        Y(I)=(((1.0D0+GM1_2*M(I)*M(I))/GP1_2)**GP1_GM1_2)/M(I)
      ENDDO
      END


      SUBROUTINE MASH(Y,M,G,N)
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
      REAL*8 GM1_2
      REAL*8 MSQ
      GM1_2 = (G-1.0D0)/2.0D0
      DO I=1,N
        MSQ = M(I)*M(I)
        IF (MSQ.lt.1.0D-3) THEN
            Y(I)=SQRT((1.0D0+GM1_2*MSQ)/(G*MSQ-GM1_2))
        ELSE
            Y(I)=SQRT((1.0D0/MSQ+GM1_2)/(G-GM1_2/MSQ))
        ENDIF
      ENDDO
      END

      SUBROUTINE POSH(Y,M,G,N)
      INTEGER N
      REAL*8 G
      REAL*8 M(N)
      REAL*8 Y(N)
Cf2py intent(in) m
Cf2py intent(in) g
Cf2py intent(out) y
Cf2py intent(hide) n
      REAL*8 GM1_2
      REAL*8 G_GM1
      REAL*8 M_GM1
      REAL*8 MSQ
      REAL*8 GP1_2
      REAL*8 G_GP1
      REAL*8 GM1_GP1
      REAL*8 A
      REAL*8 B
      GM1_2 = (G-1.0D0)/2.0D0
      G_GM1 = G/(G-1.0D0)
      G_GP1 = 2.0D0*G/(G+1.0D0)
      M_GM1 = (-1.0D0)/(G-1.0D0)
      GP1_2 = (G+1.0D0)/2.0D0
      GM1_GP1 = (G-1.0D0)/(G+1.0D0)
!     A = 0.5 * (ga + 1.) * Ma[iiv] ** 2. / To_T[iiv]
!     B = 2. * ga / (ga + 1.0) * Ma[iiv] ** 2. - 1. / (ga + 1.0) * (ga - 1.0)
!     Posh_Po[iiv] = (A ** (ga / (ga - 1.0)) * B ** (-1. / (ga - 1.)))
      DO I=1,N
        MSQ = M(I)*M(I)
        A = GP1_2/(1.0D0/MSQ+GM1_2)
        B = G_GP1 * MSQ - GM1_GP1
        Y(I)=A ** G_GM1 * B ** M_GM1
      ENDDO
      END



