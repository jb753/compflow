      SUBROUTINE TO_T(Y,M,G,N)
C	  1-(G-1)/2 * M^2
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

      SUBROUTINE V_CTPO(Y,M,G,N)
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
            Y(I)=SQRT(GM1*M(I)*M(I)/(1.0D0+GM1_2*M(I)*M(I)))
      ENDDO
      END



