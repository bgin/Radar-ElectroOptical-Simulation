c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c $Rev: 55 $ $Date: 2014-12-31 12:16:59 -0500 (Wed, 31 Dec 2014) $
c FORTRAN 77
c 
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      REAL(kind=sp) FUNCTION BDREF(MU, MUP, DPHI,
     &                     BRDF_TYPE, BRDF_ARG)
       !dir$ optimize:3
       !dir$ attributes forceinline :: BDREF
       !dir$ attributes code_align : 32 :: BDREF
       use mod_kinds, only : i4, sp
c     Supplies surface bi-directional reflectivity.
c
c     NOTE 1: Bidirectional reflectivity in DISORT is defined
c             by Eq. 39 in STWL.
c     NOTE 2: Both MU and MU0 (cosines of reflection and incidence
c             angles) are positive.
c
c  INPUT:
c
c    MU     : Cosine of angle of reflection (positive)
c
c    MUP    : Cosine of angle of incidence (positive)
c
c    DPHI   : Difference of azimuth angles of incidence and reflection
c                (radians)
c
c  LOCAL VARIABLES:
c
c    IREF   : bidirectional reflectance options
c             1 - Hapke's BDR model
c             2 - Cox-Munk BDR model
c             3 - RPV BDR model
c             4 - Ross-Li BDR model
c
c    B0     : empirical factor to account for the finite size of
c             particles in Hapke's BDR model
c
c    B      : term that accounts for the opposition effect
c             (retroreflectance, hot spot) in Hapke's BDR model
c
c    CTHETA : cosine of phase angle in Hapke's BDR model
c
c    GAMMA  : albedo factor in Hapke's BDR model
c
c    H0     : H( mu0 ) in Hapke's BDR model
c
c    H      : H( mu ) in Hapke's BDR model
c
c    HH     : angular width parameter of opposition effect in Hapke's
c             BDR model
c
c    P      : scattering phase function in Hapke's BDR model
c
c    THETA  : phase angle (radians); the angle between incidence and
c             reflection directions in Hapke's BDR model
c
c    W      : single scattering albedo in Hapke's BDR model
c
c
c   Called by- DREF, SURFAC
c +-------------------------------------------------------------------+
c     .. Scalar Arguments ..
      REAL(kind=sp)      DPHI, MU, MUP, BRDF_ARG(4)
      INTEGER(kind=i4)   BRDF_TYPE
c     ..
c     .. Local Scalars ..
      INTEGER(kind=i4)   IREF
      REAL(kind=sp)      B0, H0, HH, W
      REAL(kind=sp)      PWS, REFRAC_INDEX, BDREF_F
      REAL(kind=sp)      PI
      REAL(kind=sp)      RHO0, KAPPA, G  
      REAL(kind=sp)      K_ISO, K_VOL, K_GEO, ALPHA0 
      LOGICAL   DO_SHADOW
c     ..
c     .. External Subroutines ..
      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC COS, SQRT
c     ..

      PI   = 2.*ASIN(1.)

      IREF = BRDF_TYPE

c     ** 1. Hapke BRDF
      IF ( IREF.EQ.1 ) THEN

c       ** Hapke's BRDF model (times Pi/Mu0) (Hapke, B., Theory of reflectance
c       ** and emittance spectroscopy, Cambridge University Press, 1993, Eq.
c       ** 8.89 on page 233. Parameters are from Fig. 8.15 on page 231, expect
c       ** for w.)

        B0 = BRDF_ARG(1) !1.0
        HH = BRDF_ARG(2) !0.06
        W  = BRDF_ARG(3) !0.6

        CALL BRDF_HAPKE(MUP, MU, DPHI,
     &                  B0, HH, W, PI,
     &                  BDREF)

c     ** 2. Cox-Munk BRDF
      ELSEIF(IREF.EQ.2) THEN

c        PRINT *, "Calling oceabrdf"

        PWS          =  BRDF_ARG(1)
        REFRAC_INDEX =  BRDF_ARG(2)

        IF(BRDF_ARG(3) .EQ. 1) THEN
          DO_SHADOW = .TRUE.
        ELSEIF(BRDF_ARG(3) .EQ. 0) THEN
          DO_SHADOW = .FALSE.
        ELSE
          PRINT *, "ERROR SHADOW ARGUMENTS"
        ENDIF

        CALL OCEABRDF2(DO_SHADOW,
     &                 REFRAC_INDEX, PWS, 
     &                 MUP, MU, DPHI,
     &                 BDREF_F)

        BDREF = BDREF_F

c     ** 3. RPV BRDF
      ELSEIF(IREF .EQ. 3) THEN

        RHO0  =  BRDF_ARG(1) !0.027
        KAPPA =  BRDF_ARG(2) !0.647
        G     =  BRDF_ARG(3) !-0.169   !asymmetry factor for HG
        H0    =  BRDF_ARG(4) !0.100

        CALL BRDF_RPV(MUP, MU, DPHI,
     &                RHO0, KAPPA, G, H0,
     &                BDREF_F)

        BDREF = BDREF_F

c     ** 4. Ross-Li BRDF
      ELSEIF(IREF .EQ. 4) THEN
        
        K_ISO  = BRDF_ARG(1)   !0.200
        K_VOL  = BRDF_ARG(2)   !0.020
        K_GEO  = BRDF_ARG(3)   !0.300
        ALPHA0 = 1.5*pi/180.

        CALL BRDF_ROSSLI(MUP, MU, DPHI,
     &                   K_ISO, K_VOL, K_GEO,
     &                   ALPHA0,
     &                   BDREF_F)

        BDREF = BDREF_F

        IF(BDREF .LT. 0.00) THEN
          BDREF = 0.00
        ENDIF

      ELSE

        CALL ERRMSG( 'BDREF--Need to supply surface BDRF model',
     &                 .TRUE.)

      ENDIF

      RETURN
      END FUNCTION
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c +--------------------------------------------------------------------
      SUBROUTINE BRDF_HAPKE ( MUP, MU, DPHI,
     &                        B0, HH, W, PI,
     &                        BRDF )
       !dir$ optimize:3
       !dir$ attributes forceinline :: BRDF_HAPKE
       !dir$ attributes code_align : 32 :: BRDF_HAPKE
       use mod_kinds, only : i4, sp
c +--------------------------------------------------------------------
c Hapke "Theory of Reflectance and Emittance Spectroscopy" Chapter 10, Page 262
c Eq. (10.2).
c Version 3 fix: definition of phase angle / scattering angle see DISORT3
c paper Eqs. (25-26).
c +--------------------------------------------------------------------
      IMPLICIT NONE
      REAL(kind=sp) MUP, MU, DPHI
      REAL(kind=sp) B0, HH, W, PI
      REAL(kind=sp) BRDF
      REAL(kind=sp) CALPHA, ALPHA, P, B, H0, GAMMA, H

      CALPHA = MU * MUP - (1.-MU**2)**.5 * (1.-MUP**2)**.5
     &         * COS( DPHI )

      ALPHA = ACOS( CALPHA )

      P     = 1. + 0.5 * CALPHA

      B     = B0 * HH / ( HH + TAN( ALPHA/2.) )

      GAMMA = SQRT( 1. - W )
      H0   = ( 1. + 2.*MUP ) / ( 1. + 2.*MUP * GAMMA )
      H    = ( 1. + 2.*MU ) / ( 1. + 2.*MU * GAMMA )

c     ** Version 3: add factor PI
      BRDF = W / (4.*PI) / (MU+MUP) * ( (1.+B)* P + H0 * H - 1.0 )
c     BRDF = W / 4. / (MU+MUP) * ( (1.+B)* P + H0 * H - 1.0 )

      END
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c +--------------------------------------------------------------------
      SUBROUTINE BRDF_RPV(MU_I, MU_R, DPHI,
     &                    RHO0, KAPPA, G_HG, H0,
     &                    BRDF)
       !dir$ optimize:3
       !dir$ attributes forceinline :: BRDF_RPV
       !dir$ attributes code_align : 32 :: BRDF_RPV
       use mod_kinds, only : i4, sp
c +--------------------------------------------------------------------
c DISORT Version 3: RPV BRDF
c   Input:
c
c   MU_I:  absolute cosine of incident polar angle (positive)
c   MU_R:  absolute cosine of reflected polar angle (positive)
c   DPHI:  relative azimuth to incident vector; (pi - dphi), sun-view relative
c          azimuth sun located at phi = 180, while incident solar beam located
c          at phi = 0
c   RHO0:  RPV BRDF parameter, control reflectance
c   KAPPA: PRV BRDF parameter, control anisotropy
c   G:     RPV BRDF parameter, H-G asymmetry factor
c   H0:    RPV BRDF parameter, control hot spot (back scattering direction)
c
c   Output:
c
c   BRDF:  RPV BRDF
c +--------------------------------------------------------------------
      IMPLICIT NONE
      REAL(kind=sp) :: MU_I, MU_R, DPHI
      REAL(kind=sp) :: RHO0, KAPPA, G_HG, H0
      REAL(kind=sp) :: BRDF
      REAL(kind=sp) :: PI
      REAL(kind=sp) :: COS_ALPHA
      REAL(kind=sp) :: SIN_I, SIN_R, TAN_I, TAN_R
      REAL(kind=sp) :: G_SQ, G, F

      PI    = 2.*ASIN(1.)

      SIN_I = SQRT(1. - MU_I*MU_I)
      SIN_R = SQRT(1. - MU_R*MU_R)
      TAN_I = SIN_I/MU_I
      TAN_R = SIN_R/MU_R

      COS_ALPHA = MU_I*MU_R - SIN_I*SIN_R
     & *COS(DPHI)

      G_SQ = TAN_I*TAN_I + TAN_R*TAN_R 
     &    + 2.*TAN_I*TAN_R*COS(DPHI)

c     ** hot spot
      G = SQRT(G_SQ)

c     ** HG phase function
      F = (1. - G_HG*G_HG)/
     &     (1+G_HG*G_HG+2.*G_HG*COS_ALPHA)**1.5


c     ** BRDF semiempirical function
      BRDF = RHO0 
     &      * (MU_I*MU_R*(MU_I+MU_R))**(KAPPA-1.)
     &      * F
     &      * (1. + ((1.-H0)/(1.+G)))

      END
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c +--------------------------------------------------------------------
      SUBROUTINE BRDF_ROSSLI(MU_I, MU_R, DPHI,
     &                       K_ISO, K_VOL, K_GEO,
     &                       ALPHA0,
     &                       BRDF)
       !dir$ optimize:3
       !dir$ attributes forceinline :: BRDF_ROSSLI	
       !dir$ attributes code_align : 32 :: BRDF_ROSSLI
       use mod_kinds, only : i4, sp
c +--------------------------------------------------------------------
c Version 3: Ross-Li BRDF
c   Input:
c
c   MU_I:    absolute cosine of incident polar angle (positive)
c   MU_R:    absolute cosine of reflected polar angle (positive)
c   DPHI:  relative azimuth to incident vector; (pi - dphi), sun-view relative
c          azimuth sun located at phi = 180, while incident solar beam located
c          at phi = 0
c   K_ISO:   BRDF parameter, isotropic scattering kernel
c   K_VOL:   BRDF parameter, volume scattering kernel
c   K_GEO:   BRDF parameter, geometry scattering kernel
c   ALPHA0:  BRDF parameter, control hot spot (back scattering direction)
c
c   Output:
c   BRDF:  Ross-Li BRDF
c
c +--------------------------------------------------------------------
      IMPLICIT NONE
      REAL(kind=sp) :: MU_I, MU_R, DPHI
      REAL(kind=sp) :: F_GEO, F_VOL
      REAL(kind=sp) :: K_ISO, K_GEO, K_VOL
      REAL(kind=sp) :: RATIO_HB, RATIO_BR
      REAL(kind=sp) :: BRDF
      REAL(kind=sp) :: PI
      REAL(kind=sp) :: COS_ALPHA, SIN_ALPHA
      REAL(kind=sp) :: COS_ALPHA1
      REAL(kind=sp) :: ALPHA
      REAL(kind=sp) :: SIN_I, SIN_R, TAN_I, TAN_R
      REAL(kind=sp) :: SIN_I1, SIN_R1, COS_I1, COS_R1, TAN_I1, TAN_R1
      REAL(kind=sp) :: G_SQ, COS_T, T       
      REAL(kind=sp) :: C, ALPHA0
c +--------------------------------------------------------------------

c      PRINT *, MU_I, MU_R, DPHI,
c     &        K_ISO, K_GEO, K_VOL,
c     &        THETA0
c      PRINT *,

      RATIO_HB = 2.
      RATIO_BR = 1.
      PI       = 2.*ASIN(1.)

      SIN_I = SQRT(1. - MU_I*MU_I)
      SIN_R = SQRT(1. - MU_R*MU_R)
      TAN_I = SIN_I/MU_I
      TAN_R = SIN_R/MU_R

      COS_ALPHA = MU_I*MU_R - SIN_I*SIN_R
     & *COS(DPHI)
      SIN_ALPHA = SQRT(1. - COS_ALPHA*COS_ALPHA)
      ALPHA = ACOS(COS_ALPHA)

c     ** Compute KERNEL RossThick
      C     = 1. + 1./(1.+ALPHA/ALPHA0)
      F_VOL = 4./(3.*PI) * (1./(MU_I+MU_R))
     &       * ((PI/2. - ALPHA)*COS_ALPHA+SIN_ALPHA)*C - 1./3.

c      K1 = ((PI/2. - ALPHA)*COS_ALPHA + SIN_ALPHA)
c     &       /(MU_I + MU_R) - PI/4.


c     ** Compute KERNEL LSR
      TAN_I1 = RATIO_BR * TAN_I
      TAN_R1 = RATIO_BR * TAN_R
      SIN_I1 = TAN_I1/SQRT(1.+ TAN_I1*TAN_I1)
      SIN_R1 = TAN_R1/SQRT(1.+ TAN_R1*TAN_R1)
      COS_I1 = 1./SQRT(1.+ TAN_I1*TAN_I1)
      COS_R1 = 1./SQRT(1.+ TAN_R1*TAN_R1)

      COS_ALPHA1 = COS_I1*COS_R1 - SIN_I1*SIN_R1
     &            *COS(DPHI)

      G_SQ = TAN_I1*TAN_I1 + TAN_R1*TAN_R1 
     &      + 2.*TAN_I1*TAN_R1*COS(DPHI)

c      M = 1./COS_I1 + 1./COS_R1

      COS_T = RATIO_HB *(COS_I1*COS_R1)/(COS_I1+COS_R1)
     &       *SQRT(G_SQ + (TAN_I1*TAN_R1*SIN(DPHI))**2)
  
      IF(COS_T .LE. 1. .AND. COS_T .GE. -1.) THEN
        T = ACOS(COS_T)
      ELSE
        T = 0.
      ENDIF

      F_GEO = (COS_I1+COS_R1)/(PI*COS_I1*COS_R1)*(T-SIN(T)*COS(T)-PI)   
     &       + (1.+ COS_ALPHA1)/(2.*COS_I1*COS_R1)

c     Compute BRDF

c      PRINT *, RATIO_HB, D_SQ, 
c     &    TAN_I1*TAN_R1*SIN(DPHI),
c     &    M, COS_T

c      BRDF = K1
      BRDF = K_ISO + K_GEO*F_GEO + K_VOL*F_VOL

      END
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c +--------------------------------------------------------------------
      SUBROUTINE OCEABRDF2
     &       ( DO_SHADOW, 
     &         REFRAC_INDEX, WS,
     &         MU_I, MU_R, DPHI,
     &         BRDF)
       !dir$ optimize:3
       !dir$ attributes forceinline :: OCEABRDF2
       !dir$ attributes code_align : 32 :: OCEABRDF2
       use mod_kinds, only : i4, sp
c +--------------------------------------------------------------------
c Version 3: 1D Gaussian Rough Ocean BRDF
c   Input:
c
c   mu_i:         absolute cosine of incident polar angle (positive)
c   mu_r:         absolute cosine of reflected polar angle (positive)
c   dphi:         relative azimuth (radians) 
c   do_shadow:    BRDF parameter, open/close shadow effect 
c   refrac_index: BRDF parameter, refractive index of boundary media (water)
c   ws:           BRDF parameter, wind speed (m/s)
c
c   Output:
c
c   brdf:         1D Gaussian Rough Ocean BRDF
c          
c +--------------------------------------------------------------------
      LOGICAL :: DO_SHADOW
      REAL(kind=sp) ::    REFRAC_INDEX, WS
      REAL(kind=sp) ::    SIN_I, SIN_R, MU_I, MU_R, DPHI, BRDF
      REAL(kind=sp) ::    COS_THETA, SIGMA_SQ, MU_N_SQ, P
      REAL(kind=sp) ::    N_I, N_T, COS_LI, COS_LT, SIN_LI, SIN_LT
      REAL(kind=sp) ::    R_S, R_P, R
      REAL(kind=sp) ::    SHADOW
      REAL(kind=sp) ::    PI

      PI = 2.*ASIN(1.)

c     ** Cox Munk slope distribution
      SIN_I = SQRT(1. - MU_I*MU_I)
      SIN_R = SQRT(1. - MU_R*MU_R)

      COS_THETA = -MU_I*MU_R + SIN_I*SIN_R*COS(DPHI)
      MU_N_SQ   = (MU_I + MU_R)*(MU_I + MU_R)/(2.*(1.-COS_THETA))   

      SIGMA_SQ  = 0.003 + 0.00512*WS

      P = 1./(PI*SIGMA_SQ) * EXP( -(1-MU_N_SQ)/(SIGMA_SQ*MU_N_SQ) )

c     ** Fresnel reflectance

      N_I = 1.0
      N_T = REFRAC_INDEX

      SIN_LI = SQRT( 1.-0.5*(1.-COS_THETA) ) 
      COS_LI = SQRT( 0.5*(1.-COS_THETA) ) 
      SIN_LT = N_I*SIN_LI/N_T
      COS_LT = SQRT(1. - SIN_LT*SIN_LT)

      R_S = (N_I*COS_LI-N_T*COS_LT)/(N_I*COS_LI+N_T*COS_LT)
      R_P = (N_T*COS_LI-N_I*COS_LT)/(N_I*COS_LT+N_T*COS_LI)

      R = 0.5*(R_S*R_S + R_P*R_P)

c     ** Rough surface BRDF
      BRDF = (P*R)/(4.*MU_I*MU_R*MU_N_SQ*MU_N_SQ)

c     Shadowing effect (see Tsang, Kong, Shin, Theory of Microwave Remote
c     Sensing, Wiley-Interscience, 1985) 
      IF(DO_SHADOW) THEN
        SHADOW = 1./( SHADOW_ETA(MU_I, SIGMA_SQ, PI) 
     &          + SHADOW_ETA(MU_R, SIGMA_SQ, PI) + 1. )
        BRDF = BRDF*SHADOW
      ENDIF

      END
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c +--------------------------------------------------------------------
      REAL(kind=sp) FUNCTION SHADOW_ETA(COS_THETA, SIGMA_SQ, PI)
       !dir$ optimize:3
       !dir$ attributes forceinline :: SHADOW_ETA
       !dir$ attributes code_align : 32 :: SHADOW_ETA
c +--------------------------------------------------------------------
c Version 3: shadow effect function
c            called by OCEABRDF2
c   Input:
c
c   COS_THETA     absolute cosine of incident/reflected polar angle (positive)
c   SIGMA_SQ      slope variance 
c   PI            3.141592653... constant
c
c   Output:
c
c   SHADOW_ETA:   shadow function
c +--------------------------------------------------------------------
      REAL(kind=sp) :: COS_THETA, SIN_THETA
      REAL(kind=sp) :: MU, SIGMA_SQ, PI
      REAL(kind=sp) :: TERM1, TERM2

      SIN_THETA = SQRT(1.-COS_THETA*COS_THETA)
      MU = COS_THETA/SIN_THETA

      TERM1 = SQRT(SIGMA_SQ/PI)/MU*EXP( -MU*MU/(SIGMA_SQ) )
      TERM2 = ERFC( MU/SQRT(SIGMA_SQ) )

      SHADOW_ETA = 0.5*(TERM1 - TERM2)

      END
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
