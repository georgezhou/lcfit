!=======================================================================
!     PROGRAM JKTEBOP             John Taylor   j.k.taylor@warwick.ac.uk
!                                 Astrophysics Group  Warwick University
!-----------------------------------------------------------------------
! This is the library version of jktebop
!-----------------------------------------------------------------------
! V(1) = surface brightness ratio    V(15) = third light
! V(2) = sum of fractional radii     V(16) = phase correction
! V(3) = ratio of stellar radii      V(17) = light scaling factor
! V(4) = linear LD for star A        V(18) = integration ring size (deg)
! V(5) = linear LD for star B        V(19) = orbital period (days)
! V(6) = orbital inclination         V(20) = ephemeris timebase (days)
! V(7) = e cos(omega) OR ecentricity V(21) = nonlinear LD for star A
! V(8) = e sin(omega) OR omega       V(22) = nonlinear LD for star B
! V(9) = gravity darkening 1         VEXTRA(1) = primary star radius
! V(10)= gravity darkening 2         VEXTRA(2) = secondary star radius
! V(11) = primary reflected light    VEXTRA(3) = stellar light ratio
! V(12) = secondary reflected light  VEXTRA(4) = eccentricity
! V(13) = stellar mass ratio         VEXTRA(5) = periastron longitude
! V(14) = tidal lead/lag angle (deg) VEXTRA(6) = reduced chi-squared
!-----------------------------------------------------------------------
! Version 1: Simplex minimisation algorithm  and new input / output used
! Version 2: Monte Carlo simulation and parameter perturbation algorithm
! Version 3: Adjustments to Monte Carlo LD coeffs and input/output files
! Version 4: Now solves for sum of radii, convergence criterion modified
! Version 5: Added TASK0 to find LD and GD coeffs.   Minor modifications
! Version 6: Reflection and  scale  factor  can all be fixed or adjusted
! Version 7: Can use either (e,w) or (ecosw,esinw); SFACT modified to be
!            in magnitudes; observ'l errors found; spherical star option
! Version 8: Bootstrapping error analysis algorithm added and output mod
! Version 9: Command-line arguments allowed, Monte Carlo without parame-
!            eter kicking option, and fitting for period and Tzero added
! Version 10: Can now have 9999 datapoints. Whole code now in magnitudes
! Version 11: Bug fixes, tasks renumbered,  added sigma clipping and the
!             global fit procedures, but not thoroughly tested these yet
! Version 12: Nonlinear limb darkening law, fitting for times of minimum
!             light, and FMAX corrections included (from Alvaro Gimenez)
! Version 13: Removed  BILINEAR  and modified  TASK1  to just call JKTLD
!             Modified input file and arrays  for the diff types of data
! Version 14: Fixed the requirement for inputtd INTRING to be an integer
!             Fixed formal errors (when observational ones not supplied)
! Last modified: 10th February 2006  (working on beta Aurigae)
!-----------------------------------------------------------------------
! Possible modifications in future:
! 1) extend to WD2003 and WINK
! 4) port to F90 or F95 to use long lines, modules, improved output
! 5) JVC suggests that formal fitting errors are x2 too small
! 6) incorporate change of omega (apsidal motion)
! 8) go through subroutine LIGHT to understand it all
!-----------------------------------------------------------------------
! Miscellaneous notes:
! 1) Phase shift has been redefined compared to original EBOP so that
!    it corresponds directly to the phase of the primary minimum.
! 2) MRQMIN adjusts coeffs only if VARY (called 'ia' in MRQMIN) is 1.
! 3) If VARY=2 then the parameter is fixed during the initial fit but is
!    perturbed by a set amount (flat distribution) for later analyses.
! 4) If VARY(11) and/or  VARY(12) are -1 then V(11) and/or V(12) are
!    calculated from the system geometry; if 0 they are fixed at the
!    input value and if 1 are freely adjusted to best fit.
! 5) If the mass ratio is <= 0 then both stars are assumed to be spheres
! 6) If ecosw > 5.0 then (ecosw,esinw) will be taken to be  (10+e,omega)
!    and fitting will occur using e and omega as parameters. e and omega
!    can be strongly correlated,  but this is useful if  e  is known but
!    omnega is not;   this can happen for EBs exhibiting apsidal motion.
! 7) Observational errors are looked for in the input light curve file.
!    If they are not found then equal weight is given to each point.
! 8) Nonlinear LD is now supported for the two-coefficient logarithmic,
!    quadratic and square-root laws. The type of law must be specified.
!    on input. Star B can also be forced to the same coeffs as star A.
!    BUT: normalisation for logarithmic not possible (not got equations)
! 9) Fitting for times of minimum light is directly possible. The cycle
!    numbers and times are inputted on lines immediately below all the
!    parameter lines in the input file.
!-----------------------------------------------------------------------
! Task numbers and purposes:
! (1) This outputs LD coefficients for a given Teff, logg, [M/H], Vmicro
! (2) This outputs a model light curve for fixed input parameters.
! (3) This fits a model to an observed light curve and outputs results.
! (4) This fits a model, rejects discrepant observations, and refits.
! (5) This does a pseudo-global minimisation by perturbing input params.
! (6) This investigates how different parameters vary around best fit.
! (7) This conducts bootstrapping simulations to find robust errors.
! (8) This conducts Monte Carlo simulations to find robust errors.
!=======================================================================
! Language: JKTEBOP is written in FORTRAN 77 using standard F77 syntax
! with a few straightforward extras, and compiled using g77.
! Possible non-F77 bits: == <= < > >= /= ! and part-line comments

!=======================================================================
!=======================================================================
      DOUBLEPRECISION FUNCTION GETPHASE (HJD,PERIOD,TZERO)
            ! Returns phase from given time and orbital ephemeris
      implicit none
      real*8 HJD,PERIOD,TZERO

      GETPHASE = (HJD - TZERO) / PERIOD
      GETPHASE = GETPHASE - int(GETPHASE)
      if ( GETPHASE < 0.0d0 ) GETPHASE = GETPHASE + 1.0d0

      END FUNCTION GETPHASE
!=======================================================================
      DOUBLEPRECISION FUNCTION GETMIN (TZERO,PERIOD,ECC,OMEGA,CYCLE)
            ! Returns time of minimum for given cycle and ephemeris.  If
            ! the orbit is circular thenthe cycle number is used without
            ! restriction so can refer to any phase. If the orbit is ec-
            ! centric then the cycle number should be integer (indicates
            ! primary minimum) or half-integer (secondary minimum).
      implicit none
      real*8 TZERO,PERIOD           ! IN: reference time, orbital period
      real*8 ECC,OMEGA              ! IN: orbital (e,w) or (ecosw,esinw)
      real*8 CYCLE                  ! IN: cycle number of minimum to use
      real*8 E,W                    ! LOCAL: eccentricity and peri.long.
      real*8 THETA,EE,TANEE         ! LOCAL: true and eccentric anomaly
      real*8 PSEP                   ! LOCAL: phase diff between minima
      real*8 PI,DEG2RAD,TINY        ! LOCAL: useful variables

      PI = atan(1.0d0) * 4.0d0
      DEG2RAD = 45.0d0 / atan(1.0d0)
      TINY = 1.0d-6

            ! First must deal with the possibility that e and omega are
            ! actually e*cos(omega) and e*sin(omega)

      if ( ECC > 9.0d0 ) then
        E = ECC - 10.0d0
        W = OMEGA
      else
        E = sqrt( OMEGA*OMEGA + ECC*ECC )
        W = atan2( ECC,OMEGA )
      end if

            ! If orbit is circular then simply use the orbital ephemeris
            ! If orbit is eccentric then must calculate the phase diffe-
            ! rence between two successive primary and secondary minima.

      if ( abs(E) < TINY ) then
        GETMIN = TZERO  +  PERIOD * CYCLE
      else
        THETA = 3.0d0*PI/2.0d0 - W
        TANEE = sqrt( (1.0d0-E) / (1.0d0+E) ) * tan(THETA/2.0d0)
        EE = atan(TANEE) * 2.0d0
        PSEP = (EE - E*sin(EE)) / (2.0d0 * PI)
        if ( PSEP < 0.0d0 ) PSEP = PSEP + 1.0d0

        if ( mod(abs(CYCLE),1.0d0) < TINY ) then       ! primary minimum
          GETMIN = TZERO + PERIOD * CYCLE
        else                                         ! secondary minimum
          GETMIN = TZERO + PERIOD * int(CYCLE)
          if ( CYCLE < 1.0d0 ) GETMIN = GETMIN - PERIOD
          GETMIN = GETMIN + PERIOD*PSEP
        end if
      end if

      END FUNCTION GETMIN

!=======================================================================
      !DOUBLEPRECISION FUNCTION GETMODEL (V,LDTYPE,TIME,TYPE,LA,LB)
      SUBROUTINE GETMODEL (V,LDTYPE,TIME,TYPE,LA,LB,FLUX)
            ! Output a predicted model value according to the parameters
            ! in array V. Precise meaning of the value depends on DTYPE.
            ! DTYPE=1  it outputs an EBOP magnitude for given time
            ! DTYPE=2  it outputs a light ratio for the given time
            ! DTYPE=3  outputs a time of eclipse for the given =CYCLE=
      implicit none
      real*8 V(22)                  ! IN: Photometric parameters
      integer LDTYPE(2)             ! IN: LD law type for the two stars
      real*8 TIME                   ! IN: The given TIME, PHASE or CYCLE
      integer TYPE                  ! IN: 1, 2 or 3 dep on wanted result
      real*8 LA,LB                  ! OUT: Light produced by each star
      real FMAG,LP,LS               ! LOCAL: LIGHT subroutine output
      real*8 ECC,OMEGA              ! LOCAL: eccentricity, perilongitude
      real*8 GETMIN,GETPHASE          ! FUNCTIONS
      real*8 FLUX

      if ( TYPE == 1 ) then
        CALL LIGHT(V,LDTYPE,real(GETPHASE(TIME,V(19),V(20))),FMAG,LP,LS)
        LA = dble(LP)
        LB = dble(LS)
        !GETMODEL = dble(FMAG)
        FLUX = dble(FMAG)

      else if ( TYPE == 2 ) then
        CALL LIGHT(V,LDTYPE,real(GETPHASE(TIME,V(19),V(20))),FMAG,LP,LS)
        !GETMODEL = dble(LS) / dble(LP)
        FLUX = dble(LS)/dble(LP)

      else if ( TYPE == 3 ) then
        if ( V(7) > 5.0d0 ) then
          ECC = V(7) - 10.0d0
          OMEGA = V(8)
        else
          ECC = sqrt(V(7)**2 + V(8)**2)
          OMEGA = atan2(V(8),V(7)) * 45.0d0 / atan(1.0d0)
          if ( OMEGA < 0.0d0 ) OMEGA = OMEGA + 360.0d0
        end if
        !GETMODEL = GETMIN (V(20),V(19),ECC,OMEGA,TIME)
        FLUX = GETMIN (V(20),V(19),ECC,OMEGA,TIME)

      else
        !GETMODEL = -100.0d0
        FLUX = -100.0d0
        print*,"### ERROR: wrong datatype asked for in GETMODEL:",TYPE
        STOP
      end if

      !END FUNCTION GETMODEL
      END SUBROUTINE GETMODEL
!=======================================================================

      SUBROUTINE BIAX (R,Q,A,B,EPS)
            ! EBOP subroutine to calculate biaxial ellipsoid dimensions
            ! and oblateness for each star after Chandrasekhar (1933).

      if ( Q <= 0.0 )  then
        A = R
        B = R
        EPS = 0.0
      else
        A = R * ( 1.0 + (1.0 + 7.0*Q)/6.0 * R**3.0)
        B = R * ( 1.0 + (1.0 - 2.0*Q)/6.0 * R**3.0)
        EPS = (A - B) / A
        B=( (1.0 - EPS) * R**3.0) ** (1.0/3.0)
        A = B / (1.0 - EPS)
      end if

      END SUBROUTINE BIAX
!=======================================================================
      SUBROUTINE GETEW (ECOSW,ESINW,E,W)
            ! EBOP subroutine to calculate e and w from e(cos)w e(sin)w

      if ( ECOSW == 0.0  .and.  ESINW == 0.0 ) then
        E = 0.0
        W = 0.0
      else
        W = atan2( ESINW,ECOSW )
        E = sqrt( ESINW*ESINW + ECOSW*ECOSW )
        W = W * 180.0 / 3.1415926536
      end if

      END SUBROUTINE GETEW
!=======================================================================
      SUBROUTINE LIGHT (V,LDTYPE,PHASE,FMAG,LP,LS)
      real*8 V(22)
      REAL LP,LS,LECL,LE
      real LD1U,LD2U          ! linear LD coeff for each star
      real LD1Q,LD2Q          ! quadratic LD coeff for each star
      real LD1S,LD2S          ! square-root LD coeff for each star
      real LD1L,LD2L          ! logarithmic LD coeff for each star
      real LDU,LDQ,LDS,LDL    ! LD coeffs for the star in question
      integer LDTYPE(2)       ! LD law type for both stars
      integer GIMENEZ         ! 1 to use original FMAX calculations
                              ! 2 to use Gimenez' modified calcs
                              ! 3 to use Gimenez' calcs for nonlinear LD
      data GIMENEZ / 3 /

      data PI,TWOPI,RAD / 3.1415926536E0,6.28318531E0,0.0174532925E0 /
!
!        DETERMINE PRIMARY AND SECONDARY BIAXIAL DIMENSIONS
!        USE SPERICAL RADII FOR THE COMPUTATION OF ECLIPSE FUNCTIONS
!        USE OBLATENESSES FOR THE COMPUTATION OF THE OUTSIDE ECLIPSE
!        PHOTOMETRIC VARIATIONS WITH LIMB AND GRAVITY DARKENING
!
      BS     = real(V(1))
      RP     = real(V(2)/(1.0d0+V(3)))
      RATIO  = real(V(3))
      FI     = real(V(6))

      if ( V(7) > 5.0d0 ) then
        ECOSW = real( (V(7)-10.0d0) * cos(V(8)/57.2957795d0) )
        ESINW = real( (V(7)-10.0d0) * sin(V(8)/57.2957795d0) )
      else
        ECOSW  = real(V( 7))
        ESINW  = real(V( 8))
      end if

      YP     = real(V( 9))
      YS     = real(V(10))
      SP     = real(V(11))
      SS     = real(V(12))
      Q      = real(V(13))
      TANGL  = real(V(14))
      EL     = real(V(15))
      DPH    = 1.0 - real(V(16))
      SFACT  = real(V(17))!10.0d0**(-0.4*V(17)))
      DGAM   = real(V(18))

      if ( Q <= 0.0 ) then
        CALL BIAX (RP,0.0,RPA,RPB,EP)
        RS=RP*RATIO
        CALL BIAX (RS,0.0,RSA,RSB,ES)
      else
        CALL BIAX (RP,Q,RPA,RPB,EP)
        RS=RP*RATIO
        CALL BIAX (RS,1.0/Q,RSA,RSB,ES)
      end if

      LD1U = real(V(4))         ! linear term
      LD2U = real(V(5))
      LD1L = 0.0                ! log term
      LD2L = 0.0
      LD1S = 0.0                ! sqrt term
      LD2S = 0.0
      LD1Q = 0.                 ! quadratic term
      LD2Q = 0.0
      if ( LDTYPE(1) == 2 ) LD1L = real(V(21))
      if ( LDTYPE(1) == 3 ) LD1S = real(V(21))
      if ( LDTYPE(1) == 4 ) LD1Q = real(V(21))
      if ( LDTYPE(2) == 2 ) LD2L = real(V(22))
      if ( LDTYPE(2) == 3 ) LD2S = real(V(22))
      if ( LDTYPE(2) == 4 ) LD2Q = real(V(22))
      if ( LDTYPE(2) == 0 ) then
        LD2U = LD1U
        LD2L = LD1L
        LD2S = LD1S
        LD2Q = LD1Q
      end if

!        write(23,*),ldtype
!        write(23,*),ld1u,ld2u
!        write(23,*),ld1l,ld2l
!        write(23,*),ld1s,ld2s
!        write(23,*),ld1q,ld2q
!        write(23,*)," "

!
!        CORRECT THE OBSERVED PHASE FOR ANY EPOCH ERROR IN EPHEMERIS
!
      THETA=PHASE+DPH
!
      SINI  = SIN(FI*RAD)
      SINI2 = SINI*SINI
      COSI2 = 1.0E0  - SINI2
!
!        TRANSLATE TIDAL LEAD/LAG ANGLE TO RADIANS
      TANGR=TANGL*RAD
!
!     EQUATION 9
!        CONVERT PHASE TO RADIANS
      FMN=THETA*TWOPI
!
!        GET CURRENT VALUES OF E, AND W
      CALL GETEW (ECOSW,ESINW,E,W)
!
!        TEST FOR CIRCULAR ORBIT
      IF (E)   17,20,17
   20 COSVW=COS(FMN)
      SINVW=SIN(FMN)
      RV=1.0E0
      GO TO 25
!
!        SOLUTION OF KEPLER'S EQUATION BY DIFFERENTIAL CORRECTIONS
!        (NON-ZERO ECCENTRICITY ONLY . . . )
!
!     EQUATION 6
!
   17 OMEGA = 450.0E0  - W
   23 IF (OMEGA - 360.0E0)         22,21,21
   21 OMEGA = OMEGA - 360.0E0
      GO TO 23
   22 OMEGA = OMEGA*RAD
!        SINE AND COSINE OF OMEGA
      COSW=COS(OMEGA)
      SINW=SIN(OMEGA)
!
!        COMPUTE MEAN ANOMALY CORRECTION TO PHASE
!        CORRESPONDING TO V=OMEGA=90-W
!        AT WHICH PHASE COS(V-OMEGA)=1
      E0=ATAN2(SQRT(1.0E0-E*E)*SINW,COSW+E)
!
!        MEAN ANOMALY OF MID-PRIMARY ECLIPSE
      FMA0=E0-E*SIN(E0)
!
!        MEAN ANOMALY
      FMA=FMN+FMA0
!     FIRST APPROXIMATION OF ECCENTRIC ANOMALY
      EA=FMA+E*SIN(FMA)
!
      DO 10 J=1,15
!        EVALUATE SINE AND COSINE OF ECCENTRIC ANOMALY
      SINE=SIN(EA)
      COSE=COS(EA)
      DENOM=1.0E0-E*COSE
      DISC=FMA-EA+E*SINE
      EA=EA+DISC/DENOM
!        TEST FOR CONVERGENCE
      IF (ABS(DISC) - 2.0E-05)     15,15,10
   10 CONTINUE
!
!
!        EVALUATE SINE AND COSINE OF TRUE ANOMALY
   15 COSV=(COSE-E)/DENOM
      SINV=SINE*SQRT(1.0E0-E*E)/DENOM
!
!        RADIUS VECTOR
      RV = (1.0E0-E*E)/(1.0E0+E*COSV)
!
!        THE PHOTOMETRIC PHASE ARGUMENT IN TERMS OF ORBIT PARAMETERS
!        VW = V-OMEGA
      COSVW=COSV*COSW+SINV*SINW
      SINVW=SINV*COSW-COSV*SINW
!
   25 COS2=COSVW*COSVW
      SIN2=1.0E0-COS2
!
      CSVWT=COS(TANGR)*COSVW-SIN(TANGR)*SINVW
!
!
!        PHOTOMETRIC EFFECTS
!
!

!-----------------------------------------------------------------------
! Alvaro Gimenez has corrected the treatment of stellar shapes.
! This newer treatment can be used by putting GIMENEZ=2
! His teatment for nonlinear LD can be used by putting GIMENEZ=3
!-----------------------------------------------------------------------
! This whole thing affects only the brightness normalisation of the two
! eclipsing stars: any problems here wil affect the radiative parameters
! but not the geometric parameters (radii, inclination etc).
!-----------------------------------------------------------------------
      FMAXP = 0.0
      FMAXS = 0.0
      DELTP = 0.0
      DELTS = 0.0
      SHORT = 0.0
      if ( GIMENEZ == 1 ) then

!       FMAXP=((1.0E0-UP)+0.666666667E0*UP*(1.0E0+0.2E0*EP))  ! Original
!      1      *(1.0E0+3.0E0*YP*EP)/(1.0E0-EP)                 ! lines
!       FMAXS=((1.0E0-US)+0.666666667E0*US*(1.0E0+0.2E0*ES))  ! if
!      1      *(1.0E0+3.0E0*YS*ES)/(1.0E0-ES)                 ! stars
!       DELTP=(15.0E0+UP)/(15.0E0-5.0E0*UP)*(1.0E0+YP)*EP     ! are
!       DELTS=(15.0E0+US)/(15.0E0-5.0E0*US)*(1.0E0+YS)*ES     ! oblate
!       SHORT=SINI2*CSVWT*CSVWT

!    26 FMAXP=1.0E0-UP/3.0E0              ! Original
!       FMAXS=1.0E0-US/3.0E0              ! lines if
!       DELTP=0.0E0                       ! stars
!       DELTS=0.0E0                       ! are
!       SHORT=0.0                         ! spherical

        if ( Q >= 0.0 ) then
          FMAXP=((1.0E0-LD1U)+0.666666667E0*LD1U*(1.0E0+0.2E0*EP))*(1.0E0+3.0E0*YP*EP)/(1.0E0-EP)
          FMAXS=((1.0E0-LD2U)+0.666666667E0*LD2U*(1.0E0+0.2E0*ES))*(1.0E0+3.0E0*YS*ES)/(1.0E0-ES)
          DELTP=(15.0E0+LD1U)/(15.0E0-5.0E0*LD1U)*(1.0E0+YP)*EP
          DELTS=(15.0E0+LD2U)/(15.0E0-5.0E0*LD2U)*(1.0E0+YS)*ES
          SHORT=SINI2*CSVWT*CSVWT
        else
          FMAXP=1.0E0-LD1U/3.0E0
          FMAXS=1.0E0-LD2U/3.0E0
          DELTP=0.0E0
          DELTS=0.0E0
          SHORT=0.0
        end if
!-----------------------------------------------------------------------

      else if ( GIMENEZ == 2 ) then

!       FMAXP=(1.0E0-UP*(1.0E0-2.0E0/5.0E0*EP)/3.0E0+YP*EP
!      1      *(3.0E0-13.0E0/15.0E0*UP))/(1.0E0-EP)
!       FMAXS=(1.0E0-US*(1.0E0-2.0E0/5.0E0*ES)/3.0E0+YS*ES
!      1      *(3.0E0-13.0E0/15.0E0*US))/(1.0E0-ES)
!       DELTP=(15.0E0+UP)/(15.0E0-5.0E0*UP)*(1.0E0+YP)*EP
!       DELTS=(15.0E0+US)/(15.0E0-5.0E0*US)*(1.0E0+YS)*ES
!       SHORT=SINI2*CSVWT*CSVWT

!    26 FMAXP=1.0E0-UP/3.0E0              ! Original
!       FMAXS=1.0E0-US/3.0E0              ! lines if
!       DELTP=0.0E0                       ! stars
!       DELTS=0.0E0                       ! are
!       SHORT=0.0                         ! spherical

        if ( Q >= 0.0 ) then
          FMAXP=(1.0E0-LD1U*(1.0E0-2.0E0/5.0E0*EP)/3.0E0+YP*EP*(3.0E0-13.0E0/15.0E0*LD1U))/(1.0E0-EP)
          FMAXS=(1.0E0-LD2U*(1.0E0-2.0E0/5.0E0*ES)/3.0E0+YS*ES*(3.0E0-13.0E0/15.0E0*LD2U))/(1.0E0-ES)
          DELTP=(15.0E0+LD1U)/(15.0E0-5.0E0*LD1U)*(1.0E0+YP)*EP
          DELTS=(15.0E0+LD2U)/(15.0E0-5.0E0*LD2U)*(1.0E0+YS)*ES
          SHORT=SINI2*CSVWT*CSVWT
        else
          FMAXP=1.0E0-LD1U/3.0E0
          FMAXS=1.0E0-LD2U/3.0E0
          DELTP=0.0E0
          DELTS=0.0E0
          SHORT=0.0
        end if
!-----------------------------------------------------------------------

      else if ( GIMENEZ == 3 ) then
! And this is Gimenez's code for including nonlinear LD. He includes
! the linear (UP), quadratic (UP, U2P) and square-root (UP, U3P) laws.

!      FMAXP=1.0E0-UP*(1.0E0-2.0E0*EP/5.0E0)/3.0E0-
!     1      U2P*(1.0E0-3.0E0*EP/5.0E0)/6.0E0-
!     1      U3P*(1.0E0-4.0E0*EP/9.0E0)/5.0E0+2.0E0*YP*EP
!     1      *(1.5E0-13.0E0*UP/30.0E0-U2P/5.0E0-23.0E0*U3P/90.0E0)
!      FMAXP=FMAXP/(1.0E0-EP)
!      FMINP=1.0E0-UP*(1.0E0+4.0E0*EP/5.0E0)/3.0E0-
!     1      U2P*(1.0E0+6.0E0*EP/5.0E0)/6.0E0-
!     1      U3P*(1.0E0+8.0E0*EP/9.0E0)/5.0E0+2.0E0*YP*EP
!     1      *(1.0E0-7.0E0*UP/15.0E0-4.0E0*U2P/15.0E0-13.0E0*U3P/45.0E0)
!      FMINS=1.0E0-US*(1.0E0+4.0E0*ES/5.0E0)/3.0E0-
!     1      U2S*(1.0E0+6.0E0*ES/5.0E0)/6.0E0-
!     1      U3S*(1.0E0+8.0E0*ES/9.0E0)/5.0E0+2.0E0*YS*ES
!     1      *(1.0E0-7.0E0*US/15.0E0-4.0E0*U2S/15.0E0-13.0E0*U3S/45.0E0)
!      FMAXS=1.0E0-US*(1.0E0-2.0E0*ES/5.0E0)/3.0E0-
!     1      U2S*(1.0E0-3.0E0*ES/5.0E0)/6.0E0-
!     1      U3S*(1.0E0-4.0E0*ES/9.0E0)/5.0E0+2.0E0*YS*ES
!     1      *(1.5E0-13.0E0*US/30.0E0-U2S/5.0E0-23.0E0*U3S/90.0E0)
!      FMAXS=FMAXS/(1.0E0-ES)
!      DELTP=1.0E0-FMINP/FMAXP
!      DELTS=1.0E0-FMINS/FMAXS
!      SHORT=SINI2*CSVWT*CSVWT

!   26 FMAXP=1.0E0-UP/3.0E0-U2P/6.0E0-U3P/5.0E0
!      FMAXS=1.0E0-US/3.0E0-U2S/6.0E0-U3S/5.0E0
!      DELTP=0.0E0
!      DELTS=0.0E0
!      SHORT=0.0

        if ( Q >= 0.0d0 ) then
          FMAXP=1.0E0-LD1U*(1.0E0-2.0E0*EP/5.0E0)/3.0E0-LD1Q*(1.0E0-3.0E0*EP/5.0E0)/6.0E0-LD1S*(1.0E0-4.0E0*EP/9.0E0)/5.0E0+2.0E0*YP*EP*(1.5E0-13.0E0*LD1U/30.0E0-LD1Q/5.0E0-23.0E0*LD1S/90.0E0)
          FMAXP=FMAXP/(1.0E0-EP)
          FMINP=1.0E0-LD1U*(1.0E0+4.0E0*EP/5.0E0)/3.0E0-LD1Q*(1.0E0+6.0E0*EP/5.0E0)/6.0E0-LD1S*(1.0E0+8.0E0*EP/9.0E0)/5.0E0+2.0E0*YP*EP*(1.0E0-7.0E0*LD1U/15.0E0-4.0E0*LD1Q/15.0E0-13.0E0*LD1S/45.0E0)
          FMINS=1.0E0-LD2U*(1.0E0+4.0E0*ES/5.0E0)/3.0E0-LD2Q*(1.0E0+6.0E0*ES/5.0E0)/6.0E0-LD2S*(1.0E0+8.0E0*ES/9.0E0)/5.0E0+2.0E0*YS*ES*(1.0E0-7.0E0*LD2U/15.0E0-4.0E0*LD2Q/15.0E0-13.0E0*LD2S/45.0E0)
          FMAXS=1.0E0-LD2U*(1.0E0-2.0E0*ES/5.0E0)/3.0E0-LD2Q*(1.0E0-3.0E0*ES/5.0E0)/6.0E0-LD2S*(1.0E0-4.0E0*ES/9.0E0)/5.0E0+2.0E0*YS*ES*(1.5E0-13.0E0*LD2U/30.0E0-LD2Q/5.0E0-23.0E0*LD2S/90.0E0)
          FMAXS=FMAXS/(1.0E0-ES)
          DELTP=1.0E0-FMINP/FMAXP
          DELTS=1.0E0-FMINS/FMAXS
          SHORT=SINI2*CSVWT*CSVWT
        else
          FMAXP=1.0E0-LD1U/3.0E0-LD1Q/6.0E0-LD1S/5.0E0
          FMAXS=1.0E0-LD2U/3.0E0-LD2Q/6.0E0-LD2S/5.0E0
          DELTP=0.0E0
          DELTS=0.0E0
          SHORT=0.0
        end if
!----------------------------------------------------------------------
      end if
!----------------------------------------------------------------------
! Complete original code before the above messing:
! C
! C
! C        PHOTOMETRIC EFFECTS
! C
! C
! C        TEST FOR SIMPLE CASE OF TWO SPHERICAL STARS
!       IF (EP .EQ. 0.  .AND.  ES .EQ. 0.)   GO TO 26
! C
! C        EITHER OR BOTH STARS ARE OBLATE
! C
!       FMAXP=((1.0E0-UP)+0.666666667E0*UP*(1.0E0+0.2E0*EP))
!      1      *(1.0E0+3.0E0*YP*EP)/(1.0E0-EP)
!       FMAXS=((1.0E0-US)+0.666666667E0*US*(1.0E0+0.2E0*ES))
!      1      *(1.0E0+3.0E0*YS*ES)/(1.0E0-ES)
! C        CHANGE IN INTENSITY RATIO DUE TO OBLATENESS RELATED VARIABLES
! C        FROM QUADRATURE TO MINIMUM
! C        FACE ON TO END ON
!       DELTP=(15.0E0+UP)/(15.0E0-5.0E0*UP)*(1.0E0+YP)*EP
!       DELTS=(15.0E0+US)/(15.0E0-5.0E0*US)*(1.0E0+YS)*ES
! C        FORE-SHORTENING FUNCTION OF OBLATENESS
!       SHORT=SINI2*CSVWT*CSVWT
!       GO TO 27
! C
! C        BOTH STARS ARE SPHERICAL
! C
!    26 FMAXP=1.0E0-UP/3.0E0
!       FMAXS=1.0E0-US/3.0E0
!       DELTP=0.0E0
!       DELTS=0.0E0
!       SHORT=0.0
!----------------------------------------------------------------------

!
!        UN-NORMALIZED BRIGHTNESS OF STELLAR COMPONENTS AT QUADRATURE
   27 OP=PI*RPB*RPB*FMAXP
      OS=PI*RSB*RSB*FMAXS*BS
!        THE NORMALIZING FACTOR
      OTOT=OP+OS
!        BRIGHTNESS CONTRIBUTION FROM EACH COMPONENT
      LP=OP/OTOT*(1.0E0-DELTP*SHORT)
      LS=OS/OTOT*(1.0E0-DELTS*SHORT)
!
!        REFLECTION AND RERADIATION EQUATION
      IF (SP .EQ. 0.0E0  .AND.  SS .EQ. 0.0E0)   GO TO 28
      HEAT=SINI*COSVW
      HEAT2=0.5E0+0.5E0*HEAT*HEAT
      DLP=SP*(HEAT2+HEAT)
      DLS=SS*(HEAT2-HEAT)
      GO TO 29
   28 DLP=0.0E0
      DLS=0.0E0
!
!        WHICH ECLIPSE COULD THIS BE
   29 IF (COSVW)         40,40,30
!
!     PRIMARY ECLIPSE
!
   30 R1 = RP
      R2 = RS
!----------------------------------------------------------------------!
! JKT mod (10/8/2006): the line these replaced was      UU = UP        !
!----------------------------------------------------------------------!
      LDU = LD1U                                                       !
      LDL = LD1L                                                       !
      LDS = LD1S                                                       !
      LDQ = LD1Q                                                       !
!----------------------------------------------------------------------!
      LE=LP
      DLE=DLP
      GO TO 60
!
!
!     SECONDARY ECLIPSE
!
   40 R1 = RS
      R2 = RP
!-----------------------------------------------------------------------
! JKT mod (10/8/2006): the line these replaced was      UU = US        !
!----------------------------------------------------------------------!
      LDU = LD2U                                                       !
      LDL = LD2L                                                       !
      LDS = LD2S                                                       !
      LDQ = LD2Q                                                       !
!----------------------------------------------------------------------!
      LE=LS
      DLE=DLS
!
   60 SUM = 0.0E0
      ALAST = 0.0E0
      AREA=0.0E0
!
!     EQUATION  5
!
      DD = SINVW*SINVW + COSVW*COSVW*COSI2
      IF (DD .LE. 1.0E-06)  DD=0.0
      DD = DD*RV*RV
      D = SQRT(ABS(DD))
      R22 = R2*R2
!
!     EQUATION 17
!
      GAMN = 90.01E0*RAD
      DGAMA = DGAM*RAD
      DGM = DGAMA/2.0E0
      RK = 0.0E0
      GAM = 0.0E0
   50 GAM = GAM + DGAMA
!        HAS LIMIT OF INTEGRATION BEEN REACHED
      IF (GAM - GAMN)              48,48,49
!
   48 RR = R1*SIN(GAM)
      R12 = RR*RR
!
      AA = 0.0E0
!        ARE THE PROJECTED DISKS CONCENTRIC
      IF (D)                       405,406,405
  406 IF (RR - R2)                 230,230,403
  403 IF (RK - R2)                 404, 49, 49
  404 AA = PI*R22
      GO TO 215
!        TEST FOR NO ECLIPSE
  405 IF (D-R1-R2)                 240,216,216
  216 SUM = 0.0E0
      GO TO 49
!        DECIDE WHICH AREA EQUATIONS FOR NON-CONCENTRIC ECLIPSE
  240 IF (D-RR-R2)                 245,215,215
  245 IF (D-R2+RR)                 230,230,250
  250 IF (R1-R2)                   255,255,280
  255 IF (DD-R22+R12)              205,210,210
  280 IF (D-RR+R2)                 290,260,260
  260 IF (RR-R2)                   255,255,265
  265 IF (DD-R12+R22)              270,210,210
!
!     EQUATION 12
!
  270 S1 = ABS((R12 - R22 - DD)*0.5E0/D)
      A1 = ABS(R2-S1)
      B2 = ABS(RR-S1-D  )
      AA=PI*R22-(R22*ACOS((R2-A1)/R2)- (R2-A1)*SQRT(2.0E0*R2*A1-A1*A1))+R12*ACOS((RR-B2)/RR)-(RR-B2)*SQRT(2.0E0*RR*B2-B2*B2)
      GO TO 215
!
  290 IF (R1 - R2 - D)             260,260,295
  295 IF (RK - R2 - D)             300,215,215
  300 RR = R2 + D
      R12 = RR*RR
      GAMN = 0.0E0
      GO TO 260
!
  230 AA = PI*R12
      GO TO 215
!
!     EQUATION 10
!
  205 S = ABS((R12 - R22 + DD)*0.5E0/D)
      A = ABS(RR-S)
      B1 = ABS(R2-S-D)
      A1 = R12*ACOS((RR-A)/RR) - (RR-A)*SQRT(2.0E0*RR*A - A*A)
      AB1 = R22*ACOS((R2-B1)/R2) - (R2-B1)*SQRT(2.0E0*R2*B1-B1*B1)
      AA = PI*R12 - A1 + AB1
      GO TO 215
!
!     EQUATION 1
!
  210 S = ABS((R12 - R22 + DD)*0.5E0/D)
      A = ABS(RR-S)
      B = ABS(S-D+R2)
      A1 = R12*ACOS((RR-A)/RR) - (RR-A)*SQRT(2.0E0*RR*A - A*A)
      AA1 = R22*ACOS((R2-B)/R2) - (R2-B)*SQRT(2.0E0*R2*B - B*B)
      AA = A1 + AA1
!
  215 DAREA = AA - ALAST
!----------------------------------------------------------------------!
! JKT modification (10/9/2006). The removed line was:                  !
!     SUM = SUM + DAREA*(1.0E0  - UU + UU*COS(GAM-DGM))                !
!----------------------------------------------------------------------!
      COSGAM = cos(GAM-DGM)                                            !
      SUM = SUM + DAREA*(1.0 - LDU*(1.0-COSGAM) - LDL*COSGAM*log(COSGAM) - LDS*(1.0-sqrt(COSGAM))  - LDQ*(1.0-COSGAM)**2)       !
!----------------------------------------------------------------------!
      ALAST = AA
      AREA = AREA + DAREA
!
      RK = RR
      GO TO 50
!
!        LIGHT LOSS FROM ECLIPSE
!
   49 ADISK = PI*R1*R1
!----------------------------------------------------------------------!
! JKT modification (10/9/2006).  See 1992A+A...259..227D for more info.!
! The factor 10.3616329185 was calculated independently by JKT         !
! The removed line was:           ALPHA = SUM/(ADISK*(1.0E0-UU/3.0E0)) !
!----------------------------------------------------------------------!
      ALPHA = 1.0 - LDU/3.0 + LDL/10.3616329185 - LDS/5.0 - LDQ/6.0    !
      ALPHA = SUM/(ADISK*ALPHA)                                        !
!----------------------------------------------------------------------!
      LECL = ALPHA*LE
      AREA = AREA/ADISK
      REFL=DLP+DLS-AREA*DLE
!
!        THEORETICAL INTENSITY WITH THIRD LIGHT AND QUADRATURE
!        SCALE FACTOR APPLIED
!
!      FLITE = ((LP+LS-LECL+REFL)*(1.0E0-EL)+EL)*SFACT

      FLITE = ((LP+LS-LECL+REFL)*(1.0E0-EL)+EL)
      FMAG = -2.5d0 * log10(FLITE) + SFACT
!
!      RETURN
!      END
      END SUBROUTINE LIGHT
