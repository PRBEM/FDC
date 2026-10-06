! ==============================================================================
! SUBROUTINE: Freq_Spectral
! ==============================================================================
!
! PURPOSE:
!   Returns the wave frequency band parameters (mean frequency, bandwidth,
!   and cutoff frequencies) at a given location and equatorial gyrofrequency,
!   for the spectral model specified by FreqUnit in module constants.
!
!   These parameters define the frequency band over which wave power is
!   integrated in the diffusion coefficient calculation. They are used by:
!     - Roots_finder: to check validity of resonant frequency roots
!     - PowerSpectralDensity: as normalisation bounds (GAUSSIAN case)
!     - WeightFunction: stored in module variables omega_m, d_omega,
!       omega_lc, omega_uc after this call
!
!   Supported models and their parameter dependencies:
!
!     Model         L    lambda   MLT   Kp   Sets omega_m/d_omega?
!     HZ                                     Yes (from input file, in Hz)
!     RATIO                                  Yes (as fraction of omega_e_eq)
!     2016HISS      X                        Yes (L-dependent)
!     NBBLBC             X              X    Yes (Kp and latitude dependent)
!     NBBUBC             X              X    Yes (Kp and latitude dependent)
!     2017CHORUS                   X         No  (omega_lc, omega_uc only)
!     2018LBCHORUS                 X         No  (omega_lc, omega_uc only)
!     2018UBCHORUS                 X         No  (omega_lc, omega_uc only)
!     2026LBCHORUS                           No  (omega_lc, omega_uc only)
!     2026UBCHORUS                           No  (omega_lc, omega_uc only)
!     EMIC_HBAND                             Yes (H+ band, normalized to Omega_p)
!     EMIC_HEBAND                            Yes (He+ band, normalized to Omega_O+)
!     EMIC_OBAND                             Yes (O+ band, normalized to Omega_O+)
!
! INPUTS:
!   L           - Magnetic L-shell [Earth radii, dimensionless]
!   Lambda      - Magnetic latitude [radians]
!   omega_e_eq  - Electron gyrofrequency at the equator [rad/s]
!
! OUTPUTS:
!   omega_m     - Mean wave frequency [rad/s]
!   d_omega     - Wave frequency bandwidth [rad/s]
!   omega_lc    - Lower frequency cutoff [rad/s]
!   omega_uc    - Upper frequency cutoff [rad/s]
!
! MODULE DEPENDENCIES:
!   constants - FreqUnit, omega_m_per, d_omega_per, o_lc, o_uc, Kp, MLT,
!               PI, RAD
!
! NOTES:
!   - For models that set only omega_lc and omega_uc (2017CHORUS,
!     2018LBCHORUS, 2018UBCHORUS, 2026LBCHORUS, 2026UBCHORUS) the spectral
!     shape is defined entirely within PowerSpectralDensity via its
!     polynomial fit, so omega_m and d_omega are not needed. They are set
!     to -1 here rather than left undefined, so that an accidental use
!     shows up as a negative frequency instead of as whatever happened to
!     be in memory.
!   - All frequencies are in rad/s. The HZ case converts from Hz by
!     multiplying by 2*PI. The 2016HISS case uses Hz * 2*PI directly.
!   - The EMIC bands follow Albert (2003), Sec. 4. The He+ band is
!     normalised to the O+ gyrofrequency (3*Omega_O+ lies between Omega_O+
!     and Omega_He+ = 4*Omega_O+), as is the O+ band.
!   - omega_uc = omega_e_eq (2018UBCHORUS) puts the upper band edge at the
!     electron gyrofrequency, where the cold-plasma resonance cone closes
!     and N(omega) cannot be evaluated.
!
! ==============================================================================

subroutine Freq_Spectral(L, Lambda, omega_e_eq, omega_m, d_omega, omega_lc, omega_uc)
use constants
implicit none

    real*8, intent(in)  :: L, Lambda, omega_e_eq
    real*8, intent(out) :: omega_m, d_omega, omega_lc, omega_uc

    real*8 :: lambda_deg

    lambda_deg = Lambda / RAD

    select case(FreqUnit)

        ! ----------------------------------------------------------------------
        ! HZ: fixed Gaussian spectrum, parameters given in Hz in input file.
        ! Converted to rad/s by multiplying by 2*PI.
        ! ----------------------------------------------------------------------
        case('HZ')
            d_omega  = d_omega_per * 2.0 * PI
            omega_m  = omega_m_per * 2.0 * PI
            omega_uc = o_uc        * 2.0 * PI
            omega_lc = o_lc        * 2.0 * PI

        ! ----------------------------------------------------------------------
        ! RATIO: fixed Gaussian spectrum, parameters given as fractions of
        ! omega_e_eq in input file.
        ! ----------------------------------------------------------------------
        case('RATIO')
            d_omega  = d_omega_per * omega_e_eq
            omega_m  = omega_m_per * omega_e_eq
            omega_uc = o_uc        * omega_e_eq
            omega_lc = o_lc        * omega_e_eq

        ! ----------------------------------------------------------------------
        ! 2016HISS: L-dependent hiss frequency band.
        ! All frequencies in Hz * 2*PI = rad/s.
        ! ----------------------------------------------------------------------
        case('2016HISS')
            if (L <= 2.5) then
                d_omega = 290.0 * 2.0*PI; omega_m = 282.0 * 2.0*PI
                omega_uc = 4000.0 * 2.0*PI; omega_lc = 20.0 * 2.0*PI
            else if (L <= 3.5) then
                d_omega = 160.0 * 2.0*PI; omega_m = 252.0 * 2.0*PI
                omega_uc = 4000.0 * 2.0*PI; omega_lc = 20.0 * 2.0*PI
            else if (L <= 4.5) then
                d_omega = 110.0 * 2.0*PI; omega_m = 178.0 * 2.0*PI
                omega_uc = 4000.0 * 2.0*PI; omega_lc = 20.0 * 2.0*PI
            else
                d_omega = 90.0  * 2.0*PI; omega_m = 140.0 * 2.0*PI
                omega_uc = 4000.0 * 2.0*PI; omega_lc = 20.0 * 2.0*PI
            end if

        ! ----------------------------------------------------------------------
        ! NBBLBC: lower-band chorus frequency band.
        ! Kp and latitude dependent. omega_m and d_omega set as fractions of
        ! omega_e_eq.
        ! ----------------------------------------------------------------------
        case('NBBLBC')
            if (Kp == 0.0) then
                if (lambda_deg <= 5.0) then
                    omega_m = 0.28*omega_e_eq; d_omega = 0.07*omega_e_eq
                else if (lambda_deg <= 10.0) then
                    omega_m = 0.25*omega_e_eq; d_omega = 0.13*omega_e_eq
                else
                    omega_m = 0.22*omega_e_eq; d_omega = 0.10*omega_e_eq
                end if
            else if (Kp == 1.0) then
                if (lambda_deg <= 5.0) then
                    omega_m = 0.28*omega_e_eq; d_omega = 0.12*omega_e_eq
                else if (lambda_deg <= 10.0) then
                    omega_m = 0.30*omega_e_eq; d_omega = 0.08*omega_e_eq
                else
                    omega_m = 0.26*omega_e_eq; d_omega = 0.07*omega_e_eq
                end if
            else
                if (lambda_deg <= 5.0) then
                    omega_m = 0.35*omega_e_eq; d_omega = 0.15*omega_e_eq
                else if (lambda_deg <= 10.0) then
                    omega_m = 0.30*omega_e_eq; d_omega = 0.15*omega_e_eq
                else
                    omega_m = 0.25*omega_e_eq; d_omega = 0.10*omega_e_eq
                end if
            end if
            omega_lc = 0.10*omega_e_eq
            omega_uc = 0.50*omega_e_eq

        ! ----------------------------------------------------------------------
        ! NBBUBC: upper-band chorus frequency band.
        ! Kp and latitude dependent.
        ! ----------------------------------------------------------------------
        case('NBBUBC')
            if (Kp == 0.0) then
                if (lambda_deg <= 5.0) then
                    omega_m = 0.54*omega_e_eq; d_omega = 0.08*omega_e_eq
                else
                    omega_m = 0.56*omega_e_eq; d_omega = 0.06*omega_e_eq
                end if
            else if (Kp == 1.0) then
                if (lambda_deg <= 5.0) then
                    omega_m = 0.54*omega_e_eq; d_omega = 0.10*omega_e_eq
                else
                    omega_m = 0.57*omega_e_eq; d_omega = 0.04*omega_e_eq
                end if
            else
                if (lambda_deg <= 5.0) then
                    omega_m = 0.53*omega_e_eq; d_omega = 0.10*omega_e_eq
                else
                    omega_m = 0.58*omega_e_eq; d_omega = 0.04*omega_e_eq
                end if
            end if
            omega_lc = 0.50*omega_e_eq
            omega_uc = 0.65*omega_e_eq

        ! ----------------------------------------------------------------------
        ! 2017CHORUS: only sets band edges; omega_m and d_omega undefined.
        ! ----------------------------------------------------------------------
        case('2017CHORUS')
            omega_lc = 0.05*omega_e_eq
            omega_uc = 0.50*omega_e_eq
            omega_m  = -1.0d0      ! unused: spectral shape from PowerSpectralDensity
            d_omega  = -1.0d0

        ! ----------------------------------------------------------------------
        ! 2018LBCHORUS: lower-band chorus band edges
        ! ----------------------------------------------------------------------
        case('2018LBCHORUS')
            omega_lc = 0.05*omega_e_eq
            omega_uc = 0.50*omega_e_eq
            omega_m  = -1.0d0      ! unused: spectral shape from PowerSpectralDensity
            d_omega  = -1.0d0

        ! ----------------------------------------------------------------------
        ! 2018UBCHORUS: upper-band chorus band edges
        ! ----------------------------------------------------------------------
        case('2018UBCHORUS')
            omega_lc = 0.50*omega_e_eq
            omega_uc = 1.00*omega_e_eq
            omega_m  = -1.0d0      ! unused: spectral shape from PowerSpectralDensity
            d_omega  = -1.0d0

        ! ----------------------------------------------------------------------
        ! 2026LBCHORUS: lower-band chorus band edges
        ! ----------------------------------------------------------------------
        case('2026LBCHORUS')
            omega_lc = 0.05*omega_e_eq
            omega_uc = 0.50*omega_e_eq
            omega_m  = -1.0d0      ! unused: spectral shape from PowerSpectralDensity
            d_omega  = -1.0d0

        ! ----------------------------------------------------------------------
        ! 2026UBCHORUS: upper-band chorus band edges
        ! ----------------------------------------------------------------------
        case('2026UBCHORUS')
            omega_lc = 0.50*omega_e_eq
            omega_uc = 0.80*omega_e_eq
            omega_m  = -1.0d0      ! unused: spectral shape from PowerSpectralDensity
            d_omega  = -1.0d0

        ! ----------------------------------------------------------------------
        ! EMIC_HBAND: hydrogen-band EMIC waves Albert (2003) Sec. 4. omega_m, d_omega normalized to the
        ! proton gyrofrequency
        ! ----------------------------------------------------------------------
        case('EMIC_HBAND')
            block
                real*8 :: Omega_p_eq
                Omega_p_eq = omega_e_eq * me / mi
                omega_m = 0.6d0 * Omega_p_eq
                d_omega = 0.1d0 * Omega_p_eq
                omega_lc = omega_m - d_omega
                omega_uc = omega_m + d_omega
            end block

        ! ----------------------------------------------------------------------
        ! EMIC_HEBAND: helium-band EMIC waves Albert (2003)
        ! Sec. 4. omega_m, d_omega normalized to the oxygen gyrofrequency
        ! ----------------------------------------------------------------------
        case('EMIC_HEBAND')
            block
                real*8 :: Omega_O_eq
                Omega_O_eq = omega_e_eq * me / (16.0d0 * mi)
                omega_m = 3.0d0 * Omega_O_eq
                d_omega = 0.5d0 * Omega_O_eq
                omega_lc = omega_m - d_omega
                omega_uc = omega_m + d_omega
            end block

        ! ----------------------------------------------------------------------
        ! EMIC_OBAND: oxygen-band EMIC waves, Albert (2003)
        ! Sec. 4. omega_m, d_omega normalized to the oxygen gyrofrequency.
        ! ----------------------------------------------------------------------
        case('EMIC_OBAND')
            block
                real*8 :: Omega_O_eq
                Omega_O_eq = omega_e_eq * me / (16.0d0 * mi)
                omega_m = 0.9d0 * Omega_O_eq
                d_omega = 0.05d0 * Omega_O_eq
                omega_lc = omega_m - d_omega
                omega_uc = omega_m + d_omega
            end block

        case default
            write(*,*) 'ERROR: Invalid FreqUnit = ', trim(FreqUnit)
            stop

    end select

    return

end subroutine Freq_Spectral