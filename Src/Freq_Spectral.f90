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
!     2018LBCHORUS, 2018UBCHORUS, 2026LBCHORUS, 2026UBCHORUS), omega_m and
!     d_omega are left undefined on output. For these models, the spectral
!     shape is defined entirely within PowerSpectralDensity via its
!     polynomial fit, and omega_m/d_omega are not needed.
!   - The 2017CHORUS, 2018LBCHORUS, and 2018UBCHORUS cases have identical
!     day/night MLT branches. The MLT branching is retained in case future
!     updates introduce MLT dependence, but currently has no effect.
!   - All frequencies are in rad/s. The HZ case converts from Hz by
!     multiplying by 2*PI. The 2016HISS case uses Hz * 2*PI directly.
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

        ! ----------------------------------------------------------------------
        ! 2018LBCHORUS: lower-band chorus band edges only.
        ! ----------------------------------------------------------------------
        case('2018LBCHORUS')
            omega_lc = 0.05*omega_e_eq
            omega_uc = 0.50*omega_e_eq

        ! ----------------------------------------------------------------------
        ! 2018UBCHORUS: upper-band chorus band edges only.
        ! ----------------------------------------------------------------------
        case('2018UBCHORUS')
            omega_lc = 0.50*omega_e_eq
            omega_uc = 1.00*omega_e_eq

        ! ----------------------------------------------------------------------
        ! 2026LBCHORUS: lower-band chorus band edges only.
        ! ----------------------------------------------------------------------
        case('2026LBCHORUS')
            omega_lc = 0.05*omega_e_eq
            omega_uc = 0.50*omega_e_eq

        ! ----------------------------------------------------------------------
        ! 2026UBCHORUS: upper-band chorus band edges only.
        ! ----------------------------------------------------------------------
        case('2026UBCHORUS')
            omega_lc = 0.50*omega_e_eq
            omega_uc = 0.80*omega_e_eq

        case default
            write(*,*) 'ERROR: Invalid FreqUnit = ', trim(FreqUnit)
            stop

    end select

    return

end subroutine Freq_Spectral