! ==============================================================================
! SUBROUTINE: PowerSpectralDensity
! ==============================================================================
!
! PURPOSE:
!   Returns the normalised wave power spectral density (PSD) at a given
!   wave frequency (Freq), for the spectral model specified by thw frequency
!   distribution in the input file.
!
!   The PSD is normalised so that:
!     integral[ PSD(omega) * domega ] from omega_lc to omega_uc = 1
!
!   It appears as part of G1 in WeightFunction:
!     G1 = omega_e * PSD(omega_resonant)
!
!   Supported models and their parameter dependencies:
!
!     Model          L    lambda   MLT   Kp   Description
!     GAUSSIAN              X            Gaussian spectrum centred on omega_m
!                                       with width d_omega
!     2016HISS       X                  Hiss: L-dependent kappa spectrum
!     2017CHORUS            X      X    Chorus: 7th-degree polynomial in
!                                       Freq/omega_e_eq, MLT and latitude
!                                       dependent
!     2018LBCHORUS          X      X    Lower-band chorus: polynomial,
!                                       MLT and latitude dependent
!     2018UBCHORUS                 X    Upper-band chorus: polynomial,
!                                       MLT dependent only
!     2026LBCHORUS          X      X    Lower-band chorus (2026 update):
!                                       polynomial, MLT and latitude dependent
!     2026UBCHORUS                 X    Upper-band chorus (2026 update):
!                                       polynomial, MLT dependent only
!
! INPUTS:
!   L          - Magnetic L-shell [Earth radii, dimensionless]
!   Lambda     - Magnetic latitude [radians]
!   omega_e_eq - Electron gyrofrequency at the equator [rad/s]
!   Freq       - Wave frequency at which PSD is evaluated [rad/s]
!
! OUTPUTS:
!   PSD        - Normalised power spectral density [s/rad]
!
! MODULE DEPENDENCIES:
!   constants  - FreqDistribution, MLT, RAD, PI
!
! CALLS:
!   errorfunction - Error function erf(x), used in GAUSSIAN normalisation
!
! NOTES:
!   - For the polynomial models (2017CHORUS, 2018LBCHORUS, 2018UBCHORUS,
!     2026LBCHORUS, 2026UBCHORUS), the PSD is a 7th-degree polynomial in
!     Freq/omega_e_eq. The normalisation divides by the definite integral
!     of the polynomial from omega_lc to omega_uc, computed analytically
!     as the antiderivative evaluated at the band edges. The variables
!     norm_lower and norm_upper hold these antiderivative values (renamed
!     from erf1/erf2 in the original, which was misleading since they are
!     not error functions in the polynomial cases).
!   - A PSD floor of 1e-20 is applied in some polynomial models to avoid
!     numerical issues with near-zero or negative PSD values at the band
!     edges. It is not applied in the GAUSSIAN case.
!   - The 2016HISS model uses a kappa distribution rather than a Gaussian
!     or polynomial, with L-dependent parameters.
!   - The 2026UBCHORUS and 2018UBCHORUS models have no latitude dependence
!     (unlike the lower-band models).
!
! ==============================================================================

subroutine PowerSpectralDensity(L, Lambda, omega_e_eq, Freq, &
                                omega_m, omega_lc, omega_uc, d_omega, &
                                PSD)

use constants
implicit none

    real*8, intent(in)  :: L, Lambda, Freq, omega_e_eq
    real*8, intent(in)  :: omega_m, omega_lc, omega_uc, d_omega
    real*8, intent(out) :: PSD

    real*8 :: lambda_deg
    real*8 :: mag1, mag2, erf1, erf2        ! error function args/values (GAUSSIAN only)
    real*8 :: norm_lower, norm_upper         ! polynomial antiderivative at band edges
    real*8 :: c1, c2, kappa, tau, Bw_norm   ! 2016HISS parameters
    real*8 :: beta0, beta1, beta2, beta3, beta4, beta5, beta6, beta7  ! polynomial coefficients
    real*8 :: x_lc, x_uc, x_freq            ! normalised frequency variables

    lambda_deg = Lambda / RAD

    select case(FreqDistribution)

        ! ----------------------------------------------------------------------
        ! GAUSSIAN: Gaussian spectrum centred on omega_m with width d_omega.
        ! Normalised using error functions at the band edges.
        ! ----------------------------------------------------------------------
        case('GAUSSIAN')
            mag1 = (omega_m - omega_lc) / d_omega
            mag2 = (omega_uc - omega_m)  / d_omega

            call errorfunction(mag1, erf1)
            call errorfunction(mag2, erf2)

            PSD = 2.0 * exp(-((Freq - omega_m)/d_omega)**2.0) &
                  / (sqrt(PI) * d_omega * (erf1 + erf2))

        ! ----------------------------------------------------------------------
        ! 2016HISS: L-dependent kappa distribution for hiss waves.
        ! Parameters from Orlova et al. [2016].
        ! ----------------------------------------------------------------------
        case('2016HISS')
            if (L <= 2.5) then
                kappa   = 0.495;  tau = 5.2e4
                c1      = 7.1e-7; c2  = 2.7e-6
                Bw_norm = 9.2446
            else if (L > 2.5 .and. L <= 3.5) then
                kappa   = 0.915;  tau = 1.6e5
                c1      = 4.6e-7; c2  = 5.6e-7
                Bw_norm = 9.5703
            else if (L > 3.5 .and. L <= 4.5) then
                kappa   = 0.624;  tau = 1.6e5
                c1      = 4.8e-7; c2  = 4.5e-7
                Bw_norm = 9.2774
            else
                kappa   = 0.574;  tau = 2e5
                c1      = 5.5e-7; c2  = 3.8e-7
                Bw_norm = 9.6841
            end if

            if (Freq > omega_uc .or. Freq < omega_lc) then
                PSD = 0.0
            else if (Freq <= omega_m .and. Freq >= omega_lc) then
                PSD = c1 * exp(-d_omega/Freq) / (2.0*PI * Bw_norm**2.0) * 1e6
            else
                PSD = c2 * (1.0 + (Freq/(2.0*PI))**2.0 / (kappa*tau))**(-kappa-1.0) &
                      / (2.0*PI * Bw_norm**2.0) * 1e6
            end if

        ! ----------------------------------------------------------------------
        ! 2017CHORUS: 7th-degree polynomial in Freq/omega_e_eq.
        ! MLT-dependent (day/night) and latitude-dependent (< or >= 7 degrees).
        ! Normalised by analytic integral of polynomial over [omega_lc, omega_uc].
        ! ----------------------------------------------------------------------
        case('2017CHORUS')
            if (MLT >= 6.0 .and. MLT < 18.0) then   ! Dayside
                if (lambda_deg <= 7.0) then
                    beta0 =  1.2608e-6;  beta1 = -4.6384e-5
                    beta2 =  6.2732e-4;  beta3 = -3.5303e-3
                    beta4 =  9.3343e-3;  beta5 = -1.0391e-2
                    beta6 =  5.7168e-4;  beta7 =  4.973e-3
                else
                    beta0 =  5.974e-7;   beta1 = -2.2599e-5
                    beta2 =  3.3958e-4;  beta3 = -2.226e-3
                    beta4 =  7.985e-3;   beta5 = -1.6349e-2
                    beta6 =  1.7687e-2;  beta7 = -7.7332e-3
                end if
            else                                      ! Nightside
                if (lambda_deg <= 7.0) then
                    beta0 = -1.2233e-6;  beta1 =  5.1851e-5
                    beta2 = -8.3298e-4;  beta3 =  7.2506e-3
                    beta4 = -3.3065e-2;  beta5 =  8.1568e-2
                    beta6 = -1.0443e-1;  beta7 =  5.4606e-2
                else
                    beta0 =  1.2777e-7;  beta1 = -9.7504e-7
                    beta2 = -2.4561e-5;  beta3 =  6.7819e-4
                    beta4 = -4.5157e-3;  beta5 =  1.3066e-2
                    beta6 = -1.7702e-2;  beta7 =  9.2489e-3
                end if
            end if

            x_freq = Freq / omega_e_eq
            x_lc   = omega_lc / omega_e_eq
            x_uc   = omega_uc / omega_e_eq

            PSD = beta0 + beta1*x_freq + beta2*x_freq**2 + beta3*x_freq**3 &
                + beta4*x_freq**4 + beta5*x_freq**5 + beta6*x_freq**6 &
                + beta7*x_freq**7

            norm_lower = beta0*x_lc     + beta1/2.0*x_lc**2 + beta2/3.0*x_lc**3 &
                       + beta3/4.0*x_lc**4 + beta4/5.0*x_lc**5 + beta5/6.0*x_lc**6 &
                       + beta6/7.0*x_lc**7 + beta7/8.0*x_lc**8

            norm_upper = beta0*x_uc     + beta1/2.0*x_uc**2 + beta2/3.0*x_uc**3 &
                       + beta3/4.0*x_uc**4 + beta4/5.0*x_uc**5 + beta5/6.0*x_uc**6 &
                       + beta6/7.0*x_uc**7 + beta7/8.0*x_uc**8

            PSD = PSD / ((norm_upper - norm_lower) * omega_e_eq)

        ! ----------------------------------------------------------------------
        ! 2018LBCHORUS: lower-band chorus, polynomial PSD.
        ! MLT-dependent (day/night) and latitude-dependent (< or >= 10 degrees).
        ! ----------------------------------------------------------------------
        case('2018LBCHORUS')
            if (MLT >= 6.0 .and. MLT < 18.0) then   ! Dayside
                if (lambda_deg <= 10.0) then
                    beta0 =  3.6057e-7;  beta1 = -1.3014e-5
                    beta2 =  2.1507e-4;  beta3 = -1.4922e-3
                    beta4 =  5.5155e-3;  beta5 = -1.1273e-2
                    beta6 =  1.1739e-2;  beta7 = -4.7206e-3
                else
                    beta0 =  4.4083e-7;  beta1 = -1.4193e-5
                    beta2 =  2.0794e-4;  beta3 = -1.0750e-3
                    beta4 =  2.4121e-3;  beta5 = -1.6367e-3
                    beta6 = -2.1965e-3;  beta7 =  3.0513e-3
                end if
            else                                      ! Nightside
                if (lambda_deg <= 10.0) then
                    beta0 = -5.3446e-7;  beta1 =  2.4890e-5
                    beta2 = -4.2571e-4;  beta3 =  4.1014e-3
                    beta4 = -2.0035e-2;  beta5 =  5.2005e-2
                    beta6 = -6.9501e-2;  beta7 =  3.7751e-2
                else
                    beta0 = -1.2997e-7;  beta1 =  9.0475e-6
                    beta2 = -1.8619e-4;  beta3 =  1.9125e-3
                    beta4 = -9.8212e-3;  beta5 =  2.6168e-2
                    beta6 = -3.4857e-2;  beta7 =  1.8411e-2
                end if
            end if

            x_freq = Freq / omega_e_eq
            x_lc   = omega_lc / omega_e_eq
            x_uc   = omega_uc / omega_e_eq

            PSD = beta0 + beta1*x_freq + beta2*x_freq**2 + beta3*x_freq**3 &
                + beta4*x_freq**4 + beta5*x_freq**5 + beta6*x_freq**6 &
                + beta7*x_freq**7

            norm_lower = beta0*x_lc     + beta1/2.0*x_lc**2 + beta2/3.0*x_lc**3 &
                       + beta3/4.0*x_lc**4 + beta4/5.0*x_lc**5 + beta5/6.0*x_lc**6 &
                       + beta6/7.0*x_lc**7 + beta7/8.0*x_lc**8

            norm_upper = beta0*x_uc     + beta1/2.0*x_uc**2 + beta2/3.0*x_uc**3 &
                       + beta3/4.0*x_uc**4 + beta4/5.0*x_uc**5 + beta5/6.0*x_uc**6 &
                       + beta6/7.0*x_uc**7 + beta7/8.0*x_uc**8

            PSD = PSD / ((norm_upper - norm_lower) * omega_e_eq)
            if (PSD < 1.0e-20) PSD = 1.0e-20

        ! ----------------------------------------------------------------------
        ! 2018UBCHORUS: upper-band chorus, polynomial PSD.
        ! MLT-dependent (day/night) only, no latitude dependence.
        ! ----------------------------------------------------------------------
        case('2018UBCHORUS')
            if (MLT >= 6.0 .and. MLT < 18.0) then   ! Dayside
                beta0 =  2.0543522894746959e-4;  beta1 = -2.1251727655271302e-3
                beta2 =  9.2309482352891624e-3;  beta3 = -2.1851729739887232e-2
                beta4 =  3.0496655461723709e-2;  beta5 = -2.5133977878939284e-2
                beta6 =  1.1343438413049637e-2;  beta7 = -2.1655977959280860e-3
            else                                      ! Nightside
                beta0 = -1.4010228888340406e-4;  beta1 =  8.2669571248420001e-4
                beta2 = -1.3655153785780715e-3;  beta3 = -1.1351830489998828e-3
                beta4 =  6.6938370901288491e-3;  beta5 = -9.0877873665073227e-3
                beta6 =  5.4834375943244367e-3;  beta7 = -1.2753845539271731e-3
            end if

            x_freq = Freq / omega_e_eq
            x_lc   = omega_lc / omega_e_eq
            x_uc   = omega_uc / omega_e_eq

            PSD = beta0 + beta1*x_freq + beta2*x_freq**2 + beta3*x_freq**3 &
                + beta4*x_freq**4 + beta5*x_freq**5 + beta6*x_freq**6 &
                + beta7*x_freq**7

            norm_lower = beta0*x_lc     + beta1/2.0*x_lc**2 + beta2/3.0*x_lc**3 &
                       + beta3/4.0*x_lc**4 + beta4/5.0*x_lc**5 + beta5/6.0*x_lc**6 &
                       + beta6/7.0*x_lc**7 + beta7/8.0*x_lc**8

            norm_upper = beta0*x_uc     + beta1/2.0*x_uc**2 + beta2/3.0*x_uc**3 &
                       + beta3/4.0*x_uc**4 + beta4/5.0*x_uc**5 + beta5/6.0*x_uc**6 &
                       + beta6/7.0*x_uc**7 + beta7/8.0*x_uc**8

            PSD = PSD / ((norm_upper - norm_lower) * omega_e_eq)
            if (PSD < 1.0e-20) PSD = 1.0e-20

        ! ----------------------------------------------------------------------
        ! 2026LBCHORUS: lower-band chorus (2026 update), polynomial PSD.
        ! MLT-dependent (day/night) and latitude-dependent (low/mid/high).
        ! Latitude bands: <= 10, 10-20, > 20 degrees.
        ! ----------------------------------------------------------------------
        case('2026LBCHORUS')
            if (MLT >= 6.0 .and. MLT < 18.0) then   ! Dayside
                if (lambda_deg <= 10.0) then
                    beta0 =  1.11336540917184e-7;   beta1 = -3.72912612422969e-7
                    beta2 = -9.16958991615852e-6;   beta3 =  2.15723693886208e-4
                    beta4 = -1.26806058032446e-3;   beta5 =  3.59190154486275e-3
                    beta6 = -5.35827750817583e-3;   beta7 =  3.35331572494158e-3
                else if (lambda_deg > 10.0 .and. lambda_deg <= 20.0) then
                    beta0 =  4.30576556905264e-7;   beta1 = -1.93877149456913e-5
                    beta2 =  3.5456331476089e-4;    beta3 = -2.47357343147156e-3
                    beta4 =  9.07129101805928e-3;   beta5 = -1.89341440278942e-2
                    beta6 =  2.11954576442875e-2;   beta7 = -9.80281543417116e-3
                else
                    beta0 =  1.01790061236938e-6;   beta1 = -4.82797987996845e-5
                    beta2 =  9.29651327598179e-4;   beta3 = -7.71658195124899e-3
                    beta4 =  3.38877002399438e-2;   beta5 = -8.25886472235263e-2
                    beta6 =  1.05301327489477e-1;   beta7 = -5.46805009108624e-2
                end if
            else                                      ! Nightside
                if (lambda_deg <= 10.0) then
                    beta0 = -2.97874170950016e-7;   beta1 =  1.4322737443405e-5
                    beta2 = -2.55070014505957e-4;   beta3 =  2.66338023050004e-3
                    beta4 = -1.44207301553904e-2;   beta5 =  4.18529434765739e-2
                    beta6 = -6.20248242140477e-2;   beta7 =  3.66919071518278e-2
                else if (lambda_deg > 10.0 .and. lambda_deg <= 20.0) then
                    beta0 =  1.87395736698176e-7;   beta1 = -7.4583549236649e-6
                    beta2 =  1.15119097024855e-4;   beta3 = -7.12856986278894e-4
                    beta4 =  2.19744709855392e-3;   beta5 = -3.58152843005299e-3
                    beta6 =  2.90037927467941e-3;   beta7 = -8.78285849663509e-4
                else
                    beta0 =  1.41774960798836e-7;   beta1 = -4.18554459149626e-6
                    beta2 =  9.27712392661608e-5;   beta3 = -8.44912962502075e-4
                    beta4 =  3.8042235566409e-3;    beta5 = -9.12206447572251e-3
                    beta6 =  1.12105133295341e-2;   beta7 = -5.566214822435e-3
                end if
            end if

            x_freq = Freq / omega_e_eq
            x_lc   = omega_lc / omega_e_eq
            x_uc   = omega_uc / omega_e_eq

            PSD = beta0 + beta1*x_freq + beta2*x_freq**2 + beta3*x_freq**3 &
                + beta4*x_freq**4 + beta5*x_freq**5 + beta6*x_freq**6 &
                + beta7*x_freq**7

            norm_lower = beta0*x_lc     + beta1/2.0*x_lc**2 + beta2/3.0*x_lc**3 &
                       + beta3/4.0*x_lc**4 + beta4/5.0*x_lc**5 + beta5/6.0*x_lc**6 &
                       + beta6/7.0*x_lc**7 + beta7/8.0*x_lc**8

            norm_upper = beta0*x_uc     + beta1/2.0*x_uc**2 + beta2/3.0*x_uc**3 &
                       + beta3/4.0*x_uc**4 + beta4/5.0*x_uc**5 + beta5/6.0*x_uc**6 &
                       + beta6/7.0*x_uc**7 + beta7/8.0*x_uc**8

            PSD = PSD / ((norm_upper - norm_lower) * omega_e_eq)
            if (PSD < 1.0e-20) PSD = 1.0e-20

        ! ----------------------------------------------------------------------
        ! 2026UBCHORUS: upper-band chorus (2026 update), polynomial PSD.
        ! MLT-dependent (day/night) only, no latitude dependence.
        ! ----------------------------------------------------------------------
        case('2026UBCHORUS')
            if (MLT >= 6.0 .and. MLT < 18.0) then   ! Dayside
                beta0 =  1.4947221225619e-3;    beta1 = -1.5803770764091e-2
                beta2 =  7.11096808978984e-2;   beta3 = -1.76587151860897e-1
                beta4 =  2.61515137534286e-1;   beta5 = -2.31080159701795e-1
                beta6 =  1.12859681459824e-1;   beta7 = -2.35129368141491e-2
            else                                      ! Nightside
                beta0 =  1.61805052165801e-2;   beta1 = -1.77099131141267e-1
                beta2 =  8.25782709253675e-1;   beta3 = -2.12703533294286e0
                beta4 =  3.26969853152609e0;    beta5 = -3.0005692071183e0
                beta6 =  1.52254312687049e0;    beta7 = -3.29626122022566e-1
            end if

            x_freq = Freq / omega_e_eq
            x_lc   = omega_lc / omega_e_eq
            x_uc   = omega_uc / omega_e_eq

            PSD = beta0 + beta1*x_freq + beta2*x_freq**2 + beta3*x_freq**3 &
                + beta4*x_freq**4 + beta5*x_freq**5 + beta6*x_freq**6 &
                + beta7*x_freq**7

            norm_lower = beta0*x_lc     + beta1/2.0*x_lc**2 + beta2/3.0*x_lc**3 &
                       + beta3/4.0*x_lc**4 + beta4/5.0*x_lc**5 + beta5/6.0*x_lc**6 &
                       + beta6/7.0*x_lc**7 + beta7/8.0*x_lc**8

            norm_upper = beta0*x_uc     + beta1/2.0*x_uc**2 + beta2/3.0*x_uc**3 &
                       + beta3/4.0*x_uc**4 + beta4/5.0*x_uc**5 + beta5/6.0*x_uc**6 &
                       + beta6/7.0*x_uc**7 + beta7/8.0*x_uc**8

            PSD = PSD / ((norm_upper - norm_lower) * omega_e_eq)
            if (PSD < 1.0e-20) PSD = 1.0e-20

        case default
            write(*,*) 'ERROR: Invalid FreqDistribution = ', trim(FreqDistribution)
            stop

    end select

    return

end subroutine PowerSpectralDensity