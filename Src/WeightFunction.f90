! ==============================================================================
! SUBROUTINE: WeightFunction
! ==============================================================================
!
! PURPOSE:
!   Computes the integrand of the quasi-linear diffusion coefficients
!   (D_aa, D_ap, D_pp) at a single wave normal angle (theta) and magnetic
!   latitude (lambda), for a given resonance harmonic number (nres).
!
!   This subroutine is called by D_local_n, which integrates the output
!   over theta to produce the local diffusion coefficients. WeightFunction
!   itself returns values that are NOT yet integrated — they are the
!   integrand at one theta point, to be weighted by sin(theta)*dtheta
!   in the calling routine.
!
!   For each resonance harmonic, the resonant wave frequency is found from
!   the resonance condition:
!       omega - k_par * v_par - n * Omega_e / gamma = 0
!   The dispersion properties (refractive index, group velocity, polarisation
!   factor) are obtained via Get_dispersion_properties and
!   Get_polarisation_factor from module dispersion, which currently implements
!   the cold plasma (Stix) formulation and is designed to accept a hot plasma
!   (WHAMP) implementation in future.
!
!   The integrand is assembled from:
!     G1 : frequency-weighted wave power spectral density at resonant omega
!     G2 : wave normal angle distribution (sum of two Gaussians in tan(theta),
!          normalised by int_N_omega)
!     G3 : normalised wave amplitude factor (dBw/B_local)^2 * Omega_e / gamma^2
!     G4 : pitch-angle coupling factor, relating D_aa to D_ap and D_pp
!     Delta_n : resonance broadening term, weighted by Bessel functions
!               J_{n-1}, J_n, J_{n+1} and the wave polarisation factor Fsq
!
!   The three outputs are:
!     WeightValue(1) = G3 * Delta_n * G1 * G2              (D_aa integrand)
!     WeightValue(2) = WeightValue(1) * G4                 (D_ap integrand)
!     WeightValue(3) = WeightValue(2) * G4                 (D_pp integrand)
!
! INPUTS:
!   alpha          - Equatorial electron pitch angle [radians]
!   E              - Electron kinetic energy normalised to rest-mass energy
!                    (E = gamma - 1, dimensionless)
!   L              - Magnetic L-shell [Earth radii, dimensionless]
!   lambda         - Magnetic latitude [radians]
!   Bw             - RMS wave magnetic field amplitude [T]
!   Theta_function - Array(7) defining the wave normal angle distribution:
!                      (1) theta_m1  : mean angle of 1st Gaussian [degrees]
!                      (2) theta_w1  : width of 1st Gaussian [degrees]
!                      (3) theta_m2  : mean angle of 2nd Gaussian [degrees]
!                      (4) theta_w2  : width of 2nd Gaussian [degrees]
!                      (5) theta_min : lower cutoff angle [degrees]
!                      (6) theta_max : upper cutoff angle [degrees]
!                      (7) theta_q   : relative amplitude of 2nd Gaussian
!   nres           - Resonance harmonic number n (integer, e.g. -1, 0, 1, ...)
!   theta          - Wave normal angle at this evaluation point [radians]
!
! OUTPUTS:
!   WeightValue(1) - D_aa integrand at this theta [s^-1 rad^-1]
!   WeightValue(2) - D_ap integrand at this theta [s^-1 rad^-1]
!   WeightValue(3) - D_pp integrand at this theta [s^-1 rad^-1]
!
! CALLS:
!   Local_alpha              - Converts equatorial pitch angle to local pitch angle
!   ElectronDensity          - Returns electron number density at (L, lambda)
!   B_field                  - Returns magnetic field magnitude at (L, lambda)
!   Freq_Spectral            - Returns wave band parameters (omega_m, d_omega, etc.)
!   Roots_finder             - Finds resonant frequencies satisfying resonance condition
!   Get_dispersion_properties - Computes refractive index, group velocity (module dispersion)
!   Get_polarisation_factor  - Computes Fsq (module dispersion)
!   PowerSpectralDensity     - Returns wave PSD at the resonant frequency
!   int_N_omega              - Normalisation integral of wave normal angle distribution
!   Bessel_FirstKind         - Bessel function of the first kind J_n(x)
!
! MODULE DEPENDENCIES:
!   constants  - Physical constants (c, qe, me, mi, PI, ...)
!   variables  - Shared plasma frequency variables written as side effects:
!                omega_e, omega_pe, omega_i, omega_pi
!   dispersion - Dispersion relation interface (Module_Dispersion.f90)
!
! NOTES:
!   - The resonance cone angle is checked; if it falls within the specified
!     wave normal angle range, theta_max is reduced to 0.99 * theta_cone to
!     avoid the singularity in the dispersion relation.
!   - Roots that are zero or lie below the lower-hybrid frequency are silently
!     discarded (WeightValue contributions set to zero).
!   - The sign check on k_vec for nres > 0 will halt execution if violated;
!     this is a consistency check inherited from the original implementation.
! ==============================================================================

subroutine WeightFunction(alpha, E, L, lambda, Bw, Theta_function, nres, theta, WeightValue)
use constants
use variables
use dispersion
implicit none

    real*8, intent(in)  :: alpha, E, L, lambda, Bw, Theta_function(7), theta
    integer, intent(in) :: nres
    real*8, intent(out) :: WeightValue(3)

    ! --- Local diffusion coefficient accumulator ---
    real*8 :: D_tmp(3)

    ! --- Wave normal angle distribution parameters ---
    real*8 :: theta_min, theta_max
    real*8 :: theta_m1, theta_w1, theta_m2, theta_w2, theta_q
    real*8 :: theta_cone, theta_max_eff
    character :: res_cone

    ! --- Particle state ---
    real*8 :: gamma1, betta, alpha_ne
    real*8 :: p, vpar

    ! --- Magnetic field and plasma ---
    real*8 :: lambda_eq
    real*8 :: B_local, B_eq
    real*8 :: omega_e_eq
    real*8 :: Ne_local

    ! --- Wave normal angle in tan form ---
    real*8 :: X, X_m1, X_w1, X_m2, X_w2

    ! --- Resonant frequency roots ---
    real*8 :: roots(4)
    integer :: i

    ! --- Dispersion quantities (from module dispersion) ---
    type(DispersionProperties) :: DP
    real*8 :: n_ref_ind2
    real*8 :: dw_dk, dw_dk_par
    real*8 :: Fsq

    ! --- Resonant frequency and wave vector ---
    real*8 :: omega, k_vec

    ! --- Bessel function values ---
    real*8 :: argJ
    real*8 :: BesselF1, BesselF2, BesselF3

    ! --- Diffusion coefficient factors ---
    real*8 :: Delta_n
    real*8 :: G1, G2, G3, G4

    ! --- Wave amplitude and PSD ---
    real*8 :: dB_local, PSD

    ! --- Wave normal angle normalisation ---
    real*8 :: N_of_omega

    ! -------------------------------------------------------------------------
    ! Initialise outputs
    ! -------------------------------------------------------------------------
    D_tmp(1)    = 0.0
    D_tmp(2)    = 0.0
    D_tmp(3)    = 0.0
    WeightValue(1) = 0.0
    WeightValue(2) = 0.0
    WeightValue(3) = 0.0

    ! -------------------------------------------------------------------------
    ! Unpack wave normal angle distribution parameters
    ! -------------------------------------------------------------------------
    theta_m1  = Theta_function(1)
    theta_w1  = Theta_function(2)
    theta_m2  = Theta_function(3)
    theta_w2  = Theta_function(4)
    theta_min = Theta_function(5)
    theta_max = Theta_function(6)
    theta_q   = Theta_function(7)

    ! -------------------------------------------------------------------------
    ! Local plasma state
    ! -------------------------------------------------------------------------
    call Local_alpha(alpha, lambda, alpha_ne)
    call ElectronDensity(L, lambda, Ne_local)

    X     = tan(theta)
    X_m1  = tan(theta_m1 / 180.0 * PI)
    X_w1  = tan(theta_w1 / 180.0 * PI)
    X_m2  = tan(theta_m2 / 180.0 * PI)
    X_w2  = tan(theta_w2 / 180.0 * PI)

    betta  = sqrt(E*(E + 2.0)) / (E + 1.0)
    gamma1 = E + 1.0

    lambda_eq = 0.0
    call B_field(L, lambda,     B_local)
    call B_field(L, lambda_eq,  B_eq)

    ! Write plasma frequencies to module variables (used by Roots_finder).
    ! int_N_omega now receives these as explicit arguments.
    ! Note: the remaining shared state (for Roots_finder) is a known issue
    ! to be addressed in a future refactoring step.
    omega_e  = qe * B_local / me
    omega_e_eq = qe * B_eq / me
    omega_i  = omega_e * me / mi
    omega_pe = sqrt(Ne_local * qe * qe / (C0 * me))
    omega_pi = omega_pe * sqrt(me / mi)

    call Freq_Spectral(L, lambda, omega_e_eq, omega_m, d_omega, omega_lc, omega_uc)

    ! -------------------------------------------------------------------------
    ! Particle momentum and parallel velocity
    ! -------------------------------------------------------------------------
    p    = me * c * (sqrt((E + 1.0)**2.0 - 1.0))
    vpar = p / gamma1 / me * cos(alpha_ne)

    ! -------------------------------------------------------------------------
    ! Find resonant frequencies
    ! -------------------------------------------------------------------------
    call Roots_finder(alpha, E, L, lambda, nres, theta, roots)

    ! -------------------------------------------------------------------------
    ! Loop over resonant frequency roots
    ! -------------------------------------------------------------------------
    do i = 1, 4

        if (roots(i) == 0.0 .or. roots(i) <= sqrt((omega_pe/omega_e)**2 + 1.0)) then

            ! Root is invalid — contribute nothing
            WeightValue(1) = WeightValue(1)
            WeightValue(2) = WeightValue(2)
            WeightValue(3) = WeightValue(3)

        else

            omega = roots(i)

            ! Wave vector from resonance condition
            k_vec = (omega - (real(nres) * real(sqe) * omega_e) / gamma1) &
                    / (c * betta * cos(alpha_ne) * cos(theta))

            ! Consistency check: k_vec should be positive for nres > 0
            if ((nres > 0) .and. (k_vec < 0)) then
                write(*,*) 'ERROR: negative k_vec for positive resonance harmonic'
                write(*,*) 'N = ', nres, '  Root = ', omega
                stop
            end if

            ! -----------------------------------------------------------------
            ! Dispersion properties via module dispersion
            ! (replaces inline Stix parameter, n^2, dw/dk, and Fsq blocks)
            ! -----------------------------------------------------------------
            call Get_dispersion_properties(omega, k_vec, theta, &
                                           omega_e, omega_pe, omega_i, omega_pi, &
                                           DP)

            n_ref_ind2 = DP%n_ref_ind2
            dw_dk      = DP%dw_dk

            ! -----------------------------------------------------------------
            ! Resonance cone check
            ! -----------------------------------------------------------------
            res_cone     = 'N'
            theta_max_eff = theta_max   ! degrees

            if (DP%P / DP%S < 0) then
                theta_cone = atan(sqrt(-DP%P / DP%S)) * 180.0 / PI   ! degrees

                if (theta_cone*0.9 < theta_max .and. theta_cone*0.99 > theta_max) then
                    res_cone = 'Y'

                else if (theta_cone*0.99 <= theta_max) then
                    theta_max_eff = theta_cone * 0.99
                    res_cone = 'Y'
                end if
            end if

            ! -----------------------------------------------------------------
            ! Bessel functions
            ! -----------------------------------------------------------------
            argJ = abs(k_vec) * c * (sqrt((E+1.0)**2.0 - 1.0)) &
                   * sin(theta) * sin(alpha_ne) / (real(sqe) * omega_e)

            call Bessel_FirstKind(argJ, nres-1, BesselF1)
            call Bessel_FirstKind(argJ, nres,   BesselF2)
            call Bessel_FirstKind(argJ, nres+1, BesselF3)

            ! -----------------------------------------------------------------
            ! Polarisation factor via module dispersion
            ! -----------------------------------------------------------------
            call Get_polarisation_factor(alpha_ne, nres, theta, k_vec, omega, &
                                         BesselF1, BesselF2, BesselF3, &
                                         DP, Fsq)

            ! -----------------------------------------------------------------
            ! Group velocity component along field line
            ! -----------------------------------------------------------------
            dw_dk_par = dw_dk / cos(theta)

            ! -----------------------------------------------------------------
            ! Diffusion coefficient integrand
            ! -----------------------------------------------------------------
            if (vpar == dw_dk_par &
                .or. (theta * 180.0/PI) > theta_max_eff &
                .or. theta_min > theta_max_eff) then

                D_tmp(1) = 0.0
                D_tmp(2) = 0.0
                D_tmp(3) = 0.0

            else

                ! Resonance broadening factor
                Delta_n = PI/2.0 * 1.0/cos(theta) / ((vpar/c)**3.0) * Fsq &
                          * ((real(nres)*real(sqe)*omega_e/(gamma1*omega) &
                             - (sin(alpha_ne))**2.0)**2.0) &
                          / abs(1.0 - dw_dk_par/vpar)

                ! G1: frequency-weighted wave power spectral density
                call PowerSpectralDensity(L, lambda, omega_e_eq, omega, PSD)
                G1 = omega_e * PSD

                ! G2: normalised wave normal angle distribution
                call int_N_omega(alpha, E, L, lambda, Theta_function, &
                                 theta_max_eff, nres, omega, &
                                 omega_e, omega_pe, omega_i, omega_pi, &
                                 res_cone, N_of_omega)
                G2 = (exp(-((X - X_m1)/X_w1)**2.0) &
                     + theta_q * exp(-((X - X_m2)/X_w2)**2.0)) / N_of_omega

                ! G3: normalised wave amplitude factor
                dB_local = Bw
                G3 = omega_e / (gamma1**2.0) * (dB_local / B_local)**2.0

                ! G4: pitch-angle coupling factor (links D_aa to D_ap and D_pp)
                G4 = sin(alpha_ne) * cos(alpha_ne) &
                     / (real(nres)*real(sqe)*omega_e/(gamma1*omega) &
                        - (sin(alpha_ne))**2.0)

                D_tmp(1) = G3 * Delta_n * G1 * G2
                D_tmp(2) = D_tmp(1) * G4
                D_tmp(3) = D_tmp(2) * G4

            end if

            WeightValue(1) = WeightValue(1) + D_tmp(1)
            WeightValue(2) = WeightValue(2) + D_tmp(2)
            WeightValue(3) = WeightValue(3) + D_tmp(3)

        end if

    end do

    return

end subroutine WeightFunction