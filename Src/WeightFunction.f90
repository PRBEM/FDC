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
!   itself returns values that are NOT yet integrated. They are the
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
!                              (DiffusionFormulation = 'LYONS' only)
!   cunningham_jacobian      - Local wave-normal-angle normalisation factor
!                              (DiffusionFormulation = 'CUNNINGHAM' only;
!                              module dispersion)
!   Bessel_FirstKind         - Bessel function of the first kind J_n(x)
!
! MODULE DEPENDENCIES:
!   constants  - Physical constants (c, qe, me, mi, PI, ...)
!   dispersion - Dispersion relation interface (Module_Dispersion.f90)
!
! NOTES:
!   - The resonance cone angle is checked; if it falls within the specified
!     wave normal angle range, theta_max is reduced to 0.99 * theta_cone to
!     avoid the singularity in the dispersion relation. This guard is
!     primarily needed by the 'LYONS' formulation's int_N_omega integral,
!     which (per Cunningham 2023) is not integrable if the resonance cone
!     angle is included as an integration endpoint and g(X) > 0 there. The
!     'CUNNINGHAM' formulation's local Jacobian factor does not have this
!     non-integrability problem (no integral over theta is performed), but
!     the guard is left in place for both formulations for now since it
!     also affects the validity of n_ref_ind2 itself near the cone.
!   - Roots that are zero or lie below the lower-hybrid frequency are silently
!     discarded (WeightValue contributions set to zero).
!   - The sign check on k_vec for nres > 0 will halt execution if violated;
!     this is a consistency check inherited from the original implementation.
!   - DiffusionFormulation ('LYONS' or 'CUNNINGHAM', from module constants)
!     selects which normalisation is used for G2, as described above. The
!     'LYONS' path reproduces the code's original behaviour exactly. 
!     See Module_Dispersion.f90 for further documentation of the two formulations.
! ==============================================================================

subroutine WeightFunction(alpha, E, L, lambda, Bw, Theta_function, nres, theta, WeightValue, &
                          whamp_curve)

use constants
use dispersion
use accounting
implicit none
    include "mpif.h"

    real*8, intent(in)  :: alpha, E, L, lambda, Bw, Theta_function(7), theta
    integer, intent(in) :: nres
    real*8, intent(out) :: WeightValue(3)
    type(WHAMPCurve), intent(in) :: whamp_curve

    ! --- Local diffusion coefficient accumulator ---
    real*8 :: D_tmp(3)

    ! --- Wave normal angle distribution parameters ---
    real*8 :: theta_min, theta_max
    real*8 :: theta_m1, theta_w1, theta_m2, theta_w2, theta_q
    real*8 :: theta_cone, theta_max_eff
    integer          :: acct_status
    real*8           :: acct_N
    character :: res_cone
    character(len=8) :: whamp_path

    ! --- Particle state ---
    real*8 :: gamma1, betta, alpha_ne
    real*8 :: p, vpar

    ! --- Magnetic field and plasma ---
    real*8 :: lambda_eq
    real*8 :: B_local, B_eq
    real*8 :: omega_e_eq
    real*8 :: Ne_local

    ! --- Local plasma frequencies ---
    real*8 :: omega_e, omega_pe, omega_i, omega_pi
    real*8 :: omega_lc, omega_uc, omega_m, d_omega

    ! --- Wave normal angle in tan form ---
    real*8 :: X, X_m1, X_w1, X_m2, X_w2

    ! --- Resonant frequency roots ---
    real*8 :: roots(4)
    real*8 :: k_sign(4)   ! unused here — WeightFunction derives k_vec's sign
                          ! directly from the resonance condition below
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

    integer(8) :: t_prof_1, t_prof_2, prof_rate

    integer :: i_theta_lo_cached, i_theta_hi_cached
    real*8  :: w_theta_cached

    integer, save :: my_rank_cached = -1
    integer :: mpi_ierr

#ifdef PROFILE
    call system_clock(t_prof_1)
#endif

    ! -------------------------------------------------------------------------
    ! Initialise outputs
    ! -------------------------------------------------------------------------
    if (my_rank_cached < 0) then
        call MPI_Comm_rank(MPI_COMM_WORLD, my_rank_cached, mpi_ierr)
    end if
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

    ! Compute local plasma frequencies
    omega_e  = qe * B_local / me
    omega_e_eq = qe * B_eq / me
    omega_i  = omega_e * me / mi
    omega_pe = sqrt(Ne_local * qe * qe / (C0 * me))
    omega_pi = omega_pe * sqrt(me / mi)

    ! Initialise frequency parameters to zero before Freq_Spectral call.
    ! Some FreqUnit models (e.g. 2018UBCHORUS, 2017CHORUS) only set
    ! omega_lc and omega_uc, leaving omega_m and d_omega unchanged.
    ! Explicit initialisation here ensures no garbage values are used.
    omega_m   = 0.0d0
    d_omega   = 0.0d0
    omega_lc  = 0.0d0
    omega_uc  = 0.0d0

    call Freq_Spectral(L, lambda, omega_e_eq, omega_m, d_omega, omega_lc, omega_uc)

    ! -------------------------------------------------------------------------
    ! Particle momentum and parallel velocity
    ! -------------------------------------------------------------------------
    p    = me * c * (sqrt((E + 1.0)**2.0 - 1.0))
    vpar = p / gamma1 / me * cos(alpha_ne)

    ! -------------------------------------------------------------------------
    ! Find resonant frequencies
    ! -------------------------------------------------------------------------
    call Roots_finder(alpha, E, L, lambda, nres, theta, &
                        omega_e, omega_pe, omega_i, omega_pi, &
                        omega_lc, omega_uc, &
                        roots, whamp_curve, k_sign)

    ! -------------------------------------------------------------------------
    ! Loop over resonant frequency roots
    ! -------------------------------------------------------------------------
    do i = 1, 4

        if (roots(i) == 0.0 .or. roots(i) <= sqrt((omega_pe/omega_e)**2 + 1.0)) then

            !call acct_write(DispersionType, my_rank_cached, 1, i, nres, &
            !                L, E, alpha*180.0d0/PI, lambda*180.0d0/PI, &
            !                theta*180.0d0/PI, -1.0d0, &
            !                roots(i), -1.0d0, -1.0d0, -1.0d0, -1.0d0, 0.0d0)

        else

            omega = roots(i)
            ! Wave vector from resonance condition
            k_vec = (omega - (real(nres) * real(sqe) * omega_e) / gamma1) &
                    / (c * betta * cos(alpha_ne) * cos(theta))

            ! Consistency check: k_vec should be positive for nres > 0
            !if ((nres > 0) .and. (k_vec < 0)) then
            !    write(*,*) 'ERROR: negative k_vec for positive resonance harmonic'
            !    write(*,*) 'N = ', nres, '  Root = ', omega
            !    stop
            !end if

            ! -----------------------------------------------------------------
            ! Dispersion properties via module dispersion
            ! (replaces inline Stix parameter, n^2, dw/dk, and Fsq blocks)
            ! -----------------------------------------------------------------
            call Get_dispersion_properties(omega, k_vec, theta, &
                                           omega_e, omega_pe, omega_i, omega_pi, &
                                           DP)

            n_ref_ind2 = DP%n_ref_ind2
            dw_dk      = DP%dw_dk
            DP%EFL(1) = (0.0d0, 0.0d0)
            DP%EFL(2) = (0.0d0, 0.0d0)
            DP%EFL(3) = (0.0d0, 0.0d0)
            DP%BFL(1) = (0.0d0, 0.0d0)
            DP%BFL(2) = (0.0d0, 0.0d0)
            DP%BFL(3) = (0.0d0, 0.0d0)

            whamp_path = 'NONE'

            ! -----------------------------------------------------------------
            ! HOT case: overwrite n_ref_ind2 and dw_dk from WHAMP curve
            ! -----------------------------------------------------------------
            if ((DispersionType == 'HOT' .or. DispersionType == 'EMIC') .and. whamp_curve%n_k > 0) then
                block
                    integer :: i_k, i_theta_lo, i_theta_hi
                    real*8  :: w_theta, omega_lo_t, omega_hi_t
                    real*8  :: frac, n2_lo, n2_hi, dw_lo, dw_hi
                    logical :: valid_lo, valid_hi

                    ! Find bracketing theta indices and interpolation weight
                    i_theta_lo = 1
                    i_theta_hi = 1
                    w_theta    = 0.0d0
                    if (whamp_curve%n_theta > 1) then
                        if (theta <= whamp_curve%theta(1)) then
                            i_theta_lo = 1
                            i_theta_hi = 1
                        else if (theta >= whamp_curve%theta(whamp_curve%n_theta)) then
                            i_theta_lo = whamp_curve%n_theta
                            i_theta_hi = whamp_curve%n_theta
                        else
                            do i_k = 2, whamp_curve%n_theta
                                if (whamp_curve%theta(i_k-1) <= theta .and. &
                                    whamp_curve%theta(i_k)   >= theta) then
                                    i_theta_lo = i_k - 1
                                    i_theta_hi = i_k
                                    w_theta = (theta - whamp_curve%theta(i_theta_lo)) / &
                                              (whamp_curve%theta(i_theta_hi) - &
                                               whamp_curve%theta(i_theta_lo))
                                    exit
                                end if
                            end do
                        end if
                    end if

                    ! Find bracketing k points by interpolated omega
                    block
                        logical :: found
                        found = .false.
                        do i_k = 2, whamp_curve%n_k
                            valid_lo = whamp_curve%valid(i_theta_lo, i_k) .and. &
                                    whamp_curve%valid(i_theta_lo, i_k-1)
                            valid_hi = whamp_curve%valid(i_theta_hi, i_k) .and. &
                                    whamp_curve%valid(i_theta_hi, i_k-1)
                            if (.not. valid_lo .or. .not. valid_hi) cycle

                            omega_lo_t = (1.0d0-w_theta)*whamp_curve%omega(i_theta_lo,i_k-1) + &
                                        w_theta*whamp_curve%omega(i_theta_hi,i_k-1)
                            omega_hi_t = (1.0d0-w_theta)*whamp_curve%omega(i_theta_lo,i_k)   + &
                                        w_theta*whamp_curve%omega(i_theta_hi,i_k)

                            !if (E > 0.18d0 .and. E < 0.2d0 .and. abs(alpha_ne*180.0d0/PI - 56.67d0) < 1.0d0 &
                            !    .and. nres == 1 .and. omega_hi_t > 20000.0d0 .and. omega_lo_t < 25000.0d0) then
                            !    open(86, file='bracket_a0005_near22k.txt', status='unknown', position='append')
                            !    write(86,*) 'omega_target=', omega, ' i_k=', i_k, &
                            !                ' omega_lo_t=', omega_lo_t, ' omega_hi_t=', omega_hi_t, &
                            !                ' valid_lo=', valid_lo, ' valid_hi=', valid_hi
                            !    close(86)
                            !end if

                            if (omega_lo_t <= omega .and. omega_hi_t >= omega) then
                                frac  = (omega - omega_lo_t) / (omega_hi_t - omega_lo_t)
                                n2_lo = (whamp_curve%k(i_k-1) * c / omega_lo_t)**2.0d0
                                n2_hi = (whamp_curve%k(i_k)   * c / omega_hi_t)**2.0d0
                                dw_lo = (1.0d0-w_theta)*whamp_curve%dw_dk(i_theta_lo,i_k-1) + &
                                        w_theta*whamp_curve%dw_dk(i_theta_hi,i_k-1)
                                dw_hi = (1.0d0-w_theta)*whamp_curve%dw_dk(i_theta_lo,i_k)   + &
                                        w_theta*whamp_curve%dw_dk(i_theta_hi,i_k)
                                n_ref_ind2 = n2_lo + frac * (n2_hi - n2_lo)
                                dw_dk      = dw_lo + frac * (dw_hi - dw_lo)
                                if (k_vec < 0.0d0) then
                                    dw_dk = -dw_dk
                                end if
                                if (abs(dw_dk) < 1.0d-6 * c) then
                                    n_ref_ind2 = 0.0d0  
                                end if
                                DP%EFL(1) = (1.0d0-w_theta)*whamp_curve%EFL(1,i_theta_lo,i_k-1) + &
                                            w_theta*whamp_curve%EFL(1,i_theta_hi,i_k-1) + &
                                            frac*( ((1.0d0-w_theta)*whamp_curve%EFL(1,i_theta_lo,i_k) + &
                                                    w_theta*whamp_curve%EFL(1,i_theta_hi,i_k)) - &
                                                   ((1.0d0-w_theta)*whamp_curve%EFL(1,i_theta_lo,i_k-1) + &
                                                    w_theta*whamp_curve%EFL(1,i_theta_hi,i_k-1)) )
                                DP%EFL(2) = (1.0d0-w_theta)*whamp_curve%EFL(2,i_theta_lo,i_k-1) + &
                                            w_theta*whamp_curve%EFL(2,i_theta_hi,i_k-1) + &
                                            frac*( ((1.0d0-w_theta)*whamp_curve%EFL(2,i_theta_lo,i_k) + &
                                                    w_theta*whamp_curve%EFL(2,i_theta_hi,i_k)) - &
                                                   ((1.0d0-w_theta)*whamp_curve%EFL(2,i_theta_lo,i_k-1) + &
                                                    w_theta*whamp_curve%EFL(2,i_theta_hi,i_k-1)) )
                                DP%EFL(3) = (1.0d0-w_theta)*whamp_curve%EFL(3,i_theta_lo,i_k-1) + &
                                            w_theta*whamp_curve%EFL(3,i_theta_hi,i_k-1) + &
                                            frac*( ((1.0d0-w_theta)*whamp_curve%EFL(3,i_theta_lo,i_k) + &
                                                    w_theta*whamp_curve%EFL(3,i_theta_hi,i_k)) - &
                                                   ((1.0d0-w_theta)*whamp_curve%EFL(3,i_theta_lo,i_k-1) + &
                                                    w_theta*whamp_curve%EFL(3,i_theta_hi,i_k-1)) )
                                DP%BFL(1) = (1.0d0-w_theta)*whamp_curve%BFL(1,i_theta_lo,i_k-1) + &
                                            w_theta*whamp_curve%BFL(1,i_theta_hi,i_k-1) + &
                                            frac*( ((1.0d0-w_theta)*whamp_curve%BFL(1,i_theta_lo,i_k) + &
                                                    w_theta*whamp_curve%BFL(1,i_theta_hi,i_k)) - &
                                                   ((1.0d0-w_theta)*whamp_curve%BFL(1,i_theta_lo,i_k-1) + &
                                                    w_theta*whamp_curve%BFL(1,i_theta_hi,i_k-1)) )
                                DP%BFL(2) = (1.0d0-w_theta)*whamp_curve%BFL(2,i_theta_lo,i_k-1) + &
                                            w_theta*whamp_curve%BFL(2,i_theta_hi,i_k-1) + &
                                            frac*( ((1.0d0-w_theta)*whamp_curve%BFL(2,i_theta_lo,i_k) + &
                                                    w_theta*whamp_curve%BFL(2,i_theta_hi,i_k)) - &
                                                   ((1.0d0-w_theta)*whamp_curve%BFL(2,i_theta_lo,i_k-1) + &
                                                    w_theta*whamp_curve%BFL(2,i_theta_hi,i_k-1)) )
                                DP%BFL(3) = (1.0d0-w_theta)*whamp_curve%BFL(3,i_theta_lo,i_k-1) + &
                                            w_theta*whamp_curve%BFL(3,i_theta_hi,i_k-1) + &
                                            frac*( ((1.0d0-w_theta)*whamp_curve%BFL(3,i_theta_lo,i_k) + &
                                                    w_theta*whamp_curve%BFL(3,i_theta_hi,i_k)) - &
                                                   ((1.0d0-w_theta)*whamp_curve%BFL(3,i_theta_lo,i_k-1) + &
                                                    w_theta*whamp_curve%BFL(3,i_theta_hi,i_k-1)) )
                                found = .true.
                                whamp_path = 'INTERP'
                                exit
                            end if
                        end do

                        ! Tolerance fallback: if omega is just barely outside
                        ! the curve's valid range at this theta (a precision
                        ! edge case near the whistler cutoff), use the
                        ! nearest valid endpoint's values instead of leaving
                        ! n_ref_ind2/dw_dk at their stale HOT-placeholder
                        ! zero values.
                        if (.not. found) then
                            block
                                integer :: i_first_valid, i_last_valid, i_s2
                                real*8  :: omega_min_t, omega_max_t, tol
                                logical :: vlo2, vhi2

                                i_first_valid = 0
                                i_last_valid  = 0
                                do i_s2 = 1, whamp_curve%n_k
                                    vlo2 = whamp_curve%valid(i_theta_lo, i_s2)
                                    vhi2 = whamp_curve%valid(i_theta_hi, i_s2)
                                    if (vlo2 .and. vhi2) then
                                        if (i_first_valid == 0) i_first_valid = i_s2
                                        i_last_valid = i_s2
                                    end if
                                end do

                                if (i_first_valid > 0) then
                                    omega_min_t = (1.0d0-w_theta)*whamp_curve%omega(i_theta_lo,i_first_valid) + &
                                                  w_theta*whamp_curve%omega(i_theta_hi,i_first_valid)
                                    omega_max_t = (1.0d0-w_theta)*whamp_curve%omega(i_theta_lo,i_last_valid) + &
                                                  w_theta*whamp_curve%omega(i_theta_hi,i_last_valid)
                                    tol = 0.01d0

                                    if (omega > omega_max_t .and. abs(omega-omega_max_t) < tol*omega_max_t) then
                                        n_ref_ind2 = (whamp_curve%k(i_last_valid)*c/omega_max_t)**2.0d0
                                        dw_dk = (1.0d0-w_theta)*whamp_curve%dw_dk(i_theta_lo,i_last_valid) + &
                                                w_theta*whamp_curve%dw_dk(i_theta_hi,i_last_valid)
                                        if (k_vec < 0.0d0) dw_dk = -dw_dk
                                        if (abs(dw_dk) < 1.0d-6 * c) found = .false.
                                        DP%EFL(1) = (1.0d0-w_theta)*whamp_curve%EFL(1,i_theta_lo,i_last_valid) + &
                                                    w_theta*whamp_curve%EFL(1,i_theta_hi,i_last_valid)
                                        DP%EFL(2) = (1.0d0-w_theta)*whamp_curve%EFL(2,i_theta_lo,i_last_valid) + &
                                                    w_theta*whamp_curve%EFL(2,i_theta_hi,i_last_valid)
                                        DP%EFL(3) = (1.0d0-w_theta)*whamp_curve%EFL(3,i_theta_lo,i_last_valid) + &
                                                    w_theta*whamp_curve%EFL(3,i_theta_hi,i_last_valid)
                                        DP%BFL(1) = (1.0d0-w_theta)*whamp_curve%BFL(1,i_theta_lo,i_last_valid) + &
                                                    w_theta*whamp_curve%BFL(1,i_theta_hi,i_last_valid)
                                        DP%BFL(2) = (1.0d0-w_theta)*whamp_curve%BFL(2,i_theta_lo,i_last_valid) + &
                                                    w_theta*whamp_curve%BFL(2,i_theta_hi,i_last_valid)
                                        DP%BFL(3) = (1.0d0-w_theta)*whamp_curve%BFL(3,i_theta_lo,i_last_valid) + &
                                                    w_theta*whamp_curve%BFL(3,i_theta_hi,i_last_valid)
                                        whamp_path = 'TOLFB'
                                    else if (omega < omega_min_t .and. abs(omega-omega_min_t) < tol*omega_min_t) then
                                        n_ref_ind2 = (whamp_curve%k(i_first_valid)*c/omega_min_t)**2.0d0
                                        dw_dk = (1.0d0-w_theta)*whamp_curve%dw_dk(i_theta_lo,i_first_valid) + &
                                                w_theta*whamp_curve%dw_dk(i_theta_hi,i_first_valid)
                                        if (k_vec < 0.0d0) dw_dk = -dw_dk
                                        if (abs(dw_dk) < 1.0d-6 * c) found = .false.
                                        DP%EFL(1) = (1.0d0-w_theta)*whamp_curve%EFL(1,i_theta_lo,i_first_valid) + &
                                                    w_theta*whamp_curve%EFL(1,i_theta_hi,i_first_valid)
                                        DP%EFL(2) = (1.0d0-w_theta)*whamp_curve%EFL(2,i_theta_lo,i_first_valid) + &
                                                    w_theta*whamp_curve%EFL(2,i_theta_hi,i_first_valid)
                                        DP%EFL(3) = (1.0d0-w_theta)*whamp_curve%EFL(3,i_theta_lo,i_first_valid) + &
                                                    w_theta*whamp_curve%EFL(3,i_theta_hi,i_first_valid)
                                        DP%BFL(1) = (1.0d0-w_theta)*whamp_curve%BFL(1,i_theta_lo,i_first_valid) + &
                                                    w_theta*whamp_curve%BFL(1,i_theta_hi,i_first_valid)
                                        DP%BFL(2) = (1.0d0-w_theta)*whamp_curve%BFL(2,i_theta_lo,i_first_valid) + &
                                                    w_theta*whamp_curve%BFL(2,i_theta_hi,i_first_valid)
                                        DP%BFL(3) = (1.0d0-w_theta)*whamp_curve%BFL(3,i_theta_lo,i_first_valid) + &
                                                    w_theta*whamp_curve%BFL(3,i_theta_hi,i_first_valid)
                                        whamp_path = 'TOLFB'
                                    end if
                                end if
                            end block
                        end if
                        if (.not. found) then
                            block
                                real*8 :: n2_targeted, dw_targeted, omega_found
                                complex(kind=8) :: EFL_targeted(3)
                                complex(kind=8) :: BFL_targeted(3)
                                logical :: valid_targeted
                                call get_whamp_point(whamp_curve%fce_khz, whamp_curve%fpe_khz, theta, &
                                                    abs(k_vec), omega, &
                                                    n2_targeted, dw_targeted, EFL_targeted, BFL_targeted, valid_targeted)

                                !if (E > 0.055d0 .and. E < 0.062d0 .and. abs(alpha_ne*180.0d0/PI - 10.0d0) < 1.0d0 &
                                !    .and. nres == 1) then
                                !    open(81, file='rescue_debug.txt', status='unknown', position='append')
                                !    write(81,*) 'rescue attempted: valid_targeted=', valid_targeted, &
                                !                ' k_vec=', k_vec, ' omega=', omega, &
                                !                ' fce=', whamp_curve%fce_khz, ' fpe=', whamp_curve%fpe_khz
                                !    close(81)
                                !end if

                                if (valid_targeted .and. n2_targeted > 0.0d0) then
                                    omega_found = abs(k_vec) * c / sqrt(n2_targeted)
                                else
                                    omega_found = -1.0d0
                                end if

                                if (valid_targeted .and. n2_targeted > 0.0d0 .and. &
                                    abs(omega_found - omega) <= 1.0d-3 * omega) then
                                    n_ref_ind2 = n2_targeted
                                    dw_dk      = dw_targeted
                                    if (k_vec < 0.0d0) dw_dk = -dw_dk
                                    DP%EFL(1) = EFL_targeted(1)
                                    DP%EFL(2) = EFL_targeted(2)
                                    DP%EFL(3) = EFL_targeted(3)
                                    DP%BFL(1) = BFL_targeted(1)
                                    DP%BFL(2) = BFL_targeted(2)
                                    DP%BFL(3) = BFL_targeted(3)
                                    found = .true.
                                    whamp_path = 'RESCUE'
                                end if
                            end block
                        end if
                        !if (E > 0.18d0 .and. E < 0.2d0 .and. abs(alpha_ne*180.0d0/PI - 56.67d0) < 1.0d0 &
                        !    .and. nres == 1) then
                        !    open(87, file='interp_a0005_debug.txt', status='unknown', position='append')
                        !    write(87,*) 'omega=', omega, ' found=', found, &
                        !                ' n_ref_ind2=', n_ref_ind2, ' dw_dk=', dw_dk
                        !    close(87)
                        !end if
                        if (.not. found) then
                            !call acct_write(DispersionType, my_rank_cached, 2, i, nres, &
                            !                L, E, alpha*180.0d0/PI, lambda*180.0d0/PI, &
                            !                theta*180.0d0/PI, -1.0d0, &
                            !                omega, -1.0d0, -1.0d0, -1.0d0, -1.0d0, 0.0d0)
                            cycle
                        end if
                    end block
                    block
                    real*8 :: n2_cold_check, reldiff_interp
                    character(len=48) :: interp_debug_filename

                    call cold_plasma_refractive_index(omega, theta, DP%P, DP%S, DP%R, DP%L, &
                                                    n2_cold_check)
                    if (n2_cold_check /= 0.0d0) then
                        reldiff_interp = abs(n_ref_ind2 - n2_cold_check) / abs(n2_cold_check)
                    else
                        reldiff_interp = -1.0d0
                    end if

                end block

                ! Cache theta-bracket indices/weight before they go out of scope,
                ! so the Fsq diagnostic further down can reference them
                i_theta_lo_cached = i_theta_lo
                i_theta_hi_cached = i_theta_hi
                w_theta_cached    = w_theta
                end block

            end if   ! closes the outer HOT/EMIC interpolation block

            DP%n_ref_ind2 = n_ref_ind2

            !if (trim(DispersionType) == 'HOT' .or. trim(DispersionType) == 'EMIC') then
            !    if (n_ref_ind2 > N_max_hot**2) then
            !        cycle
            !    end if
            !end if

            ! -----------------------------------------------------------------
            ! Resonance cone check
            ! -----------------------------------------------------------------
            res_cone      = 'N'
            theta_max_eff = theta_max   ! degrees
            theta_cone    = -1.0d0

            if (trim(DispersionType) == 'COLD' .or. trim(DispersionType) == 'HOT') then
                if (DP%P / DP%S < 0) then
                    theta_cone = atan(sqrt(-DP%P / DP%S)) * 180.0 / PI

                    if (theta_cone*0.9 < theta_max .and. theta_cone*cone_backoff > theta_max) then
                        res_cone = 'Y'
                    else if (theta_cone*cone_backoff <= theta_max) then
                        theta_max_eff = theta_cone * cone_backoff
                        res_cone = 'Y'
                    end if
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

            if (trim(DispersionType) == 'HOT' .or. trim(DispersionType) == 'EMIC') then
                block
                    real*8 :: n2_cold_local, Fsq_cold_equiv
                    real*8 :: E_norm_sq_i, B_norm_sq_i
                    real*8 :: sin2_bfl_i, sin2_geom_i, n2_from_BFL
                    real*8 :: Fsq_raw_i, Fsq_lyons_cold
                    real*8 :: sin2_cold, k_sgn_d
                    complex(kind=8) :: E_cold_vec(3)
                    complex(kind=8) :: E_R_i, E_L_i, Theta_i, E_par_k_i
                    logical :: ok_cold
                    real*8, parameter :: c_undo = 299.7925d0
                    character(len=64) :: diagA_file

                    ! ---- signed k, as used in Theta_nk ----
                    if (k_vec >= 0.0d0) then
                        k_sgn_d =  1.0d0
                    else
                        k_sgn_d = -1.0d0
                    end if

                    ! ---- cold reference quantities at the SAME (omega, theta) ----
                    call cold_plasma_refractive_index(omega, theta, DP%P, DP%S, DP%R, DP%L, &
                                                      n2_cold_local)
                    call cold_plasma_polarisation(alpha_ne, nres, theta, k_vec, omega, &
                                                  BesselF1, BesselF2, BesselF3, &
                                                  n2_cold_local, DP%R, DP%L, DP%P, DP%S, &
                                                  Fsq_cold_equiv)
                    call cold_plasma_efield_diag(theta, n2_cold_local, DP%R, DP%L, DP%P, &
                                                 DP%S, DP%D, E_cold_vec, sin2_cold, ok_cold)

                    ! ---- quantities from the EFL/BFL actually used this step ----
                    E_norm_sq_i = real(DP%EFL(1)*conjg(DP%EFL(1)) &
                                     + DP%EFL(2)*conjg(DP%EFL(2)) &
                                     + DP%EFL(3)*conjg(DP%EFL(3)))
                    B_norm_sq_i = real(DP%BFL(1)*conjg(DP%BFL(1)) &
                                     + DP%BFL(2)*conjg(DP%BFL(2)) &
                                     + DP%BFL(3)*conjg(DP%BFL(3))) * c_undo**2

                    ! OLD route: sin^2 inferred from BFL and the interpolated n^2.
                    ! Kept for comparison — this is the quantity that carried the
                    ! n^2 inconsistency.
                    if (E_norm_sq_i > 0.0d0 .and. n_ref_ind2 > 0.0d0) then
                        sin2_bfl_i = B_norm_sq_i / (n_ref_ind2 * E_norm_sq_i)
                    else
                        sin2_bfl_i = -1.0d0
                    end if

                    ! NEW route: sin^2 from EFL geometry alone (khat = (sin,0,cos)
                    ! in WHAMP's basis, since e1 lies along k_perp).
                    sin2_geom_i = -1.0d0
                    if (E_norm_sq_i > 0.0d0) then
                        E_par_k_i = DP%EFL(1)*sin(theta) + DP%EFL(3)*cos(theta)
                        sin2_geom_i = 1.0d0 - real(E_par_k_i*conjg(E_par_k_i)) / E_norm_sq_i
                        sin2_geom_i = max(0.0d0, min(1.0d0, sin2_geom_i))
                    end if

                    ! Refractive index implied by BFL itself, for direct comparison
                    ! against the interpolated n_ref_ind2 and the cold value.
                    n2_from_BFL = -1.0d0
                    if (E_norm_sq_i > 0.0d0 .and. sin2_geom_i > 0.0d0) then
                        n2_from_BFL = B_norm_sq_i / (E_norm_sq_i * sin2_geom_i)
                    end if

                    ! ---- Fsq from WHAMP fields, new formula, gate bypassed ----
                    E_R_i = (DP%EFL(1) - (0.0d0,1.0d0)*DP%EFL(2)) / sqrt(2.0d0)
                    E_L_i = (DP%EFL(1) + (0.0d0,1.0d0)*DP%EFL(2)) / sqrt(2.0d0)
                    Theta_i = (E_R_i*BesselF3 + E_L_i*BesselF1) / sqrt(2.0d0) &
                            + k_sgn_d * DP%EFL(3)*BesselF2*cos(alpha_ne)/sin(alpha_ne)
                    if (E_norm_sq_i > 0.0d0 .and. sin2_geom_i > 0.0d0) then
                        Fsq_raw_i = real(Theta_i*conjg(Theta_i)) / (E_norm_sq_i * sin2_geom_i)
                    else
                        Fsq_raw_i = -1.0d0
                    end if

                    ! ---- TEST C: same Lyons formula, but on the COLD eigenvector.
                    ! Isolates formula error from field error. If this matches
                    ! Fsq_cold, the two formulations agree and any residual is in
                    ! WHAMP's EFL. If it does not, the discrepancy is in the
                    ! formula itself and has nothing to do with WHAMP. ----
                    Fsq_lyons_cold = -1.0d0
                    if (ok_cold .and. sin2_cold > 0.0d0) then
                        block
                            complex(kind=8) :: E_Rc, E_Lc, Theta_c
                            real*8 :: E_norm_c
                            E_norm_c = real(E_cold_vec(1)*conjg(E_cold_vec(1)) &
                                          + E_cold_vec(2)*conjg(E_cold_vec(2)) &
                                          + E_cold_vec(3)*conjg(E_cold_vec(3)))
                            E_Rc = (E_cold_vec(1) - (0.0d0,1.0d0)*E_cold_vec(2)) / sqrt(2.0d0)
                            E_Lc = (E_cold_vec(1) + (0.0d0,1.0d0)*E_cold_vec(2)) / sqrt(2.0d0)
                            Theta_c = (E_Rc*BesselF3 + E_Lc*BesselF1) / sqrt(2.0d0) &
                                    + k_sgn_d * E_cold_vec(3)*BesselF2 &
                                      *cos(alpha_ne)/sin(alpha_ne)
                            if (E_norm_c > 0.0d0) &
                                Fsq_lyons_cold = real(Theta_c*conjg(Theta_c)) &
                                                 / (E_norm_c * sin2_cold)
                        end block
                    end if

                    ! ---- log: equatorial slice only, to bound file size ----
                    !if (abs(lambda*180.0d0/PI) < 0.5d0) then
                    !    write(diagA_file,'(A,I0,A)') 'testA_sin2_rank', my_rank_cached, '.txt'
                    !    open(94, file=trim(diagA_file), status='unknown', position='append')
                    !    write(94,*) &
                    !        'E=', E, ' alpha=', alpha_ne*180.0d0/PI, &
                    !        ' nres=', nres, ' theta=', theta*180.0d0/PI, &
                    !        ' omega=', omega, ' ksgn=', k_sgn_d, &
                    !        ' n2_hot=', n_ref_ind2, ' n2_cold=', n2_cold_local, &
                    !        ' n2_bfl=', n2_from_BFL, &
                    !        ' sin2_bfl=', sin2_bfl_i, ' sin2_geom=', sin2_geom_i, &
                    !        ' sin2_cold=', sin2_cold, &
                    !        ' Fsq_raw=', Fsq_raw_i, ' Fsq_lyonscold=', Fsq_lyons_cold, &
                    !        ' Fsq_cold=', Fsq_cold_equiv, &
                    !        ' okcold=', ok_cold, ' path=', trim(whamp_path), &
                    !        ' EFL=', DP%EFL, ' Ecold=', E_cold_vec
                    !    close(94)
                    !end if
                end block
            end if

            !if (E > 1.9d0 .and. E < 2.0d0 .and. abs(alpha_ne*180.0d0/PI - 10.0d0) < 1.0d0 &
            !    .and. abs(lambda*180.0d0/PI) < 0.5d0 .and. nres == -1) then
            !    open(91, file='fsq_nresmin1.txt', status='unknown', position='append')
            !    write(91,*) 'DispType=', trim(DP%disp_type), ' nres=', nres, &
            !                ' argJ=', argJ, ' BesselF1=', BesselF1, ' BesselF2=', BesselF2, &
            !                ' BesselF3=', BesselF3, ' EFL=', DP%EFL, ' Fsq=', Fsq
            !    close(91)
            !end if

            ! -----------------------------------------------------------------
            ! Group velocity component along field line
            ! -----------------------------------------------------------------
            dw_dk_par = dw_dk / cos(theta)

            ! -----------------------------------------------------------------
            ! Diffusion coefficient integrand
            ! -----------------------------------------------------------------
            acct_status = 0
            acct_N      = -1.0d0

            if (theta_min > theta_max_eff) then
                acct_status = 5
            else if ((theta * 180.0/PI) > theta_max_eff) then
                acct_status = 4
            else if (abs(vpar - dw_dk_par) < 1.0d-3 * abs(vpar)) then
                acct_status = 3
            else if (Fsq <= 0.0d0) then
                acct_status = 6
            end if

            if (acct_status /= 0) then

                D_tmp(1) = 0.0
                D_tmp(2) = 0.0
                D_tmp(3) = 0.0

            else

                ! Resonance broadening factor
                Delta_n = PI/2.0 * 1.0/cos(theta) / ((vpar/c)**3.0) * Fsq &
                          * ((real(nres)*real(sqe)*omega_e/(gamma1*omega) &
                             - (sin(alpha_ne))**2.0)**2.0) &
                          / abs(1.0 - dw_dk_par/vpar)

                call PowerSpectralDensity(L, lambda, omega_e_eq, omega, &
                          omega_m, omega_lc, omega_uc, d_omega, PSD)
                G1 = omega_e * PSD

                block
                    real*8 :: theta_max_N
                    if (ConeTruncateN) then
                        theta_max_N = theta_max_eff
                    else
                        theta_max_N = theta_max
                    end if
                    acct_iroot     = i
                    acct_nres_ctx  = nres
                    acct_theta_res = theta * 180.0d0/PI
                    call int_N_omega(alpha, E, L, lambda, Theta_function, &
                                    theta_max_N, nres, omega, k_vec, &
                                    omega_e, omega_pe, omega_i, omega_pi, &
                                    res_cone, N_of_omega, whamp_curve)
                end block

                acct_N = N_of_omega

                if (N_of_omega > 0.0d0) then
                    G2 = (exp(-((X - X_m1)/X_w1)**2.0) &
                        + theta_q * exp(-((X - X_m2)/X_w2)**2.0)) / N_of_omega
                else
                    ! No wave power in the permitted theta range at this omega:
                    ! the resonance cannot contribute.
                    G2 = 0.0d0
                end if

                dB_local = Bw
                G3 = omega_e / (gamma1**2.0) * (dB_local / B_local)**2.0

                G4 = sin(alpha_ne) * cos(alpha_ne) &
                     / (real(nres)*real(sqe)*omega_e/(gamma1*omega) &
                        - (sin(alpha_ne))**2.0)

                D_tmp(1) = G3 * Delta_n * G1 * G2
                D_tmp(2) = D_tmp(1) * G4
                D_tmp(3) = D_tmp(2) * G4

            end if

            !call acct_write(DispersionType, my_rank_cached, acct_status, i, nres, &
            !                L, E, alpha*180.0d0/PI, lambda*180.0d0/PI, &
            !                theta*180.0d0/PI, theta_max_eff, &
            !                omega, n_ref_ind2, Fsq, acct_N, theta_cone, D_tmp(1))

            WeightValue(1) = WeightValue(1) + D_tmp(1)
            WeightValue(2) = WeightValue(2) + D_tmp(2)
            WeightValue(3) = WeightValue(3) + D_tmp(3)

        end if

    end do
#ifdef PROFILE
    call system_clock(t_prof_2, prof_rate)
    t_weight_total = t_weight_total + real(t_prof_2 - t_prof_1, 8) / real(prof_rate, 8)
    n_weight_calls = n_weight_calls + 1
#endif
    return

end subroutine WeightFunction