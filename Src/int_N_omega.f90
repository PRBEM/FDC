! ==============================================================================
! SUBROUTINE: int_N_omega
! ==============================================================================
!
! PURPOSE:
!   Computes the normalisation factor N(omega) for the wave normal angle
!   distribution at a given resonant frequency omega. This factor appears
!   in the denominator of the G2 term in WeightFunction, ensuring the
!   wave normal angle distribution integrates to unity.
!
!   N(omega) is defined as:
!     N(omega) = integral[ sin(theta) * g(theta) * n^2 * c / |dw/dk| * dtheta ]
!   where g(theta) is the two-Gaussian wave normal angle distribution in
!   tan(theta), n^2 is the refractive index squared, and dw/dk is the
!   group velocity along the wave vector. The integral runs from theta_min
!   to theta_max_eff using the trapezoidal rule.
!
!   Dispersion properties (n^2 and dw/dk) are obtained via
!   Get_dispersion_properties from module dispersion, which currently
!   implements the cold plasma (Stix) formulation.
!
! INPUTS:
!   alpha          - Equatorial electron pitch angle [radians]
!   E              - Electron kinetic energy normalised to rest-mass energy
!                    (E = gamma - 1, dimensionless)
!   L              - Magnetic L-shell [Earth radii, dimensionless]
!   lambda         - Magnetic latitude [radians]
!   Theta_function - Array(7) defining the wave normal angle distribution:
!                      (1) theta_m1  : mean angle of 1st Gaussian [degrees]
!                      (2) theta_w1  : width of 1st Gaussian [degrees]
!                      (3) theta_m2  : mean angle of 2nd Gaussian [degrees]
!                      (4) theta_w2  : width of 2nd Gaussian [degrees]
!                      (5) theta_min : lower cutoff angle [degrees]
!                      (6) theta_max : upper cutoff angle [degrees]
!                      (7) theta_q   : relative amplitude of 2nd Gaussian
!   theta_max_eff  - Effective upper cutoff angle, potentially reduced to
!                    0.99 * resonance cone angle [degrees]
!   nres           - Resonance harmonic number [-]
!   omega          - Resonant wave frequency [rad/s]
!   omega_e        - Electron gyrofrequency [rad/s]
!   omega_pe       - Electron plasma frequency [rad/s]
!   omega_i        - Ion gyrofrequency [rad/s]
!   omega_pi       - Ion plasma frequency [rad/s]
!   res_cone       - 'Y' if theta_max_eff is near the resonance cone,
!                    'N' otherwise. Controls integration resolution.
!
! OUTPUTS:
!   N_of_omega     - Normalisation factor N(omega) [rad]
!
! CALLS:
!   Get_dispersion_properties - Computes n^2 and dw/dk (module dispersion)
!
! MODULE DEPENDENCIES:
!   constants  - PI, nint_N, nint_N_res_cone
!   dispersion - Get_dispersion_properties, DispersionProperties
!
! NOTES:
!   - If n^2 < 0 at a given theta point (evanescent wave), that point
!     contributes zero to the integral.
!   - A finer integration grid (nint_N_res_cone points) is used when
!     res_cone = 'Y', to handle the near-singularity of the dispersion
!     relation close to the resonance cone angle.
!   - omega_e, omega_pe, omega_i, omega_pi are now passed explicitly as
!     arguments rather than read from module variables. This replaces the
!     previous implicit dependency on module variables being pre-populated
!     by WeightFunction before this routine is called.
!
! MODIFIED: 2011-10-12 Xudong - revised refractive index and dw/dk calculation
!           2011-10-24 Xudong - switched to lower-sign branch of dispersion
! ==============================================================================

subroutine int_N_omega(alpha, E, L, lambda, Theta_function, theta_max_eff, &
                       nres, omega, omega_e, omega_pe, omega_i, omega_pi, &
                       res_cone, N_of_omega, whamp_curve)
use constants
use dispersion
implicit none

    real*8, intent(in)  :: alpha, E, L, lambda, Theta_function(7)
    real*8, intent(in)  :: theta_max_eff, omega
    real*8, intent(in)  :: omega_e, omega_pe, omega_i, omega_pi
    integer, intent(in) :: nres
    character, intent(in) :: res_cone
    real*8, intent(out) :: N_of_omega
    type(WHAMPCurve), intent(in) :: whamp_curve

    ! --- Wave normal angle distribution parameters ---
    real*8 :: theta_min, theta_max
    real*8 :: theta_m1, theta_w1, theta_m2, theta_w2, theta_q
    real*8 :: X_m1, X_w1, X_m2, X_w2

    ! --- Integration ---
    real*8 :: theta, X, step, Ntmp
    integer :: nint_N_tmp, i

    ! --- Dispersion ---
    type(DispersionProperties) :: DP
    real*8 :: k_vec

    ! -------------------------------------------------------------------------
    ! Integration resolution: finer near resonance cone
    ! -------------------------------------------------------------------------
    if (res_cone == 'Y') then
        nint_N_tmp = nint_N_res_cone
    else
        nint_N_tmp = nint_N
    end if

    ! -------------------------------------------------------------------------
    ! Unpack wave normal angle distribution parameters
    ! -------------------------------------------------------------------------
    theta_m1  = Theta_function(1)
    theta_w1  = Theta_function(2)
    theta_m2  = Theta_function(3)
    theta_w2  = Theta_function(4)
    theta_min = Theta_function(5)
    theta_max = theta_max_eff
    theta_q   = Theta_function(7)

    X_m1 = tan(theta_m1 / 180.0 * PI)
    X_w1 = tan(theta_w1 / 180.0 * PI)
    X_m2 = tan(theta_m2 / 180.0 * PI)
    X_w2 = tan(theta_w2 / 180.0 * PI)

    step = (theta_max - theta_min) / 180.0 * PI / (real(nint_N_tmp) - 1.0)

    ! -------------------------------------------------------------------------
    ! Trapezoidal integration over theta
    ! -------------------------------------------------------------------------
    N_of_omega = 0.0

    do i = 1, nint_N_tmp
        theta = theta_min / 180.0 * PI + step * (real(i) - 1.0)
        X     = tan(theta)

        ! Step 1: get n_ref_ind2 with k_vec=0
        call Get_dispersion_properties(omega, 0.0d0, theta, &
                                       omega_e, omega_pe, omega_i, omega_pi, &
                                       DP)

        ! Skip evanescent points (n^2 < 0)
        if (DP%n_ref_ind2 < 0.0) then
            Ntmp = 0.0
        else
            ! Step 2: compute actual k_vec, then get dw_dk with correct k_vec
            k_vec = sqrt(DP%n_ref_ind2) * omega / c
            call Get_dispersion_properties(omega, k_vec, theta, &
                                           omega_e, omega_pe, omega_i, omega_pi, &
                                           DP)

            Ntmp = sin(theta) &
                   * (exp(-((X - X_m1)/X_w1)**2.0) &
                      + theta_q * exp(-((X - X_m2)/X_w2)**2.0)) &
                   * DP%n_ref_ind2 * c / abs(DP%dw_dk)
        end if

        ! Trapezoidal weights
        if (i == 1 .or. i == nint_N_tmp) then
            N_of_omega = N_of_omega + Ntmp / 2.0
        else
            N_of_omega = N_of_omega + Ntmp
        end if
    end do

    N_of_omega = N_of_omega * step

    return

end subroutine int_N_omega