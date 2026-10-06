! ==============================================================================
! SUBROUTINE: D_ba
! ==============================================================================
!
! PURPOSE:
!   Computes the fully bounce-averaged quasi-linear diffusion coefficients
!   (D_aa, D_ap, D_pp) for a single resonance harmonic (nres), for an
!   electron with equatorial pitch angle alpha, energy E, at L-shell L.
!
!   This is the top-level routine for a single harmonic. The caller is
!   responsible for summing over resonance harmonics to get the total
!   diffusion coefficients.
!
!   The full call chain is:
!     D_ba
!       -> Integration_lambda  (integration over lambda)
!          -> int_D_loc        (field-line geometry factors at each lambda)
!             -> D_local_n     (integration over wave normal angle theta)
!                -> WeightFunction  (integrand at a single theta point)
!
!   The bounce-averaging follows the standard quasi-linear formulation:
!     <D> = A * integral[ D_local * geometry_factors * dlambda ]
!   where A = 1/S(alpha) is the inverse of a function proportional to the
!   bounce period, approximated as:
!     S(alpha) ~ 1.30 - 0.56 * sin(alpha)
!
!   The formulation follows:
!     Albert, J. M. (2007), Simple approximations of quasi-linear diffusion
!     coefficients, J. Geophys. Res., 112, A12202, doi:10.1029/2007JA012551
!   which builds on the PADIE framework of:
!     Glauert, S. A., and R. B. Horne (2005), Calculation of pitch angle and
!     energy diffusion coefficients with the PADIE code, J. Geophys. Res.,
!     110, A04206, doi:10.1029/2004JA010851
!
! INPUTS:
!   alpha        - Equatorial electron pitch angle [radians]
!   E            - Electron kinetic energy normalised to rest-mass energy
!                  (E = gamma - 1, dimensionless)
!   L            - Magnetic L-shell [Earth radii, dimensionless]
!   nres         - Resonance harmonic number n (integer, e.g. -1, 0, 1, ...)
!   whamp_curves - Pre-computed WHAMP dispersion curves, one per latitude of
!                  lambda_grid (HOT/EMIC; a dummy entry for COLD)
!   lambda_grid  - Latitudes at which whamp_curves were computed [rad]
!   n_lambda_in  - Number of entries in whamp_curves / lambda_grid
!
! OUTPUTS:
!   BavRate(1)  - Bounce-averaged D_aa [s^-1] (pitch-angle diffusion)
!   BavRate(2)  - Bounce-averaged D_ap [s^-1] (mixed diffusion)
!   BavRate(3)  - Bounce-averaged D_pp [s^-1] (momentum diffusion)
!   In local mode (see below) these are the LOCAL coefficients at the
!   equator, without bounce averaging.
!
! CALLS:
!   Lambda_mirror_new   - Computes mirror latitude from equatorial pitch angle
!   B_field             - Returns magnetic field magnitude at (L, lambda)
!   Integration_lambda  - Integrates local diffusion rates over lambda
!   Bw_wave, Theta_wave, D_local_n  - local mode only
!
! MODULE DEPENDENCIES:
!   constants   - Physical constants and unit conversion factors (RAD, PI, ...),
!                 and integration parameters:
!                   lam_min, lam_max  : latitude integration bounds [degrees]
!                   nint              : number of lambda steps for nres /= 0
!                   nint_zero_res     : number of lambda steps for nres = 0
!                   local_mode        : .true. when lam_min = lam_max = 0
!
! NOTES:
!   - LOCAL MODE: lam_min = lam_max = 0 (set in FDC_MPI) returns the local
!     equatorial coefficients, i.e. D_local_n's output without bounce
!     averaging and without the 1/S(alpha) factor. Used for comparisons with
!     codes that compute local coefficients, and for convergence tests that
!     should not involve the latitude integration.
!   - The mirror latitude lambda_m is clamped to [lam_min, lam_max] if it
!     falls outside that range. This handles particles whose mirror point
!     lies outside the region where wave data is available.
!   - B_m is the field at the mirror point, which follows from conservation
!     of the magnetic moment: B_m = B_eq/sin^2(alpha). It must NOT be taken
!     as B_field(L, lambda_m): Lambda_mirror_new shortens lambda_m by a
!     factor 0.999 to keep the integration away from the singular endpoint,
!     which would make B_m too small. B_m enters the variable transformation
!     used for D_pp in Integration_lambda (Orlova & Shprits 2011), where the
!     value at the equator is set analytically assuming an exact B_m, so an
!     inconsistent B_m biases D_pp. -DOLD_BM restores the previous behaviour
!     for comparison runs.
!   - A finer lambda grid (nint_zero_res) is used for the nres = 0 resonance,
!     which typically requires higher resolution due to its different frequency
!     dependence.
!   - The bounce period factor S(alpha) ~ 1.30 - 0.56*sin(alpha) is an
!     analytic approximation for a dipole field (Hamlin et al. 1961; see e.g.
!     Schulz & Lanzerotti 1974). It is the least accurate of the standard
!     parametrisations; Orlova & Shprits (2011), Table I, give three-term
!     forms with maximum deviations below 0.1% against roughly 6% for this
!     one. Worth revisiting if bounce-averaged absolute values matter.
!
! ==============================================================================

subroutine D_ba(alpha, E, L, nres, BavRate, whamp_curves, lambda_grid, n_lambda_in)
use constants
use dispersion
implicit none

    real*8,  intent(in)  :: alpha, E, L
    integer, intent(in)  :: nres
    real*8,  intent(out) :: BavRate(3)

    integer,          intent(in) :: n_lambda_in
    type(WHAMPCurve), intent(in) :: whamp_curves(n_lambda_in)
    real*8,           intent(in) :: lambda_grid(n_lambda_in)

    real*8  :: B_m                  ! field at the mirror point [T]
    real*8  :: lambda_min           ! lower integration bound [rad]
    real*8  :: lambda_m             ! mirror latitude, clamped [rad]
    real*8  :: lambda_mdeg          ! mirror latitude before clamping [deg]
    real*8  :: lambda               ! latitude (local mode) [rad]
    real*8  :: A                    ! 1/S(alpha), S proportional to bounce period
    real*8  :: Bw                   ! wave amplitude at lambda [T] (local mode)
    real*8  :: Theta_function(7)    ! wave normal angle parameters (local mode)
    integer :: nint_tmp             ! number of lambda steps for this harmonic

    ! --------------------------------------------------------------------
    ! LOCAL MODE: local D at the equator, no bounce averaging.
    ! whamp_curves(1) is at lambda = lam_min (COLD: dummy; HOT/EMIC: all
    ! shared curves lie at lam_min when lam_min = lam_max).
    ! --------------------------------------------------------------------
    if (local_mode) then
        lambda = lam_min * RAD
        call Bw_wave(L, lambda, Bw)
        call Theta_wave(L, lambda, Theta_function)
        call D_local_n(alpha, E, L, nres, lambda, Bw, Theta_function, &
                       BavRate, whamp_curves(1))
        return
    end if

    A = 1.0 / (1.30 - 0.56*sin(alpha))        ! 1/S(alpha)

    call Lambda_mirror_new(alpha, lambda_m)   ! mirror latitude (shortened by 0.999)

#ifdef OLD_BM
    call B_field(L, lambda_m, B_m)            ! field at the shortened mirror latitude
#else
    call B_field(L, 0.0d0, B_m)
    B_m = B_m / sin(alpha)**2                 ! exact mirror field (magnetic moment)
#endif

    ! latitude integration range
    lambda_min  = lam_min * RAD               ! [rad]
    lambda_mdeg = lambda_m / RAD              ! [deg]

    if (lambda_mdeg >= lam_max) then          ! mirror point beyond the wave model
        lambda_m = lam_max * RAD
    else if (lambda_mdeg <= lam_min) then     ! mirror point below the lower bound
        lambda_m = lam_min * RAD
    end if

    if (nres == 0) then
        nint_tmp = nint_zero_res
    else
        nint_tmp = nint
    end if

    call Integration_lambda(lambda_m, lambda_min, nint_tmp, B_m, alpha, E, L, nres, &
                            BavRate, whamp_curves, lambda_grid, n_lambda_in)

    BavRate(1) = BavRate(1) * A
    BavRate(2) = BavRate(2) * A
    BavRate(3) = BavRate(3) * A

    return

end subroutine D_ba