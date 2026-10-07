! ==============================================================================
! SUBROUTINE: Integration_lambda
! ==============================================================================
!
! PURPOSE:
!   Computes the bounce-averaged quasi-linear diffusion coefficients
!   (D_aa, D_ap, D_pp) for a single resonance harmonic by integrating
!   the local diffusion rate integrand over magnetic latitude (lambda),
!   from lambda_min to lambda_m (the mirror latitude).
!
!   For DispersionType = 'HOT', the hot plasma dispersion curve is
!   needed at every lambda point. Where the shared curves passed in from
!   sum_n match the lambda grid, they are reused; otherwise the curve is
!   recomputed with get_whamp_dispersion and released again.
!
! INPUTS:
!   alpha        - Equatorial electron pitch angle [radians]
!   E            - Electron kinetic energy normalised to rest-mass energy
!                  (E = gamma - 1, dimensionless)
!   L            - Magnetic L-shell [Earth radii, dimensionless]
!   B_m          - Magnetic field strength at the mirror point [T], equal to
!                  B_eq/sin^2(alpha) by conservation of the magnetic moment
!                  (see D_ba)
!   lambda_min   - Lower bound of latitude integration [radians]
!   lambda_m     - Upper bound of latitude integration (mirror latitude,
!                  shortened by 0.999 in Lambda_mirror_new) [radians]
!   nint_tmp     - Number of integration steps over lambda
!   nres         - Resonance harmonic number n (integer, e.g. -1, 0, 1, ...)
!   whamp_curves - Shared WHAMP dispersion curves, one per lambda grid point
!   lambda_grid  - Latitudes at which whamp_curves were computed [radians]
!   n_lambda_in  - Number of entries in whamp_curves / lambda_grid
!
! OUTPUTS:
!   BavRate(1)  - Bounce-averaged D_aa [s^-1] (pitch-angle diffusion)
!   BavRate(2)  - Bounce-averaged D_ap [s^-1] (mixed diffusion)
!   BavRate(3)  - Bounce-averaged D_pp [s^-1] (momentum diffusion)
!   (all still to be multiplied by 1/S(alpha) in D_ba)
!
! CALLS:
!   B_field              - Magnetic field magnitude at (L, lambda)
!   ElectronDensity      - Electron density at (L, lambda)       (HOT)
!   Theta_wave           - Wave normal angle parameters           (HOT)
!   get_whamp_dispersion - Hot plasma dispersion curve            (HOT)
!   int_D_loc            - Lambda-integrand of the diffusion coefficients
!   dump_whamp_surface   - Diagnostic dump of a shared curve      (-DDIAG)
!
! MODULE DEPENDENCIES:
!   constants   - PI, qe, me, c, C0, DispersionType
!   dispersion  - WHAMPCurve type, get_whamp_dispersion
!
! INTEGRATION METHOD:
!   D_aa, D_ap: composite trapezoidal rule with step delta_lambda over
!               nint_tmp+1 points from lambda_min to lambda_m. Both
!               integrands are finite at the mirror point (see int_D_loc),
!               so no special treatment is needed.
!
!   D_pp:       the integrand diverges like 1/cos(alpha_local) at the mirror
!               point, so the integration variable is changed to the
!               bounce-phase variable psi of Orlova & Shprits (2011),
!               Phys. Plasmas 18, 092904, doi:10.1063/1.3638137:
!
!                 1 - B(lambda)/B_m = cos^2(alpha_eq) * sin^2(psi)   [their Eq. 3]
!
! NOTES:
!   - Returns immediately with BavRate = 0 if delta_lambda <= 0, i.e. if
!     lambda_m <= lambda_min.
!   - psi(nint_tmp+1) is set to 0 if B_local > B_m at the last point, i.e.
!     if the grid has run past the mirror point.
!   - B_field is evaluated one step beyond lambda_m and stored in
!     Bfield(nint_tmp+2), so that the last interior point can use a centred
!     difference for dB.
!   - The argument of the asin for psi is clamped to [0, 1]: with the exact
!     B_m it is exactly 1 at the equator, where rounding can push it just
!     above 1 and produce a NaN.
!
! ==============================================================================

subroutine Integration_lambda(lambda_m, lambda_min, nint_tmp, B_m, alpha, E, L, nres, BavRate, &
                              whamp_curves, lambda_grid, n_lambda_in)
use constants
use dispersion
implicit none

    real*8,  intent(in)  :: alpha, E, L, B_m, lambda_min, lambda_m
    integer, intent(in)  :: nint_tmp, nres
    real*8,  intent(out) :: BavRate(3)

    integer,          intent(in) :: n_lambda_in
    type(WHAMPCurve), intent(in) :: whamp_curves(n_lambda_in)
    real*8,           intent(in) :: lambda_grid(n_lambda_in)

    real*8  :: psi(nint_tmp+1)      ! bounce-phase variable, pi/2 at the equator
    real*8  :: Bfield(nint_tmp+2)   ! B at the grid points, plus one beyond
    real*8  :: dB, B_local
    real*8  :: delta_lambda         ! step in magnetic latitude [rad]
    real*8  :: lambda, lambda_tmp
    real*8  :: sin_psi              ! argument of the asin, clamped to [0,1]
    real*8  :: D_integrad(3)
    integer :: i, j

    type(WHAMPCurve) :: curve       ! freshly computed curve (HOT/EMIC)
    logical :: use_shared           ! shared curves match this lambda grid
    integer :: jmatch
    real*8, parameter :: LAMBDA_MATCH_TOL = 1.0d-8

    integer(8) :: t_prof_1, t_prof_2, prof_rate
    logical, save :: surf_dumped = .false.

    delta_lambda = (lambda_m - lambda_min) / real(nint_tmp)

    BavRate(1) = 0.0
    BavRate(2) = 0.0
    BavRate(3) = 0.0

    ! Shared curves are only valid for nres /= 0: nres = 0 uses the finer
    ! nint_zero_res grid, for which no curves were pre-computed.
    use_shared = ((DispersionType == 'HOT' .or. DispersionType == 'EMIC') .and. nres /= 0 .and. &
                  n_lambda_in == nint_tmp+1 .and. size(whamp_curves) == nint_tmp+1)

    if (delta_lambda > 0.0) then
#ifdef PROFILE
        call system_clock(t_prof_1)
#endif

        ! ---- field and bounce-phase variable at every grid point ----
        do i = 1, nint_tmp+1
            lambda_tmp = lambda_min + delta_lambda * real(i-1)
            call B_field(L, lambda_tmp, B_local)
            Bfield(i) = B_local

            if (i == 1) then
                psi(i) = PI/2.0                       ! equator
            else if (i == nint_tmp+1 .and. B_local > B_m) then
                psi(i) = 0.0                          ! past the mirror point
            else
                sin_psi = sqrt(max(0.0d0, 1.0d0 - B_local/B_m)) / cos(alpha)
                psi(i)  = asin(min(1.0d0, sin_psi))
            end if
        end do

        ! one step beyond the mirror point, for the centred difference
        call B_field(L, lambda_m + delta_lambda, B_local)
        Bfield(nint_tmp+2) = B_local

        curve%n_k     = 0
        curve%n_theta = 0
        curve%fce_khz = 0.0d0
        curve%fpe_khz = 0.0d0

        ! ---- main lambda loop ----
        do j = 1, nint_tmp+1
            lambda = lambda_min + delta_lambda * real(j-1)

            if (DispersionType == 'HOT' .or. DispersionType == 'EMIC') then
                if (use_shared) then
                    jmatch = j      ! lambda grids match by construction
                    if (abs(lambda_grid(jmatch) - lambda) > LAMBDA_MATCH_TOL) then
                        ! safety fallback: grids did not match, recompute
                        call whamp_curve_at(lambda, curve)
                        call int_D_loc(alpha, E, L, nres, lambda, D_integrad, curve)
                        call release_curve(curve)
                    else
#ifdef DIAG
                        ! Diagnostic: writes the shared WHAMP dispersion surface
                        ! once per process (whampsurf_*.txt).
                        if (.not. surf_dumped) then
                            call dump_whamp_surface(L, lambda, whamp_curves(jmatch), 'shared')
                            surf_dumped = .true.
                        end if
#endif
                        call int_D_loc(alpha, E, L, nres, lambda, D_integrad, &
                                       whamp_curves(jmatch))
                    end if
                else
                    call whamp_curve_at(lambda, curve)
                    call int_D_loc(alpha, E, L, nres, lambda, D_integrad, curve)
                    call release_curve(curve)
                end if
            else
                ! COLD: curve is a dummy argument, not used downstream
                call int_D_loc(alpha, E, L, nres, lambda, D_integrad, curve)
            end if

            ! ---- accumulate ----
            if (j == 1) then
                ! equator: l'Hopital limit of the D_pp integrand (O&S Eq. 7)
                BavRate(1) = BavRate(1) + D_integrad(1) / 2.0
                BavRate(2) = BavRate(2) + D_integrad(2) / 2.0
                BavRate(3) = BavRate(3) + D_integrad(3) * B_m/Bfield(j) * tan(alpha)**2

            else if (j == nint_tmp+1) then
                ! last point: one-sided difference in psi
                dB = (Bfield(j+1) - Bfield(j-1)) / 2.0
                BavRate(1) = BavRate(1) + D_integrad(1) / 2.0
                BavRate(2) = BavRate(2) + D_integrad(2) / 2.0
                BavRate(3) = BavRate(3) + D_integrad(3) * 2.0*B_m/abs(dB) &
                             * cos(psi(j)) * abs(psi(j-1) - psi(j))

            else
                dB = (Bfield(j+1) - Bfield(j-1)) / 2.0
                BavRate(1) = BavRate(1) + D_integrad(1)
                BavRate(2) = BavRate(2) + D_integrad(2)
                BavRate(3) = BavRate(3) + D_integrad(3) * 2.0*B_m/abs(dB) &
                             * cos(psi(j)) * abs(psi(j-1) - psi(j+1))
            end if
        end do

        BavRate(1) = BavRate(1) * delta_lambda
        BavRate(2) = BavRate(2) * delta_lambda
        BavRate(3) = BavRate(3) * delta_lambda / 2.0

#ifdef PROFILE
        call system_clock(t_prof_2, prof_rate)
        t_intlam_total = t_intlam_total + real(t_prof_2 - t_prof_1, 8) / real(prof_rate, 8)
        n_intlam_calls = n_intlam_calls + 1
#endif
    end if

contains

    ! Computes the WHAMP dispersion curve for the local plasma at lambda.
    subroutine whamp_curve_at(lam, c_out)
        real*8,           intent(in)  :: lam
        type(WHAMPCurve), intent(out) :: c_out
        real*8 :: B_loc, Ne_loc, fce_khz, fpe_khz, Theta_function(7)

        call B_field(L, lam, B_loc)
        call ElectronDensity(L, lam, Ne_loc)
        fce_khz = qe * B_loc / me / (2.0d0 * PI * 1.0d3)
        fpe_khz = sqrt(Ne_loc * qe * qe / (C0 * me)) / (2.0d0 * PI * 1.0d3)
        call Theta_wave(L, lam, Theta_function)
        call get_whamp_dispersion(fce_khz, fpe_khz, &
                                  Theta_function(5)*PI/180.0d0, &
                                  Theta_function(6)*PI/180.0d0, &
                                  c_out)
    end subroutine whamp_curve_at

    subroutine release_curve(c_in)
        type(WHAMPCurve), intent(inout) :: c_in
        deallocate(c_in%theta, c_in%k, c_in%k_whamp, &
                   c_in%omega, c_in%dw_dk, c_in%valid, c_in%EFL)
    end subroutine release_curve

end subroutine Integration_lambda