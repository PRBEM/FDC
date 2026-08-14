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
!   For DispersionType='HOT', calls get_whamp_dispersion once per lambda
!   point to obtain the hot plasma dispersion curve, which is then passed
!   down the call chain to WeightFunction and int_N_omega.
!
!
! INPUTS:
!   alpha       - Equatorial electron pitch angle [radians]
!   E           - Electron kinetic energy normalised to rest-mass energy
!                 (E = gamma - 1, dimensionless)
!   L           - Magnetic L-shell [Earth radii, dimensionless]
!   B_m         - Magnetic field strength at the mirror point [T]
!   lambda_min  - Lower bound of latitude integration (equator or loss cone
!                 boundary) [radians]
!   lambda_m    - Upper bound of latitude integration (mirror latitude) [radians]
!   nint_tmp    - Number of integration steps over lambda
!   nres        - Resonance harmonic number n (integer, e.g. -1, 0, 1, ...)
!
! OUTPUTS:
!   BavRate(1)  - Bounce-averaged D_aa [s^-1] (pitch-angle diffusion)
!   BavRate(2)  - Bounce-averaged D_ap [s^-1] (mixed diffusion)
!   BavRate(3)  - Bounce-averaged D_pp [s^-1] (momentum diffusion)
!
! CALLS:
!   B_field     - Returns magnetic field magnitude at (L, lambda)
!   ElectronDensity      - Returns electron density at (L, lambda)
!   get_whamp_dispersion - Returns hot plasma dispersion curve (HOT case only)
!   int_D_loc   - Returns the lambda-integrand of the diffusion coefficients
!
! MODULE DEPENDENCIES:
!   constants   - PI, qe, me, c, DispersionType
!   dispersion  - WHAMPCurve type, get_whamp_dispersion
!
! INTEGRATION METHOD:
!   D_aa, D_ap: Standard trapezoidal rule with step delta_lambda over
!               nint_tmp+1 points from lambda_min to lambda_m.
!
!   D_pp:       Uses a different scheme at each point, weighted by
!               2*B_m/|dB/dlambda| * cos(psi) * |delta_psi|, where psi is
!               the local pitch angle. This accounts for the mirror boundary
!               condition more carefully, particularly near lambda_m where
!               psi -> 0. The endpoint (j = nint_tmp+1) is treated separately
!               using a one-sided finite difference for dB/dlambda.
!               The final result is divided by 2.
!
! NOTES:
!   - The subroutine returns immediately with BavRate = 0 if delta_lambda <= 0,
!     i.e. if lambda_m <= lambda_min.
!   - psi(1) is set to PI/2 (equator, field-aligned particle).
!   - psi(nint_tmp+1) is set to 0 if B_local > B_m at the mirror point
!     (particle is at or beyond its mirror point).
!   - B_field is evaluated at lambda_m + delta_lambda (one step beyond the
!     mirror point) to allow a centred finite difference for dB at the last
!     interior point. This point is stored in Bfield(nint_tmp+2).
!   - nres is correctly passed as intent(in) here. Note that int_D_loc
!     contains a redundant bare integer declaration of nres alongside its
!     intent(in) declaration — this is a known bug to be fixed.
!
! ==============================================================================

subroutine Integration_lambda(lambda_m, lambda_min, nint_tmp, B_m, alpha, E, L, nres, BavRate, &
                               whamp_curves, lambda_grid, n_lambda_in)
use constants
use dispersion
implicit none

  real*8,intent(in)::alpha, E, L, B_m, lambda_min, lambda_m
  integer,intent(in)::nint_tmp, nres
  real*8,intent(out)::BavRate(3)

  integer, intent(in) :: n_lambda_in
  type(WHAMPCurve), intent(in) :: whamp_curves(n_lambda_in)
  real*8, intent(in) :: lambda_grid(n_lambda_in)  

  real*8::psi(nint_tmp+1)
  real*8::Bfield(nint_tmp+2)
  real*8::dB, B_local

  real*8::lambda_tmp
  real*8::delta_lambda   ! step of magnetic latitude in rad
  real*8::lambda

  integer::i, j

  real*8::D_integrad(3)

  ! HOT plasma: local plasma frequencies and WHAMP dispersion curve
  real*8 :: Ne_local, fce_khz, fpe_khz
  type(WHAMPCurve) :: curve

  real*8 :: t_prof_1, t_prof_2

  real*8 :: Theta_function(7)

  logical :: use_shared
  integer :: jmatch
  real*8, parameter :: LAMBDA_MATCH_TOL = 1.0d-8

  delta_lambda=(lambda_m-lambda_min)/real(nint_tmp)     ! step of magnetic latitude in rad

  BavRate(1) = 0.0
  BavRate(2) = 0.0
  BavRate(3) = 0.0

  ! Shared curves only valid for nres/=0 (nint grid); nres=0 uses nint_zero_res
  use_shared = (DispersionType == 'HOT' .and. nres /= 0 .and. &
                n_lambda_in == nint_tmp+1 .and. size(whamp_curves) == nint_tmp+1)

  if (delta_lambda > 0.0) then
    call cpu_time(t_prof_1)
  
    do i=1, nint_tmp+1

        lambda_tmp = lambda_min + delta_lambda * real(i-1)

        call b_field(L, lambda_tmp, B_local)
        Bfield(i) = B_local

        if (i == 1) then
            psi(i) = PI/2.0
        else if (i == nint_tmp+1 .and. B_local > B_m) then
            psi(i) = 0.0
        else
            psi(i) = asin(sqrt(1 - B_local/B_m) / cos(alpha))
        end if
    end do

	call b_field(L, lambda_m + delta_lambda, B_local)
    Bfield(nint_tmp+2) = B_local

    curve%n_k     = 0
    curve%n_theta = 0
    curve%fce_khz = 0.0d0
    curve%fpe_khz = 0.0d0

    ! -------------------------------------------------------------------------
    ! Main lambda loop: reuse cached curves (HOT) or do nothing extra (COLD)
    ! -------------------------------------------------------------------------
    do j = 1, nint_tmp+1
        lambda = lambda_min + delta_lambda * real(j-1)
        if (DispersionType == 'HOT') then
            if (use_shared) then
                ! Direct lookup - lambda grids match exactly by construction
                jmatch = j
                if (abs(lambda_grid(jmatch) - lambda) > LAMBDA_MATCH_TOL) then
                    ! Safety fallback: grids didn't match as expected, recompute fresh
                    call B_field(L, lambda, B_local)
                    call ElectronDensity(L, lambda, Ne_local)
                    fce_khz = qe * B_local / me / (2.0d0 * PI * 1.0d3)
                    fpe_khz = sqrt(Ne_local * qe * qe / (C0 * me)) / (2.0d0 * PI * 1.0d3)
                    call Theta_wave(L, lambda, Theta_function)
                    call get_whamp_dispersion(fce_khz, fpe_khz, &
                                              Theta_function(5)*PI/180.0d0, &
                                              Theta_function(6)*PI/180.0d0, &
                                              curve)
                    call int_D_loc(alpha, E, L, nres, lambda, D_integrad, curve)
                    !if (E > 1.9d0 .and. E < 2.0d0 .and. (abs(alpha*180.0d0/PI - 10.0d0) < 0.5d0 .or. &
                    !                                    abs(alpha*180.0d0/PI - 80.0d0) < 0.5d0)) then
                    !    open(92, file='dintegrad_debug.txt', status='unknown', position='append')
                    !    write(92,*) 'alpha(deg)=', alpha*180.0d0/PI, ' lambda(deg)=', lambda*180.0d0/PI, &
                    !                ' D_integrad(1)=', D_integrad(1), ' D_integrad(2)=', D_integrad(2)
                    !    close(92)
                    !end if
                    deallocate(curve%theta, curve%k, curve%k_whamp, &
                               curve%omega, curve%dw_dk, curve%valid, curve%EFL)
                else
                    call int_D_loc(alpha, E, L, nres, lambda, D_integrad, whamp_curves(jmatch))
                    !if (E > 1.9d0 .and. E < 2.0d0 .and. (abs(alpha*180.0d0/PI - 10.0d0) < 0.5d0 .or. &
                    !                                    abs(alpha*180.0d0/PI - 80.0d0) < 0.5d0)) then
                    !    open(92, file='dintegrad_debug.txt', status='unknown', position='append')
                    !    write(92,*) 'alpha(deg)=', alpha*180.0d0/PI, ' lambda(deg)=', lambda*180.0d0/PI, &
                    !                ' D_integrad(1)=', D_integrad(1), ' D_integrad(2)=', D_integrad(2)
                    !    close(92)
                    !end if
                end if
            else
                call B_field(L, lambda, B_local)
                call ElectronDensity(L, lambda, Ne_local)
                fce_khz = qe * B_local / me / (2.0d0 * PI * 1.0d3)
                fpe_khz = sqrt(Ne_local * qe * qe / (C0 * me)) / (2.0d0 * PI * 1.0d3)
                call Theta_wave(L, lambda, Theta_function)
                call get_whamp_dispersion(fce_khz, fpe_khz, &
                                          Theta_function(5)*PI/180.0d0, &
                                          Theta_function(6)*PI/180.0d0, &
                                          curve)
                call int_D_loc(alpha, E, L, nres, lambda, D_integrad, curve)
                    !if (E > 1.9d0 .and. E < 2.0d0 .and. (abs(alpha*180.0d0/PI - 10.0d0) < 0.5d0 .or. &
                    !                                    abs(alpha*180.0d0/PI - 80.0d0) < 0.5d0)) then
                    !    open(92, file='dintegrad_debug.txt', status='unknown', position='append')
                    !    write(92,*) 'alpha(deg)=', alpha*180.0d0/PI, ' lambda(deg)=', lambda*180.0d0/PI, &
                    !                ' D_integrad(1)=', D_integrad(1), ' D_integrad(2)=', D_integrad(2)
                    !    close(92)
                    !end if
                deallocate(curve%theta, curve%k, curve%k_whamp, &
                           curve%omega, curve%dw_dk, curve%valid, curve%EFL)
            end if
        else
            call int_D_loc(alpha, E, L, nres, lambda, D_integrad, curve)
            !if (E > 1.9d0 .and. E < 2.0d0 .and. (abs(alpha*180.0d0/PI - 10.0d0) < 0.5d0 .or. &
            !                                    abs(alpha*180.0d0/PI - 80.0d0) < 0.5d0)) then
            !    open(92, file='dintegrad_debug.txt', status='unknown', position='append')
            !    write(92,*) 'alpha(deg)=', alpha*180.0d0/PI, ' lambda(deg)=', lambda*180.0d0/PI, &
            !                ' D_integrad(1)=', D_integrad(1), ' D_integrad(2)=', D_integrad(2)
            !    close(92)
            !end if
        end if

        if (j == 1) then
             BavRate(1) = BavRate(1) + D_integrad(1) / 2.0
             BavRate(2) = BavRate(2) + D_integrad(2) / 2.0
             BavRate(3) = BavRate(3) + D_integrad(3) * B_m/Bfield(j) * tan(alpha)**2

        else if (j == nint_tmp+1) then

              dB = (Bfield(j+1) - Bfield(j-1)) / 2
              BavRate(1) = BavRate(1) + D_integrad(1) / 2.0
              BavRate(2) = BavRate(2) + D_integrad(2) / 2.0
              BavRate(3) = BavRate(3) + D_integrad(3) * 2*B_m/abs(dB) * cos(psi(j)) * abs(psi(j-1)-psi(j))


        else
              dB = (Bfield(j+1) - Bfield(j-1)) / 2
              BavRate(1) = BavRate(1) + D_integrad(1)
              BavRate(2) = BavRate(2) + D_integrad(2)
              BavRate(3) = BavRate(3) + D_integrad(3) * 2*B_m/abs(dB) * cos(psi(j)) * abs(psi(j-1)-psi(j+1))
        end if
    end do

    BavRate(1) = BavRate(1) * delta_lambda
    BavRate(2) = BavRate(2) * delta_lambda
    BavRate(3) = BavRate(3) * delta_lambda / 2
    call cpu_time(t_prof_2)
    t_intlam_total = t_intlam_total + (t_prof_2 - t_prof_1)
    n_intlam_calls = n_intlam_calls + 1
end if

end subroutine Integration_lambda
