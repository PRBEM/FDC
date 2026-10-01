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
!       -> Integration_lambda  (trapezoidal integration over lambda)
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
!   alpha   - Equatorial electron pitch angle [radians]
!   E       - Electron kinetic energy normalised to rest-mass energy
!             (E = gamma - 1, dimensionless)
!   L       - Magnetic L-shell [Earth radii, dimensionless]
!   nres    - Resonance harmonic number n (integer, e.g. -1, 0, 1, ...)
!
! OUTPUTS:
!   BavRate(1)  - Bounce-averaged D_aa [s^-1] (pitch-angle diffusion)
!   BavRate(2)  - Bounce-averaged D_ap [s^-1] (mixed diffusion)
!   BavRate(3)  - Bounce-averaged D_pp [s^-1] (momentum diffusion)
!
! CALLS:
!   Lambda_mirror_new   - Computes mirror latitude from equatorial pitch angle
!   B_field             - Returns magnetic field magnitude at (L, lambda)
!   Integration_lambda  - Integrates local diffusion rates over lambda
!
! MODULE DEPENDENCIES:
!   constants   - Physical constants and unit conversion factors (RAD, PI, ...),
!                 and integration parameters:
!                   lam_min, lam_max  : latitude integration bounds [degrees]
!                   nint              : number of lambda steps for nres /= 0
!                   nint_zero_res     : number of lambda steps for nres = 0
!
! NOTES:
!   - The mirror latitude lambda_m is clamped to [lam_min, lam_max] if it
!     falls outside that range. This handles particles whose mirror point
!     lies outside the region where wave data is available.
!   - A finer lambda grid (nint_zero_res) is used for the nres = 0 resonance,
!     which typically requires higher resolution due to its different frequency
!     dependence.
!   - The bounce period factor S(alpha) ~ 1.30 - 0.56*sin(alpha) is an
!     analytic approximation for a dipole field (see e.g. Schulz & Lanzerotti
!     1974). A = 1/S(alpha) normalises the latitude integral to a proper
!     bounce average.
!   - B_m (mirror field strength) is evaluated at the original lambda_m before
!     clamping. If lambda_m is clamped, B_m may not correspond to the actual
!     integration endpoint — this is a known approximation.
!
! ==============================================================================

subroutine D_ba(alpha,E,L,nres,BavRate,whamp_curves,lambda_grid,n_lambda_in)
use constants
use dispersion
implicit none

  real*8,intent(in)::alpha,E,L
  integer,intent(in)::nres
  real*8,intent(out)::BavRate(3)

  integer, intent(in) :: n_lambda_in
  type(WHAMPCurve), intent(in) :: whamp_curves(n_lambda_in)
  real*8, intent(in) :: lambda_grid(n_lambda_in)

  real*8::B_m
  real*8::B_local
  real*8::lambda_tmp
  
  integer::nint_tmp
  real*8::lambda_min
  real*8::lambda_max
  real*8::lambda_m
  real*8::lambda_mdeg
  integer::i,j

  real*8::lambda
  real*8::D_integrad(3)
  real*8::A   ! 1/S(alpha), where S(alpha) is proportional to bounce period
  real*8::delta_lambda   ! step of magnetic latitude in rad
  real*8::Bw                  ! wave amplitude at lambda [T] (local mode)
  real*8::Theta_function(7)   ! wave normal angle parameters (local mode)

  ! LOCAL MODE: return the local D at lambda = lam_min (equator), no bounce
  ! averaging and no 1/S(alpha) factor. D_local_n returns the theta-integrated
  ! local coefficients directly, so it writes straight into BavRate.
  ! whamp_curves(1) sits at lambda = lam_min (COLD: dummy; HOT/EMIC: all
  ! shared curves are at lam_min when lam_min = lam_max).
  if (local_mode) then
    lambda = lam_min * RAD
    call Bw_wave(L, lambda, Bw)
    call Theta_wave(L, lambda, Theta_function)
    call D_local_n(alpha, E, L, nres, lambda, Bw, Theta_function, &
                   BavRate, whamp_curves(1))
    return
  end if

  A=1.0/(1.30-0.56*sin(alpha))          ! 1/S(alpha)

  call Lambda_mirror_new(alpha,lambda_m)    ! determine the mirror latitude
	call B_field(L, 0.0d0, B_m)
  B_m = B_m / sin(alpha)**2

  ! determine the range of magnetic latitude
  lambda_min=lam_min*RAD                 ! in rad
  lambda_mdeg=lambda_m/RAD                ! in degree

  if(lambda_mdeg>=lam_max) then                ! mirror latitude larger than maximum latitude
    lambda_m=lam_max*RAD
  else if(lambda_mdeg<=lam_min) then           ! mirror latitude smaller than minimum latitude
         lambda_m=lam_min*RAD
  else
	     lambda_m=lambda_m
  end if
  

  if(nres==0) then
  	nint_tmp=nint_zero_res
  else
        nint_tmp=nint
  end if


    call Integration_lambda(lambda_m,lambda_min,nint_tmp,B_m,alpha,E,L,nres,BavRate, &
                           whamp_curves,lambda_grid,n_lambda_in)
    BavRate(1)=BavRate(1)*A
    BavRate(2)=BavRate(2)*A
    BavRate(3)=BavRate(3)*A

  return

end subroutine D_ba
