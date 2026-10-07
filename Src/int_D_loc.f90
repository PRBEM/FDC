! ==============================================================================
! SUBROUTINE: int_D_loc
! ==============================================================================
!
! PURPOSE:
!   Computes the integrand for the bounce-averaging integral over magnetic
!   latitude (lambda), by taking the local diffusion coefficients from
!   D_local_n and applying the field-line geometry factors appropriate for
!   a dipole magnetic field.
!
!   This routine sits at the third level of the call chain:
!     sum_n -> D_ba -> Integration_lambda -> int_D_loc -> D_local_n -> WeightFunction
!
!   The geometric factors applied are:
!     D_aa integrand: D_aa_local / cos(alpha_ne) * (tan(alpha)/tan(alpha_ne))^2
!                     * cos(lambda) * sqrt(1 + 3*sin^2(lambda))
!     D_ap integrand: D_ap_local * tan(alpha)/sin(alpha_ne)
!                     * cos(lambda) * sqrt(1 + 3*sin^2(lambda))
!     D_pp integrand: D_pp_local * cos(alpha)
!                     * cos(lambda) * sqrt(1 + 3*sin^2(lambda))
!
!   where cos(lambda) * sqrt(1 + 3*sin^2(lambda)) is the field-line arc
!   length element ds/dlambda for a dipole field, and the remaining factors
!   transform the diffusion coefficients between the local and equatorial
!   pitch angle frames (the Jacobians (d alpha_eq/d alpha)^2 and
!   d alpha_eq/d alpha, together with the 1/cos(alpha_local) from ds/v_par).
!
!   Behaviour at the mirror point:
!     D_aa integrand -> 0
!     D_ap integrand -> finite
!     D_pp integrand -> diverges as 1/cos(alpha_local)
!   Only D_pp therefore needs the change of variables applied in
!   Integration_lambda; D_aa and D_ap are integrated directly over lambda.
!
! INPUTS:
!   alpha        - Equatorial electron pitch angle [radians]
!   E            - Electron kinetic energy normalised to rest-mass energy
!                  (E = gamma - 1, dimensionless)
!   L            - Magnetic L-shell [Earth radii, dimensionless]
!   nres         - Resonance harmonic number n (integer, e.g. -1, 0, 1, ...)
!   lambda       - Magnetic latitude at this evaluation point [radians]
!   whamp_curve  - Hot plasma dispersion curve from WHAMP, passed through to
!                  D_local_n (a dummy entry for COLD)
!
! OUTPUTS:
!   D_integrad(1)  - D_aa integrand at this lambda [s^-1 rad^-1]
!   D_integrad(2)  - D_ap integrand at this lambda [s^-1 rad^-1]
!   D_integrad(3)  - D_pp integrand at this lambda [s^-1 rad^-1]
!
! CALLS:
!   Local_alpha  - Converts equatorial pitch angle to local pitch angle
!   Bw_wave      - Returns RMS wave amplitude at (L, lambda)
!   Theta_wave   - Returns wave normal angle distribution parameters at (L, lambda)
!   D_local_n    - Returns local diffusion coefficients integrated over theta
!
! MODULE DEPENDENCIES:
!   dispersion   - WHAMPCurve type
!
! NOTES:
!   - Bw and Theta_function are obtained here by calling Bw_wave and
!     Theta_wave, rather than being passed in from Integration_lambda, so the
!     wave parameters are re-evaluated at every lambda point. They depend
!     only on (L, lambda), so they could be computed once per lambda in
!     Integration_lambda and passed down.
!
! ==============================================================================

subroutine int_D_loc(alpha, E, L, nres, lambda, D_integrad, whamp_curve)
use dispersion
implicit none

    real*8,  intent(in)  :: alpha, E, L, lambda
    integer, intent(in)  :: nres
    real*8,  intent(out) :: D_integrad(3)
    type(WHAMPCurve), intent(in) :: whamp_curve

    real*8 :: Bw                  ! wave amplitude at lambda [T]
    real*8 :: Theta_function(7)   ! wave normal angle distribution parameters
    real*8 :: alpha_ne            ! local pitch angle at lambda [rad]
    real*8 :: D_local(3)          ! local D_aa, D_ap, D_pp [s^-1]
    real*8 :: arc_element         ! ds/dlambda for a dipole field

    call Local_alpha(alpha, lambda, alpha_ne)
    call Bw_wave(L, lambda, Bw)
    call Theta_wave(L, lambda, Theta_function)
    call D_local_n(alpha, E, L, nres, lambda, Bw, Theta_function, D_local, whamp_curve)

    arc_element = cos(lambda) * sqrt(1.0 + 3.0*sin(lambda)**2)

    D_integrad(1) = D_local(1) / cos(alpha_ne) &
                    * (tan(alpha)/tan(alpha_ne))**2 &
                    * arc_element

    D_integrad(2) = D_local(2) * tan(alpha) / sin(alpha_ne) &
                    * arc_element

    D_integrad(3) = D_local(3) * cos(alpha) &
                    * arc_element

    return

end subroutine int_D_loc