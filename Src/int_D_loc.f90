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
!   pitch angle frames.
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
!   lambda  - Magnetic latitude at this evaluation point [radians]
!   nres    - Resonance harmonic number n (integer, e.g. -1, 0, 1, ...)
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
!   (none explicitly declared; constants and variables accessed via called
!    routines)
!
! NOTES:
!   - Bw and Theta_function are obtained internally by calling Bw_wave and
!     Theta_wave; they are not passed in from Integration_lambda. This means
!     wave parameters are re-fetched at every lambda point.
!
! AUTHOR: [original author]
! MODIFIED: [date] - Step 5 refactoring: removed redundant bare integer
!                    declaration of nres and old commented-out geometric
!                    factors block
! ==============================================================================

subroutine int_D_loc(alpha, E, L, nres, lambda, D_integrad)
implicit none
    real*8, intent(in)  :: alpha, E, L, lambda
    integer, intent(in) :: nres
    real*8, intent(out) :: D_integrad(3)

    real*8  :: Bw
    real*8  :: Theta_function(7)
    real*8  :: alpha_ne
    real*8  :: D_local(3)

    ! Arc length element ds/dlambda for dipole field
    real*8  :: arc_element

    call Local_alpha(alpha, lambda, alpha_ne)
    call Bw_wave(L, lambda, Bw)
    call Theta_wave(L, lambda, Theta_function)
    call D_local_n(alpha, E, L, nres, lambda, Bw, Theta_function, D_local)

    arc_element = cos(lambda) * sqrt(1.0 + 3.0*(sin(lambda))**2.0)

    D_integrad(1) = D_local(1) / cos(alpha_ne) &
                    * (tan(alpha)/tan(alpha_ne))**2.0 &
                    * arc_element

    D_integrad(2) = D_local(2) * tan(alpha) / sin(alpha_ne) &
                    * arc_element

    D_integrad(3) = D_local(3) * cos(alpha) &
                    * arc_element

    return

end subroutine int_D_loc