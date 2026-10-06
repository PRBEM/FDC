! ==============================================================================
! SUBROUTINE: Local_alpha
! ==============================================================================
!
! PURPOSE:
!   Converts the equatorial electron pitch angle (alpha) to the local pitch
!   angle (alpha_ne) at a given magnetic latitude (lambda), using
!   conservation of the first adiabatic invariant (magnetic moment):
!
!     sin^2(alpha_ne) / B_local = sin^2(alpha) / B_eq
!
!   For a dipole field this gives:
!     sin(alpha_ne) = sin(alpha) * sqrt(f(lambda))
!
!   where f(lambda) = B_local / B_eq is computed by f_function.
!
! INPUTS:
!   alpha    - Equatorial electron pitch angle [radians]
!   lambda   - Magnetic latitude [radians]
!
! OUTPUTS:
!   alpha_ne - Local electron pitch angle at lambda [radians]
!
! CALLS:
!   f_function - Returns f(lambda) = B_local/B_eq for a dipole field
!
! NOTES:
!   - The argument of the asin, sqrt(f(lambda))*sin(alpha), equals 1 exactly
!     at the mirror point and would exceed 1 beyond it. It is clamped to 1
!     so that rounding at the mirror point cannot produce a NaN; the callers
!     (Integration_lambda via int_D_loc, and D_ba) keep lambda below the
!     mirror latitude, which Lambda_mirror_new additionally shortens by
!     0.999.
!   - This routine assumes the dipole geometry encoded in f_function. If
!     B_field is ever changed to a non-dipolar model, alpha_ne must be
!     computed as asin(sin(alpha)*sqrt(B_field(L,lambda)/B_field(L,0)))
!     instead (see the note in B_field).
!
! ==============================================================================

subroutine Local_alpha(alpha, lambda, alpha_ne)
implicit none

    real*8, intent(in)  :: alpha, lambda
    real*8, intent(out) :: alpha_ne

    real*8 :: f_value, sin_alpha_ne

    call f_function(lambda, f_value)

    sin_alpha_ne = sqrt(f_value) * sin(alpha)
    alpha_ne     = asin(min(1.0d0, sin_alpha_ne))

    return

end subroutine Local_alpha