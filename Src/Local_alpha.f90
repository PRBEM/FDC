! ==============================================================================
! SUBROUTINE: Local_alpha
! ==============================================================================
!
! PURPOSE:
!   Converts the equatorial electron pitch angle (alpha) to the local pitch
!   angle (alpha_ne) at a given magnetic latitude (lambda), using the
!   conservation of the first adiabatic invariant (magnetic moment):
!
!     sin^2(alpha_ne) / B_local = sin^2(alpha) / B_eq
!
!   For a dipole field this gives:
!     sin(alpha_ne) = sin(alpha) * sqrt(f(lambda))
!
!   where f(lambda) = B_local / B_eq is the ratio of the local to equatorial
!   magnetic field strength, computed by f_function.
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
!   - The argument of asin is sqrt(f(lambda))*sin(alpha). For this to be
!     valid, sqrt(f(lambda))*sin(alpha) must not exceed 1.0, i.e. the
!     particle must not have reached its mirror point. The calling routine
!     (Integration_lambda) is responsible for ensuring lambda does not
!     exceed the mirror latitude.
!
! ==============================================================================

subroutine Local_alpha(alpha, lambda, alpha_ne)
implicit none

    real*8, intent(in)  :: alpha, lambda
    real*8, intent(out) :: alpha_ne

    real*8 :: f_value, temp

    call f_function(lambda, f_value)

    temp     = sqrt(f_value) * sin(alpha)
    alpha_ne = asin(temp)

    return

end subroutine Local_alpha