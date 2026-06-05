! ==============================================================================
! SUBROUTINE: f_function
! ==============================================================================
!
! PURPOSE:
!   Computes the ratio of the local magnetic field strength to the equatorial
!   field strength at a given magnetic latitude lambda, for a dipole field:
!
!     f(lambda) = B(lambda) / B_eq
!               = sqrt(4 - 3*cos^2(lambda)) / cos^6(lambda)
!
!   This follows from the dipole field magnitude (Baumjohann & Treumann,
!   Basic Space Plasma Physics, eq. 3.8):
!     B(lambda, L) = B_E/L^3 * sqrt(1 + 3*sin^2(lambda)) / cos^6(lambda)
!   where B_E is the equatorial surface field. The L^3 and B_E factors are
!   handled by B_field; f_function returns only the latitude-dependent part.
!
!   f(lambda) appears in:
!     - Local_alpha: sin(alpha_ne) = sin(alpha) * sqrt(f(lambda))
!     - B_field:     B_local = B0 * f(lambda) / L^3
!
!   Reference: Baumjohann & Treumann, Basic Space Plasma Physics,
!              equation (3.9)
!
! INPUTS:
!   lambda   - Magnetic latitude [radians]
!
! OUTPUTS:
!   f_value  - B(lambda)/B_eq, dimensionless ratio [-]
!
! NOTES:
!   - f(0) = 1 at the equator (lambda = 0), as expected.
!   - f(lambda) diverges as lambda -> 90 degrees (cos(lambda) -> 0),
!     corresponding to the magnetic pole. The calling routines are
!     responsible for ensuring lambda stays within physically valid bounds.
!
! ==============================================================================
 
subroutine f_function(lambda, f_value)
implicit none
 
    real*8, intent(in)  :: lambda
    real*8, intent(out) :: f_value
 
    f_value = sqrt(4.0 - 3.0*(cos(lambda))**2.0) / (cos(lambda))**6.0
 
    return
 
end subroutine f_function
