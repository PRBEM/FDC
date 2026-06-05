! ==============================================================================
! SUBROUTINE: B_field
! ==============================================================================
!
! PURPOSE:
!   Returns the local magnetic field magnitude B_local [T] at a given
!   L-shell and magnetic latitude, using the dipole field formula:
!
!     B(L, lambda) = B0 * f(lambda) / L^3
!
!   where f(lambda) = sqrt(4 - 3*cos^2(lambda)) / cos^6(lambda) is the
!   latitude variation factor (see f_function), and B0 is the equatorial
!   surface field strength set in module constants.
!
!   This implements equation (3.8) of Baumjohann & Treumann, Basic Space
!   Plasma Physics, split across f_function (latitude factor) and this
!   routine (L-shell scaling and absolute field strength).
!
! INPUTS:
!   L        - Magnetic L-shell [Earth radii, dimensionless]
!   lambda   - Magnetic latitude [radians]
!
! OUTPUTS:
!   Blocal   - Local magnetic field magnitude [T]
!
! MODULE DEPENDENCIES:
!   constants - B0 (equatorial surface field [T]), PT (planet name)
!
! CALLS:
!   f_function - Returns f(lambda) = B(lambda)/B_eq for a dipole field
!
! NOTES:
!   - For PT = 'EARTH', 'JUPITER', 'SATURN': uses the full dipole formula
!     B0 * f(lambda) / L^3, where B0 is the planet's equatorial surface
!     field (set in main based on PT).
!   - For PT = 'DIPOLE': omits the L^3 factor, returning B0 * f(lambda).
!     This means the field does not vary with L-shell, only with latitude.
!     This may be intentional for uniform-field test cases, or a bug.
!     Verify before using PT='DIPOLE' in production runs.
!   - KNOWN LIMITATION: f_function encodes the dipole geometry. If B_field
!     is modified to return non-dipolar field values (e.g. from satellite
!     measurements), Local_alpha must also be updated to use
!     B_field(L,lambda)/B_field(L,0) rather than f_function directly.
!
! ==============================================================================

subroutine B_field(L, lambda, Blocal)
use constants, only: B0, PT
implicit none

    real*8, intent(in)  :: L, lambda
    real*8, intent(out) :: Blocal

    real*8 :: f_value

    call f_function(lambda, f_value)

    select case(PT)
        case('JUPITER', 'EARTH', 'SATURN')
            Blocal = B0 * f_value / (L**3.0)
        case('DIPOLE')
            ! Note: L^3 factor omitted, so the field is independent of L-shell.
            Blocal = B0 * f_value
        case default
            write(*,*) 'ERROR: Invalid Planet = ', trim(PT)
            stop
    end select

    return

end subroutine B_field