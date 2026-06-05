! ==============================================================================
! SUBROUTINE: resonanceorder
! ==============================================================================
!
! PURPOSE:
!   Estimates the maximum resonance harmonic number (resonmax) needed for
!   convergence of the bounce-averaged diffusion coefficients, following
!   the empirical formula of Mourenas & Ripoll (2012).
!
!   The formula estimates the highest harmonic n for which resonance can
!   occur, based on the electron momentum, wave frequency, and local plasma
!   parameters at the equator:
!
!     resonmax = int( p * sqrt(omega_m * omega_pe^2 / omega_e^3)
!                      * sin(85 deg) ) + 1
!
!   where:
!     p        = normalised electron momentum = sqrt((E+1)^2 - 1) [m_e * c]
!     omega_m  = mean wave frequency [rad/s]
!     omega_pe = equatorial electron plasma frequency [rad/s]
!     omega_e  = equatorial electron gyrofrequency [rad/s]
!     sin(85 deg) ~ 0.996 represents the maximum pitch angle considered
!
!   All plasma parameters are evaluated at the equator (lambda = 0).
!
!   Reference: Mourenas, D., and J.-F. Ripoll (2012), Analytical estimates
!   of quasi-linear diffusion coefficients and electron lifetimes in the
!   inner radiation belt, J. Geophys. Res., 117, A01204,
!   doi:10.1029/2011JA016985
!
! INPUTS:
!   E  - Electron kinetic energy normalised to rest-mass energy
!        (E = gamma - 1, dimensionless)
!   L  - Magnetic L-shell [Earth radii, dimensionless]
!
! OUTPUTS:
!   resonmax - Maximum resonance harmonic number (positive integer)
!              The summation in sum_n runs from -resonmax to +resonmax.
!
! CALLS:
!   B_field         - Returns equatorial magnetic field at (L, lambda=0)
!   ElectronDensity - Returns equatorial electron density at (L, lambda=0)
!
! MODULE DEPENDENCIES:
!   constants - qe, me, C0, RAD, omega_m_per
!
! NOTES:
!   - omega_m_per is used directly in the formula. If FreqUnit='HZ',
!     omega_m_per is in Hz and should be multiplied by 2*PI before use.
!     If FreqUnit='RATIO', omega_m_per is a fraction of omega_e_eq and
!     the formula is dimensionally consistent. For other FreqUnit types
!     (data-driven models), omega_m_per may not be set and this formula
!     may give incorrect results. Verify FreqUnit before using MR12 mode.
!   - The +1 ensures resonmax >= 1 even for very low energy electrons.
!
! ==============================================================================

subroutine resonanceorder(E, L, resonmax)
use constants
implicit none

    real*8,  intent(in)  :: E, L
    integer, intent(out) :: resonmax

    real*8 :: Ne_local, p, B_eq, omega_pe, omega_e_eq

    ! Equatorial plasma parameters
    call B_field(L, 0.0d0, B_eq)
    call ElectronDensity(L, 0.0d0, Ne_local)

    omega_e_eq = qe * B_eq / me
    omega_pe   = sqrt(Ne_local * qe * qe / (C0 * me))

    ! Normalised electron momentum p = sqrt((gamma)^2 - 1) in units of m_e*c
    p = sqrt((E + 1.0)**2.0 - 1.0)

    ! Mourenas & Ripoll [2012] estimate of maximum resonance harmonic
    resonmax = int(p * sqrt(omega_m_per * omega_e_eq * omega_pe**2.0 &
                            / omega_e_eq**3.0) * sin(85.0*RAD)) + 1

    return

end subroutine resonanceorder