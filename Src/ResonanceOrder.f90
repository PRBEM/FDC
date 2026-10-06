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
!   - Used only when NresType = 'MR12'. With NresType = 'CONSTANT', the
!     harmonic range comes from nres_min / nres_max in the input file.
!   - omega_m_per is taken directly from the input file, in whatever unit
!     FreqUnit specifies. The formula is dimensionally consistent only for
!     FreqUnit = 'RATIO', where omega_m_per is a fraction of omega_e_eq.
!     For 'HZ' it would have to be multiplied by 2*PI, and for the
!     data-driven models (2016HISS, the chorus fits, EMIC) omega_m_per is
!     not set at all, so the result is meaningless. Check FreqUnit before
!     using MR12 mode; this routine does not.
!   - Albert (2005), Sec. 5.1.2, gives a different and more rigorous route
!     to the same quantity: since V_max decreases with |n|, there is an |n|
!     beyond which no resonance exists at any theta, and that bound can be
!     computed from the refractive index. That would remove the dependence
!     on omega_m_per and the fixed 85 degrees, and cover the data-driven
!     spectra as well.
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

    ! Normalised electron momentum p = sqrt(gamma^2 - 1) in units of m_e*c
    p = sqrt((E + 1.0)**2 - 1.0)

    ! Mourenas & Ripoll [2012] estimate of the maximum resonance harmonic
    resonmax = int(p * sqrt(omega_m_per * omega_e_eq * omega_pe**2 &
                            / omega_e_eq**3) * sin(85.0*RAD)) + 1

    return

end subroutine resonanceorder