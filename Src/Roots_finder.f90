! ==============================================================================
! SUBROUTINE: Roots_finder
! ==============================================================================
!
! PURPOSE:
!   Finds the resonant wave frequencies (roots of the combined resonance
!   condition and cold plasma dispersion relation) for a given electron
!   state (alpha, E), field line location (L, lambda), resonance harmonic
!   (nres), and wave normal angle (theta).
!
!   The resonance condition is:
!     omega - k_par * v_par - n * Omega_e / gamma = 0
!   Combined with the cold plasma dispersion relation, this yields a
!   polynomial in omega whose real roots within the wave frequency band
!   [omega_lc, omega_uc] are the resonant frequencies.
!
!   For nres /= 0: a 10th-degree polynomial in omega is solved (general case).
!   For nres  = 0: symmetry reduces the problem to a cubic in omega^2,
!                  solved as a 3rd-degree polynomial.
!
!   Polynomial roots are found using LAPACK via zrhqr2 (nres /= 0) or
!   zrhqr3 (nres = 0).
!
! INPUTS:
!   alpha   - Equatorial electron pitch angle [radians]
!   E       - Electron kinetic energy normalised to rest-mass energy
!             (E = gamma - 1, dimensionless)
!   L       - Magnetic L-shell [Earth radii, dimensionless]
!   lambda  - Magnetic latitude [radians]
!   nres    - Resonance harmonic number n (integer, e.g. -1, 0, 1, ...)
!   theta   - Wave normal angle [radians]
!
! OUTPUTS:
!   roots(4) - Up to 4 resonant frequencies [rad/s]. Unused slots are set
!              to 0.0. Valid roots are real, lie within [omega_lc, omega_uc],
!              and satisfy the group velocity condition |b2| > 1.
!
! MODULE DEPENDENCIES:
!   constants  - sqe, and indirectly via called routines
!   variables  - READS:  omega_e, omega_pe, omega_i, omega_pi,
!                        omega_lc, omega_uc (set by Freq_Spectral/WeightFunction)
!                WRITES: a1, a2, a3, a4, a5, a6 (polynomial coefficients,
!                        also used by WeightFunction after this call returns)
!
! CALLS:
!   Local_alpha  - Converts equatorial pitch angle to local pitch angle
!   zrhqr2       - LAPACK-based polynomial root finder (degree 10, nres /= 0)
!   zrhqr3       - LAPACK-based polynomial root finder (degree 3,  nres = 0)
!
! SIDE EFFECTS:
!   Writes a1-a6 to module variables. These are combinations of the plasma
!   frequencies used both here (as polynomial coefficients a7-a16) and
!   subsequently in WeightFunction for the analytic dS/dw, dR/dw, dL/dw
!   dispersion derivatives. WeightFunction silently depends on Roots_finder
!   having been called first to populate a1-a6. This implicit dependency
!   is a known issue to be addressed in a future refactoring step.
!
!   Definitions:
!     a1 = omega_i^2 + omega_e^2 + omega_pe^2 + omega_pi^2
!     a2 = (omega_i*omega_e)^2 + (omega_pe*omega_i)^2 + (omega_pi*omega_e)^2
!     a3 = omega_i + sqe*omega_e
!     a4 = omega_i*sqe*omega_e - omega_pe^2 - omega_pi^2
!     a5 = 0.0  (term analytically zero; retained for structural completeness)
!     a6 = omega_pi^2 + omega_pe^2
!
! NOTES:
!   - a5 is hardcoded to 0.0. The original commented-out definition was
!     -(omega_pi^2 * sqe*omega_e + omega_pe^2 * omega_i). Whether the zero
!     value is a valid simplification or a bug should be verified against
!     Albert (2007).
!   - For nres = 0, only even powers of omega appear, so the polynomial is
!     solved in omega^2 and roots are recovered by taking the square root.
!     Coefficients are taken from alternating entries of the general coeff()
!     array (indices 5, 7, 9, 11).
!   - Root validity checks:
!       nres /= 0: root must be real (rti = 0), lie in [omega_lc, omega_uc],
!                  and satisfy |(omega - b1)*b2/omega| > 1
!       nres  = 0: root^2 must be real, lie in [omega_lc^2, omega_uc^2],
!                  and |b2| > 1
!
! ==============================================================================

subroutine Roots_finder(alpha, E, L, lambda, nres, theta, roots)
use constants
use variables
implicit none

    real*8, intent(in)  :: alpha, E, L, lambda, theta
    integer, intent(in) :: nres
    real*8, intent(out) :: roots(4)

    real*8 :: coeff(11)
    real*8 :: rtr(10), rti(10)

    real*8 :: coeff1(4)
    real*8 :: rtr1(3), rti1(3)

    ! Plasma frequency combination coefficients
    real*8 :: a1, a2, a3, a4, a5, a6

    ! Intermediate polynomial coefficients (local to this routine)
    real*8 :: a7, a8, a9, a10, a11, a12, a13, a14, a15, a16

    ! Resonance condition parameters
    real*8 :: b1, b2

    ! Particle state
    real*8 :: gamma1, betta, alpha_ne

    ! Wave normal angle
    real*8 :: X

    integer :: i, NN
    real*8  :: Temp2

    betta  = sqrt(E*(E + 2.0)) / (E + 1.0)
    gamma1 = E + 1.0

    call Local_alpha(alpha, lambda, alpha_ne)

    X = tan(theta)

    ! -------------------------------------------------------------------------
    ! Plasma frequency combination coefficients
    ! Written to module variables (a1-a6) for subsequent use in WeightFunction.
    ! Also used locally below as inputs to the polynomial coefficients a7-a16.
    ! -------------------------------------------------------------------------
    a1 = omega_i**2.0 + omega_e**2.0 + omega_pe**2.0 + omega_pi**2.0
    a2 = (omega_i*omega_e)**2.0 + (omega_pe*omega_i)**2.0 + (omega_pi*omega_e)**2.0
    a3 = omega_i + real(sqe)*omega_e
    a4 = omega_i*real(sqe)*omega_e - omega_pe**2.0 - omega_pi**2.0
    a5 = 0.0   ! analytically zero; see note in header regarding Albert (2007)
    a6 = omega_pi**2.0 + omega_pe**2.0

    ! Intermediate polynomial coefficients (local only)
    a7  = (X**2.0)*a1 + omega_i**2.0 + omega_e**2.0 + a6
    a8  = (X**2.0)*a2 + (omega_i*omega_e)**2.0 + a6*(omega_i**2.0 + omega_e**2.0)
    a9  = a6*((omega_i*omega_e)**2.0)
    a10 = 2.0*a4 - a3**2.0 - a6
    a11 = a4**2.0 - 2.0*a3*a5 - 2.0*a4*a6 + (a3**2.0)*a6
    a12 = a5**2.0 + (a4**2.0)*a6 - 2.0*a3*a5*a6
    a13 = (a5**2.0)*a6
    a14 = (X**2.0)*(2.0*a4 - a3**2.0) - (2.0 + X**2.0)*(a1 + a6)
    a15 = (X**2.0)*(a4**2.0 - 2.0*a3*a5) + (2.0 + X**2.0)*(a2 + a1*a6)
    a16 = (X**2.0)*(a5**2.0) + (2.0 + X**2.0)*a2*a6

    ! Resonance condition parameters
    b1 = real(nres)*real(sqe)*omega_e / gamma1
    b2 = 1.0 / (betta*cos(alpha_ne)*cos(theta))

    ! -------------------------------------------------------------------------
    ! Polynomial coefficients
    ! -------------------------------------------------------------------------
    coeff(11) = (1 + X**2.0)*((1 - b2**2.0)**2.0)
    coeff(10) = 4.0*(1 + X**2.0)*b1*(b2**2.0)*(1 - b2**2.0)
    coeff(9)  = (6.0*(1 + X**2.0)*(b1**2.0) - a7)*(b2**4.0) &
                - (a14 + 2.0*(1 + X**2.0)*(b1**2.0))*(b2**2.0) &
                + (1 + X**2.0)*a10
    coeff(8)  = (4.0*a7*b1 - 4.0*(1 + X**2.0)*(b1**3.0))*(b2**4.0) &
                + 2.0*a14*b1*(b2**2.0)
    coeff(7)  = (a8 - 6.0*a7*(b1**2.0) + (1 + X**2.0)*(b1**4.0))*(b2**4.0) &
                - (a15 + a14*(b1**2.0))*(b2**2.0) &
                + (1 + X**2.0)*a11
    coeff(6)  = (4.0*a7*(b1**3.0) - 4.0*a8*b1)*(b2**4.0) &
                + 2.0*a15*b1*(b2**2.0)
    coeff(5)  = (6.0*a8*(b1**2.0) - a7*(b1**4.0) - a9)*(b2**4.0) &
                - (a15*(b1**2.0) - a16)*(b2**2.0) &
                - (1 + X**2.0)*a12
    coeff(4)  = (4.0*a9*b1 - 4.0*a8*(b1**3.0))*(b2**4.0) &
                - 2.0*a16*b1*(b2**2.0)
    coeff(3)  = (a8*(b1**4.0) - 6.0*a9*(b1**2.0))*(b2**4.0) &
                + a16*(b1**2.0)*(b2**2.0) &
                + (1 + X**2.0)*a13
    coeff(2)  = 4.0*a9*(b1**3.0)*(b2**4.0)
    coeff(1)  = -a9*(b1**4.0)*(b2**4.0)

    ! -------------------------------------------------------------------------
    ! Root finding
    ! -------------------------------------------------------------------------
    roots(1) = 0.0
    roots(2) = 0.0
    roots(3) = 0.0
    roots(4) = 0.0
    NN = 0

    if (nres /= 0) then

        call zrhqr2(coeff, rtr, rti)

        do i = 1, 10
            if (rti(i) == 0) then
                if ((rtr(i) > omega_lc) .and. (rtr(i) < omega_uc)) then
                    Temp2 = abs((rtr(i) - b1)*b2 / rtr(i))
                    if (Temp2 > 1.0) then
                        NN = NN + 1
                        roots(NN) = rtr(i)
                    end if
                end if
            end if
        end do

    else

        ! nres = 0: polynomial is even in omega; solve as cubic in omega^2
        coeff1(4) = coeff(11)
        coeff1(3) = coeff(9)
        coeff1(2) = coeff(7)
        coeff1(1) = coeff(5)

        call zrhqr3(coeff1, rtr1, rti1)

        do i = 1, 3
            if (rti1(i) == 0) then
                if ((rtr1(i) > (omega_lc**2.0)) .and. &
                    (rtr1(i) < (omega_uc**2.0)) .and. &
                    (abs(b2) > 1.0)) then
                    NN = NN + 1
                    roots(NN) = sqrt(rtr1(i))
                end if
            end if
        end do

    end if

    return

end subroutine Roots_finder