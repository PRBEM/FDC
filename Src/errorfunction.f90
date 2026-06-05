! ==============================================================================
! SUBROUTINE: errorfunction
! ==============================================================================
!
! PURPOSE:
!   Computes the error function erf(x) for real argument x, defined as:
!
!     erf(x) = (2/sqrt(pi)) * integral[ exp(-t^2) dt ] from 0 to x
!
!   Two different series expansions are used depending on |x|:
!
!   |x| < 3.5 (small argument): Taylor series expansion
!     erf(x) = (2/sqrt(pi)) * x * exp(-x^2) * sum[ x^(2k) / prod(i+0.5) ]
!     Converges rapidly for small x; iterated until relative error < 1e-15
!     or 50 terms reached.
!
!   |x| >= 3.5 (large argument): asymptotic expansion of erfc(x)
!     erf(x) = 1 - erfc(x), where erfc(x) ~ exp(-x^2)/(|x|*sqrt(pi))
!     * sum[ (-1)^k * (k-0.5)!! / x^(2k) ]
!     Fixed 12-term expansion; sign corrected for negative x.
!
! INPUTS:
!   x    - Argument of erf(x) [dimensionless]
!
! OUTPUTS:
!   ERR  - erf(x) [dimensionless], in range (-1, 1)
!
! MODULE DEPENDENCIES:
!   constants - PI
!
! NOTES:
!   - Convergence tolerance for the small-argument series is 1e-15,
!     giving approximately 15 significant digits.
!   - The large-argument asymptotic series uses a fixed 12 terms, which
!     is accurate to machine precision for |x| >= 3.5.
!   - This subroutine is called by PowerSpectralDensity (GAUSSIAN case)
!     to compute the normalisation of the truncated Gaussian spectrum.
!
! ==============================================================================

subroutine errorfunction(x, ERR)
use constants, only: PI
implicit none

    real*8,  intent(in)  :: x
    real*8,  intent(out) :: ERR

    real*8  :: c0, eps, er, r, x2
    integer :: k

    eps = 1.0d-15
    x2  = x * x

    if (abs(x) < 3.5) then

        ! Small argument: Taylor series expansion
        er = 1.0
        r  = 1.0
        do k = 1, 50
            r  = r * x2 / (k + 0.5)
            er = er + r
            if (abs(r) <= abs(er)*eps) exit
        end do
        c0  = 2.0 / sqrt(PI) * x * exp(-x2)
        ERR = c0 * er

    else

        ! Large argument: asymptotic expansion of erfc(x)
        er = 1.0
        r  = 1.0
        do k = 1, 12
            r  = -r * (k - 0.5) / x2
            er = er + r
        end do
        c0  = exp(-x2) / (abs(x) * sqrt(PI))
        ERR = 1.0 - c0 * er
        if (x < 0.0) ERR = -ERR

    end if

    return

end subroutine errorfunction