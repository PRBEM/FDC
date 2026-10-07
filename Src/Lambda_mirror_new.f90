! ==============================================================================
! SUBROUTINE: Lambda_mirror_new
! ==============================================================================
!
! PURPOSE:
!   Computes the mirror latitude (lambda_m) for an electron with equatorial
!   pitch angle alpha, in a dipole magnetic field.
!
!   The mirror point is where the local pitch angle reaches 90 degrees:
!     sin^2(alpha) = B_eq / B_m = cos^6(lambda_m) / sqrt(1 + 3*sin^2(lambda_m))
!
!   With x = cos^2(lambda_m) and s = sin(alpha), and using
!   sin^2(lambda) = 1 - x, this becomes
!     x^3 = s^2 * sqrt(4 - 3x)
!   and, after squaring,
!     x^6 + 3*s^4*x - 4*s^4 = 0
!   which is solved for x with the companion-matrix root finder zrhqr1.
!
!   The physically valid root satisfies 0 <= x <= 1, and
!     lambda_m = acos(sqrt(x)) * 0.999
!
! INPUTS:
!   alpha    - Equatorial electron pitch angle [radians]
!
! OUTPUTS:
!   lambda_m - Mirror latitude [radians], shortened by the factor 0.999
!
! CALLS:
!   zrhqr1   - Polynomial root finder (degree 6, via LAPACK dgeev)
!
! NOTES:
!   - Only real roots in [0, 1] are accepted. If several qualify, the last
!     one found is used; in practice there is exactly one physical root.
!   - THE 0.999 factor shortens the mirror latitude so that the local
!     diffusion coefficient is never evaluated exactly at the turning point,
!     where the local pitch angle is 90 degrees, v_parallel = 0, and
!     quantities such as 1/(beta*cos(alpha_local)*cos(theta)) diverge. It is
!     a property of the integration endpoint only.
!
! ==============================================================================

subroutine Lambda_mirror_new(alpha, lambda_m)
implicit none

    real*8, intent(in)  :: alpha
    real*8, intent(out) :: lambda_m

    real*8  :: s4             ! sin^4(alpha)
    real*8  :: P(7)           ! polynomial coefficients, ascending powers
    real*8  :: rtr(6), rti(6) ! real and imaginary parts of the roots
    real*8  :: root           ! x = cos^2(lambda_m)
    integer :: i

    real*8, parameter :: TOL_IMAG = 1.0d-8

    s4 = sin(alpha)**4

    ! x^6 + 3*s^4*x - 4*s^4 = 0
    P(7) =  1.0
    P(6) =  0.0
    P(5) =  0.0
    P(4) =  0.0
    P(3) =  0.0
    P(2) =  3.0 * s4
    P(1) = -4.0 * s4

    call zrhqr1(P, rtr, rti)

    ! physically valid root: real and in [0, 1]
    root = 0.0
    do i = 1, 6
        if (abs(rti(i)) < TOL_IMAG) then
            if (rtr(i) >= 0.0 .and. rtr(i) <= 1.0) root = rtr(i)
        end if
    end do

    if (root == 0.0) then
        write(*,*) 'WARNING: Lambda_mirror_new found no valid root for alpha=', alpha
    end if

    ! shortened by 0.999: integration endpoint only, see the notes above
    lambda_m = acos(sqrt(root)) * 0.999

    return

end subroutine Lambda_mirror_new