! ==============================================================================
! SUBROUTINE: Lambda_mirror_new
! ==============================================================================
!
! PURPOSE:
!   Computes the mirror latitude (lambda_m) for an electron with equatorial
!   pitch angle alpha, in a dipole magnetic field.
!
!   The mirror point is defined by the condition that the local pitch angle
!   reaches 90 degrees, i.e.:
!     sin^2(alpha) = B_eq / B_mirror = cos^6(lambda_m) / f(lambda_m)
!
!   For a dipole field, f(lambda) = (1 + 3*sin^2(lambda))^0.5 / cos^6(lambda),
!   so the mirror condition reduces to:
!     cos^6(lambda_m) = sin^2(alpha) * (1 + 3*sin^2(lambda_m))^0.5
!
!   This is solved as a 6th-degree polynomial in cos^2(lambda_m) = x:
!     x^3 - 3*sin^4(alpha)*x - 4*sin^4(alpha) = 0   (... wait -- x^6 form)
!   using the LAPACK-based root finder zrhqr1.
!
!   The physically valid root satisfies 0 <= x <= 1 (since cos^2(lambda_m)
!   must be in [0,1]), and the mirror latitude is:
!     lambda_m = acos(sqrt(x)) * 0.999
!
!   The factor 0.999 keeps the integration endpoint just below the exact
!   mirror point, avoiding the singularity where sin(alpha_ne) = 1 in
!   Local_alpha.
!
! INPUTS:
!   alpha    - Equatorial electron pitch angle [radians]
!
! OUTPUTS:
!   lambda_m - Mirror latitude [radians], reduced by factor 0.999 to avoid
!              the mirror point singularity
!
! CALLS:
!   zrhqr1   - LAPACK-based polynomial root finder (degree 6)
!
! NOTES:
!   - The polynomial coefficients correspond to:
!       P(7)*x^6 + P(2)*x + P(1) = 0
!     where P(7)=1, P(2)=3*sin^4(alpha), P(1)=-4*sin^4(alpha),
!     and all other coefficients are zero.
!   - Only real roots in [0, 1] are accepted. If multiple such roots exist,
!     the last one found is used (in practice only one physical root exists).
!   - The 0.999 reduction factor is applied to lambda_m (not to x), so the
!     integration upper bound is 99.9% of the mirror latitude in radians.
!
! ==============================================================================

subroutine Lambda_mirror_new(alpha, lambda_m)
implicit none

    real*8, intent(in)  :: alpha
    real*8, intent(out) :: lambda_m

    real*8  :: BB
    real*8  :: P(7)
    real*8  :: rtr(6), rti(6)
    real*8  :: root
    integer :: i

    BB = sin(alpha)

    ! Polynomial coefficients for cos^6(lambda_m) = sin^2(alpha)*(1+3sin^2(lambda_m))^0.5
    ! expressed as a degree-6 polynomial in x = cos^2(lambda_m)
    P(7) =  1.0
    P(6) =  0.0
    P(5) =  0.0
    P(4) =  0.0
    P(3) =  0.0
    P(2) =  3.0 * (BB**4.0)
    P(1) = -4.0 * (BB**4.0)

    call zrhqr1(P, rtr, rti)

    ! Find the physically valid root: real, in [0, 1]
    root = 0.0
    do i = 1, 6
        if (rti(i) == 0.0) then
            if (rtr(i) >= 0.0 .and. rtr(i) <= 1.0) then
                root = rtr(i)
            end if
        end if
    end do

    ! Mirror latitude, reduced slightly to avoid the mirror point singularity
    lambda_m = acos(sqrt(root)) * 0.999

    return

end subroutine Lambda_mirror_new
