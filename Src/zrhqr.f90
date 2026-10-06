! ==============================================================================
! SUBROUTINES: zrhqr_roots, zrhqr1, zrhqr2, zrhqr3
! ==============================================================================
!
! PURPOSE:
!   Find all roots of a real polynomial using LAPACK's dgeev eigenvalue
!   solver. The polynomial is converted to companion matrix form, whose
!   eigenvalues are exactly the roots of the polynomial.
!
!   zrhqr_roots does the work for any degree; zrhqr1, zrhqr2 and zrhqr3 are
!   thin wrappers kept for the existing call sites:
!     zrhqr1 - degree 6  (Lambda_mirror_new)
!     zrhqr2 - degree 10 (Roots_finder, nres /= 0, fallback and -DROOTS_POLY)
!     zrhqr3 - degree 3  (Roots_finder, nres = 0,  fallback and -DROOTS_POLY)
!
! ALGORITHM:
!   1. Construct the companion (Frobenius) matrix from the polynomial
!      coefficients A(1..N+1), where A(N+1) is the leading coefficient:
!        hess(1,j) = -A(N+1-j) / A(N+1)   for j = 1..N
!        hess(k+1,k) = 1                    for k = 1..N-1
!        all other entries = 0
!   2. Compute all eigenvalues of hess via LAPACK dgeev.
!   3. Sort eigenvalues by ascending real part using insertion sort.
!
! INPUTS:
!   N       - Polynomial degree
!   A(N+1)  - Polynomial coefficients in ascending order:
!             A(1) + A(2)*x + ... + A(N+1)*x^N
!             A(N+1) must be non-zero (leading coefficient).
!
! OUTPUTS:
!   WR(N)   - Real parts of the roots, sorted by ascending real part
!   WI(N)   - Imaginary parts of the roots (corresponding to WR)
!
! CALLS:
!   dgeev   - LAPACK double precision eigenvalue solver
!
! NOTES:
!   - Callers rely on WI being EXACTLY zero for real roots, which is how
!     dgeev reports eigenvalues it resolved as 1x1 blocks of the real Schur
!     form. Any replacement must preserve that property.
!   - dgeev is called with JOBVL = JOBVR = 'N', so no eigenvectors are
!     computed and VL/VR are not referenced; they are passed as 1x1 dummies.
!   - dgeev balances the matrix before reducing it to Hessenberg form. The
!     companion matrix is already Hessenberg, so the reduction is redundant,
!     but balancing is not: the coefficients span many orders of magnitude
!     (omega in rad/s to the 10th power), and skipping it loses accuracy.
!     Calling dhseqr directly without balancing was measured to change the
!     number of exactly-real roots in about 6% of cases; dgebal + dhseqr is
!     accurate but only ~10% faster, so dgeev is kept.
!   - If dgeev fails (info /= 0), a one-time warning is issued; the first
!     info-1 eigenvalues it did compute are still returned.
!
! ==============================================================================


! ------------------------------------------------------------------------------
! zrhqr_roots: roots of a real polynomial of degree N
! ------------------------------------------------------------------------------
subroutine zrhqr_roots(N, A, WR, WI)
implicit none

    integer, intent(in)  :: N
    real*8,  intent(in)  :: A(N+1)
    real*8,  intent(out) :: WR(N), WI(N)

    real*8  :: hess(N,N)
    real*8  :: WORK(4*N)
    real*8  :: VL(1,1), VR(1,1)      ! not referenced for JOBVL = JOBVR = 'N'
    real*8  :: xr, xi
    integer :: info, j, k

    logical, save :: warned = .false.

    ! Companion matrix
    hess(1,:)   = -A(N:1:-1) / A(N+1)
    hess(2:N,:) = 0.0
    do k = 1, N-1
        hess(k+1,k) = 1.0
    end do

    call dgeev('N', 'N', N, hess, N, WR, WI, VL, 1, VR, 1, WORK, 4*N, info)

    if (info /= 0 .and. .not. warned) then
        write(*,*) 'WARNING: dgeev returned info = ', info, ' for a degree-', N, &
                   ' polynomial; roots 1..', info-1, ' are valid'
        warned = .true.
    end if

    ! Sort roots by ascending real part (insertion sort)
    do j = 2, N
        xr = WR(j); xi = WI(j)
        k = j - 1
        do while (k >= 1)
            if (WR(k) <= xr) exit
            WR(k+1) = WR(k); WI(k+1) = WI(k)
            k = k - 1
        end do
        WR(k+1) = xr; WI(k+1) = xi
    end do

    return
end subroutine zrhqr_roots


! ------------------------------------------------------------------------------
! zrhqr1: degree 6   (Lambda_mirror_new)
! ------------------------------------------------------------------------------
subroutine zrhqr1(A, WR, WI)
implicit none
    real*8, intent(in)  :: A(7)
    real*8, intent(out) :: WR(6), WI(6)
    call zrhqr_roots(6, A, WR, WI)
end subroutine zrhqr1


! ------------------------------------------------------------------------------
! zrhqr2: degree 10  (Roots_finder, nres /= 0)
! ------------------------------------------------------------------------------
subroutine zrhqr2(A, WR, WI)
implicit none
    real*8, intent(in)  :: A(11)
    real*8, intent(out) :: WR(10), WI(10)
    call zrhqr_roots(10, A, WR, WI)
end subroutine zrhqr2


! ------------------------------------------------------------------------------
! zrhqr3: degree 3   (Roots_finder, nres = 0)
! ------------------------------------------------------------------------------
subroutine zrhqr3(A, WR, WI)
implicit none
    real*8, intent(in)  :: A(4)
    real*8, intent(out) :: WR(3), WI(3)
    call zrhqr_roots(3, A, WR, WI)
end subroutine zrhqr3