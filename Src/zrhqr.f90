! ==============================================================================
! SUBROUTINES: zrhqr1, zrhqr2, zrhqr3
! ==============================================================================
!
! PURPOSE:
!   Find all roots of a real polynomial using LAPACK's dgeev eigenvalue
!   solver. The polynomial is converted to companion matrix form, whose
!   eigenvalues are exactly the roots of the polynomial.
!
!   Three variants exist for different polynomial degrees:
!     zrhqr1 - degree 6  (used by Lambda_mirror_new)
!     zrhqr2 - degree 10 (used by Roots_finder, nres /= 0 case)
!     zrhqr3 - degree 3  (used by Roots_finder, nres = 0 case)
!
!   All three are identical in structure — only the polynomial degree N2
!   and the corresponding array sizes differ. This duplication is a known
!   issue; a future refactoring could replace all three with a single
!   general subroutine taking N2 as an argument.
!
! ALGORITHM:
!   1. Construct the companion (Frobenius) matrix from the polynomial
!      coefficients A(1..N2+1), where A(N2+1) is the leading coefficient:
!        hess(1,j) = -A(N2+1-j) / A(N2+1)   for j = 1..N2
!        hess(k+1,k) = 1                      for k = 1..N2-1
!        all other entries = 0
!   2. Compute all eigenvalues of hess via LAPACK dgeev.
!   3. Sort eigenvalues by ascending real part using insertion sort.
!
! INPUTS:
!   A(N2+1) - Polynomial coefficients in ascending order:
!             A(1) + A(2)*x + ... + A(N2+1)*x^N2
!             A(N2+1) must be non-zero (leading coefficient).
!
! OUTPUTS:
!   WR(N2)  - Real parts of the roots, sorted by ascending real part
!   WI(N2)  - Imaginary parts of the roots (corresponding to WR)
!
! CALLS:
!   dgeev   - LAPACK double precision eigenvalue solver
!
! NOTES:
!   - WORK is declared as integer but dgeev expects a real*8 workspace.
!     This type mismatch is silently accepted by most compilers but is
!     technically non-conforming. A future fix should declare WORK as
!     real*8.
!   - The LWORK values (30, 45, 15) are set larger than the minimum
!     required (4*N2) to provide workspace margin. The optimal value
!     can be queried by calling dgeev with LWORK=-1.
!
! ==============================================================================


! ------------------------------------------------------------------------------
! zrhqr1: polynomial root finder, degree 6
! Used by: Lambda_mirror_new
! ------------------------------------------------------------------------------
subroutine zrhqr1(A, WR, WI)
implicit none
 
    integer, parameter :: N2    = 6
    integer, parameter :: LDA   = 6
    integer, parameter :: LDVL  = 6
    integer, parameter :: LDVR  = 6
    integer, parameter :: LWORK = 30
 
    real*8,  intent(in)  :: A(N2+1)
    real*8,  intent(out) :: WR(N2), WI(N2)
 
    real*8  :: hess(N2,N2), xr, xi
    real*8  :: VL(LDVL,N2), VR(LDVR,N2)
    integer :: WORK(LWORK)
    integer :: info, j, k
 
    WORK(:) = 0
 
    ! Construct companion matrix
    hess(1,:)    = -A(N2:1:-1) / A(N2+1)
    hess(2:N2,:) = 0.0
    do k = 1, N2-1
        hess(k+1,k) = 1.0
    end do
 
    call dgeev('N', 'N', N2, hess, LDA, WR, WI, VL, LDVL, VR, LDVR, &
               WORK, LWORK, INFO)
 
    ! Sort roots by ascending real part (insertion sort)
    do j = 2, N2
        xr = WR(j); xi = WI(j)
        k = j - 1
        do while (k >= 1 .and. WR(k) > xr)
            WR(k+1) = WR(k); WI(k+1) = WI(k)
            k = k - 1
        end do
        WR(k+1) = xr; WI(k+1) = xi
    end do
 
    return
end subroutine zrhqr1
 
 
! ------------------------------------------------------------------------------
! zrhqr2: polynomial root finder, degree 10
! Used by: Roots_finder (nres /= 0 case)
! ------------------------------------------------------------------------------
subroutine zrhqr2(A, WR, WI)
implicit none
 
    integer, parameter :: N2    = 10
    integer, parameter :: LDA   = 10
    integer, parameter :: LDVL  = 10
    integer, parameter :: LDVR  = 10
    integer, parameter :: LWORK = 45
 
    real*8,  intent(in)  :: A(N2+1)
    real*8,  intent(out) :: WR(N2), WI(N2)
 
    real*8  :: hess(N2,N2), xr, xi
    real*8  :: VL(LDVL,N2), VR(LDVR,N2)
    integer :: WORK(LWORK)
    integer :: info, j, k
 
    WORK(:) = 0
 
    ! Construct companion matrix
    hess(1,:)    = -A(N2:1:-1) / A(N2+1)
    hess(2:N2,:) = 0.0
    do k = 1, N2-1
        hess(k+1,k) = 1.0
    end do
 
    call dgeev('N', 'N', N2, hess, LDA, WR, WI, VL, LDVL, VR, LDVR, &
               WORK, LWORK, INFO)
 
    ! Sort roots by ascending real part (insertion sort)
    do j = 2, N2
        xr = WR(j); xi = WI(j)
        k = j - 1
        do while (k >= 1 .and. WR(k) > xr)
            WR(k+1) = WR(k); WI(k+1) = WI(k)
            k = k - 1
        end do
        WR(k+1) = xr; WI(k+1) = xi
    end do
 
    return
end subroutine zrhqr2
 
 
! ------------------------------------------------------------------------------
! zrhqr3: polynomial root finder, degree 3
! Used by: Roots_finder (nres = 0 case)
! ------------------------------------------------------------------------------
subroutine zrhqr3(A, WR, WI)
implicit none
 
    integer, parameter :: N2    = 3
    integer, parameter :: LDA   = 3
    integer, parameter :: LDVL  = 3
    integer, parameter :: LDVR  = 3
    integer, parameter :: LWORK = 15
 
    real*8,  intent(in)  :: A(N2+1)
    real*8,  intent(out) :: WR(N2), WI(N2)
 
    real*8  :: hess(N2,N2), xr, xi
    real*8  :: VL(LDVL,N2), VR(LDVR,N2)
    integer :: WORK(LWORK)
    integer :: info, j, k
 
    WORK(:) = 0
 
    ! Construct companion matrix
    hess(1,:)    = -A(N2:1:-1) / A(N2+1)
    hess(2:N2,:) = 0.0
    do k = 1, N2-1
        hess(k+1,k) = 1.0
    end do
 
    call dgeev('N', 'N', N2, hess, LDA, WR, WI, VL, LDVL, VR, LDVR, &
               WORK, LWORK, INFO)
 
    ! Sort roots by ascending real part (insertion sort)
    do j = 2, N2
        xr = WR(j); xi = WI(j)
        k = j - 1
        do while (k >= 1 .and. WR(k) > xr)
            WR(k+1) = WR(k); WI(k+1) = WI(k)
            k = k - 1
        end do
        WR(k+1) = xr; WI(k+1) = xi
    end do
 
    return
end subroutine zrhqr3
