! ==============================================================================
! bench_roots: standalone benchmark for the degree-10 polynomial root finder
! used in Roots_finder (zrhqr2). Answers two questions without a profiler:
!
!   1. How much does one zrhqr2 call cost on its own? Comparing this with the
!      ~19 us measured per Roots_finder call tells how much of Roots_finder is
!      the LAPACK eigensolve and how much is the coefficient construction.
!   2. What do the cheaper alternatives cost, and do they give the same roots?
!        A: zrhqr2 as in FDC (dgeev: balancing + Hessenberg reduction + QR)
!        B: dhseqr directly on the companion matrix (no balancing)
!        C: dgebal (scaling only) + dhseqr (keeps balancing, skips dgehrd)
!
! The companion matrix built by zrhqr2 (full first row + unit subdiagonal) is
! already upper Hessenberg, which is why B and C are possible at all.
!
! Test polynomials: roots with magnitudes 1e3-1e6 (rad/s scale, as in FDC),
! mixed signs, 6 real roots + 2 complex pairs, random per set.
!
! Build (same libraries as FDC):
!   gfortran -O2 bench_roots.f90 zrhqr.f90 -o bench_roots \
!            -L/usr/lib/ -l:liblapack.so.3 -l:libblas.so.3
! Run:
!   ./bench_roots
! ==============================================================================
program bench_roots
    implicit none
    integer, parameter :: N = 10, NSETS = 500, NREP = 400
    real*8     :: coeff(N+1, NSETS)
    real*8     :: wr(N), wi(N), wrA(N, NSETS), wiA(N, NSETS)
    real*8     :: tsec(3), maxrel(3), rel
    integer    :: nreal_mismatch(3), s, r, i, v
    integer(8) :: t0, t1, rate

    call make_polynomials(coeff)

    ! Reference roots from the current code (variant A)
    do s = 1, NSETS
        call zrhqr2(coeff(:, s), wrA(:, s), wiA(:, s))
    end do

    do v = 1, 3
        ! --- timing ---
        call system_clock(t0, rate)
        do r = 1, NREP
            do s = 1, NSETS
                call solve(v, coeff(:, s), wr, wi)
            end do
        end do
        call system_clock(t1)
        tsec(v) = real(t1 - t0, 8) / real(rate, 8)

        ! --- agreement with variant A ---
        maxrel(v) = 0.0d0
        nreal_mismatch(v) = 0
        do s = 1, NSETS
            call solve(v, coeff(:, s), wr, wi)
            ! FDC accepts a root only if its imaginary part is exactly zero,
            ! so the number of exactly-real roots must match
            if (count(wi == 0.0d0) /= count(wiA(:, s) == 0.0d0)) then
                nreal_mismatch(v) = nreal_mismatch(v) + 1
                cycle
            end if
            do i = 1, N
                if (wi(i) == 0.0d0 .and. wiA(i, s) == 0.0d0) then
                    rel = abs(wr(i) - wrA(i, s)) / abs(wrA(i, s))
                    maxrel(v) = max(maxrel(v), rel)
                end if
            end do
        end do
    end do

    write(*,'(a)') 'Variant                         us/call   speed-up   ' // &
                   'max rel diff   real-count mismatches'
    write(*,'(a,f9.2,f10.2,es15.2,i10)') 'A  dgeev (current zrhqr2)    ', &
        1.0d6*tsec(1)/(NREP*NSETS), 1.0d0, maxrel(1), nreal_mismatch(1)
    write(*,'(a,f9.2,f10.2,es15.2,i10)') 'B  dhseqr, no balancing      ', &
        1.0d6*tsec(2)/(NREP*NSETS), tsec(1)/tsec(2), maxrel(2), nreal_mismatch(2)
    write(*,'(a,f9.2,f10.2,es15.2,i10)') 'C  dgebal(S) + dhseqr        ', &
        1.0d6*tsec(3)/(NREP*NSETS), tsec(1)/tsec(3), maxrel(3), nreal_mismatch(3)
    write(*,'(a,i0,a)') '(', NSETS, ' polynomials; mismatch = number of polynomials ' // &
        'where the count of exactly-real roots differs from A)'

contains

    subroutine solve(v, a, wr, wi)
        integer, intent(in)  :: v
        real*8,  intent(in)  :: a(N+1)
        real*8,  intent(out) :: wr(N), wi(N)
        select case (v)
            case (1); call zrhqr2(a, wr, wi)
            case (2); call roots_hseqr(a, wr, wi, .false.)
            case (3); call roots_hseqr(a, wr, wi, .true.)
        end select
    end subroutine solve

    subroutine roots_hseqr(a, wr, wi, balance)
        real*8,  intent(in)  :: a(N+1)
        real*8,  intent(out) :: wr(N), wi(N)
        logical, intent(in)  :: balance
        real*8  :: h(N, N), scal(N), z(1, 1), work(N)
        integer :: k, ilo, ihi, info

        ! Same companion matrix as zrhqr2
        h(1, :)   = -a(N:1:-1) / a(N+1)
        h(2:N, :) = 0.0d0
        do k = 1, N-1
            h(k+1, k) = 1.0d0
        end do

        ilo = 1; ihi = N
        ! 'S' = scaling only: D^-1 H D stays upper Hessenberg
        if (balance) call dgebal('S', N, h, N, ilo, ihi, scal, info)

        call dhseqr('E', 'N', N, ilo, ihi, h, N, wr, wi, z, 1, work, N, info)
        call sort_roots(wr, wi)
    end subroutine roots_hseqr

    subroutine sort_roots(wr, wi)
        ! Same insertion sort as zrhqr2 (ascending real part)
        real*8, intent(inout) :: wr(N), wi(N)
        real*8  :: xr, xi
        integer :: j, k
        do j = 2, N
            xr = wr(j); xi = wi(j)
            k = j - 1
            do while (k >= 1)
                if (wr(k) <= xr) exit
                wr(k+1) = wr(k); wi(k+1) = wi(k)
                k = k - 1
            end do
            wr(k+1) = xr; wi(k+1) = xi
        end do
    end subroutine sort_roots

    subroutine make_polynomials(c)
        ! Coefficients c(1..N+1) in ascending powers, from random roots
        real*8, intent(out) :: c(N+1, NSETS)
        complex*16 :: rt(N), p(0:N)
        real*8     :: u(4)
        integer    :: s, k, j
        call random_seed()
        do s = 1, NSETS
            do k = 1, 6
                call random_number(u)
                rt(k) = cmplx(sign(10.0d0**(3.0d0 + 3.0d0*u(1)), u(2) - 0.5d0), &
                              0.0d0, kind=8)
            end do
            do k = 7, 9, 2
                call random_number(u)
                rt(k)   = cmplx(10.0d0**(3.0d0 + 3.0d0*u(1)) * sign(1.0d0, u(2) - 0.5d0), &
                                10.0d0**(2.0d0 + 3.0d0*u(3)), kind=8)
                rt(k+1) = conjg(rt(k))
            end do
            p = (0.0d0, 0.0d0); p(0) = (1.0d0, 0.0d0)
            do k = 1, N
                do j = k, 1, -1
                    p(j) = p(j-1) - rt(k) * p(j)
                end do
                p(0) = -rt(k) * p(0)
            end do
            c(:, s) = real(p(0:N), 8)
        end do
    end subroutine make_polynomials

end program bench_roots
