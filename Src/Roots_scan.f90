! ==============================================================================
! FILE: Roots_scan.f90
! ==============================================================================
!
! Cold-plasma resonant frequency finder that replaces the degree-10
! polynomial + companion-matrix eigensolve (zrhqr2/zrhqr3) in Roots_finder.
!
! CONTENTS
!   module rootcheck_stats  - counters for the -DROOTCHECK comparison mode
!   roots_scan_cold         - the new root finder
!   rootcheck_compare       - compares polynomial and scan roots (ROOTCHECK)
!   rootcheck_summary       - per-rank summary at the end of a run (ROOTCHECK)
!
! METHOD
!   Roots_finder's polynomial encodes, for both branches of the cold
!   dispersion relation, the condition
!
!       n^2_branch(omega, theta) = n_res(omega)^2,
!       n_res(omega) = (omega - b1) * b2 / omega,
!
!   where n_res is the refractive index demanded by the resonance condition
!   (the same quantity as Temp2 in Roots_finder). Instead of expanding this
!   into a polynomial, roots_scan_cold looks for the zeros of
!
!       H(omega) = n^2_branch(omega) - n_res(omega)^2
!
!   directly, and only inside [omega_lc, omega_uc]:
!     1. evaluate H on NSCAN log-spaced frequencies across the band,
!     2. refine every sign change with Brent's method,
!     3. split cells where A = S*tan^2(theta) + P changes sign (resonance
!        cone: one branch has a pole there, which would otherwise hide a
!        root next to it or be mistaken for one),
!     4. at local extrema of H that approach zero without crossing it,
!        search for a hidden pair of close roots (golden section),
!     5. apply the same acceptance tests as Roots_finder: omega strictly
!        inside the band and |n_res| > 1; poles are rejected by requiring
!        |H| to be small at the converged point.
!
!   n^2 uses the same R, L, P, S as cold_plasma_parameters and the same two
!   branches as the polynomial. The branches are evaluated in a
!   cancellation-free form that is algebraically identical to
!   cold_plasma_refractive_index:
!     minus branch = (B - sqrt(disc)) / (2A)   (the one WeightFunction uses)
!     plus  branch = (B + sqrt(disc)) / (2A)
!   with disc written as (RL - PS)^2 X^4 + 4 P^2 D^2 (1 + X^2), which equals
!   B^2 - 4AC exactly and is never negative.
!
! VALIDITY GUARD
!   If the band contains the electron or ion gyrofrequency (poles of R, L,
!   S), or the inputs are not finite, roots_scan_cold returns ok = .false.
!   and Roots_finder falls back to the polynomial method.
!
! ==============================================================================

module rootcheck_stats
    implicit none
    integer(8), save :: rc_calls    = 0   ! calls compared
    integer(8), save :: rc_differ   = 0   ! calls where the root sets differ
    integer(8), save :: rc_missing  = 0   ! polynomial roots not found by scan
    integer(8), save :: rc_extra    = 0   ! scan roots not found by polynomial
    integer(8), save :: rc_matched  = 0   ! roots found by both
    integer(8), save :: rc_fallback = 0   ! calls where scan declined (ok=.false.)
    integer(8), save :: rc_overflow = 0   ! calls with more than 4 roots
    integer(8), save :: rc_bins(4)  = 0   ! matched rel. diff > 1e-12, 1e-10, 1e-8, 1e-6
    real*8,     save :: rc_maxrel   = 0.0d0
    integer,    save :: rc_unit     = -1
    logical,    save :: rc_open     = .false.
    integer,    save :: rc_lines    = 0
    integer, parameter :: RC_MAX_LINES = 5000
    real*8,  parameter :: RC_MATCH_TOL = 1.0d-6   ! roots closer than this are "the same"
end module rootcheck_stats


! ==============================================================================
! SUBROUTINE: roots_scan_cold
! ==============================================================================
! INPUTS:
!   b1, b2     - resonance parameters exactly as computed in Roots_finder:
!                  b1 = nres * sqe * omega_e / gamma
!                  b2 = 1 / (beta * cos(alpha_local) * cos(theta))
!   theta      - wave normal angle [rad]
!   omega_e, omega_pe, omega_i, omega_pi - as in Roots_finder [rad/s]
!   omega_lc, omega_uc - frequency band [rad/s]
! OUTPUTS:
!   roots(4)   - resonant frequencies, ascending, unused slots = 0
!   nroots     - number of roots found (may exceed 4; only 4 are returned)
!   ok         - .false. if the method declined; caller must fall back
! ==============================================================================
subroutine roots_scan_cold(b1, b2, theta, omega_e, omega_pe, omega_i, omega_pi, &
                           omega_lc, omega_uc, roots, nroots, ok)
implicit none
    real*8,  intent(in)  :: b1, b2, theta
    real*8,  intent(in)  :: omega_e, omega_pe, omega_i, omega_pi
    real*8,  intent(in)  :: omega_lc, omega_uc
    real*8,  intent(out) :: roots(4)
    integer, intent(out) :: nroots
    logical, intent(out) :: ok

    integer, parameter :: NSCAN    = 32      ! scan points across the band
    integer, parameter :: NBRANCH  = 2       ! 1 = minus branch, 2 = plus branch
    integer, parameter :: MAXFOUND = 16
    real*8,  parameter :: SQE      = -1.0d0  ! sign of electron charge (as in constants)
    real*8,  parameter :: POLE_EPS = 1.0d-12 ! relative gap left around a pole
    real*8,  parameter :: ACCEPT   = 1.0d-6  ! |H| <= ACCEPT*(|n^2| + n_res^2) at a root

    real*8  :: wpe2, wpi2, X2, X4, onepX2, twopX2
    real*8  :: w(NSCAN), Hs(NSCAN, NBRANCH), As(NSCAN)
    real*8  :: n2m, n2p, Aval, ratio, wp
    real*8  :: found(MAXFOUND)
    integer :: nfound, k, br, i, j
    real*8  :: tmp, sg

    ok     = .true.
    nroots = 0
    nfound = 0
    roots  = 0.0d0

    ! ---- validity guard: the scan assumes R, L, S are smooth in the band ----
    if (.not. (omega_lc > 0.0d0 .and. omega_uc > omega_lc)) then
        ok = .false.; return
    end if
    if ((omega_lc <= omega_e .and. omega_e <= omega_uc) .or. &
        (omega_lc <= omega_i .and. omega_i <= omega_uc)) then
        ok = .false.; return
    end if
    if (b2 /= b2 .or. abs(b2) > huge(1.0d0) .or. b1 /= b1) then
        ok = .false.; return
    end if

    ! ---- per-call constants ----
    wpe2   = omega_pe * omega_pe
    wpi2   = omega_pi * omega_pi
    X2     = tan(theta)**2
    X4     = X2 * X2
    onepX2 = 1.0d0 + X2
    twopX2 = 2.0d0 + X2

    ! ---- scan: log-spaced samples, endpoints exactly at the band edges ----
    ratio = exp(log(omega_uc / omega_lc) / real(NSCAN - 1, 8))
    w(1)  = omega_lc
    do k = 2, NSCAN - 1
        w(k) = w(k-1) * ratio
    end do
    w(NSCAN) = omega_uc

    do k = 1, NSCAN
        call eval_n2(w(k), n2m, n2p, Aval)
        tmp = nres(w(k))**2
        Hs(k, 1) = n2m - tmp
        Hs(k, 2) = n2p - tmp
        As(k)    = Aval
        if (Hs(k,1) /= Hs(k,1) .or. As(k) /= As(k)) then   ! NaN guard
            ok = .false.; return
        end if
    end do

    do br = 1, NBRANCH

        ! ---- sign changes, with cells containing a pole split in two ----
        do k = 1, NSCAN - 1
            if (As(k) * As(k+1) < 0.0d0) then
                wp = brent_zero(w(k), w(k+1), As(k), As(k+1), 0)
                call try_bracket(w(k), wp * (1.0d0 - POLE_EPS), br)
                call try_bracket(wp * (1.0d0 + POLE_EPS), w(k+1), br)
            else
                call try_bracket_known(w(k), w(k+1), Hs(k,br), Hs(k+1,br), br)
            end if
        end do

        ! ---- hidden pairs: H approaches zero between samples ----
        ! Checked on every triple of neighbouring samples with no sign change
        ! and no pole. A dip is suspected if the middle sample is a local
        ! extremum, or if the parabola through the three samples has its
        ! extremum inside the triple and approaches zero. The parabola test
        ! also covers the first and last cell of the band, where a dip
        ! next to the band edge leaves no interior extremum among the samples.
        do k = 2, NSCAN - 1
            if (As(k-1)*As(k) > 0.0d0 .and. As(k)*As(k+1) > 0.0d0 .and. &
                Hs(k-1,br)*Hs(k,br) > 0.0d0 .and. Hs(k,br)*Hs(k+1,br) > 0.0d0) then
                sg = sign(1.0d0, Hs(k,br))
                if (dip_possible(w(k-1), w(k), w(k+1), &
                                 sg*Hs(k-1,br), sg*Hs(k,br), sg*Hs(k+1,br))) then
                    call search_hidden_pair(w(k-1), w(k+1), sg, br)
                end if
            end if
        end do

    end do

    ! ---- sort ascending, return up to 4 ----
    do i = 2, nfound
        tmp = found(i)
        j = i - 1
        do while (j >= 1)
            if (found(j) <= tmp) exit
            found(j+1) = found(j)
            j = j - 1
        end do
        found(j+1) = tmp
    end do
    nroots = nfound
    do i = 1, min(nfound, 4)
        roots(i) = found(i)
    end do

contains

    ! Refractive index squared on both branches, and A (zero at the cone)
    subroutine eval_n2(om, n2minus, n2plus, Aout)
        real*8, intent(in)  :: om
        real*8, intent(out) :: n2minus, n2plus, Aout
        real*8 :: R, L, P, S, D, B, C, disc, q, r1, r2

        ! identical to cold_plasma_parameters (sqe = -1)
        R = 1.0d0 - wpe2 / (om * (om + SQE * omega_e)) - wpi2 / (om * (om + omega_i))
        L = 1.0d0 - wpe2 / (om * (om - SQE * omega_e)) - wpi2 / (om * (om - omega_i))
        P = 1.0d0 - wpe2 / (om * om) - wpi2 / (om * om)
        S = 0.5d0 * (R + L)
        D = 0.5d0 * (R - L)

        Aout = S * X2 + P
        B    = R * L * X2 + twopX2 * P * S
        C    = P * R * L * onepX2
        disc = (R*L - P*S)**2 * X4 + 4.0d0 * P*P * D*D * onepX2

        ! cancellation-free roots of A n^4 - B n^2 + C = 0
        if (B >= 0.0d0) then
            q = B + sqrt(disc)
        else
            q = B - sqrt(disc)
        end if
        if (Aout /= 0.0d0) then
            r1 = q / (2.0d0 * Aout)
        else
            r1 = huge(1.0d0)
        end if
        if (q /= 0.0d0) then
            r2 = 2.0d0 * C / q
        else
            r2 = 0.0d0
        end if
        if (B >= 0.0d0) then
            n2minus = r2; n2plus = r1
        else
            n2minus = r1; n2plus = r2
        end if
    end subroutine eval_n2

    real*8 function nres(om)
        real*8, intent(in) :: om
        nres = b2 * (om - b1) / om
    end function nres

    ! br = 0: A(omega); br = 1, 2: H on the minus / plus branch
    real*8 function hfun(om, br)
        real*8,  intent(in) :: om
        integer, intent(in) :: br
        real*8 :: a, b, c
        call eval_n2(om, a, b, c)
        select case (br)
            case (0);    hfun = c
            case (1);    hfun = a - nres(om)**2
            case default; hfun = b - nres(om)**2
        end select
    end function hfun

    subroutine try_bracket(a, b, br)
        real*8,  intent(in) :: a, b
        integer, intent(in) :: br
        if (b <= a) return
        call try_bracket_known(a, b, hfun(a, br), hfun(b, br), br)
    end subroutine try_bracket

    subroutine try_bracket_known(a, b, fa, fb, br)
        real*8,  intent(in) :: a, b, fa, fb
        integer, intent(in) :: br
        if (fa == 0.0d0) then
            call accept_root(a, br)
        else if (fa * fb < 0.0d0) then
            call accept_root(brent_zero(a, b, fa, fb, br), br)
        end if
    end subroutine try_bracket_known

    ! Same acceptance as Roots_finder, plus rejection of poles
    subroutine accept_root(x, br)
        real*8,  intent(in) :: x
        integer, intent(in) :: br
        real*8  :: a, b, c, n2, nr2
        integer :: i
        if (.not. (x > omega_lc .and. x < omega_uc)) return
        nr2 = nres(x)**2
        if (nr2 <= 1.0d0) return                          ! |n_res| > 1
        call eval_n2(x, a, b, c)
        if (br == 1) then
            n2 = a
        else
            n2 = b
        end if
        if (abs(n2 - nr2) > ACCEPT * (abs(n2) + nr2)) return   ! pole, not a root
        do i = 1, nfound                                  ! duplicate guard
            if (abs(found(i) - x) <= 1.0d-12 * x) return
        end do
        if (nfound < MAXFOUND) then
            nfound = nfound + 1
            found(nfound) = x
        end if
    end subroutine accept_root

    ! f1, f2, f3 > 0 are sgn*H at x1 < x2 < x3. True if H may approach
    ! zero between x1 and x3 without a sign change at the samples.
    logical function dip_possible(x1, x2, x3, f1, f2, f3)
        real*8, intent(in) :: x1, x2, x3, f1, f2, f3
        real*8 :: d1, d2, cc, xv, fv
        dip_possible = .false.
        if (f2 < f1 .and. f2 < f3) then          ! sampled local minimum
            dip_possible = .true.; return
        end if
        d1 = (f2 - f1) / (x2 - x1)
        d2 = (f3 - f2) / (x3 - x2)
        cc = (d2 - d1) / (x3 - x1)               ! half the second derivative
        if (cc <= 0.0d0) return                  ! no minimum of the parabola
        xv = 0.5d0 * (x1 + x2) - d1 / (2.0d0 * cc)
        if (xv <= x1 .or. xv >= x3) return
        fv = f1 + d1 * (xv - x1) + cc * (xv - x1) * (xv - x2)
        dip_possible = (fv < 0.5d0 * min(f1, f2, f3))
    end function dip_possible

    ! Golden-section search on [lo, hi] for a point where H changes sign
    ! relative to sgn; if found, the two sub-intervals are refined.
    subroutine search_hidden_pair(lo_in, hi_in, sgn, br)
        real*8,  intent(in) :: lo_in, hi_in, sgn
        integer, intent(in) :: br
        real*8, parameter :: GR = 0.6180339887498949d0
        real*8  :: lo, hi, x1, x2, f1, f2, xs
        integer :: it
        lo = lo_in; hi = hi_in
        x1 = hi - GR * (hi - lo); f1 = sgn * hfun(x1, br)
        x2 = lo + GR * (hi - lo); f2 = sgn * hfun(x2, br)
        xs = -1.0d0
        do it = 1, 40
            if (f1 < 0.0d0) then; xs = x1; exit; end if
            if (f2 < 0.0d0) then; xs = x2; exit; end if
            if (f1 < f2) then
                hi = x2; x2 = x1; f2 = f1
                x1 = hi - GR * (hi - lo); f1 = sgn * hfun(x1, br)
            else
                lo = x1; x1 = x2; f1 = f2
                x2 = lo + GR * (hi - lo); f2 = sgn * hfun(x2, br)
            end if
        end do
        if (xs > 0.0d0) then
            call try_bracket(lo_in, xs, br)
            call try_bracket(xs, hi_in, br)
        end if
    end subroutine search_hidden_pair

    ! Brent's method (Brent 1973) on [a0, b0] with fa0*fb0 < 0
    real*8 function brent_zero(a0, b0, fa0, fb0, br)
        real*8,  intent(in) :: a0, b0, fa0, fb0
        integer, intent(in) :: br
        real*8  :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm
        integer :: iter
        a = a0; b = b0; fa = fa0; fb = fb0
        c = b;  fc = fb
        d = b - a; e = d
        do iter = 1, 200
            if (fb * fc > 0.0d0) then
                c = a; fc = fa; d = b - a; e = d
            end if
            if (abs(fc) < abs(fb)) then
                a = b; b = c; c = a
                fa = fb; fb = fc; fc = fa
            end if
            tol1 = 2.0d0 * epsilon(1.0d0) * abs(b)
            xm   = 0.5d0 * (c - b)
            if (abs(xm) <= tol1 .or. fb == 0.0d0) exit
            if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
                s = fb / fa
                if (a == c) then
                    p = 2.0d0 * xm * s
                    q = 1.0d0 - s
                else
                    q = fa / fc
                    r = fb / fc
                    p = s * (2.0d0 * xm * q * (q - r) - (b - a) * (r - 1.0d0))
                    q = (q - 1.0d0) * (r - 1.0d0) * (s - 1.0d0)
                end if
                if (p > 0.0d0) q = -q
                p = abs(p)
                if (2.0d0 * p < min(3.0d0 * xm * q - abs(tol1 * q), abs(e * q))) then
                    e = d; d = p / q
                else
                    d = xm; e = d
                end if
            else
                d = xm; e = d
            end if
            a = b; fa = fb
            if (abs(d) > tol1) then
                b = b + d
            else
                b = b + sign(tol1, xm)
            end if
            fb = hfun(b, br)
        end do
        brent_zero = b
    end function brent_zero

end subroutine roots_scan_cold


! ==============================================================================
! SUBROUTINE: rootcheck_compare   (used only with -DROOTCHECK)
! ==============================================================================
! Compares the polynomial roots (roots_poly, npoly <= 4, ascending) with the
! scan roots and records differences. Discrepancies are written to a
! rank-specific file rootcheck_rank_NNN.log (at most RC_MAX_LINES lines).
! ==============================================================================
subroutine rootcheck_compare(E, alpha, lambda, nres, theta, &
                             roots_poly, npoly, roots_scan, nscan, scan_ok)
use rootcheck_stats
implicit none
    include "mpif.h"
    real*8,  intent(in) :: E, alpha, lambda, theta
    integer, intent(in) :: nres, npoly, nscan
    real*8,  intent(in) :: roots_poly(4), roots_scan(4)
    logical, intent(in) :: scan_ok

    integer :: i, j, ns, nmiss, nextra, rank, ierr
    real*8  :: best, rel
    logical :: used(4)
    character(len=64) :: fname

    rc_calls = rc_calls + 1
    if (.not. scan_ok) then
        rc_fallback = rc_fallback + 1
        return
    end if
    if (nscan > 4 .or. npoly > 4) rc_overflow = rc_overflow + 1

    used = .false.
    nmiss = 0
    do i = 1, npoly
        best = huge(1.0d0); j = 0
        do ns = 1, min(nscan, 4)
            if (used(ns)) cycle
            rel = abs(roots_scan(ns) - roots_poly(i)) / abs(roots_poly(i))
            if (rel < best) then
                best = rel; j = ns
            end if
        end do
        if (j > 0 .and. best <= RC_MATCH_TOL) then
            used(j) = .true.
            rc_matched = rc_matched + 1
            rc_maxrel  = max(rc_maxrel, best)
            if (best > 1.0d-12) rc_bins(1) = rc_bins(1) + 1
            if (best > 1.0d-10) rc_bins(2) = rc_bins(2) + 1
            if (best > 1.0d-8)  rc_bins(3) = rc_bins(3) + 1
            if (best > 1.0d-6)  rc_bins(4) = rc_bins(4) + 1
        else
            nmiss = nmiss + 1
        end if
    end do
    nextra = count(.not. used(1:min(nscan, 4)))

    if (nmiss > 0 .or. nextra > 0) then
        rc_differ  = rc_differ + 1
        rc_missing = rc_missing + nmiss
        rc_extra   = rc_extra + nextra
        if (rc_lines < RC_MAX_LINES) then
            if (.not. rc_open) then
                call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
                write(fname, '(a,i3.3,a)') 'rootcheck_rank_', rank, '.log'
                open(newunit=rc_unit, file=trim(fname), status='replace', action='write')
                rc_open = .true.
                write(rc_unit, '(a)') '# E[mc^2] alpha_eq[deg] lambda[deg] nres theta[deg] ' // &
                                      '| poly roots [rad/s] | scan roots [rad/s]'
            end if
            write(rc_unit, '(es12.5,2f9.4,i5,f10.5,a,4es22.14,a,4es22.14)') &
                E, alpha*57.29577951308232d0, lambda*57.29577951308232d0, nres, &
                theta*57.29577951308232d0, ' |', roots_poly, ' |', roots_scan
            rc_lines = rc_lines + 1
        end if
    end if
end subroutine rootcheck_compare


! ==============================================================================
! SUBROUTINE: rootcheck_summary   (used only with -DROOTCHECK)
! ==============================================================================
subroutine rootcheck_summary()
use rootcheck_stats
implicit none
    include "mpif.h"
    integer :: rank, ierr
    if (rc_calls == 0) return
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    write(*,'(a,i4,a)') ' === ROOTCHECK SUMMARY (rank', rank, ') ==='
    write(*,'(a,i12)')  '   calls compared:           ', rc_calls
    write(*,'(a,i12)')  '   calls with differences:   ', rc_differ
    write(*,'(a,i12)')  '   poly roots missed by scan:', rc_missing
    write(*,'(a,i12)')  '   extra roots from scan:    ', rc_extra
    write(*,'(a,i12)')  '   roots matched:            ', rc_matched
    write(*,'(a,4i10)') '   matched rel.diff > 1e-12, 1e-10, 1e-8, 1e-6:', rc_bins
    write(*,'(a,es10.2)') '   max rel.diff (matched):   ', rc_maxrel
    write(*,'(a,i12)')  '   scan declined (fallback): ', rc_fallback
    write(*,'(a,i12)')  '   calls with > 4 roots:     ', rc_overflow
    if (rc_open) close(rc_unit)
end subroutine rootcheck_summary