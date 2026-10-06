! ==============================================================================
! SUBROUTINE: dump_dispersion_curve   (DIAGNOSTIC ONLY)
! ==============================================================================
!
! PURPOSE:
!   Writes the hot (WHAMP) and the corresponding cold dispersion curve for one
!   field line, so the two can be plotted against each other: f(k) at the first
!   wave normal angle of the curve, in units of k/k_h with k_h = omega_pe/c.
!
!   Use this to see where the kinetic branch departs from the cold branch, and
!   whether the hot curve covers the frequency band of interest at all. For the
!   full omega(theta, k) surface rather than one slice, use dump_whamp_surface.
!
! USAGE:
!   Called from sum_n (the call is commented out there); only compiled with
!   -DDIAG. Both files are rewritten on every call, so with several MPI ranks
!   the content is whichever rank wrote last. Use it in a single-rank run.
!
! INPUTS:
!   L_shell - Magnetic L-shell, used for the plasma parameters and the file name
!   curve   - WHAMP dispersion curve for this field line
!
! OUTPUT FILES:
!   dispersion_hot_L<L>.txt   k/k_h, f [Hz]  from the WHAMP curve (valid points)
!   dispersion_cold_L<L>.txt  k/k_h, f [Hz]  from the cold dispersion relation,
!                             swept over 200 logarithmically spaced frequencies
!                             from 1e-3 to 1 times omega_e, at the same theta
!
! NOTES:
!   - Both curves are evaluated at curve%theta(1), i.e. the first wave normal
!     angle of the WHAMP grid, and at the equator (lambda = 0) for the plasma
!     parameters, whereas the WHAMP curve itself was built for the latitude it
!     was computed at. The two therefore only agree when the curve passed in is
!     the equatorial one.
!   - The cold sweep skips frequencies above omega_e*cos(theta), where the
!     whistler branch has its resonance.
!
! ==============================================================================

subroutine dump_dispersion_curve(L_shell, curve)
use constants
use dispersion
implicit none

    real*8,           intent(in) :: L_shell
    type(WHAMPCurve), intent(in) :: curve

    real*8  :: B_local, Ne_local, omega_e, omega_pe, omega_i, omega_pi
    real*8  :: R_stix, L_stix, P_stix, S_stix, D_stix
    real*8  :: kh, omega_sweep, n2_cold, k_cold
    real*8  :: theta_match
    integer :: i_k, n_sweep, unit_hot, unit_cold
    character(len=64) :: fname

    ! Equatorial plasma parameters
    call B_field(L_shell, 0.0d0, B_local)
    call ElectronDensity(L_shell, 0.0d0, Ne_local)
    omega_e  = qe * B_local / me
    omega_pe = sqrt(Ne_local * qe * qe / (C0 * me))
    omega_i  = omega_e * me / mi
    omega_pi = omega_pe * sqrt(me / mi)

    kh = omega_pe / c

    ! Both curves at the first wave normal angle of the WHAMP grid
    theta_match = curve%theta(1)

    ! ---- hot (WHAMP) curve ----
    write(fname, '(A,I0,A)') 'dispersion_hot_L', int(L_shell), '.txt'
    open(newunit=unit_hot, file=trim(fname), status='replace')
    write(unit_hot,*) '# k_over_kh   f_Hz   theta_deg=', theta_match*180.0d0/PI
    do i_k = 1, curve%n_k
        if (curve%valid(1, i_k)) then
            write(unit_hot,*) curve%k(i_k)/kh, curve%omega(1,i_k)/(2.0d0*PI)
        end if
    end do
    close(unit_hot)

    ! ---- cold curve at the same theta ----
    write(fname, '(A,I0,A)') 'dispersion_cold_L', int(L_shell), '.txt'
    open(newunit=unit_cold, file=trim(fname), status='replace')
    write(unit_cold,*) '# k_over_kh   f_Hz   theta_deg=', theta_match*180.0d0/PI
    n_sweep = 200
    do i_k = 1, n_sweep
        omega_sweep = omega_e * 10.0d0**(-3.0d0 + 3.0d0*real(i_k-1,8)/real(n_sweep-1,8))
        if (omega_sweep >= omega_e * cos(theta_match)) cycle    ! whistler resonance
        call cold_plasma_parameters(omega_sweep, omega_e, omega_pe, omega_i, omega_pi, &
                                    R_stix, L_stix, P_stix, S_stix, D_stix)
        call cold_plasma_refractive_index(omega_sweep, theta_match, P_stix, S_stix, &
                                          R_stix, L_stix, n2_cold)
        if (n2_cold > 0.0d0) then
            k_cold = sqrt(n2_cold) * omega_sweep / c
            write(unit_cold,*) k_cold/kh, omega_sweep/(2.0d0*PI)
        end if
    end do
    close(unit_cold)

end subroutine dump_dispersion_curve