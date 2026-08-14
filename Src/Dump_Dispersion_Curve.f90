subroutine dump_dispersion_curve(L_shell, curve)
use constants
use dispersion
implicit none
    real*8, intent(in) :: L_shell
    type(WHAMPCurve), intent(in) :: curve

    real*8 :: B_local, Ne_local, omega_e, omega_pe, omega_i, omega_pi
    real*8 :: R_stix, L_stix, P_stix, S_stix, D_stix
    real*8 :: kh, omega_sweep, n2_cold, k_cold
    real*8 :: theta_match
    integer :: i_k, n_sweep
    character(len=64) :: fname

    integer :: i_k_dbg

    !open(76, file='dump_debug_check.txt', status='replace')
    !write(76,*) 'n_k=', curve%n_k, ' n_theta=', curve%n_theta
    !do i_k_dbg = 1, 20
    !    write(76,*) 'i_k=', i_k_dbg, ' valid=', curve%valid(1,i_k_dbg), &
    !                ' k=', curve%k(i_k_dbg), ' omega=', curve%omega(1,i_k_dbg)
    !end do
    !close(76)

    ! Equatorial plasma parameters
    call B_field(L_shell, 0.0d0, B_local)
    call ElectronDensity(L_shell, 0.0d0, Ne_local)
    omega_e  = qe * B_local / me
    omega_pe = sqrt(Ne_local * qe * qe / (C0 * me))
    omega_i  = omega_e * me / mi
    omega_pi = omega_pe * sqrt(me / mi)

    kh = omega_pe / c

    ! Use the same wave normal angle as the hot curve first theta point
    theta_match = curve%theta(1)

    ! Hot curve 
    write(fname, '(A,I0,A)') 'dispersion_hot_L', int(L_shell), '.txt'
    open(77, file=trim(fname), status='replace')
    write(77,*) '# k_over_kh   f_Hz   theta_deg=', theta_match*180.0d0/PI
    do i_k = 1, curve%n_k
        if (curve%valid(1, i_k)) then
            write(77,*) curve%k(i_k)/kh, curve%omega(1,i_k)/(2.0d0*PI)
        end if
    end do
    close(77)

    ! Cold curve
    write(fname, '(A,I0,A)') 'dispersion_cold_L', int(L_shell), '.txt'
    open(78, file=trim(fname), status='replace')
    write(78,*) '# k_over_kh   f_Hz   theta_deg=', theta_match*180.0d0/PI
    n_sweep = 200
    do i_k = 1, n_sweep
        omega_sweep = omega_e * 10.0d0**(-3.0d0 + 3.0d0*real(i_k-1,8)/real(n_sweep-1,8))
        if (omega_sweep >= omega_e * cos(theta_match)) cycle
        call cold_plasma_parameters(omega_sweep, omega_e, omega_pe, omega_i, omega_pi, &
                                    R_stix, L_stix, P_stix, S_stix, D_stix)
        call cold_plasma_refractive_index(omega_sweep, theta_match, P_stix, S_stix, &
                                          R_stix, L_stix, n2_cold)
        if (n2_cold > 0.0d0) then
            k_cold = sqrt(n2_cold) * omega_sweep / c
            write(78,*) k_cold/kh, omega_sweep/(2.0d0*PI)
        end if
    end do
    close(78)

end subroutine dump_dispersion_curve