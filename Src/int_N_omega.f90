! ==============================================================================
! SUBROUTINE: int_N_omega
! ==============================================================================
!
! PURPOSE:
!   Computes the normalisation factor N(omega) for the wave normal angle
!   distribution at a given resonant frequency omega. This factor appears
!   in the denominator of the G2 term in WeightFunction, ensuring the
!   wave normal angle distribution integrates to unity.
!
!   N(omega) is defined as:
!     N(omega) = integral[ sin(theta) * g(theta) * n^2 * c / |dw/dk| * dtheta ]
!   where g(theta) is the two-Gaussian wave normal angle distribution in
!   tan(theta), n^2 is the refractive index squared, and dw/dk is the
!   group velocity along the wave vector. The integral runs from theta_min
!   to theta_max_eff using the trapezoidal rule.
!
!   Dispersion properties (n^2 and dw/dk) are obtained via
!   Get_dispersion_properties from module dispersion, which currently
!   implements the cold plasma (Stix) formulation.
!
! INPUTS:
!   alpha          - Equatorial electron pitch angle [radians]
!   E              - Electron kinetic energy normalised to rest-mass energy
!                    (E = gamma - 1, dimensionless)
!   L              - Magnetic L-shell [Earth radii, dimensionless]
!   lambda         - Magnetic latitude [radians]
!   Theta_function - Array(7) defining the wave normal angle distribution:
!                      (1) theta_m1  : mean angle of 1st Gaussian [degrees]
!                      (2) theta_w1  : width of 1st Gaussian [degrees]
!                      (3) theta_m2  : mean angle of 2nd Gaussian [degrees]
!                      (4) theta_w2  : width of 2nd Gaussian [degrees]
!                      (5) theta_min : lower cutoff angle [degrees]
!                      (6) theta_max : upper cutoff angle [degrees]
!                      (7) theta_q   : relative amplitude of 2nd Gaussian
!   theta_max_eff  - Effective upper cutoff angle, potentially reduced to
!                    0.99 * resonance cone angle [degrees]
!   nres           - Resonance harmonic number [-]
!   omega          - Resonant wave frequency [rad/s]
!   omega_e        - Electron gyrofrequency [rad/s]
!   omega_pe       - Electron plasma frequency [rad/s]
!   omega_i        - Ion gyrofrequency [rad/s]
!   omega_pi       - Ion plasma frequency [rad/s]
!   res_cone       - 'Y' if theta_max_eff is near the resonance cone,
!                    'N' otherwise. Controls integration resolution.
!
! OUTPUTS:
!   N_of_omega     - Normalisation factor N(omega) [rad]
!
! CALLS:
!   Get_dispersion_properties - Computes n^2 and dw/dk (module dispersion)
!
! MODULE DEPENDENCIES:
!   constants  - PI, nint_N, nint_N_res_cone
!   dispersion - Get_dispersion_properties, DispersionProperties
!
! NOTES:
!   - If n^2 < 0 at a given theta point (evanescent wave), that point
!     contributes zero to the integral.
!   - A finer integration grid (nint_N_res_cone points) is used when
!     res_cone = 'Y', to handle the near-singularity of the dispersion
!     relation close to the resonance cone angle.
!   - omega_e, omega_pe, omega_i, omega_pi are now passed explicitly as
!     arguments rather than read from module variables. This replaces the
!     previous implicit dependency on module variables being pre-populated
!     by WeightFunction before this routine is called.
!
! MODIFIED: 2011-10-12 Xudong - revised refractive index and dw/dk calculation
!           2011-10-24 Xudong - switched to lower-sign branch of dispersion
! ==============================================================================

subroutine int_N_omega(alpha, E, L, lambda, Theta_function, theta_max_eff, &
                       nres, omega, k_vec_res, omega_e, omega_pe, omega_i, omega_pi, &
                       res_cone, N_of_omega, whamp_curve)
use constants
use dispersion
implicit none

    real*8, intent(in)  :: alpha, E, L, lambda, Theta_function(7)
    real*8, intent(in)  :: theta_max_eff, omega, k_vec_res
    real*8, intent(in)  :: omega_e, omega_pe, omega_i, omega_pi
    integer, intent(in) :: nres
    character, intent(in) :: res_cone
    real*8, intent(out) :: N_of_omega
    type(WHAMPCurve), intent(in) :: whamp_curve

    ! --- Wave normal angle distribution parameters ---
    real*8 :: theta_min, theta_max
    real*8 :: theta_m1, theta_w1, theta_m2, theta_w2, theta_q
    real*8 :: X_m1, X_w1, X_m2, X_w2

    ! --- Integration ---
    real*8 :: theta, X, step, Ntmp
    integer :: nint_N_tmp, i

    ! --- Dispersion ---
    type(DispersionProperties) :: DP
    real*8 :: k_vec

    real*8 :: R_stix, L_stix, P_stix, S_stix, D_stix
    real*8 :: dPdw_stix, dSdw_stix, dRdw_stix, dLdw_stix

    real*8 :: t_prof_1, t_prof_2

    integer, save :: call_count = 0
    integer, save :: n_debug_writes = 0

    call cpu_time(t_prof_1)
    call_count = call_count + 1

    ! -------------------------------------------------------------------------
    ! Integration resolution: finer near resonance cone
    ! -------------------------------------------------------------------------
    if (res_cone == 'Y') then
        nint_N_tmp = nint_N_res_cone
    else
        nint_N_tmp = nint_N
    end if

    ! -------------------------------------------------------------------------
    ! Unpack wave normal angle distribution parameters
    ! -------------------------------------------------------------------------
    theta_m1  = Theta_function(1)
    theta_w1  = Theta_function(2)
    theta_m2  = Theta_function(3)
    theta_w2  = Theta_function(4)
    theta_min = Theta_function(5)
    theta_max = theta_max_eff
    theta_q   = Theta_function(7)

    X_m1 = tan(theta_m1 / 180.0 * PI)
    X_w1 = tan(theta_w1 / 180.0 * PI)
    X_m2 = tan(theta_m2 / 180.0 * PI)
    X_w2 = tan(theta_w2 / 180.0 * PI)

    step = (theta_max - theta_min) / 180.0 * PI / (real(nint_N_tmp) - 1.0)

    ! Parallel-only (delta-function angle distribution)
    if ((theta_max - theta_min) < 1.0d-10) then
        N_of_omega = 1.0d0
        return
    end if

    ! -------------------------------------------------------------------------
    ! Trapezoidal integration over theta
    ! -------------------------------------------------------------------------
    N_of_omega = 0.0

    if (DispersionType == 'COLD') then
    call cold_plasma_parameters(omega, omega_e, omega_pe, omega_i, omega_pi, &
                                R_stix, L_stix, P_stix, S_stix, D_stix)

    call cold_plasma_stix_derivatives(omega, omega_e, omega_pe, omega_i, omega_pi, &
                                          dPdw_stix, dSdw_stix, dRdw_stix, dLdw_stix)
    end if

    !if (n_debug_writes < 20 .and. E > 1.0d0) then   ! E>1.0 → >511 keV
    !    n_debug_writes = n_debug_writes + 1
    !    open(97, file='resonance_debug_highE.txt', status='unknown', position='append')
    !    write(97,*) 'call=', call_count, &
    !                ' E_keV=', E*511.0d0, &
    !                ' alpha(deg)=', alpha*180.0d0/PI, &
    !                ' omega(Hz)=', omega/(2.0d0*PI), &
    !                ' n_k=', whamp_curve%n_k
    !    close(97)
    !end if

    do i = 1, nint_N_tmp
        theta = theta_min / 180.0 * PI + step * (real(i) - 1.0)
        X     = tan(theta)

        ! HOT case: interpolate n_ref_ind2 and dw_dk from WHAMP curve
        if ((DispersionType == 'HOT' .or. DispersionType == 'EMIC') .and. whamp_curve%n_k > 0) then
            block
                integer :: i_k, i_theta_lo, i_theta_hi, i_s
                real*8  :: w_theta, omega_lo_t, omega_hi_t, frac
                logical :: found, valid_lo, valid_hi
                found = .false.

                ! Find bracketing theta indices and interpolation weight
                i_theta_lo = 1
                i_theta_hi = 1
                w_theta    = 0.0d0
                if (whamp_curve%n_theta > 1) then
                    if (theta <= whamp_curve%theta(1)) then
                        i_theta_lo = 1
                        i_theta_hi = 1
                    else if (theta >= whamp_curve%theta(whamp_curve%n_theta)) then
                        i_theta_lo = whamp_curve%n_theta
                        i_theta_hi = whamp_curve%n_theta
                    else
                        do i_s = 2, whamp_curve%n_theta
                            if (whamp_curve%theta(i_s-1) <= theta .and. &
                                whamp_curve%theta(i_s)   >= theta) then
                                i_theta_lo = i_s - 1
                                i_theta_hi = i_s
                                w_theta = (theta - whamp_curve%theta(i_theta_lo)) / &
                                          (whamp_curve%theta(i_theta_hi) - &
                                           whamp_curve%theta(i_theta_lo))
                                exit
                            end if
                        end do
                    end if
                end if

                do i_k = 2, whamp_curve%n_k
                    valid_lo = whamp_curve%valid(i_theta_lo, i_k) .and. &
                               whamp_curve%valid(i_theta_lo, i_k-1)
                    valid_hi = whamp_curve%valid(i_theta_hi, i_k) .and. &
                               whamp_curve%valid(i_theta_hi, i_k-1)
                    if (.not. valid_lo .or. .not. valid_hi) cycle
                    omega_lo_t = (1.0d0-w_theta)*whamp_curve%omega(i_theta_lo,i_k-1) + &
                                 w_theta*whamp_curve%omega(i_theta_hi,i_k-1)
                    omega_hi_t = (1.0d0-w_theta)*whamp_curve%omega(i_theta_lo,i_k)   + &
                                 w_theta*whamp_curve%omega(i_theta_hi,i_k)
                    if (omega_lo_t <= omega .and. omega_hi_t >= omega) then
                        frac = (omega - omega_lo_t) / (omega_hi_t - omega_lo_t)
                        DP%n_ref_ind2 = (whamp_curve%k(i_k-1)*c/omega_lo_t)**2.0d0 + &
                            frac*((whamp_curve%k(i_k)*c/omega_hi_t)**2.0d0 - &
                                (whamp_curve%k(i_k-1)*c/omega_lo_t)**2.0d0)
                        DP%dw_dk = ((1.0d0-w_theta)*whamp_curve%dw_dk(i_theta_lo,i_k-1) + &
                                    w_theta*whamp_curve%dw_dk(i_theta_hi,i_k-1)) + &
                            frac*( ((1.0d0-w_theta)*whamp_curve%dw_dk(i_theta_lo,i_k) + &
                                    w_theta*whamp_curve%dw_dk(i_theta_hi,i_k)) - &
                                ((1.0d0-w_theta)*whamp_curve%dw_dk(i_theta_lo,i_k-1) + &
                                    w_theta*whamp_curve%dw_dk(i_theta_hi,i_k-1)) )
                        if (k_vec_res < 0.0d0) then
                            DP%dw_dk = -DP%dw_dk
                        end if
                        found = .true.
                        if (abs(DP%dw_dk) < 1.0d-6 * c) then
                            Ntmp = 0.0d0
                            found = .false.
                        end if
                        exit
                    end if
                end do

                ! Tolerance fallback: if omega is just barely outside the
                ! curve's valid range at this theta (a precision edge case
                ! near the whistler cutoff, not a real physics gap), use the
                ! nearest valid endpoint's values rather than returning zero.
                if (.not. found) then
                    block
                        integer :: i_first_valid, i_last_valid, i_s2
                        real*8  :: omega_min_t, omega_max_t, tol
                        logical :: vlo, vhi

                        i_first_valid = 0
                        i_last_valid  = 0
                        do i_s2 = 1, whamp_curve%n_k
                            vlo = whamp_curve%valid(i_theta_lo, i_s2)
                            vhi = whamp_curve%valid(i_theta_hi, i_s2)
                            if (vlo .and. vhi) then
                                if (i_first_valid == 0) i_first_valid = i_s2
                                i_last_valid = i_s2
                            end if
                        end do

                        if (i_first_valid > 0) then
                            omega_min_t = (1.0d0-w_theta)*whamp_curve%omega(i_theta_lo,i_first_valid) + &
                                          w_theta*whamp_curve%omega(i_theta_hi,i_first_valid)
                            omega_max_t = (1.0d0-w_theta)*whamp_curve%omega(i_theta_lo,i_last_valid) + &
                                          w_theta*whamp_curve%omega(i_theta_hi,i_last_valid)
                            tol = 0.001d0

                            if (omega > omega_max_t .and. abs(omega-omega_max_t) < tol*omega_max_t) then
                                DP%n_ref_ind2 = (whamp_curve%k(i_last_valid)*c/omega_max_t)**2.0d0
                                DP%dw_dk = (1.0d0-w_theta)*whamp_curve%dw_dk(i_theta_lo,i_last_valid) + &
                                           w_theta*whamp_curve%dw_dk(i_theta_hi,i_last_valid)
                                if (k_vec_res < 0.0d0) DP%dw_dk = -DP%dw_dk
                                found = .true.
                                if (abs(DP%dw_dk) < 1.0d-6 * c) then
                                    found = .false.
                                end if
                            else if (omega < omega_min_t .and. abs(omega-omega_min_t) < tol*omega_min_t) then
                                DP%n_ref_ind2 = (whamp_curve%k(i_first_valid)*c/omega_min_t)**2.0d0
                                DP%dw_dk = (1.0d0-w_theta)*whamp_curve%dw_dk(i_theta_lo,i_first_valid) + &
                                           w_theta*whamp_curve%dw_dk(i_theta_hi,i_first_valid)
                                if (k_vec_res < 0.0d0) DP%dw_dk = -DP%dw_dk
                                found = .true.
                                if (abs(DP%dw_dk) < 1.0d-6 * c) then
                                    found = .false.
                                end if
                            end if
                        end if
                    end block
                end if

                if (.not. found) then
                    block
                        real*8          :: n2_targeted, dw_targeted
                        complex(kind=8) :: EFL_targeted(3)
                        complex(kind=8) :: BFL_targeted(3)
                        logical         :: valid_targeted
                        call get_whamp_point(whamp_curve%fce_khz, whamp_curve%fpe_khz, &
                                            theta, abs(k_vec_res), omega, &
                                            n2_targeted, dw_targeted, EFL_targeted, &
                                            BFL_targeted, valid_targeted)
                        if (valid_targeted) then
                            DP%n_ref_ind2 = n2_targeted
                            DP%dw_dk      = dw_targeted
                            if (k_vec_res < 0.0d0) DP%dw_dk = -DP%dw_dk
                            found = .true.
                        end if
                    end block
                end if

                if (.not. found) then
                    Ntmp = 0.0d0
                else if (DP%n_ref_ind2 <= 0.0d0 .or. abs(DP%dw_dk) < 1.0d-6 * c) then
                    Ntmp = 0.0d0
                else
                    Ntmp = sin(theta) &
                        * (exp(-((X - X_m1)/X_w1)**2.0) &
                            + theta_q * exp(-((X - X_m2)/X_w2)**2.0)) &
                        * DP%n_ref_ind2 * c / abs(DP%dw_dk)
                end if
            end block

        else
            ! COLD case
            call cold_plasma_refractive_index(omega, theta, &
                                              P_stix, S_stix, R_stix, L_stix, &
                                              DP%n_ref_ind2)
            if (DP%n_ref_ind2 < 0.0) then
                Ntmp = 0.0
            else
                k_vec = sqrt(DP%n_ref_ind2) * omega / c
                block
                    real*8 :: dF1dw, dF2dw, dF3dw, dDdw, dDdk, X2
                    X2 = X*X   ! X = tan(theta), already computed above
                    dDdk = (S_stix*X2 + P_stix) * (c**4.0d0) / (omega**4.0d0) * 4.0d0*(k_vec**3.0d0) &
                           - (R_stix*L_stix*X2 + (2.0d0+X2)*P_stix*S_stix) * (c**2.0d0)/(omega**2.0d0) &
                             * 2.0d0*k_vec
                    dF1dw = ((k_vec*c)**4.0d0) / (omega**5.0d0) * &
                            (omega*X2*dSdw_stix + omega*dPdw_stix - 4.0d0*(S_stix*X2 + P_stix))
                    dF2dw = ((k_vec*c)**2.0d0) / (omega**3.0d0) * &
                            (omega*(X2*(L_stix*dRdw_stix + R_stix*dLdw_stix) + &
                                    (2.0d0+X2)*(P_stix*dSdw_stix + S_stix*dPdw_stix)) - &
                             2.0d0*(R_stix*L_stix*X2 + (2.0d0+X2)*P_stix*S_stix))
                    dF3dw = (1.0d0 + X2) * (R_stix*L_stix*dPdw_stix + P_stix*L_stix*dRdw_stix + &
                                             P_stix*R_stix*dLdw_stix)
                    dDdw = dF1dw - dF2dw + dF3dw
                    DP%dw_dk = -(dDdk / dDdw)
                end block
                Ntmp = sin(theta) &
                       * (exp(-((X - X_m1)/X_w1)**2.0) &
                          + theta_q * exp(-((X - X_m2)/X_w2)**2.0)) &
                       * DP%n_ref_ind2 * c / abs(DP%dw_dk)
            end if
        end if

        ! Trapezoidal weights
        if (i == 1 .or. i == nint_N_tmp) then
            N_of_omega = N_of_omega + Ntmp / 2.0
        else
            N_of_omega = N_of_omega + Ntmp
        end if
    end do

    N_of_omega = N_of_omega * step

    call cpu_time(t_prof_2)
    t_intn_total = t_intn_total + (t_prof_2 - t_prof_1)
    n_intn_calls = n_intn_calls + 1

    return

end subroutine int_N_omega