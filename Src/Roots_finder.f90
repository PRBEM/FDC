! ==============================================================================
! SUBROUTINE: Roots_finder
! ==============================================================================
!
! PURPOSE:
!   Finds the resonant wave frequencies (roots of the combined resonance
!   condition and cold plasma dispersion relation) for a given electron
!   state (alpha, E), field line location (L, lambda), resonance harmonic
!   (nres), and wave normal angle (theta).
!
!   The resonance condition is:
!     omega - k_par * v_par - n * Omega_e / gamma = 0
!   Combined with the cold plasma dispersion relation, this yields a
!   polynomial in omega whose real roots within the wave frequency band
!   [omega_lc, omega_uc] are the resonant frequencies.
!
!   For nres /= 0: a 10th-degree polynomial in omega is solved (general case).
!   For nres  = 0: symmetry reduces the problem to a cubic in omega^2,
!                  solved as a 3rd-degree polynomial.
!
!   Polynomial roots are found using LAPACK via zrhqr2 (nres /= 0) or
!   zrhqr3 (nres = 0).
!
! INPUTS:
!   alpha   - Equatorial electron pitch angle [radians]
!   E       - Electron kinetic energy normalised to rest-mass energy
!             (E = gamma - 1, dimensionless)
!   L       - Magnetic L-shell [Earth radii, dimensionless]
!   lambda  - Magnetic latitude [radians]
!   nres    - Resonance harmonic number n (integer, e.g. -1, 0, 1, ...)
!   theta   - Wave normal angle [radians]
!   omega_e    - Electron gyrofrequency [rad/s]
!   omega_pe   - Electron plasma frequency [rad/s]
!   omega_i    - Ion gyrofrequency [rad/s]
!   omega_pi   - Ion plasma frequency [rad/s]
!   omega_lc   - Lower frequency cutoff [rad/s]
!   omega_uc   - Upper frequency cutoff [rad/s]
!
! OUTPUTS:
!   roots(4) - Up to 4 resonant frequencies [rad/s]. Unused slots are set
!              to 0.0. Valid roots are real, lie within [omega_lc, omega_uc],
!              and satisfy the group velocity condition |b2| > 1.
!
! MODULE DEPENDENCIES:
!   constants  - sqe, and indirectly via called routines
!
! CALLS:
!   Local_alpha  - Converts equatorial pitch angle to local pitch angle
!   zrhqr2       - LAPACK-based polynomial root finder (degree 10, nres /= 0)
!   zrhqr3       - LAPACK-based polynomial root finder (degree 3,  nres = 0)
!
! NOTES:
!   - a5 is hardcoded to 0.0. The original commented-out definition was
!     -(omega_pi^2 * sqe*omega_e + omega_pe^2 * omega_i). Whether the zero
!     value is a valid simplification or a bug should be verified.
!   - For nres = 0, only even powers of omega appear, so the polynomial is
!     solved in omega^2 and roots are recovered by taking the square root.
!     Coefficients are taken from alternating entries of the general coeff()
!     array (indices 5, 7, 9, 11).
!   - Root validity checks:
!       nres /= 0: root must be real (rti = 0), lie in [omega_lc, omega_uc],
!                  and satisfy |(omega - b1)*b2/omega| > 1
!       nres  = 0: root^2 must be real, lie in [omega_lc^2, omega_uc^2],
!                  and |b2| > 1
!
! ==============================================================================

subroutine Roots_finder(alpha, E, L, lambda, nres, theta, &
                        omega_e, omega_pe, omega_i, omega_pi, &
                        omega_lc, omega_uc, &
                        roots, whamp_curve)
use constants
use dispersion
implicit none

    real*8, intent(in)  :: alpha, E, L, lambda, theta
    integer, intent(in) :: nres
    real*8, intent(in)  :: omega_e, omega_pe, omega_i, omega_pi
    real*8, intent(in)  :: omega_lc, omega_uc
    real*8, intent(out) :: roots(4)
    type(WHAMPCurve), intent(in) :: whamp_curve

    ! Cold plasma variables
    real*8 :: coeff(11)
    real*8 :: rtr(10), rti(10)

    real*8 :: coeff1(4)
    real*8 :: rtr1(3), rti1(3)

    ! Plasma frequency combination coefficients
    real*8 :: a1, a2, a3, a4, a5, a6

    ! Intermediate polynomial coefficients (local to this routine)
    real*8 :: a7, a8, a9, a10, a11, a12, a13, a14, a15, a16

    ! Resonance condition parameters
    real*8 :: b1, b2

    ! Particle state
    real*8 :: gamma1, betta, alpha_ne

    ! Wave normal angle
    real*8 :: X

    integer :: i, NN
    real*8  :: Temp2

    ! HOT plasma variables
    real*8  :: p, vpar
    real*8  :: f_prev_pos, f_curr_pos, f_prev_neg, f_curr_neg
    real*8  :: k_prev, k_curr, omega_prev, omega_curr
    real*8  :: k_res, omega_res, frac
    integer :: i_k, i_theta_lo, i_theta_hi
    real*8  :: w_theta
    logical :: root_found

    real*8 :: t_prof_1, t_prof_2

    real*8 :: f_at_zero_pos, f_at_zero_neg
    real*8 :: omega_first, k_first
    real*8 :: f_first_pos, f_first_neg

    call cpu_time(t_prof_1)

    betta  = sqrt(E*(E + 2.0)) / (E + 1.0)
    gamma1 = E + 1.0

    call Local_alpha(alpha, lambda, alpha_ne)

    X = tan(theta)

    ! -------------------------------------------------------------------------
    ! HOT PLASMA CASE: find resonant k from WHAMP dispersion curve
    ! -------------------------------------------------------------------------
    if (DispersionType == 'HOT' .and. whamp_curve%n_k > 0) then
        p    = me * c * sqrt((E + 1.0d0)**2.0d0 - 1.0d0)
        vpar = p / gamma1 / me * cos(alpha_ne)
 
        roots(1) = 0.0d0
        roots(2) = 0.0d0
        roots(3) = 0.0d0
        roots(4) = 0.0d0
        root_found = .false.
 
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
                do i_k = 2, whamp_curve%n_theta
                    if (whamp_curve%theta(i_k-1) <= theta .and. &
                        whamp_curve%theta(i_k)   >= theta) then
                        i_theta_lo = i_k - 1
                        i_theta_hi = i_k
                        w_theta = (theta - whamp_curve%theta(i_theta_lo)) / &
                                  (whamp_curve%theta(i_theta_hi) - &
                                   whamp_curve%theta(i_theta_lo))
                        exit
                    end if
                end do
            end if
        end if
 
        ! Scan k points for sign change in resonance condition
        ! f(k) = omega(k) - k*cos(theta)*vpar - nres*sqe*omega_e/gamma1
        ! Scan both +k and -k since the resonant wavevector sign depends on
        ! the harmonic number and propagation geometry.
        k_prev      = 0.0d0
        omega_prev  = 0.0d0
        f_prev_pos  = 0.0d0
        f_prev_neg  = 0.0d0

        ! DEBUG - verify whether the resonance crossing happens below k_min
        !if (E > 1.9d0 .and. E < 2.0d0 .and. abs(alpha_ne*180.0d0/PI - 10.0d0) < 1.0d0 &
        !    .and. nres == -1) then
        !    open(89, file='rootscan_debug2.txt', status='unknown', position='append')
        !    write(89,*) 'theta=', theta, ' vpar=', vpar, &
        !                ' f_neg_at_kmin=', ( ((1.0d0-w_theta)*whamp_curve%omega(i_theta_lo,1) + &
        !                                    w_theta*whamp_curve%omega(i_theta_hi,1)) &
        !                                    + whamp_curve%k(1)*cos(theta)*vpar &
        !                                    - real(nres)*real(sqe)*omega_e/gamma1 ), &
        !                ' f_neg_at_kmax=', ( ((1.0d0-w_theta)*whamp_curve%omega(i_theta_lo,whamp_curve%n_k) + &
        !                                    w_theta*whamp_curve%omega(i_theta_hi,whamp_curve%n_k)) &
        !                                    + whamp_curve%k(whamp_curve%n_k)*cos(theta)*vpar &
        !                                    - real(nres)*real(sqe)*omega_e/gamma1 )
        !    close(89)
        !end if

        if (whamp_curve%valid(i_theta_lo, 1) .and. whamp_curve%valid(i_theta_hi, 1)) then
            f_at_zero_pos = 0.0d0 - real(nres)*real(sqe)*omega_e/gamma1
            f_at_zero_neg = 0.0d0 - real(nres)*real(sqe)*omega_e/gamma1  ! same at k=0 since k*vpar term vanishes

            omega_first = (1.0d0-w_theta)*whamp_curve%omega(i_theta_lo,1) + w_theta*whamp_curve%omega(i_theta_hi,1)
            k_first     = whamp_curve%k(1)

            f_first_pos = omega_first - k_first*cos(theta)*vpar - real(nres)*real(sqe)*omega_e/gamma1
            f_first_neg = omega_first + k_first*cos(theta)*vpar - real(nres)*real(sqe)*omega_e/gamma1

            if (f_at_zero_pos * f_first_pos < 0.0d0) then
                frac = f_at_zero_pos / (f_at_zero_pos - f_first_pos)
                omega_res = 0.0d0 + frac * (omega_first - 0.0d0)
                if (omega_res > omega_lc .and. omega_res < omega_uc) then
                    roots(1)   = omega_res
                    root_found = .true.
                end if
            else if (.not. root_found .and. f_at_zero_neg * f_first_neg < 0.0d0) then
                frac = f_at_zero_neg / (f_at_zero_neg - f_first_neg)
                omega_res = 0.0d0 + frac * (omega_first - 0.0d0)
                if (omega_res > omega_lc .and. omega_res < omega_uc) then
                    roots(1)   = omega_res
                    root_found = .true.
                end if
            end if
        end if

        do i_k = 1, whamp_curve%n_k
            if (.not. whamp_curve%valid(i_theta_lo, i_k)) cycle
            if (.not. whamp_curve%valid(i_theta_hi, i_k)) cycle
            k_curr     = whamp_curve%k(i_k)
            omega_curr = (1.0d0 - w_theta) * whamp_curve%omega(i_theta_lo, i_k) + &
                         w_theta * whamp_curve%omega(i_theta_hi, i_k)
            f_curr_pos = omega_curr - k_curr*cos(theta)*vpar &
                         - real(nres)*real(sqe)*omega_e/gamma1
            f_curr_neg = omega_curr + k_curr*cos(theta)*vpar &
                         - real(nres)*real(sqe)*omega_e/gamma1

            if (i_k > 1 .and. k_prev > 0.0d0) then
                if (f_prev_pos * f_curr_pos < 0.0d0) then
                    frac      = f_prev_pos / (f_prev_pos - f_curr_pos)
                    omega_res = omega_prev + frac * (omega_curr - omega_prev)
                    if (omega_res > omega_lc .and. omega_res < omega_uc) then
                        roots(1)   = omega_res
                        root_found = .true.
                        exit
                    end if
                end if
                if (f_prev_neg * f_curr_neg < 0.0d0) then
                    frac      = f_prev_neg / (f_prev_neg - f_curr_neg)
                    omega_res = omega_prev + frac * (omega_curr - omega_prev)
                    if (omega_res > omega_lc .and. omega_res < omega_uc) then
                        roots(1)   = omega_res
                        root_found = .true.
                        exit
                    end if
                end if
            end if

            k_prev     = k_curr
            omega_prev = omega_curr
            f_prev_pos = f_curr_pos
            f_prev_neg = f_curr_neg
        end do

        if (E > 0.18d0 .and. E < 0.2d0 .and. abs(alpha_ne*180.0d0/PI - 56.67d0) < 1.0d0 &
            .and. (nres == 1 .or. nres == -1)) then
            open(88, file='root_a0005_debug.txt', status='unknown', position='append')
            write(88,*) 'nres=', nres, ' root_found=', root_found, ' omega=', roots(1)
            close(88)
        end if
        return
    end if

    ! -------------------------------------------------------------------------
    ! COLD PLASMA CASE: polynomial root finding
    ! -------------------------------------------------------------------------
    ! Plasma frequency combination coefficients
    ! Written to module variables (a1-a6) for subsequent use in WeightFunction.
    ! Also used locally below as inputs to the polynomial coefficients a7-a16.
    ! -------------------------------------------------------------------------
    a1 = omega_i**2.0 + omega_e**2.0 + omega_pe**2.0 + omega_pi**2.0
    a2 = (omega_i*omega_e)**2.0 + (omega_pe*omega_i)**2.0 + (omega_pi*omega_e)**2.0
    a3 = omega_i + real(sqe)*omega_e
    a4 = omega_i*real(sqe)*omega_e - omega_pe**2.0 - omega_pi**2.0
    a5 = 0.0   ! analytically zero; see note in header
    a6 = omega_pi**2.0 + omega_pe**2.0

    ! Intermediate polynomial coefficients (local only)
    a7  = (X**2.0)*a1 + omega_i**2.0 + omega_e**2.0 + a6
    a8  = (X**2.0)*a2 + (omega_i*omega_e)**2.0 + a6*(omega_i**2.0 + omega_e**2.0)
    a9  = a6*((omega_i*omega_e)**2.0)
    a10 = 2.0*a4 - a3**2.0 - a6
    a11 = a4**2.0 - 2.0*a3*a5 - 2.0*a4*a6 + (a3**2.0)*a6
    a12 = a5**2.0 + (a4**2.0)*a6 - 2.0*a3*a5*a6
    a13 = (a5**2.0)*a6
    a14 = (X**2.0)*(2.0*a4 - a3**2.0) - (2.0 + X**2.0)*(a1 + a6)
    a15 = (X**2.0)*(a4**2.0 - 2.0*a3*a5) + (2.0 + X**2.0)*(a2 + a1*a6)
    a16 = (X**2.0)*(a5**2.0) + (2.0 + X**2.0)*a2*a6

    ! Resonance condition parameters
    b1 = real(nres)*real(sqe)*omega_e / gamma1
    b2 = 1.0 / (betta*cos(alpha_ne)*cos(theta))

    ! -------------------------------------------------------------------------
    ! Polynomial coefficients
    ! -------------------------------------------------------------------------
    coeff(11) = (1 + X**2.0)*((1 - b2**2.0)**2.0)
    coeff(10) = 4.0*(1 + X**2.0)*b1*(b2**2.0)*(1 - b2**2.0)
    coeff(9)  = (6.0*(1 + X**2.0)*(b1**2.0) - a7)*(b2**4.0) &
                - (a14 + 2.0*(1 + X**2.0)*(b1**2.0))*(b2**2.0) &
                + (1 + X**2.0)*a10
    coeff(8)  = (4.0*a7*b1 - 4.0*(1 + X**2.0)*(b1**3.0))*(b2**4.0) &
                + 2.0*a14*b1*(b2**2.0)
    coeff(7)  = (a8 - 6.0*a7*(b1**2.0) + (1 + X**2.0)*(b1**4.0))*(b2**4.0) &
                - (a15 + a14*(b1**2.0))*(b2**2.0) &
                + (1 + X**2.0)*a11
    coeff(6)  = (4.0*a7*(b1**3.0) - 4.0*a8*b1)*(b2**4.0) &
                + 2.0*a15*b1*(b2**2.0)
    coeff(5)  = (6.0*a8*(b1**2.0) - a7*(b1**4.0) - a9)*(b2**4.0) &
                - (a15*(b1**2.0) - a16)*(b2**2.0) &
                - (1 + X**2.0)*a12
    coeff(4)  = (4.0*a9*b1 - 4.0*a8*(b1**3.0))*(b2**4.0) &
                - 2.0*a16*b1*(b2**2.0)
    coeff(3)  = (a8*(b1**4.0) - 6.0*a9*(b1**2.0))*(b2**4.0) &
                + a16*(b1**2.0)*(b2**2.0) &
                + (1 + X**2.0)*a13
    coeff(2)  = 4.0*a9*(b1**3.0)*(b2**4.0)
    coeff(1)  = -a9*(b1**4.0)*(b2**4.0)

    ! -------------------------------------------------------------------------
    ! Root finding
    ! -------------------------------------------------------------------------
    roots(1) = 0.0
    roots(2) = 0.0
    roots(3) = 0.0
    roots(4) = 0.0
    NN = 0

    if (nres /= 0) then

        call zrhqr2(coeff, rtr, rti)

        do i = 1, 10
            if (rti(i) == 0) then
                if ((rtr(i) > omega_lc) .and. (rtr(i) < omega_uc)) then
                    Temp2 = abs((rtr(i) - b1)*b2 / rtr(i))
                    if (Temp2 > 1.0) then
                        NN = NN + 1
                        roots(NN) = rtr(i)
                    end if
                end if
            end if
        end do

    else

        ! nres = 0: polynomial is even in omega; solve as cubic in omega^2
        coeff1(4) = coeff(11)
        coeff1(3) = coeff(9)
        coeff1(2) = coeff(7)
        coeff1(1) = coeff(5)

        call zrhqr3(coeff1, rtr1, rti1)

        do i = 1, 3
            if (rti1(i) == 0) then
                if ((rtr1(i) > (omega_lc**2.0)) .and. &
                    (rtr1(i) < (omega_uc**2.0)) .and. &
                    (abs(b2) > 1.0)) then
                    NN = NN + 1
                    roots(NN) = sqrt(rtr1(i))
                end if
            end if
        end do

    end if

    if (E > 1.9d0 .and. E < 2.0d0 .and. abs(alpha_ne*180.0d0/PI - 10.0d0) < 1.0d0 &
        .and. nres == -1 .and. trim(DispersionType) == 'COLD') then
        open(88, file='cold_omega_nresmin1.txt', status='unknown', position='append')
        write(88,*) 'theta=', theta*180.0d0/PI, ' NN=', NN, &
                    ' root1(Hz)=', roots(1)/(2.0d0*PI), &
                    ' root2(Hz)=', roots(2)/(2.0d0*PI), &
                    ' root3(Hz)=', roots(3)/(2.0d0*PI), &
                    ' root4(Hz)=', roots(4)/(2.0d0*PI)
        close(88)
    end if

    n_roots_calls = n_roots_calls + 1

    call cpu_time(t_prof_2)
    t_roots_total = t_roots_total + (t_prof_2 - t_prof_1)
    n_roots_calls = n_roots_calls + 1

    return

end subroutine Roots_finder