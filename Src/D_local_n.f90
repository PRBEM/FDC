! ==============================================================================
! SUBROUTINE: D_local_n
! ==============================================================================
!
! PURPOSE:
!   Computes the local quasi-linear diffusion coefficients (D_aa, D_ap, D_pp)
!   for a single resonance harmonic at a fixed magnetic latitude (lambda), by
!   integrating the diffusion coefficient integrand over wave normal angle theta.
!
!   The integration uses the trapezoidal rule over a uniform grid in theta,
!   from theta_min to theta_max (taken from Theta_function). The integrand at
!   each theta point is provided by WeightFunction.
!
!   The result is the local (non-bounce-averaged) diffusion rate. Bounce-
!   averaging over latitude and summation over resonance harmonics are
!   performed by the calling routine.
!
!
! INPUTS:
!   alpha          - Equatorial electron pitch angle [radians]
!   E              - Electron kinetic energy normalised to rest-mass energy
!                    (E = gamma - 1, dimensionless)
!   L              - Magnetic L-shell [Earth radii, dimensionless]
!   lambda         - Magnetic latitude [radians]
!   Bw             - RMS wave magnetic field amplitude [T]
!   Theta_function - Array(7) defining the wave normal angle distribution:
!                      (1) theta_m1  : mean angle of 1st Gaussian [degrees]
!                      (2) theta_w1  : width of 1st Gaussian [degrees]
!                      (3) theta_m2  : mean angle of 2nd Gaussian [degrees]
!                      (4) theta_w2  : width of 2nd Gaussian [degrees]
!                      (5) theta_min : lower cutoff of integration range [degrees]
!                      (6) theta_max : upper cutoff of integration range [degrees]
!                      (7) theta_q   : relative amplitude of 2nd Gaussian
!   nres           - Resonance harmonic number n (integer, e.g. -1, 0, 1, ...)
!   whamp_curve    - (optional) Hot plasma dispersion curve from WHAMP.
!                    Passed through to WeightFunction. If absent, cold plasma
!                    dispersion is used.
!
! OUTPUTS:
!   DDN(1)         - Local D_aa at this latitude [s^-1] (pitch-angle diffusion)
!   DDN(2)         - Local D_ap at this latitude [s^-1] (mixed diffusion)
!   DDN(3)         - Local D_pp at this latitude [s^-1] (momentum diffusion)
!
! CALLS:
!   WeightFunction - Diffusion coefficient integrand at a single theta value
!
! MODULE DEPENDENCIES:
!   constants      - Physical constants and unit conversion factors (RAD, PI, ...)
!   (nint_X is read from a module: number of theta integration points)
!   dispersion     - WHAMPCurve type
!
! INTEGRATION METHOD:
!   Trapezoidal rule on a uniform grid of nint_X points in theta [radians].
!   The integrand is weighted by sin(theta) before integration, consistent
!   with integration over wave normal angle in spherical geometry.
!
! NOTES:
!   - theta_min and theta_max are read in degrees from Theta_function(5:6)
!     and converted to radians internally.
!   - nint_X is a global parameter from a module rather than an argument;
!     this is a known issue to be addressed in future refactoring.

subroutine D_local_n(alpha, E, L, nres, lambda, Bw, Theta_function, DDN, whamp_curve)
use constants
use dispersion
implicit none

    real*8,intent(in) :: alpha, E, L, lambda, Bw, Theta_function(7)
    integer,intent(in) :: nres
    real*8,intent(out) :: DDN(3)
    type(WHAMPCurve), intent(in) :: whamp_curve

    real*8 :: WeightValue(3)
    real*8 :: theta_min, theta_max, theta
    real*8 :: BB1, BB2, BB3
    integer :: i

    ! --- CUNNINGHAM path: pre-loop setup variables ---
    real*8  :: Ne_local, B_local, B_eq, lambda_eq
    real*8  :: omega_e, omega_e_eq, omega_i, omega_pe, omega_pi
    real*8  :: omega_m, d_omega, omega_lc, omega_uc
    real*8  :: betta, gamma1, p, vpar, p_perp, alpha_ne
    real*8  :: b1, integral_gX
    real*8  :: theta_m1, theta_w1, theta_m2, theta_w2, theta_q
    real*8  :: X_m1, X_w1, X_m2, X_w2
 
    ! --- CUNNINGHAM path: per-theta / per-root variables ---
    integer :: j
    real*8  :: X, gX_raw, roots(4), k_sign(4), root_omega
    real*8  :: stix_R, stix_L, stix_P, stix_S, stix_D
    real*8  :: mu2, k_vec, k_par, k_perp
    real*8  :: dD_dk, dD_dw, sterm, phi_sq, norm_factor, PSD, norm_intensity
    real*8  :: DnXaa, DnXap, DnXpp
    real*8  :: term1, term2, term3
    real*8  :: n_ref_ind2, dw_dk
    real*8  :: argJ, BesselF1, BesselF2, BesselF3
    real*8  :: hot_jacobian
    type(DispersionProperties) :: DP


    theta_min = Theta_function(5)
    theta_max = Theta_function(6)

    BB1 = (theta_max  -theta_min) * RAD
    BB2 = 1.0 / (real(nint_X) - 1.0)
    BB3 = BB1 * BB2

    DDN(1) = 0.0
    DDN(2) = 0.0
    DDN(3) = 0.0

    ! Parallel-only (delta-function angle distribution)
    if ((theta_max - theta_min) < 1.0d-10) then
        block
            real*8 :: WeightValue_parallel(3)
            call WeightFunction(alpha, E, L, lambda, Bw, Theta_function, nres, &
                                theta_min/180.0d0*PI, WeightValue_parallel, whamp_curve)
            DDN(1) = WeightValue_parallel(1)
            DDN(2) = WeightValue_parallel(2)
            DDN(3) = WeightValue_parallel(3)
        end block
        return
    end if

    if (trim(DiffusionFormulation) == 'CUNNINGHAM') then

      theta_m1 = Theta_function(1)
      theta_w1 = Theta_function(2)
      theta_m2 = Theta_function(3)
      theta_w2 = Theta_function(4)
      theta_q  = Theta_function(7)

      X_m1 = tan(theta_m1/180.0d0*PI)
      X_w1 = tan(theta_w1/180.0d0*PI)
      X_m2 = tan(theta_m2/180.0d0*PI)
      X_w2 = tan(theta_w2/180.0d0*PI)


      ! Local plasma state (same calls as WeightFunction)
      call Local_alpha(alpha, lambda, alpha_ne)
      call ElectronDensity(L, lambda, Ne_local)

      lambda_eq = 0.0d0
      call B_field(L, lambda,    B_local)
      call B_field(L, lambda_eq, B_eq)

      omega_e    = qe * B_local / me
      omega_e_eq = qe * B_eq    / me
      omega_i    = omega_e * me / mi
      omega_pe   = sqrt(Ne_local * qe * qe / (C0 * me))
      omega_pi   = omega_pe * sqrt(me / mi)

      omega_m  = 0.0d0
      d_omega  = 0.0d0
      omega_lc = 0.0d0
      omega_uc = 0.0d0
      call Freq_Spectral(L, lambda, omega_e_eq, omega_m, d_omega, omega_lc, omega_uc)

      ! Particle state
      betta  = sqrt(E*(E + 2.0d0)) / (E + 1.0d0)
      gamma1 = E + 1.0d0
      p      = me * c * sqrt((E + 1.0d0)**2.0d0 - 1.0d0)
      vpar   = p / gamma1 / me * cos(alpha_ne)
      p_perp = p * sin(alpha_ne)

      ! Resonance condition term b1 = nres*sqe*omega_e/gamma1
      b1 = real(nres, 8) * real(sqe, 8) * omega_e / gamma1

      ! g(X) normalization integral (once for this lambda, not per theta/root)
      call gX_normalization_integral(X_m1, X_w1, X_m2, X_w2, theta_q, &
                                      tan(theta_min/180.0d0*PI), &
                                      tan(theta_max/180.0d0*PI), &
                                      nint_X, integral_gX)

    end if  

    do i = 1, nint_X
      theta = theta_min / 180.0 * PI + (real(i) - 1.0) * BB3
      
      WeightValue(1) = 0.0
      WeightValue(2) = 0.0
      WeightValue(3) = 0.0

      if (trim(DiffusionFormulation) /= 'CUNNINGHAM') then
        call WeightFunction(alpha, E, L, lambda, Bw, Theta_function, nres, theta, &
                            WeightValue, whamp_curve)

      else
        X = tan(theta)
        ! Unnormalised g(X) at this theta (same definition as in int_N_omega's Ntmp)
        gX_raw = exp(-((X - X_m1)/X_w1)**2.0d0) &
                  + theta_q * exp(-((X - X_m2)/X_w2)**2.0d0)


        ! Resonant frequencies at this theta
        roots = 0.0d0
        k_sign = 0.0d0
        call Roots_finder(alpha, E, L, lambda, nres, theta, &
                          omega_e, omega_pe, omega_i, omega_pi, &
                          omega_lc, omega_uc, &
                          roots, whamp_curve, k_sign)

        DnXaa = 0.0d0
        DnXap = 0.0d0
        DnXpp = 0.0d0

        do j = 1, 4

            root_omega = roots(j)

            if (root_omega == 0.0d0) cycle
            ! Skip roots below lower hybrid (same check as WeightFunction)
            if (root_omega <= sqrt((omega_pe/omega_e)**2.0d0 + 1.0d0)) cycle

            if (trim(DispersionType) == 'HOT' .or. trim(DispersionType) == 'EMIC') then

                ! Initialize
                n_ref_ind2 = 0.0d0
                dw_dk      = 0.0d0
                k_vec = 0.0d0
                DP%EFL(1)  = (0.0d0, 0.0d0)
                DP%EFL(2)  = (0.0d0, 0.0d0)
                DP%EFL(3)  = (0.0d0, 0.0d0)
                DP%BFL(1) = (0.0d0, 0.0d0)
                DP%BFL(2) = (0.0d0, 0.0d0)
                DP%BFL(3) = (0.0d0, 0.0d0)


              if (whamp_curve%n_k > 0) then
                block
                    integer :: i_k, i_theta_lo, i_theta_hi
                    real*8  :: w_theta, omega_lo_t, omega_hi_t
                    real*8  :: frac, n2_lo, n2_hi, dw_lo, dw_hi
                    logical :: valid_lo, valid_hi, found

                    ! Find bracketing theta indices
                    i_theta_lo = 1; i_theta_hi = 1; w_theta = 0.0d0
                    if (whamp_curve%n_theta > 1) then
                        if (theta <= whamp_curve%theta(1)) then
                            i_theta_lo = 1; i_theta_hi = 1
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

                    ! Find bracketing k by interpolated omega
                    found = .false.
                    do i_k = 2, whamp_curve%n_k
                        valid_lo = whamp_curve%valid(i_theta_lo, i_k) .and. &
                                  whamp_curve%valid(i_theta_lo, i_k-1)
                        valid_hi = whamp_curve%valid(i_theta_hi, i_k) .and. &
                                  whamp_curve%valid(i_theta_hi, i_k-1)
                        if (.not. valid_lo .or. .not. valid_hi) cycle

                        omega_lo_t = (1.0d0-w_theta)*whamp_curve%omega(i_theta_lo,i_k-1) + &
                                      w_theta*whamp_curve%omega(i_theta_hi,i_k-1)
                        omega_hi_t = (1.0d0-w_theta)*whamp_curve%omega(i_theta_lo,i_k) + &
                                      w_theta*whamp_curve%omega(i_theta_hi,i_k)

                        if (omega_lo_t <= root_omega .and. omega_hi_t >= root_omega) then
                            frac  = (root_omega - omega_lo_t) / (omega_hi_t - omega_lo_t)
                            n2_lo = (whamp_curve%k(i_k-1) * c / omega_lo_t)**2.0d0
                            n2_hi = (whamp_curve%k(i_k)   * c / omega_hi_t)**2.0d0
                            dw_lo = (1.0d0-w_theta)*whamp_curve%dw_dk(i_theta_lo,i_k-1) + &
                                    w_theta*whamp_curve%dw_dk(i_theta_hi,i_k-1)
                            dw_hi = (1.0d0-w_theta)*whamp_curve%dw_dk(i_theta_lo,i_k) + &
                                    w_theta*whamp_curve%dw_dk(i_theta_hi,i_k)
                            n_ref_ind2 = n2_lo + frac * (n2_hi - n2_lo)
                            dw_dk      = dw_lo + frac * (dw_hi - dw_lo)
                            if (abs(dw_dk) < 1.0d-6 * c) n_ref_ind2 = 0.0d0
                            DP%EFL(1) = (1.0d0-w_theta)*whamp_curve%EFL(1,i_theta_lo,i_k-1) + &
                                        w_theta*whamp_curve%EFL(1,i_theta_hi,i_k-1) + &
                                        frac*( ((1.0d0-w_theta)*whamp_curve%EFL(1,i_theta_lo,i_k) + &
                                                w_theta*whamp_curve%EFL(1,i_theta_hi,i_k)) - &
                                                ((1.0d0-w_theta)*whamp_curve%EFL(1,i_theta_lo,i_k-1) + &
                                                w_theta*whamp_curve%EFL(1,i_theta_hi,i_k-1)) )
                            DP%EFL(2) = (1.0d0-w_theta)*whamp_curve%EFL(2,i_theta_lo,i_k-1) + &
                                        w_theta*whamp_curve%EFL(2,i_theta_hi,i_k-1) + &
                                        frac*( ((1.0d0-w_theta)*whamp_curve%EFL(2,i_theta_lo,i_k) + &
                                                w_theta*whamp_curve%EFL(2,i_theta_hi,i_k)) - &
                                                ((1.0d0-w_theta)*whamp_curve%EFL(2,i_theta_lo,i_k-1) + &
                                                w_theta*whamp_curve%EFL(2,i_theta_hi,i_k-1)) )
                            DP%EFL(3) = (1.0d0-w_theta)*whamp_curve%EFL(3,i_theta_lo,i_k-1) + &
                                        w_theta*whamp_curve%EFL(3,i_theta_hi,i_k-1) + &
                                        frac*( ((1.0d0-w_theta)*whamp_curve%EFL(3,i_theta_lo,i_k) + &
                                                w_theta*whamp_curve%EFL(3,i_theta_hi,i_k)) - &
                                                ((1.0d0-w_theta)*whamp_curve%EFL(3,i_theta_lo,i_k-1) + &
                                                w_theta*whamp_curve%EFL(3,i_theta_hi,i_k-1)) )
                            DP%BFL(1) = (1.0d0-w_theta)*whamp_curve%BFL(1,i_theta_lo,i_k-1) + &
                                        w_theta*whamp_curve%BFL(1,i_theta_hi,i_k-1) + &
                                        frac*( ((1.0d0-w_theta)*whamp_curve%BFL(1,i_theta_lo,i_k) + &
                                                w_theta*whamp_curve%BFL(1,i_theta_hi,i_k)) - &
                                                ((1.0d0-w_theta)*whamp_curve%BFL(1,i_theta_lo,i_k-1) + &
                                                w_theta*whamp_curve%BFL(1,i_theta_hi,i_k-1)) )
                            DP%BFL(2) = (1.0d0-w_theta)*whamp_curve%BFL(2,i_theta_lo,i_k-1) + &
                                        w_theta*whamp_curve%BFL(2,i_theta_hi,i_k-1) + &
                                        frac*( ((1.0d0-w_theta)*whamp_curve%BFL(2,i_theta_lo,i_k) + &
                                                w_theta*whamp_curve%BFL(2,i_theta_hi,i_k)) - &
                                                ((1.0d0-w_theta)*whamp_curve%BFL(2,i_theta_lo,i_k-1) + &
                                                w_theta*whamp_curve%BFL(2,i_theta_hi,i_k-1)) )
                            DP%BFL(3) = (1.0d0-w_theta)*whamp_curve%BFL(3,i_theta_lo,i_k-1) + &
                                        w_theta*whamp_curve%BFL(3,i_theta_hi,i_k-1) + &
                                        frac*( ((1.0d0-w_theta)*whamp_curve%BFL(3,i_theta_lo,i_k) + &
                                                w_theta*whamp_curve%BFL(3,i_theta_hi,i_k)) - &
                                                ((1.0d0-w_theta)*whamp_curve%BFL(3,i_theta_lo,i_k-1) + &
                                                w_theta*whamp_curve%BFL(3,i_theta_hi,i_k-1)) )
                            found = .true.
                            exit
                        end if
                    end do

                    ! Tolerance fallback: nearest valid endpoint
                    if (.not. found) then
                        block
                            integer :: i_first_valid, i_last_valid, i_s2
                            real*8  :: omega_min_t, omega_max_t, tol
                            logical :: vlo2, vhi2
                            i_first_valid = 0; i_last_valid = 0
                            do i_s2 = 1, whamp_curve%n_k
                                vlo2 = whamp_curve%valid(i_theta_lo, i_s2)
                                vhi2 = whamp_curve%valid(i_theta_hi, i_s2)
                                if (vlo2 .and. vhi2) then
                                    if (i_first_valid == 0) i_first_valid = i_s2
                                    i_last_valid = i_s2
                                end if
                            end do
                            if (i_first_valid > 0) then
                                omega_min_t = (1.0d0-w_theta)*whamp_curve%omega(i_theta_lo,i_first_valid) + &
                                              w_theta*whamp_curve%omega(i_theta_hi,i_first_valid)
                                omega_max_t = (1.0d0-w_theta)*whamp_curve%omega(i_theta_lo,i_last_valid) + &
                                              w_theta*whamp_curve%omega(i_theta_hi,i_last_valid)
                                tol = 0.01d0
                                if (root_omega > omega_max_t .and. &
                                    abs(root_omega-omega_max_t) < tol*omega_max_t) then
                                    n_ref_ind2 = (whamp_curve%k(i_last_valid)*c/omega_max_t)**2.0d0
                                    dw_dk = (1.0d0-w_theta)*whamp_curve%dw_dk(i_theta_lo,i_last_valid) + &
                                            w_theta*whamp_curve%dw_dk(i_theta_hi,i_last_valid)
                                    if (abs(dw_dk) < 1.0d-6*c) n_ref_ind2 = 0.0d0
                                    DP%EFL(1) = (1.0d0-w_theta)*whamp_curve%EFL(1,i_theta_lo,i_last_valid) + &
                                                w_theta*whamp_curve%EFL(1,i_theta_hi,i_last_valid)
                                    DP%EFL(2) = (1.0d0-w_theta)*whamp_curve%EFL(2,i_theta_lo,i_last_valid) + &
                                                w_theta*whamp_curve%EFL(2,i_theta_hi,i_last_valid)
                                    DP%EFL(3) = (1.0d0-w_theta)*whamp_curve%EFL(3,i_theta_lo,i_last_valid) + &
                                                w_theta*whamp_curve%EFL(3,i_theta_hi,i_last_valid)
                                    DP%BFL(1) = (1.0d0-w_theta)*whamp_curve%BFL(1,i_theta_lo,i_last_valid) + &
                                                w_theta*whamp_curve%BFL(1,i_theta_hi,i_last_valid)
                                    DP%BFL(2) = (1.0d0-w_theta)*whamp_curve%BFL(2,i_theta_lo,i_last_valid) + &
                                                w_theta*whamp_curve%BFL(2,i_theta_hi,i_last_valid)
                                    DP%BFL(3) = (1.0d0-w_theta)*whamp_curve%BFL(3,i_theta_lo,i_last_valid) + &
                                                w_theta*whamp_curve%BFL(3,i_theta_hi,i_last_valid)
                                    found = .true.
                                else if (root_omega < omega_min_t .and. &
                                        abs(root_omega-omega_min_t) < tol*omega_min_t) then
                                    n_ref_ind2 = (whamp_curve%k(i_first_valid)*c/omega_min_t)**2.0d0
                                    dw_dk = (1.0d0-w_theta)*whamp_curve%dw_dk(i_theta_lo,i_first_valid) + &
                                            w_theta*whamp_curve%dw_dk(i_theta_hi,i_first_valid)
                                    if (abs(dw_dk) < 1.0d-6*c) n_ref_ind2 = 0.0d0
                                    DP%EFL(1) = (1.0d0-w_theta)*whamp_curve%EFL(1,i_theta_lo,i_first_valid) + &
                                                w_theta*whamp_curve%EFL(1,i_theta_hi,i_first_valid)
                                    DP%EFL(2) = (1.0d0-w_theta)*whamp_curve%EFL(2,i_theta_lo,i_first_valid) + &
                                                w_theta*whamp_curve%EFL(2,i_theta_hi,i_first_valid)
                                    DP%EFL(3) = (1.0d0-w_theta)*whamp_curve%EFL(3,i_theta_lo,i_first_valid) + &
                                                w_theta*whamp_curve%EFL(3,i_theta_hi,i_first_valid)
                                    DP%BFL(1) = (1.0d0-w_theta)*whamp_curve%BFL(1,i_theta_lo,i_first_valid) + &
                                                w_theta*whamp_curve%BFL(1,i_theta_hi,i_first_valid)
                                    DP%BFL(2) = (1.0d0-w_theta)*whamp_curve%BFL(2,i_theta_lo,i_first_valid) + &
                                                w_theta*whamp_curve%BFL(2,i_theta_hi,i_first_valid)
                                    DP%BFL(3) = (1.0d0-w_theta)*whamp_curve%BFL(3,i_theta_lo,i_first_valid) + &
                                                w_theta*whamp_curve%BFL(3,i_theta_hi,i_first_valid)
                                    found = .true.
                                end if
                            end if
                        end block
                    end if

                    ! Targeted rescue via get_whamp_point
                    if (.not. found) then
                        block
                            real*8          :: n2_targeted, dw_targeted
                            complex(kind=8) :: EFL_targeted(3)
                            complex(kind=8) :: BFL_targeted(3)
                            logical         :: valid_targeted
                            call get_whamp_point(whamp_curve%fce_khz, whamp_curve%fpe_khz, &
                                                theta, abs(k_vec), root_omega, &
                                                n2_targeted, dw_targeted, EFL_targeted, &
                                                BFL_targeted, valid_targeted)
                            if (valid_targeted) then
                                n_ref_ind2 = n2_targeted
                                dw_dk      = dw_targeted
                                if (k_vec < 0.0d0) dw_dk = -dw_dk
                                DP%EFL(1)  = EFL_targeted(1)
                                DP%EFL(2)  = EFL_targeted(2)
                                DP%EFL(3)  = EFL_targeted(3)
                                DP%BFL(1) = BFL_targeted(1)
                                DP%BFL(2) = BFL_targeted(2)
                                DP%BFL(3) = BFL_targeted(3)
                                found = .true.
                            end if
                        end block
                    end if

                    if (.not. found) cycle   ! no valid WHAMP data at this point

                end block
              end if   ! whamp_curve%n_k > 0

              if (n_ref_ind2 <= 0.0d0) cycle
              if (abs(dw_dk)  < 1.0d-99) cycle

              k_vec  = sqrt(n_ref_ind2) * root_omega / c
              k_par  = k_sign(j) * k_vec * cos(theta)
              k_perp = k_vec * sin(theta)

              if (k_perp < 1.0d-30 .and. abs(nres) > 1) cycle

              if (k_par >= 0.0d0) then
                  sterm = vpar - dw_dk / cos(theta)
              else
                  sterm = vpar + dw_dk / cos(theta)
              end if

              if (abs(sterm) < 1.0d-99) cycle

              argJ = abs(k_perp) * p_perp / (me * abs(omega_e))
              call Bessel_FirstKind(argJ, nres-1, BesselF1)
              call Bessel_FirstKind(argJ, nres,   BesselF2)
              call Bessel_FirstKind(argJ, nres+1, BesselF3)

              call hot_plasma_polarisation(alpha_ne, nres, &
                                          BesselF1, BesselF2, BesselF3, &
                                          DP%EFL, DP%BFL, n_ref_ind2, phi_sq)

              hot_jacobian = (k_vec**2 / (1.0d0 + X**2)) * &
                            (1.0d0/(k_vec * dw_dk) + 1.0d0/root_omega)
              norm_factor = k_vec * X * abs(hot_jacobian) / (2.0d0*PI**2 * sqrt(1.0d0 + X**2))


            else

              ! Stix parameters at this root
              call cold_plasma_parameters(root_omega, omega_e, omega_pe, &
                                          omega_i, omega_pi, &
                                          stix_R, stix_L, stix_P, stix_S, stix_D)

              ! Refractive index squared (k_vec=0 since unused internally)
              call cold_plasma_refractive_index(root_omega, theta, &
                                                stix_P, stix_S, stix_R, stix_L, mu2)

              if (mu2 <= 0.0d0) cycle

              k_vec = sqrt(mu2) * root_omega / c
              k_par = (root_omega - b1) / vpar
              k_perp = k_vec * sin(theta)

              if (k_perp < 1.0d-30 .and. abs(nres) > 1) cycle
              ! Dispersion derivatives
              call dD_dk_dD_dw_cunningham(root_omega, theta, k_vec, mu2, &
                                          stix_R, stix_L, stix_P, stix_S, &
                                          omega_e, omega_pe, omega_i, omega_pi, &
                                          dD_dk, dD_dw)

              if (abs(dD_dw) < 1.0d-99) cycle

              if (k_par >= 0.0d0) then
                  sterm = vpar + (dD_dk/dD_dw)*sqrt(1.0d0 + X**2.0d0)
              else
                  sterm = vpar - (dD_dk/dD_dw)*sqrt(1.0d0 + X**2.0d0)
              end if

              if (abs(sterm) < 1.0d-99) cycle

              ! Phi squared (gyroresonance coupling factor)
              call phi_squared_cunningham(theta, k_par, k_perp, mu2, &
                                          stix_R, stix_L, stix_P, stix_S, &
                                          p_perp, alpha_ne, omega_e, nres, &
                                          phi_sq)

              call cunningham_jacobian(root_omega, k_vec, theta, mu2, &
                                      stix_R, stix_L, stix_P, stix_S, &
                                      omega_e, omega_pe, omega_i, omega_pi, &
                                      norm_factor)

              norm_factor = k_vec * X * abs(norm_factor) / (2.0d0*PI**2 * sqrt(1.0d0 + X**2))

            end if

            if (norm_factor < 1.0d-99) cycle

            ! Normalised power spectral density
            call PowerSpectralDensity(L, lambda, omega_e_eq, root_omega, &
                                      omega_m, omega_lc, omega_uc, d_omega, &
                                      PSD)

            ! Normalised intensity = Bw^2 * PSD * g(X)/integral_gX / norm_factor
            call get_normalised_intensity_cunningham(Bw, PSD, gX_raw, &
                                                      integral_gX, norm_factor, &
                                                      norm_intensity)

            ! DnX at this root
            term1 = (qe**2.0 * root_omega**2.0) / (4.0*PI*(1.0 + X**2.0) * p**2.0)
            term2 = (b1/root_omega - sin(alpha_ne)**2.0) / cos(alpha_ne)
            term3 = norm_intensity * phi_sq / abs(sterm)

            DnXaa = DnXaa + term1 * term2**2.0 * term3
            DnXap = DnXap + term1 * term2**2.0 * term3 * sin(alpha_ne)/term2
            DnXpp = DnXpp + term1 * term2**2.0 * term3 * (sin(alpha_ne)/term2)**2.0

        end do   ! root loop

        WeightValue(1) = DnXaa
        WeightValue(2) = DnXap
        WeightValue(3) = DnXpp

      end if  
                          
      WeightValue(1:3) = WeightValue(1:3) * sin(theta)
      
      if(i == 1 .or. i == nint_X) then
        DDN(1) = DDN(1) + WeightValue(1) / 2.0
        DDN(2) = DDN(2) + WeightValue(2) / 2.0
        DDN(3) = DDN(3) + WeightValue(3) / 2.0
      else
        DDN(1) = DDN(1) + WeightValue(1)
        DDN(2) = DDN(2) + WeightValue(2)
        DDN(3) = DDN(3) + WeightValue(3)
      end if
    end do
    
    DDN(1:3)=DDN(1:3)*BB3

end subroutine D_local_n
