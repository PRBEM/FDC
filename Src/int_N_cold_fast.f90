! ==============================================================================
! SUBROUTINE: int_N_cold_fast
! ==============================================================================
!
! PURPOSE:
!   Fast path for the COLD branch of int_N_omega: the normalisation integral
!
!     N(omega) = integral over theta of
!                sin(theta) * g(X) * n^2 * c / |d omega / dk|,   X = tan(theta)
!
! with the same trapezoidal rule, the same theta grid and the same formulas
! as the COLD branch of int_N_omega, reorganised for speed:
!
!   1. The theta grid and everything that depends only on it - sin(theta),
!      X^2 = tan^2(theta) and the wave-normal distribution g(X) (two exp()
!      calls) - is cached between calls. The grid depends only on theta_min,
!      the step (i.e. theta_max_eff) and the number of points, so the cache is
!      reused as long as the resonance cone does not truncate the domain
!      (res_cone = 'N'). With cone truncation the grid changes with omega and
!      the cache is simply rebuilt, at the old cost.
!   2. Everything that depends only on omega (powers of omega and c, products
!      of the Stix parameters and their derivatives) is computed once per
!      call, and all real-literal exponents (**4.0d0, **5.0d0, ...) are
!      replaced by multiplications.
!   3. The refractive index is evaluated inline with the formula of
!      cold_plasma_refractive_index (same branch: (B - sqrt(disc)) / 2A),
!      using the cached X^2 instead of calling tan() again.
!
!   -DINTN_ORIG switches int_N_omega back to its own loop and leaves this
!   routine unused, which is the way to check that the two still agree.
!
! CONSISTENCY:
!   The refractive index here must stay identical to
!   cold_plasma_refractive_index in Module_Dispersion (same branch of the
!   quadratic, same expression for the discriminant), and the trapezoidal
!   rule and the theta grid must stay identical to int_N_omega. If either
!   is changed, change both.
!
! REFERENCE:
!   The normalisation N(omega) is the one of Lyons (1974b), as used in the
!   LYONS formulation; see WeightFunction and int_N_omega. Note that it does
!   not converge when the resonance cone lies inside the integration range:
!   the integrand diverges there like (theta_cone - theta)^(-5/2), so the
!   value depends on where the integral is truncated (ConeTruncateN,
!   cone_backoff).
!
! INPUTS
!   omega                   wave frequency [rad/s]
!   theta0, step            first theta [rad] and grid step [rad]. The caller
!                           must compute these exactly as int_N_omega does
!                           (theta_min/180.0*PI and the step from
!                           theta_max_eff), so that the grids coincide.
!   nint                    number of theta points (nint_N_tmp)
!   X_m1, X_w1, X_m2, X_w2, theta_q   wave-normal distribution parameters
!   R, L, P, S              Stix parameters at omega
!   dPdw, dSdw, dRdw, dLdw  their omega-derivatives
!   c                       speed of light [m/s]
! OUTPUT
!   N_of_omega              the integral, already multiplied by step
! ==============================================================================
subroutine int_N_cold_fast(omega, theta0, step, nint, X_m1, X_w1, X_m2, X_w2, theta_q, &
                           R, L, P, S, dPdw, dSdw, dRdw, dLdw, c, N_of_omega)
implicit none
    real*8,  intent(in)  :: omega, theta0, step, X_m1, X_w1, X_m2, X_w2, theta_q
    real*8,  intent(in)  :: R, L, P, S, dPdw, dSdw, dRdw, dLdw, c
    integer, intent(in)  :: nint
    real*8,  intent(out) :: N_of_omega

    ! ---- theta-grid cache (kept between calls) ----
    real*8, allocatable, save :: sin_t(:), X2_t(:), g_t(:)
    integer, save :: ck_nint = -1
    real*8,  save :: ck_theta0 = 0.0d0, ck_step = 0.0d0
    real*8,  save :: ck_Xm1 = 0.0d0, ck_Xw1 = 0.0d0, ck_Xm2 = 0.0d0, ck_Xw2 = 0.0d0, ck_q = 0.0d0

    real*8  :: theta, X
    real*8  :: RL, PS, PRL, T3, LdR_RdL, PdS_SdP
    real*8  :: om2, om3, om5, c2, c4_om4, c2_om2
    real*8  :: X2, A, B, disc, n2, k, k3, kc2, dDdk, dF1, dF2, dF3, dwdk, Ntmp
    integer :: i

    ! ---- rebuild the cache if the grid or g(X) parameters changed ----
    if (nint /= ck_nint .or. theta0 /= ck_theta0 .or. step /= ck_step .or. &
        X_m1 /= ck_Xm1 .or. X_w1 /= ck_Xw1 .or. X_m2 /= ck_Xm2 .or. &
        X_w2 /= ck_Xw2 .or. theta_q /= ck_q) then
        if (allocated(sin_t)) then
            if (size(sin_t) < nint) deallocate(sin_t, X2_t, g_t)
        end if
        if (.not. allocated(sin_t)) allocate(sin_t(nint), X2_t(nint), g_t(nint))
        do i = 1, nint
            ! same expression as int_N_omega
            theta    = theta0 + step * (real(i) - 1.0)
            X        = tan(theta)
            sin_t(i) = sin(theta)
            X2_t(i)  = X * X
            g_t(i)   = exp(-((X - X_m1)/X_w1)**2) + theta_q * exp(-((X - X_m2)/X_w2)**2)
        end do
        ck_nint = nint; ck_theta0 = theta0; ck_step = step
        ck_Xm1 = X_m1; ck_Xw1 = X_w1; ck_Xm2 = X_m2; ck_Xw2 = X_w2; ck_q = theta_q
    end if

    ! ---- per-call (omega-only) constants ----
    RL      = R * L
    PS      = P * S
    PRL     = P * R * L
    LdR_RdL = L * dRdw + R * dLdw
    PdS_SdP = P * dSdw + S * dPdw
    T3      = R * L * dPdw + P * L * dRdw + P * R * dLdw
    om2     = omega * omega
    om3     = om2 * omega
    om5     = om3 * om2
    c2      = c * c
    c2_om2  = c2 / om2
    c4_om4  = c2_om2 * c2_om2

    ! ---- trapezoidal integration over the cached grid ----
    N_of_omega = 0.0d0
    do i = 1, nint
        X2 = X2_t(i)

        ! refractive index, as in cold_plasma_refractive_index
        A    = S * X2 + P
        B    = RL * X2 + (2.0d0 + X2) * PS
        disc = B * B - 4.0d0 * A * PRL * (1.0d0 + X2)
        n2   = (B - sqrt(disc)) / (2.0d0 * A)

        if (n2 < 0.0d0) then
            Ntmp = 0.0d0
        else
            k   = sqrt(n2) * omega / c
            k3  = k * k * k
            kc2 = (k * c) * (k * c)
            ! d(dispersion)/dk and d(dispersion)/d(omega), as in int_N_omega
            dDdk = A * c4_om4 * 4.0d0 * k3 - B * c2_om2 * 2.0d0 * k
            dF1  = kc2 * kc2 / om5 * (omega * X2 * dSdw + omega * dPdw - 4.0d0 * A)
            dF2  = kc2 / om3 * (omega * (X2 * LdR_RdL + (2.0d0 + X2) * PdS_SdP) - 2.0d0 * B)
            dF3  = (1.0d0 + X2) * T3
            dwdk = -(dDdk / (dF1 - dF2 + dF3))
            Ntmp = sin_t(i) * g_t(i) * n2 * c / abs(dwdk)
        end if

        if (i == 1 .or. i == nint) then
            N_of_omega = N_of_omega + 0.5d0 * Ntmp
        else
            N_of_omega = N_of_omega + Ntmp
        end if
    end do

    N_of_omega = N_of_omega * step

end subroutine int_N_cold_fast