! ==============================================================================
! FILE: Cunningham_normalization.f90
! ==============================================================================
!
! PURPOSE:
!   Collects the building-block subroutines needed to implement Cunningham
!   (2023)'s diffusion-coefficient formula (equivalently, Glauert & Horne
!   2005 Eq. 11-13 with Cunningham's local wave-power normalisation, Eq. 4),
!   as a self-contained alternative to the existing Albert (2007)-based
!   pipeline in WeightFunction.f90 / int_N_omega.f90.
!
!   These routines are direct Fortran transcriptions
!   of the reference Python implementation in RB-ENVIRONMENT/PIRAN.
!
! CONTENTS:
!   - dD_dk_dD_dw_cunningham            dispersion relation derivatives in
!                                       PIRAN's literal sin/cos Stix convention
!   - phi_squared_cunningham            Glauert & Horne (2005) Eq. A13
!                                       gyroresonance coupling factor
!   - gX_normalization_integral         integral of g(X) over the configured
!                                       X range (Simpson)
!   - get_normalised_intensity_cunningham   Cunningham (2023) Eq. 4b
!                                       normalised wave intensity
!
!   The remaining pieces of the pipeline live elsewhere: cunningham_jacobian
!   in Module_Dispersion.f90, the resonant roots from Roots_finder, and the
!   assembly of the kernel and the theta integration in the CUNNINGHAM branch
!   of D_local_n.f90, which calls the routines above.
!
! REFERENCES:
!   Cunningham, G. (2023), Resolution of a Few Problems in the Application of Quasilinear
!     Theory to Calculating Diffusion Coefficients in Heliophysics,
!     J. Geophys. Res. Space Physics.
!   Glauert, S. A., and R. B. Horne (2005), Calculation of pitch angle and
!     energy diffusion coefficients with the PADIE code, J. Geophys. Res.,
!     110, A04206, doi:10.1029/2004JA010851.
!   RB-ENVIRONMENT/PIRAN, https://github.com/RB-ENVIRONMENT/PIRAN
!     (Allanson, O., Kappas, T., Tyrrell, J., Garcia, A.)
!
! ==============================================================================

! ==============================================================================
! SUBROUTINE: dD_dk_dD_dw_cunningham
! ==============================================================================
!
! PURPOSE:
!   Computes dD/dk and dD/domega for the cold plasma dispersion relation
!   D = A*mu^4 - B*mu^2 + C = 0, where mu = c*k/omega, using the SAME
!   A, B, C convention as PIRAN:
!
!     A = S*sin^2(theta) + P*cos^2(theta)
!     B = R*L*sin^2(theta) + P*S*(2 + cos^2(theta))
!     C = P*R*L*(1 + cos^2(theta))
!
!   Formulas:
!
!     dD_dk = (2/k) * (2*A*mu^4 - B*mu^2)
!     dD_dw = (dA - 4*A/omega)*mu^4 - (dB - 2*B/omega)*mu^2 + dC
!
! INPUTS:
!   omega   - Wave frequency [rad/s]
!   theta   - Wave normal angle [radians]
!   k       - Wavenumber magnitude [rad/m]
!   mu2     - Refractive index squared at this (omega, theta) [-]
!   R,L,P,S - Stix parameters at this omega [-]
!   omega_e, omega_pe, omega_i, omega_pi - Local plasma frequencies [rad/s]
!
! OUTPUTS:
!   dD_dk - d(D)/d(k) at this point [units: 1/(rad/m), consistent with the
!           dimensionless D and k in rad/m]
!   dD_dw - d(D)/d(omega) at this point
!
! CALLS:
!   cold_plasma_stix_derivatives - existing FDC routine for dP/domega,
!                                  dS/domega, dR/domega, dL/domega
!
! ==============================================================================
subroutine dD_dk_dD_dw_cunningham(omega, theta, k, mu2, R, L, P, S, &
                                  omega_e, omega_pe, omega_i, omega_pi, &
                                  dD_dk, dD_dw)
use dispersion
implicit none
    real*8, intent(in)  :: omega, theta, k, mu2
    real*8, intent(in)  :: R, L, P, S
    real*8, intent(in)  :: omega_e, omega_pe, omega_i, omega_pi
    real*8, intent(out) :: dD_dk, dD_dw
 
    real*8 :: sin2, cos2
    real*8 :: A_, B_
    real*8 :: dPdw, dSdw, dRdw, dLdw
    real*8 :: dA_, dB_, dC_
    real*8 :: mu4
 
    sin2 = sin(theta)**2.0
    cos2 = cos(theta)**2.0
    mu4  = mu2**2.0
 
    ! A, B in piran's literal sin/cos convention (NOT FDC's X^2 convention)
    A_ = S*sin2 + P*cos2
    B_ = R*L*sin2 + P*S*(2.0 + cos2)
 
    call cold_plasma_stix_derivatives(omega, omega_e, omega_pe, omega_i, omega_pi, &
                                      dPdw, dSdw, dRdw, dLdw)
 
    dA_ = dSdw*sin2 + dPdw*cos2
    dB_ = (dRdw*L + R*dLdw)*sin2 + (dPdw*S + P*dSdw)*(2.0 + cos2)
    dC_ = dPdw*R*L + P*dRdw*L + P*R*dLdw
 
    dD_dk = (2.0/k) * (2.0*A_*mu4 - B_*mu2)
    dD_dw = (dA_ - 4.0*A_/omega)*mu4 - (dB_ - 2.0*B_/omega)*mu2 + dC_
 
end subroutine dD_dk_dD_dw_cunningham


! ==============================================================================
! SUBROUTINE: phi_squared_cunningham
! ==============================================================================
!
! PURPOSE:
!   Computes Phi_{n,k}^2 (Glauert & Horne 2005, Eq. A13).
!   Could possibly be merged with cold_plasma_polarisation in Module_Dispersion
!   in the future.
!
!     point = k_perp * p_perp / (mass * gyrofreq)     [gyrofreq SIGNED]
!     J(-1+nres), J(nres), J(1+nres) at this point
!     term1 = ((mu2-L)/(mu2-S))*J(nres+1) + ((mu2-R)/(mu2-S))*J(nres-1)
!     term2 = (mu2*sin(psi)^2 - P) / (2*mu2)
!     term3 = (1/tan(alpha)) * sin(psi)*cos(psi) * J(nres)
!     term4 = ((R-L)/(2*(mu2-S)))^2
!     term5 = ((P - mu2*sin(psi)^2)/mu2)^2
!     term6 = (P*cos(psi)/mu2)^2
!     phi_squared = (term1*term2 + term3)^2 / (term4*term5 + term6)
!
! INPUTS:
!   theta    - Wave normal angle [radians]
!   k_par    - Signed parallel wavenumber [rad/m] (determines psi branch)
!   k_perp   - Perpendicular wavenumber [rad/m] (non-negative)
!   mu2      - Refractive index squared [-]
!   R,L,P,S  - Stix parameters at this omega [-]
!   p_perp   - Perpendicular momentum of the resonant electron [kg m/s]
!   alpha_ne - Local pitch angle [radians]
!   sqe      - Sign of the electron charge (FDC module constant, -1)
!   omega_e  - Electron gyrofrequency magnitude [rad/s]
!   me       - Electron rest mass [kg]
!   nres     - Resonance harmonic number (integer)
!
! OUTPUTS:
!   phi_sq   - Phi_{n,k}^2 [dimensionless]
!
! CALLS:
!   Bessel_FirstKind
!
! ==============================================================================
subroutine phi_squared_cunningham(theta, k_par, k_perp, mu2, R, L, P, S, &
                                  p_perp, alpha_ne, omega_e, nres, &
                                  phi_sq)
use constants   
implicit none
    real*8, intent(in)  :: theta, k_par, k_perp, mu2
    real*8, intent(in)  :: R, L, P, S
    real*8, intent(in)  :: p_perp, alpha_ne
    real*8, intent(in)  :: omega_e
    integer, intent(in) :: nres
    real*8, intent(out) :: phi_sq
 
    real*8 :: sin_psi, cos_psi
    real*8 :: point, point_nm1, point_n, point_np1
    real*8 :: J_nm1, J_n, J_np1
    real*8 :: term1, term2, term3, term4, term5, term6
    real*8 :: gyrofreq_signed
 
    ! psi differs from theta by sign of cos() when k_par < 0
    sin_psi = sin(theta)
    if (k_par >= 0.0d0) then
        cos_psi = cos(theta)
    else
        cos_psi = -cos(theta)
    end if
 
    gyrofreq_signed = real(sqe, 8) * omega_e
 
    point = k_perp * p_perp / (me * gyrofreq_signed)
 
    point_nm1 = point
    point_n   = point
    point_np1 = point
 
    call Bessel_FirstKind(point_nm1, nres - 1, J_nm1)
    call Bessel_FirstKind(point_n,   nres,     J_n)
    call Bessel_FirstKind(point_np1, nres + 1, J_np1)
 
    term1 = ((mu2 - L)/(mu2 - S))*J_np1 + ((mu2 - R)/(mu2 - S))*J_nm1
    term2 = (mu2*sin_psi**2.0 - P) / (2.0*mu2)
    term3 = (cos(alpha_ne)/sin(alpha_ne)) * sin_psi*cos_psi * J_n
    term4 = ((R - L)/(2.0*(mu2 - S)))**2.0
    term5 = ((P - mu2*sin_psi**2.0)/mu2)**2.0
    term6 = (P*cos_psi/mu2)**2.0
 
    phi_sq = (term1*term2 + term3)**2.0 / (term4*term5 + term6)

end subroutine phi_squared_cunningham

! ==============================================================================
! SUBROUTINE: gX_normalization_integral
! ==============================================================================
!
! PURPOSE:
!   Computes integral(g(X) dX) over the configured wave-normal-angle range
!   [theta_min, theta_max] (converted to X = tan(theta)), where g(X) is
!   the wave-normal-angle distribution:
!
!     g(X) = X * [ exp(-((X-X_m1)/X_w1)^2)
!                  + theta_q*exp(-((X-X_m2)/X_w2)^2) ]
!
!   The leading factor X enforces the boundary condition g(X=0) = 0; the
!   same form is used for gX_raw in the CUNNINGHAM branch of D_local_n.
!
!   This integral is required by Cunningham (2023)'s formulation: under
!   this method, g(X) must be PRE-NORMALISED by its own integral over the
!   full configured X range BEFORE being used in the diffusion coefficient
!   formula (i.e. the quantity actually used is g(X)/integral(g(X)dX), not
!   raw g(X) as the existing LYONS/Albert path uses).
!
!   Theta_function (and hence g(X)) CAN vary with magnetic
!   latitude lambda via Theta_wave(L, lambda, Theta_function), so this
!   integral must be recomputed once per (alpha, E, L, lambda) evaluation
!   (i.e. once per call to the Cunningham equivalent of D_local_n), NOT
!   cached globally for the whole run -- but it does NOT depend on omega,
!   theta (the integration variable), or the resonant root, so it should
!   be computed ONCE per such call and reused across all theta points and
!   resonant roots within that call, not recomputed per-point.
!
! INPUTS:
!   X_m1, X_w1    - Centre and width of the first Gaussian lobe, in
!                   X = tan(theta) (the caller converts from degrees)
!   X_m2, X_w2    - Centre and width of the second Gaussian lobe
!   theta_q       - Relative amplitude of the second Gaussian lobe [-]
!   X_min, X_max  - Integration bounds in X
!   nx_pts        - Number of integration points (>= 2); Simpson's rule, with
!                   a trapezoidal last interval if the count is even
!
! OUTPUTS:
!   integral_gX - integral of g(X) dX over [X_min, X_max] [-]
!
! ==============================================================================
subroutine gX_normalization_integral(X_m1, X_w1, X_m2, X_w2, theta_q, &
                                     X_min, X_max, nx_pts, integral_gX)
implicit none
    real*8, intent(in)  :: X_m1, X_w1, X_m2, X_w2, theta_q
    real*8, intent(in)  :: X_min, X_max
    integer, intent(in) :: nx_pts
    real*8, intent(out) :: integral_gX
 
    real*8 :: dX, X
    real*8, allocatable :: gX(:)
    integer :: i, n_intervals
 
    n_intervals = nx_pts - 1
    dX = (X_max - X_min) / real(n_intervals, 8)
 
    allocate(gX(nx_pts))
    do i = 1, nx_pts
        X = X_min + real(i-1, 8)*dX
        gX(i) = X * ( exp(-((X - X_m1)/X_w1)**2) &
                + theta_q * exp(-((X - X_m2)/X_w2)**2) )
    end do
 
    if (mod(n_intervals, 2) == 0) then
        integral_gX = gX(1) + gX(nx_pts)
        do i = 2, nx_pts-1, 2
            integral_gX = integral_gX + 4.0d0*gX(i)
        end do
        do i = 3, nx_pts-2, 2
            integral_gX = integral_gX + 2.0d0*gX(i)
        end do
        integral_gX = integral_gX * dX / 3.0d0
    else
        integral_gX = gX(1) + gX(nx_pts-1)
        do i = 2, nx_pts-2, 2
            integral_gX = integral_gX + 4.0d0*gX(i)
        end do
        do i = 3, nx_pts-3, 2
            integral_gX = integral_gX + 2.0d0*gX(i)
        end do
        integral_gX = integral_gX * dX / 3.0d0
        integral_gX = integral_gX + 0.5d0*dX*(gX(nx_pts-1) + gX(nx_pts))
    end if
 
    deallocate(gX)
 
end subroutine gX_normalization_integral

! ==============================================================================
! SUBROUTINE: get_normalised_intensity_cunningham
! ==============================================================================
!
! PURPOSE:
!   Computes the normalised wave intensity |B_k^norm|^2 (Cunningham 2023
!   Eq. 4b), used inside DnX_single_root_cunningham's term3. 
!
!     normalised_intensity = power_spectral_density * eval_gx / norm_factor
!   where (for Cunningham's method) eval_gx = g(X) / integral(g(X)dX), and
!   power_spectral_density = wave_amplitude^2 * [normalised lineshape].
!
! INPUTS:
!   Bw        - RMS wave magnetic field amplitude [T] (obtained from Bw_wave)
!   PSD       - FDC's normalised power spectral density at this omega,
!               from PowerSpectralDensity [s/rad]
!   gX_raw    - Unnormalised g(X) at this X, i.e.
!               X*[exp(-((X-X_m1)/X_w1)^2) + theta_q*exp(-((X-X_m2)/X_w2)^2)]
!   integral_gX - integral(g(X)dX) over the configured X range, from
!               gX_normalization_integral
!   norm_factor - Cunningham (2023) Eq. 4b local Jacobian normalisation,
!               from cunningham_jacobian (Module_Dispersion.f90)
!               [rad*s/m^2]
!
! OUTPUTS:
!   normalised_intensity - |B_k^norm|^2 [T^2*m^2/rad^2]
!
! ==============================================================================
subroutine get_normalised_intensity_cunningham(Bw, PSD, gX_raw, integral_gX, &
                                               norm_factor, normalised_intensity)
implicit none
    real*8, intent(in)  :: Bw, PSD, gX_raw, integral_gX, norm_factor
    real*8, intent(out) :: normalised_intensity
 
    real*8 :: power_spectral_density, eval_gX
 
    power_spectral_density = (Bw**2) * PSD
    eval_gX = gX_raw / integral_gX
 
    normalised_intensity = power_spectral_density * eval_gX / abs(norm_factor)
 
end subroutine get_normalised_intensity_cunningham

