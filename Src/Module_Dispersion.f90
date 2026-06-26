! ==============================================================================
! MODULE: dispersion
! ==============================================================================
!
! PURPOSE:
!   Provides a unified interface for evaluating cold or hot plasma dispersion
!   properties at a given wave frequency, wave normal angle, and local plasma
!   state. The dispersion model is selected at runtime via the DispersionType
!   flag in module constants ('COLD' or 'HOT').
!
!   Currently implemented:
!     'COLD' - Cold plasma dispersion relation (Stix formulation)
!
!   Planned but not yet implemented:
!     'HOT'  - Hot plasma dispersion relation via WHAMP
!
!   This module is intended to be the single point of change when adding
!   new dispersion models. WeightFunction and any other caller should obtain
!   all dispersion quantities exclusively through Get_dispersion_properties
!   and Get_polarisation_factor, and should have no direct knowledge of
!   which model is active.
!
! CONTENTS:
!   DispersionProperties        - Derived type carrying dispersion outputs
!   Get_dispersion_properties   - Dispatcher: selects cold or hot plasma path
!   Get_polarisation_factor     - Computes Fsq from a DispersionProperties instance
!   cold_plasma_parameters      - Computes Stix parameters R, L, P, S, D
!   cold_plasma_group_velocity  - Computes dw/dk from analytic Stix derivatives
!   cold_plasma_polarisation    - Computes Fsq from Stix parameters and Bessel functions
!
! MODULE DEPENDENCIES:
!   constants  - Physical constants and DispersionType flag
!
! NOTES:
!   - This module replaces the inline Stix parameter, group velocity, and
!     polarisation factor calculations that were previously embedded in
!     WeightFunction, and the a1-a6 side effects previously written to
!     module variables by Roots_finder.
!   - module variables (R_stix, L_stix, P_stix, S_stix, D_stix) are no
!     longer written as side effects after this refactoring. Any subroutine
!     that previously read these from module variables (e.g. int_N_omega)
!     will need to be updated in a future refactoring step.
!
! ==============================================================================
 
module dispersion
use constants
implicit none
 
! ==============================================================================
! DERIVED TYPE: DispersionProperties
! ==============================================================================
!
! PURPOSE:
!   Carries all dispersion relation outputs for a single evaluation point
!   (omega, k, theta, local plasma state). Passed from Get_dispersion_properties
!   to Get_polarisation_factor, avoiding repeated computation.
!
! FIELDS:
!   disp_type   - Which model populated this instance ('COLD' or 'HOT')
!   n_ref_ind2  - Refractive index squared (n^2) [-]
!   dw_dk       - Group velocity along wave vector dw/dk [m/s]
!
!   Cold plasma fields (populated when disp_type = 'COLD'):
!     R, L, P, S, D  - Stix cold plasma parameters [-]
!
!   Hot plasma fields (populated when disp_type = 'HOT'):
!     Fields to be determined once WHAMP interface is finalised.
!     Commented out below to document intent without committing to an
!     interface that may change.
!
! ==============================================================================
type :: DispersionProperties
    character(len=10) :: disp_type
 
    ! Common outputs
    real*8 :: n_ref_ind2    ! refractive index squared
    real*8 :: dw_dk         ! group velocity along wave vector [m/s]
 
    ! Cold plasma (Stix) outputs
    real*8 :: R             ! Stix R parameter
    real*8 :: L             ! Stix L parameter
    real*8 :: P             ! Stix P parameter
    real*8 :: S             ! Stix S = (R+L)/2
    real*8 :: D             ! Stix D = (R-L)/2
 
    ! Hot plasma (WHAMP) outputs (to be implemented)
    complex(kind=8) :: EFL(3)   ! E field components from WHAMP
    ! real*8 :: bx, by, bz  ! normalised magnetic field components
 
end type DispersionProperties

! ==============================================================================
! DERIVED TYPE: WHAMPCurve
! ==============================================================================
!
! PURPOSE:
!   Holds the hot plasma dispersion surface computed by WHAMP for a single
!   local plasma condition (fce, fpe). Computed once per (L, lambda) point
!   in Integration_lambda and passed down the call chain to WeightFunction
!   and int_N_omega.
!
!   The surface represents omega(k, theta), the dispersion relation as a
!   function of total wave vector magnitude k and wave normal angle theta.
!   WHAMP is called once per theta point, scanning k at fixed theta by
!   setting k_perp = k*sin(theta) and k_par = k*cos(theta).
!
!   At runtime, for each theta point in the FDC integration, the code
!   searches the 1D slice omega(k) at that theta to find the resonant k
!   satisfying the cyclotron resonance condition:
!     omega(k, theta) - k*cos(theta)*v_par - n*Omega_e/gamma = 0
!
! FIELDS:
!   n_k       - Number of k points per theta slice
!   n_theta   - Number of theta points
!   fce_khz   - Electron gyrofrequency [kHz]
!   fpe_khz   - Electron plasma frequency [kHz]
!   theta     - Wave normal angle grid [radians], size (n_theta)
!   k         - Total wave vector magnitude in SI units [rad/m], size (n_k)
!   k_whamp   - k in WHAMP normalised units (k*c/Omega_e), size (n_k)
!   omega     - Wave frequency [rad/s], size (n_theta, n_k)
!   dw_dk     - Group velocity dω/dk [m/s], size (n_theta, n_k)
!   valid     - True if WHAMP converged at this (theta, k) point,
!               size (n_theta, n_k)
!
! NOTES:
!   - n_k and n_theta are set by WHAMP_N_K and WHAMP_N_THETA in module
!     constants
!   - Only points where valid=.true. should be used for interpolation
!   - For parallel propagation only (theta=0), n_theta=1
!   - fce_khz and fpe_khz stored for reference/debugging
!
! ==============================================================================
type :: WHAMPCurve
    integer :: n_k              ! number of k points per theta slice
    integer :: n_theta          ! number of theta points
    real*8  :: fce_khz          ! electron gyrofrequency [kHz]
    real*8  :: fpe_khz          ! electron plasma frequency [kHz]
 
    real*8,  allocatable :: theta(:)      ! wave normal angle [rad], (n_theta)
    real*8,  allocatable :: k(:)          ! wave vector SI [rad/m], (n_k)
    real*8,  allocatable :: k_whamp(:)    ! wave vector WHAMP units, (n_k)
    real*8,  allocatable :: omega(:,:)    ! frequency [rad/s], (n_theta, n_k)
    real*8,  allocatable :: dw_dk(:,:)    ! group velocity [m/s], (n_theta, n_k)
    logical, allocatable :: valid(:,:)    ! convergence flag, (n_theta, n_k)
    complex(kind=8), allocatable :: EFL(:,:,:) ! E field components (3, n_theta, n_k)
end type WHAMPCurve

 
contains
 
 
! ==============================================================================
! SUBROUTINE: Get_dispersion_properties
! ==============================================================================
!
! PURPOSE:
!   Dispatcher routine. Computes the refractive index, group velocity, and
!   polarisation-related dispersion properties at a given wave frequency
!   (omega), wave vector magnitude (k_vec), wave normal angle (theta), and
!   local plasma state. Branches on DispersionType from module constants.
!
! INPUTS:
!   omega       - Wave frequency [rad/s]
!   k_vec       - Wave vector magnitude [rad/m]
!   theta       - Wave normal angle [radians]
!   omega_e     - Electron gyrofrequency [rad/s]
!   omega_pe    - Electron plasma frequency [rad/s]
!   omega_i     - Ion gyrofrequency [rad/s]
!   omega_pi    - Ion plasma frequency [rad/s]
!
! OUTPUTS:
!   DP          - DispersionProperties instance populated by the chosen model
!
! NOTES:
!   - 'HOT' case currently halts with an informative error. It will be
!     implemented once the WHAMP interface is finalised.
!
! ==============================================================================
subroutine Get_dispersion_properties(omega, k_vec, theta, &
                                     omega_e, omega_pe, omega_i, omega_pi, &
                                     DP)
use constants
implicit none
    real*8, intent(in)  :: omega, k_vec, theta
    real*8, intent(in)  :: omega_e, omega_pe, omega_i, omega_pi
    type(DispersionProperties), intent(out) :: DP
 
    select case(trim(DispersionType))
 
        case('COLD')
            DP%disp_type = 'COLD'
            call cold_plasma_parameters(omega, omega_e, omega_pe, omega_i, omega_pi, &
                                        DP%R, DP%L, DP%P, DP%S, DP%D)
            call cold_plasma_refractive_index(omega, k_vec, theta, &
                                              DP%P, DP%S, DP%R, DP%L, &
                                              DP%n_ref_ind2)
            call cold_plasma_group_velocity(omega, k_vec, theta, &
                                            DP%R, DP%L, DP%P, DP%S, &
                                            omega_e, omega_pe, omega_i, omega_pi, &
                                            DP%dw_dk)
 
        case('HOT')
            ! HOT plasma dispersion via WHAMP.
            ! Stix parameters (R, L, P, S, D) are computed from cold plasma
            ! as an approximation for the polarisation factor and resonance
            ! cone check. n_ref_ind2 and dw_dk are supplied by the WHAMP
            ! curve via WeightFunction and are not computed here.
            DP%disp_type = 'HOT'
            call cold_plasma_parameters(omega, omega_e, omega_pe, omega_i, omega_pi, &
                                        DP%R, DP%L, DP%P, DP%S, DP%D)
            ! n_ref_ind2 and dw_dk are set in WeightFunction from the WHAMP
            ! curve interpolation, not here.
            DP%n_ref_ind2 = 0.0d0
            DP%dw_dk      = 0.0d0
 
        case default
            write(*,*) 'ERROR: Unknown DispersionType = ', trim(DispersionType)
            write(*,*) 'Valid options are: COLD, HOT'
            stop
 
    end select
 
end subroutine Get_dispersion_properties
 
 
! ==============================================================================
! SUBROUTINE: Get_polarisation_factor
! ==============================================================================
!
! PURPOSE:
!   Computes the wave polarisation factor Fsq for a given dispersion model,
!   electron state, and Bessel function values. Branches on DP%disp_type.
!
!   Fsq weights the resonance broadening term Delta_n by the relative
!   amplitude of the wave field components driving the resonance. In the
!   cold plasma case it is derived from the Stix parameters and Bessel
!   functions J_{n-1}, J_n, J_{n+1}. In the hot plasma case it will be
!   derived directly from the WHAMP field components.
!
! INPUTS:
!   alpha_ne    - Local electron pitch angle [radians]
!   nres        - Resonance harmonic number [-]
!   theta       - Wave normal angle [radians]
!   k_vec       - Wave vector magnitude [rad/m]
!   omega       - Wave frequency [rad/s]
!   BesselF1    - Bessel function J_{n-1}(argJ) [-]
!   BesselF2    - Bessel function J_n(argJ)   [-]
!   BesselF3    - Bessel function J_{n+1}(argJ) [-]
!   DP          - DispersionProperties instance from Get_dispersion_properties
!
! OUTPUTS:
!   Fsq         - Polarisation factor [-]
!
! ==============================================================================
subroutine Get_polarisation_factor(alpha_ne, nres, theta, k_vec, omega, &
                                   BesselF1, BesselF2, BesselF3, &
                                   DP, Fsq)
implicit none
    real*8, intent(in)  :: alpha_ne, theta, k_vec, omega
    integer, intent(in) :: nres
    real*8, intent(in)  :: BesselF1, BesselF2, BesselF3
    type(DispersionProperties), intent(in) :: DP
    real*8, intent(out) :: Fsq
 
    select case(trim(DP%disp_type))
 
        case('COLD')
            call cold_plasma_polarisation(alpha_ne, nres, theta, k_vec, omega, &
                                          BesselF1, BesselF2, BesselF3, &
                                          DP%n_ref_ind2, DP%R, DP%L, DP%P, DP%S, &
                                          Fsq)
        case('HOT')
            call hot_plasma_polarisation(alpha_ne, nres, &
                                         BesselF1, BesselF2, BesselF3, &
                                         DP%EFL, Fsq)
 
        case default
            write(*,*) 'ERROR: Unknown disp_type in DispersionProperties: ', &
                        trim(DP%disp_type)
            stop
 
    end select
 
end subroutine Get_polarisation_factor
 
 
! ==============================================================================
! SUBROUTINE: cold_plasma_parameters
! ==============================================================================
!
! PURPOSE:
!   Computes the Stix cold plasma parameters R, L, P, S, D for a two-component
!   plasma (electrons + one ion species) at a given wave frequency omega.
!
!   Definitions (Stix 1992):
!     R = 1 - sum[ omega_ps^2 / (omega * (omega + sqe_s * omega_cs)) ]
!     L = 1 - sum[ omega_ps^2 / (omega * (omega - sqe_s * omega_cs)) ]
!     P = 1 - sum[ omega_ps^2 / omega^2 ]
!     S = (R + L) / 2
!     D = (R - L) / 2
!
! INPUTS:
!   omega       - Wave frequency [rad/s]
!   omega_e     - Electron gyrofrequency [rad/s]
!   omega_pe    - Electron plasma frequency [rad/s]
!   omega_i     - Ion gyrofrequency [rad/s]
!   omega_pi    - Ion plasma frequency [rad/s]
!
! OUTPUTS:
!   R, L, P, S, D  - Stix parameters [-]
!
! ==============================================================================
subroutine cold_plasma_parameters(omega, omega_e, omega_pe, omega_i, omega_pi, &
                                  R, L, P, S, D)
implicit none
    real*8, intent(in)  :: omega, omega_e, omega_pe, omega_i, omega_pi
    real*8, intent(out) :: R, L, P, S, D
 
    R = 1.0 - omega_pe**2.0 / (omega * (omega + real(sqe) * omega_e)) &
            - omega_pi**2.0 / (omega * (omega + omega_i))
 
    L = 1.0 - omega_pe**2.0 / (omega * (omega - real(sqe) * omega_e)) &
            - omega_pi**2.0 / (omega * (omega - omega_i))
 
    P = 1.0 - (omega_pe / omega)**2.0 &
            - (omega_pi / omega)**2.0
 
    S = 0.5 * (R + L)
    D = 0.5 * (R - L)
 
end subroutine cold_plasma_parameters
 
 
! ==============================================================================
! SUBROUTINE: cold_plasma_refractive_index
! ==============================================================================
!
! PURPOSE:
!   Computes the cold plasma refractive index squared (n^2) from the Stix
!   parameters and the wave normal angle theta, using the Appleton-Hartree
!   form of the cold plasma dispersion relation.
!
!   The lower-sign (whistler-mode) branch is used:
!     n^2 = (R*L*X^2 + (2+X^2)*P*S - sqrt(discriminant)) / (2*(S*X^2 + P))
!   where X = tan(theta).
!
!   This corresponds to the branch selected in the original WeightFunction
!   (Xudong revision, 2011-10-24).
!
! INPUTS:
!   omega       - Wave frequency [rad/s]  (unused directly; for interface
!                 consistency and future use)
!   k_vec       - Wave vector magnitude [rad/m]
!   theta       - Wave normal angle [radians]
!   P, S, R, L  - Stix parameters [-]
!
! OUTPUTS:
!   n_ref_ind2  - Refractive index squared [-]
!
! NOTES:
!   - Only the lower-sign branch is implemented, appropriate for whistler-mode
!     waves. If other wave modes are needed in future, this subroutine should
!     be extended with a branch selector argument.
!
! ==============================================================================
subroutine cold_plasma_refractive_index(omega, k_vec, theta, P, S, R, L, &
                                        n_ref_ind2)
implicit none
    real*8, intent(in)  :: omega, k_vec, theta, P, S, R, L
    real*8, intent(out) :: n_ref_ind2
 
    real*8 :: X
    real*8 :: numerator_linear, discriminant
 
    X = tan(theta)
 
    numerator_linear = R*L*(X**2.0) + (2.0 + X**2.0)*P*S
    discriminant     = numerator_linear**2.0 &
                       - 4.0*(S*(X**2.0) + P)*P*R*L*(1.0 + X**2.0)
 
    ! Lower-sign (whistler-mode) branch
    n_ref_ind2 = (numerator_linear - sqrt(discriminant)) &
                 / (2.0 * (S*(X**2.0) + P))
 
end subroutine cold_plasma_refractive_index
 
 
! ==============================================================================
! SUBROUTINE: cold_plasma_group_velocity
! ==============================================================================
!
! PURPOSE:
!   Computes the cold plasma group velocity dw/dk along the wave vector
!   direction, using analytic derivatives of the Stix parameters with
!   respect to omega.
!
!   The group velocity is derived from the implicit dispersion relation D=0:
!     dw/dk = -(dD/dk) / (dD/dw)
!   where D is the cold plasma dispersion function expressed in terms of
!   k and omega.
!
!   The analytic expressions for dP/dw, dS/dw, dR/dw, dL/dw follow
!   Albert (2007), equations for the derivative of the dispersion tensor.
!
! INPUTS:
!   omega       - Wave frequency [rad/s]
!   k_vec       - Wave vector magnitude [rad/m]
!   theta       - Wave normal angle [radians]
!   R, L, P, S  - Stix parameters [-]
!   omega_e     - Electron gyrofrequency [rad/s]
!   omega_i     - Ion gyrofrequency [rad/s]
!
! OUTPUTS:
!   dw_dk       - Group velocity along wave vector [m/s]
!
! NOTES:
!   - The coefficients a1-a6 are computed locally here from the plasma
!     frequencies. They are no longer written to module variables (this
!     replaces the side effect previously set in Roots_finder).
!   - a5 is set to 0.0, consistent with the original implementation.
!     See note in Roots_finder header regarding whether this is a valid
!     simplification.
!
! ==============================================================================
subroutine cold_plasma_group_velocity(omega, k_vec, theta, R, L, P, S, &
                                      omega_e, omega_pe, omega_i, omega_pi, &
                                      dw_dk)
implicit none
    real*8, intent(in)  :: omega, k_vec, theta, R, L, P, S
    real*8, intent(in)  :: omega_e, omega_pe, omega_i, omega_pi
    real*8, intent(out) :: dw_dk
 
    real*8 :: X
 
    ! Local polynomial coefficients (previously a1-a6 in module variables)
    real*8 :: a1, a2, a3, a4, a5, a6
 
    ! Analytic Stix parameter derivatives w.r.t. omega
    real*8 :: dPdw, dSdw, dRdw, dLdw
 
    ! Dispersion function derivatives
    real*8 :: dF1dw, dF2dw, dF3dw, dDdw, dDdk
 
    X = tan(theta)
 
    ! Plasma frequency combination coefficients
    ! (definitions match those set in Roots_finder)
    a1 = omega_i**2.0 + omega_e**2.0 + omega_pe**2.0 + omega_pi**2.0
    a2 = (omega_i*omega_e)**2.0 + (omega_pe*omega_i)**2.0 + (omega_pi*omega_e)**2.0
    a3 = omega_i + real(sqe)*omega_e
    a4 = omega_i*real(sqe)*omega_e - omega_pe**2.0 - omega_pi**2.0
    a5 = 0.0   ! analytically zero; see note in header regarding Albert (2007)
    a6 = omega_pi**2.0 + omega_pe**2.0
 
    ! Analytic derivatives of Stix parameters w.r.t. omega
    dPdw = 2.0 * a6 / (omega**3.0)
 
    dSdw = (2.0*(a1 - omega_i**2.0 - omega_e**2.0)*(omega**5.0) + &
             4.0*((omega_i**2.0)*(omega_e**2.0) - a2)*(omega**3.0) + &
             2.0*(a2*(omega_i**2.0) + a2*(omega_e**2.0) - &
                  a1*(omega_i**2.0)*(omega_e**2.0))*omega) / &
            (((omega**2.0 - omega_i**2.0)**2.0) * &
             ((omega**2.0 - omega_e**2.0)**2.0))
 
    dRdw = (2.0*(omega_i*real(sqe)*omega_e - a4)*(omega**3.0) + &
             (a3*omega_i*real(sqe)*omega_e - a4*omega_i - &
              a4*real(sqe)*omega_e - 3.0*a5)*(omega**2.0) - &
             2.0*a5*(omega_i + real(sqe)*omega_e)*omega - &
             a5*omega_i*real(sqe)*omega_e) / &
            ((omega**2.0) * ((omega + omega_i)**2.0) * &
             ((omega + real(sqe)*omega_e)**2.0))
 
    dLdw = (2.0*(omega_i*real(sqe)*omega_e - a4)*omega**3.0 - &
             (a3*omega_i*real(sqe)*omega_e - a4*omega_i - &
              a4*real(sqe)*omega_e - 3.0*a5)*(omega**2.0) - &
             2.0*a5*(omega_i + real(sqe)*omega_e)*omega + &
             a5*omega_i*real(sqe)*omega_e) / &
            ((omega**2.0) * ((omega - omega_i)**2.0) * &
             ((omega - real(sqe)*omega_e)**2.0))
 
    ! dD/dk: derivative of dispersion function w.r.t. k
    dDdk = (S*(X**2.0) + P) * (c**4.0) / (omega**4.0) * 4.0*(k_vec**3.0) &
           - (R*L*(X**2.0) + (2.0 + X**2.0)*P*S) * (c**2.0)/(omega**2.0) &
             * 2.0*k_vec
 
    ! dD/dw: derivative of dispersion function w.r.t. omega
    dF1dw = ((k_vec*c)**4.0) / (omega**5.0) * &
            (omega*(X**2.0)*dSdw + omega*dPdw - 4.0*(S*(X**2.0) + P))
 
    dF2dw = ((k_vec*c)**2.0) / (omega**3.0) * &
            (omega*((X**2.0)*(L*dRdw + R*dLdw) + &
                    (2.0 + X**2.0)*(P*dSdw + S*dPdw)) - &
             2.0*(R*L*(X**2.0) + (2.0 + X**2.0)*P*S))
 
    dF3dw = (1.0 + X**2.0) * (R*L*dPdw + P*L*dRdw + P*R*dLdw)
 
    dDdw = dF1dw - dF2dw + dF3dw
 
    dw_dk = -dDdk / dDdw
 
end subroutine cold_plasma_group_velocity
 
 
! ==============================================================================
! SUBROUTINE: cold_plasma_polarisation
! ==============================================================================
!
! PURPOSE:
!   Computes the cold plasma wave polarisation factor Fsq, which appears in
!   the quasi-linear diffusion coefficient integrand as a weight on the
!   resonance broadening term Delta_n.
!
!   Fsq encodes how strongly each wave field component drives the resonant
!   interaction for a given electron pitch angle, resonance harmonic, and
!   wave normal angle. It is constructed from Bessel functions J_{n-1}, J_n,
!   J_{n+1} and combinations of Stix parameters following Albert (2007).
!
!   The factor is:
!     Fsq = | (F1*F2 + F3)^2 / ((F4*F5)^2 + F6^2) |
!
!   where F1-F6 are intermediate combinations of Stix parameters, Bessel
!   functions, and wave normal angle (defined in code below and in
!   Albert 2007, equations [A1]-[A6]).
!
! INPUTS:
!   alpha_ne    - Local electron pitch angle [radians]
!   nres        - Resonance harmonic number [-]
!   theta       - Wave normal angle [radians]
!   k_vec       - Wave vector magnitude [rad/m]
!   omega       - Wave frequency [rad/s]
!   BesselF1    - J_{n-1}(argJ) [-]
!   BesselF2    - J_n(argJ)     [-]
!   BesselF3    - J_{n+1}(argJ) [-]
!   n_ref_ind2  - Refractive index squared [-]
!   R, L, P, S  - Stix parameters [-]
!
! OUTPUTS:
!   Fsq         - Polarisation factor [-]
!
! ==============================================================================
subroutine cold_plasma_polarisation(alpha_ne, nres, theta, k_vec, omega, &
                                    BesselF1, BesselF2, BesselF3, &
                                    n_ref_ind2, R, L, P, S, &
                                    Fsq)
implicit none
    real*8, intent(in)  :: alpha_ne, theta, k_vec, omega
    integer, intent(in) :: nres
    real*8, intent(in)  :: BesselF1, BesselF2, BesselF3
    real*8, intent(in)  :: n_ref_ind2, R, L, P, S
    real*8, intent(out) :: Fsq
 
    real*8 :: F1, F2, F3, F4, F5, F6
 
    F1 = BesselF3 * (n_ref_ind2 - L) / (n_ref_ind2 - S) &
       + BesselF1 * (n_ref_ind2 - R) / (n_ref_ind2 - S)
 
    F2 = (n_ref_ind2 * (sin(theta))**2.0 - P) / (2.0 * n_ref_ind2)
 
    F3 = k_vec / abs(k_vec) * (1.0/tan(alpha_ne)) &
         * sin(theta) * cos(theta) * BesselF2
 
    F4 = (R - L) / (2.0 * (n_ref_ind2 - S))
 
    F5 = (P - n_ref_ind2 * (sin(theta))**2.0) / n_ref_ind2
 
    F6 = P * cos(theta) / n_ref_ind2
 
    Fsq = abs(((F1*F2 + F3)**2.0) / ((F4*F5)**2.0 + F6**2.0))
 
end subroutine cold_plasma_polarisation
 
! ==============================================================================
! SUBROUTINE: get_whamp_dispersion
! ==============================================================================
!
! PURPOSE:
!   Computes the hot plasma dispersion surface omega(k, theta) for a given
!   local plasma condition (fce, fpe) by calling the WHAMP engine library.
!   Called once per (L, lambda) point in Integration_lambda.
!
!   For each theta point, WHAMP is called with:
!     k_perp = k * sin(theta)  (WHAMP P parameter)
!     k_par  = k * cos(theta)  (WHAMP Z parameter)
!   where k is scanned over WHAMP_N_K log-spaced points.
!
!   Results are stored in a WHAMPCurve instance and passed down the FDC
!   call chain to WeightFunction and int_N_omega.
!
! INPUTS:
!   fce_khz   - Electron gyrofrequency [kHz]
!   fpe_khz   - Electron plasma frequency [kHz]
!   theta_min - Minimum wave normal angle [radians]
!   theta_max - Maximum wave normal angle [radians]
!
! OUTPUTS:
!   curve     - WHAMPCurve instance with dispersion surface omega(k, theta)
!
! MODULE DEPENDENCIES:
!   constants   - WHAMP_N_K, WHAMP_N_THETA, hot plasma parameters
!   comin       - WHAMP input module (DN, TA, AA, ASS, VD, XC, PM, ZM, XOI, PZL)
!   comoutput   - WHAMP output module (fOUT, VGZOUT, kparOUT, kperpOUT,
!                                      flagSolutionFoundOUT)
!
! CALLS:
!   whamp_engine  - WHAMP dispersion relation solver (from libwhamp.a)
!
! NOTES:
!   - Hot electron population is fixed: 10% abundance, T_par=15 keV,
!     T_perp/T_par=2, proton ions. These are hardcoded for now and should
!     be moved to module constants or input file in a future version.
!   - WHAMP k normalisation: k_WHAMP = k * c / Omega_e
!     Convert to SI: k_SI = k_WHAMP * Omega_e / c = k_WHAMP * 2*pi*fce/c
!   - Starting frequency XOI set to 0.05 * fce_khz (middle of hiss band)
!   - Log-spaced k grid from 10^WHAMP_K_MIN to 10^WHAMP_K_MAX
!   - For parallel propagation (theta=0), k_perp=0 exactly
!
! ==============================================================================
subroutine get_whamp_dispersion(fce_khz, fpe_khz, theta_min, theta_max, curve)
use comin
use comoutput
use comcout
implicit none
 
    real*8, intent(in)               :: fce_khz, fpe_khz
    real*8, intent(in)               :: theta_min, theta_max
    type(WHAMPCurve), intent(out)    :: curve
 
    ! Local variables
    integer :: i_theta, i_k
    real*8  :: theta, k_whamp, k_si, omega_e_rads
    real*8  :: k_min_exp, k_max_exp, k_step_exp
    real*8  :: ne_total, ne_cold, ne_hot
 
    ! Hot plasma parameters (fixed for now)
    real*8, parameter :: hot_fraction  = 0.0d0    ! cold-only validation test
    real*8, parameter :: T_cold_keV    = 0.001d0  ! cold species temperature
    real*8, parameter :: T_hot_par_keV = 15.0d0   ! hot electron T_par
    real*8, parameter :: T_hot_anis    = 2.0d0    ! T_perp/T_par for hot electrons
    real*8, parameter :: k_min_exp_val = -3.0d0   ! k grid: 10^-3 to 10^0
    real*8, parameter :: k_max_exp_val =  1.5d0
    integer, parameter :: n_k          = WHAMP_N_K
    integer, parameter :: n_theta      = WHAMP_N_THETA

    real*8 :: CV_local, TR_local
 
    ! Allocate curve arrays
    curve%n_k     = n_k
    curve%n_theta = n_theta
    curve%fce_khz = fce_khz
    curve%fpe_khz = fpe_khz
 
    allocate(curve%theta(n_theta))
    allocate(curve%k(n_k))
    allocate(curve%k_whamp(n_k))
    allocate(curve%omega(n_theta, n_k))
    allocate(curve%dw_dk(n_theta, n_k))
    allocate(curve%valid(n_theta, n_k))
    allocate(curve%EFL(3, n_theta, n_k))
    
    curve%EFL = (0.0d0, 0.0d0)
    curve%omega  = 0.0d0
    curve%dw_dk  = 0.0d0
    curve%valid  = .false.
 
    ! Omega_e in rad/s (for unit conversions)
    omega_e_rads = 2.0d0 * PI * fce_khz * 1.0d3
 
    ! Compute plasma densities from fce and fpe
    ! ne = (fpe/fce)^2 * fce^2 * 4*pi^2 * eps0 * me / qe^2
    ! Note: fpe^2 = ne * qe^2 / (eps0 * me) => ne = (2*pi*fpe)^2 * eps0 * me / qe^2
    ne_total = (2.0d0 * PI * fpe_khz * 1.0d3)**2 * C0 * me / (qe * qe)
    ne_cold  = (1.0d0 - hot_fraction) * ne_total
    ne_hot   = hot_fraction * ne_total

    TR_local = T_cold_keV   ! since REN(1)=1, TR=TA(1)/1=TA(1)=T_cold_keV
    CV_local = TR_local*(1022.0d0 + TR_local)/(511.0d0 + TR_local)**2.0d0
    CV_local = 1.0d0 / sqrt(CV_local)
 
    ! Set up WHAMP k grid (log-spaced in WHAMP units)
    ! k_whamp stores log10(k) — actual k = 10^k_whamp
    k_step_exp = (k_max_exp_val - k_min_exp_val) / real(n_k - 1, 8)
    do i_k = 1, n_k
        curve%k_whamp(i_k) = k_min_exp_val + (i_k - 1) * k_step_exp
        curve%k(i_k) = 10.0d0**curve%k_whamp(i_k) * omega_e_rads * CV_local / c 
    end do
 
    ! Set up theta grid
    do i_theta = 1, n_theta
        if (n_theta == 1) then
            curve%theta(i_theta) = theta_min
        else
            curve%theta(i_theta) = theta_min + &
                (theta_max - theta_min) / real(n_theta - 1, 8) * (i_theta - 1)
        end if
    end do
 
    ! -------------------------------------------------------------------------
    ! Set up WHAMP comin module variables (plasma parameters)
    ! Species: cold electrons (90%), hot anisotropic electrons (10%), protons
    ! -------------------------------------------------------------------------

    DN(1)   = ne_cold           ! cold electrons [m^-3]
    DN(2)   = ne_total          ! cold protons   [m^-3] (charge neutralizing)
    DN(3)   = ne_hot            ! hot electrons  [m^-3] (0 => cold-only test)


    TA(1)   = T_cold_keV        ! cold electrons T_par [keV]
    TA(2)   = T_cold_keV        ! cold protons   T_par [keV]
    TA(3)   = T_hot_par_keV     ! hot electrons  T_par [keV]

    AA(1,1) = 1.0d0             ! cold electrons T_perp/T_par = 1 (isotropic)
    AA(2,1) = 1.0d0             ! cold protons   T_perp/T_par = 1 (isotropic)
    AA(3,1) = T_hot_anis        ! hot electrons  T_perp/T_par anisotropy

    ASS(1)  = 0.0d0             ! electrons (mass=0 means electron mass)
    ASS(2)  = 1.0d0             ! protons (mass=1 proton mass)
    ASS(3)  = 0.0d0             ! electrons (hot)
 
    XC      = fce_khz           ! electron gyrofrequency [kHz]
    PZL     = 1.0d0             ! log scale for k grid
    XOI     = 0.05d0 / fce_khz      ! starting frequency: 0.5 kHz in omega/Omega_e units
 
    ! -------------------------------------------------------------------------
    ! Call WHAMP for each theta point
    ! -------------------------------------------------------------------------
    do i_theta = 1, n_theta
        theta = curve%theta(i_theta)
        XOI = 0.05d0 / fce_khz   ! reset starting guess at the start of each theta

        do i_k = 1, n_k
            k_whamp = curve%k_whamp(i_k)

            PZL = 1.0d0
            if (abs(theta) < 1.0d-6) then
                PM(1) = -10.0d0
                PM(2) = -10.0d0
                PM(3) = 1.0d0
                ZM(1) = k_whamp
                ZM(2) = k_whamp
                ZM(3) = 1.0d0
            else
                PM(1) = k_whamp + log10(sin(theta))
                PM(2) = PM(1)
                PM(3) = 1.0d0
                ZM(1) = k_whamp + log10(cos(theta))
                ZM(2) = ZM(1)
                ZM(3) = 1.0d0
            end if

            call whamp_engine

            if (flagSolutionFoundOUT(1, 1) == 1) then
                block
                    real*8 :: omega_candidate
                    omega_candidate = real(fOUT(1, 1)) * omega_e_rads
                    if (omega_candidate > 0.0d0 .and. omega_candidate < omega_e_rads) then
                        curve%valid(i_theta, i_k) = .true.
                        curve%omega(i_theta, i_k) = omega_candidate
                        curve%dw_dk(i_theta, i_k) = VGZOUT(1, 1) * c / CV
                        curve%EFL(1, i_theta, i_k) = EFL(1)
                        curve%EFL(2, i_theta, i_k) = EFL(2)
                        curve%EFL(3, i_theta, i_k) = EFL(3)
                        XOI = real(fOUT(1, 1))
                    end if
                end block
            end if

        end do
    end do
 
end subroutine get_whamp_dispersion

subroutine hot_plasma_polarisation(alpha_ne, nres, &
                                   BesselF1, BesselF2, BesselF3, &
                                   EFL, Fsq)
use constants
implicit none
    real*8, intent(in)          :: alpha_ne
    integer, intent(in)         :: nres
    real*8, intent(in)          :: BesselF1, BesselF2, BesselF3
    complex(kind=8), intent(in) :: EFL(3)
    real*8, intent(out)         :: Fsq

    complex(kind=8) :: E_R, E_L, numerator
    real*8          :: E_norm_sq

    ! Right- and left-hand circular components
    ! EFL(1)=Ex, EFL(2)=Ey, EFL(3)=Ez in field-aligned coords (B along z)
    E_R = 0.5d0 * (EFL(1) + (0d0,1d0)*EFL(2))
    E_L = 0.5d0 * (EFL(1) - (0d0,1d0)*EFL(2))

    ! Resonance interaction: Jn-1*E_R + Jn+1*E_L + Jn*Ez*cos(alpha)/sin(alpha)
    numerator = BesselF1 * E_R + BesselF3 * E_L &
                + BesselF2 * EFL(3) * cos(alpha_ne) / sin(alpha_ne)

    ! Normalise by wave energy |E|^2
    E_norm_sq = real(EFL(1)*conjg(EFL(1)) + EFL(2)*conjg(EFL(2)) + &
                     EFL(3)*conjg(EFL(3)))

    if (E_norm_sq > 0.0d0) then
        Fsq = real(numerator * conjg(numerator)) / E_norm_sq
    else
        Fsq = 0.0d0
    end if

end subroutine hot_plasma_polarisation

end module dispersion
