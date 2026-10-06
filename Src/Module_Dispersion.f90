! ==============================================================================
! MODULE: dispersion
! ==============================================================================
!
! PURPOSE:
!   Provides the dispersion physics used throughout the code: the cold plasma
!   (Stix) relation, the kinetic relation via WHAMP, the EMIC band logic, and
!   the derived quantities built on them (group velocity, polarisation,
!   Jacobians). The model in use is selected by DispersionType in module
!   constants ('COLD', 'HOT' or 'EMIC').
!
!   STATUS: COLD and HOT are in production use and validated. EMIC is work
!   in progress: unvalidated and incomplete, so EMIC results should not be
!   used without checking the implementation first.
!
! CONTENTS:
!   Derived types
!     DispersionProperties   - dispersion outputs for one evaluation point
!     WHAMPCurve             - tabulated omega(theta, k) surface for one
!                              local plasma state
!     WHAMPCurveCache        - set of curves along a field line, keyed by
!                              (L, alpha); WHAMP_CACHE_SIZE entries, LRU
!
!   Dispatch
!     Get_dispersion_properties - cold or hot path for one (omega, k, theta)
!     Get_polarisation_factor   - Fsq from a DispersionProperties instance
!
!   Cold plasma
!     cold_plasma_parameters       - Stix R, L, P, S, D
!     cold_plasma_parameters_multi - the same for several ion species
!     cold_plasma_refractive_index - n^2 from the Stix parameters
!     cold_plasma_stix_derivatives - d(R,L,P,S)/d omega
!     cold_plasma_group_velocity   - dw/dk from the analytic derivatives
!     cold_plasma_polarisation     - Fsq from Stix parameters and Bessel terms
!     cold_plasma_efield_diag      - E-field eigenvector (-DDIAG only)
!
!   Hot plasma (WHAMP)
!     get_whamp_dispersion   - builds a WHAMPCurve: WHAMP_N_THETA x WHAMP_N_K
!                              kinetic solutions for one (fce, fpe)
!     get_whamp_point        - single targeted solution, used as a rescue
!                              when interpolation on the curve fails
!     hot_plasma_polarisation- Fsq from WHAMP's E and B field components
!
!   EMIC  (WORK IN PROGRESS - unvalidated, incomplete)
!     emic_band_bounds, emic_find_cutoff, emic_find_seed_omega
!
!   Cunningham formulation
!     cunningham_jacobian    - Jacobian of the Cunningham (2023) kernel
!
!   Profiling accumulators (-DPROFILE): t_*_total, n_*_calls, summed per rank
!   and printed by mpi_parallel.
!
! MODULE DEPENDENCIES:
!   constants  - physical constants, DispersionType, WHAMP_N_K,
!                WHAMP_N_THETA, N_max_hot, k_min/max_exp_val_input,
!                hot_species_* configuration
!
! NOTES:
!   - The original intent was that callers obtain every dispersion quantity
!     through Get_dispersion_properties and Get_polarisation_factor, with no
!     knowledge of the active model. That holds for the cold path; for
!     HOT/EMIC, WeightFunction, int_N_omega and D_local_n interpolate the
!     WHAMPCurve inline (bracketing in theta and k, tolerance fallback,
!     targeted rescue), in three near-identical copies. Moving that into a
!     routine here would restore the intent and remove several hundred lines.
!   - WHAMP curves are expensive: one call to get_whamp_dispersion solves
!     WHAMP_N_THETA x WHAMP_N_K kinetic dispersion problems. They depend only
!     on the local plasma state, i.e. on (L, lambda), not on energy or pitch
!     angle, which is what the WHAMPCurveCache exploits.
!   - The cold Stix parameters are no longer written to module variables as
!     a side effect; they are returned through the argument list.
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
!     EFL, BFL       - wave E and B field components from WHAMP, used by
!                      hot_plasma_polarisation
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
 
    ! Hot plasma (WHAMP) outputs
    complex(kind=8) :: EFL(3)   ! E field components from WHAMP
    complex(kind=8) :: BFL(3)   ! B field components from WHAMP
 
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
!   - fce_khz and fpe_khz identify the plasma state the curve was built for;
!     they are also the cache key used by WHAMPCurveCache
!   - EFL and BFL carry the wave field components per (theta, k), used by
!     hot_plasma_polarisation
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
    complex(kind=8), allocatable :: BFL(:,:,:) ! B field components (3, n_theta, n_k)
end type WHAMPCurve

integer, parameter :: WHAMP_CACHE_SIZE = 4   ! number of (L,alpha) curve-sets to hold

type :: WHAMPCurveCache
    logical :: valid = .false.
    real*8  :: L_cached = -1.0d0
    real*8  :: alpha_cached = -1.0d0
    integer :: n_lambda_cached = 0
    real*8, allocatable :: lambda_cached(:)
    type(WHAMPCurve), allocatable :: curves(:)
end type WHAMPCurveCache

type(WHAMPCurveCache), save :: whamp_cache(WHAMP_CACHE_SIZE)
integer, save :: whamp_cache_lru_counter = 0
integer, save :: whamp_cache_last_used(WHAMP_CACHE_SIZE) = 0

! Profiling accumulators
real*8,  save :: t_roots_total  = 0.0d0
real*8,  save :: t_weight_total = 0.0d0
real*8,  save :: t_intn_total   = 0.0d0
real*8,  save :: t_intlam_total = 0.0d0
integer, save :: n_roots_calls  = 0
integer, save :: n_weight_calls = 0
integer, save :: n_intn_calls   = 0
integer, save :: n_intlam_calls = 0

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
            call cold_plasma_refractive_index(omega, theta, &
                                              DP%P, DP%S, DP%R, DP%L, &
                                              DP%n_ref_ind2)
            call cold_plasma_group_velocity(omega, k_vec, theta, &
                                            DP%R, DP%L, DP%P, DP%S, &
                                            omega_e, omega_pe, omega_i, omega_pi, &
                                            DP%dw_dk)
 
        case('HOT', 'EMIC')
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
        case('HOT', 'EMIC')
            call hot_plasma_polarisation(alpha_ne, nres, k_vec, theta, &
                                         BesselF1, BesselF2, BesselF3, &
                                         DP%EFL, DP%BFL, DP%n_ref_ind2, Fsq)
 
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
subroutine cold_plasma_refractive_index(omega, theta, P, S, R, L, &
                                        n_ref_ind2)
implicit none
    real*8, intent(in)  :: omega, theta, P, S, R, L
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
! SUBROUTINE: cold_plasma_stix_derivatives
! ==============================================================================
!
! PURPOSE:
!   Computes the analytic derivatives of the Stix cold plasma parameters
!   R, L, P, S with respect to omega, at fixed plasma state. Factored out
!   of cold_plasma_group_velocity so the same derivatives can be reused by
!   cunningham_jacobian (Cunningham 2023, Eq. 4) without duplicating the
!   formulas.
!
!   These are the same dPdw, dSdw, dRdw, dLdw expressions previously
!   computed inline in cold_plasma_group_velocity (Albert 2007 derivative
!   formulas). Behaviour and numerical results are UNCHANGED from the
!   original implementation; only the computation has been moved here.
!
! INPUTS:
!   omega       - Wave frequency [rad/s]
!   omega_e     - Electron gyrofrequency [rad/s]
!   omega_pe    - Electron plasma frequency [rad/s]
!   omega_i     - Ion gyrofrequency [rad/s]
!   omega_pi    - Ion plasma frequency [rad/s]
!
! OUTPUTS:
!   dPdw, dSdw, dRdw, dLdw - d/domega of the Stix P, S, R, L parameters [s/rad]
!
! NOTES:
!   - a1-a6 are the same local polynomial coefficients previously computed
!     inline in cold_plasma_group_velocity. a5 is analytically zero; see
!     the original note regarding Albert (2007).
!
! ==============================================================================
subroutine cold_plasma_stix_derivatives(omega, omega_e, omega_pe, omega_i, omega_pi, &
                                        dPdw, dSdw, dRdw, dLdw)
use constants
implicit none
    real*8, intent(in)  :: omega, omega_e, omega_pe, omega_i, omega_pi
    real*8, intent(out) :: dPdw, dSdw, dRdw, dLdw
 
    real*8 :: a1, a2, a3, a4, a5, a6
 
    ! Plasma frequency combination coefficients
    ! (definitions match those previously set in Roots_finder /
    !  cold_plasma_group_velocity)
    a1 = omega_i**2.0 + omega_e**2.0 + omega_pe**2.0 + omega_pi**2.0
    a2 = (omega_i*omega_e)**2.0 + (omega_pe*omega_i)**2.0 + (omega_pi*omega_e)**2.0
    a3 = omega_i + real(sqe)*omega_e
    a4 = omega_i*real(sqe)*omega_e - omega_pe**2.0 - omega_pi**2.0
    a5 = 0.0   ! analytically zero; see note in original cold_plasma_group_velocity
    a6 = omega_pi**2.0 + omega_pe**2.0
 
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
 
end subroutine cold_plasma_stix_derivatives

 
! ==============================================================================
! SUBROUTINE: cold_plasma_group_velocity  [REFACTORED]
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
!   Albert (2007), equations for the derivative of the dispersion tensor,
!   and are now obtained from the shared helper cold_plasma_stix_derivatives
!   rather than being computed inline. This is a pure refactoring: the
!   numerical result of this subroutine is UNCHANGED from the original
!   implementation.
!
! INPUTS:
!   omega       - Wave frequency [rad/s]
!   k_vec       - Wave vector magnitude [rad/m]
!   theta       - Wave normal angle [radians]
!   R, L, P, S  - Stix parameters [-]
!   omega_e     - Electron gyrofrequency [rad/s]
!   omega_pe    - Electron plasma frequency [rad/s]
!   omega_i     - Ion gyrofrequency [rad/s]
!   omega_pi    - Ion plasma frequency [rad/s]
!
! OUTPUTS:
!   dw_dk       - Group velocity along wave vector [m/s]
!
! NOTES:
!   - a5 (inside cold_plasma_stix_derivatives) is set to 0.0, consistent
!     with the original implementation. See note there regarding whether
!     this is a valid simplification (inherited, unchanged).
!
! ==============================================================================
subroutine cold_plasma_group_velocity(omega, k_vec, theta, R, L, P, S, &
                                      omega_e, omega_pe, omega_i, omega_pi, &
                                      dw_dk)
use constants
implicit none
    real*8, intent(in)  :: omega, k_vec, theta, R, L, P, S
    real*8, intent(in)  :: omega_e, omega_pe, omega_i, omega_pi
    real*8, intent(out) :: dw_dk
 
    real*8 :: X
 
    ! Analytic Stix parameter derivatives w.r.t. omega (now from shared helper)
    real*8 :: dPdw, dSdw, dRdw, dLdw
 
    ! Dispersion function derivatives
    real*8 :: dF1dw, dF2dw, dF3dw, dDdw, dDdk
 
    X = tan(theta)
 
    call cold_plasma_stix_derivatives(omega, omega_e, omega_pe, omega_i, omega_pi, &
                                      dPdw, dSdw, dRdw, dLdw)
 
    ! dD/dk: derivative of dispersion function w.r.t. k
    dDdk = (S*(X**2.0) + P) * (c**4.0) / (omega**4.0) * 4.0*(k_vec**3.0) &
           - (R*L*(X**2.0) + (2.0 + X**2.0)*P*S) * (c**2.0)/(omega**2.0) &
             * 2.0*k_vec
 
    ! dD/dw: derivative of dispersion function w.r.t. omega
    dF1dw = ((k_vec*c)**4.0) / (omega**5.0) * &
            (omega*(X**2.0)*dSdw + omega*dPdw - 4.0*(S*(X**2.0) + P)) !((k_vec*c)**4.0) / (omega**5.0) * (omega*dA_-4*A_)
 
    dF2dw = ((k_vec*c)**2.0) / (omega**3.0) * &
            (omega*((X**2.0)*(L*dRdw + R*dLdw) + &
                    (2.0 + X**2.0)*(P*dSdw + S*dPdw)) - &
             2.0*(R*L*(X**2.0) + (2.0 + X**2.0)*P*S)) !((k_vec*c)**2.0) / (omega**3.0) *(omega*dB_-2*B_)
 
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
! SUBROUTINE: cold_plasma_efield_diag   (DIAGNOSTIC ONLY)
! ==============================================================================
!   Returns the cold-plasma E-field eigenvector in WHAMP's basis
!   (e1 along k_perp, e2 = e3 x e1, e3 along B0) and the transverse
!   fraction sin^2(theta_kE) = 1 - |E.khat|^2/|E|^2.
!
!   At T->0 this MUST agree with the same quantity derived from WHAMP's
!   EFL/BFL. Used to test whether the WHAMP eigenvector pipeline is
!   returning physically plausible polarisation at modest theta.
!
!   Null vector obtained as the cross product of two rows of the wave
!   matrix M (valid because M is singular at the dispersion root).
!   Overall complex phase is arbitrary; sin^2(theta_kE) is phase-invariant.
! ==============================================================================
subroutine cold_plasma_efield_diag(theta, n2_cold, R, L, P, S, D, &
                                   E_cold, sin2_theta_kE_cold, ok)
implicit none
    real*8, intent(in)  :: theta, n2_cold, R, L, P, S, D
    complex(kind=8), intent(out) :: E_cold(3)
    real*8, intent(out) :: sin2_theta_kE_cold
    logical, intent(out) :: ok

    complex(kind=8) :: r1(3), r2(3), r3(3), v(3), vb(3)
    complex(kind=8) :: E_par_k
    real*8 :: st, ct, nrm, nrm_best, E_norm_sq
    complex(kind=8), parameter :: ii = (0.0d0, 1.0d0)

    ok = .false.
    E_cold = (0.0d0, 0.0d0)
    sin2_theta_kE_cold = -1.0d0
    if (n2_cold <= 0.0d0) return

    st = sin(theta)
    ct = cos(theta)

    ! Wave matrix rows: (n^2 dij - ni nj - eps_ij) E = 0,  k = n(sin,0,cos)
    r1(1) = cmplx(S - n2_cold*ct*ct, 0.0d0, kind=8)
    r1(2) = -ii*D
    r1(3) = cmplx(n2_cold*st*ct, 0.0d0, kind=8)

    r2(1) = ii*D
    r2(2) = cmplx(S - n2_cold, 0.0d0, kind=8)
    r2(3) = (0.0d0, 0.0d0)

    r3(1) = cmplx(n2_cold*st*ct, 0.0d0, kind=8)
    r3(2) = (0.0d0, 0.0d0)
    r3(3) = cmplx(P - n2_cold*st*st, 0.0d0, kind=8)

    ! Pick the row pair giving the best-conditioned null vector
    nrm_best = 0.0d0
    call cross3(r2, r3, v)
    nrm = cnorm(v)
    if (nrm > nrm_best) then
        nrm_best = nrm
        vb = v
    end if
    call cross3(r1, r3, v)
    nrm = cnorm(v)
    if (nrm > nrm_best) then
        nrm_best = nrm
        vb = v
    end if
    call cross3(r1, r2, v)
    nrm = cnorm(v)
    if (nrm > nrm_best) then
        nrm_best = nrm
        vb = v
    end if

    if (nrm_best <= 0.0d0) return

    E_cold = vb / nrm_best
    E_norm_sq = real(E_cold(1)*conjg(E_cold(1)) + E_cold(2)*conjg(E_cold(2)) &
                   + E_cold(3)*conjg(E_cold(3)))
    if (E_norm_sq <= 0.0d0) return

    E_par_k = E_cold(1)*st + E_cold(3)*ct
    sin2_theta_kE_cold = 1.0d0 - real(E_par_k*conjg(E_par_k)) / E_norm_sq
    sin2_theta_kE_cold = max(0.0d0, min(1.0d0, sin2_theta_kE_cold))
    ok = .true.

contains

    subroutine cross3(a, b, cr)
        complex(kind=8), intent(in)  :: a(3), b(3)
        complex(kind=8), intent(out) :: cr(3)
        cr(1) = a(2)*b(3) - a(3)*b(2)
        cr(2) = a(3)*b(1) - a(1)*b(3)
        cr(3) = a(1)*b(2) - a(2)*b(1)
    end subroutine cross3

    real*8 function cnorm(a)
        complex(kind=8), intent(in) :: a(3)
        cnorm = sqrt(real(a(1)*conjg(a(1)) + a(2)*conjg(a(2)) + a(3)*conjg(a(3))))
    end function cnorm

end subroutine cold_plasma_efield_diag
 
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
!   - Plasma species composition (mass, density fraction, T_par, anisotropy)
!     is fully configurable via module constants (n_hot_species,
!     hot_species_mass, hot_species_density_frac, hot_species_T_par_keV,
!     hot_species_T_anis), populated from the input file in FDC_MPI.f90.
!     No default composition is assumed; DispersionType='HOT' requires
!     this block to be present in the input file.
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


    integer :: my_rank_local, mpi_ierr_local
 
    ! Local variables
    integer :: i_theta, i_k, i_sp
    real*8  :: theta, k_whamp, k_si, omega_e_rads
    real*8  :: k_step_exp
    real*8  :: ne_total
 
    integer, parameter :: n_k          = WHAMP_N_K
    integer             :: n_theta

    real*8 :: CV_local, TR_local
    real*8  :: omega_target_ion, omega_next_below_ion, omega_cutoff_ion
    logical :: cutoff_success
 
    ! Parallel-only (delta-function angle distribution) callers pass
    ! theta_min=theta_max=0; skip the redundant oblique theta sweep.
    if ((theta_max - theta_min) < 1.0d-10) then
        n_theta = 1
    else
        n_theta = WHAMP_N_THETA
    end if

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
    allocate(curve%BFL(3, n_theta, n_k))
    
    curve%BFL = (0.0d0, 0.0d0)
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

    TR_local = hot_species_T_par_keV(1)   ! WHAMP normalizes CV by species 1's thermal velocity
    CV_local = TR_local*(1022.0d0 + TR_local)/(511.0d0 + TR_local)**2.0d0
    CV_local = 1.0d0 / sqrt(CV_local)

    ! For EMIC, precompute the cold N-ion L-mode pass-band bracket ONCE per curve
    ! (composition/fce/fpe are fixed here -- independent of theta and k). This
    ! replaces the old single-ion whistler-shaped XOI formula, which had no
    ! lower cutoff and so had no way to distinguish the H+ band from the
    ! neighboring He+/O+ bands.
    cutoff_success   = .false.
    omega_target_ion = 0.0d0
    omega_cutoff_ion = 0.0d0
    if (trim(DispersionType) == 'EMIC') then
        omega_target_ion = omega_e_rads * me / (hot_species_mass(1) * mi)
        call emic_band_bounds(omega_e_rads, omega_target_ion, omega_next_below_ion)
        call emic_find_cutoff(omega_e_rads, ne_total, omega_next_below_ion, &
                            omega_target_ion, omega_cutoff_ion, cutoff_success)
    end if
 
    ! Set up WHAMP k grid (log-spaced in WHAMP units)
    ! k_whamp stores log10(k) — actual k = 10^k_whamp
    k_step_exp = (k_max_exp_val_input - k_min_exp_val_input) / real(n_k - 1, 8)
    do i_k = 1, n_k
        curve%k_whamp(i_k) = k_min_exp_val_input + (i_k - 1) * k_step_exp
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

    ! Set up WHAMP comin module variables (plasma parameters)
    do i_sp = 1, n_hot_species
        DN(i_sp)   = hot_species_density_frac(i_sp) * ne_total
        TA(i_sp)   = hot_species_T_par_keV(i_sp)
        AA(i_sp,1) = hot_species_T_anis(i_sp)
        AA(i_sp,2) = hot_species_T_anis(i_sp)
        DD(i_sp)   = 1.0d0
        ASS(i_sp)  = hot_species_mass(i_sp)
    end do
 
    XC      = fce_khz           ! electron gyrofrequency [kHz]
    PZL     = 1.0d0             ! log scale for k grid
    !XOI = 0.05d0 / fce_khz      ! starting frequency: 0.5 kHz in omega/Omega_e units
 
! -------------------------------------------------------------------------
    ! Call WHAMP for each theta point
    ! -------------------------------------------------------------------------
    do i_theta = 1, n_theta
        theta = curve%theta(i_theta)

        if (trim(DispersionType) == 'EMIC' .and. cutoff_success) then

            block
                logical :: have_prev, seed_ok, valid_tmp
                real*8  :: prev_omega, cold_seed_omega
                real*8  :: n_ref_ind2_tmp, dw_dk_tmp, omega_tmp
                complex(kind=8) :: EFL_tmp(3), BFL_tmp(3)

                logical :: walk_ok
                real*8  :: walk_omega, dw_dk_walk
                complex(kind=8) :: EFL_walk(3), BFL_walk(3)
                real*8  :: k_whamp_start, k_whamp_end, k_whamp_sub, k_sub_si
                integer :: j_sub
                integer, parameter :: n_sub_steps = 20

                have_prev  = .false.
                prev_omega = 0.0d0

                do i_k = 1, n_k

                    if (i_k == 1) then
                        ! First grid point: cold seed only, single attempt.
                        valid_tmp = .false.
                        call emic_find_seed_omega(curve%k(i_k), omega_e_rads, ne_total, &
                                                   omega_cutoff_ion, omega_target_ion, &
                                                   cold_seed_omega, seed_ok)
                        if (seed_ok) then
                            call get_whamp_point(fce_khz, fpe_khz, theta, curve%k(i_k), &
                                                  cold_seed_omega, n_ref_ind2_tmp, dw_dk_tmp, &
                                                  EFL_tmp, BFL_tmp, valid_tmp)
                            if (valid_tmp) omega_tmp = curve%k(i_k) * c / sqrt(n_ref_ind2_tmp)
                        end if

                        if (valid_tmp .and. omega_tmp > omega_cutoff_ion * 0.98d0 .and. &
                            omega_tmp < omega_target_ion * 1.02d0) then
                            curve%valid(i_theta, i_k) = .true.
                            curve%omega(i_theta, i_k) = omega_tmp
                            curve%dw_dk(i_theta, i_k) = dw_dk_tmp
                            curve%EFL(:, i_theta, i_k) = EFL_tmp
                            curve%BFL(:, i_theta, i_k) = BFL_tmp
                            have_prev  = .true.
                            prev_omega = omega_tmp
                        end if

                        cycle
                    end if

                    ! ------------------------------------------------------
                    ! i_k >= 2: sub-step from grid point i_k-1 to grid point
                    ! i_k in n_sub_steps small increments of log10(k).
                    ! ------------------------------------------------------
                    k_whamp_start = curve%k_whamp(i_k - 1)
                    k_whamp_end   = curve%k_whamp(i_k)

                    walk_ok    = have_prev
                    walk_omega = prev_omega

                    do j_sub = 1, n_sub_steps
                        k_whamp_sub = k_whamp_start + &
                            (k_whamp_end - k_whamp_start) * real(j_sub, 8) / real(n_sub_steps, 8)
                        k_sub_si = 10.0d0**k_whamp_sub * omega_e_rads * CV_local / c

                        valid_tmp = .false.

                        if (walk_ok) then
                            call get_whamp_point(fce_khz, fpe_khz, theta, k_sub_si, &
                                                  walk_omega, n_ref_ind2_tmp, dw_dk_tmp, &
                                                  EFL_tmp, BFL_tmp, valid_tmp)
                            if (valid_tmp) then
                                omega_tmp = k_sub_si * c / sqrt(n_ref_ind2_tmp)
                                if (.not. (omega_tmp > omega_cutoff_ion * 0.98d0 .and. &
                                           omega_tmp < omega_target_ion * 1.02d0 .and. &
                                           omega_tmp >= walk_omega * 0.999d0)) then
                                    valid_tmp = .false.   ! converged, but wrong branch -- retry cold below
                                end if
                            end if
                        end if

                        if (.not. valid_tmp) then
                            call emic_find_seed_omega(k_sub_si, omega_e_rads, ne_total, &
                                                       omega_cutoff_ion, omega_target_ion, &
                                                       cold_seed_omega, seed_ok)
                            if (seed_ok) then
                                call get_whamp_point(fce_khz, fpe_khz, theta, k_sub_si, &
                                                      cold_seed_omega, n_ref_ind2_tmp, dw_dk_tmp, &
                                                      EFL_tmp, BFL_tmp, valid_tmp)
                                if (valid_tmp) then
                                    omega_tmp = k_sub_si * c / sqrt(n_ref_ind2_tmp)
                                    if (.not. (omega_tmp > omega_cutoff_ion * 0.98d0 .and. &
                                               omega_tmp < omega_target_ion * 1.02d0 .and. &
                                               (.not. walk_ok .or. omega_tmp >= walk_omega * 0.999d0))) then
                                        valid_tmp = .false.
                                    end if
                                end if
                            end if
                        end if

                        if (valid_tmp) then
                            walk_omega = omega_tmp
                            dw_dk_walk = dw_dk_tmp
                            EFL_walk   = EFL_tmp
                            BFL_walk   = BFL_tmp
                            walk_ok    = .true.
                        else
                            walk_ok = .false.
                            ! subsequent sub-steps in this interval will each
                            ! retry fresh from a cold seed (walk_ok=.false.
                            ! skips the continuation attempt above), so the
                            ! walk can still recover before reaching i_k
                        end if
                    end do   ! j_sub

                    if (walk_ok) then
                        curve%valid(i_theta, i_k) = .true.
                        curve%omega(i_theta, i_k) = walk_omega
                        curve%dw_dk(i_theta, i_k) = dw_dk_walk
                        curve%EFL(:, i_theta, i_k) = EFL_walk
                        curve%BFL(:, i_theta, i_k) = BFL_walk
                        have_prev  = .true.
                        prev_omega = walk_omega
                    else
                        have_prev = .false.   ! next interval starts fresh, not from a stale/far seed
                    end if

                end do   ! i_k
            end block

        else
            ! ---- existing per-point independent-seed logic (HOT
            !      chorus/hiss, or EMIC without a valid cutoff bracket) --
            !      UNCHANGED from current file ----
            do i_k = 1, n_k
                k_whamp = curve%k_whamp(i_k)

                block
                    real*8 :: k_si_guess, kc_over_wp_sq
                    real*8 :: omega_ref_rads, wp_ref_rads

                    k_si_guess = 10.0d0**k_whamp * omega_e_rads * CV_local / c

                    if (trim(DispersionType) == 'EMIC') then
                        if (cutoff_success) then
                            block
                                real*8  :: omega_seed_local
                                logical :: seed_success
                                call emic_find_seed_omega(k_si_guess, omega_e_rads, ne_total, &
                                                        omega_cutoff_ion, omega_target_ion, &
                                                        omega_seed_local, seed_success)
                                if (seed_success) then
                                    XOI = omega_seed_local / omega_e_rads
                                else
                                    XOI = omega_cutoff_ion / omega_e_rads
                                end if
                            end block
                        else
                            omega_ref_rads = omega_e_rads * me / (hot_species_mass(1) * mi)
                            wp_ref_rads    = 2.0d0*PI*fpe_khz*1.0d3 * &
                                            sqrt(hot_species_density_frac(1) * me / &
                                                (hot_species_mass(1) * mi))
                            kc_over_wp_sq  = (k_si_guess * c / wp_ref_rads)**2.0d0
                            XOI = cos(theta) * kc_over_wp_sq / (1.0d0 + kc_over_wp_sq) &
                                * (omega_ref_rads / omega_e_rads)
                        end if
                    else
                        omega_ref_rads = omega_e_rads
                        wp_ref_rads    = 2.0d0*PI*fpe_khz*1.0d3
                        kc_over_wp_sq  = (k_si_guess * c / wp_ref_rads)**2.0d0
                        XOI = cos(theta) * kc_over_wp_sq / (1.0d0 + kc_over_wp_sq)
                    end if
                    if (XOI < 1.0d-6) XOI = 1.0d-6
                end block

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
                        real*8  :: omega_candidate, ratio_check
                        logical :: accept

                        omega_candidate = real(fOUT(1, 1)) * omega_e_rads

                        if (XOI > 1.0d-10) then
                            ratio_check = real(fOUT(1,1)) / XOI
                        else
                            ratio_check = 1.0d0
                        end if

                        if (trim(DispersionType) == 'EMIC' .and. cutoff_success) then
                            accept = (omega_candidate > omega_cutoff_ion * 0.98d0) .and. &
                                    (omega_candidate < omega_target_ion * 1.02d0) .and. &
                                    (ratio_check > 0.05d0) .and. (ratio_check < 20.0d0)
                        else
                            accept = (omega_candidate > 0.0d0) .and. (omega_candidate < omega_e_rads) &
                                    .and. (ratio_check > 0.2d0) .and. (ratio_check < 5.0d0)
                        end if

                        if (accept) then
                            curve%valid(i_theta, i_k) = .true.
                            curve%omega(i_theta, i_k) = omega_candidate
                            curve%dw_dk(i_theta, i_k) = VGZOUT(1, 1) * c / CV
                            curve%EFL(1, i_theta, i_k) = ExOUT(1, 1)
                            curve%EFL(2, i_theta, i_k) = EyOUT(1, 1)
                            curve%EFL(3, i_theta, i_k) = EzOUT(1, 1)
                            curve%BFL(1, i_theta, i_k) = BxOUT(1, 1)
                            curve%BFL(2, i_theta, i_k) = ByOUT(1, 1)
                            curve%BFL(3, i_theta, i_k) = BzOUT(1, 1)
                        end if
                    end block
                end if
            end do
        end if
    end do

    !open(99, file='whamp_debug_check_a0005.txt', status='unknown', position='append')
    !write(99,*) 'fce=', fce_khz, ' fpe=', fpe_khz, &
    !            ' omega_min(Hz)=', minval(curve%omega, mask=curve%valid)/(2.0d0*PI), &
    !            ' omega_max(Hz)=', maxval(curve%omega, mask=curve%valid)/(2.0d0*PI), &
    !            ' valid_count=', count(curve%valid), '/', size(curve%valid)
    !close(99)

    block
        real*8, parameter :: WHAMP_DEBUG_FCE_MIN = 0.0d0
        real*8, parameter :: WHAMP_DEBUG_FCE_MAX = 1.0d20
        real*8, parameter :: WHAMP_DEBUG_THETA_TARGET = 8.9795916527509672d0  ! degrees
        integer :: i_theta_dbg, i_k_dbg, i_theta_closest
        real*8  :: min_theta_diff, theta_diff

        if (fce_khz >= WHAMP_DEBUG_FCE_MIN .and. fce_khz <= WHAMP_DEBUG_FCE_MAX) then
            min_theta_diff = 1.0d20
            i_theta_closest = 1
            do i_theta_dbg = 1, n_theta
                theta_diff = abs(curve%theta(i_theta_dbg)*180.0d0/PI - WHAMP_DEBUG_THETA_TARGET)
                if (theta_diff < min_theta_diff) then
                    min_theta_diff = theta_diff
                    i_theta_closest = i_theta_dbg
                end if
            end do

            !open(97, file='whamp_nres1_target_debug.txt', status='unknown', position='append')
            !write(97,*) '=== fce_khz=', fce_khz, ' fpe_khz=', fpe_khz, &
            !            ' target_theta_found(deg)=', curve%theta(i_theta_closest)*180.0d0/PI, ' ==='
            !do i_k_dbg = 1, n_k
            !    write(97,*) '  i_k=', i_k_dbg, ' k=', curve%k(i_k_dbg), &
            !                ' omega=', curve%omega(i_theta_closest,i_k_dbg), &
            !                ' valid=', curve%valid(i_theta_closest,i_k_dbg)
            !end do
            !close(97)
        end if
    end block

end subroutine get_whamp_dispersion
 
subroutine get_whamp_point(fce_khz, fpe_khz, theta, k_target, omega_guess, &
                            n_ref_ind2_out, dw_dk_out, EFL_out, BFL_out, valid_out)
    use constants
    use comin
    use comoutput
    use comcout
    implicit none
    real*8, intent(in)  :: fce_khz, fpe_khz, theta, k_target, omega_guess
    real*8, intent(out) :: n_ref_ind2_out, dw_dk_out
    complex(kind=8), intent(out) :: EFL_out(3)
    complex(kind=8), intent(out) :: BFL_out(3)
    logical, intent(out) :: valid_out
 
    real*8 :: omega_e_rads, CV_local, TR_local, k_whamp
    real*8 :: ratio_check
 
    valid_out = .false.
    n_ref_ind2_out = 0.0d0
    dw_dk_out = 0.0d0
    EFL_out = (0.0d0, 0.0d0)
    BFL_out = (0.0d0, 0.0d0)
 
    if (k_target <= 0.0d0 .or. omega_guess <= 0.0d0) return
 
    omega_e_rads = 2.0d0 * PI * fce_khz * 1.0d3
    TR_local = hot_species_T_par_keV(1)
    CV_local = TR_local*(1022.0d0 + TR_local)/(511.0d0 + TR_local)**2.0d0
    CV_local = 1.0d0 / sqrt(CV_local)
 
    k_whamp = log10(k_target / (omega_e_rads * CV_local / c))
 
    XC  = fce_khz
    PZL = 1.0d0
    XOI = omega_guess / omega_e_rads   ! excellent seed - we already know the answer closely
 
    if (abs(theta) < 1.0d-6) then
        PM(1) = -10.0d0; PM(2) = -10.0d0; PM(3) = 1.0d0
        ZM(1) = k_whamp; ZM(2) = k_whamp; ZM(3) = 1.0d0
    else
        PM(1) = k_whamp + log10(sin(theta)); PM(2) = PM(1); PM(3) = 1.0d0
        ZM(1) = k_whamp + log10(cos(theta)); ZM(2) = ZM(1); ZM(3) = 1.0d0
    end if
 
    call whamp_engine
 
    !open(80, file='whamp_point_internal.txt', status='unknown', position='append')
    !write(80,*) 'flag=', flagSolutionFoundOUT(1,1), ' fOUT=', real(fOUT(1,1)), ' XOI_in=', XOI
    !close(80)
 
    if (flagSolutionFoundOUT(1, 1) == 1) then
        block
            real*8 :: omega_candidate
            omega_candidate = real(fOUT(1, 1)) * omega_e_rads
            ratio_check = real(fOUT(1,1)) / XOI
            if (omega_candidate > 0.0d0 .and. omega_candidate < omega_e_rads &
                .and. ratio_check > 0.5d0 .and. ratio_check < 2.0d0) then
                valid_out = .true.
                n_ref_ind2_out = (k_target * c / omega_candidate)**2.0d0
                dw_dk_out = VGZOUT(1, 1) * c / CV_local
                EFL_out(1) = ExOUT(1, 1)
                EFL_out(2) = EyOUT(1, 1)
                EFL_out(3) = EzOUT(1, 1)
                BFL_out(1) = BxOUT(1, 1)
                BFL_out(2) = ByOUT(1, 1)
                BFL_out(3) = BzOUT(1, 1)
            end if
        end block
    end if
end subroutine get_whamp_point


subroutine hot_plasma_polarisation(alpha_ne, nres, theta, k_vec, &
                                   BesselF1, BesselF2, BesselF3, &
                                   EFL, BFL, n_ref_ind2, Fsq)
use constants
implicit none
    real*8, intent(in)          :: alpha_ne, theta, k_vec
    integer, intent(in)         :: nres
    real*8, intent(in)          :: BesselF1, BesselF2, BesselF3
    complex(kind=8), intent(in) :: EFL(3)
    complex(kind=8), intent(in) :: BFL(3)
    real*8, intent(in)          :: n_ref_ind2
    real*8, intent(out)         :: Fsq

    complex(kind=8) :: E_R, E_L, Theta_nk, E_par_k
    real*8          :: E_norm_sq, sin2_theta_kE, k_sgn

    real*8, parameter :: sin2_theta_kE_min = 1.0d-12

    Fsq = 0.0d0

    E_norm_sq = real(EFL(1)*conjg(EFL(1)) + EFL(2)*conjg(EFL(2)) &
                   + EFL(3)*conjg(EFL(3)))
    if (E_norm_sq <= 0.0d0) return

    ! transverse fraction of E relative to k, from EFL alone
    E_par_k = EFL(1)*sin(theta) + EFL(3)*cos(theta)
    sin2_theta_kE = 1.0d0 - real(E_par_k*conjg(E_par_k)) / E_norm_sq
    sin2_theta_kE = max(0.0d0, min(1.0d0, sin2_theta_kE))   ! guard roundoff
    if (sin2_theta_kE < sin2_theta_kE_min) return

    ! E_R, E_L per Lyons (1974), text following Eq. (3):
    !    "E_k,L^R = (Ex -+ i Ey)/sqrt(2)", R takes the upper (minus) sign.
    E_R = (EFL(1) - (0.0d0,1.0d0)*EFL(2)) / sqrt(2.0d0)
    E_L = (EFL(1) + (0.0d0,1.0d0)*EFL(2)) / sqrt(2.0d0)

    !  Theta_nk, Lyons (1974) Eq. (3) 
    if (k_vec >= 0.0d0) then
        k_sgn =  1.0d0
    else
        k_sgn = -1.0d0
    end if

    Theta_nk = (E_R * BesselF3 + E_L * BesselF1) / sqrt(2.0d0) &
             + k_sgn * EFL(3) * BesselF2 * cos(alpha_ne) / sin(alpha_ne)

    ! Phi_n^2 = |Theta_nk|^2 / (|E|^2 sin^2(theta_kE))
    Fsq = real(Theta_nk * conjg(Theta_nk)) / (E_norm_sq * sin2_theta_kE)

end subroutine hot_plasma_polarisation

 
! ==============================================================================
! SUBROUTINE: cunningham_jacobian
! ==============================================================================
!
! PURPOSE:
!   Computes the Jacobian used in the denominator of Cunningham (2023) Eq. 4b
!
!   Verbatim formula, following the PIRAN code (with mu = c*k/omega):
!     A  = S*X^2 + P
!     B  = R*L*X^2 + P*S*(2+X^2)
!     C  = P*R*L*(1+X^2)
!     dA = dS*X^2 + dP
!     dB = (dR*L + R*dL)*X^2 + (dP*S + P*dS)*(2+X^2)
!     dC = (dP*R*L + P*dR*L + P*R*dL)*(1+X^2)
!     jacobian = (k^2/(1+X^2)) * [ (dA*mu^4 - dB*mu^2 + dC)
!                                    / (2*(2*A*mu^4 - B*mu^2)) - 1/omega ]
!
!
! IMPORTANT WARNING -- DO NOT WIRE THIS DIRECTLY INTO G2 = g(X)/N:
!   Investigation during development of the 'CUNNINGHAM' DiffusionFormulation
!   option established that this Jacobian is NOT a drop-in replacement for
!   int_N_omega's N(omega) inside WeightFunction's existing
!   G2 = g(X)/N(omega) formula. Two structural facts must both be respected
!   together by any future correct implementation, established by reviewing
!   the PIRAN reference implementation (RB-ENVIRONMENT/PIRAN):
!     (1) Under Cunningham's method, g(X) must be pre-normalised by its own
!         integral over the configured X range (integral_gx in PIRAN),
!         BEFORE being divided by this Jacobian -- i.e. the quantity used
!         is g(X)/integral(g(X)dX), NOT raw g(X) as int_N_omega's path uses.
!     (2) This Jacobian substitutes for a DIFFERENT term, in a DIFFERENTLY
!         -structured diffusion-coefficient formula (Cunningham/Glauert &
!         Horne Eq. 11-13, via PIRAN's get_DnX_single_root), than the one
!         WeightFunction currently implements (Albert 2007's Delta_n/G1/G3/
!         G4 reformulation, which is specific to the LYONS path and is NOT
!         interchangeable with Cunningham's formula at the level of a
!         single substituted term).
!   A correct 'CUNNINGHAM' implementation requires re-deriving
!   WeightFunction's diffusion-coefficient kernel from Cunningham/Glauert &
!   Horne Eq. 11-13 directly (see PIRAN's get_phi_squared, get_singular_term,
!   get_DnX_single_root for the reference formulas), not substituting this
!   Jacobian into the existing Albert-based pipeline. This subroutine is
!   retained as a verified building block for that future work, but is NOT
!   currently called from WeightFunction.f90.
!
! INPUTS:
!   omega       - Resonant wave frequency [rad/s]
!   k_vec       - Wave vector magnitude [rad/m]
!   theta       - Wave normal angle [radians]
!   n_ref_ind2  - Refractive index squared (mu^2) at this (omega, theta) [-]
!   R, L, P, S  - Stix parameters at this omega [-]
!   omega_e, omega_pe, omega_i, omega_pi - Local plasma frequencies [rad/s]
!
! OUTPUTS:
!   jacobian_out - PIRAN-equivalent Jacobian value [rad*s/m^2]
!
! ==============================================================================
subroutine cunningham_jacobian(omega, k_vec, theta, n_ref_ind2, R, L, P, S, &
                               omega_e, omega_pe, omega_i, omega_pi, &
                               jacobian_out)
use constants
implicit none
    real*8, intent(in)  :: omega, k_vec, theta, n_ref_ind2
    real*8, intent(in)  :: R, L, P, S
    real*8, intent(in)  :: omega_e, omega_pe, omega_i, omega_pi
    real*8, intent(out) :: jacobian_out
 
    real*8 :: X, mu, mu2
    real*8 :: A_, B_
    real*8 :: dPdw, dSdw, dRdw, dLdw
    real*8 :: dA_, dB_, dC_
 
    X   = tan(theta)
    mu2 = n_ref_ind2
    mu  = sqrt(mu2)
 
    ! A, B (note: C is not needed by jacobian() itself, only dC)
    A_ = S*(X**2.0) + P
    B_ = R*L*(X**2.0) + P*S*(2.0 + X**2.0)
 
    call cold_plasma_stix_derivatives(omega, omega_e, omega_pe, omega_i, omega_pi, &
                                      dPdw, dSdw, dRdw, dLdw)
 
    dA_ = dSdw*(X**2.0) + dPdw
    dB_ = (dRdw*L + R*dLdw)*(X**2.0) + (dPdw*S + P*dSdw)*(2.0 + X**2.0)
    dC_ = (dPdw*R*L + P*dRdw*L + P*R*dLdw)*(1.0 + X**2.0) !same as dF3dw in cold_plasma_stix_derivatives
 
    jacobian_out = (k_vec**2.0 / (1.0 + X**2.0)) * ( &
        (dA_*mu**4.0 - dB_*mu**2.0 + dC_) &
        / (2.0*(2.0*A_*mu**4.0 - B_*mu**2.0)) &
        - 1.0/omega )
 
end subroutine cunningham_jacobian

! ==============================================================================
! SUBROUTINE: cold_plasma_parameters_multi
! ==============================================================================
! PURPOSE:
!   N-species generalization of cold_plasma_parameters, using the
!   hot_species_mass/hot_species_density_frac composition list (mass=0 taken
!   to be the electron species). Used only to seed the EMIC WHAMP root-finder
!   at theta=0 (n^2 = L exactly there)
! ==============================================================================
subroutine cold_plasma_parameters_multi(omega, omega_e_rads, ne_total, R, L, P, S, D)
use constants
implicit none
    real*8, intent(in)  :: omega          ! trial wave frequency [rad/s]
    real*8, intent(in)  :: omega_e_rads   ! electron gyrofrequency [rad/s]
    real*8, intent(in)  :: ne_total       ! total electron number density [m^-3]
    real*8, intent(out) :: R, L, P, S, D

    integer :: i_sp
    real*8  :: mass_kg, n_sp, wp_sp2, wc_sp

    R = 1.0d0
    L = 1.0d0
    P = 1.0d0

    do i_sp = 1, n_hot_species
        n_sp = hot_species_density_frac(i_sp) * ne_total
        if (hot_species_mass(i_sp) <= 0.0d0) then
            mass_kg = me
            wc_sp   = real(sqe, 8) * omega_e_rads         ! signed: -Omega_e
        else
            mass_kg = hot_species_mass(i_sp) * mi
            wc_sp   = omega_e_rads * me / mass_kg          ! signed: +Omega_ion
        end if

        wp_sp2 = n_sp * qe * qe / (C0 * mass_kg)

        R = R - wp_sp2 / (omega * (omega + wc_sp))
        L = L - wp_sp2 / (omega * (omega - wc_sp))
        P = P - wp_sp2 / (omega**2.0d0)
    end do

    S = 0.5d0 * (R + L)
    D = 0.5d0 * (R - L)

end subroutine cold_plasma_parameters_multi

! Highest ion gyrofrequency strictly below the target species -- delimits
! the target pass band from the stop band/next band below it.
subroutine emic_band_bounds(omega_e_rads, omega_target, omega_next_below)
use constants
implicit none
    real*8, intent(in)  :: omega_e_rads, omega_target
    real*8, intent(out) :: omega_next_below
    integer :: i_sp
    real*8  :: omega_sp

    omega_next_below = 0.0d0
    do i_sp = 1, n_hot_species
        if (hot_species_mass(i_sp) <= 0.0d0) cycle   ! skip electrons
        omega_sp = omega_e_rads * me / (hot_species_mass(i_sp) * mi)
        if (omega_sp < omega_target .and. omega_sp > omega_next_below) then
            omega_next_below = omega_sp
        end if
    end do
end subroutine emic_band_bounds


! Bisect L(omega)=0 between the next resonance below and the target
! resonance, to find this pass band's lower cutoff (its k=0 edge).
subroutine emic_find_cutoff(omega_e_rads, ne_total, omega_next_below, &
                             omega_target, omega_cutoff, success)
use constants
implicit none
    real*8, intent(in)  :: omega_e_rads, ne_total, omega_next_below, omega_target
    real*8, intent(out) :: omega_cutoff
    logical, intent(out) :: success

    real*8 :: w_lo, w_hi, w_mid, L_lo, L_hi, L_mid
    real*8 :: R_dum, P_dum, S_dum, D_dum
    integer :: iter
    integer, parameter :: max_iter = 60
    real*8, parameter :: eps_frac = 1.0d-4

    success = .false.
    omega_cutoff = 0.0d0

    w_lo = omega_next_below + eps_frac * (omega_target - omega_next_below)
    w_hi = omega_target - eps_frac * (omega_target - omega_next_below)

    call cold_plasma_parameters_multi(w_lo, omega_e_rads, ne_total, R_dum, L_lo, P_dum, S_dum, D_dum)
    call cold_plasma_parameters_multi(w_hi, omega_e_rads, ne_total, R_dum, L_hi, P_dum, S_dum, D_dum)

    ! Expect L<0 just above the lower resonance (evanescent) and L>0 just
    ! below the target resonance. If not, this composition doesn't give a
    ! clean single pass band here -- caller falls back.
    if (.not. (L_lo < 0.0d0 .and. L_hi > 0.0d0)) return

    do iter = 1, max_iter
        w_mid = 0.5d0 * (w_lo + w_hi)
        call cold_plasma_parameters_multi(w_mid, omega_e_rads, ne_total, R_dum, L_mid, P_dum, S_dum, D_dum)
        if (L_mid < 0.0d0) then
            w_lo = w_mid
        else
            w_hi = w_mid
        end if
    end do

    omega_cutoff = 0.5d0 * (w_lo + w_hi)
    success = .true.
end subroutine emic_find_cutoff


! Bisect L(omega) = (k*c/omega)^2 between the cutoff and the target
! resonance to get the exact cold-plasma seed frequency for a given k.
subroutine emic_find_seed_omega(k_si, omega_e_rads, ne_total, omega_cutoff, &
                                 omega_target, omega_seed, success)
use constants
implicit none
    real*8, intent(in)  :: k_si, omega_e_rads, ne_total, omega_cutoff, omega_target
    real*8, intent(out) :: omega_seed
    logical, intent(out) :: success

    real*8 :: w_lo, w_hi, w_mid, f_lo, f_hi, f_mid
    real*8 :: R_dum, L_val, P_dum, S_dum, D_dum
    integer :: iter
    integer, parameter :: max_iter = 60
    real*8, parameter :: eps_frac = 1.0d-4

    success = .false.
    omega_seed = 0.0d0

    if (k_si <= 0.0d0) then
        omega_seed = omega_cutoff
        success = .true.
        return
    end if

    w_lo = omega_cutoff
    w_hi = omega_target - eps_frac * (omega_target - omega_cutoff)

    call cold_plasma_parameters_multi(w_lo, omega_e_rads, ne_total, R_dum, L_val, P_dum, S_dum, D_dum)
    f_lo = L_val - (k_si*c/w_lo)**2.0d0
    call cold_plasma_parameters_multi(w_hi, omega_e_rads, ne_total, R_dum, L_val, P_dum, S_dum, D_dum)
    f_hi = L_val - (k_si*c/w_hi)**2.0d0

    if (.not. (f_lo <= 0.0d0 .and. f_hi >= 0.0d0)) return   ! k outside this band's range

    do iter = 1, max_iter
        w_mid = 0.5d0 * (w_lo + w_hi)
        call cold_plasma_parameters_multi(w_mid, omega_e_rads, ne_total, R_dum, L_val, P_dum, S_dum, D_dum)
        f_mid = L_val - (k_si*c/w_mid)**2.0d0
        if (f_mid < 0.0d0) then
            w_lo = w_mid
        else
            w_hi = w_mid
        end if
    end do

    omega_seed = 0.5d0 * (w_lo + w_hi)
    success = .true.
end subroutine emic_find_seed_omega
 
end module dispersion