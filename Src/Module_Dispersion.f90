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
    ! real*8 :: ex, ey, ez  ! normalised electric field components
    ! real*8 :: bx, by, bz  ! normalised magnetic field components
 
end type DispersionProperties
 
 
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
            ! HOT plasma dispersion via WHAMP (not yet implemented).
            ! When implemented, this case should populate DP%n_ref_ind2,
            ! DP%dw_dk, and the WHAMP field component fields (ex, ey, ez,
            ! bx, by, bz), then set DP%disp_type = 'HOT'.
            write(*,*) 'ERROR: DispersionType = HOT is not yet implemented.'
            write(*,*) 'Set DispersionType = COLD in the input file, or '
            write(*,*) 'implement the WHAMP interface in Module_Dispersion.f90.'
            stop
 
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
            ! To be implemented using WHAMP field components (ex,ey,ez,bx,by,bz)
            write(*,*) 'ERROR: HOT plasma polarisation factor not yet implemented.'
            stop
 
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
 
 
end module dispersion
