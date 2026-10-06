! ==============================================================================
! MODULE: constants
! ==============================================================================
!
! PURPOSE:
!   Provides two categories of data used throughout the code:
!
!   1. TRUE PHYSICAL CONSTANTS (parameter, never changed at runtime):
!        PI, RAD        - mathematical constants
!        me, mi         - electron and proton rest mass [kg]
!        c              - speed of light [m/s]
!        C0             - permittivity of free space [F/m]
!        qe             - elementary charge [C]
!        sqe            - sign of electron charge (-1)
!
!   2. RUNTIME CONFIGURATION PARAMETERS (set at startup from the input file
!      in FDC_MPI.f90, then treated as read-only during computation):
!
!      Run control:
!        PT             - Planet name
!        B0             - Equatorial field strength [T]; for PT = 'DIPOLE'
!                         this is the field AT the L-shell of interest, not
!                         at the planetary surface (see B_field)
!        Output_file    - Output filename
!        log_file       - Log filename
!        unit_input/log/output - Fortran I/O unit numbers
!
!      Particle grid:
!        NL, NE, NPA    - Number of L, energy, pitch angle grid points
!        L_min, L_max   - L-shell range [Earth radii]
!        E_min, E_max   - Energy range [keV]
!        PA_min, PA_max - Pitch angle range [degrees]
!
!      Resonance control:
!        NresType       - Resonance order mode: 'CONSTANT' or 'MR12'
!        nres_min       - Minimum resonance harmonic (if NresType='CONSTANT')
!        nres_max       - Maximum resonance harmonic (if NresType='CONSTANT')
!
!      Plasma/wave environment:
!        Kp, Kp10       - Geomagnetic activity index
!        MLT            - Magnetic local time [hours]
!        DenModel       - Electron density model name
!        DenValue       - Density model parameter
!        DenLatIndex    - Density latitude dependence index
!
!      Hot plasma (WHAMP) species (required when DispersionType='HOT'):
!        n_hot_species            - Number of plasma species configured
!        hot_species_mass         - ASS value per species (0=e-, 1=H+, 4=He+, 16=O+)
!        hot_species_density_frac - Density fraction of ne_total, per species
!        hot_species_T_par_keV    - Parallel temperature [keV], per species
!        hot_species_T_anis       - T_perp/T_par anisotropy, per species
!        k_min_exp_val_input      - Lower bound of WHAMP log10(k) grid
!        k_max_exp_val_input      - Upper bound of WHAMP log10(k) grid
!        WHAMP_N_K, WHAMP_N_THETA - size of the WHAMP dispersion surface
!        N_max_hot                - hot electron density [cm^-3] used in the
!                                   WHAMP setup
!
!      Wave frequency spectrum:
!        FreqUnit         - Frequency unit ('Hz' or 'rad/s')
!        FreqDistribution - Spectral shape ('Gaussian', etc.)
!        omega_m_per      - Mean wave frequency [rad/s or Hz, see FreqUnit]
!        d_omega_per      - Frequency bandwidth [rad/s or Hz]
!        o_lc, o_uc       - Lower/upper frequency cutoffs [rad/s or Hz]
!
!      Wave amplitude:
!        BwType         - Wave amplitude model ('CONSTANT' or data-driven)
!        Bw_Input       - Wave amplitude if BwType='CONSTANT' [pT]
!
!      Wave normal angle distribution:
!        ThetaType      - Wave normal angle model ('CONSTANT' or data-driven)
!        Theta_Input(7) - Wave normal angle parameters if ThetaType='CONSTANT'
!                         (same layout as Theta_function throughout the code)
!
!      Integration resolution and mode:
!        nint            - Lambda integration points (nres /= 0)
!        nint_zero_res   - Lambda integration points (nres = 0)
!        nint_X          - Theta integration points (in D_local_n)
!        nint_N          - Points for the int_N_omega normalisation integral
!        nint_N_res_cone - Points for int_N_omega near the resonance cone
!        lam_min, lam_max - Latitude integration bounds [degrees]
!        local_mode      - local (non-bounce-averaged) mode, see below
!        cone_backoff, ConeTruncateN - resonance cone handling, see below
!
! NOTES:
!   - Despite its name, this module mixes true compile-time constants with
!     runtime configuration parameters. A future refactoring should separate
!     these into two modules, e.g. 'constants' (parameters only) and
!     'config' (runtime settings read from the input file).
!   - sqe = -1 encodes the sign of the electron charge, used in the resonance
!     condition and in the Stix parameters throughout the code.
!   - The physical constants are written with d0 exponents. Without them the
!     literals would be default (single precision) reals, so a real*8
!     parameter would carry only about seven significant digits.
!
! ==============================================================================

module constants
implicit none

real*8,parameter::PI=3.14159265358979323846d0
real*8,parameter::RAD=PI/180.0d0
real*8,parameter::me=9.1093897d-31,mi=1.6726231d-27,c=2.99792458d8, &
                  C0=8.8541878188d-12,qe=1.60217733d-19,sqe=-1
integer,parameter::unit_input=10,unit_log=20,unit_output=114

character(len=30)::PT              ! planet
real*8::B0                         ! equatorial field [T]; see note above

character(len=100)::Output_file,log_file

integer::NL,NE,NPA
real*8::PA_min,PA_max,E_min,E_max,L_min,L_max

real*8::Kp10,Kp,MLT

character(len=30)::NresType
integer::nres_min,nres_max

character(len=30)::DenModel
real*8::DenValue,DenLatIndex

character(len=30)::DispersionType   ! COLD, HOT or EMIC

! Hot plasma (WHAMP) species configuration. Required when DispersionType='HOT';
! populated from the input file in FDC_MPI.f90 (no default -- see there).
integer, parameter :: max_hot_species = 10   ! matches comin's ASS/TA/AA/DN/DD sizing
integer :: n_hot_species
real*8  :: hot_species_mass(max_hot_species)         ! ASS value: 0=electron, 1=H+, 4=He+, 16=O+
real*8  :: hot_species_density_frac(max_hot_species) ! fraction of ne_total, per species
real*8  :: hot_species_T_par_keV(max_hot_species)
real*8  :: hot_species_T_anis(max_hot_species)       ! T_perp/T_par
real*8  :: k_min_exp_val_input
real*8  :: k_max_exp_val_input

character(len=30)::FreqUnit
character(len=30)::FreqDistribution
real*8::omega_m_per,d_omega_per,o_lc,o_uc

integer::nint,nint_zero_res,nint_X,nint_N,nint_N_res_cone

! Resonance cone handling in the N(omega) normalisation integral
! (WeightFunction / int_N_omega):
!   ConeTruncateN = .true.  : the theta integration stops at
!                             cone_backoff * theta_cone when the resonance
!                             cone falls inside [theta_min, theta_max], and
!                             the finer grid nint_N_res_cone is used there.
!                 = .false. : the integral runs to theta_max regardless.
real*8,  parameter :: cone_backoff   = 0.99d0
logical, parameter :: ConeTruncateN  = .true.

integer, parameter :: WHAMP_N_K = 80       ! k points in the WHAMP dispersion curve
integer, parameter :: WHAMP_N_THETA = 40   ! theta points in the dispersion surface
real*8,  parameter :: N_max_hot = 150.0d0  ! test value; Li et al. 2014 range: 130-280

real*8::lam_min,lam_max
! Set in FDC_MPI when lam_min = lam_max = 0: D_ba then returns the LOCAL
! equatorial coefficients, without bounce averaging and without the
! 1/S(alpha) factor. Any other lam_min = lam_max is rejected there.
logical::local_mode = .false.

character(len=30)::BwType
real*8::Bw_Input                   ! wave amplitude [pT] if BwType = 'CONSTANT'

character(len=30)::ThetaType
real*8::Theta_Input(7)

character(len=30)::DiffusionFormulation   ! LYONS or CUNNINGHAM

real::time_start,time_finish

end module constants