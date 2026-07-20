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
!   2. RUNTIME CONFIGURATION PARAMETERS (set at startup from input file,
!      then treated as read-only during computation):
!
!      Run control:
!        PT             - Planet name
!        B0             - Equatorial surface field strength [T]
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
!        nres_min       - Minimum resonance harmonic (used if NresType='CONSTANT')
!        nres_max       - Maximum resonance harmonic (used if NresType='CONSTANT')
!
!      Plasma/wave environment:
!        Kp, Kp10       - Geomagnetic activity index
!        MLT            - Magnetic local time [hours]
!        DenModel       - Electron density model name
!        DenValue       - Density model parameter
!        DenLatIndex    - Density latitude dependence index
!
!      Wave frequency spectrum:
!        FreqUnit       - Frequency unit ('Hz' or 'rad/s')
!        FreqDistribution - Spectral shape ('Gaussian', etc.)
!        omega_m_per    - Mean wave frequency [rad/s or Hz, see FreqUnit]
!        d_omega_per    - Frequency bandwidth [rad/s or Hz]
!        o_lc, o_uc     - Lower/upper frequency cutoffs [rad/s or Hz]
!
!      Wave amplitude:
!        BwType         - Wave amplitude model ('CONSTANT' or data-driven)
!        Bw_Input       - Wave amplitude if BwType='CONSTANT' [T]
!
!      Wave normal angle distribution:
!        ThetaType      - Wave normal angle model ('CONSTANT' or data-driven)
!        Theta_Input(7) - Wave normal angle parameters if ThetaType='CONSTANT'
!                         (same layout as Theta_function throughout the code)
!
!      Integration resolution:
!        nint           - Lambda integration points (nres /= 0)
!        nint_zero_res  - Lambda integration points (nres = 0)
!        nint_X         - Theta integration points (in D_local_n)
!        nint_N         - Points for int_N_omega normalisation integral
!        nint_N_res_cone- Points for int_N_omega near resonance cone
!        lam_min, lam_max - Latitude integration bounds [degrees]
!
! NOTES:
!   - Despite its name, this module mixes true compile-time constants with
!     runtime configuration parameters. A future refactoring should separate
!     these into two modules, e.g. 'constants' (parameters only) and
!     'config' (runtime settings read from input file).
!   - sqe = -1 encodes the sign of the electron charge, used in resonance
!     condition and Stix parameter calculations throughout the code.
!
! ==============================================================================

module constants
implicit none

real*8,parameter::PI=3.1415926
real*8,parameter::RAD=PI/180.0
real*8,parameter::me=9.1093897e-31,mi=1.6726231e-27,c=2.99792458e8,C0=8.84194128288e-12,qe=1.60217733e-19,sqe=-1
integer,parameter::unit_input=10,unit_log=20,unit_output=114

character(len=30)::PT  !!planet
REAL*8::B0

character(len=100)::Output_file,log_file

integer::NL,NE,NPA
real*8::PA_min,PA_max,E_min,E_max,L_min,L_max

real*8::Kp10,Kp,MLT

character(len=30)::NresType
integer::nres_min,nres_max

character(len=30)::DenModel
real*8::DenValue,DenLatIndex

character(len=30)::DispersionType

character(len=30)::FreqUnit
character(len=30)::FreqDistribution
real*8::omega_m_per,d_omega_per,o_lc,o_uc

integer::nint,nint_zero_res,nint_X,nint_N,nint_N_res_cone

integer, parameter :: WHAMP_N_K = 80  ! number of k points in WHAMP dispersion curve
integer, parameter :: WHAMP_N_THETA = 20  ! theta points in WHAMP dispersion surface

real*8::lam_min,lam_max

character(len=30)::BwType
real*8::Bw_Input

character(len=30)::ThetaType
real*8::Theta_Input(7)

character(len=30)::DiffusionFormulation

real::time_start,time_finish

end module constants
