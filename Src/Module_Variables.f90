! ==============================================================================
! MODULE: variables
! ==============================================================================
!
! PURPOSE:
!   Provides shared mutable state for plasma frequencies and wave band
!   parameters that are computed in WeightFunction and read by Roots_finder.
!
!   PLASMA FREQUENCIES (written by WeightFunction, read by Roots_finder):
!     omega_e    - Electron gyrofrequency at local point [rad/s]
!     omega_pe   - Electron plasma frequency at local point [rad/s]
!     omega_i    - Ion gyrofrequency at local point [rad/s]
!     omega_pi   - Ion plasma frequency at local point [rad/s]
!
!   WAVE BAND PARAMETERS (written by WeightFunction via Freq_Spectral,
!   read by Roots_finder for frequency validity checks):
!     omega_m    - Mean wave frequency [rad/s]
!     d_omega    - Wave frequency bandwidth [rad/s]
!     omega_lc   - Lower frequency cutoff [rad/s]
!     omega_uc   - Upper frequency cutoff [rad/s]
!
! NOTES:
!   - The Stix parameters (R_stix, L_stix, P_stix, S_stix, D_stix) and
!     polynomial coefficients (a1-a6) have been removed from this module. 
!     They are now computed locally within
!     Module_Dispersion.f90 and no longer need to be shared state.
!   - The remaining variables (plasma frequencies and wave band parameters)
!     are still shared state because Roots_finder reads them implicitly
!     rather than receiving them as arguments. Removing this remaining
!     dependency would require refactoring Roots_finder to accept these
!     as explicit arguments.
!
! ==============================================================================

module variables
implicit none

    ! Plasma frequencies. Written by WeightFunction, read by Roots_finder
    real*8 :: omega_e, omega_pe, omega_i, omega_pi

    ! Wave band parameters. Written by WeightFunction via Freq_Spectral,
    ! read by Roots_finder for frequency validity checks
    real*8 :: omega_lc, omega_uc

    ! omega_m and d_omega are returned by Freq_Spectral as explicit arguments
    ! in WeightFunction and are no longer needed as shared module state.
    ! They are kept here temporarily. Remove once confirmed unused.
    real*8 :: omega_m, d_omega

end module variables