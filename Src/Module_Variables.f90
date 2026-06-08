! ==============================================================================
! MODULE: variables
! ==============================================================================
!
! PURPOSE:
!   Previously provided shared mutable state for plasma frequencies and wave
!   band parameters that were computed in WeightFunction and read by
!   Roots_finder.
!
!   As of the refactoring for hot plasma, all variables have been removed:
!     - Stix parameters (R_stix etc.) and a1-a6: removed,
!       now computed locally in Module_Dispersion.f90
!     - omega_e, omega_pe, omega_i, omega_pi, omega_lc, omega_uc: removed,
!       now passed as explicit arguments to Roots_finder
!     - omega_m, d_omega: now local in WeightFunction
!
!   This module is now empty and is retained only to avoid breaking any
!   subroutines that still have 'use variables' in their source. It can
!   be deleted once all 'use variables' statements have been removed from
!   the codebase.
!
! ==============================================================================

module variables
implicit none

    ! This module is intentionally empty.
    ! See header for history of removed variables.

end module variables