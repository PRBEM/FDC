! ==============================================================================
! PROGRAM: test_whamp_emic_bands
! ==============================================================================
! PURPOSE:
!   Extends test_whamp_emic_hband to all three EMIC bands (H+, He+, O+) for
!   the same L=4 plasma point and composition (H+/He+/O+ = 85/10/5%). For
!   each band, reorders the species array so the target ion sits at index 1
!   (the convention emic_band_bounds/emic_find_cutoff/emic_find_seed_omega
!   use), then runs get_whamp_dispersion and reports the curve normalized
!   to that ion's own gyrofrequency.
!
!   Also independently prints the cold cutoff/resonance bracket for each
!   band (via emic_band_bounds + emic_find_cutoff directly) so the WHAMP
!   valid/invalid boundary can be checked against the analytic prediction,
!   the same cross-check done by hand for the H+ band.
!
!   O+ band specifically exercises the omega_next_below_ion=0 fallback in
!   emic_band_bounds (O+ is the lightest-gyrofrequency ion in this
!   composition, so there is no ion resonance below it -- the band should
!   extend down toward omega=0 rather than a finite lower cutoff).
!
! BUILD:
!   gfortran Module_Constants.f90 Module_Dispersion.f90 \
!       ./tests/test_whamp_emic_bands.f90 \
!       -I../whamp/src ../whamp/src/libwhamp.a \
!       -llapack -lblas -o test_whamp_emic_bands
!   ./test_whamp_emic_bands
!
! ==============================================================================
program test_whamp_emic_bands
use constants
use dispersion
implicit none

    type(WHAMPCurve) :: curve
    integer :: i_k, i_band
    real*8  :: fce_khz, fpe_khz, omega_e_rads, ne_total
    real*8  :: omega_target, omega_next_below, omega_cutoff
    logical :: cutoff_ok
    real*8  :: freq_norm
    character(len=4) :: band_name(3)

    ! Base composition, physical fractions fixed (H+/He+/O+/e- = 85/10/5/100%)
    real*8 :: base_mass(4) = (/ 1.0d0, 4.0d0, 16.0d0, 0.0d0 /)
    real*8 :: base_frac(4) = (/ 0.85d0, 0.10d0, 0.05d0, 1.00d0 /)
    integer :: target_idx(3) = (/ 1, 2, 3 /)   ! index into base_mass/base_frac for H+, He+, O+

    band_name = (/ 'H+  ', 'He+ ', 'O+  ' /)

    ! L=4 equatorial values, same as the H+ band test
    fce_khz = 13.1133d0
    fpe_khz = 429.47d0
    omega_e_rads = 2.0d0*PI*fce_khz*1.0d3
    ne_total = (2.0d0*PI*fpe_khz*1.0d3)**2 * C0 * me / (qe*qe)

    do i_band = 1, 3
        write(*,*) '================================================'
        write(*,*) 'Band: ', trim(band_name(i_band))
        write(*,*) '================================================'

        ! Reorder so the target ion is species 1; everything else keeps its
        ! own mass/fraction, just shifted to fill the remaining slots.
        call build_species_list(target_idx(i_band), base_mass, base_frac)

        ! Independent check of the analytic band bounds, same computation
        ! get_whamp_dispersion will do internally.
        omega_target = omega_e_rads * me / (hot_species_mass(1) * mi)
        call emic_band_bounds(omega_e_rads, omega_target, omega_next_below)
        call emic_find_cutoff(omega_e_rads, ne_total, omega_next_below, &
                               omega_target, omega_cutoff, cutoff_ok)

        write(*,'(A,F10.6,A)') '  Omega_target      = ', omega_target/(2.0d0*PI), ' Hz'
        write(*,'(A,F10.6,A)') '  Omega_next_below  = ', omega_next_below/(2.0d0*PI), ' Hz'
        if (cutoff_ok) then
            write(*,'(A,F10.6,A)') '  cold cutoff        = ', omega_cutoff/omega_target, &
                                    ' (fraction of Omega_target)'
        else
            write(*,'(A)') '  cold cutoff bisection FAILED (unexpected bracket)'
        end if
        if (omega_next_below == 0.0d0) then
            write(*,'(A)') '  (no ion resonance below target -- fallback path exercised)'
        end if

        DispersionType = 'EMIC'
        k_min_exp_val_input = -6.0d0
        k_max_exp_val_input =  0.0d0

        call get_whamp_dispersion(fce_khz, fpe_khz, 0.0d0, 0.0d0, curve)

        write(*,'(A)') '  i_k    k_WHAMP    omega [Hz]   omega/Omega_target'
        write(*,'(A)') '  ' // repeat('-', 50)
        do i_k = 1, curve%n_k
            if (curve%valid(1, i_k)) then
                freq_norm = curve%omega(1, i_k) / omega_target
                write(*,'(I5, F11.5, F14.4, F16.5)') i_k, 10.0d0**curve%k_whamp(i_k), &
                    curve%omega(1, i_k) / (2.0d0*PI), freq_norm
            end if
        end do
        write(*,'(A,I4,A,I4)') '  Valid points: ', count(curve%valid), ' / ', curve%n_k
        write(*,*)

        deallocate(curve%theta, curve%k, curve%k_whamp, curve%omega, &
                   curve%dw_dk, curve%valid, curve%EFL)
    end do

contains

    subroutine build_species_list(target_i, mass_arr, frac_arr)
        implicit none
        integer, intent(in) :: target_i
        real*8, intent(in)  :: mass_arr(4), frac_arr(4)
        integer :: j, slot

        n_hot_species = 4
        hot_species_mass(1) = mass_arr(target_i)
        hot_species_density_frac(1) = frac_arr(target_i)
        hot_species_T_par_keV(1) = 0.001d0
        hot_species_T_anis(1) = 1.0d0

        slot = 2
        do j = 1, 4
            if (j == target_i) cycle
            hot_species_mass(slot) = mass_arr(j)
            hot_species_density_frac(slot) = frac_arr(j)
            hot_species_T_par_keV(slot) = 0.001d0
            hot_species_T_anis(slot) = 1.0d0
            slot = slot + 1
        end do
    end subroutine build_species_list

end program test_whamp_emic_bands