! ==============================================================================
! PROGRAM: test_whamp_emic_hband
! ==============================================================================
!
! PURPOSE:
!   Standalone test for get_whamp_dispersion with an EMIC H+ band multi-ion
!   species configuration, parallel propagation only (theta=0). Prints the
!   resulting omega(k) curve so it can be checked against Albert (2003)
!   Figure 1's pass-band structure BEFORE running the full diffusion
!   coefficient grid (which is much more expensive to iterate on if the
!   k_min_exp/k_max_exp guess turns out to be wrong).
!
!   Plasma point: L=4, n_e=2285 cm^-3 (paper's stated value for
!   omega_pe^2/Omega_e^2 = 10^3 at L=4), H+/He+/O+ = 85%/10%/5%.
!
! BUILD:
!   gfortran Module_Constants.f90 Module_Dispersion.f90 \
!       test_whamp_emic_hband.f90 \
!       -I../whamp/src ../whamp/src/libwhamp.a \
!       -llapack -lblas -o test_whamp_emic
!   ./test_whamp_emic
!
! ==============================================================================
program test_whamp_emic_hband
use constants
use dispersion
implicit none
    type(WHAMPCurve) :: curve
    integer :: i_k
    real*8  :: fce_khz, fpe_khz
    real*8  :: Omega_H_plus_khz, freq_norm_Hp

    ! L=4 equatorial values matching the input file
    fce_khz = 13.1133d0
    fpe_khz = 429.47d0

    ! H+ band composition (Albert 2003 Sec. 4, sudden commencement phase)
    n_hot_species = 4
    hot_species_mass(1) = 1.0d0;  hot_species_density_frac(1) = 0.85d0
    hot_species_T_par_keV(1) = 0.001d0; hot_species_T_anis(1) = 1.0d0
    hot_species_mass(2) = 4.0d0;  hot_species_density_frac(2) = 0.10d0
    hot_species_T_par_keV(2) = 0.001d0; hot_species_T_anis(2) = 1.0d0
    hot_species_mass(3) = 16.0d0; hot_species_density_frac(3) = 0.05d0
    hot_species_T_par_keV(3) = 0.001d0; hot_species_T_anis(3) = 1.0d0
    hot_species_mass(4) = 0.0d0;  hot_species_density_frac(4) = 1.00d0   ! electrons
    hot_species_T_par_keV(4) = 0.001d0; hot_species_T_anis(4) = 1.0d0

    DispersionType = 'EMIC'

    k_min_exp_val_input = -6.0d0
    k_max_exp_val_input =  0.0d0

    write(*,*) '================================================'
    write(*,*) 'Test: get_whamp_dispersion (EMIC H+ band, parallel)'
    write(*,*) 'fce =', fce_khz, 'kHz'
    write(*,*) 'fpe =', fpe_khz, 'kHz'
    write(*,*) '(fpe/fce)^2 =', (fpe_khz/fce_khz)**2, ' (paper states ~1000 at L=4)'
    write(*,*) '================================================'
    write(*,*)

    call get_whamp_dispersion(fce_khz, fpe_khz, 0.0d0, 0.0d0, curve)

    Omega_H_plus_khz = fce_khz * me / mi   ! H+ mass = 1 proton mass

    write(*,'(A)') 'Dispersion curve (parallel propagation, theta=0):'
    write(*,'(A6, A14, A16, A14)') 'i_k', 'k_WHAMP', 'omega [kHz]', 'omega/Omega_H+'
    write(*,'(A)') repeat('-', 55)
    do i_k = 1, curve%n_k
        if (curve%valid(1, i_k)) then
            freq_norm_Hp = curve%omega(1, i_k) / (2.0d0*PI*1.0d3) / Omega_H_plus_khz
            write(*,'(I6, F14.5, F16.4, F14.5)') i_k, 10.0d0**curve%k_whamp(i_k), &
                curve%omega(1, i_k) / (2.0d0*PI*1.0d3), freq_norm_Hp
        end if
    end do

    write(*,*)
    write(*,'(A,I4,A,I4,A)') 'Valid points: ', count(curve%valid), ' / ', curve%n_k, ' total'
    write(*,'(A)') 'Expect a pass band clustering near omega/Omega_H+ ~ 0.5-0.7 ' // &
                   '(paper H-band: omega_m=0.6, d_omega=0.1 * Omega_H+)'

    deallocate(curve%theta, curve%k, curve%k_whamp, curve%omega, &
               curve%dw_dk, curve%valid, curve%EFL)

end program test_whamp_emic_hband
