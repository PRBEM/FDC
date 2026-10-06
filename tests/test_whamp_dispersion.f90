! ==============================================================================
! PROGRAM: test_whamp_dispersion
! ==============================================================================
!
! PURPOSE:
!   Standalone test for get_whamp_dispersion. Calls WHAMP for alpha*=0.01
!   at L=3 (fce=32.4 kHz, fpe=101.2 kHz) and prints the dispersion curve
!   for parallel propagation (theta=0).
!
!   Expected output should match interactive WHAMP run with:
!     plasma_hot_01.txt, f=0.5, p=0,0,1, z=-3,1,0.2
!
! BUILD:
!   gfortran Module_Constants.f90 Module_Variables.f90 Module_Dispersion.f90 \
!       test_whamp_dispersion.f90 \
!       -I../whamp/src ../whamp/src/libwhamp.a \
!       -llapack -lblas -o test_whamp
!   ./test_whamp
!
! ==============================================================================

program test_whamp_dispersion
use constants
use dispersion
implicit none

    type(WHAMPCurve) :: curve
    integer :: i_k
    real*8  :: fce_khz, fpe_khz
    real*8  :: k_norm, freq_norm

    ! alpha*=0.01 at L=3: fce=32.4 kHz, fpe = fce/sqrt(alpha*) = 32.4/sqrt(0.01) = 324 kHz
    fce_khz = 32.4d0
    fpe_khz = 324.0d0

    write(*,*) '================================================'
    write(*,*) 'Test: get_whamp_dispersion'
    write(*,*) 'fce =', fce_khz, 'kHz'
    write(*,*) 'fpe =', fpe_khz, 'kHz'
    write(*,*) 'alpha* = fce^2/fpe^2 =', (fce_khz/fpe_khz)**2
    write(*,*) '================================================'
    write(*,*)

    ! Call with parallel propagation only (theta_min=theta_max=0)
    call get_whamp_dispersion(fce_khz, fpe_khz, 0.0d0, 0.0d0, curve)

    write(*,'(A)') 'Dispersion curve (parallel propagation, theta=0):'
    write(*,'(A6, A12, A12, A14, A12, A8)') &
        'i_k', 'k_WHAMP', 'k*c/wpe', 'omega [kHz]', 'w/Omega_e', 'valid'
    write(*,'(A)') repeat('-', 65)

    do i_k = 1, curve%n_k
        ! k normalised to c/omega_pe (Cao et al. convention)
        k_norm  = 10.0d0**curve%k_whamp(i_k) * (fce_khz/fpe_khz)
        ! frequency normalised to Omega_e
        freq_norm = curve%omega(1, i_k) / (2.0d0 * PI * fce_khz * 1.0d3)

        write(*,'(I6, F12.5, F12.5, F14.4, F12.5, L8)') &
            i_k, &
            10.0d0**curve%k_whamp(i_k), &
            k_norm, &
            curve%omega(1, i_k) / (2.0d0 * PI * 1.0d3), &
            freq_norm, &
            curve%valid(1, i_k)
    end do

    write(*,*)
    write(*,'(A,I4,A,I4,A)') 'Valid points: ', &
        count(curve%valid), ' / ', curve%n_k*curve%n_theta, ' total'

    ! Cleanup
    deallocate(curve%theta, curve%k, curve%k_whamp, &
               curve%omega, curve%dw_dk, curve%valid)

end program test_whamp_dispersion