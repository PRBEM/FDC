! ==============================================================================
! SUBROUTINE: dump_whamp_surface   (DIAGNOSTIC ONLY)
! ==============================================================================
!
! PURPOSE:
!   Writes the whole tabulated WHAMP dispersion surface omega(theta, k) for one
!   cached curve, so the branch structure can be inspected directly rather than
!   inferred from the interpolated results in WeightFunction / int_N_omega.
!
!   Use this to test whether a saturation or kink in n^2(theta) is a genuine
!   turnover of the kinetic branch or an artefact of the k grid running out:
!   if the root at the target frequency sits at the last valid k index, the
!   plateau is the grid edge, not physics.
!
! USAGE:
!   Called from Integration_lambda, once per process, only in builds with
!   -DDIAG. Each rank writes its own file.
!
! OUTPUT FILE:
!   whampsurf_<tag>_L<L>_lam<lambda*1000 deg>_rank<rank>.txt
!
! OUTPUT COLUMNS:
!   i_theta  i_k  theta_deg  k_SI  k_whamp(log10)  omega  n2  dw_dk  valid
!
!   n2 = (k*c/omega)^2, set to -1 where the point did not converge.
!
! ==============================================================================
subroutine dump_whamp_surface(L_shell, lambda, curve, tag)
use constants
use dispersion
implicit none
    include "mpif.h"

    real*8,           intent(in) :: L_shell, lambda
    type(WHAMPCurve), intent(in) :: curve
    character(len=*), intent(in) :: tag

    integer            :: i_t, i_k, lam_tag, unit_surf, my_rank, ierr
    real*8             :: n2
    character(len=128) :: fname

    ! NOTE: the intrinsic NINT is shadowed by the integer variable `nint`
    ! in Module_Constants, so round explicitly instead.
    lam_tag = int(lambda * 180.0d0 / PI * 1000.0d0 + 0.5d0)

    ! One file per rank: without the rank in the name, all processes would
    ! write to the same file (see the concurrent-write problem in the other
    ! per-rank diagnostics).
    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

    write(fname,'(A,A,A,I0,A,I0,A,I3.3,A)') 'whampsurf_', trim(tag), '_L', &
         int(L_shell), '_lam', lam_tag, '_rank', my_rank, '.txt'

    open(newunit=unit_surf, file=trim(fname), status='replace')
    write(unit_surf,'(A)') '# i_theta i_k theta_deg k_SI k_whamp omega n2 dw_dk valid'
    write(unit_surf,'(A,I0,A,I0,A,ES14.6,A,ES14.6)') '# n_theta=', curve%n_theta, &
         ' n_k=', curve%n_k, ' fce_kHz=', curve%fce_khz, ' fpe_kHz=', curve%fpe_khz

    do i_t = 1, curve%n_theta
        do i_k = 1, curve%n_k
            if (curve%valid(i_t, i_k) .and. curve%omega(i_t, i_k) > 0.0d0) then
                n2 = (curve%k(i_k) * c / curve%omega(i_t, i_k))**2
            else
                n2 = -1.0d0
            end if
            write(unit_surf,'(I4,1X,I4,1X,F10.5,4(1X,ES16.8),1X,ES14.6,1X,L2)') &
                i_t, i_k, curve%theta(i_t)*180.0d0/PI, &
                curve%k(i_k), curve%k_whamp(i_k), &
                curve%omega(i_t, i_k), n2, curve%dw_dk(i_t, i_k), &
                curve%valid(i_t, i_k)
        end do
    end do
    close(unit_surf)

end subroutine dump_whamp_surface