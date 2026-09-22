! ==============================================================================
! MODULE: accounting
!
! Diagnostic-only. Writes one record per (theta, nres, root) evaluation in
! WeightFunction, recording which guard (if any) suppressed the contribution.
!
!   status: 0 accepted             | 1 invalid root (Roots_finder cutoff)
!           2 WHAMP lookup failed  | 3 dw_dk_par ~ vpar
!           4 theta > theta_max_eff| 5 theta_min > theta_max_eff
!           6 Fsq == 0
!
! Compare HOT vs COLD tallies to locate where the two branches diverge.
! Remove or comment out the calls once the question is settled -- the
! open/close per record is slow.
! ==============================================================================
module accounting
implicit none

integer :: acct_iroot     = -1
integer :: acct_nres_ctx  = -999
real*8  :: acct_theta_res = -1.0d0
contains

    subroutine acct_write(disp, rank, st, iroot, nres, &
                          L, E, alpha_deg, lambda_deg, theta_deg, &
                          theta_max_eff, om, n2, fs, Nom, thc, dt1)
    implicit none
        character(len=*), intent(in) :: disp
        integer, intent(in) :: rank, st, iroot, nres
        real*8,  intent(in) :: L, E, alpha_deg, lambda_deg, theta_deg
        real*8,  intent(in) :: theta_max_eff, om, n2, fs, Nom, thc, dt1

        character(len=80) :: acct_file

        write(acct_file,'(A,A,A,I0,A)') 'acct_', trim(disp), '_rank', rank, '.txt'
        open(96, file=trim(acct_file), status='unknown', position='append')
        write(96,'(I2,1X,I2,1X,I4,6(1X,F12.5),6(1X,ES15.7))') &
            st, iroot, nres, &
            L, E, alpha_deg, lambda_deg, theta_deg, theta_max_eff, &
            om, n2, fs, Nom, thc, dt1
        close(96)

    end subroutine acct_write

end module accounting