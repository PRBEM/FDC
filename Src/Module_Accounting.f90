! ==============================================================================
! MODULE: accounting                                        (DIAGNOSTIC ONLY)
! ==============================================================================
!
! PURPOSE:
!   Writes one record per (theta, nres, root) evaluation in WeightFunction,
!   recording which guard (if any) suppressed the contribution. Comparing the
!   HOT and COLD tallies shows where the two branches start to diverge.
!
!   status: 0 accepted              | 1 invalid root (Roots_finder cutoff)
!           2 WHAMP lookup failed   | 3 dw_dk_par ~ vpar
!           4 theta > theta_max_eff | 5 theta_min > theta_max_eff
!           6 Fsq == 0
!
! USAGE:
!   The acct_write calls in WeightFunction are commented out; uncomment them
!   (and build with -DDIAG) when the tallies are needed. One file per rank:
!
!     acct_<DispersionType>_rank<rank>.txt
!
!   The context variables acct_iroot, acct_nres_ctx and acct_theta_res are
!   set by the caller before acct_write, for quantities that are not part of
!   the argument list.
!
! NOTES:
!   - The output file is opened once per rank and left open; it used to be
!     opened and closed for every record, which made the diagnostic far
!     slower than the calculation it was measuring.
!   - The file is replaced at the first write of a run. Previously it was
!     appended to, so records from earlier runs accumulated in the same file
!     and had to be separated by hand.
!   - Output is not flushed after every record, so a run that aborts may
!     lose the last few lines.
!
! ==============================================================================
module accounting
implicit none

    ! context set by the caller before acct_write
    integer :: acct_iroot     = -1
    integer :: acct_nres_ctx  = -999
    real*8  :: acct_theta_res = -1.0d0

    integer, save, private :: acct_unit = -1
    logical, save, private :: acct_open = .false.

contains

    subroutine acct_write(disp, rank, st, iroot, nres, &
                          L, E, alpha_deg, lambda_deg, theta_deg, &
                          theta_max_eff, om, n2, fs, Nom, thc, dt1)
    implicit none
        character(len=*), intent(in) :: disp
        integer,          intent(in) :: rank, st, iroot, nres
        real*8,           intent(in) :: L, E, alpha_deg, lambda_deg, theta_deg
        real*8,           intent(in) :: theta_max_eff, om, n2, fs, Nom, thc, dt1

        character(len=80) :: acct_file

        if (.not. acct_open) then
            write(acct_file,'(A,A,A,I0,A)') 'acct_', trim(disp), '_rank', rank, '.txt'
            open(newunit=acct_unit, file=trim(acct_file), status='replace', action='write')
            write(acct_unit,'(A)') '# status iroot nres  L  E  alpha_deg  lambda_deg  ' // &
                                   'theta_deg  theta_max_eff  omega  n2  Fsq  N_omega  ' // &
                                   'theta_cone  Delta_n'
            acct_open = .true.
        end if

        write(acct_unit,'(I2,1X,I2,1X,I4,6(1X,F12.5),6(1X,ES15.7))') &
            st, iroot, nres, &
            L, E, alpha_deg, lambda_deg, theta_deg, theta_max_eff, &
            om, n2, fs, Nom, thc, dt1

    end subroutine acct_write

end module accounting