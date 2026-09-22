! ==============================================================================
! SUBROUTINE: mpi_parallel
! ==============================================================================
!
! PURPOSE:
!   Top-level MPI parallelisation routine. Constructs the full (L, E, alpha)
!   computation grid, distributes grid points across MPI processes, collects
!   results, and writes the output file.
!
!   Uses a master/slave pattern:
!     Master (rank 0):
!       1. Constructs the full input grid of (L, E, alpha) points
!       2. Opens the output file and writes the header
!       3. Sends grid points one at a time to slave processes
!       4. Receives computed diffusion coefficients from slaves
!       5. Writes results to the output file
!       6. Closes the output file
!
!     Slaves (ranks 1 to numprocs-1):
!       1. Receive a grid point (L, E, alpha) from master
!       2. Call sum_n to compute bounce-averaged diffusion coefficients
!       3. Send results back to master
!       4. Repeat for all assigned grid points
!
!   Grid construction:
!     L     : uniform spacing from L_min to L_max, NL points
!     E     : logarithmic spacing from E_min to E_max, NE points [keV]
!     alpha : uniform spacing from PA_min to PA_max, NPA points [degrees]
!
!   Output file format: Tecplot ASCII with columns:
!     L [Earth radii], E [keV], alpha [degrees],
!     Daa [day^-1], Dpp [day^-1], Dap [day^-1]
!
! CALLS:
!   sum_n       - Computes bounce-averaged diffusion coefficients for one
!                 (L, E, alpha) point, summing over resonance harmonics
!
! MODULE DEPENDENCIES:
!   constants   - NL, NE, NPA, L_min, L_max, E_min, E_max, PA_min, PA_max,
!                 unit_output, unit_log, Output_file
!
! NOTES:
!   - PERFORMANCE: currently sends one grid point per MPI message
!     (noutput_rows=1). For large grids (~1.1M points) this generates
!     ~2.2M messages, which has significant communication overhead.
!     A future improvement would be to batch multiple points per message.
!     The commented-out MPI I/O block at the top of the original was an
!     early attempt to address this.
!   - The grid size is hardcoded to a maximum of 1,200,000 points
!     (InputData(1200000,3)), sufficient for 89*211*61 ~ 1.1M points.
!     If NL*NE*NPA exceeds this, the code will silently overflow.
!     A future improvement would use allocatable arrays sized from NL*NE*NPA.
!   - The load is distributed so that remainder rows (mod(rows, numprocs-1))
!     are assigned to the lowest-ranked slaves (ranks 1 to row_remainder).
!
! ==============================================================================

subroutine mpi_parallel
use constants
use dispersion
implicit none

    include "mpif.h"

    ! Grid dimensions
    integer :: rows
    integer, parameter :: columns = 3   ! L, Energy, pitch angle

    ! Grid point temporaries
    real*8 :: Ltemp, Etemp, PAtemp

    ! MPI partitioning
    integer   :: local_rows, start_row, end_row
    integer   :: noutput_rows, noutput_columns
    integer*8 :: row_number, row_remainder

    ! Data arrays
    real*8, target      :: InputData(1200000, columns)
    real*8, allocatable :: LocalData(:,:)
    real*8, pointer     :: SendData(:,:)
    real*8, allocatable :: OutputData(:,:), RecvData(:)

    ! MPI variables
    integer :: ierr, numprocs, mype, master
    integer :: status(MPI_STATUS_SIZE)
    integer :: i, j, k

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

    rows   = NL * NE * NPA
    master = 0

    ! -------------------------------------------------------------------------
    ! MASTER: build grid, distribute work, collect and write results
    ! -------------------------------------------------------------------------
    if (mype == master) then

        ! Build input grid
        ! L: uniform spacing; E: logarithmic spacing; alpha: uniform spacing
        do i = 1, NL
            if (NL == 1) then
                Ltemp = L_min
            else
                Ltemp = (L_max - L_min) / (NL - 1.0) * (i - 1.0) + L_min
            end if

            do j = 1, NE
                if (NE == 1) then
                    Etemp = E_min
                else
                    Etemp = 10**((log10(E_max) - log10(E_min)) &
                                 / (NE - 1.0) * (j - 1.0) + log10(E_min))
                end if

                do k = 1, NPA
                    PAtemp = (PA_max - PA_min) / (NPA - 1.0) * (k - 1.0) + PA_min
                    InputData(((i-1)*NE + (j-1))*NPA + k, 1) = Ltemp
                    InputData(((i-1)*NE + (j-1))*NPA + k, 2) = Etemp
                    InputData(((i-1)*NE + (j-1))*NPA + k, 3) = PAtemp
                end do
            end do
        end do

        ! Open output file and write Tecplot header
        open(unit=unit_output, file=Output_file, status='replace')
        write(unit_output, "('VARIABLES = ""L"", ""Energy"", ""Pitch-angle"", " // &
              """Daa Day^-1"", ""Dpp Day^-1"",""Dap Day^-1""')")
        write(unit_output, "(/)")
        write(unit_output, "('Zone T=""some name"", I=' I3 ',' 1X 'J=' I3, ',' 1X 'K=' I3)") &
              NPA, NE, NL
        write(unit_output, "(/)")

        ! Partition and distribute grid points to slaves
        row_remainder = mod(rows, numprocs-1)
        row_number    = (rows - row_remainder) / (numprocs - 1d0)
        noutput_rows    = 1
        noutput_columns = 6

        allocate(RecvData(noutput_rows * noutput_columns))

        write(unit_log, "(I6 4x I6 4X I6 4x I6)") &
              rows, numprocs-1, row_number, row_remainder

        ! Send/receive main block of rows
        do k = 1, row_number
            do i = 1, numprocs-1
                start_row = i + (k-1.0) * (numprocs-1d0)
                end_row   = start_row
                SendData => InputData(start_row:end_row, 1:columns)
                call MPI_SEND(SendData, (end_row-start_row+1)*columns, &
                              MPI_DOUBLE_PRECISION, i, 110+i, MPI_COMM_WORLD, ierr)
            end do

            do i = 1, numprocs-1
                call MPI_RECV(RecvData, noutput_columns*noutput_rows, &
                              MPI_DOUBLE_PRECISION, i, 110+i, MPI_COMM_WORLD, status, ierr)
                do j = 1, noutput_rows
                    write(unit_output, &
                          "(f6.3,5X,E15.6E4,5X,f5.2,5X,E15.6E4,5X,E15.6E4,5X,E15.6E4)") &
                          RecvData(j + noutput_rows*0), &
                          RecvData(j + noutput_rows*1), &
                          RecvData(j + noutput_rows*2), &
                          RecvData(j + noutput_rows*3), &
                          RecvData(j + noutput_rows*4), &
                          RecvData(j + noutput_rows*5)
                end do
            end do
        end do

        ! Send/receive remainder rows
        do i = 1, row_remainder
            start_row = i + row_number * (numprocs-1d0)
            end_row   = start_row
            SendData => InputData(start_row:end_row, 1:columns)
            call MPI_SEND(SendData, (end_row-start_row+1)*columns, &
                          MPI_DOUBLE_PRECISION, i, 110+i, MPI_COMM_WORLD, ierr)
        end do

        do i = 1, row_remainder
            call MPI_RECV(RecvData, noutput_columns*noutput_rows, &
                          MPI_DOUBLE_PRECISION, i, 110+i, MPI_COMM_WORLD, status, ierr)
            do j = 1, noutput_rows
                write(unit_output, &
                      "(f6.3,5X,E15.6E4,5X,f5.2,5X,E15.6E4,5X,E15.6E4,5X,E15.6E4)") &
                      RecvData(j + noutput_rows*0), &
                      RecvData(j + noutput_rows*1), &
                      RecvData(j + noutput_rows*2), &
                      RecvData(j + noutput_rows*3), &
                      RecvData(j + noutput_rows*4), &
                      RecvData(j + noutput_rows*5)
            end do
        end do

        close(unit_output)

    ! -------------------------------------------------------------------------
    ! SLAVES: receive grid points, compute diffusion coefficients, send back
    ! -------------------------------------------------------------------------
    else

        row_remainder = mod(rows, numprocs-1)
        row_number    = (rows - row_remainder) / (numprocs - 1d0)
        noutput_rows    = 1
        noutput_columns = 6

        ! Slaves with rank <= row_remainder get one extra point
        if (mype <= row_remainder) row_number = row_number + 1

        allocate(LocalData(1, columns))
        allocate(OutputData(noutput_rows, noutput_columns))

        do k = 1, row_number
            call MPI_RECV(LocalData, columns, MPI_DOUBLE_PRECISION, master, &
                          110+mype, MPI_COMM_WORLD, status, ierr)

            call sum_n(LocalData, 1, columns, OutputData, noutput_rows, noutput_columns)

            call MPI_SEND(OutputData, noutput_rows*noutput_columns, &
                          MPI_DOUBLE_PRECISION, master, 110+mype, MPI_COMM_WORLD, ierr)
        end do

    end if
#ifdef PROFILE
    write(*,*) '=== PROFILING SUMMARY (rank', mype, ') ==='
    write(*,*) 'Integration_lambda: ', t_intlam_total, 's (', n_intlam_calls, ' calls)'
    write(*,*) 'Roots_finder:       ', t_roots_total,  's (', n_roots_calls,  ' calls)'
    write(*,*) 'WeightFunction:     ', t_weight_total, 's (', n_weight_calls, ' calls)'
    write(*,*) 'int_N_omega:        ', t_intn_total,   's (', n_intn_calls,   ' calls)'
    write(*,*) 'Total accounted:    ', &
        t_intlam_total+t_roots_total+t_weight_total+t_intn_total, 's'
#endif


    call MPI_FINALIZE(ierr)

end subroutine mpi_parallel