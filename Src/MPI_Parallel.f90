! ==============================================================================
! SUBROUTINE: mpi_parallel
! ==============================================================================
!
! PURPOSE:
!   Top-level MPI parallelisation routine. Constructs the full (L, E, alpha)
!   computation grid, distributes grid points across MPI processes in chunks,
!   collects results, and writes the output file.
!
!   Uses a master/slave pattern with chunked communication:
!     Master (rank 0):
!       1. Constructs the full input grid of (L, E, alpha) points
!       2. Opens the output file and writes the header
!       3. Sends each slave its entire chunk of grid points in one message
!       4. Receives each slave's results in one message and writes to file
!
!     Slaves (ranks 1 to numprocs-1):
!       1. Receive entire chunk of grid points
!       2. Loop over chunk, calling sum_n for each point
!       3. Send all results back to master in one message
!
!   Chunked communication reduces MPI messages from ~2*rows to
!   2*(numprocs-1), significantly reducing communication overhead.
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
!   Output row order is preserved: results are written in the same order
!   as the input grid (L outer loop, E middle loop, alpha inner loop).
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
!   - The grid is hardcoded to a maximum of 1,200,000 points.
!     If NL*NE*NPA exceeds this, the code will overflow.
!   - chunk_size and all array indices are plain integer (32-bit) throughout
!     to avoid MPI type mismatch issues.
!   - Load is distributed so that remainder rows go to the lowest-ranked
!     slaves (ranks 1 to row_remainder).
!   - Results are received and written in slave rank order (1, 2, ...,
!     numprocs-1), preserving the input grid row order.
!
! ==============================================================================

subroutine mpi_parallel
use constants
implicit none

    include "mpif.h"

    integer, parameter :: columns         = 3   ! L, Energy, pitch angle
    integer, parameter :: noutput_columns = 6   ! Daa, Dpp, Dap, L, E, alpha
    integer, parameter :: max_rows        = 1200000

    ! Grid dimensions and partitioning — all plain integer
    integer :: rows, row_number, row_remainder
    integer :: chunk_size, start_row

    ! Grid point temporaries
    real*8 :: Ltemp, Etemp, PAtemp

    ! Data arrays
    real*8 :: InputData(max_rows, columns)
    real*8, allocatable :: LocalData(:,:)
    real*8, allocatable :: OutputData(:,:)
    real*8, allocatable :: RecvData(:,:)

    ! MPI variables
    integer :: ierr, numprocs, mype, master
    integer :: status(MPI_STATUS_SIZE)
    integer :: i, j, k

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

    rows   = NL * NE * NPA
    master = 0

    ! Partition rows across slaves
    row_remainder = mod(rows, numprocs-1)
    row_number    = (rows - row_remainder) / (numprocs - 1)

    ! -------------------------------------------------------------------------
    ! MASTER: build grid, send chunks, receive and write results
    ! -------------------------------------------------------------------------
    if (mype == master) then

        ! Build input grid
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

        write(unit_log, "(I6 4x I6 4X I6 4x I6)") &
              rows, numprocs-1, row_number, row_remainder

        ! Send each slave its chunk in one message
        ! Use a contiguous temporary buffer to avoid memory layout mismatch
        ! between InputData(max_rows, columns) and LocalData(chunk_size, columns)
        start_row = 1
        do i = 1, numprocs-1
            if (i <= row_remainder) then
                chunk_size = row_number + 1
            else
                chunk_size = row_number
            end if

            block
                real*8 :: SendBuf(chunk_size, columns)
                integer :: m
                do m = 1, chunk_size
                    SendBuf(m,:) = InputData(start_row + m - 1, :)
                end do
                call MPI_SEND(SendBuf(1,1), chunk_size*columns, &
                              MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, ierr)
            end block

            start_row = start_row + chunk_size
        end do

        ! Receive results from each slave in order and write to file
        do i = 1, numprocs-1
            if (i <= row_remainder) then
                chunk_size = row_number + 1
            else
                chunk_size = row_number
            end if

            allocate(RecvData(chunk_size, noutput_columns))

            call MPI_RECV(RecvData(1,1), chunk_size*noutput_columns, &
                          MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)

            do j = 1, chunk_size
                write(unit_output, &
                      "(f6.3,5X,E15.6E4,5X,f5.2,5X,E15.6E4,5X,E15.6E4,5X,E15.6E4)") &
                      RecvData(j,1), RecvData(j,2), RecvData(j,3), &
                      RecvData(j,4), RecvData(j,5), RecvData(j,6)
            end do

            deallocate(RecvData)
        end do

        close(unit_output)

    ! -------------------------------------------------------------------------
    ! SLAVES: receive chunk, compute all points, send results back
    ! -------------------------------------------------------------------------
    else

        if (mype <= row_remainder) then
            chunk_size = row_number + 1
        else
            chunk_size = row_number
        end if

        allocate(LocalData(chunk_size, columns))
        allocate(OutputData(chunk_size, noutput_columns))

        call MPI_RECV(LocalData(1,1), chunk_size*columns, &
                      MPI_DOUBLE_PRECISION, master, mype, MPI_COMM_WORLD, status, ierr)

        ! Use explicit temporaries to avoid non-contiguous array slice issues
        block
            real*8 :: input_point(1, columns)
            real*8 :: output_point(1, noutput_columns)
            do k = 1, chunk_size
                input_point(1,:) = LocalData(k,:)
                call sum_n(input_point, 1, columns, &
                           output_point, 1, noutput_columns)
                OutputData(k,:) = output_point(1,:)
            end do
        end block

        call MPI_SEND(OutputData(1,1), chunk_size*noutput_columns, &
                      MPI_DOUBLE_PRECISION, master, mype, MPI_COMM_WORLD, ierr)

        deallocate(LocalData)
        deallocate(OutputData)

    end if

    call MPI_FINALIZE(ierr)

end subroutine mpi_parallel