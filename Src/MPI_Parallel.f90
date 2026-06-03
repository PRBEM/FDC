subroutine mpi_parallel
  use constants
  IMPLICIT NONE

  include "mpif.h"

!  *********
  ! the size of the 2-D input data file
  integer::rows
  integer,parameter::columns=3      !!! L -- Energy -- equatorial pitch angle -- equatorial electron density  
  real*8::Ltemp,Etemp,PAtemp

  integer::local_rows,start_row,end_row,noutput_rows,noutput_columns,&
           noutput_rows_m,noutput_columns_m
  integer*8::row_number,row_remainder
 
  real*8,target::InputData(1200000,columns)  !1200000~89*211*61
  !real*8,allocatable::InputData(:,:)

  real*8,allocatable::LocalData(:,:)
  real*8,pointer::SendData(:,:)
  real*8,allocatable::OutputData(:,:),RecvData(:)
!  *********


!  *********     MPI stuff
  integer::ierr,numprocs,mype,fh,master
  integer::status(MPI_STATUS_SIZE)

  integer::i,j,k,rank
  INTEGER(KIND=MPI_OFFSET_KIND)::size

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  rows=NL*NE*NPA
!allocate(InputData(rows,columns))
  master = 0
  
!  ******** THIS CODE WILL USE MPI I/O EVENTUALLY BUT IS NOT WORKING YET **********
!!$  call MPI_FILE_OPEN(MPI_COMM_WORLD,'CRRES_data.txt',
!!$  $     MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
!!$  call MPI_FILE_GET_SIZE(fh,size,ierr)
!!$  write(*,*) rank,size
!!$  call MPI_FILE_SET_VIEW(fh, 0, MPI_REAL8, MPI_REAL, 'native',
!!$  $     MPI_INFO_NULL, ierr ) 
!!$  call MPI_FILE_READ_ALL(fh,testarr,10,MPI_REAL8,
!!$  $     status,ierr)
!!$  call MPI_FILE_READ(fh,testarr,10,MPI_REAL8,
!!$  $     status,ierr)
!!$  call MPI_FILE_CLOSE(fh,ierr)
!  ********************************************************************************


  if (mype.eq.master) then  ! START MASTER OPERATIONS

	do i=1,NL
	if(NL.eq.1)then
	Ltemp=L_min
	else 
	Ltemp=(L_max-L_min)/(NL-1.0)*(i-1.0)+L_min
	endif
	
        do j=1,NE
	if(NE.eq.1)then
	Etemp=E_min
	else 
	Etemp=10**((log10(E_max)-log10(E_min))/(NE-1.0)*(j-1.0)+log10(E_min))
	endif
	    do k=1,NPA
		PAtemp=(PA_max-PA_min)/(NPA-1.0)*(k-1.0)+PA_min
		InputData(((i-1)*NE+(j-1))*NPA+k,1)=Ltemp
		InputData(((i-1)*NE+(j-1))*NPA+k,2)=Etemp
		InputData(((i-1)*NE+(j-1))*NPA+k,3)=PAtemp
	    enddo
        enddo
     enddo
! *** Open output file and start writing first our own results ***********************
     open(unit=unit_output,file=Output_file,status='replace')

     write(unit_output,"('VARIABLES = ""L"", ""Energy"", ""Pitch-angle"", ""Daa Day^-1"", ""Dpp Day^-1"",""Dap Day^-1""')")
     write(unit_output,"(/)")
     write(unit_output,"('Zone T=""some name"", I=' I3 ',' 1X 'J=' I3, ',' 1X 'K=' I3)")NPA,NE,NL
     write(unit_output,"(/)")

!  ************************************************************

!  ***** Partition crres data and send to slaves
   row_remainder=mod(rows,numprocs-1)
   row_number=(rows-row_remainder)/(numprocs-1d0)
	noutput_rows=1
        noutput_columns=6
     allocate(RecvData(noutput_rows*noutput_columns))
   write(unit_log,"(I6 4x I6 4X I6 4x I6)")rows,numprocs-1,row_number,row_remainder
   do k=1,row_number
     do i=1,numprocs-1
	start_row=i+(k-1.0)*(numprocs-1d0)
	end_row=start_row
        SendData=>InputData(start_row:end_row,1:columns)

        call MPI_SEND(SendData,(end_row-start_row+1)*columns,MPI_DOUBLE_PRECISION,i,&
             110+i,MPI_COMM_WORLD,ierr)

     enddo
        
     do i=1,numprocs-1
        ! get partitioning information again and define noutput_rows and noutput_columns based on that

!************ !!!!!!! DEFINE noutput_rows and noutput_columns correctly here *****************
        ! noutput_rows would be local_rows*points_PA, I guess
        ! noutput_columns would be 4 judging by the write statement at l.242 in lstarcalculation
        
!************ !!!!!!! DEFINE noutput_rows and noutput_columns correctly here *****************
        call MPI_RECV(RecvData,noutput_columns*noutput_rows,MPI_DOUBLE_PRECISION,i, &
             110+i,MPI_COMM_WORLD, status, ierr)

        do j=1,noutput_rows
           write(unit_output,"(f6.3,5X,E15.6E4,5X,f5.2,5X,E15.6E4,5X,E15.6E4,5X,E15.6E4)") &
                RecvData(j+noutput_rows*0),&
                RecvData(j+noutput_rows*1),&
                RecvData(j+noutput_rows*2),&
                RecvData(j+noutput_rows*3),&
		RecvData(j+noutput_rows*4),&
                RecvData(j+noutput_rows*5)
        enddo

     enddo
   enddo

     do i=1,row_remainder
	start_row=i+row_number*(numprocs-1d0)
	end_row=start_row
        SendData=>InputData(start_row:end_row,1:columns)

        call MPI_SEND(SendData,(end_row-start_row+1)*columns,MPI_DOUBLE_PRECISION,i,&
             110+i,MPI_COMM_WORLD,ierr)

     enddo

     
     do i=1,row_remainder
        ! get partitioning information again and define noutput_rows and noutput_columns based on that

!************ !!!!!!! DEFINE noutput_rows and noutput_columns correctly here *****************
        ! noutput_rows would be local_rows*points_PA, I guess
        ! noutput_columns would be 4 judging by the write statement at l.242 in lstarcalculation

!************ !!!!!!! DEFINE noutput_rows and noutput_columns correctly here *****************
        call MPI_RECV(RecvData,noutput_columns*noutput_rows,MPI_DOUBLE_PRECISION,i, &
             110+i,MPI_COMM_WORLD, status, ierr)

        do j=1,noutput_rows
           write(unit_output,"(f6.3,5X,E15.6E4,5X,f5.2,5X,E15.6E4,5X,E15.6E4,5X,E15.6E4)") &
                RecvData(j+noutput_rows*0),&
                RecvData(j+noutput_rows*1),&
                RecvData(j+noutput_rows*2),&
                RecvData(j+noutput_rows*3),&
		RecvData(j+noutput_rows*4),&
                RecvData(j+noutput_rows*5)
        enddo

     enddo

     close(unit_output)
! END OF MASTER OPERATIONS
  else ! START SLAVE OPERATIONS
     ! Allocate and receive our local chunk of crres data     
!     call row_parti(rows,numprocs,mype,start_row,end_row)
 row_remainder=mod(rows,numprocs-1)
   row_number=(rows-row_remainder)/(numprocs-1d0)
        noutput_rows=1
        noutput_columns=6
    if(mype.le.row_remainder)row_number=row_number+1
        allocate(LocalData(1,columns))
        allocate(OutputData(noutput_rows,noutput_columns))
         do k=1,row_number
     call MPI_RECV(LocalData,columns,MPI_DOUBLE_PRECISION, master, &
          110+mype,MPI_COMM_WORLD, status, ierr)
!************ !!!!!!! DEFINE noutput_rows and noutput_columns correctly here *****************
        ! noutput_rows would be local_rows*points_PA, I guess
        ! noutput_columns would be 4 judging by the write statement at l.242 in lstarcalculation
!************ !!!!!!! DEFINE noutput_rows and noutput_columns correctly here *****************
    OutputData(1,4)=1
        OutputData(1,5)=2
        OutputData(1,6)=3
     call sum_n(LocalData,1,columns,OutputData,noutput_rows,noutput_columns)
! send results back to master
     
! f = fopen('output'+ slave_number);


!
!   OutputData(1,1)=LocalData(1,1)
 !       OutputData(1,2)=LocalData(1,2)
 !       OutputData(1,3)=LocalData(1,3)
 !   OutputData(1,4)=1
 !       OutputData(1,5)=2
 !       OutputData(1,6)=3
        call MPI_SEND(OutputData,noutput_rows*noutput_columns,MPI_DOUBLE_PRECISION,master,&
          110+mype,MPI_COMM_WORLD,ierr)
        enddo
  endif

! OK, WE ARE DONE
  call MPI_FINALIZE(ierr)
end subroutine mpi_parallel

