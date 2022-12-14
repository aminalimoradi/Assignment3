! Manager routine. Broadcast matrix size, then distribute work tickets (random seeds) to workers and recieve results
! until enough eigenvalues are received.
subroutine manager
use globals
use auxiliary
implicit none
logical :: start(1:num_procs-1)
integer :: seed,ierr,recvd=0,worker,exit_tags_sent,tag,ndat,p,failed
integer, dimension(mpi_status_size) :: status
double precision :: buf
double precision, allocatable, dimension(:) :: eigs

! Set matrix size and number of eigenvalues
call init(ndat)
! Broadcast matrix size
call mpi_bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! Allocate array for received eigenvalues
allocate(eigs(ndat))

start=.true.     ! Array of Booleans, each element is True if the corresponding worker is ready to start a new computation.
exit_tags_sent=0 ! Counter for the number exit tags sent.
failed=0         ! Counter for the number of failed eigenvalue computations.
open(22,file='eigs')

! Receive data until at least ndat have been received
do while(.true.)
   ! For each slave check if a seed is needed
   do p=1,num_procs-1
      if(start(p)) then
         call get_random_seed(seed) ! Obtains random seed from the MKL library (see main)
         if(recvd.lt.ndat-num_procs+2) then ! Assuming all workers that are still busy will return a result, this aims for ndat results.
            tag=DEFAULT_TAG
         else
            tag=EXIT_TAG                    ! When enough results are in (or pending), send the exit tag.
            exit_tags_sent=exit_tags_sent+1
         end if
         call mpi_send(seed,1,MPI_INTEGER,p,tag,MPI_COMM_WORLD,ierr)
         start(p)=.false.                   ! If worker k receives a work order, set its "start" flag to False.
      end if
   end do
   ! If all workers have been told to exit, the manager exits the loop as well.
   if(exit_tags_sent.eq.num_procs-1) then
      print *,'Master exiting...'
      exit
   end if

   ! After sending out a round of work orders, wait for results. Accept messages from any worker.
   call mpi_recv(buf,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
   worker = status(MPI_SOURCE) ! Extract the source - which worker sent that result?
   tag = status(MPI_TAG)       ! Extract the tag to find out if the computation worked out.
   if(tag.eq.0) then           ! Store the result if the computation worked out.
      recvd=recvd+1
      eigs(recvd)=buf
      write(22,'(i4,e14.7)') worker,buf
   else
      failed=failed+1
      print '(i5,a7,i5,a8)',recvd,' received, ',failed,' failed'
   end if
   start(worker)=.true.        ! Set the worker's "start" flag to True.
end do

close(22)

deallocate(eigs)
return
end subroutine manager

! Initialization - now just hard-coded matrix size and number of eigenvalues to compute. Perhaps change to reading from disk..
subroutine init(ndat)
use globals
implicit none
integer :: ndat

n=100
ndat=8

return
end subroutine init
