! Worker routine. Wait for a message (random seed) and, when received, create a random matrix,
! compute its largest eigenvalue and MPI_Send it back to the manager. Exit when the exit flag is received.
subroutine worker
use globals
use auxiliary
implicit none
logical :: conv
integer :: seed,ierr,recvd_tag,tag
integer, dimension(mpi_status_size) :: status
double precision :: eig
double precision, allocatable, dimension(:,:) :: A

! Receive matrix size - a single "one-to-all" communication.
call mpi_bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
allocate(A(n,n))

do while(.true.)
   call mpi_recv(seed,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,STATUS,ierr)
   ! MPI_Recv is blocking, so the worker waits here until a message arrives from the manager.
   recvd_tag = status (MPI_TAG) ! Extract the tag.
   if(recvd_tag.eq.EXIT_TAG)  then
      print *,'Slave ',proc_num,' exiting...'
      exit ! Exit from the "inifinite loop" when the manager sends the exit tag.
   end if
!  If the tag is not the exit tag, perform the computation.
   call make_random_matrix(A,seed)
   call compute_eig(A,eig,conv)
   if(.not.(conv)) print *,'No convergence in the eigenvalue computation!'
      ! Send result, along with a tag that indicates if it is correct.
   if(conv) then
      tag=0
   else
      tag=1
   end if
   call mpi_send(eig,1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierr) ! Send result and tag to manager.
end do

deallocate(A)

return
end subroutine worker

subroutine make_random_matrix(A,seed)
use globals
use auxiliary
implicit none
integer :: seed,ierr,i,j
double precision :: A(n,n)

call open_random_stream(seed) ! Open_random_stream is included in the main code - it sets up the PRN stream using the MKL.
do i=1,n
   do j=1,i
      ierr=vdRNGGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,stream,1,A(i,j),0d0,1d0)
      if(i.ne.j) then
         A(j,i) = A(i,j)
      end if
   end do
   A(i,i) = 2d0 * A(i,i)
end do
! Return a symmetric matrix with elements drawn from a standard normal distribution - except that the diagonal elements have standard deviation of 2.
call close_random_stream()

return
end subroutine make_random_matrix

subroutine compute_eig(A,eig,conv)
use globals
use auxiliary
implicit none
logical :: conv
integer :: ierr,M,lwork,info,i
double precision :: A(n,n),eig,eps,mu,wr(n),wi(n),B(n,n),vl(1,n),vr(1,n)
double precision, allocatable, dimension(:) :: work

conv=.true.

! Like most LAPACK routines, DGEEV requires an array for temporary use.
lwork=4*n
allocate(work(lwork))

! Pre-computing the spectrum
call DGEEV('N','N',n,A,n,wr,wi,vl,1,vr,1,work,lwork,info)
if(info.ne.0) then
   print *,'Error in DGEEV ',info
   conv=.false.
end if
eig=wr(maxloc(wr,1)) ! Pull out the largest eigenvalue - note that they are real-valued.

deallocate(work)

return
end subroutine compute_eig

