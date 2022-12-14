include 'mkl_vsl.f90'                     ! For pseudo-random number (PRN) generation with the MKL VSL library
module globals                            ! For MKL and MPI stuff and truly global variables
USE MKL_VSL_TYPE                          ! Routines for quasi-random number generation 
USE MKL_VSL                                
implicit none
include 'mpif.h'                          ! Fortran MPI header file
integer :: n                              ! Size of the matrices
integer, parameter :: EXIT_TAG=1,DEFAULT_TAG=0,void=0 ! Used to signal workers to either compute or exit
integer :: proc_num,num_procs             ! MPI and marix splitting parameters
double precision, parameter :: pi=2q0*asin(1q0)
end module globals

! This module now holds the PRN MKL stuff
module auxiliary
use globals
implicit none
TYPE (VSL_STREAM_STATE) :: stream         ! Variable that designates the PRN stream
contains
subroutine get_random_seed(seed)
integer :: un,seed,istat

! Read seed integer from /dev/urandom, a file filled with PRN generated from the state of the computer -
! pretty much the most random way to set a random seed on a computer
open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
   read(un) seed
close(un)

end subroutine get_random_seed

! Start the PRN stream - called every time a matrix is filled
subroutine open_random_stream(seed)
! Setup of stream of pseudo-random numbers using MKL library.
integer :: ierr,seed

! Use a Mersenne Twister pseudorandom number generator; the PRN sequence is identified by the variable "stream"
ierr=vslnewstream(stream,VSL_BRNG_MT19937,seed)
if(ierr.ne.0) then
   print *,'Trouble with the MKL RNG, returns flag ',ierr
end if

end subroutine open_random_stream

subroutine close_random_stream
! Close (de-allocate) pseudo-random number stream
integer :: rnd_ierr

rnd_ierr=vsldeletestream( stream )

end subroutine close_random_stream
end module auxiliary

! The only thing happening here is: initialize MPI, send the manager and workers to their respective routines
! and wrap up when done
program main
  use globals
  implicit none
  logical :: ok
  real(8) :: wtime
  
  call init_MPI(ok)                   ! Start MPI, all processes let all processes know if it worked.
  if(proc_num.eq.0) wtime=MPI_wtime() ! Start the clock using a function from the MPI library

  if(ok) then
     if(proc_num.eq.0) then
        call manager()
     else
        call worker()
     end if
  end if
! Print the wall time
if(proc_num.eq.0) then
   wtime=MPI_wtime()-wtime
   print *,'wtime=',wtime,'(s)'
end if

  ! Deallocate stuff
  call finalize()
end program main

subroutine init_MPI(ok)
    use globals
implicit none
! MPI initialization
logical :: my_ok=.true.,ok
integer :: ierr
! Initialialize MPI, obtain the number of procs and proc number, return error flag in case of errors
call mpi_init(ierr)
if(ierr.ne.0) then
   print *,'MPI_init failed!'
   my_ok=.false.
end if
if(my_ok) then
   call mpi_comm_size(MPI_COMM_WORLD,num_procs,ierr)
   if(ierr.ne.0) then
      print *,'MPI_comm_size failed!'
      my_ok=.false.
   end if
   if(my_ok) then
      call mpi_comm_rank(MPI_COMM_WORLD,proc_num,ierr)
      if(ierr.ne.0) then
         print *,'MPI_comm_rank failed!'
         my_ok=.false.
      end if
   end if
end if
! Check if everyone is ok.
call mpi_allreduce(my_ok,ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
end subroutine init_MPI

subroutine finalize
use globals
implicit none
logical :: init
integer :: ierr

! Wait until everybody is done. One process "finalize"ing while others are still working can cause ugly crashes
call mpi_barrier(MPI_COMM_WORLD,ierr)
! Check if MPI has been initialized - trying to "finalize" on a proc that has not "initialized" leads to ugly crashes
call mpi_initialized(init,ierr)
! If it is, call finalize.
if(init) call mpi_finalize(ierr)

end subroutine finalize
