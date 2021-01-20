program example_mpi_complex_2d
   use, intrinsic :: iso_c_binding
   implicit none
   include 'mpif.h'
   include 'fftw3-mpi.f03'
!   include 'aslfftw3-mpi.f03'

   ! Parameter Definition
   integer(C_INTPTR_T), parameter :: NX = 4, NY = 3

   ! Variable Definition
   integer(C_INTPTR_T) :: i, ix, iy, alloc_local, local_n0, local_0_start

   complex(C_DOUBLE_COMPLEX), allocatable :: zin(:), zout(:), zin_all(:), zout_all(:)

   type(C_PTR) :: planf, planb

   integer :: my_rank, num_procs, ierr, local_n0_tmp

   integer, allocatable :: recv_counts(:), recv_offsets(:)

   ! MPI Preparation
   call MPI_Init(ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
   call fftw_mpi_init()

   ! Memory Allocation
   allocate(recv_counts(num_procs), recv_offsets(num_procs))
   alloc_local = fftw_mpi_local_size_2d(NY, NX, MPI_COMM_WORLD, local_n0, local_0_start)
   allocate(zin(alloc_local), zout(alloc_local), zin_all(NX * NY), zout_all(NX * NY))

   ! Plan Creation (out-of-place forward and backward FFT)
   planf = fftw_mpi_plan_dft_2d(NY, NX, zin, zout, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE)
   planb = fftw_mpi_plan_dft_2d(NY, NX, zout, zin, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE)
   if ((.not. c_associated(planf)) .or. (.not. c_associated(planb))) then
      write(*,*) "plan creation error!!"
      stop
   end if

   ! (just preparation for print of whole input/result data)
   local_n0_tmp = int(local_n0)
   call MPI_Allgather(local_n0_tmp, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
   do i = 1, num_procs
      recv_counts(i) = recv_counts(i) * int(NX)
   end do
   recv_offsets(1) = 0
   do i = 2, num_procs
      recv_offsets(i) = recv_offsets(i - 1) + recv_counts(i - 1)
   end do

   ! Input Initialization
   do iy = 0, local_n0 - 1
   do ix = 0, NX - 1
      zin(ix + NX * iy + 1) = dcmplx(dble(ix + (local_0_start + iy)), dble(ix * (local_0_start + iy)))
   end do
   end do

   ! (print of whole input data)
   call MPI_Allgatherv &
      & (zin, recv_counts(my_rank + 1), MPI_DOUBLE_COMPLEX, &
      & zin_all, recv_counts, recv_offsets, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
   if (my_rank == 0) then
      write(*,*) "===INPUT==="
      call print_complex(zin_all, NX, NY, .true., 1.0d0)
   end if

   ! FFT Execution (forward)
   call fftw_mpi_execute_dft(planf, zin, zout)

   ! (print of whole result data)
   call MPI_Allgatherv &
      & (zout, recv_counts(my_rank + 1), MPI_DOUBLE_COMPLEX, &
      & zout_all, recv_counts, recv_offsets, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
   if (my_rank == 0) then
      write(*,*) "===FORWARD==="
      call print_complex(zout_all, NX, NY, .false., 1.0d0)
   end if

   ! FFT Execution (backward)
   call fftw_mpi_execute_dft(planb, zout, zin)

   ! (print of whole result data)
   call MPI_Allgatherv &
      & (zin, recv_counts(my_rank + 1), MPI_DOUBLE_COMPLEX, &
      & zin_all, recv_counts, recv_offsets, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
   if (my_rank == 0) then
      write(*,*) "===BACKWARD==="
      call print_complex(zin_all, NX, NY, .true., dble(NX * NY))
   end if

   ! Plan Destruction
   call fftw_destroy_plan(planf)
   call fftw_destroy_plan(planb)

   ! Memory Deallocation
   deallocate(zin, zout, zin_all, zout_all, recv_counts, recv_offsets)

   call MPI_Finalize(ierr)
end program example_mpi_complex_2d

!
! === Print of Complex Data ===
!
subroutine print_complex(z, nx, ny, io, dscale)
   use, intrinsic :: iso_c_binding
   implicit none

   integer(C_INTPTR_T), intent(in) :: nx, ny
   complex(C_DOUBLE_COMPLEX), intent(in) :: z(nx * ny)
   logical, intent(in) :: io
   real(C_DOUBLE), intent(in) :: dscale

   character(27), parameter :: c1 = '############ zin ##########'
   character(27), parameter :: c2 = '########### zout ##########'
   character(2), parameter :: c3 = 'ix', c4 = 'iy'
   character(6), parameter :: c6 = '-real-', c7 = '-imag-'
   integer(C_INTPTR_T) :: ix, iy, i

   if (io) then
      write(*,'(A42)') c1
   else
      write(*,'(A42)') c2
   end if
   write(*,'(A5,A5,2X,A13,A13)') c3, c4, c6, c7
   do iy = 0, ny - 1
   do ix = 0, nx - 1
      i = ix + nx * iy + 1
      write(*,'(2I5,A6,E12.3,A,E12.3,A)') ix, iy, &
         & '(', dble(z(i)) / dscale, ',', dimag(z(i)) / dscale, ')'
   end do
   end do
end subroutine print_complex
