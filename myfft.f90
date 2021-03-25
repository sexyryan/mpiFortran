MODULE myfft

!   Because all sizes in the MPI FFTW interface are declared as ptrdiff_t in C,
!   you should use integer(C_INTPTR_T) in Fortran
    use, intrinsic :: iso_c_binding
    IMPLICIT NONE
!   INCLUDE 'fftw3.f'           ! needed for defining the plan
!   Instead of including fftw3.f03 as in Overview of Fortran interface,
!   you should include 'fftw3-mpi.f03' (after use, intrinsic :: iso_c_binding as before).
!   The fftw3-mpi.f03 file includes fftw3.f03, so you should not include them both yourself
    include 'mpif.h'
    include 'fftw3-mpi.f03'
    integer :: ierr, myid, nproc
!   REAL(KIND(0.D0)), ALLOCATABLE :: viscfac(:),vfac(:),con1(:),con2(:)
!   REAL(KIND(0.D0)), ALLOCATABLE :: lam(:,:), in(:,:), laminv(:,:), in_ddti(:,:), pre_ainv(:,:), laminv_ddti(:,:)
!   REAL(KIND(0.D0)), ALLOCATABLE :: lam1(:,:,:), lam1i(:,:,:), lam1inv(:,:,:)

!   You may be wondering if you need to search-and-replace real(kind(0.0d0)) (or whatever your favorite Fortran spelling of “double precision” is)
!   with real(C_DOUBLE) everywhere in your program, and similarly for complex and integer types.
!   The answer is no; you can still use your existing types.
    REAL(C_DOUBLE), ALLOCATABLE :: viscfac(:),vfac(:),con1(:),con2(:)
    REAL(C_DOUBLE), ALLOCATABLE :: lam(:,:), laminv(:,:), pre_ainv(:,:), laminv_ddti(:,:)
    REAL(C_DOUBLE), ALLOCATABLE :: lam1(:,:,:), lam1i(:,:,:), lam1inv(:,:,:)
!   Declare a pointer, arr, to your array of the desired type and dimensions.
!   For example, real(C_DOUBLE), pointer :: a(:,:) for a 2d real array, or complex(C_DOUBLE_COMPLEX), pointer :: a(:,:,:) for a 3d complex array.
    REAL(C_DOUBLE), POINTER :: in(:,:), in_ddti(:,:)
    REAL(C_DOUBLE), POINTER :: out(:,:), out_ddti(:,:)
    REAL(C_DOUBLE), POINTER :: in_all(:,:), in_all_ddti(:,:)
    REAL(C_DOUBLE), POINTER :: out_all(:,:), out_all_ddti(:,:)
!   REAL(KIND(0.D0)), DIMENSION(:,:,:), ALLOCATABLE :: intfac1,intfac2,intfac3
!   INTEGER*8 :: forward, inverse, forward_ddti
!   INTEGER :: mm, nn
!   REAL*8  :: normalize

!   Declare a type(C_PTR) :: p to hold the return value from FFTW’s allocation routine.
!   Set p = fftw_alloc_real(sz) for a real array, or p = fftw_alloc_complex(sz) for a complex array.
    TYPE(C_PTR) :: plan, plan_ddti ,inverse
    TYPE(C_PTR) :: cdata, cdata_ddti


!   Because all sizes in the MPI FFTW interface are declared as ptrdiff_t in C, you should use integer(C_INTPTR_T) in Fortran (see FFTW Fortran type reference).
    INTEGER(C_INTPTR_T) :: id, ix, iy
    INTEGER(C_INTPTR_T) :: mm, nn
    INTEGER(C_INTPTR_T) :: alloc_local, local_M, local_M_n, local_j_offset
    INTEGER(C_INTPTR_T) :: alloc_local_ddti, local_M_ddti,local_M_ddti_n, local_j_offset_ddti


!   In writing the transform pair, we have used the fact that the sine transform can be normalized so that it is identical to its inverse.
    REAL(C_DOUBLE) :: normalize


    integer,allocatable:: rcounts(:), rcounts_n(:) ! array of local_M's (for mpi_gatrherv)
    integer,allocatable:: displs(:), displs_n(:)  ! array of local_j_offset (for mpi_gatherv)

    integer,allocatable:: rcounts_ddti(:), rcounts_ddti_n(:) ! array of local_M's (for mpi_gatrherv)
    integer,allocatable:: displs_ddti(:), displs_ddti_n(:)  ! array of local_j_offset (for mpi_gatherv)
CONTAINS

! *****************************************************************************************
  SUBROUTINE myfft_setup_fft

    USE grid
    USE parameters
    IMPLICIT NONE

    INTEGER  :: i,j,k
    REAL(C_DOUBLE)   :: del2, del22, normalize_ddti
    REAL(C_DOUBLE)   :: lam_ddti(m,n)

    mm = m
    nn = n

    ALLOCATE( viscfac(mgridlev), con1(mgridlev), con2(mgridlev) )
    ALLOCATE( vfac(mgridlev) )
    ALLOCATE( laminv(mm-1,nn-1) )
    ALLOCATE( lam(mm-1,nn-1) )
    ALLOCATE( lam1(mm-1,nn-1,mgridlev), &
              lam1i(mm-1,nn-1,mgridlev), &
              lam1inv(mm-1,nn-1,mgridlev) )
!    ALLOCATE( intfac1(mm-1,nn-1,mgridlev), intfac2(mm-1,nn-1,mgridlev), intfac3(mm-1,nn-1,mgridlev)  )
    ALLOCATE( laminv_ddti(mm,nn)  )

    allocate(in_all(mm-1 ,nn-1), in_all_ddti(mm , nn))
    allocate(out_all(mm-1 ,nn-1), out_all_ddti(mm , nn))




    normalize_ddti = 4.d0*REAL( (mm+1)*(nn+1) )
    normalize      = 4.d0*REAL(     mm*nn     )

    ! eigenvalues for inverse of C^T C
    del2 = delta*delta
    DO k=1,mgridlev
       del22      =  del2*4.d0**(REAL(k)-1)
       viscfac(k) =  dt/Re/del22
       vfac(k)    =  0.5d0*dt/Re/del22/normalize
       con1(k)    =  1.5d0*dt/del22/normalize
       con2(k)    = -0.5d0*dt/del22/normalize
    ENDDO

    DO j=1,nn-1
       DO i=1,mm-1

          lam(i,j) = -2.d0*( COS( pi*REAL(i)/REAL(mm) ) + &
                             COS( pi*REAL(j)/REAL(nn) ) - 2.d0 )
          laminv(i,j) = 1.d0/lam(i,j)/normalize
          ! integrating factors for viscous terms
          DO k=1,mgridlev
             del22 = del2* 4.d0**(REAL(k)-1.D0)
             lam1(i,j,k)    =      (1.d0 - 0.5d0*dt*lam(i,j)/del22/Re)/normalize ! original
             lam1i(i,j,k)   = 1.d0/(1.d0 + 0.5d0*dt*lam(i,j)/del22/Re) ! original

! it appears that this HUGE array is NEVER used anywhere, so why compute/store it????
!             lam1inv(i,j,k) = 1.d0/(1.d0 - 0.5d0*dt*lam(i,j)/del22/Re)

!computed, but never USED !!????
!             intfac2(i,j,k) =  1.5*dt* EXP(      -dt*lam(i,j) / del22 / Re )/del22/normalize
!             intfac3(i,j,k) = -0.5*dt* EXP( -2.d0*dt*lam(i,j) / del22 / Re )/del22/normalize
!             intfac1(i,j,k) =          EXP(      -dt*lam(i,j) / del22 / Re )      /normalize
          END DO
      END DO
    END DO


    DO j=1,nn
       DO i=1,mm
          lam_ddti(i,j)    = 2.d0*( COS( pi*REAL(i)/REAL(mm+1) ) + &
                                    COS( pi*REAL(j)/REAL(nn+1) ) - 2.d0 )
          laminv_ddti(i,j) = 1.d0/lam_ddti(i,j)/normalize_ddti
       END DO
    END DO




  END SUBROUTINE myfft_setup_fft
! *****************************************************************************************


! *****************************************************************************************
  SUBROUTINE myfft_destroy_fft
    IMPLICIT NONE

    DEALLOCATE( lam, laminv, viscfac,vfac,con1,con2 )
!    deallocate( intfac1,intfac2,intfac3 )
    DEALLOCATE( lam1, lam1i, lam1inv )

    DEALLOCATE( laminv_ddti )

  END SUBROUTINE myfft_destroy_fft
! *****************************************************************************************
  function dst(psi)

    REAL(C_DOUBLE) :: psi(:,:)
    REAL(C_DOUBLE) :: dst(2:mm,2:nn)
    integer :: ix, iy, j


    alloc_local = fftw_mpi_local_size_2d(nn-1, mm-1,&
    & MPI_COMM_WORLD, local_M, local_j_offset)

    cdata = fftw_alloc_real(alloc_local)
    CALL c_f_pointer(cdata, in, [mm-1,local_M])

    plan = fftw_mpi_plan_r2r_2d &
    & (nn-1, mm-1, in, out, MPI_COMM_WORLD, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE)


    if ((.not. c_associated(plan))) then
       write(*,*) "plan creation error!!"
       stop
    end if

    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
    CALL MPI_ALLGATHER (local_M, 1, MPI_INTEGER, rcounts, 1, MPI_INTEGER,&
     & MPI_COMM_WORLD, IERR)

    displs(1) = 0
    do j=1,nproc
      if((j-1).ne.0)displs(j) = displs(j-1) + rcounts(j-1)
    enddo

    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

    !Input Initialization
   do ix = 1, local_M
     do iy = 1, mm-1
     in(iy,ix) = psi(iy, ix + local_j_offset)
     end do
   end do

   local_M = mm-1 * local_M
   rcounts = mm-1 * rcounts
   displs  = mm-1 * displs

   CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
   call fftw_mpi_execute_r2r(plan, in, out)
   CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

   call MPI_gatherv (out,local_M, MPI_REAL8,&
   & out_all, rcounts, displs,MPI_REAL8,&
   & 0, MPI_COMM_WORLD, ierr)

    dst = out_all

    call fftw_destroy_plan(plan)
    call fftw_mpi_cleanup()
    call fftw_free(cdata)

  END function dst
! *****************************************************************************************
!FUNCTION idst( psi )
!    IMPLICIT NONE
!    REAL(KIND(0.D0)) :: psi(:,:)
!    REAL(KIND(0.D0)) :: idst(2:mm,2:nn)


!    in  = psi / normalize
!    CALL fftw_mpi_execute_r2r(plan,in,in)
!    idst = in

!END FUNCTION idst
! *****************************************************************************************
function ctci( omega )
    IMPLICIT NONE
    REAL(C_DOUBLE) :: omega(:,:)
    REAL(C_DOUBLE) :: ctci(2:mm,2:nn)
    integer :: ix, iy, j

    alloc_local = fftw_mpi_local_size_2d(nn-1, mm-1,&
    & MPI_COMM_WORLD, local_M, local_j_offset)
    cdata = fftw_alloc_real(alloc_local)
    CALL c_f_pointer(cdata, in, [mm-1,local_M])

    plan = fftw_mpi_plan_r2r_2d &
    & (nn-1, mm-1, in, out, MPI_COMM_WORLD, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE)

    if ((.not. c_associated(plan))) then
       write(*,*) "plan creation error!!"
       stop
    end if

    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
    CALL MPI_ALLGATHER (local_M, 1, MPI_INTEGER, rcounts, 1, MPI_INTEGER,&
     & MPI_COMM_WORLD, IERR)

     displs(1) = 0
     do j=1,nproc
       if((j-1).ne.0)displs(j) = displs(j-1) + rcounts(j-1)
     enddo

     CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
    !   in =  omega
    do ix = 1, local_M
      do iy = 1, mm-1
      in(iy,ix) = omega(iy, ix + local_j_offset)
      end do
    end do

    local_M_n = mm-1 * local_M
    rcounts_n = mm-1 * rcounts
    displs_n  = mm-1 * displs

    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
    CALL fftw_mpi_execute_r2r(plan,in,out)
    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)


    call MPI_gatherv (out,local_M_n, MPI_REAL8,&
    & out_all, rcounts_n, displs_n, MPI_REAL8,&
    & 0, MPI_COMM_WORLD, ierr)

!    in = laminv * in
    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
    out_all = laminv * out_all
    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

    do ix = 1, local_M
      do iy = 1, mm-1
      in(iy,ix) = out_all(iy, ix + local_j_offset)
      end do
    end do

    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
    CALL fftw_mpi_execute_r2r(plan,in,out)
    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

    call MPI_gatherv (out,local_M_n, MPI_REAL8,&
    & out_all, rcounts_n, displs_n, MPI_REAL8,&
    & 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

    ctci =  out_all

    call fftw_destroy_plan(plan)
    call fftw_mpi_cleanup()
    call fftw_free(cdata)

  END function ctci
! *****************************************************************************************
  FUNCTION ainv( omega )
    IMPLICIT NONE
    REAL(C_DOUBLE) :: omega(:,:)
    REAL(C_DOUBLE) :: ainv(2:mm,2:nn)
    integer :: ix, iy, j

    alloc_local = fftw_mpi_local_size_2d(nn-1, mm-1,&
    & MPI_COMM_WORLD, local_M, local_j_offset)
    cdata = fftw_alloc_real(alloc_local)
    CALL c_f_pointer(cdata, in, [mm-1,local_M])


    plan = fftw_mpi_plan_r2r_2d &
    & (nn-1, mm-1, in, out, MPI_COMM_WORLD, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE)

    if ((.not. c_associated(plan))) then
       write(*,*) "plan creation error!!"
       stop
    end if

    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
    CALL MPI_ALLGATHER (local_M, 1, MPI_INTEGER, rcounts, 1, MPI_INTEGER,&
     & MPI_COMM_WORLD, IERR)

     displs(1) = 0
     do j=1,nproc
       if((j-1).ne.0)displs(j) = displs(j-1) + rcounts(j-1)
     enddo

     CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

!   in = omega
    do ix = 1, local_M
      do iy = 1, mm-1
      in(iy,ix) = omega(iy, ix + local_j_offset)
      end do
    end do

    local_M_n = mm-1 * local_M
    rcounts_n = mm-1 * rcounts
    displs_n  = mm-1 * displs

    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
    CALL fftw_mpi_execute_r2r(plan,in,out)
    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

    call MPI_gatherv (out,local_M_n, MPI_REAL8,&
    & out_all, rcounts_n, displs_n, MPI_REAL8,&
    & 0, MPI_COMM_WORLD, ierr)

    out_all = lam1i(:,:,1) * out_all / normalize
!    in = lam1i(:,:,1) * in / normalize

    do ix = 1, local_M
      do iy = 1, mm-1
      in(iy,ix) = out_all(iy, ix + local_j_offset)
      end do
    end do

    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
    CALL fftw_mpi_execute_r2r(plan,in,out)
    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

    call MPI_gatherv (out,local_M_n, MPI_REAL8,&
    & out_all, rcounts_n, displs_n,MPI_REAL8,&
    & 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

     ainv = out_all
!    in = omega
!    pre_ainv = dst (in)
!    pre_ainv = pre_ainv * lam1i(:,:,1)
!    ainv = idst (pre_ainv)

    call fftw_destroy_plan(plan)
    call fftw_mpi_cleanup()
    call fftw_free(cdata)

  END FUNCTION ainv
! *****************************************************************************************
function  ddti(phi)
    IMPLICIT NONE
    REAL(C_DOUBLE) :: phi(:,:)
    REAL(C_DOUBLE) :: ddti(1:mm,1:nn)
    integer :: ix, iy, j

    alloc_local_ddti = fftw_mpi_local_size_2d(nn, mm,&
    & MPI_COMM_WORLD, local_M_ddti, local_j_offset_ddti)
    cdata_ddti = fftw_alloc_real(alloc_local_ddti)
    CALL c_f_pointer(cdata_ddti, in_ddti, [mm,local_M_ddti])

    plan_ddti = fftw_mpi_plan_r2r_2d &
    & (nn, mm, in_ddti, out_ddti, MPI_COMM_WORLD, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE)

    if ((.not. c_associated(plan)) .or. (.not. c_associated(plan_ddti))) then
       write(*,*) "plan creation error!!"
       stop
    end if

    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
    CALL MPI_ALLGATHER (local_M_ddti, 1, MPI_INTEGER, rcounts_ddti, 1, MPI_INTEGER,&
     & MPI_COMM_WORLD, IERR)

     displs(1) = 0
     do j=1,nproc
       if((j-1).ne.0)displs(j) = displs_ddti(j-1) + rcounts_ddti(j-1)
     enddo
    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

    local_M_ddti_n = mm-1 * local_M_ddti
    rcounts_ddti_n = mm-1 * rcounts_ddti
    displs_ddti_n  = mm-1 * displs_ddti

!    in_ddti = phi
    do ix = 1, local_M_ddti
      do iy = 1, mm
      in_ddti(iy,ix) = out_all_ddti(iy, ix + local_j_offset_ddti)
      end do
    end do

    local_M_ddti_n = mm-1 * local_M_ddti
    rcounts_ddti_n = mm-1 * rcounts_ddti
    displs_ddti_n  = mm-1 * displs_ddti

    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
    call fftw_mpi_execute_r2r(plan_ddti, in_ddti, out_ddti)
    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

    call MPI_gatherv (out_ddti,local_M_ddti_n, MPI_REAL8,&
    & out_all_ddti, rcounts_ddti_n, displs_ddti_n, MPI_REAL8,&
    & 0, MPI_COMM_WORLD, ierr)

    out_all_ddti = laminv_ddti * out_all_ddti

    do ix = 1, local_M_ddti
      do iy = 1, mm
       in_ddti(iy,ix) = out_all_ddti(iy, ix + local_j_offset_ddti)
      end do
    end do

    CALL fftw_mpi_execute_r2r(plan_ddti,in_ddti,out_ddti)

    call MPI_gatherv (out_ddti,local_M_ddti_n, MPI_REAL8,&
    & out_all_ddti, rcounts_ddti_n, displs_ddti_n, MPI_REAL8,&
    & 0, MPI_COMM_WORLD, ierr)

    ddti = out_all_ddti

    call fftw_destroy_plan(plan_ddti)
    call fftw_mpi_cleanup()
    call fftw_free(cdata_ddti)

  END function ddti
! *****************************************************************************************

END MODULE myfft
