MODULE myfft
    use, intrinsic :: iso_c_binding
    IMPLICIT NONE
    include 'mpif.h'
    include 'fftw3-mpi.f03'

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
    TYPE(C_PTR) :: cdataf, cdataf_ddti
    TYPE(C_PTR) :: cdatab, cdatab_ddti

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
    integer :: ierr, myid, nprocs
contains

! *****************************************************************************************
  SUBROUTINE myfft_setup_fft

    USE grid
    USE parameters
    IMPLICIT NONE

    INTEGER  :: i,j,k,l
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
    allocate(rcounts(nprocs), displs(nprocs))
    allocate(rcounts_ddti(nprocs), displs_ddti(nprocs))

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


!****************************************************************

    displs(1) = 0
    displs_ddti(1) = 0

    alloc_local = fftw_mpi_local_size_2d(mm-1, nn-1, &
                & MPI_COMM_WORLD, local_M, local_j_offset)

    alloc_local_ddti = fftw_mpi_local_size_2d(mm, nn, &
                & MPI_COMM_WORLD, local_M_ddti, local_j_offset_ddti)

    cdataf = fftw_alloc_real(alloc_local)
    call c_f_pointer(cdataf, in, [mm-1, local_M])

    cdatab = fftw_alloc_real(alloc_local)
    call c_f_pointer(cdatab, out, [mm-1, local_M])

    cdataf_ddti = fftw_alloc_real(alloc_local_ddti)
    call c_f_pointer(cdataf_ddti, in_ddti, [mm, local_M_ddti])

    cdatab_ddti = fftw_alloc_real(alloc_local_ddti)
    call c_f_pointer(cdatab_ddti, out_ddti, [mm, local_M_ddti])
   ! Plan Creation
   plan = fftw_mpi_plan_r2r_2d &
      & (nn-1, mm-1, in, out, MPI_COMM_WORLD, FFTW_RODFT00, FFTW_REDFT00, FFTW_ESTIMATE)

   plan_ddti = fftw_mpi_plan_r2r_2d &
      & (nn, mm, in_ddti, out_ddti, MPI_COMM_WORLD, FFTW_RODFT00, FFTW_REDFT00, FFTW_ESTIMATE)

   if ((.not. c_associated(plan)) .or. (.not. c_associated(plan_ddti))) then
      write(*,*) "plan creation error!!"
      stop
   end if

   CALL MPI_ALLGATHER (local_M, 1, MPI_INTEGER, rcounts, 1, MPI_INTEGER,&
    & MPI_COMM_WORLD, IERR)

   CALL MPI_ALLGATHER (local_M_ddti, 1, MPI_INTEGER, rcounts_ddti, 1, MPI_INTEGER,&
    & MPI_COMM_WORLD, IERR)

    do k=1,nprocs
      if((k-1).ne.0)displs(k) = displs(k-1) + rcounts(k-1)
    enddo

    do k=1,nprocs
      if((k-1).ne.0)displs_ddti(k) = displs_ddti(k-1) + rcounts_ddti(k-1)
    enddo

    local_M_n = nn-1 * local_M
    rcounts_n = nn-1 * rcounts
    displs_n  = nn-1 * displs

    local_M_ddti_n = nn * local_M_ddti
    rcounts_ddti_n = nn * rcounts_ddti
    displs_ddti_n  = nn * displs_ddti



! ****************************************************************

  END SUBROUTINE myfft_setup_fft
! *****************************************************************************************


! *****************************************************************************************
  SUBROUTINE myfft_destroy_fft
    IMPLICIT NONE

    deallocate(in, out, in_all, out_all)
    deallocate(in_ddti, out_ddti, in_all_ddti, out_all_ddti)
    DEALLOCATE( lam, laminv, viscfac,vfac,con1,con2 )
!    deallocate( intfac1,intfac2,intfac3 )
    DEALLOCATE( lam1, lam1i, lam1inv )
    DEALLOCATE( laminv_ddti )

    call fftw_destroy_plan(plan)
    call fftw_destroy_plan(plan_ddti)

  END SUBROUTINE myfft_destroy_fft
! *****************************************************************************************
  function dst(psi)

    REAL(C_DOUBLE) :: psi(:,:)
    REAL(C_DOUBLE) :: dst(2:mm,2:nn)

    in_all = psi

    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
    CALL MPI_SCATTERV( in_all, rcounts_n, displs_n, MPI_REAL8, &
     & in, local_M_n, MPI_REAL8, &
     & 0, MPI_COMM_WORLD, IERR)

    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
    CALL fftw_mpi_execute_r2r(plan, in, out)
    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

    CALL MPI_gatherv (out,local_M_n, MPI_REAL8,&
    & out_all, rcounts_n, displs_n,MPI_REAL8,&
    & 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

    dst(2:mm,2:nn) = out_all(1:mm-1,1:nn-1)

  END function dst
! *****************************************************************************************
function dst_ddti(psi)

  REAL(C_DOUBLE) :: psi(:,:)
  REAL(C_DOUBLE) :: dst_ddti(1:mm,1:nn)

  in_all_ddti = psi

  CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
  CALL MPI_SCATTERV( in_all_ddti, rcounts_ddti_n, displs_ddti_n, MPI_REAL8, &
   & in_ddti, local_M_ddti_n, MPI_REAL8, &
   & 0, MPI_COMM_WORLD, IERR)

  CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
  CALL fftw_mpi_execute_r2r(plan_ddti, in_ddti, out_ddti)
  CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

  CALL MPI_gatherv (out_ddti,local_M_ddti_n, MPI_REAL8,&
  & out_all_ddti, rcounts_ddti_n, displs_ddti_n,MPI_REAL8,&
  & 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

  dst_ddti(1:mm,1:nn) = out_all_ddti(1:mm,1:nn)

END function dst_ddti
! *****************************************************************************************
function ctci( omega )
    IMPLICIT NONE
    REAL(C_DOUBLE) :: omega(:,:)
    REAL(C_DOUBLE) :: ctci(2:mm,2:nn)
    integer ::i, j

    in_all = omega
    in_all = dst( laminv * dst ( in_all ) )
    ctci = in_all

  END function ctci
! *****************************************************************************************
  FUNCTION ainv( omega )
    IMPLICIT NONE
    REAL(C_DOUBLE) :: omega(:,:)
    REAL(C_DOUBLE) :: ainv(2:mm,2:nn)

    in_all = omega

    in_all = dst ( ( dst(in_all) * lam1i(:,:,1) ) / normalize )

    ainv = in_all

  END FUNCTION ainv
! *****************************************************************************************
function  ddti(phi)
    IMPLICIT NONE
    REAL(C_DOUBLE) :: phi(:,:)
    REAL(C_DOUBLE) :: ddti(1:mm,1:nn)

    in_all_ddti = phi

    in_all_ddti = dst_ddti  ( dst_ddti(in_all_ddti) * laminv_ddti )

    ddti = in_all_ddti
  END function ddti
! *****************************************************************************************

END MODULE myfft
