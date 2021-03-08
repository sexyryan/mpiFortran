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



!   REAL(KIND(0.D0)), DIMENSION(:,:,:), ALLOCATABLE :: intfac1,intfac2,intfac3
!   INTEGER*8 :: forward, inverse, forward_ddti
!   INTEGER :: mm, nn
!   REAL*8  :: normalize

!   Declare a type(C_PTR) :: p to hold the return value from FFTW’s allocation routine.
!   Set p = fftw_alloc_real(sz) for a real array, or p = fftw_alloc_complex(sz) for a complex array.
    TYPE(C_PTR) :: forward, inverse, forward_ddti
    TYPE(C_PTR) :: indata, indata_ddti

!   Because all sizes in the MPI FFTW interface are declared as ptrdiff_t in C, you should use integer(C_INTPTR_T) in Fortran (see FFTW Fortran type reference).
    INTEGER(C_INTPTR_T) :: mm, nn
    INTEGER(C_INTPTR_T) :: alloc_local, local_M, local_j_offset
    INTEGER(C_INTPTR_T) :: alloc_local_ddti, local_M_ddti, local_j_offset_ddti
!   In writing the transform pair, we have used the fact that the sine transform can be normalized so that it is identical to its inverse.
    REAL(C_DOUBLE) :: normalize

    INTEGER :: ierr, myid, nproc

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
    ALLOCATE( in(mm-1,nn-1), laminv(mm-1,nn-1) )
    ALLOCATE( lam(mm-1,nn-1) )
    ALLOCATE( lam1(mm-1,nn-1,mgridlev), &
              lam1i(mm-1,nn-1,mgridlev), &
              lam1inv(mm-1,nn-1,mgridlev) )
!    ALLOCATE( intfac1(mm-1,nn-1,mgridlev), intfac2(mm-1,nn-1,mgridlev), intfac3(mm-1,nn-1,mgridlev)  )
    ALLOCATE( in_ddti(mm,nn), laminv_ddti(mm,nn)  )


!   Initialize for MPI
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    CALL FFTW_MPI_INIT()



!   Get local data size and allocate (note dimensin reversal)
!   Because you need to call the ‘local_size’ function to find out how much space to allocate, and this may be larger than the local portion of the array (see MPI Data Distribution),
!   you should always allocate your arrays dynamically using FFTW’s allocation routines as described in Allocating aligned memory in Fortran.
!   (Coincidentally, this also provides the best performance by guaranteeding proper data alignment.)

    alloc_local = fftw_mpi_local_size_2d(mm-1, nn-1, MPI_COMM_WORLD, local_M, local_j_offset)
    forward = fftw_alloc_real(alloc_local)
!   Associate your pointer arr with the allocated memory p using the standard c_f_pointer subroutine: call c_f_pointer(p, arr, [...dimensions...]),
!   where [...dimensions...]) are an array of the dimensions of the array (in the usual Fortran order).
!   e.g. call c_f_pointer(p, arr, [L,M,N]) for an L × M × N array.
!   (Alternatively, you can omit the dimensions argument if you specified the shape explicitly when declaring arr.) You can now use arr as a usual multidimensional array.
    CALL c_f_pointer(indata, in, [nn-1,local_M])

    alloc_local_ddti = fftw_mpi_local_size_2d(mm, nn, MPI_COMM_WORLD, local_M_ddti, local_j_offset_ddti)
    forward_ddti = fftw_alloc_real(alloc_local_ddti)
    CALL c_f_pointer(indata_ddti, in_ddti, [nn,local_M_ddti])



!   Create MPI plan for in-place forward DFT (note dimension reversal)
!   For serial
!   FFTW_EXHAUSTIVE is like FFTW_PATIENT, but considers an even wider range of algorithms,
!   including many that we think are unlikely to be fast, to produce the most optimal plan but with a substantially increased planning time.
!   CALL dfftw_plan_r2r_2d(forward, mm-1,nn-1,in,in, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE)
!   CALL dfftw_plan_r2r_2d(forward_ddti, mm,nn,in_ddti,in_ddti, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE)

!   For MPI
!   Reference from 6.6 Other multi-dimensional Real-Data MPI Transforms in FFTW manual
!   To perform a two-dimensional L × M that is an REDFT10 (DCT-II) in the first dimension and an RODFT10 (DST-II) in the second dimension
!   plan = fftw_mpi_plan_r2r_2d(L, M, data, data, MPI_COMM_WORLD, FFTW_REDFT10, FFTW_RODFT10, FFTW_MEASURE);

    forward = fftw_mpi_plan_r2r_2d(mm-1, nn-1, in, in, MPI_COMM_WORLD, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE)
    forward_ddti = fftw_mpi_plan_r2r_2d(mm, nn, in_ddti, in_ddti, MPI_COMM_WORLD, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE)



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
  SUBROUTINE myfft_destroy_fft
    IMPLICIT NONE


    DEALLOCATE( lam,in, laminv, viscfac,vfac,con1,con2 )
!    deallocate( intfac1,intfac2,intfac3 )
    DEALLOCATE( lam1, lam1i, lam1inv )
    DEALLOCATE( in_ddti, laminv_ddti )

    CALL fftw_destroy_plan(forward)
    CALL fftw_destroy_plan(forward_ddti)
    CALL fftw_mpi_cleanup()
!   When you are done using the array, deallocate the memory by call fftw_free(p) on p.
    CALL fftw_free(indata)
    CALL fftw_free(indata_ddti)
    CALL mpi_finalize(ierr)


  END SUBROUTINE myfft_destroy_fft
! *****************************************************************************************
  FUNCTION dst( psi )
    IMPLICIT NONE
    REAL(KIND(0.D0)) :: psi(:,:)
    REAL(KIND(0.D0)) :: dst(2:mm,2:nn)

    ! discrete sine transform
    ! careful...two calls of dst need to be divided by "normalize"
    ! to return original vector
    !You must use the correct type of execute function, matching the way the plan was created.
    !Complex DFT plans should use fftw_execute_dft, Real-input (r2c) DFT plans should use use fftw_execute_dft_r2c, and real-output (c2r) DFT plans should use fftw_execute_dft_c2r.
    !The various r2r plans should use fftw_execute_r2r. Fortunately, if you use the wrong one you will get a compile-time type-mismatch error (unlike legacy Fortran).
    !However, note that in the MPI interface these functions are changed: fftw_execute_dft becomes fftw_mpi_execute_dft, etcetera. See Using MPI Plans.
    !void fftw_mpi_execute_r2r(fftw_plan p, double *in, double *out);
    in  = psi
    CALL fftw_mpi_execute_r2r(forward,in,in)
    dst = in

  END FUNCTION dst
! *****************************************************************************************
FUNCTION idst( psi )
    IMPLICIT NONE
    REAL(KIND(0.D0)) :: psi(:,:)
    REAL(KIND(0.D0)) :: idst(2:mm,2:nn)


    in  = psi / normalize
    CALL fftw_mpi_execute_r2r(forward,in,in)
    idst = in

END FUNCTION idst
! *****************************************************************************************
  FUNCTION ctci( omega )
    IMPLICIT NONE
    REAL(KIND(0.D0)) :: omega(:,:)
    REAL(KIND(0.D0)) :: ctci(2:mm,2:nn)

    in =  omega
    !CALL dfftw_execute_r2r(forward,in,in)
    in = dst(in)
    in = laminv * in
    !CALL dfftw_execute_r2r(forward,in,in)
    in = dst(in)
    ctci =  in

  END FUNCTION ctci
! *****************************************************************************************
  FUNCTION ainv( omega )
    IMPLICIT NONE
    REAL(KIND(0.D0)) :: omega(:,:)
    REAL(KIND(0.D0)) :: ainv(2:mm,2:nn)

    in = omega
    pre_ainv = dst (in)
    pre_ainv = pre_ainv * lam1i(:,:,1)
    ainv = idst (pre_ainv)

  END FUNCTION ainv
! *****************************************************************************************
  FUNCTION ddti( phi )
    IMPLICIT NONE
    REAL(KIND(0.D0)) :: phi(:,:)
    REAL(KIND(0.D0)) :: ddti(1:mm,1:nn)

    in_ddti = phi

    !CALL dfftw_execute_r2r(forward_ddti,in_ddti,in_ddti)
    in_ddti = dst(in_ddti)
    in_ddti = laminv_ddti * in_ddti
    !CALL dfftw_execute_r2r(forward_ddti,in_ddti,in_ddti)
    in_ddti = dst(in_ddti)

    ddti = in_ddti

  END FUNCTION ddti
! *****************************************************************************************

END MODULE myfft
