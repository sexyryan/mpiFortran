# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <mpi.h>
# include <fftw3-mpi.h>

     int main(int argc, char **argv)
     {
         const ptrdiff_t L = 8192, M =8192;
         fftw_plan plan;
         fftw_complex *data ;
         ptrdiff_t alloc_local, local_L, local_L_start, i, j, ii;
         int proc_id, nproc;
         double xx, yy, rr, r2, t0, t1, t2, t3, tplan, texec,tmid;
         const double amp = 0.25;
         /* Initialize */
         MPI_Init(&argc, &argv);
         fftw_mpi_init();
         MPI_Comm_size(MPI_COMM_WORLD,&nproc);
         MPI_Comm_rank(MPI_COMM_WORLD,&proc_id);

         /* get local data size and allocate */
         alloc_local = fftw_mpi_local_size_2d(L, M, MPI_COMM_WORLD, &local_L, &local_L_start);
         data = fftw_alloc_complex(alloc_local);
         /* create plan for in-place forward DFT */
         t0 = MPI_Wtime();
         plan = fftw_mpi_plan_dft_2d(L, M, data, data, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);
         t1 = MPI_Wtime();
         /* initialize data to some function my_function(x,y) */
         for (i = 0; i < local_L; ++i) for (j = 0; j < M; ++j)
         {
	    ii = i + local_L_start;
            xx = ( (double) ii - (double)L/2 )/(double)L;
            yy = ( (double)j - (double)M/2 )/(double)M;
            r2 = pow(xx, 2) + pow(yy, 2);
            rr = sqrt(r2);
            if (rr <= amp)
            {
              data[i*M + j][0] = 1.;
              data[i*M + j][1] = 0.;
            }
            else
            {
              data[i*local_L + j][0] = 0.;
              data[i*local_L + j][1] = 1.;
            }
         }
         /* compute transforms, in-place, as many times as desired */
         t2 = MPI_Wtime();
         fftw_execute(plan);
         t3 = MPI_Wtime();
         /* Print results */
         tplan = t1 - t0;
         tmid = t2 - t1;
         texec = t3 - t2;
         if (proc_id == 0) printf(" T_plan = %f, T_mid = %f, T_exec = %f \n",tplan,tmid,texec);
         /* deallocate and destroy plans */
         fftw_destroy_plan(plan);
         fftw_mpi_cleanup();
         fftw_free ( data );
         MPI_Finalize();
     }
/*to compile mpicc -o ex2_c ex2_c.c -lfftw3_mpi -lfftw3 -lm*/
