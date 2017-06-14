/*
                             Developed by
                 Laboratory of Computational Astrophysics
                University of Illinois at Urbana-Champaign
*/
/* By: P.S. Li, 12/30/99 */

#include <stdio.h>
#include <math.h>
#include <rfftw_mpi.h>

const double pi2=2.*3.14159265358979324;

void poisson_solver_(int *myid_w, int *nprocs_w, int *nx, int *ny,
                     int *nz, double data[*nz][*ny][*nx])
{ 
     int i,j,k,jindex,lny;
     int local_nz,local_z_start,local_ny_after_transpose,
         local_y_start_after_transpose,total_local_size;
     fftw_real *work;
     fftw_real *local_data;
     fftw_real coeff;
     fftw_complex *cdata;
     rfftwnd_mpi_plan plan, iplan;
     FILE *input_file;

/*   Scaling factor for FFTW algorithm and negative sign for ZeusMP */

     double scale=-1./((*nx)*(*ny)*(*nz)*(*nprocs_w));

/*   Create plans */

     input_file=fopen("plan_wisdom","r");

     plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,*nz*(*nprocs_w),
            *ny,*nx,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE|FFTW_USE_WISDOM);

     fclose(input_file);

     input_file=fopen("iplan_wisdom","r");

     iplan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,*nz*(*nprocs_w),
             *ny,*nx,FFTW_COMPLEX_TO_REAL,FFTW_MEASURE|FFTW_USE_WISDOM);

     fclose(input_file);

/*   Dim.'s of REAL mydata */

     rfftwnd_mpi_local_sizes(plan,&local_nz,&local_z_start,
			     &local_ny_after_transpose,
			     &local_y_start_after_transpose,
			     &total_local_size);

     work = (fftw_real*)malloc(total_local_size*sizeof(fftw_real));

     local_data = (fftw_real*)malloc(total_local_size*sizeof(fftw_real));

     for(i=0; i < local_nz; i++)
        {
        for(j=0; j < *ny; j++)
          {
           for(k=0; k < *nx; k++)
            {
	      local_data[(i*(*ny)+j)*(2*(*nx/2+1))+k] =
		   (fftw_real) data[i][j][k];
	    }

	  }
	}

/*   Compute the forward transform */

     rfftwnd_mpi(plan,1,local_data,work,FFTW_TRANSPOSED_ORDER);

     cdata=(fftw_complex*) local_data;

/*   Mulitply tranformed coeff. */

     lny=*ny/(*nprocs_w);

     for(j=0; j < local_ny_after_transpose; j++)
       {
	 for(i=0; i < *nz*(*nprocs_w); i++)
	   {
	     for(k=0; k < (*nx/2+1); k++)
	       {
		 jindex=j+*myid_w*lny;
		 if (*myid_w == 0 && i == 0 && jindex == 0 && k == 0)
		   coeff=-1.;
		 else
		   coeff=1./(cos(pi2*k/(*nx))+cos(pi2*jindex/
			 (*ny))+cos(pi2*i/(*nz*(*nprocs_w)))-3.);
		 cdata[(j*(*nz*(*nprocs_w))+i)*(*nx/2+1)+k].im=
		   coeff*cdata[(j*(*nz*(*nprocs_w))+i)*
			      (*nx/2+1)+k].im;
		 cdata[(j*(*nz*(*nprocs_w))+i)*(*nx/2+1)+k].re=
		   coeff*cdata[(j*(*nz*(*nprocs_w))+i)*
			      (*nx/2+1)+k].re;
	       }
	   }
       }

/*   Compute inverse transform for potential */

     rfftwnd_mpi(iplan,1,local_data,work,FFTW_TRANSPOSED_ORDER);

/*   Scaling and send data back */

     for(i=0; i < local_nz; i++)
       {
	 for(j=0; j < *ny; j++)
	   {
	     for(k=0; k < *nx; k++)
	       {
		 data[i][j][k] = (fftw_real) local_data[(i*(*ny)+j)*
				 (2*(*nx/2+1))+k]*scale;
	       }

           }
       }

     free(work);
     free(local_data);
     rfftwnd_mpi_destroy_plan(plan);
     rfftwnd_mpi_destroy_plan(iplan);
     
}





