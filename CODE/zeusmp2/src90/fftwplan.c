#include <rfftw_mpi.h>
/*
                             Developed by
                 Laboratory of Computational Astrophysics
                University of Illinois at Urbana-Champaign
*/
/* By: PSLi, 12/30/99 */

void create_plan_(int *nx, int *ny, int *nz)
{ 
     FILE *output_file;
     rfftwnd_mpi_plan plan, iplan;

/*   Create plans */

     output_file=fopen("plan_wisdom","w");

     plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,*nz,*ny,*nx,
            FFTW_REAL_TO_COMPLEX,FFTW_MEASURE|FFTW_USE_WISDOM);

     fftw_export_wisdom_to_file(output_file);

     fclose(output_file);

     output_file=fopen("iplan_wisdom","w");
     
     iplan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,*nz,*ny,*nx,
             FFTW_COMPLEX_TO_REAL,FFTW_MEASURE|FFTW_USE_WISDOM);     

     fftw_export_wisdom_to_file(output_file);

     fclose(output_file);
}











