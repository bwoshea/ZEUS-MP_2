c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              R A D S O L V E              \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c=======================================================================
      subroutine rad_solve
c
      use real_prec
      use config
      use root
      use field
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
c-----------------------------------------------------------------------
c     radiation solver menu:
c    
c     LRAD = 0:  No radiation transport
c     LRAD = 1:  Grey, flux-limited (non-equilibrium) diffusion
c     LRAD > 1:  User-supplied options
c-----------------------------------------------------------------------
c
      if(lrad .eq. 1) then
       call grey_fld(w3dd, e, w3de, w3dh, er)
       call rad_imp_dt(w3dh,er)
       dtnri2 = dtnri2 / (1.26*1.26)
      endif
      if(lrad .gt. 1) then
       if(myid .eq. 0) 
     .    write(*,"('RAD_SOLVE: Unsupported value of LRAD')")
#ifdef MPI_USED
       call mpi_finalize(ierr)
#endif
       stop
      endif
c
      return
      end
