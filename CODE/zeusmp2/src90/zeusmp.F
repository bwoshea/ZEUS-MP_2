c=======================================================================
c
c    \\\\\\\\\\        B E G I N   P R O G R A M          //////////
c    //////////               Z E U S M P                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
      program zeusmp
c
c PURPOSE
c   Main program for 3-D MPI version of ZEUS.
c
c AUTHOR
c   Robert A. Fiedler
c
c LAST MODIFIED by JCH, for F90
c   6/26/02.
c.......................................................................
c
c DECLARATIONS
c
      use real_prec
      use config
      use param
      use field
      use grid
      use root
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use gravmod
      use clockmod
      use impsoln
c
      implicit NONE
c
      real(rl4) :: cputime, wclock
c
      real(rl) :: zcs, etot, etot_glb
      real(rl) :: cpuall
c
      integer :: i , j , k, nwrite
      integer :: nx, ny, nz, snz, maxn
c
      ifsen(1) = 0
      ifsen(2) = 0
      ifsen(3) = 1
      ifsen(4) = 1
      ifsen(5) = 1
c
      myid_w    = 0
      myid      = 0
      nprocs_w  = 1
      nprocs    = 1
      coords(1) = 0
      coords(2) = 0
      coords(3) = 0
      reorder   = .false.
      totcgit = 0
      ncgcall = 0
c
c  call CONFIGURE
c
      call configure
c
c  Master writes greeting.
c
      if (myid_w .eq. 0) then
       call options
      endif
c
c Set up the problem: read input deck and possibly a restart file.
c
      call mstart
c
c Write out initial data dumps.
c
      call dataio( ifsen(2), ifsen(3), ifsen(4), ifsen(5), ifsen(6))
c
c Create best calculating plans for FFTW to perform forward
c and backward FFT.
c 
#ifdef FFT
      if(xgrvfft) then
        nx=ie-is+1
        ny=je-js+1
        nz=ke-ks+1
        call create_plan(ntiles(1)*nx,ntiles(2)*ny,ntiles(3)*nz)
      endif
#endif
c
c  Initialize cpu and wall clocks.  The quantities "cputime" and
c  "wclock" are the CPU and wall clock times (in seconds) in the main 
c  loop.
c
        wclock0 = 0.0
        cputime0 = 0.0
        call clocks (cputime, wclock)
        wclock0 = wclock
        cputime0 = cputime
c
      if (myid .eq. 0)
     &  write(6,"(/,' Set-up complete with ',i2,' warning(s):'
     &             ,' entering main loop...')") nwarn
c
      if(ldimen .eq. 1) nwrite = 1000
      if(ldimen .eq. 2) nwrite = 100
      if(ldimen .eq. 3) nwrite = 10
c
c--------------------------  start of main loop  -----------------------
c
c Execution ends when INTCHK returns a value of 1 for ifsen(1).
c
1000  continue

      nsub = 1

c
c
c Solve Poisson's equation for the gravitational potential.
c
      if(xgrav .or. xgrvfft) then
       call gravity
      endif ! xgrav
C
c
c Evaluate all non-advective terms in the evolution equations.
c
      if(lrad .eq. 0) then
       if(myid_w .eq. 0) then
        if(mod(nhy,nwrite) .eq. 0) then
         write(*,"('nhy, time, dt = ',i6,1p2d12.4)")
     .              nhy, time, dt
        endif ! mod
       endif ! myid_w
      endif ! lrad
c
      call srcstep
c
c Compute the advection of all field variables.
c
      if(ldimen .eq. 3) call transprt
      if(ldimen .eq. 2) call transprt_2D
      if(ldimen .eq. 1) call transprt_1D
c
c Update the step counter and current time value.
c
      nhy   = nhy   + 1
      time  = time  + dt
c
c Check the CPU time, number of steps, output times to determine if
c a stopping criterion has been met or output is desired.
c Also check for keyboard input, depending on the value of mbatch.
c
      call intchk( ifsen(2), ifsen(3), ifsen(4), ifsen(5), ifsen(6) )
c
c Compute new timestep
c
      call nudt
c
c Update the grid and related quantites.
c
      if(xvgrid) call newgrid
c
c Write out any desired output files now that everything has been
c updated.
c Skip dataio if the run is being terminated to avoid duplicate output.
c
      if (ifsen(1) .eq. 1) goto 2000
      call dataio( ifsen(2), ifsen(3), ifsen(4), ifsen(5), ifsen(6))
c
      goto 1000  !  Loop back to begin the next time step.
c
c--------------------------  end of main loop  -------------------------
c
c Terminate the run by making final dumps, write goodbyes
c
2000  continue
      call clocks (cputime, wclock)
#ifndef ARCH_CRAY
      tused = real(cputime)
#else
      tused = wclock
#endif
      ifsen(2) = 1
      ifsen(3) = 1
      ifsen(4) = 1
      ifsen(5) = 1
      ifsen(6) = 1
      call dataio( ifsen(2), ifsen(3), ifsen(4), ifsen(5) , ifsen(6))
c
#ifdef MPI_USED
c
c Sum everyone's cputime (stored in tused) to get CPU time for all 
c processes.
c
      call MPI_REDUCE(tused, cpuall, 1, MPI_FLOAT,
     &                MPI_SUM, 0, comm3d, ierr )
#else /* MPI */
       cpuall = tused
#endif /* MPI */
      if (myid .eq. 0) then      
c
c Let's assume tused is user + system time on the master thread.
c One would also like to keep track of wall-clock time and the sum
c of CPU times used by each processor.
c
        zcs = real(nprocs_w*nhy*nx1z*nx2z*nx3z)/(tused+tiny)
        write(6,"(/' Execution terminated with ',i4,' warning(s)')") 
     &     nwarn
        write(6,"(/' Performance summary:')")
        write(6,"('  zone-cycles per cpu second =',1pe12.5)") zcs
        write(6,"('  Master CPU                 =',1pe12.5, ' sec')") 
     &     tused
        write(6,"('  Average CPU/node           =',1pe12.5, ' sec')") 
     &     cpuall/dble(nprocs_w)
        write(6,"('  Wall Clock                 =',1pe12.5, ' sec')") 
     &     wclock
        write(6,"()")
c
        if(xgrav .and. (.not. xsphgrv)) then
         if(lgeom .ne. 1) then
          write(6,"(/' GRAVITY SUMMARY:')")
          write(6,"('  Percentage of cycles with Phi updates: ',i3,
     .              '%')")int(100*float(ncgcall)/float(nhy))
          write(6,"('  Average number of iterations/update  : ',
     .              1pd12.4)")float(totcgit)/float(ncgcall)
         endif ! lgeom
        endif ! xgrav
        if(lrad .ne. 0) then
         write(6,"(/' RADIATION SUMMARY:')")
         write(6,"('  Average number of N-R iterations/cycle: ',
     .              1pd12.4)")totnrit/float(nhy)
         write(6,"('  Average number of CG  iterations/cycle: ',
     .              1pd12.4)")totlsit/float(nhy)
        endif ! lrad
c
        close(unit=2)
C        close(unit=3)
        close(unit=30)
       if(xtsl) then
        close(unit=31)
       endif ! xtsl
      endif
c
c Turn off MPI
c
#ifdef MPI_USED
      call MPI_FINALIZE ( ierr )
#endif
c
c=======================================================================
c
c    \\\\\\\\\\          E N D  P R O G R A M             //////////
c    //////////               Z E U S M P                 \\\\\\\\\\
c
c=======================================================================
c
      end
