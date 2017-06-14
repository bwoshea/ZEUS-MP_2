c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////             C O N F I G U R E           \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c     PURPOSE:  Reads the code configuration namelists that debuted
c               with ZEUS-MP Version 2.  These namelists serve the
c               functions formerly performed by CPP macros in the
c               discarded "zeusmp.def" file used to configure
c               ZEUS-MP Version 1.
c
c     Written by:  John Hayes; way back in '03.
c     Modified by: John Hayes; repeatedly since then.
c     Modified 09/01/2006: moved file OPEN for zmp_inp after the MPI
c     initialization and added "if(myid_w .eq. 0)" condition to the
c     OPEN statement.
c
c=======================================================================
c
      subroutine configure
c
      use real_prec
      use param
      use config
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer ::  confi_buf(10)
      logical ::  confl_buf(20)
      real(rl)::  confr_buf(2)
c
      namelist /geomconf/  lgeom, ldimen
      namelist /physconf/  lrad   , xhydro , xgrav, xmhd    , xgrvfft,
     .                     xptmass, xtotnrg, xiso , xvgrid  , xsubav,
     .                     xforce , xsphgrv, leos , nspec   , lopac
      namelist /ioconf/    xascii , xhdf, xrestart, xtsl
      namelist /preconf/   small_no, large_no
      namelist /arrayconf/ izones, jzones, kzones, maxijk
c
c----------------------------------------------------------------------
c     If parallel execution, start up MPI
c----------------------------------------------------------------------
c
#ifdef MPI_USED
       call MPI_INIT( ierr )
       call MPI_COMM_RANK( MPI_COMM_WORLD, myid_w  , ierr )
       call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs_w, ierr )
#else
       myid_w = 0
       myid   = 0
#endif /* MPI_USED */
c
c----------------------------------------------------------------------
c     Open zmp_conf run configuration file
c----------------------------------------------------------------------
c
      if(myid_w .eq. 0) open(unit=1,file='zmp_inp',status='old')
c
c----------------------------------------------------------------------
c     Initialize all parameters to default values before continuing
c     with configuration READ
c----------------------------------------------------------------------
c
      lgeom    = 1
      ldimen   = 3
      lrad     = 0
      leos     = 1
      lopac    = 0
      nspec    = 1
      izones   = 32
      jzones   = 32
      kzones   = 32
      maxijk   = 32
      xvgrid   = .false.
      xhydro   = .false.
      xforce   = .false.
      xgrav    = .false.
      xgrvfft  = .false.
      xsphgrv  = .false.
      xptmass  = .false.
      xmhd     = .false.
      xtotnrg  = .false.
      xiso     = .false.
      xascii   = .false.
      xhdf    = .false.
      xrestart = .false.
      xtsl     = .false.
      large_no = 1.0D99
      small_no = 1.0D-99
c
c----------------------------------------------------------------------
c     Read remaining namelists
c----------------------------------------------------------------------
c
      if(myid_w .eq. 0) then
       read(1,geomconf)
       confi_buf(1) = lgeom
       confi_buf(2) = ldimen
c
       read(1,physconf)
       confi_buf(3) = lrad
       confi_buf(4) = lopac
       confi_buf(5) = leos
       confi_buf(6) = nspec
       confl_buf(1) = xhydro
       confl_buf(2) = xforce
       confl_buf(3) = xmhd
       confl_buf(4) = xgrav
       confl_buf(5) = xgrvfft
       confl_buf(6) = xptmass
       confl_buf(7) = xtotnrg
       confl_buf(8) = xiso
       confl_buf(9) = xvgrid
       confl_buf(10) = xsubav
       confl_buf(11) = xsphgrv
c
       read(1,ioconf)
       confl_buf(12) = xascii
       confl_buf(13) = xhdf
       confl_buf(14) = xrestart
       confl_buf(15) = xtsl
c
c
       read(1,preconf)
       confr_buf(1) = small_no
       confr_buf(2) = large_no
c
       read(1,arrayconf)
       confi_buf( 7) = izones
       confi_buf( 8) = jzones
       confi_buf( 9) = kzones
       confi_buf(10) = maxijk
      endif ! myid_w
c
#ifdef MPI_USED
       call MPI_BCAST(confi_buf,10,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(confl_buf,15,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(confr_buf,2,MPI_DOUBLE_PRECISION,0,
     .                            MPI_COMM_WORLD,ierr)
       if(myid_w .ne. 0) then
        lgeom  = confi_buf(1)
        ldimen = confi_buf(2)
        lrad   = confi_buf(3)
        lopac  = confi_buf(4)
        leos   = confi_buf(5)
        nspec  = confi_buf(6)
        izones = confi_buf(7)
        jzones = confi_buf(8)
        kzones = confi_buf(9)
        maxijk = confi_buf(10)
        xhydro  = confl_buf(1)
        xforce  = confl_buf(2)
        xmhd    = confl_buf(3)
        xgrav   = confl_buf(4)
        xgrvfft = confl_buf(5)
        xptmass = confl_buf(6)
        xtotnrg = confl_buf(7)
        xiso    = confl_buf(8)
        xvgrid  = confl_buf(9)
        xsubav  = confl_buf(10)
        xsphgrv = confl_buf(11)
        xascii   = confl_buf(12)
        xhdf     = confl_buf(13)
        xrestart = confl_buf(14)
        xtsl     = confl_buf(15)
        small_no = confr_buf(1)
        large_no = confr_buf(2)
       endif ! myid_w
#endif /* MPI_USED */
c
#ifdef MPI_USED 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
#if defined USE_HDF4 && defined USE_HDF5
c
c-----------------------------------------------------------------------
c     Check to see that only one HDF viz file format is enabled
c     (HDF4/HDF5)
c-----------------------------------------------------------------------
c
      if(myid_w.eq. 0) then
       write(*,"(/'Build Error: you may not select both HDF4 and')")
       write(*,"('HDF5 formats for visulization data output!')")
       write(*,"(/'The format is selected in the Makefile by the')")
       write(*,"( 'numeric value of N in the -DUSE_HDFN entry in')")
       write(*,"( 'the ZMP_CPP definition line:')")
       write(*,"(/'   N = 4: HDF4 format')")
       write(*,"( '   N = 5: HDF5 format')")
       write(*,"(/'ZEUS-MP terminating...')")
      endif
#ifdef MPI_USED 
      call mpi_finalize(ierr)
#endif
      stop
#endif /* USE_HDF4 and USE_HDF5 */
#ifndef USE_HDF4 
#ifndef USE_HDF5
c
c-----------------------------------------------------------------------
c     if "xhdf = .true.", make sure that at least one HDF CPP macro
c     was defined in the Makefile
c-----------------------------------------------------------------------
c
      if(xhdf) then
       if(myid_w .eq. 0) then
        write(*,"(/'Build Error: since XHDF=.TRUE. in zmp_inp,')")
        write(*,"( 'you must also include -DUSE_HDFN in the')")
        write(*,"( 'ZMP_CPP definition line in the Makefile, where:')")
        write(*,"(/'   N = 4 for HDF4 format')")
        write(*,"( '       ** OR **         ')")
        write(*,"( '   N = 5 for HDF5 format')")
        write(*,"(/'NOTE: you may choose one format or the other,')")
        write(*,"( '      but NOT BOTH!')")
        write(*,"(/'ZEUS-MP terminating...')")
       endif ! myid
#ifdef MPI_USED 
       call mpi_finalize(ierr)
#endif
       stop
      endif ! xhdf
#endif /* NOT USE_HDF5 */
#endif /* NOT USE_HDF4 */
c
      call alloc_arrays
c
      return
      end
