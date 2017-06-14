c=======================================================================
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
      subroutine mstart
c
c  PURPOSE:  Starts a run.
c
c  EXTERNALS: SETUP, MGET, RESTART
c
c  LOCALS:
c
c  MODIFIED: 7/20/01 by PSLi
c  Modified  09/01/2006 by John Hayes; added "if(myid_w .eq. 0) condition
c  to CLOSE statement for unit 1 (zmp_inp).
c-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use root
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer  :: irestart
c
      integer incr,strtoi
c
      namelist /mpitop/ ntiles, periodic
      namelist /rescon/ irestart,tdump,dtdump,id,resfile
      namelist /iocon/ thdf,dthdf,thist,dthist,tusr,dtusr,t_out
     .                ,ttsl,dttsl
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c
c Open input and log files -- master thread only.
c
      if (myid_w .eq. 0) then
C        open(unit=1,file='zmp_inp',status='old')
        open(unit=2,file='zmp_log',status='unknown')
      endif
CJH
c------------------------  MPI TOPOLOGY  -------------------------------c
c  ntiles:   elements equal the number of tiles in each direction.
c  periodic: elements are true if grid is periodic in that direction --
c            we check hydro BC flags for defaults, but can override.
c
      if (myid_w .eq. 0) then
        ntiles(1) = 1
        ntiles(2) = 1
        ntiles(3) = 1
        periodic(1) = .false.
        periodic(2) = .false.
        periodic(3) = .false.
        read (1,mpitop)
        write(2,mpitop)
      endif
c
#ifdef MPI_USED
c
c Tell the others what the master has read.
c
      call MPI_BCAST(ntiles  , 3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(periodic, 3,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c
c Quit if the number of processors indicated on the command line
c differs from the number of tiles specified in the input file.
c
      if (nprocs_w .ne. ntiles(1)*ntiles(2)*ntiles(3)) then
        if (myid_w .eq. 0)
     1  write(*,"(/'MSTART: The number of threads ',i3,' does not match',
     2            /'MSTART: the number input via ntiles ',i3,' in',
     3            /'MSTART: input file zmp_inp; aborting the run...')")
     4    nprocs_w, ntiles(1)*ntiles(2)*ntiles(3)
        call MPI_FINALIZE( ierr )
        stop
      endif
c
c Create a virtual Cartesian topology for the domain decomposition.
c
      call MPI_CART_CREATE( MPI_COMM_WORLD, 3, ntiles, periodic
     &                    , reorder, comm3d, ierr )
      call MPI_COMM_RANK( comm3d, myid,     ierr )
      call MPI_COMM_SIZE( comm3d, nprocs,   ierr )
c
c Find the ranks of my neighbors; find my virtual Cartesian coords.
c
      call MPI_CART_SHIFT( comm3d, 0, 1, n1m, n1p, ierr )
      call MPI_CART_SHIFT( comm3d, 1, 1, n2m, n2p, ierr )
      call MPI_CART_SHIFT( comm3d, 2, 1, n3m, n3p, ierr )
c
      call MPI_CART_COORDS( comm3d, myid, 3, coords, ierr )
c
#endif /* MPI */
c
c------------------------  RESTART CONTROL  ----------------------------
c
c  irestart: set to one for calculation restarted from restart dump
c  tdump: time of last restart dump
c  dtdump: time between restart dumps
c  id: character*2 tag attended to filenames to identify run
c  resfile: name of restart file to restart from
c
      if (myid_w .eq. 0) then
        irestart = 0
        tdump = 0.0
        dtdump = 0.0
        id     = 'aa'
        resfile= 'resaa000000.000'
        read (1,rescon)
        write(2,rescon)
      endif
#ifdef MPI_USED
c
c Tell the others what the master has read (do tdump, dtdump later). 
c
      call MPI_BCAST(irestart, 1,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(id      , 2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(resfile ,15,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif /* MPI */
c
c Read the remaining namelists (except for iocon) in setup/restart.
c
      if (irestart .eq. 0) then
        call setup
        time = 0.0D0
      else
#ifdef MPI_USED
        incr = strtoi(resfile,13,15)
        write(resfile,"(a3,a2,3i2.2,'.',i3.3)") 'res',id,coords(1)
     &                                      ,coords(2),coords(3),incr
#endif
        call mget(resfile)
        nwarn = 0
        ifsen(1) = 0
        call restart
      endif
c
c------------------------  I/O CONTROL ---------------------------------
c
c  thdf: time of last HDF dump
c  dthdf: time between HDF dumps
c  ttsl: time of last tslice dump
c  dttsl: time between tslice dumps
c  thist: time of last history dump
c  dthist: time between history dumps
c  tusr: time of last user dump
c  dtusr: time between user dumps
c
      if (myid_w .eq. 0) then
        if (irestart .eq. 0) then
           thdf  = 0.0
          dthdf  = 0.0
           ttsl  = 0.0
          dttsl  = 0.0
           thist = 0.0
          dthist = 0.0
           tusr  = 0.0
          dtusr  = 0.0
          do incr=1,nbuff-8
            t_out(incr) = 0.0
          enddo ! incr
        endif
        read (1,iocon)
        write(2,iocon)
      endif
#ifdef MPI_USED
      if (myid_w .eq. 0) then
        buf_in(1) = thdf
        buf_in(2) = dthdf
        buf_in(3) = thist
        buf_in(4) = dthist
        buf_in(5) = tusr
        buf_in(6) = dtusr
        buf_in(7) = tdump
        buf_in(8) = dtdump
        buf_in(9) = 0.0
        buf_in(10) = 0.0
        buf_in(11) = ttsl
        buf_in(12) = dttsl
        do incr=13,nbuff
          buf_in(incr) = t_out(incr-12)
        enddo ! incr
      endif
c
c Tell the others what the master has read.  An array
c is used to pass just one message.
c
      call MPI_BCAST(buf_in,nbuff,MPI_FLOAT,0,MPI_COMM_WORLD
     &              ,ierr)
      if (myid_w .ne. 0) then
          thdf  = buf_in(1)
         dthdf  = buf_in(2)
          thist = buf_in(3)
         dthist = buf_in(4)
          tusr  = buf_in(5)
         dtusr  = buf_in(6)
         tdump  = buf_in(7)
        dtdump  = buf_in(8)
        ttsl   = buf_in(11)
        dttsl  = buf_in(12)
        do incr=1,nbuff-12
          t_out(incr) = buf_in(incr+12)
        enddo ! incr
      endif
#endif /* MPI */
c
c Output file names are of the form "hdfidccc.n" for ease of
c use with graphics packages that can process a series of files.
c
      if (irestart .eq. 0) then
        incr = 0
        ifsen(2) = 1
        ifsen(3) = 1
        ifsen(4) = 1
        ifsen(5) = 1
        ifsen(6) = 1

        write(tslfile,"(a3,i3.3,a2)") 'tsl',incr,id

        write(resfile,"(a3,a2,3i2.2,'.',i3.3)") 'res',id,coords(1)
     &                                      ,coords(2),coords(3),incr
        write(hdffile,"(a3,a2,3i2.2,'.',i3.3)") 'hdf',id,coords(1)
     &                                      ,coords(2),coords(3),incr
        write(hstfile,"(a3,a2,3i2.2,'.',i3.3)") 'hst',id,coords(1)
     &                                      ,coords(2),coords(3),incr
        write(usrfile,"(a3,a2,3i2.2,'.',i3.3)") 'usr',id,coords(1)
     &                                      ,coords(2),coords(3),incr

      else
        ifsen(2) = 0
        ifsen(3) = 0
        ifsen(4) = 0
        ifsen(5) = 0
        ifsen(6) = 0

        incr = strtoi(tslfile,4,6) + 1
        write(tslfile,"(a3,i3.3,a2)") 'tsl',incr,id 

        incr = strtoi(resfile,13,15) + 1
        write(resfile,"(a3,a2,3i2.2,'.',i3.3)") 'res',id,coords(1)
     &                                      ,coords(2),coords(3),incr
        incr = strtoi(hdffile,13,15) + 1
        write(hdffile,"(a3,a2,3i2.2,'.',i3.3)") 'hdf',id,coords(1)
     &                                      ,coords(2),coords(3),incr
        incr = strtoi(usrfile,13,15) + 1
        write(usrfile,"(a3,a2,3i2.2,'.',i3.3)") 'usr',id,coords(1)
     &                                      ,coords(2),coords(3),incr
      endif
      if(xtsl) then
       if (myid_w .eq. 0) then
         open (unit=31,file=tslfile,status='unknown')
       endif
      endif ! xtsl
      if (incr+1 .le. nbuff-10) then
c
c Adjust dump times to hit the next t_out, if defined non-zero.
c
        if (t_out(incr+1) .gt. 0.0) then
          tusr = t_out(incr+1) - dtusr
          thdf = t_out(incr+1) - dthdf
        endif
      endif
c
c  Close unit=1 (input deck).  Unit=2 (zeuslp) is kept open throughout
c  entire run to accept warning messages.  It is closed in zeusmp.
c
      if(myid_w .eq. 0) close(unit=1)
c
c Write out initial data dumps.
c
      if(irestart .eq. 0) then
c      call dataio( iswres, iswhdf, iswhst, iswusr
c    .            )
       call dataio( ifsen(2), ifsen(3), ifsen(4), ifsen(5), ifsen(6))
       nhy = 0
      else
       nhy = 0
      endif
c
      return
      end
