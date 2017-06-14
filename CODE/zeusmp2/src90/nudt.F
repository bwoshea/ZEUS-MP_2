c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                  N U D T                  \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine nudt
c
c    mln:zeus3d.nudt <------------------------ mhd time step calculation
c                                                              may, 1986
c
c    written by: Mike Norman
c    modified 1: May, 1986 by David Clarke; adapted for mhd.
c    modified 2: October, 1987 by Mike Norman; reworked for covariant
c                formulation.
c    modified 3: June, 1988 by Jim Stone; incorporated into ZEUS2D.
c    modified 4: February, 1990 by David Clarke; incorporated into
c                ZEUS3D.
c    modified 5: September, 1990 by David Clarke; moved magnetic fields
c                to face-centres.
c    modified 6: March 5, 1996 by Robert Fiedler; completely rewritten
c                for ZEUS-MP.
c    modified 7: Aug. 23, 1996 by Robert Fiedler; minor change to 
c                substep counter.
c    modified 8: Dec. 19, 1996 by Robert Fiedler; added radiation 
c                diffusion time step.
c    modified 9: Feb. 14, 1997 by RAF; corrected for non-hydro and
c                lower dimensionality options.
c    modified 10: kluge for Garching T3E (search on M-MML) 26 Feb 98
c    modified 11: December, 1999 by PSLi; add in dt calculation in case
c                 of subcycle of artificial viscosity.
c
c  PURPOSE:  Computes the new timestep for explicit calculations from
c  the values of the field variables updated by the source and transport
c  steps.
c
c  In explicit calculations, the timestep is given by:
c
c     dt = courno * sqrt [ min ( dtcs**2 + dtv1**2 + dtv2**2 + dtv3**2
c                              + dtal**2 + dtqq**2 + dtrd**2 ) ]
c
c  where the variable names are described below.  The timestep can be
c  reduced in size by any amount, but can be larger than the old timstep
c  by no more than a factor of 1.26.
c
c  LOCAL VARIABLES:
c
c  i-sweep
c  dr*i      inverse of distance across zone in 1-, 2-, or 3-direction
c  drimax    maximum of dr1i, dr2i, and dr3i
c  dt**i2i   square of the inverse time step of the ** physical process
c            gathered during the i-sweep.  Possible values of ** are:
c                cs = sound speed
c                v1 = fluid motion in x1 direction
c                v2 = fluid motion in x2 direction
c                v3 = fluid motion in x2 direction
c                al = Alfven speed
c                qq = artificial viscosity
c                rd = radiation diffusion
c            The first five are vectors; the next to last is a scalar
c            which has been computed in ARTIFICIALVISC (passed in root).
c            Last is a scalar computed in rad_pfl (passed in radexp).
c  dttoi2i   sum of dt**i2i (without artificial viscosity contribution)
c
c  j-sweep
c  dttoi2j   vector of maximum values of dttoi2i from each i-sweep.
c            This vector is filled during a j-sweep.
c  imaxj     i-index of each dttoi2j
c  dt**i2j   values of above at zone i=imaxj
c
c  k-sweep
c  dttoi2k   vector of maximum values of dttoi2j from each j-sweep.
c  imaxk     i-index of each dttoi2k.
c  jmaxk     j-index of each dttoi2k.
c  dt**i2j   values of above at zone i=imaxk, j=jmaxk
c
c  grand maximum inverse time
c  imin      i-index where dttoi2k is a maximum
c  jmin      j-index where dttoi2k is a maximum
c  kmin      k-index where dttoi2k is a maximum
c  dt**      time step of the ** physical process at i=imax, j=jmax,
c            and k=kmax
c  dtnew     new timestep, which is limited to be no greater than 1.26
c            times the previous timestep.
c
c  EXTERNALS:
c
c-----------------------------------------------------------------------
c
      use real_prec
      use config
      use param
      use field
      use bndry
      use root
      use scratch
      use grid
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer  :: i, j, k,
     .            imax, jmax, kmax,
     .            imin, jmin, kmin,
     .            k0  , k1  , k2  , j0
c
      real(rl) :: dtcsm, dtv1m, dtv2m, dtv3m, dtalm, dttoi2m
c
c-----------------------------------------------------------------------
c
c Find the minimum time step required by the Courant condition for
c this tile.  We will first compute 1/dt**2 required by each of the
c various physical processes in the calculation, and save the maximum
c value of their sum at each zone.  In the process, we will compute
c the updated velocity.
c
       dttoi2m = 0.0
       dtcsm   = 0.0
       dtv1m   = 0.0
       dtv2m   = 0.0
       dtv3m   = 0.0
       dtalm   = 0.0
       imin    = is
       jmin    = js
       kmin    = ks
       if (nx2z .eq. 1) then
c
c For one-dimensional problems only
c
         j0 = js
       else
         j0 = js + 1
       endif
       if (nx3z .eq. 1) then
c
c For one-dimensional or two-dimensional problems only
c
         k0 = ks
       else
         k0 = ks + 1
       endif
c
c Divide the computational volume into three equal pieces.  We must
c have at least 3 active zones in the 3-direction.
c
      if(ldimen .eq. 3) then
       k1 = int( real( ke - ks + 1 ) / 3.0 ) + ks
       k2 = int( real( ke - ks + 1 ) / 3.0 ) + k1
      else
       k1 = ks
       k2 = ks
      endif
c
      do i = 1, 6
       bvstat(i,1) = 0
       bvstat(i,3) = 0
       bvstat(i,4) = 0
       bvstat(i,5) = 0
       bvstat(i,7) = 0
       bvstat(i,8) = 0
      enddo
c......................................................................
c
c i boundaries
c
c    1) Post sends and receives. 
c
       nreq = 0
       nsub = nsub + 1
      if(xtotnrg) then
       imax = ie - 1
       if (nx2z .eq. 1) then
c
c One-dimensional problems only
c
         jmax = js
       else
         jmax = je - 1
       endif
       kmax = ke - 1
       call bvald  (1,1,0,0,0,0,d   )
       call bvalv1 (0,1,0,0,0,0,w3da)
      else ! xtotnrg
       imax = ie
CJH       jmax = je
CJH       kmax = ke
       if(ldimen .gt. 1) then
        jmax = je
       else
        jmax = js
       endif
       if(ldimen .gt. 2) then
        kmax = ke
       else
        kmax = ks
       endif
       call bvald  (1,0,0,0,0,0,d   )
      endif ! xtotnrg
c
c    2) Do first 1/3 of the interior points.
c
       call newdt (is+1,imax,j0  ,jmax,k0  ,k1
     &                ,  imin,jmin,kmin
     &                ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
c
c       subroutine newdt (ibeg,iend,jbeg,jend,kbeg,kend
c     &                ,  imin,jmin,kmin
c     &                ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
c


c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
       if(nreq .ne. 0) 
     .    call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
c j boundaries
c
c    1) Post sends and receives.
c
      if(ldimen .gt. 1) then
       nreq = 0
       nsub = nsub + 1
       if(xtotnrg) then
        call bvald  (0,0,1,1,0,0,d   )
        call bvalv2 (0,0,0,1,0,0,w3db)
       else ! xtotnrg
        call bvald  (0,0,1,0,0,0,d   )
       endif ! xtotnrg
      endif ! ldimen

c
c    2) Do middle 1/3 of the interior points, plus some on borders.
c
       call newdt (is  ,is  ,j0  ,jmax,k0  ,k1
     &                ,  imin,jmin,kmin
     &                ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
      if(xtotnrg) then
       call newdt (ie  ,ie  ,j0  ,jmax,k0  ,k1
     &                ,  imin,jmin,kmin
     &                ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
      endif ! xtotnrg
      if (ldimen .eq. 3) then
         call newdt (is  ,ie  ,j0  ,jmax,k1+1,k2
     &                  ,  imin,jmin,kmin
     &                  ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
      endif


c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
       if(nreq .ne. 0) 
     .    call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
c k boundaries
c
c    1) Post sends and receives.
c
      if(ldimen .eq. 3) then
       nreq = 0
       nsub = nsub + 1
       if(xtotnrg) then
        call bvald  (0,0,0,0,1,1,d   )
        call bvalv3 (0,0,0,0,0,1,w3dc)
       else ! xtotnrg
        call bvald  (0,0,0,0,1,0,d   )
       endif ! xtotnrg
      endif ! ldimen

c
c    2) Do last 1/3 of the interior points, plus some on borders.
c
CJH       if (nx2z .gt. 1) then
       if (ldimen .gt. 1) then
         call newdt (is  ,ie  ,js  ,js  ,k0  ,k2
     &                  ,  imin,jmin,kmin
     &                  ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
C        if(ldimen .eq. 2) then
C         call newdt (is  ,ie  ,je  ,je  ,k0  ,k2
C     &                  ,  imin,jmin,kmin
C     &                  ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
C        endif ! ldimen
        if(xtotnrg) then
         call newdt (is  ,ie  ,je  ,je  ,k0  ,k2
     &                  ,  imin,jmin,kmin
     &                  ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
        endif ! xtotnrg
       endif
CJH       if (nx3z .gt. 1) then
       if (ldimen .gt. 2) then
         call newdt (is  ,ie  ,js  ,je  ,k2+1,kmax
     &                  ,  imin,jmin,kmin
     &                  ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
       endif
c
c      Mark the velocity boundary values out of date.
c
       do 10 i = 1,6
         bvstat(i,3) = 0      !  v1
         bvstat(i,4) = 0      !  v2
         bvstat(i,5) = 0      !  v3
10     continue


c
c    3) Wait for communications to complete. 
c
#ifdef MPI_USED
       if(nreq .ne. 0) 
     .    call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
c Finally, do the remaining border zones.
c
       if (ke .gt. ks) then
         call newdt (is  ,ie  ,js  ,je  ,ks  ,ks  
     &                  ,  imin,jmin,kmin
     &                  ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
        if(xtotnrg) then
         call newdt (is  ,ie  ,js  ,je  ,ke  ,ke  
     &                  ,  imin,jmin,kmin
     &                  ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
        endif ! xtotnrg
       endif


c
c-----------------------------------------------------------------------
c
c Compute preliminary new time step.
c
#ifdef MPI_USED
c
c Now find the smallest dtnew among all tiles, and send the result
c to all in buf_out.  We need the info for the print statement below
c if the time step is too small, so use the MPI_MINLOC operation
c to pass values AND RANKS.
c
      if(xhydro) then
       buf_in(1) = dttoi2m
       buf_in(2) = real( myid )
       call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1
     &                   , MPI_2DOUBLE_PRECISION
     &                   , MPI_MAXLOC, comm3d, ierr)
       dttoi2m   =   buf_out(1)
c
       buf_in(1) = dtqqi2
       buf_in(3) = real( myid )
       call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1
     &                   , MPI_2DOUBLE_PRECISION
     &                   , MPI_MAXLOC, comm3d, ierr)
       dtqqi2    =   buf_out(1)
      endif
      if(lrad .ne. 0) then
       buf_in(1) = dtnri2
       buf_in(4) = real( myid )
       call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1
     &                   , MPI_2DOUBLE_PRECISION
     &                   , MPI_MAXLOC, comm3d, ierr)
       dtnri2    =   buf_out(1)
c
       buf_in(1) = dtimrdi2
       buf_in(5) = real( myid )
       call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1
     &                   , MPI_2DOUBLE_PRECISION
     &                   , MPI_MAXLOC, comm3d, ierr)
       dtimrdi2  =   buf_out(1)
      endif
#endif /* MPI */
c
      if(lrad .ne. 0) then
       if(xforce .eqv. .false.) then
        dtnew = courno / (sqrt( dtnri2 + dtimrdi2 ) + tiny )
       else ! xforce
        dtnew = courno / ( sqrt ( dttoi2m + dtqqi2 + dtnri2 +
     .                            dtimrdi2) + tiny )
       endif ! xforce
      endif ! lrad
c
      if(lrad .eq. 0) then
CPS
       if(xsubav) then
        dtnew  = courno / ( sqrt ( dttoi2m ) + tiny )
       else
        dtnew  = courno / ( sqrt ( dttoi2m + dtqqi2 ) + tiny )
       endif
C
      endif ! lrad
c
C#ifdef MPI_USED
Cc
Cc Now find the smallest dtnew among all tiles, and send the result
Cc to all in buf_out.  We need the info for the print statement below 
Cc if the time step is too small, so use the MPI_MINLOC operation
Cc to pass values AND RANKS.
Cc
C       buf_in(1) = dtnew
C       buf_in(2) = real( myid )
Cc  M-MML: on advice of Jakob Pichlmeier (Cray Munich) changed MPI_2FLOAT
Cc  to MPI_2DOUBLE_PRECISION...  26.2.98
C       call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1
C     &                   , MPI_2DOUBLE_PRECISION
C     &                   , MPI_MINLOC, comm3d, ierr)
C       dtnew  =   buf_out(1)
C#endif /* MPI */
       dt     =   min ( dtnew, 1.26*dt )
       if(time .ne. tlim) then
        if ((time+dt) .gt. tlim) dt = tlim-time
       endif
c
       if (dt .le. dtmin) then
c
c Determine which tile requires this short time step and report its
c coordinates and other info.  Convert 1/dt**2 to dt for each 
c physical process.
c

#ifdef MPI_USED
         if ( myid .eq. int( buf_out(2) ) ) then
c
           dtcs = 1.0 / ( sqrt ( dtcsm ) + tiny )
           dtv1 = 1.0 / ( sqrt ( dtv1m ) + tiny )
           dtv2 = 1.0 / ( sqrt ( dtv2m ) + tiny )
           dtv3 = 1.0 / ( sqrt ( dtv3m ) + tiny )
           dtal = 1.0 / ( sqrt ( dtalm ) + tiny )
           dtqq = 1.0 / ( sqrt ( dtqqi2) + tiny )
           write (2, 2010) coords(1), coords(2), coords(3)
     &                   , imin, jmin, kmin, nhy, dt, dtmin, dtcs
     &                   , dtv1, dtv2, dtv3, dtqq, dtal
         endif
         if ( myid .eq. int( buf_out(3) ) ) then
c
           dtcs = 1.0 / ( sqrt ( dtcsm ) + tiny )
           dtv1 = 1.0 / ( sqrt ( dtv1m ) + tiny )
           dtv2 = 1.0 / ( sqrt ( dtv2m ) + tiny )
           dtv3 = 1.0 / ( sqrt ( dtv3m ) + tiny )
           dtal = 1.0 / ( sqrt ( dtalm ) + tiny )
           dtqq = 1.0 / ( sqrt ( dtqqi2) + tiny )
           write (2, 2011) coords(1), coords(2), coords(3)
     &                   , imin, jmin, kmin, nhy, dt, dtmin, dtcs
     &                   , dtv1, dtv2, dtv3, dtqq, dtal
         endif
#endif /* MPI */
#ifndef MPI_USED
           dtcs = 1.0 / ( sqrt ( dtcsm ) + tiny )
           dtv1 = 1.0 / ( sqrt ( dtv1m ) + tiny )
           dtv2 = 1.0 / ( sqrt ( dtv2m ) + tiny )
           dtv3 = 1.0 / ( sqrt ( dtv3m ) + tiny )
           dtal = 1.0 / ( sqrt ( dtalm ) + tiny )
           dtqq = 1.0 / ( sqrt ( dtqqi2) + tiny )
           write (2, 2010) coords(1), coords(2), coords(3)
     &                   , imin, jmin, kmin, nhy, dt, dtmin, dtcs
     &                   , dtv1, dtv2, dtv3, dtqq, dtal
#endif /* NO MPI */
         nwarn = nwarn + 1
       endif
       return
c
c-----------------------------------------------------------------------
c----------------------- Write format statements -----------------------
c-----------------------------------------------------------------------
c
2010   format('NUDT    : **** VELOCITY **** Hot zone on tile (',i4,','
     &       ,i4,',',i4,')',/
     &       ,'NUDT    : at i=',i4,' j=',i4,' k=',i4,' (dt < dtmin)'
     &       ,', nhy  = ',i6,',',/
     2       ,'NUDT    : dt   = ',1pe12.5,', dtmin= ',1e12.5,', dtcs = '
     3       ,1e12.5,',',/
     4       ,'NUDT    : dtv1 = ',1e12.5,', dtv2 = ',1e12.5,', dtv3 = '
     5       ,1e12.5,',',/
     6       ,'NUDT    : dtqq = ',1e12.5,', dtal = '
     7       ,1e12.5,'.')
c
2011   format('NUDT    : **** VISCOSITY **** Hot zone on tile (',i4,','
     &       ,i4,',',i4,')',/
     &       ,'NUDT    : at i=',i4,' j=',i4,' k=',i4,' (dt < dtmin)'
     &       ,', nhy  = ',i6,',',/
     2       ,'NUDT    : dt   = ',1pe12.5,', dtmin= ',1e12.5,', dtcs = '
     3       ,1e12.5,',',/
     4       ,'NUDT    : dtv1 = ',1e12.5,', dtv2 = ',1e12.5,', dtv3 = '
     5       ,1e12.5,',',/
     6       ,'NUDT    : dtqq = ',1e12.5,', dtal = '
     7       ,1e12.5,'.')
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                  N U D T                  \\\\\\\\\\
c
c=======================================================================
c
c
