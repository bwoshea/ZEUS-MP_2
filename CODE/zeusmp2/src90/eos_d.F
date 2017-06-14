c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 E O S _ D                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c=======================================================================
      subroutine eos_d
c
c     Driver for non-ideal equations of state (leos > 1).  The logic
c     is cloned from the forces_d routine written by R. A. Fiedler.
c     Note that the EOS must be computed in one layer of inner ghost
c     zones on each face so that pressure gradients at the inner
c     boundaries can be computed by FORCES.
c
c     Written by J. Hayes, 5-2003
c
      use real_prec
      use config
      use param
      use grid
      use root
      use field
      use bndry
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
      integer :: i, j, k, k1, k2, imax, jmax, kmax
c
c-----------------------------------------------------------------------
c
      if(.not. xtotnrg) then
       imax = ie
       jmax = je
       kmax = ke
      else
       imax = ie-2
       jmax = je-2
       kmax = ke-2
      endif
      if(ldimen .eq. 2) go to 222
      if(ldimen .eq. 1) go to 111
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     3D Grids
c-----------------------------------------------------------------------
c
c Divide the computational volume into three equal pieces.  We must
c have at least 5 active zones in the 3-direction.
c
       k1 = int( real( ke - ks + 1 ) / 3.0 ) + ks
       k2 = int( real( ke - ks + 1 ) / 3.0 ) + k1
c
c
c    2) Do first portion of the interior points.
c
c
c i boundaries
c
      nreq = 0
      nsub = nsub + 1
      call bvald(1,0,0,0,0,0,d )
      call bvale(1,0,0,0,0,0,e )
      if(xtotnrg) then
       call bvalv1(1,1,0,0,0,0,v1)
       call bvalv2(1,1,0,0,0,0,v2)
       call bvalv3(1,1,0,0,0,0,v3)
      endif
      if(leos .ne. 1) call bvalt(1,0,0,0,0,0,tt)
      if(nspec .gt. 1) call bvalabuns(1,0,0,0,0,0,abun)
c
      call eos(is+1,imax,js+1,jmax,ks+1,k1)
c
#ifdef MPI_USED
c
c    3) Wait for communications to complete.
c
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
c j boundaries
c
      nreq = 0
      nsub = nsub + 1
      call bvald(0,0,1,0,0,0,d )
      call bvale(0,0,1,0,0,0,e )
      if(xtotnrg) then
       call bvalv1(0,0,1,1,0,0,v1)
       call bvalv2(0,0,1,1,0,0,v2)
       call bvalv3(0,0,1,1,0,0,v3)
      endif
      if(leos .ne. 1) call bvalt(0,0,1,0,0,0,tt)
      if(nspec .gt. 1) call bvalabuns(0,0,1,0,0,0,abun)
c
c    2) Do middle portion of the interior points, plus some on borders.
c
      call eos(is  ,is  ,js+1,jmax,ks+1,k1)
      call eos(is  ,ie  ,js+1,jmax,k1+1,k2)
#ifdef MPI_USED
c
c    3) Wait for communications to complete.
c
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
c k boundaries
c
c    1) Post sends and receives.  We need d(k+1) for ZRP or RTP.
c
      nreq = 0
      nsub = nsub + 1
      call bvald(0,0,0,0,1,0,d )
      call bvale(0,0,0,0,1,0,e )
      if(xtotnrg) then
       call bvalv1(0,0,0,0,1,1,v1)
       call bvalv2(0,0,0,0,1,1,v2)
       call bvalv3(0,0,0,0,1,1,v3)
      endif
      if(leos .ne. 1) call bvalt(0,0,0,0,1,0,tt)
      if(nspec .gt. 1) call bvalabuns(0,0,0,0,1,0,abun)
c
c    2) Do last portion of the interior points, plus some on borders.
c
      call eos(is  ,ie  ,js  ,js  ,ks+1,k2)
      call eos(is  ,ie  ,js  ,je  ,k2+1,kmax)
c
c      Mark the boundary values out of date.
c
      do 10 i = 1,6
       bvstat(i,7) = 0      !  v1
       bvstat(i,8) = 0      !  v2
10    continue
#ifdef MPI_USED
c
c    3) Wait for communications to complete. 
c
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
c Finally, do the remaining border zones.
c
      call eos(is-1,is  ,js  ,je  ,ks  ,ke) ! 1-face + ghost layer
      call eos(is  ,ie  ,js-1,js  ,ks  ,ke) ! 2-face + ghost layer
      call eos(is  ,ie  ,js  ,je  ,ks-1,ks) ! 3-face + ghost layer
      if(xtotnrg) then
       call eos(ie-1,ie  ,js  ,je  ,ks  ,ke-2)
       call eos(is  ,ie-2,je-1,je  ,ks  ,ke-2)
       call eos(is  ,ie  ,js  ,je  ,ke-1,ke  )
      endif
      go to 999
c
c-----------------------------------------------------------------------
c     2D Grids
c-----------------------------------------------------------------------
c
222   continue
      if(xtotnrg) then
       do i = 1, 6
        bvstat(i,3) = 0
        bvstat(i,4) = 0
       enddo
      endif
c
c i boundaries
c
      nreq = 0
      nsub = nsub + 1
      call bvald(1,0,0,0,0,0,d )
      call bvale(1,0,0,0,0,0,e )
      if(xtotnrg) then
       call bvalv1(1,1,0,0,0,0,v1)
       call bvalv2(1,1,0,0,0,0,v2)
      endif
      if(leos .ne. 1) call bvalt(1,0,0,0,0,0,tt)
      if(nspec .gt. 1) call bvalabuns(1,0,0,0,0,0,abun)
c
c    2) Do first portion of the interior points.
c
      call eos(is+1,imax,js+1,jmax,ks,ks)
c
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
c j boundaries
c
      nreq = 0
      nsub = nsub + 1
      call bvald(0,0,1,0,0,0,d )
      call bvale(0,0,1,0,0,0,e )
      if(xtotnrg) then
       call bvalv1(0,0,1,1,0,0,v1)
       call bvalv2(0,0,1,1,0,0,v2)
      endif
      if(leos .ne. 1) call bvalt    (0,0,1,0,0,0,tt)
      if(nspec .gt. 1) call bvalabuns(0,0,1,0,0,0,abun)
c
      call eos(is  ,is  ,js+1,jmax,ks,ks)
      if(xtotnrg) call eos(ie-1,ie  ,js+1,jmax,ks,ks)
c
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c
      call eos(is  ,ie  ,js  ,js  ,ks,ks)
c
      if(xtotnrg) call eos(is  ,ie  ,je-1,je  ,ks,ks)
c
      call eos(is-1,is-1,js  ,je  ,ks,ks) ! 1-dir inner ghost layer
      call eos(is  ,ie  ,js-1,js-1,ks,ks) ! 2-dir inner ghost layer
c
      if(.false.) then
       call eos(ie+1,ie+1,js  ,je  ,ks,ks) ! 1-dir outer ghost layer
       call eos(is  ,ie  ,je+1,je+1,ks,ks) ! 2-dir outer ghost layer
      endif
c
c      Mark the boundary values out of date.
c
       do i = 1,6
         bvstat(i,7) = 0      !  T
         bvstat(i,8) = 0      !  abun
       enddo
c
      go to 999
c-----------------------------------------------------------------------
c     1D Grids
c-----------------------------------------------------------------------
c
111   continue
c
c i boundaries
c
      nreq = 0
      nsub = nsub + 1
      call bvald(1,0,0,0,0,0,d )
      call bvale(1,0,0,0,0,0,e )
      if(xtotnrg) then
       call bvalv1(1,1,0,0,0,0,v1)
      endif
      if(leos .ne. 1) call bvalt    (1,0,0,0,0,0,tt)
      if(nspec .gt. 1) call bvalabuns(1,0,0,0,0,0,abun)
c
c    2) Do first portion of the interior points.
c
      call eos(is+1,imax,js,js,ks,ks)
c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c
c     do first i zone plus inner ghost zone
c
      call eos(is-1,is,js,js,ks,ks)
      if(xtotnrg) call eos(ie-1,ie,js,js,ks,ks)
c
c      Mark the boundary values out of date.
c
      do i = 1,6
       bvstat(i,7) = 0 
       bvstat(i,8) = 0 
      enddo
c
999   return
      end
