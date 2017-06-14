c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              F O R C E S _ D              \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine forces_d
     1 (v1old, v2old, v3old, v1new, v2new, v3new)
c
c  PURPOSE
c    Driver for forces module.  Updates the velocities due to the
c    following body forces:
c
c    1) Thermal pressure gradient
c    2) Self-gravity
c    3) Rotational Pseudo-Forces
c    4) Gravitational point mass
c    5) Radiation force
c    6) Magnetic pressure
c
c  AUTHOR
c    Robert A. Fiedler
c
c  LAST MODIFIED
c    02/17/97, MHD on 3 Mar 98 by M-M Mac Low
c
c  INPUT
c    v*old    velocity components before acceleration
c
c  OUTPUT
c    v*new    velocity components after update
c
c  EXTERNALS:
c    BVALV1  , BVALV2  , BVALV3
c    BVALD   , BVALE   , BVALER
c    FORCES
c
c.......................................................................
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
      real(rl) :: v1old(in,jn,kn), v2old(in,jn,kn), v3old(in,jn,kn),
     .            v1new(in,jn,kn), v2new(in,jn,kn), v3new(in,jn,kn)
c
      real(rl) :: tg
      integer  :: sf
c
      integer :: i, j, k, k1, k2, imax, jmax, kmax
c
c-----------------------------------------------------------------------
c      Not doing forces, so bail
c-----------------------------------------------------------------------
c
      if(xforce .eqv. .false.) then
c
c Just copy the old velocities to the new for all active zones.
c
       do k=ks,ke
         do j=js,je
           do i=is,ie
             v1new(i,j,k) = v1old(i,j,k)
             v2new(i,j,k) = v2old(i,j,k)
             v3new(i,j,k) = v3old(i,j,k)
           enddo ! i
         enddo ! j
       enddo ! k
c
       do i = 1, 6
        bvstat(i,3) = 0
        bvstat(i,4) = 0
        bvstat(i,5) = 0
       enddo
c
       return
      endif ! xforce
c
c-----------------------------------------------------------------------
      if(ldimen .eq. 2) go to 222
      if(ldimen .eq. 1) go to 111
c-----------------------------------------------------------------------
c
c Divide the computational volume into three equal pieces.  We must
c have at least 5 active zones in the 3-direction.
c
       k1 = int( real( ke - ks + 1 ) / 3.0 ) + ks
       k2 = int( real( ke - ks + 1 ) / 3.0 ) + k1
c
c......................................................................
c
c i boundaries
c
      nreq = 0
      nsub = nsub + 1
      call bvald  (1,0,0,0,0,0,d )
      if(xiso .eqv. .false.) call bvale  (1,0,0,0,0,0,e )
      imax = ie
      if(lgeom .eq. 3) then
       jmax = je - 1
       call bvalv2 (1,0,0,0,0,0,v2old)
       call bvalv3 (1,0,0,0,0,0,v3old)
      else
       jmax = je
      endif ! RTP
      if(lrad .ne. 0) call bvalers (1,0,0,0,0,0,er)
c
c    2) Do first portion of the interior points.
c
      call forces (is+1,imax,js+1,jmax,ks+1,k1
     .             ,v1old,v2old,v3old,v1new,v2new,v3new)
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
c    1) Post sends and receives.  We need d(j+1) only for RTP.
c
      nreq = 0
      nsub = nsub + 1
c
      if(lgeom .eq. 3) then
       call bvald  (0,0,1,1,0,0,d )
      else
       call bvald  (0,0,1,0,0,0,d )
      endif ! lgeom
c
      if(xiso .eqv. .false.) call bvale  (0,0,1,0,0,0,e )
      if(lgeom .eq. 3) then
       call bvalv2 (0,0,0,1,0,0,v2old)
      endif
      if(lgeom .ne. 1) then
       call bvalv3 (0,0,1,0,0,0,v3old)
      endif
      if(lrad .ne. 0) call bvalers (0,0,1,0,0,0,er)
c
c    2) Do middle portion of the interior points, plus some on borders.
c
      call forces (is  ,is  ,js+1,jmax,ks+1,k1
     .            ,v1old,v2old,v3old,v1new,v2new,v3new)
      call forces (is  ,ie  ,js+1,jmax,k1+1,k2
     .            ,v1old,v2old,v3old,v1new,v2new,v3new)
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
c
      if(lgeom .ne. 1) then
       call bvald(0,0,0,0,1,1,d )
      else
       call bvald(0,0,0,0,1,0,d )
      endif
c
      if(xiso .eqv. .false.) call bvale  (0,0,0,0,1,0,e )
      if(lgeom .ne. 1) then
       kmax = ke - 1
       call bvalv3 (0,0,0,0,1,1,v3old)
      else ! lgeom
       kmax = ke
      endif ! lgeom
      if(lrad .ne. 0) call bvalers (0,0,0,0,1,0,er)
c
c    2) Do last portion of the interior points, plus some on borders.
c
      call forces (is  ,ie  ,js  ,js  ,ks+1,k2
     &             ,v1old,v2old,v3old,v1new,v2new,v3new)
c
      if(lgeom .eq. 3) then
       call forces (is  ,ie  ,je  ,je  ,ks+1,k2
     &             ,v1old,v2old,v3old,v1new,v2new,v3new)
      endif ! rtp
c
      call forces (is  ,ie  ,js  ,je  ,k2+1,kmax
     &             ,v1old,v2old,v3old,v1new,v2new,v3new)
c
c      Mark the velocity boundary values out of date.
c
       do 10 i = 1,6
         bvstat(i,3) = 0      !  v1
         bvstat(i,4) = 0      !  v2
         bvstat(i,5) = 0      !  v3
10     continue
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
      call forces (is  ,ie  ,js  ,je  ,ks  ,ks  
     .             ,v1old,v2old,v3old,v1new,v2new,v3new)
c
      if(lgeom .ne. 1) then
       call forces (is  ,ie  ,js  ,je  ,ke  ,ke  
     &             ,v1old,v2old,v3old,v1new,v2new,v3new)
      endif ! xtotnrg or zrp/rtp
c
      go to 999
c
c......................................................................
c     2D TRANSPORT
c......................................................................
222   continue
c
c i boundaries
c
c    1) Post sends and receives.  We never use d(i+1) in any of the
c       source term substeps, so we need to pass only 1 "m" layer
c       of density boundary values.  We also need only one
c       "m" layer of energy values.  If we need the velocity, we
c       need just one layer, but both "m" and "p".  This should
c       minimize communication for any physics/geometry.  By
c       exchanging the i, j, and k boundary data in three separate 
c       steps, we ensure that the corner and edge data are correctly
c       updated.
c
      nreq = 0
      nsub = nsub + 1
      if(leos .ne. 2) call bvald  (1,0,0,0,0,0,d )
c
      if(xiso .eqv. .false.) call bvale  (1,0,0,0,0,0,e )
      imax = ie
c
      if(lgeom .eq. 3) then
       jmax = je - 1
       call bvalv2 (1,0,0,0,0,0,v2old)
       call bvalv3 (1,0,0,0,0,0,v3old)
      else
       jmax = je
      endif ! rtp
c
      if(lrad .ne. 0) call bvalers (1,0,0,0,0,0,er)
c
c    2) Do first portion of the interior points.
c
      call forces (is+1,imax,js+1,jmax,ks,ks
     .             ,v1old,v2old,v3old,v1new,v2new,v3new)
c
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
c j boundaries
c
c    1) Post sends and receives.  We need d(j+1) only for RTP.
c
      nreq = 0
      nsub = nsub + 1
c
      if(lgeom .eq. 3) then
       call bvald  (0,0,1,1,0,0,d )
      else
       call bvald  (0,0,1,0,0,0,d )
      endif ! lgeom
      if(xiso .eqv. .false.) call bvale  (0,0,1,0,0,0,e )
c
      if(lgeom .eq. 3) then
       call bvalv2 (0,0,0,1,0,0,v2old)
      endif
c
      if(lgeom .ne. 1) then
       call bvalv3 (0,0,1,0,0,0,v3old)
      endif
c
      if(lrad .ne. 0) call bvalers (0,0,1,0,0,0,er)
c
c    2) Do middle portion of the interior points, plus some on borders.
c
      call forces (is  ,is  ,js+1,jmax,ks,ks
     .             ,v1old,v2old,v3old,v1new,v2new,v3new)
c
c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c
c    2) Do last portion of the interior points, plus some on borders.
c
      call forces (is  ,ie  ,js  ,js  ,ks,ks
     .             ,v1old,v2old,v3old,v1new,v2new,v3new)
c
      if(lgeom .eq. 3) then
       call forces (is  ,ie  ,je  ,je  ,ks,ks
     &             ,v1old,v2old,v3old,v1new,v2new,v3new)
      endif ! rtp
c
c      Mark the velocity boundary values out of date.
c
       do i = 1,6
         bvstat(i,3) = 0      !  v1
         bvstat(i,4) = 0      !  v2
         bvstat(i,5) = 0      !  v2
       enddo
c
      go to 999
c......................................................................
c     1D TRANSPORT
c......................................................................
111   continue
c
c i boundaries
c
c    1) Post sends and receives.  We never use d(i+1) in any of the
c       source term substeps, so we need to pass only 1 "m" layer
c       of density boundary values.  We also need only one
c       "m" layer of energy values.  If we need the velocity, we
c       need just one layer, but both "m" and "p".  This should
c       minimize communication for any physics/geometry.  By
c       exchanging the i, j, and k boundary data in three separate 
c       steps, we ensure that the corner and edge data are correctly
c       updated.
c
      nreq = 0
      nsub = nsub + 1
      call bvald  (1,0,0,0,0,0,d )
c
      if(xiso .eqv. .false.) call bvale  (1,0,0,0,0,0,e )
      imax = ie
c
      if(lrad .ne. 0) call bvalers (1,0,0,0,0,0,er)
c
c    2) Do first portion of the interior points.
c
      call forces (is+1,imax,js,js,ks,ks
     .             ,v1old,v2old,v3old,v1new,v2new,v3new)
c
c       subroutine forces (ibeg,iend,jbeg,jend,kbeg,kend
c     &                   ,u1,u2,u3,w1,w2,w3)
cc
cc Arrays u1 , u2 , u3  hold the old velocity values, while
cc        w1 , w2 , w3  receive the updated values.
c
c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
      call forces (is  ,is  ,js,js,ks,ks
     .             ,v1old,v2old,v3old,v1new,v2new,v3new)
c
c      Mark the velocity boundary values out of date.
c
       do i = 1,6
         bvstat(i,3) = 0      !  v1
       enddo
c
999   return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////             F O R C E S _ D               \\\\\\\\\\
c
c=======================================================================
c
