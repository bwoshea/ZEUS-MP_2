c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              L O R E N T Z _ D            \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine lorentz_d (w1, w2, w3, u1, u2, u3)
c
c  PURPOSE
c    Driver for lorentz update. Updates the velocities due to the
c    transverse magnetic forces (the longitudinal forces are handled 
c    in FORCES).  Adapted from FORCES_D; but LORENTZ needs two slabs 
c    of velocity both above and below
c
c  AUTHOR
c    Mordecai-Mark Mac Low
c
c  LAST MODIFIED
c    3 Mar 98 by M-M Mac Low
c
c    April 29, 2006 by John Hayes:  consolidated BVAL* and LORENTZ
c    calls so that all BVAL* updates are completed before calling
c    subroutine LORENTZ.  LORENTZ is therefore called one time with
c    all zones being processed as a group.  The original form of this
c    routine attempted to subdivide the LORENTZ update and perform
c    it in stages, alternating calls to LORENTZ with calls to the
c    BVAL* routines on each tile face.  Unfortunately, this does not work
c    properly in parallel.  It may be possible to correct the problem,
c    but as of this date it remains unsolved.  This
c    driver routine, in its simplified form, does produce proper
c    results (when compare to a 1-processor run) when used in a
c    parallel calculation.  -- JCH 04/29/2006
c
c    May 15, 2006 by John Hayes: changed velocity references to point
c    to local scratch arrays w[123] and u[123].  Only density is accessed
c    via module FIELD.  Calls to BVALV[123] and LORENTZ are affected.
c    These changes, in combination with those made in CT/CT_2D/CT_1D,
c    correct an error in which two variables were assigned the same
c    space in memory.
c
c  INPUT
c    w[1,2,3]    velocity components before acceleration
c
c  OUTPUT
c    u[1,2,3]    velocity components after update
c
c  EXTERNALS:
c    BVALV1  , BVALV2  , BVALV3
c    BVALD   , BVALE  
c    LORENTZ
c
c.......................................................................
c
      use real_prec
      use config
      use param
      use grid
      use root
      use field, ONLY : d
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
      integer  :: i, j, k, k1, k2
      real(rl) :: u1(in,jn,kn), u2(in,jn,kn), u3(in,jn,kn)
      real(rl) :: w1(in,jn,kn), w2(in,jn,kn), w3(in,jn,kn)
c
      if(xforce .eqv. .false.) then
c
c Just copy the old velocities to the new for all active zones.
c
      do k=ks,ke
       do j=js,je
         do i=is,ie
          u1(i,j,k) = w1(i,j,k)
          u2(i,j,k) = w2(i,j,k)
          u3(i,j,k) = w3(i,j,k)
         enddo ! i
        enddo ! j
       enddo ! k
       return
      endif ! xforce

c-----------------------------------------------------------------------
c
c Divide the computational volume into three equal pieces.  We must
c have at least 5 active zones in the 3-direction.
c
      k1 = int( real( ke - ks + 1 ) / 3.0 ) + ks
      k2 = int( real( ke - ks + 1 ) / 3.0 ) + k1
c
c
c     asif: accomodate 2D MHD runs
c
      if(ldimen.eq.2) goto 222
c
c     JCH: accomodate 1D MHD runs
c
      if(ldimen.eq.1) goto 111
c......................................................................
c
c i boundaries
c
       nreq = 0
       nsub = nsub + 1
       call bvald  (3,1,0,0,0,0,d )
c
       call bvalv1 (3,3,0,0,0,0,w1)
       call bvalv2 (3,3,0,0,0,0,w2)
       call bvalv3 (3,3,0,0,0,0,w3)
#ifdef MPI_USED
      if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c
c j boundaries
c
       nreq = 0
       nsub = nsub + 1
       call bvald  (0,0,3,1,0,0,d )
       call bvalv1 (0,0,3,3,0,0,w1)
       call bvalv2 (0,0,3,3,0,0,w2)
       call bvalv3 (0,0,3,3,0,0,w3)
#ifdef MPI_USED
      if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c
c k boundaries
c
       nreq = 0
       nsub = nsub + 1
       call bvald  (0,0,0,0,3,1,d )
       call bvalv1 (0,0,0,0,3,3,w1)
       call bvalv2 (0,0,0,0,3,3,w2)
       call bvalv3 (0,0,0,0,3,3,w3)
#ifdef MPI_USED
      if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c
       call lorentz (is,ie,js,je,ks,ke,w1,w2,w3,u1,u2,u3)
c
       goto 999
222     continue
c	asif: 2D lorentz
c......................................................................
c
c i boundaries
c
       nreq = 0
       nsub = nsub + 1
       call bvald  (3,1,0,0,0,0,d )
       call bvalv1 (3,3,0,0,0,0,w1)
       call bvalv2 (3,3,0,0,0,0,w2)
       call bvalv3 (3,3,0,0,0,0,w3)
#ifdef MPI_USED
      if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c
c j boundaries
c
       nreq = 0
       nsub = nsub + 1
       call bvald  (0,0,3,1,0,0,d )
       call bvalv1 (0,0,3,3,0,0,w1)
       call bvalv2 (0,0,3,3,0,0,w2)
       call bvalv3 (0,0,3,3,0,0,w3)
#ifdef MPI_USED
      if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c
       call lorentz_2D(is,ie,js,je,ks,ks,w1,w2,w3,u1,u2,u3)
c
       goto 999
c
111     continue
c
c      1-D Lorentz driver (JCH)
c......................................................................
c
c i boundaries
c
       nreq = 0
       nsub = nsub + 1
       call bvald  (3,1,0,0,0,0,d )
       call bvalv1 (3,3,0,0,0,0,w1)
       call bvalv2 (3,3,0,0,0,0,w2)
       call bvalv3 (3,3,0,0,0,0,w3)
#ifdef MPI_USED
      if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c
       call lorentz_1D(is,ie,js,js,ks,ks,w1,w2,w3,u1,u2,u3)
c
999    continue
c
c      Mark the velocities as out of date.
c
      do  i = 1,6
       bvstat(i,1) = 0      !  d
       bvstat(i,3) = 0      !  w1
       bvstat(i,4) = 0      !  w2
       bvstat(i,5) = 0      !  w3
      enddo
c
      return
c
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////           L O R E N T Z _ D               \\\\\\\\\\
c
c=======================================================================
