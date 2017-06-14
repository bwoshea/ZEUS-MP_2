c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////               A V I S C _ D
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine avisc_d
     1            (v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
c
c  PURPOSE
c    Driver for artificial viscosity module.  Updates velocity
c    components and material energy.
c
c  AUTHOR
c    Robert A. Fiedler
c
c  LAST MODIFIED
c    01/21/97
c
c  INPUT
c    v*old     Velocity components before viscous update.
c
c  OUPUT
c    v*new     Velocity components after  viscous update.
c    s*        Momentum components for use in the transport step.
c
c  EXTERNALS:
c    BVALV1  , BVALV2  , BVALV3
c    BVALD
c    AVISC
c
c.......................................................................
c
      use real_prec
      use config
      use param
      use root
      use field
      use bndry
      use grid
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit none
c
      real(rl) :: v1old(in,jn,kn),v2old(in,jn,kn),v3old(in,jn,kn),
     .            v1new(in,jn,kn),v2new(in,jn,kn),v3new(in,jn,kn),
     .            s1   (in,jn,kn),s2   (in,jn,kn),s3   (in,jn,kn)
c
      integer  :: i, k1, k2, imax, jmax, kmax
      real(rl) :: dvdxmn
c
c-----------------------------------------------------------------------
c
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
c
c Save the largest velocity gradient for computing the viscous
c time step.
c
      dvdxmn = 0.0
c
      nsub = nsub + 1
      if(xtotnrg) then
       call bvald  (1,1,0,0,0,0,d    )
       call bvalv1 (3,3,0,0,0,0,v1old)
      else
       call bvald  (1,0,0,0,0,0,d    )
       call bvalv1 (1,1,0,0,0,0,v1old)
      endif
c
      call avisc (is+1,ie-2,js+1,je-2,ks+1,k1  ,dvdxmn
     &           ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
c
#ifdef MPI_USED
c
c    3) Wait for communications to complete.
c
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
       nreq = 0
       nsub = nsub + 1
       if(xtotnrg) then
        call bvald  (0,0,1,1,0,0,d    )
        call bvalv2 (0,0,3,3,0,0,v2old)
       else
        call bvald  (0,0,1,0,0,0,d    )
        call bvalv2 (0,0,1,1,0,0,v2old)
       endif
c
       call avisc (is  ,is  ,js+1,je-2,ks+1,k1  ,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (ie-1,ie  ,js+1,je-2,ks+1,k1  ,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (is  ,ie  ,js+1,je-2,k1+1,k2  ,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)


#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
       nreq = 0
       nsub = nsub + 1
       if(xtotnrg) then
        call bvald (0,0,0,0,1,1,d    )
        call bvalv3(0,0,0,0,3,3,v3old)
       else
        call bvald (0,0,0,0,1,0,d    )
        call bvalv3(0,0,0,0,1,1,v3old)
       endif
c
       call avisc (is  ,ie  ,js  ,js  ,ks+1,k2  ,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (is  ,ie  ,je-1,je  ,ks+1,k2  ,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (is  ,ie  ,js  ,je  ,k2+1,ke-2,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
c
c      Mark the boundary values out of date.
c
       do 20 i = 1,6
         if(xiso .eqv. .false.) bvstat(i,2) = 0      !  e
         bvstat(i,3) = 0      !  v1
         bvstat(i,4) = 0      !  v2
         bvstat(i,5) = 0      !  v3
20     continue
#ifdef MPI_USED
c
c    3) Wait for communications to complete.
c
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
       call avisc (is  ,ie  ,js  ,je  ,ks  ,ks  ,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (is  ,ie  ,js  ,je  ,ke-1,ke  ,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
c
      go to 999
c-----------------------------------------------------------------------
222   continue
c
c Save the largest velocity gradient for computing the viscous
c time step.
c
       dvdxmn = 0.0
c
       nreq = 0
       nsub = nsub + 1
       if(xtotnrg) then
        call bvald  (1,1,0,0,0,0,d    )
        call bvalv1 (3,3,0,0,0,0,v1old)
       else
        call bvald  (1,0,0,0,0,0,d    )
        call bvalv1 (1,1,0,0,0,0,v1old)
       endif
c
       call avisc (is+1,ie-2,js+1,je-2,ks,ks  ,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
c
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
       nreq = 0
       nsub = nsub + 1
       if(xtotnrg) then
        call bvald  (0,0,1,1,0,0,d    )
        call bvalv2 (0,0,3,3,0,0,v2old)
       else
        call bvald  (0,0,1,0,0,0,d    )
        call bvalv2 (0,0,1,1,0,0,v2old)
       endif
c
       call avisc (is  ,is  ,js+1,je-2,ks,ks  ,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (ie-1,ie  ,js+1,je-2,ks,ks  ,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
c
c      Mark the boundary values out of date.
c
       do i = 1,6
         if(xiso .eqv. .false.) bvstat(i,2) = 0      !  e
         bvstat(i,3) = 0      !  v1
         bvstat(i,4) = 0      !  v2
         bvstat(i,5) = 0      !  v3
       enddo
c
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
       call avisc (is  ,ie  ,js  ,js  ,ks  ,ks  ,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (is  ,ie  ,je-1,je  ,ks  ,ks  ,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
c
      go to 999
c-----------------------------------------------------------------------
111   continue
c
c Save the largest velocity gradient for computing the viscous
c time step.
c
       dvdxmn = 0.0
c
c i boundaries
c
c    1) Post sends and receives.
c
       nreq = 0
       nsub = nsub + 1
       if(xtotnrg) then
        call bvald  (1,1,0,0,0,0,d    )
        call bvalv1 (3,3,0,0,0,0,v1old)
       else
        call bvald  (1,0,0,0,0,0,d    )
        call bvalv1 (1,1,0,0,0,0,v1old)
       endif


c
c    2) Do first portion of the interior points.
c
       call avisc (is+1,is+3,js,js,ks,ks  ,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
c
c       subroutine avisc (ibeg,iend,jbeg,jend,kbeg,kend,dvdxmn
c     &                  ,w1,w2,w3,u1,u2,u3,s1,s2,s3)
c
cc    w1,w2,w3   velocity values prior to viscosity update.
cc    u1,u2,u3   velocity values after    viscosity update.
cc    s1,s2,s3   updated momentum densities for transport step.
c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c
c      Mark the boundary values out of date.
c
       do i = 1,6
         if(xiso .eqv. .false.) bvstat(i,2) = 0      !  e
         bvstat(i,3) = 0      !  v1
         bvstat(i,4) = 0      !  v2
         bvstat(i,5) = 0      !  v3
       enddo
c......................................................................
c
c Do the remaining border zones.
c
       call avisc (is  ,is  ,js  ,js  ,ks  ,ks  ,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
       call avisc (is+4,ie  ,js  ,js  ,ks  ,ks  ,dvdxmn
     &            ,v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
c
c......................................................................
999   continue
c......................................................................
c
c  Compute the viscous timestep.  Note that the minimum dv/dx is found
c  since it is less than 0.  Thus the minimum dv/dx gives the maximum
c  absolute value.  We'll do a global min on this later, in nudt.
c
       dtqqi2 = ( 4.0 * qcon * dvdxmn )**2
c
       return
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////               A V I S C _ D               \\\\\\\\\\
c
c=======================================================================
c
