c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 P D V _ D                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine pdv_d
     1            (dcopy, eod)
c
c  PURPOSE
c    Compute the compressional work term (PdV) in the material
c    energy equation (for no radiation).
c
c  AUTHOR
c    Robert A. Fiedler
c
c  LAST MODIFIED
c    01/21/97
c
c  INPUT
c    (none)
c
c  OUTPUT
c    dcopy     A copy of the density.
c    eod       e/d for the transport step.
c
c  We need just 1 "p" layer of updated boundary data for v1, v2, and v3,
c  but none for d and e.  
c
c  Routine pdv also saves the density and e/d in (dcopy,eod) for use
c  in the transport step.  This is why it is being called even when
c  TOTAL_ENERGY is defined.
c
c  EXTERNALS:
c    BVALV1  , BVALV2  , BVALV3
c    PDV
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
      implicit NONE
c
      real(rl) :: dcopy(in,jn,kn), eod(in,jn,kn)
c
      integer  :: i, k1, k2, imax, jmax, kmax
      integer  :: j, k
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
c i boundaries
c
c    1) Post sends and receives.
c
       nreq = 0
       nsub = nsub + 1
       call bvalv1 (0,1,0,0,0,0,v1)


c
c    2) Do first portion of the interior points.
c
c      Copy density boundary values.
c
       do k=ks-2,ks-1
         do j=js-2,je+2
           do i=is-2,ie+2
             dcopy(i,j,k) = d(i,j,k)
           enddo ! i
         enddo ! j
       enddo ! k
c
       call pdv (is,ie-1,js,je-1,ks  ,k1
     &          ,dcopy ,eod )
c
c       subroutine pdv (ibeg,iend,jbeg,jend,kbeg,kend
c     &                ,dlo,eod)
cc
cc    dlo   Mass density (copy)                      for transport step.
cc    eod   Specific energy density e/d (or (e+p)/d] for transport step.
c
c
c      Copy density boundary values.
c
       do k=ks,ke
         do j=js-2,js-1
           do i=is-2,ie+2
             dcopy(i,j,k) = d(i,j,k)
           enddo ! i
         enddo ! j
       enddo ! k


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
       nreq = 0
       nsub = nsub + 1
       call bvalv2 (0,0,0,1,0,0,v2)


c
c    2) Do middle 1/3 of the interior points, and some on borders.
c
       call pdv (ie  ,ie  ,js,je-1,ks  ,k1
     &          ,dcopy ,eod )
       call pdv (is  ,ie  ,js,je-1,k1+1,k2
     &          ,dcopy ,eod )
c
c      Copy density boundary values.
c
       do k=ks,ke
         do j=je+1,je+2
           do i=is-2,ie+2
             dcopy(i,j,k) = d(i,j,k)
           enddo ! i
         enddo ! j
       enddo ! k
       do k=ks,ke
         do j=js,je
           do i=is-2,is-1
             dcopy(i,j,k) = d(i,j,k)
           enddo ! i
         enddo ! j
       enddo ! k


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
       nreq = 0
       nsub = nsub + 1
       call bvalv3 (0,0,0,0,0,1,v3)


c
c    2) Do last 1/3 of the interior points, and some on borders.
c
       call pdv (is  ,ie  ,je  ,je  ,ks  ,k2
     &          ,dcopy ,eod )
       call pdv (is  ,ie  ,js  ,je  ,k2+1,ke-1   
     &          ,dcopy ,eod )
c
c      Copy density boundary values.
c
       do k=ke+1,ke+2
         do j=js-2,je+2
           do i=is-2,ie+2
             dcopy(i,j,k) = d(i,j,k)
           enddo ! i
         enddo ! j
       enddo ! k
       do k=ks,ke
         do j=js,je
           do i=ie+1,ie+2
             dcopy(i,j,k) = d(i,j,k)
           enddo ! i
         enddo ! j
       enddo ! k
c
c    3) Wait for communications to complete. 
c
#ifdef MPI_USED
       if(nreq .ne. 0) 
     .    call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
c Finally, do the remaining border.
c
       call pdv (is  ,ie  ,js  ,je  ,ke  ,ke
     &          ,dcopy ,eod )
c
      go to 999
c======================================================================
c     2D TRANSPORT
c======================================================================
c
222   continue
c
c i boundaries
c
c    1) Post sends and receives.
c
       nreq = 0
       nsub = nsub + 1
       call bvalv1 (0,1,0,0,0,0,v1)


c
c    2) Do first portion of the interior points.
c
c      Copy density values.
c
         do j=js-2,je+2
           do i=is-2,ie+2
             dcopy(i,j,ks) = d(i,j,ks)
           enddo ! i
         enddo ! j
c
       call pdv (is,ie-1,js,je-1,ks  ,ks
     &          ,dcopy ,eod )
c
c       subroutine pdv (ibeg,iend,jbeg,jend,kbeg,kend
c     &                ,dlo,eod)
cc
cc    dlo   Mass density (copy)                      for transport step.
cc    eod   Specific energy density e/d (or (e+p)/d] for transport step.
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
       nreq = 0
       nsub = nsub + 1
       call bvalv2 (0,0,0,1,0,0,v2)
c
       call pdv (ie  ,ie  ,js,je-1,ks  ,ks
     &          ,dcopy ,eod )
c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
       if(nreq .ne. 0) 
     .    call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c
c Finally, do the remaining border.
c
       call pdv (is  ,ie  ,je  ,je  ,ks  ,ks
     &          ,dcopy ,eod )
c
c      Mark the boundary values out of date.
c
      go to 999
c======================================================================
c     1D TRANSPORT
c======================================================================
c
111   continue
c
c i boundaries
c
c    1) Post sends and receives.
c
       nreq = 0
       nsub = nsub + 1
       call bvalv1 (0,1,0,0,0,0,v1)


c
c    2) Do first portion of the interior points.
c
c      Copy density values.
c
         do j=js,js
           do i=is-2,ie+2
             dcopy(i,j,ks) = d(i,j,ks)
           enddo ! i
         enddo ! j
c
       call pdv (is,ie-1,js,js,ks  ,ks
     &          ,dcopy ,eod )
c
c       subroutine pdv (ibeg,iend,jbeg,jend,kbeg,kend
c     &                ,dlo,eod)
cc
cc    dlo   Mass density (copy)                      for transport step.
cc    eod   Specific energy density e/d (or (e+p)/d] for transport step.
c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
       if(nreq .ne. 0) 
     .    call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
c Finally, do the remaining border.
c
       call pdv (ie  ,ie  ,js  ,js  ,ks  ,ks
     &          ,dcopy ,eod )
c
c      Mark the boundary values out of date.
c
999   continue
      do i = 1,6
       bvstat(i,2) = 0      !  e
       if(leos .ne. 1) bvstat(i,8) = 0 ! T
      enddo
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                 P D V _ D                 \\\\\\\\\\
c
c=======================================================================
c
