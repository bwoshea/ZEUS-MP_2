c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 A D V X 3                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine advx3 (dlo,den
     &                  ,eod,edn
     &                  ,ero,ern
     &                  ,abo,abn
     &                  ,mflx,s1,s2,s3)
c
c    RAF, 2/19/97
c
c  PURPOSE: 
c    Controls the update of density, energy, and momenta
c    from the advection terms in the 3-direction.
c
c  INPUT:
c    dlo         Mass            density prior to update
c    eod         Specific energy density prior to update
c    ero         Specific radiation energy density prior to update
c
c  OUTPUT:
c    den         Mass            density    after update
c    edn         Specific energy density    after update
c    ern         Specific radiation energy density    after update
c
c  I/O:
c    s1,s2,s3    Momentum density components (get updated)
c
c  LOCAL:
c    mflx        Mass flux in the 3-direction at zone faces 
c
c  EXTERNALS:
c    BVALV1  , BVALV2  , BVALV3
c    BVALD   , BVALE   , BVALER
c
c-----------------------------------------------------------------------
      use real_prec
      use param
      use config
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
      integer  :: k1,k2,i
      integer  :: kbeg, kend, krange, kblocks, kskip, ktlb, kblk, kpage
      integer  :: j, k
      real(rl) :: p3
c
      real(rl) :: dlo(in,jn,kn), den(in,jn,kn), mflx(in,jn,kn),
     .            s1 (in,jn,kn), s2 (in,jn,kn), s3  (in,jn,kn),
     .            eod(in,jn,kn), edn(in,jn,kn),
     .            ero(in,jn,kn), ern(in,jn,kn),
     .            abo(in,jn,kn,nspec), abn(in,jn,kn,nspec)
c
      real(rl) :: atwid (ijkn)
      real(rl) :: mflux (ijkn,1),
     .            dtwid (ijkn,1), dd  (ijkn,1),
     .            etwid (ijkn,1), deod(ijkn,1)
c
      real(rl) :: atwid1 (ijkn), atwid2 (ijkn), atwid3 (ijkn),
     .            atwidj1(ijkn), atwidj2(ijkn), atwidj3(ijkn)
c
      real(rl) :: sflx  (ijkn,1), dq   (ijkn,1)
c
c      Tunable data
c
      data p3 / 0.9 /   ! Fraction of interior points to do 3rd stage
c
c Set blocking factor for 3rd stage k-loops.  The number of iterations
c in a block should be small enough so that all the data fits on the
c number of pages that the TLB can hold.  Thus,
c
c iterations = <TLB entries * page size> / <2-D arrays * word size>
c
c In tranx3 and momx3, 11 different 2-D arrays of data are used.  
c Assume 16kB pages and 8B words.  The number of iterations per
c block is forced to be at least 5.
c
C#ifdef TLB
C       data ktlb / TLB /  ! Get TLB size in MB from cpp -DTLB=n
       data ktlb / 128 /  ! This is the default for compilation on 
C                           SGI O2K
C#else /* TLB */
C       data ktlb / 9999 /    ! Assume many TLB entries
C#endif /* TLB */
C#ifdef BLK_MIN
C       data kblk / BLK_MIN /  ! Get min block size from cpp -DBLK_MIN=n
C#else /* BLK_MIN */
       data kblk / 5 /        ! Default min block size is this many
C#endif /* BLK_MIN */
C#ifdef PAGE
C       data kpage / PAGE /    ! Get PAGE size in bytes from cpp -DPAGE=n
C#else /* PAGE */
       data kpage / 16384 /   ! Assume this many bytes per page
C#endif /* PAGE */
c
c-----------------------------------------------------------------------
c
c Divide up the work.  Since we must have
c
c   ks+4 < k1-1
c   k1   < k2-1
c   k2   < ke-3
c
c   ke - ks .ge. 12   --- this is the smallest allowable k range.
c
       nseq = nseq + 1        !   nseq indicates the sweep sequence.
       k2   = ke - int( p3 * nx3z )
       k2   = min( k2, ke - 3 )
       k2   = max( k2, ks + 6 )
       k1   = ( k2 + ks ) / 2
       k1   = max( k1, ks + 5 )


c......................................................................
c
c i boundaries
c
c    1) Post sends and receives. 
c       By exchanging the i, j, and k boundary data in three separate 
c       steps, we ensure that the corner and edge data are correctly
c       updated.
c
       nreq = 0
       nsub = nsub + 1
       call bvalv3 (1,0,0,0,0,0,v3 )
c
      if(xiso) then
       if (nseq .eq. 1) then
c
c We need to make a copy of the density, since we skipped pdv.
c
         do 30 k=ks-2,ke+2
           do 20 j=js-2,je+2
             do 10 i=is-2,ie+2
               dlo(i,j,k) = den(i,j,k)
10           continue
20         continue
30       continue
       endif
c
      endif ! xiso
       call bvald  (1,0,0,0,0,0,dlo)


c
c    2) Do first portion of the interior points.
c
       call tranx3 (is+1,ie,js+1,je,ks+3,k1,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod)


       call momx3  (is+2,ie,js+2,je,ks+4,k1-1,s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj1,atwidj2,atwidj3,
     .              sflx,dq)


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
c    1) Post sends and receives.
c
       nreq = 0
       nsub = nsub + 1
       call bvald  (0,0,1,0,0,0,dlo)
       call bvalv3 (0,0,1,0,0,0,v3 )


c
c    2) Do second portion of the interior points, plus some on borders.
c
       call tranx3 (is  ,is  ,js+1,je,ks+3,k1,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod)


       call momx3  (is  ,is+1,js+2,je,ks+4,k1-1,s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj1,atwidj2,atwidj3,
     .              sflx,dq)


c
       call tranx3 (is  ,ie  ,js+1,je,k1+1,k2,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod)


       call momx3  (is  ,ie  ,js+2,je,k1  ,k2-1,s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj1,atwidj2,atwidj3,
     .              sflx,dq)


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
c    1) Post sends and receives.
c
       nreq = 0
       nsub = nsub + 1
c
c We need all the density slabs.
c
       call bvald  (0,0,0,0,3,3,dlo)

       if(nspec .gt. 1) call bvalabuns(0,0,0,0,3,3,abo)
c
c
c We need all slabs of eod.
c
       if(xiso .eqv. .false.) call bvale  (0,0,0,0,3,3,eod)
c
       if(lrad .ne. 0) call bvalert(0,0,0,0,3,3,ero)
c
c TRANX3 and MOMX3 together need all 3 velocities at ks-2 through ke+2.
c
       call bvalv1 (0,0,0,0,3,3,v1 )
       call bvalv2 (0,0,0,0,3,3,v2 )
       call bvalv3 (0,0,0,0,3,3,v3 )


c
c    2) Do last portion of the interior points, plus some on borders.
c
       call tranx3 (is  ,ie  ,js  ,js  ,ks+3,k2,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod)


       call momx3  (is  ,ie  ,js  ,js+1,ks+4,k2-1,s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj1,atwidj2,atwidj3,
     .              sflx,dq)


c
c Block the k loop to reduce TLB misses; the k ranges for the earlier
c stages above should be small enough already.
c
       kskip   = max( ktlb * kpage / (11*8*in*jn), kblk)
       krange = ke-2 - (k2+1) + 1
       kblocks = krange / kskip
       if ( mod(krange,kskip) .eq. 0) kblocks = max(kblocks-1,0)
c
       do 35 kbeg = k2+1, k2+1 + kblocks*kskip, kskip 
         kend = min( kbeg + kskip - 1, ke-2 )
         call tranx3 (is  ,ie  ,js  ,je  ,kbeg,kend,dlo,den
     &               ,eod,edn
     &               ,ero,ern
     &               ,abo,abn
     &             ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod)


       call momx3  (is  ,ie  ,js  ,je  ,kbeg-1,kend-1,s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj1,atwidj2,atwidj3,
     .              sflx,dq)


35     continue
c
c Mark d and e/d (e) boundary values out of date.
c
       do 40 i=1,6
         bvstat(i,1) = 0  !  d
         if(xiso .eqv. .false.) bvstat(i,2) = 0  !  e or e/d
         if(lrad .ne. 0) bvstat(i,6) = 0  !  er
         if(nspec .gt. 1) bvstat(i,7) = 0  !  abun
40     continue
c

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
       call tranx3 (is  ,ie  ,js  ,je  ,ks, ks+2, dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod)


       call momx3  (is  ,ie  ,js  ,je  ,ks, ks+3, s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj1,atwidj2,atwidj3,
     .              sflx,dq)


c
       call tranx3 (is  ,ie  ,js  ,je  ,ke-1, ke, dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod)


       call momx3  (is  ,ie  ,js  ,je  ,ke-2, ke, s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj1,atwidj2,atwidj3,
     .              sflx,dq)


c
       return
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                 A D V X 3                 \\\\\\\\\\
c
c=======================================================================
c
c
