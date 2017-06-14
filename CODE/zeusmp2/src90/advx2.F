c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 A D V X 2                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine advx2 (dlo,den
     &                  ,eod,edn
     &                  ,ero,ern
     &                  ,abo,abn
     &                  ,mflx,s1,s2,s3)
c
c    RAF, 2/17/97
c
c  PURPOSE: 
c    Controls the update of density, energy, and momenta
c    from the advection terms in the 2-direction.
c
c  INPUT:
c    dlo         Mass            density prior to update
c    eod         Specific energy density prior to update
c
c  OUTPUT:
c    den         Mass            density    after update
c    edn         Specific energy density    after update
c
c  I/O:
c    s1,s2,s3    Momentum density components (get updated)
c
c  LOCAL:
c    mflx        Mass flux in the 2-direction at zone faces 
c
c  EXTERNALS:
c    BVALV1  , BVALV2  , BVALV3
c    BVALD   , BVALE
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
      integer  :: k1, k2, i, j, k
      real(rl) :: p2
c
      real(rl) :: dlo(in,jn,kn), den(in,jn,kn), mflx(in,jn,kn),
     .            s1 (in,jn,kn), s2 (in,jn,kn), s3  (in,jn,kn),
     .            eod(in,jn,kn), edn(in,jn,kn),
     .            ero(in,jn,kn), ern(in,jn,kn),
     .            abo(in,jn,kn,nspec), abn(in,jn,kn,nspec)
c
      real(rl) :: atwid (ijkn)
      real(rl) :: atwid1 (ijkn), atwid2 (ijkn), atwid3 (ijkn),
     .            atwidj (ijkn)
c
c      Tunable data
c
      data p2 / 0.9 /
c
c-----------------------------------------------------------------------
c
       nseq = nseq + 1        !   nseq indicates the sweep sequence.
       k1   = nint( ( real( ke + ks ) - p2 * nx3z ) * haf )
       k1   = max( k1, ks + 3)
       k2   = ke + ks - k1
       k2   = min( k2, ke - 1 )
c
c-----------------------------------------------------------------------
      if(ldimen .eq. 2) go to 222
c-----------------------------------------------------------------------
c
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
       call bvalv2 (1,0,0,0,0,0,v2 )
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
       call tranx2 (is+1,ie,js+3,je-2,ks+1,k1,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid)


       call momx2  (is+2,ie,js+4,je-3,ks+2,k1,s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj)


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
c
c We need all the density slabs.
c
       call bvald  (0,0,3,3,0,0,dlo)

       if(nspec .gt. 1) call bvalabuns(0,0,3,3,0,0,abo)
c
c
c We need all slabs of eod.
c
       if(xiso .eqv. .false.) call bvale  (0,0,3,3,0,0,eod)
c
       if(lrad .ne. 0) call bvalert(0,0,3,3,0,0,ero)
c
c TRANX2 and MOMX2 together need all 3 velocities at js-2 through je+2.
c
       call bvalv1 (0,0,3,3,0,0,v1 )
       call bvalv2 (0,0,3,3,0,0,v2 )
       call bvalv3 (0,0,3,3,0,0,v3 )


c
c    2) Do middle portion of the interior points, plus some on borders.
c
       call tranx2 (is  ,is  ,js+3,je-2,ks+1,k1,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid)


       call momx2  (is  ,is+1,js+4,je-3,ks+2,k1,s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj)


c
       call tranx2 (is  ,ie  ,js+3,je-2,k1+1,k2,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid)


       call momx2  (is  ,ie  ,js+4,je-3,k1+1,k2,s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj)


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
       call bvald  (0,0,0,0,1,0,dlo)
       call bvalv2 (0,0,0,0,1,0,v2 )


c
c    2) Do last portion of the interior points, plus some on borders.
c
       call tranx2 (is  ,ie  ,js  ,js+2,ks+1,k2,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid)


       call momx2  (is  ,ie  ,js  ,js+3,ks+2,k2,s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj)


c
       call tranx2 (is  ,ie  ,je-1,je  ,ks+1,k2,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid)


       call momx2  (is  ,ie  ,je-2,je  ,ks+2,k2,s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj)


c
       call tranx2 (is  ,ie  ,js  ,je  ,k2+1,ke,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid)


       call momx2  (is  ,ie  ,js  ,je  ,k2+1,ke,s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj)
c
c Mark d and e/d (e) boundary values out of date.
c
       do 40 i=1,7
         bvstat(i,1) = 0  !  d
         if(xiso .eqv. .false.) bvstat(i,2) = 0  !  e or e/d
         if(lrad .ne. 0) bvstat(i,6) = 0  !  er
         if(nspec .gt. 1) bvstat(i,7) = 0  !  abun
40     continue


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
       call tranx2 (is  ,ie  ,js  ,je  ,ks, ks, dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid)


c      write(*,"('ADVX2 before: s2 =',1pd13.5)")s2(3,3,3)
       call momx2  (is  ,ie  ,js  ,je  ,ks, ks+1, s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj)
c      write(*,"('ADVX2 after: s2 =',1pd13.5)")s2(3,3,3)
c
      go to 999
c
c=======================================================================
c     2D TRANSPORT
c=======================================================================
c
222   continue
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
       call bvalv2 (1,0,0,0,0,0,v2 )
c
      if(xiso) then
       if (nseq .eq. 1) then
c
c We need to make a copy of the density, since we skipped pdv.
c
           do j=js-2,je+2
             do i=is-2,ie+2
               dlo(i,j,ks) = den(i,j,ks)
             enddo
           enddo
       endif
c
      endif ! xiso
       call bvald  (1,0,0,0,0,0,dlo)


c
c    2) Do first portion of the interior points.
c
       call tranx2 (is+1,ie,js+3,je-2,ks,ks,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid)
c
       call momx2  (is+2,ie,js+4,je-3,ks,ks,s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj)
c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
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
c
c We need all the density slabs.
c
       call bvald  (0,0,3,3,0,0,dlo)

       if(nspec .gt. 1) call bvalabuns(0,0,3,3,0,0,abo)
c
c
c We need all slabs of eod.
c
       if(xiso .eqv. .false.) call bvale  (0,0,3,3,0,0,eod)
c
       if(lrad .ne. 0) call bvalert(0,0,3,3,0,0,ero)
c
c TRANX2 and MOMX2 together need all 3 velocities at js-2 through je+2.
c
       call bvalv1 (0,0,3,3,0,0,v1 )
       call bvalv2 (0,0,3,3,0,0,v2 )
       call bvalv3 (0,0,3,3,0,0,v3 )
c
c    2) Do middle portion of the interior points, plus some on borders.
c
       call tranx2 (is  ,is  ,js+3,je-2,ks,ks,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid)
c
       call momx2  (is  ,is+1,js+4,je-3,ks,ks,s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj)
c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
c Finally, do the remaining border zones.
c
       call tranx2 (is  ,ie  ,js  ,js+2,ks, ks, dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid)
c
       call momx2  (is  ,ie  ,js  ,js+3,ks,ks,s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj)
c
       call tranx2 (is  ,ie  ,je-1,je  ,ks,ks,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid)
c
       call momx2  (is  ,ie  ,je-2,je  ,ks, ks, s1,s2,s3,mflx,
     .              atwid1,atwid2,atwid3,atwidj)
c
c Mark d and e/d (e) boundary values out of date.
c
999    continue
c
       do i=1,7
         bvstat(i,1) = 0  !  d
         if(xiso .eqv. .false.) bvstat(i,2) = 0  !  e or e/d
         bvstat(i,6) = 0  !  er
         if(nspec .gt. 1) bvstat(i,7) = 0  !  abun
       enddo
c
       return
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                 A D V X 2                 \\\\\\\\\\
c
c=======================================================================
c
c
