c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 A D V X 1                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c    
c=======================================================================
c
       subroutine advx1 (dlo,den
     &                  ,eod,edn
     &                  ,ero,ern
     .                  ,abo,abn
     &                  ,mflx,s1,s2,s3)
c
c    RAF, 2/17/97
c
c  PURPOSE: 
c    Controls the update of density, energy, and momenta
c    from the advection terms in the 1-direction.
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
c    mflx        Mass flux in the 1-direction at zone faces 
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
      integer  :: k1, k2, i, j, k
c
      real(rl) :: p1
c
      real(rl) :: dlo(in,jn,kn), den(in,jn,kn), mflx(in,jn,kn),
     .            s1 (in,jn,kn), s2 (in,jn,kn), s3  (in,jn,kn),
     .            eod(in,jn,kn), edn(in,jn,kn),
     .            ero(in,jn,kn), ern(in,jn,kn),
     .            abo(in,jn,kn,nspec), abn(in,jn,kn,nspec)
c
      real(rl) :: atwid (ijkn),  mflux (ijkn),
     .            dtwid (ijkn),  dd    (ijkn),
     .            etwid (ijkn),  deod  (ijkn)
c
      real(rl) :: atwid1(ijkn)
      real(rl) :: vtwid ( ijkn ), sflx  ( ijkn ),
     .            dq ( ijkn )
c
c      Tunable data
c
      data p1 / 0.9 /
c
c-----------------------------------------------------------------------
c
       nseq = nseq + 1        !   nseq indicates the sweep sequence.
       k1   = int( real( nx3z ) * p1 ) + ks   !  The lion's share...
       k1   = min( k1, ke - 2 )               !  but not too much!
       k2   = ( ke + k1 ) / 2                 !  Half the remainder.

      if(ldimen .eq. 2) go to 222
      if(ldimen .eq. 1) go to 111

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
c
c TRANX1 and MOMX1 together need all 3 velocities at is-2 through ie+2.
c
       call bvalv1 (3,3,0,0,0,0,v1 )
       call bvalv2 (3,3,0,0,0,0,v2 )
       call bvalv3 (3,3,0,0,0,0,v3 )
c
c We need all slabs of er/d.
c
      if(lrad .ne. 0) call bvalert(3,3,0,0,0,0,ero)
      if(xiso .eqv. .false.) then
c
c We need all slabs of eod.
c
       call bvale  (3,3,0,0,0,0,eod)
c

      else ! xiso
c
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
c
c We need all the density slabs.
c
       call bvald  (3,3,0,0,0,0,dlo)

       if(nspec .gt. 1) call bvalabuns(3,3,0,0,0,0,abo)

c
c    2) Do first portion of the interior points.
c
       call tranx1 (is+3,ie-2,js+1,je  ,ks+1,k1,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,dtwid,etwid,mflux,dd,deod)


       call momx1  (is+4,ie-3,js+2,je  ,ks+2,k1,s1,s2,s3,mflx,
     .              atwid1,vtwid,sflx,dq)


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
       call bvalv1 (0,0,1,0,0,0,v1 )


c
c    2) Do second portion of the interior points, plus some on borders.
c
       call tranx1 (is  ,is+2,js+1,je  ,ks+1,k1,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,dtwid,etwid,mflux,dd,deod)


       call momx1  (is  ,is+3,js+2,je  ,ks+2,k1,s1,s2,s3,mflx,
     .              atwid1,vtwid,sflx,dq)


c
       call tranx1 (ie-1,ie  ,js+1,je  ,ks+1,k1,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,dtwid,etwid,mflux,dd,deod)


       call momx1  (ie-2,ie  ,js+2,je  ,ks+2,k1,s1,s2,s3,mflx,
     .              atwid1,vtwid,sflx,dq)


c
       call tranx1 (is  ,ie  ,js+1,je  ,k1+1,k2,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,dtwid,etwid,mflux,dd,deod)


       call momx1  (is  ,ie  ,js+2,je  ,k1+1,k2,s1,s2,s3,mflx,
     .              atwid1,vtwid,sflx,dq)


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
       call bvalv1 (0,0,0,0,1,0,v1 )
c       call bvalv3 (0,0,0,0,0,1,v3 )     !   Did this already in pdv.


c
c    2) Do last portion of the interior points, plus some on borders.
c
       call tranx1 (is  ,ie  ,js  ,js  ,ks+1,k2,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,dtwid,etwid,mflux,dd,deod)


       call momx1  (is  ,ie  ,js  ,js+1,ks+2,k2,s1,s2,s3,mflx,
     .              atwid1,vtwid,sflx,dq)


c
       call tranx1 (is  ,ie  ,js  ,je  ,k2+1,ke  ,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,dtwid,etwid,mflux,dd,deod)


       call momx1  (is  ,ie  ,js  ,je  ,k2+1,ke  ,s1,s2,s3,mflx,
     .              atwid1,vtwid,sflx,dq)
c
c Mark d and e/d (e) boundary values out of date.
c
       do 40 i=1,7
         bvstat(i,1) = 0  !  d
         if(xiso .eqv. .false.) bvstat(i,2) = 0  !  e or e/d
         if(lrad .ne. 0) bvstat(i,6) = 0  !  er
         if(nspec .gt. 1) bvstat(i,7) = 0
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
       call tranx1 (is  ,ie  ,js  ,je  ,ks, ks, dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,dtwid,etwid,mflux,dd,deod)


       call momx1  (is  ,ie  ,js  ,je  ,ks, ks+1, s1,s2,s3,mflx,
     .              atwid1,vtwid,sflx,dq)
c
       go to 999

c=======================================================================
c      2D TRANSPORT
c=======================================================================
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
c
c TRANX1 and MOMX1 together need all 3 velocities at is-2 through ie+2.
c
      call bvalv1 (3,3,0,0,0,0,v1 )
      call bvalv2 (3,3,0,0,0,0,v2 )
      call bvalv3 (3,3,0,0,0,0,v3 )
c
c We need all slabs of er/d.
c
      if(lrad .ne. 0) call bvalert(3,3,0,0,0,0,ero)
      if(xiso .eqv. .false.) then
c
c We need all slabs of eod.
c
       call bvale  (3,3,0,0,0,0,eod)
c

      else ! xiso
c
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
c
c We need all the density slabs.
c
       call bvald  (3,3,0,0,0,0,dlo)

       if(nspec .gt. 1) call bvalabuns(3,3,0,0,0,0,abo)

c
c    2) Do first portion of the interior points.
c
       call tranx1 (is+3,ie-2,js+1,je  ,ks,ks,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,dtwid,etwid,mflux,dd,deod)


       call momx1  (is+4,ie-3,js+2,je  ,ks,ks,s1,s2,s3,mflx,
     .              atwid1,vtwid,sflx,dq)


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
       call bvalv1 (0,0,1,0,0,0,v1 )


c
c    2) Do second portion of the interior points, plus some on borders.
c
       call tranx1 (is  ,is+2,js+1,je  ,ks,ks,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,dtwid,etwid,mflux,dd,deod)


       call momx1  (is  ,is+3,js+2,je  ,ks,ks,s1,s2,s3,mflx,
     .              atwid1,vtwid,sflx,dq)


c
       call tranx1 (ie-1,ie  ,js+1,je  ,ks,ks,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,dtwid,etwid,mflux,dd,deod)


       call momx1  (ie-2,ie  ,js+2,je  ,ks,ks,s1,s2,s3,mflx,
     .              atwid1,vtwid,sflx,dq)
c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c
c Finally, do the remaining border zones.
c
       call tranx1 (is  ,ie  ,js  ,js  ,ks, ks, dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
c
       call momx1  (is  ,ie  ,js  ,js+1,ks, ks, s1,s2,s3,mflx,
     .              atwid1,vtwid,sflx,dq)
c
c Mark d and e/d (e) boundary values out of date.
c
       do i=1,7
         bvstat(i,1) = 0  !  d
         if(xiso .eqv. .false.) bvstat(i,2) = 0  !  e or e/d
         if(lrad .ne. 0) bvstat(i,6) = 0  !  er
         if(nspec .gt. 1) bvstat(i,7) = 0
       enddo
c
      go to 999
c=======================================================================
c      1D TRANSPORT
c=======================================================================
111   continue

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
c
c TRANX1 and MOMX1 together need all 3 velocities at is-2 through ie+2.
c
      call bvalv1 (3,3,0,0,0,0,v1 )
      call bvalv2 (3,3,0,0,0,0,v2 )
      call bvalv3 (3,3,0,0,0,0,v3 )
c
c We need all slabs of er/d.
c
      if(lrad .ne. 0) call bvalert(3,3,0,0,0,0,ero)
      if(xiso .eqv. .false.) then
c
c We need all slabs of eod.
c
       call bvale  (3,3,0,0,0,0,eod)
c

      else ! xiso
c
       if (nseq .eq. 1) then
c
c We need to make a copy of the density, since we skipped pdv.
c
           do j=js,js
             do i=is-2,ie+2
               dlo(i,j,ks) = den(i,j,ks)
             enddo
           enddo
       endif
c
      endif ! xiso
c
c We need all the density slabs.
c
       call bvald  (3,3,0,0,0,0,dlo)

       if(nspec .gt. 1) call bvalabuns(3,3,0,0,0,0,abo)
c
c    2) Do first portion of the interior points.
c
       call tranx1 (is+3,ie-2,js,js  ,ks,ks,dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,dtwid,etwid,mflux,dd,deod)


       call momx1  (is+4,ie-3,js,js  ,ks,ks,s1,s2,s3,mflx,
     .              atwid1,vtwid,sflx,dq)


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
       call tranx1 (is  ,is+2,js  ,js,ks, ks, dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
c
       call momx1  (is  ,is+3,js  ,js,ks, ks, s1,s2,s3,mflx,
     .              atwid1,vtwid,sflx,dq)
c
       call tranx1 (ie-1,ie  ,js  ,js,ks, ks, dlo,den
     &             ,eod,edn
     &             ,ero,ern
     &             ,abo,abn
     &             ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
c
       call momx1  (ie-2,ie  ,js  ,js,ks, ks, s1,s2,s3,mflx,
     .              atwid1,vtwid,sflx,dq)
c
c Mark d and e/d (e) boundary values out of date.
c
       do i=1,7
         bvstat(i,1) = 0  !  d
         if(xiso .eqv. .false.) bvstat(i,2) = 0  !  e or e/d
         if(lrad .ne. 0) bvstat(i,6) = 0  !  er
         if(nspec .gt. 1) bvstat(7,1) = 0
       enddo
c
999    return
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                 A D V X 1                 \\\\\\\\\\
c
c=======================================================================
c
c
