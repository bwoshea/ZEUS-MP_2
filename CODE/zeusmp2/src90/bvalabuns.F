c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////             B V A L A B U N S             \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c=======================================================================
c
       subroutine bvalabuns( rl1, ru1, rl2, ru2, rl3, ru3, ab )
c
c      a clone of BVALD generated for the ABUN abundance array
c
      use real_prec
      use config
      use param
      use root
      use grid
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      real(rl) :: ab(in,jn,kn,nspec)
c
      integer  :: i,j,k,l1,l2,l3,u1,u2,u3,
     .            rl1,rl2,rl3,ru1,ru2,ru3,
     .           ls,ll,lu,us,ul,uu,n
c
c-----------------------------------------------------------------------
c------------------------  I - B O U N D A R Y  ------------------------
c-----------------------------------------------------------------------
c
       l1 = rl1 - bvstat(1,7)
       u1 = ru1 - bvstat(2,7)
c
c      Inner i boundary.
c
#ifdef MPI_USED
c
c Count slabs, compute positioning indices.
c
       ls = max (l1-1,1)  ! number of slabs to send/receive
       ll = min (l1,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l1 not 2
       us = max (u1-1,1)  ! number of slabs to send/receive
       ul = min (u1,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if u1 not 2
c
c Post a receive for a slab of data from the interior of the 
c neighboring tile to fill my ghost zones.  Initiate a send 
c to pass a slab of my interior data for my neighbor's ghost zones.
c 
       if (niis(1).eq.0 .or. niis(1).eq.4) then
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_IRECV(ab(is-ll+i-1, 1, 1, 1), 1, iab_slice, n1m
     &                     ,11100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,7) = rl1
         endif
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_ISEND(ab(is+uu+i-1, 1, 1, 1), 1, iab_slice, n1m
     &                     ,11200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,7) = ru1
         endif
       else
         if (l1 .gt. 0) then
          do n = 1, nspec
           do k=ks-1,ke+1
cdir$ ivdep
            do j=js-1,je+1
             if ( abs(niib(j,k)) .eq. 1) then
               ab(is-1,j,k,n) = ab(is  ,j,k,n)
               ab(is-2,j,k,n) = ab(is+1,j,k,n)
             endif
             if (niib(j,k) .eq. 2) then
               ab(is-1,j,k,n) = ab(is  ,j,k,n)
               ab(is-2,j,k,n) = ab(is-1,j,k,n)
             endif
             if (niib(j,k) .eq. 3) then
               ab(is-1,j,k,n) = abiib(j,k,1,n)
               ab(is-2,j,k,n) = abiib(j,k,2,n)
             endif
             if (niib(j,k) .eq. 5) then
               ab(is-1,j,k,n) = ab(is  ,j,k,n)
               ab(is-2,j,k,n) = ab(is+1,j,k,n)
             endif
            enddo 
           enddo
         enddo ! n
         bvstat(1,7) = rl1
         endif
       endif
c
c      Outer i boundary.
c
       if (nois(1).eq.0 .or. nois(1).eq.4) then
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_IRECV(ab(ie+i+uu, 1, 1, 1), 1, iab_slice, n1p
     &                     ,11200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,7) = ru1
         endif
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_ISEND(ab(ie+i-ll, 1, 1, 1), 1, iab_slice, n1p
     &                     ,11100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,7) = rl1
         endif
       else
         if (u1 .gt. 0) then
         do n = 1, nspec
          do k=ks-1,ke+1
cdir$ ivdep
           do j=js-1,je+1
             if ( abs(noib(j,k)) .eq. 1) then
               ab(ie+1,j,k,n) = ab(ie,j,k,n)
               ab(ie+2,j,k,n) = ab(ie-1,j,k,n)
             endif
             if (noib(j,k) .eq. 2) then
               ab(ie+1,j,k,n) = ab(ie,j,k,n)
               ab(ie+2,j,k,n) = ab(ie+1,j,k,n)
             endif
             if (noib(j,k) .eq. 3) then
               ab(ie+1,j,k,n) = aboib(j,k,1,n)
               ab(ie+2,j,k,n) = aboib(j,k,2,n)
             endif
             if (noib(j,k) .eq. 5) then
               ab(ie+1,j,k,n) = ab(ie  ,j,k,n)
               ab(ie+2,j,k,n) = ab(ie-1,j,k,n)
             endif
            enddo 
           enddo
         enddo
         bvstat(2,7) = ru1
         endif
       endif
#endif /* MPI_USED */
#ifndef MPI_USED
         if (l1 .gt. 0) then
          do n = 1, nspec
           do k=ks-1,ke+1
cdir$ ivdep
            do j=js-1,je+1
             if ( abs(niib(j,k)) .eq. 1) then
               ab(is-1,j,k,n) = ab(is  ,j,k,n)
               ab(is-2,j,k,n) = ab(is+1,j,k,n)
             endif
             if (niib(j,k) .eq. 2) then
               ab(is-1,j,k,n) = ab(is  ,j,k,n)
               ab(is-2,j,k,n) = ab(is-1,j,k,n)
             endif
             if (niib(j,k) .eq. 3) then
               ab(is-1,j,k,n) = abiib(j,k,1,n)
               ab(is-2,j,k,n) = abiib(j,k,2,n)
             endif
             if (niib(j,k) .eq. 4) then
               ab(is-1,j,k,n) = ab(ie  ,j,k,n)
               ab(is-2,j,k,n) = ab(ie-1,j,k,n)
             endif
             if (niib(j,k) .eq. 5) then
               ab(is-1,j,k,n) = ab(is  ,j,k,n)
               ab(is-2,j,k,n) = ab(is+1,j,k,n)
             endif
            enddo 
           enddo
         enddo
         bvstat(1,7) = rl1
         endif
c
c      Outer i boundary.
c
         if (u1 .gt. 0) then
          do n = 1, nspec
           do k=ks-1,ke+1
cdir$ ivdep
            do j=js-1,je+1
             if ( abs(noib(j,k)) .eq. 1) then
               ab(ie+1,j,k,n) = ab(ie,j,k,n)
               ab(ie+2,j,k,n) = ab(ie-1,j,k,n)
             endif
             if (noib(j,k) .eq. 2) then
               ab(ie+1,j,k,n) = ab(ie,j,k,n)
               ab(ie+2,j,k,n) = ab(ie+1,j,k,n)
             endif
             if (noib(j,k) .eq. 3) then
               ab(ie+1,j,k,n) = aboib(j,k,1,n)
               ab(ie+2,j,k,n) = aboib(j,k,2,n)
             endif
             if (noib(j,k) .eq. 4) then
               ab(ie+1,j,k,n) = ab(is  ,j,k,n)
               ab(ie+2,j,k,n) = ab(is+1,j,k,n)
             endif
             if (noib(j,k) .eq. 5) then
               ab(ie+1,j,k,n) = ab(ie  ,j,k,n)
               ab(ie+2,j,k,n) = ab(ie-1,j,k,n)
             endif
            enddo 
           enddo
         enddo
         bvstat(2,7) = ru1
         endif
#endif /* NO MPI */
c
c-----------------------------------------------------------------------
c------------------------  J - B O U N D A R Y  ------------------------
c-----------------------------------------------------------------------
c
       l2 = rl2 - bvstat(3,7)
       u2 = ru2 - bvstat(4,7)
c
c      Inner j boundary.
c
#ifdef MPI_USED
c
c Count slabs, compute positioning indices.
c
       ls = max (l2-1,1)  ! number of inner slabs to send/receive
       ll = min (l2,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l2 not 2
       us = max (u2-1,1)  ! number of outer slabs to send/receive
       ul = min (u2,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l2 not 2
c
       if (nijs(1).eq.0 .or. nijs(1).eq.4) then
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_IRECV(ab( 1,js-ll+j-1, 1, 1), 1, jab_slice, n2m
     &                     ,11300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,7) = rl2
         endif
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_ISEND(ab ( 1,js+uu+j-1, 1, 1), 1, jab_slice, n2m
     &                     ,11400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,7) = ru2
         endif
       else
         if (l2 .gt. 0) then
          do n = 1, nspec
           do k=ks-1,ke+1
cdir$ ivdep
            do i=is-1,ie+1
             if ( abs(nijb(i,k)) .eq. 1) then
               ab(i,js-1,k,n) = ab(i,js  ,k,n)
               ab(i,js-2,k,n) = ab(i,js+1,k,n)
             endif
             if (nijb(i,k) .eq. 2) then
               ab(i,js-1,k,n) = ab(i,js  ,k,n)
               ab(i,js-2,k,n) = ab(i,js-1,k,n)
             endif
             if (nijb(i,k) .eq. 3) then
               ab(i,js-1,k,n) = abijb(i,k,1,n)
               ab(i,js-2,k,n) = abijb(i,k,2,n)
             endif
             if (nijb(i,k) .eq. 5) then
               ab(i,js-1,k,n) = ab(i,js  ,k,n)
               ab(i,js-2,k,n) = ab(i,js+1,k,n)
             endif
            enddo
           enddo
         enddo
         bvstat(3,7) = rl2
         endif
       endif
c
c      Outer j boundary.
c
       if (nojs(1).eq.0 .or. nojs(1).eq.4) then
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_IRECV(ab ( 1,je+j+uu, 1, 1), 1, jab_slice, n2p
     &                     ,11400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,7) = ru2
         endif
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_ISEND(ab( 1,je+j-ll,   1, 1), 1, jab_slice, n2p
     &                     ,11300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,7) = rl2
         endif
       else
         if (u2 .gt. 0) then
          do n = 1, nspec
           do k=ks-1,ke+1
cdir$ ivdep
            do i=is-1,ie+1
             if ( abs(nojb(i,k)) .eq. 1) then
               ab(i,je+1,k,n) = ab(i,je  ,k,n)
               ab(i,je+2,k,n) = ab(i,je-1,k,n)
             endif
             if (nojb(i,k) .eq. 2) then
               ab(i,je+1,k,n) = ab(i,je  ,k,n)
               ab(i,je+2,k,n) = ab(i,je+1,k,n)
             endif
             if (nojb(i,k) .eq. 3) then
               ab(i,je+1,k,n) = abojb(i,k,1,n)
               ab(i,je+2,k,n) = abojb(i,k,2,n)
             endif
             if (nojb(i,k) .eq. 5) then
               ab(i,je+1,k,n) = ab(i,je  ,k,n)
               ab(i,je+2,k,n) = ab(i,je-1,k,n)
             endif
            enddo
           enddo
         enddo
         bvstat(4,7) = ru2
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED 
         if (l2 .gt. 0) then
          do n = 1, nspec
           do k=ks-1,ke+1
cdir$ ivdep
            do i=is-1,ie+1
             if ( abs(nijb(i,k)) .eq. 1) then
               ab(i,js-1,k,n) = ab(i,js  ,k,n)
               ab(i,js-2,k,n) = ab(i,js+1,k,n)
             endif
             if (nijb(i,k) .eq. 2) then
               ab(i,js-1,k,n) = ab(i,js  ,k,n)
               ab(i,js-2,k,n) = ab(i,js-1,k,n)
             endif
             if (nijb(i,k) .eq. 3) then
               ab(i,js-1,k,n) = abijb(i,k,1,n)
               ab(i,js-2,k,n) = abijb(i,k,2,n)
             endif
             if (nijb(i,k) .eq. 4) then
               ab(i,js-1,k,n) = ab(i,je  ,k,n)
               ab(i,js-2,k,n) = ab(i,je-1,k,n)
             endif
             if (nijb(i,k) .eq. 5) then
               ab(i,js-1,k,n) = ab(i,js  ,k,n)
               ab(i,js-2,k,n) = ab(i,js+1,k,n)
             endif
            enddo
           enddo
         enddo
         bvstat(3,7) = rl2
         endif
c
c      Outer j boundary.
c
         if (u2 .gt. 0) then
          do n = 1, nspec
           do k=ks-1,ke+1
cdir$ ivdep
            do i=is-1,ie+1
             if ( abs(nojb(i,k)) .eq. 1) then
               ab(i,je+1,k,n) = ab(i,je  ,k,n)
               ab(i,je+2,k,n) = ab(i,je-1,k,n)
             endif
             if (nojb(i,k) .eq. 2) then
               ab(i,je+1,k,n) = ab(i,je  ,k,n)
               ab(i,je+2,k,n) = ab(i,je+1,k,n)
             endif
             if (nojb(i,k) .eq. 3) then
               ab(i,je+1,k,n) = abojb(i,k,1,n)
               ab(i,je+2,k,n) = abojb(i,k,2,n)
             endif
             if (nojb(i,k) .eq. 4) then
               ab(i,je+1,k,n) = ab(i,js  ,k,n)
               ab(i,je+2,k,n) = ab(i,js+1,k,n)
             endif
             if (nojb(i,k) .eq. 5) then
               ab(i,je+1,k,n) = ab(i,je  ,k,n)
               ab(i,je+2,k,n) = ab(i,je-1,k,n)
             endif
            enddo
           enddo
         enddo
         bvstat(4,7) = ru2
         endif
#endif /* NO MPI */
c
c-----------------------------------------------------------------------
c------------------------  K - B O U N D A R Y  ------------------------
c-----------------------------------------------------------------------
c
       l3 = rl3 - bvstat(5,7)
       u3 = ru3 - bvstat(6,7)
c
c      Inner k boundary.
c
#ifdef MPI_USED
c
c Count slabs, compute positioning indices.
c
       ls = max (l3-1,1)  ! number of inner slabs to send/receive
       ll = min (l3,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l3 not 2
       us = max (u3-1,1)  ! number of outer slabs to send/receive
       ul = min (u3,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l3 not 2
c
       if (niks(1).eq.0 .or. niks(1).eq.4) then
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_IRECV(ab( 1, 1,ks-ll+k-1, 1), 1, kab_slice, n3m
     &                     ,11500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,7) = rl3
         endif
         if (u3 .gt. 0) then
         do k=1,us
             nreq = nreq + 1
             call MPI_ISEND(ab( 1, 1,ks+uu+k-1, 1), 1, kab_slice, n3m
     &                     ,11600+25*k+nsub, comm3d, req(nreq), ierr)
         enddo
         bvstat(6,7) = ru3
         endif
       else
         if (l3 .gt. 0) then
          do n = 1, nspec
           do j=js-1,je+1
cdir$ ivdep
            do i=is-1,ie+1
             if ( abs(nikb(i,j)) .eq. 1) then
               ab(i,j,ks-1,n) = ab(i,j,ks  ,n)
               ab(i,j,ks-2,n) = ab(i,j,ks+1,n)
             endif
             if (nikb(i,j) .eq. 2) then
               ab(i,j,ks-1,n) = ab(i,j,ks  ,n)
               ab(i,j,ks-2,n) = ab(i,j,ks-1,n)
             endif
             if (nikb(i,j) .eq. 3) then
               ab(i,j,ks-1,n) = abikb(i,j,1,n)
               ab(i,j,ks-2,n) = abikb(i,j,2,n)
             endif
             if (nikb(i,j) .eq. 5) then
               ab(i,j,ks-1,n) = ab(i,j,ks  ,n)
               ab(i,j,ks-2,n) = ab(i,j,ks+1,n)
             endif
            enddo
           enddo
         enddo
         bvstat(5,7) = rl3
         endif
       endif
c
c      Outer k boundary.
c
       if (noks(1).eq.0 .or. noks(1).eq.4) then
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_IRECV(ab( 1, 1,ke+k+uu, 1), 1, kab_slice, n3p
     &                     ,11600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,7) = ru3
         endif
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_ISEND(ab( 1, 1,ke+k-ll, 1), 1, kab_slice, n3p
     &                     ,11500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,7) = rl3
         endif
       else
         if (u3 .gt. 0) then
          do n = 1, nspec
           do j=js-1,je+1
cdir$ ivdep
            do i=is-1,ie+1
             if ( abs(nokb(i,j)) .eq. 1) then
               ab(i,j,ke+1,n) = ab(i,j,ke  ,n)
               ab(i,j,ke+2,n) = ab(i,j,ke-1,n)
             endif
             if (nokb(i,j) .eq. 2) then
               ab(i,j,ke+1,n) = ab(i,j,ke  ,n)
               ab(i,j,ke+2,n) = ab(i,j,ke+1,n)
             endif
             if (nokb(i,j) .eq. 3) then
               ab(i,j,ke+1,n) = abokb(i,j,1,n)
               ab(i,j,ke+2,n) = abokb(i,j,2,n)
             endif
             if (nokb(i,j) .eq. 5) then
               ab(i,j,ke+1,n) = ab(i,j,ke  ,n)
               ab(i,j,ke+2,n) = ab(i,j,ke-1,n)
             endif
            enddo 
           enddo
         enddo
         bvstat(6,7) = ru3
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED 
         if (l3 .gt. 0) then
          do n = 1, nspec
           do j=js-1,je+1
cdir$ ivdep
            do i=is-1,ie+1
             if ( abs(nikb(i,j)) .eq. 1) then
               ab(i,j,ks-1,n) = ab(i,j,ks  ,n)
               ab(i,j,ks-2,n) = ab(i,j,ks+1,n)
             endif
             if (nikb(i,j) .eq. 2) then
               ab(i,j,ks-1,n) = ab(i,j,ks  ,n)
               ab(i,j,ks-2,n) = ab(i,j,ks-1,n)
             endif
             if (nikb(i,j) .eq. 3) then
               ab(i,j,ks-1,n) = abikb(i,j,1,n)
               ab(i,j,ks-2,n) = abikb(i,j,2,n)
             endif
             if (nikb(i,j) .eq. 4) then
               ab(i,j,ks-1,n) = ab(i,j,ke  ,n)
               ab(i,j,ks-2,n) = ab(i,j,ke-1,n)
             endif
             if (nikb(i,j) .eq. 5) then
               ab(i,j,ks-1,n) = ab(i,j,ks  ,n)
               ab(i,j,ks-2,n) = ab(i,j,ks+1,n)
             endif
            enddo
           enddo
         enddo
         bvstat(5,7) = rl3
         endif
c
c      Outer k boundary.
c
         if (u3 .gt. 0) then
          do n = 1, nspec
           do j=js-1,je+1
cdir$ ivdep
            do i=is-1,ie+1
             if ( abs(nokb(i,j)) .eq. 1) then
               ab(i,j,ke+1,n) = ab(i,j,ke  ,n)
               ab(i,j,ke+2,n) = ab(i,j,ke-1,n)
             endif
             if (nokb(i,j) .eq. 2) then
               ab(i,j,ke+1,n) = ab(i,j,ke  ,n)
               ab(i,j,ke+2,n) = ab(i,j,ke+1,n)
             endif
             if (nokb(i,j) .eq. 3) then
               ab(i,j,ke+1,n) = abokb(i,j,1,n)
               ab(i,j,ke+2,n) = abokb(i,j,2,n)
             endif
             if (nokb(i,j) .eq. 4) then
               ab(i,j,ke+1,n) = ab(i,j,ks  ,n)
               ab(i,j,ke+2,n) = ab(i,j,ks+1,n)
             endif
             if (nokb(i,j) .eq. 5) then
               ab(i,j,ke+1,n) = ab(i,j,ke  ,n)
               ab(i,j,ke+2,n) = ab(i,j,ke-1,n)
             endif
            enddo
           enddo
         enddo
         bvstat(6,7) = ru3
         endif
#endif /* NO MPI */
c
      return
      end
