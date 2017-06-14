c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                   B V A L T               \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c=======================================================================
       subroutine bvalt ( rl1, ru1, rl2, ru2, rl3, ru3, t )
c
c    PURPOSE: updates temperature in boundary zones; only needed if
c             LEOS > 1.
c
c    Cloned from BVALE by John hayes
c
c    i-boundaries:                    j = js  , je     k = ks  , ke
c    j-boundaries:   i = is  , ie                      k = ks  , ke
c    k-boundaries:   i = is  , ie     j = js  , je
c
c  See comments in BVALD.
c
c  Flags l1, l2, l3, activate the 1-, 2-, and 3- loops when nonzero.
c  Their values give the number of layers to pass.
c  Flag iupper activates enables sends    in the "m" direction and
c                                receives in the "p" direction when
c  nonzero. 
c
c  EXTERNALS: [NONE]
c
c-----------------------------------------------------------------------
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
      real(rl) :: t(in,jn,kn)
c
      integer  :: i,j,k,l1,l2,l3,u1,u2,u3,
     .            rl1,rl2,rl3,ru1,ru2,ru3,
     .            ls,ll,lu,us,ul,uu
c
c-----------------------------------------------------------------------
c------------------------  I - B O U N D A R Y  ------------------------
c-----------------------------------------------------------------------
c
       l1 = rl1 - bvstat(1,8)
       u1 = ru1 - bvstat(2,8)
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
             call MPI_IRECV( t(is-ll+i-1,   1,   1), 1, i_slice, n1m
     &                     ,13100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,8) = rl1
         endif
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_ISEND( t(is+uu+i-1,   1,   1), 1, i_slice, n1m
     &                   ,13200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,8) = ru1
         endif
       else
         if (l1 .gt. 0) then
         do k=ks-1,ke+1
cdir$ ivdep
           do j=js-1,je+1
             if ( abs(niib(j,k)) .eq. 1) then
               t(is-1,j,k) = t(is  ,j,k)
               t(is-2,j,k) = t(is+1,j,k)
             endif
             if (niib(j,k) .eq. 2) then
             endif
             if (niib(j,k) .eq. 3) then
               t(is-1,j,k) = t(is  ,j,k)
               t(is-2,j,k) = t(is-1,j,k)
             endif
             if (niib(j,k) .eq. 5) then
               t(is-1,j,k) = t(is  ,j,k)
               t(is-2,j,k) = t(is+1,j,k)
             endif
           enddo
         enddo
         bvstat(1,8) = rl1
         endif
       endif
c
c      Outer i boundary.
c
       if (nois(1).eq.0 .or. nois(1).eq.4) then
         if (u1 .gt. 0) then
           do i=1,us
             nreq = nreq + 1
             call MPI_IRECV( t(ie+i+uu,   1,   1), 1, i_slice, n1p
     &                   ,13200+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(2,8) = ru1
         endif
         if (l1 .gt. 0) then
           do i=1,ls
             nreq = nreq + 1
             call MPI_ISEND( t(ie+i-ll,   1,   1), 1, i_slice, n1p
     &                     ,13100+25*i+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(1,8) = rl1
         endif
       else
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
cdir$ ivdep
           do j=js-1,je+1
             if ( abs(noib(j,k)) .eq. 1) then
               t(ie+1,j,k) = t(ie,j,k)
               t(ie+2,j,k) = t(ie-1,j,k)
             endif
             if (noib(j,k) .eq. 2) then
               t(ie+1,j,k) = t(ie,j,k)
               t(ie+2,j,k) = t(ie+1,j,k)
             endif
             if (noib(j,k) .eq. 3) then
               t(ie+1,j,k) = t(ie,j,k)
               t(ie+2,j,k) = t(ie+1,j,k)
             endif
             if (noib(j,k) .eq. 5) then
               t(ie+1,j,k) = t(ie  ,j,k)
               t(ie+2,j,k) = t(ie-1,j,k)
             endif
           enddo
         enddo
         bvstat(2,8) = ru1
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED 
         if (l1 .gt. 0) then
         do k=ks-1,ke+1
cdir$ ivdep
           do j=js-1,je+1
             if ( abs(niib(j,k)) .eq. 1) then
               t(is-1,j,k) = t(is  ,j,k)
               t(is-2,j,k) = t(is+1,j,k)
             endif
             if (niib(j,k) .eq. 2) then
               t(is-1,j,k) = t(is  ,j,k)
               t(is-2,j,k) = t(is-1,j,k)
             endif
             if (niib(j,k) .eq. 3) then
               t(is-1,j,k) = t(is  ,j,k)
               t(is-2,j,k) = t(is-1,j,k)
             endif
             if (niib(j,k) .eq. 4) then
               t(is-1,j,k) = t(ie  ,j,k)
               t(is-2,j,k) = t(ie-1,j,k)
             endif
             if (niib(j,k) .eq. 5) then
               t(is-1,j,k) = t(is  ,j,k)
               t(is-2,j,k) = t(is+1,j,k)
             endif
           enddo
         enddo
         bvstat(1,8) = rl1
         endif
c
c      Outer i boundary.
c
         if (u1 .gt. 0) then
         do k=ks-1,ke+1
cdir$ ivdep
           do j=js-1,je+1
             if ( abs(noib(j,k)) .eq. 1) then
               t(ie+1,j,k) = t(ie,j,k)
               t(ie+2,j,k) = t(ie-1,j,k)
             endif
             if (noib(j,k) .eq. 2) then
               t(ie+1,j,k) = t(ie,j,k)
               t(ie+2,j,k) = t(ie+1,j,k)
             endif
             if (noib(j,k) .eq. 3) then
               t(ie+1,j,k) = t(ie,j,k)
               t(ie+2,j,k) = t(ie+1,j,k)
             endif
             if (noib(j,k) .eq. 4) then
               t(ie+1,j,k) = t(is  ,j,k)
               t(ie+2,j,k) = t(is+1,j,k)
             endif
             if (noib(j,k) .eq. 5) then
               t(ie+1,j,k) = t(ie  ,j,k)
               t(ie+2,j,k) = t(ie-1,j,k)
             endif
           enddo
         enddo
         bvstat(2,8) = ru1
         endif
#endif /* NO MPI */
c
c-----------------------------------------------------------------------
c------------------------  J - B O U N D A R Y  ------------------------
c-----------------------------------------------------------------------
c
       l2 = rl2 - bvstat(3,8)
       u2 = ru2 - bvstat(4,8)
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
             call MPI_IRECV( t(   1,js-ll+j-1,   1), 1, j_slice, n2m
     &                   ,13300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,8) = rl2
         endif
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_ISEND( t(   1,js+uu+j-1,   1), 1, j_slice, n2m
     &                   ,13400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,8) = ru2
         endif
       else
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
cdir$ ivdep
           do i=is-1,ie+1
             if ( abs(nijb(i,k)) .eq. 1) then
               t(i,js-1,k) = t(i,js  ,k)
               t(i,js-2,k) = t(i,js+1,k)
             endif
             if (nijb(i,k) .eq. 2) then
               t(i,js-1,k) = t(i,js  ,k)
               t(i,js-2,k) = t(i,js-1,k)
             endif
             if (nijb(i,k) .eq. 3) then
               t(i,js-1,k) = t(i,js  ,k)
               t(i,js-2,k) = t(i,js-1,k)
             endif
             if (nijb(i,k) .eq. 5) then
               t(i,js-1,k) = t(i,js  ,k)
               t(i,js-2,k) = t(i,js+1,k)
             endif
           enddo
         enddo
         bvstat(3,8) = rl2
         endif
       endif
c
c      Outer j boundary.
c
       if (nojs(1).eq.0 .or. nojs(1).eq.4) then
         if (u2 .gt. 0) then
           do j=1,us
             nreq = nreq + 1
             call MPI_IRECV( t(   1,je+j+uu,   1), 1, j_slice, n2p
     &                   ,13400+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(4,8) = ru2
         endif
         if (l2 .gt. 0) then
           do j=1,ls
             nreq = nreq + 1
             call MPI_ISEND( t(   1,je+j-ll,   1), 1, j_slice, n2p
     &                   ,13300+25*j+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(3,8) = rl2
         endif
       else
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
cdir$ ivdep
           do i=is-1,ie+1
             if ( abs(nojb(i,k)) .eq. 1) then
               t(i,je+1,k) = t(i,je  ,k)
               t(i,je+2,k) = t(i,je-1,k)
             endif
             if (nojb(i,k) .eq. 2) then
               t(i,je+1,k) = t(i,je  ,k)
               t(i,je+2,k) = t(i,je+1,k)
             endif
             if (nojb(i,k) .eq. 3) then
               t(i,je+1,k) = t(i,je  ,k)
               t(i,je+2,k) = t(i,je+1,k)
             endif
             if (nojb(i,k) .eq. 5) then
               t(i,je+1,k) = t(i,je  ,k)
               t(i,je+2,k) = t(i,je-1,k)
             endif
           enddo
         enddo
         bvstat(4,8) = ru2
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED 
         if (l2 .gt. 0) then
         do k=ks-1,ke+1
cdir$ ivdep
           do i=is-1,ie+1
             if ( abs(nijb(i,k)) .eq. 1) then
               t(i,js-1,k) = t(i,js  ,k)
               t(i,js-2,k) = t(i,js+1,k)
             endif
             if (nijb(i,k) .eq. 2) then
               t(i,js-1,k) = t(i,js  ,k)
               t(i,js-2,k) = t(i,js-1,k)
             endif
             if (nijb(i,k) .eq. 3) then
               t(i,js-1,k) = t(i,js  ,k)
               t(i,js-2,k) = t(i,js-1,k)
             endif
             if (nijb(i,k) .eq. 4) then
               t(i,js-1,k) = t(i,je  ,k)
               t(i,js-2,k) = t(i,je-1,k)
             endif
             if (nijb(i,k) .eq. 5) then
               t(i,js-1,k) = t(i,js  ,k)
               t(i,js-2,k) = t(i,js+1,k)
             endif
           enddo
         enddo
         bvstat(3,8) = rl2
         endif
c
c      Outer j boundary.
c
         if (u2 .gt. 0) then
         do k=ks-1,ke+1
cdir$ ivdep
           do i=is-1,ie+1
             if ( abs(nojb(i,k)) .eq. 1) then
               t(i,je+1,k) = t(i,je  ,k)
               t(i,je+2,k) = t(i,je-1,k)
             endif
             if (nojb(i,k) .eq. 2) then
               t(i,je+1,k) = t(i,je  ,k)
               t(i,je+2,k) = t(i,je+1,k)
             endif
             if (nojb(i,k) .eq. 3) then
               t(i,je+1,k) = t(i,je  ,k)
               t(i,je+2,k) = t(i,je+1,k)
             endif
             if (nojb(i,k) .eq. 4) then
               t(i,je+1,k) = t(i,js  ,k)
               t(i,je+2,k) = t(i,js+1,k)
             endif
             if (nojb(i,k) .eq. 5) then
               t(i,je+1,k) = t(i,je  ,k)
               t(i,je+2,k) = t(i,je-1,k)
             endif
           enddo
         enddo
         bvstat(4,8) = ru2
         endif
#endif /* NO MPI */
c
c-----------------------------------------------------------------------
c------------------------  K - B O U N D A R Y  ------------------------
c-----------------------------------------------------------------------
c
       l3 = rl3 - bvstat(5,8)
       u3 = ru3 - bvstat(6,8)
c
c      Inner k boundary.
c
#ifdef MPI_USED
c
c Count slabs, compute positioning indices.
c
       ls = max (l3-1,1)  ! number of inner slabs to send/receive
       ll = min (l3,2)    ! index for lower plane; 1 or 2
       lu = ll - ls       ! index for upper plane; 0 if l3 not
       us = max (u3-1,1)  ! number of outer slabs to send/receive
       ul = min (u3,2)    ! index for lower plane; 1 or 2
       uu = ul - us       ! index for upper plane; 0 if l3 not 2
c
       if (niks(1).eq.0 .or. niks(1).eq.4) then
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_IRECV( t(   1,   1,ks-ll+k-1), 1, k_slice, n3m
     &                     ,13500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,8) = rl3
         endif
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_ISEND( t(   1,   1,ks+uu+k-1), 1, k_slice, n3m
     &                     ,13600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,8) = ru3
         endif
       else
         if (l3 .gt. 0) then
         do j=js-1,je+1
cdir$ ivdep
           do i=is-1,ie+1
             if ( abs(nikb(i,j)) .eq. 1) then
               t(i,j,ks-1) = t(i,j,ks  )
               t(i,j,ks-2) = t(i,j,ks+1)
             endif
             if (nikb(i,j) .eq. 2) then
               t(i,j,ks-1) = t(i,j,ks  )
               t(i,j,ks-2) = t(i,j,ks-1)
             endif
             if (nikb(i,j) .eq. 3) then
               t(i,j,ks-1) = t(i,j,ks  )
               t(i,j,ks-2) = t(i,j,ks-1)
             endif
             if (nikb(i,j) .eq. 5) then
               t(i,j,ks-1) = t(i,j,ks  )
               t(i,j,ks-2) = t(i,j,ks+1)
             endif
           enddo
         enddo
         bvstat(5,8) = rl3
         endif
       endif
c
c      Outer k boundary.
c
       if (noks(1).eq.0 .or. noks(1).eq.4) then
         if (u3 .gt. 0) then
           do k=1,us
             nreq = nreq + 1
             call MPI_IRECV( t(   1,   1,ke+k+uu), 1, k_slice, n3p
     &                     ,13600+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(6,8) = ru3
         endif
         if (l3 .gt. 0) then
           do k=1,ls
             nreq = nreq + 1
             call MPI_ISEND( t(   1,   1,ke+k-ll), 1, k_slice, n3p
     &                     ,13500+25*k+nsub, comm3d, req(nreq), ierr)
           enddo
           bvstat(5,8) = rl3
         endif
       else
         if (u3 .gt. 0) then
         do j=js-1,je+1
cdir$ ivdep
           do i=is-1,ie+1
             if ( abs(nokb(i,j)) .eq. 1) then
               t(i,j,ke+1) = t(i,j,ke  )
               t(i,j,ke+2) = t(i,j,ke-1)
             endif
             if (nokb(i,j) .eq. 2) then
               t(i,j,ke+1) = t(i,j,ke  )
               t(i,j,ke+2) = t(i,j,ke+1)
             endif
             if (nokb(i,j) .eq. 3) then
               t(i,j,ke+1) = t(i,j,ke  )
               t(i,j,ke+2) = t(i,j,ke+1)
             endif
             if (nokb(i,j) .eq. 5) then
               t(i,j,ke+1) = t(i,j,ke  )
               t(i,j,ke+2) = t(i,j,ke-1)
             endif
           enddo
         enddo
         bvstat(6,8) = ru3
         endif
       endif
#endif /* MPI */
#ifndef MPI_USED
         if (l3 .gt. 0) then
         do j=js-1,je+1
cdir$ ivdep
           do i=is-1,ie+1
             if ( abs(nikb(i,j)) .eq. 1) then
               t(i,j,ks-1) = t(i,j,ks  )
               t(i,j,ks-2) = t(i,j,ks+1)
             endif
             if (nikb(i,j) .eq. 2) then
               t(i,j,ks-1) = t(i,j,ks  )
               t(i,j,ks-2) = t(i,j,ks-1)
             endif
             if (nikb(i,j) .eq. 3) then
               t(i,j,ks-1) = t(i,j,ks  )
               t(i,j,ks-2) = t(i,j,ks-1)
             endif
             if (nikb(i,j) .eq. 4) then
               t(i,j,ks-1) = t(i,j,ke  )
               t(i,j,ks-2) = t(i,j,ke-1)
             endif
             if (nikb(i,j) .eq. 5) then
               t(i,j,ks-1) = t(i,j,ks  )
               t(i,j,ks-2) = t(i,j,ks+1)
             endif
           enddo
         enddo
         bvstat(5,8) = rl3
         endif
c
c      Outer k boundary.
c
         if (u3 .gt. 0) then
         do j=js-1,je+1
cdir$ ivdep
           do i=is-1,ie+1
             if ( abs(nokb(i,j)) .eq. 1) then
               t(i,j,ke+1) = t(i,j,ke  )
               t(i,j,ke+2) = t(i,j,ke-1)
             endif
             if (nokb(i,j) .eq. 2) then
               t(i,j,ke+1) = t(i,j,ke  )
               t(i,j,ke+2) = t(i,j,ke+1)
             endif
             if (nokb(i,j) .eq. 3) then
               t(i,j,ke+1) = t(i,j,ke  )
               t(i,j,ke+2) = t(i,j,ke+1)
             endif
             if (nokb(i,j) .eq. 4) then
               t(i,j,ke+1) = t(i,j,ks  )
               t(i,j,ke+2) = t(i,j,ks+1)
             endif
             if (nokb(i,j) .eq. 5) then
               t(i,j,ke+1) = t(i,j,ke  )
               t(i,j,ke+2) = t(i,j,ke-1)
             endif
           enddo
         enddo
         bvstat(6,8) = ru3
         endif
#endif /* NO MPI */
c
      return
      end
