c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              B V A L E M F 1              \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine bvalemf1 ( v2b3, v3b2 )
c
c    dac:zeus3d.bvalemf1 <-------------- boundary values for 1-emf terms
c                                                         february, 1990
c
c    written by: David Clarke, February 1990
c    modified 1: September, 1990 by David Clarke; moved magnetic fields
c                to face-centres.
c    modified 2: minimal rewrite for ZEUS-MP by M-MML 9.3.98
c    modified 3: added restrictions for symmetry about the J and K axes;
c                John Hayes, 10/2005
c
c  PURPOSE: This routine sets boundary values for the two terms in the
c  1-emf (centred on the 1-edges).  The "active" zones for "emf1" are:
c
c    i = is to ie;  j = js+1 to je;  k = ks+1 to ke
c
c  In order to update both the active and ghost zones of the 2- and
c  3-magnetic field components, all edges in the boundary regions are
c  required.  This gives a complete grid of values for "emf1".  Thus,
c  the ranges for the boundary values are:
c
c    i-boundaries:                    j = js  , je+1   k = ks  , ke+1
c    j-boundaries:   i = is-2, ie+2                    k = ks  , ke+1
c    k-boundaries:   i = is-2, ie+2   j = js-2, je+3
c
c  Note that the boundary values must be set even if is > ismn, etc.
c  because the emfs are stored in worker arrays and it is likely that
c  the boundary values have been overwritten.
c
c  See comments in BVALD.
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
      integer  :: i, j, k, ii, jj, kk
c
      real(rl) :: v2b3(  in,  jn,  kn), v3b2(  in,  jn,  kn)
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c------------------------  I - B O U N D A R Y  ------------------------
c-----------------------------------------------------------------------
c
#ifdef MPI_USED
         nreq = 0
         nsub = nsub + 1
c
c Post a receive for a slab of data from the interior of the 
c neighboring tile to fill my ghost zones.  Initiate a send 
c to pass a slab of my interior data for my neighbor's ghost zones.
c 
       if (niis(1).eq.0 .or. niis(1).eq.4) then
c Trying to send two i_slices at once did not work, so use a loop here 
c to grab both boundary slices   M-MML  
c
         do ii = 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v2b3(is  +ii, 1, 1), 1, i_slice, n1m, 
     &                   7200+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v3b2(is  +ii, 1, 1), 1, i_slice, n1m, 
     &                   7300+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v2b3(is-2+ii, 1, 1), 1, i_slice, n1m, 
     &                   7000+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v3b2(is-2+ii, 1, 1), 1, i_slice, n1m, 
     &                   7100+nsub+ii, comm3d , req(nreq), ierr)
         enddo
       else
c
c      Inner i boundary.
c
         do k=ks,ke+1
cdir$ ivdep
           do j=js,je+1
             if ( abs(niib23(j,k)) .eq. 1) then
               v2b3(is-1,j,k) = v2b3(is  ,j,k)
               v2b3(is-2,j,k) = v2b3(is+1,j,k)
               v3b2(is-1,j,k) = v3b2(is  ,j,k)
               v3b2(is-2,j,k) = v3b2(is+1,j,k)
             endif
             if (niib23(j,k) .eq. 2) then
               v2b3(is-1,j,k) = v2b3(is  ,j,k)
               v2b3(is-2,j,k) = v2b3(is  ,j,k)
               v3b2(is-1,j,k) = v3b2(is  ,j,k)
               v3b2(is-2,j,k) = v3b2(is  ,j,k)
             endif
             if (niib23(j,k) .eq. 3) then
               v2b3(is-1,j,k) = emf1iib(j,k,1)
               v2b3(is-2,j,k) = emf1iib(j,k,2)
               v3b2(is-1,j,k) = 0.0
               v3b2(is-2,j,k) = 0.0
             endif
             if (niib23(j,k) .eq. 5) then
               v2b3(is-1,j,k) =-v2b3(is  ,j,k)
               v2b3(is-2,j,k) =-v2b3(is+1,j,k)
               v3b2(is-1,j,k) =-v3b2(is  ,j,k)
               v3b2(is-2,j,k) =-v3b2(is+1,j,k)
             endif
           enddo
         enddo
       endif
c
c      Outer i boundary.
c
c
c Post a receive for a slab of data from the interior of the 
c neighboring tile to fill my ghost zones.  Initiate a send 
c to pass a slab of my interior data for my neighbor's ghost zones.
c 
       if (nois(1).eq.0 .or. nois(1).eq.4) then
         do ii= 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v2b3(ie-1+ii, 1, 1), 1, i_slice, n1p, 
     &                   7000+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v3b2(ie-1+ii, 1, 1), 1, i_slice, n1p, 
     &                   7100+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v2b3(ie+1+ii, 1, 1), 1, i_slice, n1p, 
     &                   7200+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v3b2(ie+1+ii, 1, 1), 1, i_slice, n1p, 
     &                   7300+nsub+ii, comm3d , req(nreq), ierr)
         enddo
       else
c
         do k=ks,ke+1
cdir$ ivdep
           do j=js,je+1
             if ( abs(noib23(j,k)) .eq. 1) then
               v2b3(ie+1,j,k) = v2b3(ie  ,j,k)
               v2b3(ie+2,j,k) = v2b3(ie-1,j,k)
               v3b2(ie+1,j,k) = v3b2(ie  ,j,k)
               v3b2(ie+2,j,k) = v3b2(ie-1,j,k)
             endif
             if (noib23(j,k) .eq. 2) then
               v2b3(ie+1,j,k) = v2b3(ie  ,j,k)
               v2b3(ie+2,j,k) = v2b3(ie  ,j,k)
               v3b2(ie+1,j,k) = v3b2(ie  ,j,k)
               v3b2(ie+2,j,k) = v3b2(ie  ,j,k)
             endif
             if (noib23(j,k) .eq. 3) then
               v2b3(ie+1,j,k) = emf1oib(j,k,1)
               v2b3(ie+2,j,k) = emf1oib(j,k,2)
               v3b2(ie+1,j,k) = 0.0
               v3b2(ie+2,j,k) = 0.0
             endif
             if (noib23(j,k) .eq. 5) then
               v2b3(ie+1,j,k) =-v2b3(ie  ,j,k)
               v2b3(ie+2,j,k) =-v2b3(ie-1,j,k)
               v3b2(ie+1,j,k) =-v3b2(ie  ,j,k)
               v3b2(ie+2,j,k) =-v3b2(ie-1,j,k)
             endif
c
           enddo
         enddo
       endif
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )

#endif /* MPI */
#ifndef MPI_USED
c
c      Inner i boundary.
c
         do k=ks,ke+1
cdir$ ivdep
           do j=js,je+1
             if ( abs(niib23(j,k)) .eq. 1) then
               v2b3(is-1,j,k) = v2b3(is  ,j,k)
               v2b3(is-2,j,k) = v2b3(is+1,j,k)
               v3b2(is-1,j,k) = v3b2(is  ,j,k)
               v3b2(is-2,j,k) = v3b2(is+1,j,k)
             endif
             if (niib23(j,k) .eq. 2) then
               v2b3(is-1,j,k) = v2b3(is  ,j,k)
               v2b3(is-2,j,k) = v2b3(is  ,j,k)
               v3b2(is-1,j,k) = v3b2(is  ,j,k)
               v3b2(is-2,j,k) = v3b2(is  ,j,k)
             endif
             if (niib23(j,k) .eq. 3) then
               v2b3(is-1,j,k) = emf1iib(j,k,1)
               v2b3(is-2,j,k) = emf1iib(j,k,2)
               v3b2(is-1,j,k) = 0.0
               v3b2(is-2,j,k) = 0.0
             endif
             if (niib23(j,k) .eq. 4) then
               v2b3(is-1,j,k) = v2b3(ie  ,j,k)
               v2b3(is-2,j,k) = v2b3(ie-1,j,k)
               v3b2(is-1,j,k) = v3b2(ie  ,j,k)
               v3b2(is-2,j,k) = v3b2(ie-1,j,k)
             endif
             if (niib23(j,k) .eq. 5) then
               v2b3(is-1,j,k) =-v2b3(is  ,j,k)
               v2b3(is-2,j,k) =-v2b3(is+1,j,k)
               v3b2(is-1,j,k) =-v3b2(is  ,j,k)
               v3b2(is-2,j,k) =-v3b2(is+1,j,k)
             endif
           enddo
         enddo
c
c      Outer i boundary.
c
c
         do k=ks,ke+1
cdir$ ivdep
           do j=js,je+1
             if ( abs(noib23(j,k)) .eq. 1) then
               v2b3(ie+1,j,k) = v2b3(ie  ,j,k)
               v2b3(ie+2,j,k) = v2b3(ie-1,j,k)
               v3b2(ie+1,j,k) = v3b2(ie  ,j,k)
               v3b2(ie+2,j,k) = v3b2(ie-1,j,k)
             endif
             if (noib23(j,k) .eq. 2) then
               v2b3(ie+1,j,k) = v2b3(ie  ,j,k)
               v2b3(ie+2,j,k) = v2b3(ie  ,j,k)
               v3b2(ie+1,j,k) = v3b2(ie  ,j,k)
               v3b2(ie+2,j,k) = v3b2(ie  ,j,k)
             endif
             if (noib23(j,k) .eq. 3) then
               v2b3(ie+1,j,k) = emf1oib(j,k,1)
               v2b3(ie+2,j,k) = emf1oib(j,k,2)
               v3b2(ie+1,j,k) = 0.0
               v3b2(ie+2,j,k) = 0.0
             endif
             if (noib23(j,k) .eq. 4) then
               v2b3(ie+1,j,k) = v2b3(is  ,j,k)
               v2b3(ie+2,j,k) = v2b3(is+1,j,k)
               v3b2(ie+1,j,k) = v3b2(is  ,j,k)
               v3b2(ie+2,j,k) = v3b2(is+1,j,k)
             endif
             if (noib23(j,k) .eq. 5) then
               v2b3(ie+1,j,k) =-v2b3(ie  ,j,k)
               v2b3(ie+2,j,k) =-v2b3(ie-1,j,k)
               v3b2(ie+1,j,k) =-v3b2(ie  ,j,k)
               v3b2(ie+2,j,k) =-v3b2(ie-1,j,k)
             endif
c
           enddo
         enddo
#endif /* NO MPI */
c
c
c-----------------------------------------------------------------------
c------------------------  J - B O U N D A R Y  ------------------------
c-----------------------------------------------------------------------
c
      if(ldimen .eq. 1) return
c
c      Inner j boundary.
c
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1

c  add a loop to get around j_slice problem M-MML 25.3.98
       if (nijs(1).eq.0 .or. nijs(1).eq.4) then
         do jj = 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v2b3(1, js+1+jj, 1), 1, j_slice, n2m, 
     &             7600+nsub+jj, comm3d, req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v3b2(1, js+1+jj, 1), 1, j_slice, n2m, 
     &             7700+nsub+jj, comm3d, req(nreq),ierr)
           nreq = nreq + 1
           call MPI_IRECV(v2b3(1, js-2+jj, 1), 1, j_slice, n2m, 
     &             7400+nsub+jj, comm3d, req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v3b2(1, js-2+jj, 1), 1, j_slice, n2m, 
     &             7500+nsub+jj, comm3d, req(nreq), ierr)
         enddo
       else
         do k=ks,ke+1
cdir$ ivdep
           do i=is-2,ie+2
             if (nijb3(i,k) .eq. 1) then
               v2b3(i,js  ,k) = 0.0
               v2b3(i,js-1,k) =-v2b3(i,js+1,k)
               v2b3(i,js-2,k) =-v2b3(i,js+2,k)
               v3b2(i,js  ,k) = 0.0
               v3b2(i,js-1,k) =-v3b2(i,js+1,k)
               v3b2(i,js-2,k) =-v3b2(i,js+2,k)
             endif
             if (nijb3(i,k) .eq.-1) then
              if(lgeom .eq. 1) then
               v2b3(i,js  ,k) = 0.0
               v2b3(i,js-1,k) =-v2b3(i,js+1,k)
               v2b3(i,js-2,k) =-v2b3(i,js+2,k)
               v3b2(i,js  ,k) = 0.0
               v3b2(i,js-1,k) =-v3b2(i,js+1,k)
               v3b2(i,js-2,k) =-v3b2(i,js+2,k)
              else ! lgeom
               v2b3(i,js  ,k) = 0.0
               v2b3(i,js-1,k) = v2b3(i,js+1,k)
               v2b3(i,js-2,k) = v2b3(i,js+2,k)
               v3b2(i,js  ,k) = 0.0
               v3b2(i,js-1,k) = v3b2(i,js+1,k)
               v3b2(i,js-2,k) = v3b2(i,js+2,k)
              endif ! lgeom
             endif
             if (nijb3(i,k) .eq. 2) then
               v2b3(i,js-1,k) = v2b3(i,js  ,k)
               v2b3(i,js-2,k) = v2b3(i,js  ,k)
               v3b2(i,js-1,k) = v3b2(i,js  ,k)
               v3b2(i,js-2,k) = v3b2(i,js  ,k)
             endif
             if (nijb3(i,k) .eq. 3) then
               v2b3(i,js  ,k) = emf1ijb(i,k,1)
               v2b3(i,js-1,k) = emf1ijb(i,k,2)
               v2b3(i,js-2,k) = emf1ijb(i,k,3)
               v3b2(i,js  ,k) = 0.0
               v3b2(i,js-1,k) = 0.0
               v3b2(i,js-2,k) = 0.0
             endif
             if (nijb3(i,k) .eq. 5) then
               v2b3(i,js  ,k) = 0.0
               v2b3(i,js-1,k) = v2b3(i,js+1,k)
               v2b3(i,js-2,k) = v2b3(i,js+2,k)
               v3b2(i,js-1,k) = v3b2(i,js+1,k)
               v3b2(i,js-2,k) = v3b2(i,js+2,k)
             endif
           enddo
         enddo
       endif
c
c      Outer j boundary.
c
c
c  modified to only pass 1 j_slice at a time M-MML 25.3.98
       if (nojs(1).eq.0 .or. nojs(1).eq.4) then
         do jj = 0,1
           nreq = nreq + 1
           call MPI_ISEND(v2b3(1, je-1+jj, 1), 1, j_slice, n2p, 
     &              7400+nsub+jj, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v3b2(1, je-1+jj, 1), 1, j_slice, n2p, 
     &              7500+nsub+jj, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v2b3(1, je+2+jj, 1), 1, j_slice, n2p, 
     &              7600+nsub+jj, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v3b2(1, je+2+jj, 1), 1, j_slice, n2p, 
     &              7700+nsub+jj, comm3d , req(nreq), ierr)
         enddo
       else
         do k=ks,ke+1
cdir$ ivdep
           do i=is-2,ie+2
             if (nojb3(i,k) .eq. 1) then
               v2b3(i,je+1,k) = 0.0
               v2b3(i,je+2,k) =-v2b3(i,je  ,k)
               v2b3(i,je+3,k) =-v2b3(i,je-1,k)
               v3b2(i,je+1,k) = 0.0
               v3b2(i,je+2,k) =-v3b2(i,je  ,k)
               v3b2(i,je+3,k) =-v3b2(i,je-1,k)
             endif
             if (nojb3(i,k) .eq.-1) then
              if(lgeom .eq. 3) then
               v2b3(i,je+1,k) = 0.0
               v2b3(i,je+2,k) = v2b3(i,je  ,k)
               v2b3(i,je+3,k) = v2b3(i,je-1,k)
               v3b2(i,je+1,k) = 0.0
               v3b2(i,je+2,k) = v3b2(i,je  ,k)
               v3b2(i,je+3,k) = v3b2(i,je-1,k)
              endif
              if(lgeom .ne. 3) then
               v2b3(i,je+1,k) = 0.0
               v2b3(i,je+2,k) =-v2b3(i,je  ,k)
               v2b3(i,je+3,k) =-v2b3(i,je-1,k)
               v3b2(i,je+1,k) = 0.0
               v3b2(i,je+2,k) =-v3b2(i,je  ,k)
               v3b2(i,je+3,k) =-v3b2(i,je-1,k)
              endif
             endif
             if (nojb3(i,k) .eq. 2) then
c*if alias PROBLEM.eq.advect
c               v2b3(i,je+2,k) = 3.0 * ( v2b3(i,je+1,k) - v2b3(i,je  ,k) )
c       1                      + v2b3(i,je-1,k)
c               v2b3(i,je+3,k) = 3.0 * ( v2b3(i,je+2,k) - v2b3(i,je+1,k) )
c       1                      + v2b3(i,je  ,k)
c*else
               v2b3(i,je+2,k) = v2b3(i,je+1,k)
               v2b3(i,je+3,k) = v2b3(i,je+1,k)
               v3b2(i,je+2,k) = v3b2(i,je+1,k)
               v3b2(i,je+3,k) = v3b2(i,je+1,k)
             endif
             if (nojb3(i,k) .eq. 3) then
               v2b3(i,je+1,k) = emf1ojb(i,k,1)
               v2b3(i,je+2,k) = emf1ojb(i,k,2)
               v2b3(i,je+3,k) = emf1ojb(i,k,3)
               v3b2(i,je+1,k) = 0.0
               v3b2(i,je+2,k) = 0.0
               v3b2(i,je+3,k) = 0.0
             endif
             if (nojb3(i,k) .eq. 5) then
               v2b3(i,je+1,k) = 0.0
               v2b3(i,je+2,k) = v2b3(i,je  ,k)
               v2b3(i,je+3,k) = v2b3(i,je-1,k)
               v3b2(i,je+2,k) = v3b2(i,je  ,k)
               v3b2(i,je+3,k) = v3b2(i,je-1,k)
             endif
c
           enddo
         enddo
       endif

       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )

#endif /* MPI */
#ifndef MPI_USED
         do k=ks,ke+1
cdir$ ivdep
           do i=is-2,ie+2
             if (nijb3(i,k) .eq. 1) then
               v2b3(i,js  ,k) = 0.0
               v2b3(i,js-1,k) =-v2b3(i,js+1,k)
               v2b3(i,js-2,k) =-v2b3(i,js+2,k)
               v3b2(i,js  ,k) = 0.0
               v3b2(i,js-1,k) =-v3b2(i,js+1,k)
               v3b2(i,js-2,k) =-v3b2(i,js+2,k)
             endif
             if (nijb3(i,k) .eq.-1) then
              if(lgeom .eq. 1) then
               v2b3(i,js  ,k) = 0.0
               v2b3(i,js-1,k) =-v2b3(i,js+1,k)
               v2b3(i,js-2,k) =-v2b3(i,js+2,k)
               v3b2(i,js  ,k) = 0.0
               v3b2(i,js-1,k) =-v3b2(i,js+1,k)
               v3b2(i,js-2,k) =-v3b2(i,js+2,k)
              else ! lgeom
               v2b3(i,js  ,k) = 0.0
               v2b3(i,js-1,k) = v2b3(i,js+1,k)
               v2b3(i,js-2,k) = v2b3(i,js+2,k)
               v3b2(i,js  ,k) = 0.0
               v3b2(i,js-1,k) = v3b2(i,js+1,k)
               v3b2(i,js-2,k) = v3b2(i,js+2,k)
              endif ! lgeom
             endif
             if (nijb3(i,k) .eq. 2) then
               v2b3(i,js-1,k) = v2b3(i,js  ,k)
               v2b3(i,js-2,k) = v2b3(i,js  ,k)
               v3b2(i,js-1,k) = v3b2(i,js  ,k)
               v3b2(i,js-2,k) = v3b2(i,js  ,k)
             endif
             if (nijb3(i,k) .eq. 3) then
               v2b3(i,js  ,k) = emf1ijb(i,k,1)
               v2b3(i,js-1,k) = emf1ijb(i,k,2)
               v2b3(i,js-2,k) = emf1ijb(i,k,3)
               v3b2(i,js  ,k) = 0.0
               v3b2(i,js-1,k) = 0.0
               v3b2(i,js-2,k) = 0.0
             endif
             if (nijb3(i,k) .eq. 4) then
               v2b3(i,js-1,k) = v2b3(i,je  ,k)
               v2b3(i,js-2,k) = v2b3(i,je-1,k)
               v3b2(i,js-1,k) = v3b2(i,je  ,k)
               v3b2(i,js-2,k) = v3b2(i,je-1,k)
             endif
             if (nijb3(i,k) .eq. 5) then
               v2b3(i,js  ,k) = 0.0
               v2b3(i,js-1,k) = v2b3(i,js+1,k)
               v2b3(i,js-2,k) = v2b3(i,js+2,k)
               v3b2(i,js-1,k) = v3b2(i,js+1,k)
               v3b2(i,js-2,k) = v3b2(i,js+2,k)
             endif
           enddo
         enddo
c
c      Outer j boundary.
c
         do k=ks,ke+1
cdir$ ivdep
           do i=is-2,ie+2
             if (nojb3(i,k) .eq. 1) then
               v2b3(i,je+1,k) = 0.0
               v2b3(i,je+2,k) =-v2b3(i,je  ,k)
               v2b3(i,je+3,k) =-v2b3(i,je-1,k)
               v3b2(i,je+1,k) = 0.0
               v3b2(i,je+2,k) =-v3b2(i,je  ,k)
               v3b2(i,je+3,k) =-v3b2(i,je-1,k)
             endif
             if (nojb3(i,k) .eq.-1) then
              if(lgeom .eq. 3) then
               v2b3(i,je+1,k) = 0.0
               v2b3(i,je+2,k) = v2b3(i,je  ,k)
               v2b3(i,je+3,k) = v2b3(i,je-1,k)
               v3b2(i,je+1,k) = 0.0
               v3b2(i,je+2,k) = v3b2(i,je  ,k)
               v3b2(i,je+3,k) = v3b2(i,je-1,k)
              endif
              if(lgeom .ne. 3) then
               v2b3(i,je+1,k) = 0.0
               v2b3(i,je+2,k) =-v2b3(i,je  ,k)
               v2b3(i,je+3,k) =-v2b3(i,je-1,k)
               v3b2(i,je+1,k) = 0.0
               v3b2(i,je+2,k) =-v3b2(i,je  ,k)
               v3b2(i,je+3,k) =-v3b2(i,je-1,k)
              endif
             endif
             if (nojb3(i,k) .eq. 2) then
c*if alias PROBLEM.eq.advect
c               v2b3(i,je+2,k) = 3.0 * ( v2b3(i,je+1,k) - v2b3(i,je  ,k) )
c       1                      + v2b3(i,je-1,k)
c               v2b3(i,je+3,k) = 3.0 * ( v2b3(i,je+2,k) - v2b3(i,je+1,k) )
c       1                      + v2b3(i,je  ,k)
c*else
               v2b3(i,je+2,k) = v2b3(i,je+1,k)
               v2b3(i,je+3,k) = v2b3(i,je+1,k)
               v3b2(i,je+2,k) = v3b2(i,je+1,k)
               v3b2(i,je+3,k) = v3b2(i,je+1,k)
             endif
             if (nojb3(i,k) .eq. 3) then
               v2b3(i,je+1,k) = emf1ojb(i,k,1)
               v2b3(i,je+2,k) = emf1ojb(i,k,2)
               v2b3(i,je+3,k) = emf1ojb(i,k,3)
               v3b2(i,je+1,k) = 0.0
               v3b2(i,je+2,k) = 0.0
               v3b2(i,je+3,k) = 0.0
             endif
             if (nojb3(i,k) .eq. 4) then
               v2b3(i,je+2,k) = v2b3(i,js+1,k)
               v2b3(i,je+3,k) = v2b3(i,js+2,k)
               v3b2(i,je+2,k) = v3b2(i,js+1,k)
               v3b2(i,je+3,k) = v3b2(i,js+2,k)
             endif
             if (nojb3(i,k) .eq. 5) then
               v2b3(i,je+1,k) = 0.0
               v2b3(i,je+2,k) = v2b3(i,je  ,k)
               v2b3(i,je+3,k) = v2b3(i,je-1,k)
               v3b2(i,je+2,k) = v3b2(i,je  ,k)
               v3b2(i,je+3,k) = v3b2(i,je-1,k)
             endif
c
           enddo
         enddo
#endif /* NO MPI */
c
c-----------------------------------------------------------------------
c------------------------  K - B O U N D A R Y  ------------------------
c-----------------------------------------------------------------------
c
      if(ldimen .lt. 3) return
c
c      Inner k boundary.
c
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1

       if (niks(1).eq.0 .or. niks(1).eq.4) then
         do kk = 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v2b3(1, 1, ks+1+kk), 1, k_slice, n3m, 
     &                   8000+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v3b2(1, 1, ks+1+kk), 1, k_slice, n3m, 
     &                   8100+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v2b3(1, 1, ks-2+kk), 1, k_slice, n3m, 
     &                   7800+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v3b2(1, 1, ks-2+kk), 1, k_slice, n3m, 
     &                   7900+nsub+kk, comm3d , req(nreq), ierr)
         enddo
       else
         do j=js-2,je+3
cdir$ ivdep
           do i=is-2,ie+2
c 
c      Inner k boundary.
c 
             if ( abs(nikb2(i,j)) .eq. 1) then
               v2b3(i,j,ks  ) = 0.0
               v2b3(i,j,ks-1) =-v2b3(i,j,ks+1)
               v2b3(i,j,ks-2) =-v2b3(i,j,ks+2)
               v3b2(i,j,ks  ) = 0.0
               v3b2(i,j,ks-1) =-v3b2(i,j,ks+1)
               v3b2(i,j,ks-2) =-v3b2(i,j,ks+2)
             endif
             if (nikb2(i,j) .eq. 2) then
               v2b3(i,j,ks-1) = v2b3(i,j,ks  )
               v2b3(i,j,ks-2) = v2b3(i,j,ks  )
               v3b2(i,j,ks-1) = v3b2(i,j,ks  )
               v3b2(i,j,ks-2) = v3b2(i,j,ks  )
             endif
             if (nikb2(i,j) .eq. 3) then
               v2b3(i,j,ks  ) = 0.0
               v2b3(i,j,ks-1) = 0.0
               v2b3(i,j,ks-2) = 0.0
               v3b2(i,j,ks  ) =-emf1ikb(i,j,1)
               v3b2(i,j,ks-1) =-emf1ikb(i,j,2)
               v3b2(i,j,ks-2) =-emf1ikb(i,j,3)
             endif
            if (nikb2(i,j) .eq. 5) then
               v2b3(i,j,ks-1) = v2b3(i,j,ks+1)
               v2b3(i,j,ks-2) = v2b3(i,j,ks+2)
               v3b2(i,j,ks  ) = 0.0
               v3b2(i,j,ks-1) = v3b2(i,j,ks+1)
               v3b2(i,j,ks-2) = v3b2(i,j,ks+2)
             endif
c
           enddo
         enddo
       endif
c 
c      Outer k boundary.
c 
       if (noks(1).eq.0 .or. noks(1).eq.4) then
         do kk = 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v2b3(1, 1, ke-1+kk), 1, k_slice, n3p, 
     &                   7800+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v3b2(1, 1, ke-1+kk), 1, k_slice, n3p, 
     &                   7900+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v2b3(1, 1, ke+2+kk), 1, k_slice, n3p, 
     &                   8000+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v3b2(1, 1, ke+2+kk), 1, k_slice, n3p, 
     &                   8100+nsub+kk, comm3d , req(nreq), ierr)
         enddo
       else
         do j=js-2,je+3
cdir$ ivdep
           do i=is-2,ie+2
             if ( abs(nokb2(i,j)) .eq. 1) then
               v2b3(i,j,ke+1) = 0.0
               v2b3(i,j,ke+2) =-v2b3(i,j,ke  )
               v2b3(i,j,ke+3) =-v2b3(i,j,ke-1)
               v3b2(i,j,ke+1) = 0.0
               v3b2(i,j,ke+2) =-v3b2(i,j,ke  )
               v3b2(i,j,ke+3) =-v3b2(i,j,ke-1)
             endif
             if (nokb2(i,j) .eq. 2) then
               v2b3(i,j,ke+2) = v2b3(i,j,ke+1)
               v2b3(i,j,ke+3) = v2b3(i,j,ke+1)
               v3b2(i,j,ke+2) = v3b2(i,j,ke+1)
               v3b2(i,j,ke+3) = v3b2(i,j,ke+1)
             endif
             if (nokb2(i,j) .eq. 3) then
               v2b3(i,j,ke+1) = 0.0
               v2b3(i,j,ke+2) = 0.0
               v2b3(i,j,ke+3) = 0.0
               v3b2(i,j,ke+1) =-emf1okb(i,j,1)
               v3b2(i,j,ke+2) =-emf1okb(i,j,2)
               v3b2(i,j,ke+3) =-emf1okb(i,j,3)
             endif
             if (nokb2(i,j) .eq. 5) then
               v2b3(i,j,ke+2) = v2b3(i,j,ke  )
               v2b3(i,j,ke+3) = v2b3(i,j,ke-1)
               v3b2(i,j,ke+1) = 0.0
               v3b2(i,j,ke+2) = v3b2(i,j,ke  )
               v3b2(i,j,ke+3) = v3b2(i,j,ke-1)
             endif
c 
           enddo
         enddo
       endif
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )

#endif /* MPI */
#ifndef MPI_USED
         do j=js-2,je+3
cdir$ ivdep
           do i=is-2,ie+2
c 
c      Inner k boundary.
c 
             if ( abs(nikb2(i,j)) .eq. 1) then
               v2b3(i,j,ks  ) = 0.0
               v2b3(i,j,ks-1) =-v2b3(i,j,ks+1)
               v2b3(i,j,ks-2) =-v2b3(i,j,ks+2)
               v3b2(i,j,ks  ) = 0.0
               v3b2(i,j,ks-1) =-v3b2(i,j,ks+1)
               v3b2(i,j,ks-2) =-v3b2(i,j,ks+2)
             endif
             if (nikb2(i,j) .eq. 2) then
               v2b3(i,j,ks-1) = v2b3(i,j,ks  )
               v2b3(i,j,ks-2) = v2b3(i,j,ks  )
               v3b2(i,j,ks-1) = v3b2(i,j,ks  )
               v3b2(i,j,ks-2) = v3b2(i,j,ks  )
             endif
             if (nikb2(i,j) .eq. 3) then
               v2b3(i,j,ks  ) = 0.0
               v2b3(i,j,ks-1) = 0.0
               v2b3(i,j,ks-2) = 0.0
               v3b2(i,j,ks  ) =-emf1ikb(i,j,1)
               v3b2(i,j,ks-1) =-emf1ikb(i,j,2)
               v3b2(i,j,ks-2) =-emf1ikb(i,j,3)
             endif
             if (nikb2(i,j) .eq. 4) then
               v2b3(i,j,ks-1) = v2b3(i,j,ke  )
               v2b3(i,j,ks-2) = v2b3(i,j,ke-1)
               v3b2(i,j,ks-1) = v3b2(i,j,ke  )
               v3b2(i,j,ks-2) = v3b2(i,j,ke-1)
             endif
            if (nikb2(i,j) .eq. 5) then
               v2b3(i,j,ks-1) = v2b3(i,j,ks+1)
               v2b3(i,j,ks-2) = v2b3(i,j,ks+2)
               v3b2(i,j,ks  ) = 0.0
               v3b2(i,j,ks-1) = v3b2(i,j,ks+1)
               v3b2(i,j,ks-2) = v3b2(i,j,ks+2)
             endif
c
           enddo
         enddo
c 
c      Outer k boundary.
c 
         do j=js-2,je+3
cdir$ ivdep
           do i=is-2,ie+2
             if ( abs(nokb2(i,j)) .eq. 1) then
               v2b3(i,j,ke+1) = 0.0
               v2b3(i,j,ke+2) =-v2b3(i,j,ke  )
               v2b3(i,j,ke+3) =-v2b3(i,j,ke-1)
               v3b2(i,j,ke+1) = 0.0
               v3b2(i,j,ke+2) =-v3b2(i,j,ke  )
               v3b2(i,j,ke+3) =-v3b2(i,j,ke-1)
             endif
             if (nokb2(i,j) .eq. 2) then
               v2b3(i,j,ke+2) = v2b3(i,j,ke+1)
               v2b3(i,j,ke+3) = v2b3(i,j,ke+1)
               v3b2(i,j,ke+2) = v3b2(i,j,ke+1)
               v3b2(i,j,ke+3) = v3b2(i,j,ke+1)
             endif
             if (nokb2(i,j) .eq. 3) then
               v2b3(i,j,ke+1) = 0.0
               v2b3(i,j,ke+2) = 0.0
               v2b3(i,j,ke+3) = 0.0
               v3b2(i,j,ke+1) =-emf1okb(i,j,1)
               v3b2(i,j,ke+2) =-emf1okb(i,j,2)
               v3b2(i,j,ke+3) =-emf1okb(i,j,3)
             endif
             if (nokb2(i,j) .eq. 4) then
               v2b3(i,j,ke+2) = v2b3(i,j,ks+1)
               v2b3(i,j,ke+3) = v2b3(i,j,ks+2)
               v3b2(i,j,ke+2) = v3b2(i,j,ks+1)
               v3b2(i,j,ke+3) = v3b2(i,j,ks+2)
             endif
             if (nokb2(i,j) .eq. 5) then
               v2b3(i,j,ke+2) = v2b3(i,j,ke  )
               v2b3(i,j,ke+3) = v2b3(i,j,ke-1)
               v3b2(i,j,ke+1) = 0.0
               v3b2(i,j,ke+2) = v3b2(i,j,ke  )
               v3b2(i,j,ke+3) = v3b2(i,j,ke-1)
             endif
c 
           enddo
         enddo
#endif /* NO MPI */
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////              B V A L E M F 1              \\\\\\\\\\
c
c=======================================================================
c
c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              B V A L E M F 2              \\\\\\\\\\
c
c=======================================================================
c
       subroutine bvalemf2 ( v3b1, v1b3 )
c
c    dac:zeus3d.bvalemf2 <-------------- boundary values for 2-emf terms
c                                                         february, 1990
c
c    written by: David Clarke, February 1990
c    modified 1: September, 1990 by David Clarke; moved magnetic fields
c                to face-centres.
c    modified 2: minimal rewrite for ZEUS-MP by M-MML 10.3.98
c    modified 3: added restrictions for symmetry about the J and K
c                axes; John Hayes, 10/2005.
c
c  PURPOSE: This routine sets boundary values for the two terms in the
c  2-emf (centred on the 2-edges).  The active zones for "emf2" are:
c
c    i = is to ie+1;  j = js to je;  k = ks tp ke+1
c
c  In order to update both the active and ghost zones of the 3- and
c  1-magnetic field components, all edges in the boundary regions are
c  required.  This gives a complete grid of values for "emf2".  Thus,
c  the ranges for the boundary values are:
c
c    j-boundaries:   i = is  , ie+1                    k = ks  , ke+1
c    k-boundaries:   i = is  , ie+1   j = js-2, je+2
c    i-boundaries:                    j = js-2, je+2   k = ks-2, ke+3
c
c  Note that the boundary values must be set even if js > jsmn, etc.
c  because the emfs are stored in worker arrays and it is likely that
c  the boundary values have been overwritten.

c  See comments in BVALD.
c
c  EXTERNALS: [NONE}
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
      integer  :: i, j, k, ii, jj, kk
c
      real(rl) :: v3b1(  in,  jn,  kn), v1b3(  in,  jn,  kn)
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c------------------------  J - B O U N D A R Y  ------------------------
c-----------------------------------------------------------------------
c
      if(ldimen .eq. 1) go to 111
c
c      Inner j boundary.
c
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1
c
       if (nijs(1).eq.0 .or. nijs(1).eq.4) then
         do jj= 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v3b1(1, js  +jj, 1), 1, j_slice, n2m, 
     &                   8400+nsub+jj, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v1b3(1, js  +jj, 1), 1, j_slice, n2m, 
     &                   8500+nsub+jj, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v3b1(1, js-2+jj, 1), 1, j_slice, n2m, 
     &                   8200+nsub+jj, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v1b3(1, js-2+jj, 1), 1, j_slice, n2m, 
     &                   8300+nsub+jj, comm3d , req(nreq), ierr)
         enddo
       else
         do k=ks,ke+1
cdir$ ivdep
           do i=is,ie+1
             if (nijb31(i,k) .eq. 1) then
               v3b1(i,js-1,k) = v3b1(i,js  ,k)
               v3b1(i,js-2,k) = v3b1(i,js+1,k)
               v1b3(i,js-1,k) = v1b3(i,js  ,k)
               v1b3(i,js-2,k) = v1b3(i,js+1,k)
             endif
             if (nijb31(i,k) .eq.-1) then
              if(lgeom .eq. 1) then
               v3b1(i,js-1,k) = v3b1(i,js  ,k)
               v3b1(i,js-2,k) = v3b1(i,js+1,k)
               v1b3(i,js-1,k) = v1b3(i,js  ,k)
               v1b3(i,js-2,k) = v1b3(i,js+1,k)
              else ! lgeom
               v3b1(i,js-1,k) =-v3b1(i,js  ,k)
               v3b1(i,js-2,k) =-v3b1(i,js+1,k)
               v1b3(i,js-1,k) =-v1b3(i,js  ,k)
               v1b3(i,js-2,k) =-v1b3(i,js+1,k)
              endif ! lgeom
             endif
             if (nijb31(i,k) .eq. 2) then
               v3b1(i,js-1,k) = v3b1(i,js  ,k)
               v3b1(i,js-2,k) = v3b1(i,js  ,k)
               v1b3(i,js-1,k) = v1b3(i,js  ,k)
               v1b3(i,js-2,k) = v1b3(i,js  ,k)
             endif
             if (nijb31(i,k) .eq. 3) then
               v3b1(i,js-1,k) = emf2ijb(i,k,1)
               v3b1(i,js-2,k) = emf2ijb(i,k,2)
               v1b3(i,js-1,k) = 0.0
               v1b3(i,js-2,k) = 0.0
             endif
             if (nijb31(i,k) .eq. 5) then
               v3b1(i,js-1,k) =-v3b1(i,js  ,k)
               v3b1(i,js-2,k) =-v3b1(i,js+1,k)
               v1b3(i,js-1,k) =-v1b3(i,js  ,k)
               v1b3(i,js-2,k) =-v1b3(i,js+1,k)
             endif
           enddo
         enddo
       endif
c 
c      Outer j boundary.
c 
c
       if (nojs(1).eq.0 .or. nojs(1).eq.4) then
         do jj = 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v3b1(1, je-1+jj, 1), 1, j_slice, n2p, 
     &           8200+nsub+jj, comm3d , req(nreq), ierr) 
           nreq = nreq + 1
           call MPI_ISEND(v1b3(1, je-1+jj, 1), 1, j_slice, n2p, 
     &           8300+nsub+jj, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v3b1(1, je+1+jj, 1), 1, j_slice, n2p, 
     &           8400+nsub+jj, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v1b3(1, je+1+jj, 1), 1, j_slice, n2p, 
     &           8500+nsub+jj, comm3d , req(nreq), ierr)
         enddo
       else
         do k=ks,ke+1
cdir$ ivdep
           do i=is,ie+1
             if (nojb31(i,k) .eq. 1) then
               v3b1(i,je+1,k) = v3b1(i,je  ,k)
               v3b1(i,je+2,k) = v3b1(i,je-1,k)
               v1b3(i,je+1,k) = v1b3(i,je  ,k)
               v1b3(i,je+2,k) = v1b3(i,je-1,k)
             endif
             if (nojb31(i,k) .eq.-1) then
              if(lgeom .eq. 3) then
               v3b1(i,je+1,k) =-v3b1(i,je  ,k)
               v3b1(i,je+2,k) =-v3b1(i,je-1,k)
               v1b3(i,je+1,k) =-v1b3(i,je  ,k)
               v1b3(i,je+2,k) =-v1b3(i,je-1,k)
              else
               v3b1(i,je+1,k) = v3b1(i,je  ,k)
               v3b1(i,je+2,k) = v3b1(i,je-1,k)
               v1b3(i,je+1,k) = v1b3(i,je  ,k)
               v1b3(i,je+2,k) = v1b3(i,je-1,k)
              endif ! lgeom
             endif
             if (nojb31(i,k) .eq. 2) then
               v3b1(i,je+1,k) = v3b1(i,je  ,k)
               v3b1(i,je+2,k) = v3b1(i,je  ,k)
               v1b3(i,je+1,k) = v1b3(i,je  ,k)
               v1b3(i,je+2,k) = v1b3(i,je  ,k)
             endif
             if (nojb31(i,k) .eq. 3) then
               v3b1(i,je+1,k) = emf2ojb(i,k,1)
               v3b1(i,je+2,k) = emf2ojb(i,k,2)
               v1b3(i,je+1,k) = 0.0
               v1b3(i,je+2,k) = 0.0
             endif
             if (nojb31(i,k) .eq. 5) then
               v3b1(i,je+1,k) =-v3b1(i,je  ,k)
               v3b1(i,je+2,k) =-v3b1(i,je-1,k)
               v1b3(i,je+1,k) =-v1b3(i,je  ,k)
               v1b3(i,je+2,k) =-v1b3(i,je-1,k)
             endif
c 
           enddo
          enddo
       endif
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif /* MPI */
#ifndef MPI_USED
         do k=ks,ke+1
cdir$ ivdep
           do i=is,ie+1
             if (nijb31(i,k) .eq. 1) then
               v3b1(i,js-1,k) = v3b1(i,js  ,k)
               v3b1(i,js-2,k) = v3b1(i,js+1,k)
               v1b3(i,js-1,k) = v1b3(i,js  ,k)
               v1b3(i,js-2,k) = v1b3(i,js+1,k)
             endif
             if (nijb31(i,k) .eq.-1) then
              if(lgeom .eq. 1) then
               v3b1(i,js-1,k) = v3b1(i,js  ,k)
               v3b1(i,js-2,k) = v3b1(i,js+1,k)
               v1b3(i,js-1,k) = v1b3(i,js  ,k)
               v1b3(i,js-2,k) = v1b3(i,js+1,k)
              else ! lgeom
               v3b1(i,js-1,k) =-v3b1(i,js  ,k)
               v3b1(i,js-2,k) =-v3b1(i,js+1,k)
               v1b3(i,js-1,k) =-v1b3(i,js  ,k)
               v1b3(i,js-2,k) =-v1b3(i,js+1,k)
              endif ! lgeom
             endif
             if (nijb31(i,k) .eq. 2) then
               v3b1(i,js-1,k) = v3b1(i,js  ,k)
               v3b1(i,js-2,k) = v3b1(i,js  ,k)
               v1b3(i,js-1,k) = v1b3(i,js  ,k)
               v1b3(i,js-2,k) = v1b3(i,js  ,k)
             endif
             if (nijb31(i,k) .eq. 3) then
               v3b1(i,js-1,k) = emf2ijb(i,k,1)
               v3b1(i,js-2,k) = emf2ijb(i,k,2)
               v1b3(i,js-1,k) = 0.0
               v1b3(i,js-2,k) = 0.0
             endif
             if (nijb31(i,k) .eq. 4) then
               v3b1(i,js-1,k) = v3b1(i,je  ,k)
               v3b1(i,js-2,k) = v3b1(i,je-1,k)
               v1b3(i,js-1,k) = v1b3(i,je  ,k)
               v1b3(i,js-2,k) = v1b3(i,je-1,k)
             endif
             if (nijb31(i,k) .eq. 5) then
               v3b1(i,js-1,k) =-v3b1(i,js  ,k)
               v3b1(i,js-2,k) =-v3b1(i,js+1,k)
               v1b3(i,js-1,k) =-v1b3(i,js  ,k)
               v1b3(i,js-2,k) =-v1b3(i,js+1,k)
             endif
           enddo
         enddo
c 
c      Outer j boundary.
c 
         do k=ks,ke+1
cdir$ ivdep
           do i=is,ie+1
             if (nojb31(i,k) .eq. 1) then
               v3b1(i,je+1,k) = v3b1(i,je  ,k)
               v3b1(i,je+2,k) = v3b1(i,je-1,k)
               v1b3(i,je+1,k) = v1b3(i,je  ,k)
               v1b3(i,je+2,k) = v1b3(i,je-1,k)
             endif
             if (nojb31(i,k) .eq.-1) then
              if(lgeom .eq. 3) then
               v3b1(i,je+1,k) =-v3b1(i,je  ,k)
               v3b1(i,je+2,k) =-v3b1(i,je-1,k)
               v1b3(i,je+1,k) =-v1b3(i,je  ,k)
               v1b3(i,je+2,k) =-v1b3(i,je-1,k)
              else
               v3b1(i,je+1,k) = v3b1(i,je  ,k)
               v3b1(i,je+2,k) = v3b1(i,je-1,k)
               v1b3(i,je+1,k) = v1b3(i,je  ,k)
               v1b3(i,je+2,k) = v1b3(i,je-1,k)
              endif ! lgeom
             endif
             if (nojb31(i,k) .eq. 2) then
               v3b1(i,je+1,k) = v3b1(i,je  ,k)
               v3b1(i,je+2,k) = v3b1(i,je  ,k)
               v1b3(i,je+1,k) = v1b3(i,je  ,k)
               v1b3(i,je+2,k) = v1b3(i,je  ,k)
             endif
             if (nojb31(i,k) .eq. 3) then
               v3b1(i,je+1,k) = emf2ojb(i,k,1)
               v3b1(i,je+2,k) = emf2ojb(i,k,2)
               v1b3(i,je+1,k) = 0.0
               v1b3(i,je+2,k) = 0.0
             endif
             if (nojb31(i,k) .eq. 4) then
               v3b1(i,je+1,k) = v3b1(i,js  ,k)
               v3b1(i,je+2,k) = v3b1(i,js+1,k)
               v1b3(i,je+1,k) = v1b3(i,js  ,k)
               v1b3(i,je+2,k) = v1b3(i,js+1,k)
             endif
             if (nojb31(i,k) .eq. 5) then
               v3b1(i,je+1,k) =-v3b1(i,je  ,k)
               v3b1(i,je+2,k) =-v3b1(i,je-1,k)
               v1b3(i,je+1,k) =-v1b3(i,je  ,k)
               v1b3(i,je+2,k) =-v1b3(i,je-1,k)
             endif
c 
           enddo
          enddo
#endif /* NO MPI */
c
c
c-----------------------------------------------------------------------
c------------------------  K - B O U N D A R Y  ------------------------
c-----------------------------------------------------------------------
c
       if(ldimen .ne. 3) go to 111
c
c      Inner k boundary.
c
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1
       if (niks(1).eq.0 .or. niks(1).eq.4) then
         do kk = 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v3b1(1, 1, ks+1+kk), 1, k_slice, n3m, 
     &                   8800+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v1b3(1, 1, ks+1+kk), 1, k_slice, n3m, 
     &                   8900+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v3b1(1, 1, ks-2+kk), 1, k_slice, n3m, 
     &                   8600+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v1b3(1, 1, ks-2+kk), 1, k_slice, n3m, 
     &                   8700+nsub+kk, comm3d , req(nreq), ierr)
         enddo
       else
         do j=js-2,je+2
cdir$ ivdep
           do i=is,ie+1
             if ( abs(nikb1(i,j)) .eq. 1) then
               v3b1(i,j,ks  ) = 0.0
               v3b1(i,j,ks-1) =-v3b1(i,j,ks+1)
               v3b1(i,j,ks-2) =-v3b1(i,j,ks+2)
               v1b3(i,j,ks  ) = 0.0
               v1b3(i,j,ks-1) =-v1b3(i,j,ks+1)
               v1b3(i,j,ks-2) =-v1b3(i,j,ks+2)
             endif
             if (nikb1(i,j) .eq. 2) then
               v3b1(i,j,ks-1) = v3b1(i,j,ks  )
               v3b1(i,j,ks-2) = v3b1(i,j,ks  )
               v1b3(i,j,ks-1) = v1b3(i,j,ks  )
               v1b3(i,j,ks-2) = v1b3(i,j,ks  )
             endif
             if (nikb1(i,j) .eq. 3) then
               v3b1(i,j,ks  ) = emf2ikb(i,j,1)
               v3b1(i,j,ks-1) = emf2ikb(i,j,2)
               v3b1(i,j,ks-2) = emf2ikb(i,j,3)
               v1b3(i,j,ks  ) = 0.0
               v1b3(i,j,ks-1) = 0.0
               v1b3(i,j,ks-2) = 0.0
             endif
             if (nikb1(i,j) .eq. 5) then
               v3b1(i,j,ks  ) = 0.0
               v3b1(i,j,ks-1) = v3b1(i,j,ks+1)
               v3b1(i,j,ks-2) = v3b1(i,j,ks+2)
               v1b3(i,j,ks-1) = v1b3(i,j,ks+1)
               v1b3(i,j,ks-2) = v1b3(i,j,ks+2)
             endif
           enddo
         enddo
       endif
c 
c      Outer k boundary.
c 
       if (noks(1).eq.0 .or. noks(1).eq.4) then
         do kk = 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v3b1(1, 1, ke-1+kk), 1, k_slice, n3p, 
     &                   8600+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v1b3(1, 1, ke-1+kk), 1, k_slice, n3p, 
     &                   8700+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v3b1(1, 1, ke+2+kk), 1, k_slice, n3p, 
     &                   8800+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v1b3(1, 1, ke+2+kk), 1, k_slice, n3p, 
     &                   8900+nsub+kk, comm3d , req(nreq), ierr)
         enddo
       else
         do j=js-2,je+2
cdir$ ivdep
           do i=is,ie+1
             if ( abs(nokb1(i,j)) .eq. 1) then
               v3b1(i,j,ke+1) = 0.0
               v3b1(i,j,ke+2) =-v3b1(i,j,ke  )
               v3b1(i,j,ke+3) =-v3b1(i,j,ke-1)
               v1b3(i,j,ke+1) = 0.0
               v1b3(i,j,ke+2) =-v1b3(i,j,ke  )
               v1b3(i,j,ke+3) =-v1b3(i,j,ke-1)
             endif
             if (nokb1(i,j) .eq. 2) then
               v3b1(i,j,ke+2) = v3b1(i,j,ke+1)
               v3b1(i,j,ke+3) = v3b1(i,j,ke+1)
               v1b3(i,j,ke+2) = v1b3(i,j,ke+1)
               v1b3(i,j,ke+3) = v1b3(i,j,ke+1)
             endif
             if (nokb1(i,j) .eq. 3) then
               v3b1(i,j,ke+1) = emf2okb(i,j,1)
               v3b1(i,j,ke+2) = emf2okb(i,j,2)
               v3b1(i,j,ke+3) = emf2okb(i,j,3)
               v1b3(i,j,ke+1) = 0.0
               v1b3(i,j,ke+2) = 0.0
               v1b3(i,j,ke+3) = 0.0
             endif
             if (nokb1(i,j) .eq. 5) then
               v3b1(i,j,ke+1) = 0.0
               v3b1(i,j,ke+2) = v3b1(i,j,ke  )
               v3b1(i,j,ke+3) = v3b1(i,j,ke-1)
               v1b3(i,j,ke+2) = v1b3(i,j,ke  )
               v1b3(i,j,ke+3) = v1b3(i,j,ke-1)
             endif
c 
           enddo
         enddo
       endif
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif /* MPI */
#ifndef MPI_USED
         do j=js-2,je+2
cdir$ ivdep
           do i=is,ie+1
             if ( abs(nikb1(i,j)) .eq. 1) then
               v3b1(i,j,ks  ) = 0.0
               v3b1(i,j,ks-1) =-v3b1(i,j,ks+1)
               v3b1(i,j,ks-2) =-v3b1(i,j,ks+2)
               v1b3(i,j,ks  ) = 0.0
               v1b3(i,j,ks-1) =-v1b3(i,j,ks+1)
               v1b3(i,j,ks-2) =-v1b3(i,j,ks+2)
             endif
             if (nikb1(i,j) .eq. 2) then
               v3b1(i,j,ks-1) = v3b1(i,j,ks  )
               v3b1(i,j,ks-2) = v3b1(i,j,ks  )
               v1b3(i,j,ks-1) = v1b3(i,j,ks  )
               v1b3(i,j,ks-2) = v1b3(i,j,ks  )
             endif
             if (nikb1(i,j) .eq. 3) then
               v3b1(i,j,ks  ) = emf2ikb(i,j,1)
               v3b1(i,j,ks-1) = emf2ikb(i,j,2)
               v3b1(i,j,ks-2) = emf2ikb(i,j,3)
               v1b3(i,j,ks  ) = 0.0
               v1b3(i,j,ks-1) = 0.0
               v1b3(i,j,ks-2) = 0.0
             endif
             if (nikb1(i,j) .eq. 4) then
               v3b1(i,j,ks-1) = v3b1(i,j,ke  )
               v3b1(i,j,ks-2) = v3b1(i,j,ke-1)
               v1b3(i,j,ks-1) = v1b3(i,j,ke  )
               v1b3(i,j,ks-2) = v1b3(i,j,ke-1)
             endif
             if (nikb1(i,j) .eq. 5) then
               v3b1(i,j,ks  ) = 0.0
               v3b1(i,j,ks-1) = v3b1(i,j,ks+1)
               v3b1(i,j,ks-2) = v3b1(i,j,ks+2)
               v1b3(i,j,ks-1) = v1b3(i,j,ks+1)
               v1b3(i,j,ks-2) = v1b3(i,j,ks+2)
             endif
           enddo
         enddo
c 
c      Outer k boundary.
c 
         do j=js-2,je+2
cdir$ ivdep
           do i=is,ie+1
             if ( abs(nokb1(i,j)) .eq. 1) then
               v3b1(i,j,ke+1) = 0.0
               v3b1(i,j,ke+2) =-v3b1(i,j,ke  )
               v3b1(i,j,ke+3) =-v3b1(i,j,ke-1)
               v1b3(i,j,ke+1) = 0.0
               v1b3(i,j,ke+2) =-v1b3(i,j,ke  )
               v1b3(i,j,ke+3) =-v1b3(i,j,ke-1)
             endif
             if (nokb1(i,j) .eq. 2) then
               v3b1(i,j,ke+2) = v3b1(i,j,ke+1)
               v3b1(i,j,ke+3) = v3b1(i,j,ke+1)
               v1b3(i,j,ke+2) = v1b3(i,j,ke+1)
               v1b3(i,j,ke+3) = v1b3(i,j,ke+1)
             endif
             if (nokb1(i,j) .eq. 3) then
               v3b1(i,j,ke+1) = emf2okb(i,j,1)
               v3b1(i,j,ke+2) = emf2okb(i,j,2)
               v3b1(i,j,ke+3) = emf2okb(i,j,3)
               v1b3(i,j,ke+1) = 0.0
               v1b3(i,j,ke+2) = 0.0
               v1b3(i,j,ke+3) = 0.0
             endif
             if (nokb1(i,j) .eq. 4) then
               v3b1(i,j,ke+2) = v3b1(i,j,ks+1)
               v3b1(i,j,ke+3) = v3b1(i,j,ks+2)
               v1b3(i,j,ke+2) = v1b3(i,j,ks+1)
               v1b3(i,j,ke+3) = v1b3(i,j,ks+2)
             endif
             if (nokb1(i,j) .eq. 5) then
               v3b1(i,j,ke+1) = 0.0
               v3b1(i,j,ke+2) = v3b1(i,j,ke  )
               v3b1(i,j,ke+3) = v3b1(i,j,ke-1)
               v1b3(i,j,ke+2) = v1b3(i,j,ke  )
               v1b3(i,j,ke+3) = v1b3(i,j,ke-1)
             endif
c 
           enddo
         enddo
#endif /* NO MPI */
c
c-----------------------------------------------------------------------
c------------------------  I - B O U N D A R Y  ------------------------
c-----------------------------------------------------------------------
c
111   CONTINUE
c
c      Inner i boundary.
c
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1

       if (niis(1).eq.0 .or. niis(1).eq.4) then
         do ii = 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v3b1(is+1+ii, 1, 1), 1, i_slice, n1m, 
     &                   9200+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v1b3(is+1+ii, 1, 1), 1, i_slice, n1m, 
     &                   9300+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v3b1(is-2+ii, 1, 1), 1, i_slice, n1m, 
     &                   9000+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v1b3(is-2+ii, 1, 1), 1, i_slice, n1m, 
     &                   9100+nsub+ii, comm3d , req(nreq), ierr)
         enddo
       else
         do k=ks-2,ke+3
cdir$ ivdep
           do j=js-2,je+2
             if (niib3(j,k) .eq. 1) then
               v3b1(is  ,j,k) = 0.0
               v3b1(is-1,j,k) =-v3b1(is+1,j,k)
               v3b1(is-2,j,k) =-v3b1(is+2,j,k)
               v1b3(is  ,j,k) = 0.0
               v1b3(is-1,j,k) =-v1b3(is+1,j,k)
               v1b3(is-2,j,k) =-v1b3(is+2,j,k)
             endif
             if (niib3(j,k) .eq.-1) then
              if(lgeom .eq. 3) then
               v3b1(is  ,j,k) = 0.0
               v3b1(is-1,j,k) = v3b1(is+1,j,k)
               v3b1(is-2,j,k) = v3b1(is+2,j,k)
               v1b3(is  ,j,k) = 0.0
               v1b3(is-1,j,k) = v1b3(is+1,j,k)
               v1b3(is-2,j,k) = v1b3(is+2,j,k)
              else
               v3b1(is  ,j,k) = 0.0
               v3b1(is-1,j,k) =-v3b1(is+1,j,k)
               v3b1(is-2,j,k) =-v3b1(is+2,j,k)
               v1b3(is  ,j,k) = 0.0
               v1b3(is-1,j,k) =-v1b3(is+1,j,k)
               v1b3(is-2,j,k) =-v1b3(is+2,j,k)
              endif ! lgeom
             endif
             if (niib3(j,k) .eq. 2) then
               v3b1(is-1,j,k) = v3b1(is  ,j,k)
               v3b1(is-2,j,k) = v3b1(is  ,j,k)
               v1b3(is-1,j,k) = v1b3(is  ,j,k)
               v1b3(is-2,j,k) = v1b3(is  ,j,k)
             endif
             if (niib3(j,k) .eq. 3) then
               v3b1(is  ,j,k) = 0.0
               v3b1(is-1,j,k) = 0.0
               v3b1(is-2,j,k) = 0.0
               v1b3(is  ,j,k) =-emf2iib(j,k,1)
               v1b3(is-1,j,k) =-emf2iib(j,k,2)
               v1b3(is-2,j,k) =-emf2iib(j,k,3)
             endif
             if (niib3(j,k) .eq. 5) then
               v3b1(is-1,j,k) = v3b1(is+1,j,k)
               v3b1(is-2,j,k) = v3b1(is+2,j,k)
               v1b3(is  ,j,k) = 0.0
               v1b3(is-1,j,k) = v1b3(is+1,j,k)
               v1b3(is-2,j,k) = v1b3(is+2,j,k)
             endif
           enddo
         enddo
       endif
c 
c      Outer i boundary.
c 
       if (nois(1).eq.0 .or. nois(1).eq.4) then
         do ii = 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v3b1(ie-1+ii, 1, 1), 1, i_slice, n1p, 
     &                   9000+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v1b3(ie-1+ii, 1, 1), 1, i_slice, n1p, 
     &                   9100+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v3b1(ie+2+ii, 1, 1), 1, i_slice, n1p, 
     &                   9200+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v1b3(ie+2+ii, 1, 1), 1, i_slice, n1p, 
     &                   9300+nsub+ii, comm3d , req(nreq), ierr)
         enddo
       else
         do k=ks-2,ke+3
cdir$ ivdep
           do j=js-2,je+2
             if ( abs(noib3(j,k)) .eq. 1) then
               v3b1(ie+1,j,k) = 0.0
               v3b1(ie+2,j,k) =-v3b1(ie  ,j,k)
               v3b1(ie+3,j,k) =-v3b1(ie-1,j,k)
               v1b3(ie+1,j,k) = 0.0
               v1b3(ie+2,j,k) =-v1b3(ie  ,j,k)
               v1b3(ie+3,j,k) =-v1b3(ie-1,j,k)
             endif
             if (noib3(j,k) .eq. 2) then
               v3b1(ie+2,j,k) = v3b1(ie+1,j,k)
               v3b1(ie+3,j,k) = v3b1(ie+1,j,k)
c*if alias PROBLEM.eq.advect
c               v1b3(ie+2,j,k) = 3.0 * ( v1b3(ie+1,j,k) - v1b3(ie  ,j,k) )
c       1                      + v1b3(ie-1,j,k)
c               v1b3(ie+3,j,k) = 3.0 * ( v1b3(ie+2,j,k) - v1b3(ie+1,j,k) )
c       1                      + v1b3(ie  ,j,k)
c*else
               v1b3(ie+2,j,k) = v1b3(ie+1,j,k)
               v1b3(ie+3,j,k) = v1b3(ie+1,j,k)
c*endif
             endif
             if (noib3(j,k) .eq. 3) then
               v3b1(ie+1,j,k) = 0.0
               v3b1(ie+2,j,k) = 0.0
               v3b1(ie+3,j,k) = 0.0
               v1b3(ie+1,j,k) =-emf2oib(j,k,1)
               v1b3(ie+2,j,k) =-emf2oib(j,k,2)
               v1b3(ie+3,j,k) =-emf2oib(j,k,3)
             endif
             if (noib3(j,k) .eq. 5) then
               v3b1(ie+2,j,k) = v3b1(ie  ,j,k)
               v3b1(ie+3,j,k) = v3b1(ie-1,j,k)
               v1b3(ie+1,j,k) = 0.0
               v1b3(ie+2,j,k) = v1b3(ie  ,j,k)
               v1b3(ie+3,j,k) = v1b3(ie-1,j,k)
             endif
c
           enddo
         enddo
c
       endif
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif /* MPI */
#ifndef MPI_USED
         do k=ks-2,ke+3
cdir$ ivdep
           do j=js-2,je+2
             if (niib3(j,k) .eq. 1) then
               v3b1(is  ,j,k) = 0.0
               v3b1(is-1,j,k) =-v3b1(is+1,j,k)
               v3b1(is-2,j,k) =-v3b1(is+2,j,k)
               v1b3(is  ,j,k) = 0.0
               v1b3(is-1,j,k) =-v1b3(is+1,j,k)
               v1b3(is-2,j,k) =-v1b3(is+2,j,k)
             endif
             if (niib3(j,k) .eq.-1) then
              if(lgeom .eq. 3) then
               v3b1(is  ,j,k) = 0.0
               v3b1(is-1,j,k) = v3b1(is+1,j,k)
               v3b1(is-2,j,k) = v3b1(is+2,j,k)
               v1b3(is  ,j,k) = 0.0
               v1b3(is-1,j,k) = v1b3(is+1,j,k)
               v1b3(is-2,j,k) = v1b3(is+2,j,k)
              else
               v3b1(is  ,j,k) = 0.0
               v3b1(is-1,j,k) =-v3b1(is+1,j,k)
               v3b1(is-2,j,k) =-v3b1(is+2,j,k)
               v1b3(is  ,j,k) = 0.0
               v1b3(is-1,j,k) =-v1b3(is+1,j,k)
               v1b3(is-2,j,k) =-v1b3(is+2,j,k)
              endif ! lgeom
             endif
             if (niib3(j,k) .eq. 2) then
               v3b1(is-1,j,k) = v3b1(is  ,j,k)
               v3b1(is-2,j,k) = v3b1(is  ,j,k)
               v1b3(is-1,j,k) = v1b3(is  ,j,k)
               v1b3(is-2,j,k) = v1b3(is  ,j,k)
             endif
             if (niib3(j,k) .eq. 3) then
               v3b1(is  ,j,k) = 0.0
               v3b1(is-1,j,k) = 0.0
               v3b1(is-2,j,k) = 0.0
               v1b3(is  ,j,k) =-emf2iib(j,k,1)
               v1b3(is-1,j,k) =-emf2iib(j,k,2)
               v1b3(is-2,j,k) =-emf2iib(j,k,3)
             endif
             if (niib3(j,k) .eq. 4) then
               v3b1(is-1,j,k) = v3b1(ie  ,j,k)
               v3b1(is-2,j,k) = v3b1(ie-1,j,k)
               v1b3(is-1,j,k) = v1b3(ie  ,j,k)
               v1b3(is-2,j,k) = v1b3(ie-1,j,k)
             endif
             if (niib3(j,k) .eq. 5) then
               v3b1(is-1,j,k) = v3b1(is+1,j,k)
               v3b1(is-2,j,k) = v3b1(is+2,j,k)
               v1b3(is  ,j,k) = 0.0
               v1b3(is-1,j,k) = v1b3(is+1,j,k)
               v1b3(is-2,j,k) = v1b3(is+2,j,k)
             endif
           enddo
         enddo
c 
c      Outer i boundary.
c 
         do k=ks-2,ke+3
cdir$ ivdep
           do j=js-2,je+2
             if ( abs(noib3(j,k)) .eq. 1) then
               v3b1(ie+1,j,k) = 0.0
               v3b1(ie+2,j,k) =-v3b1(ie  ,j,k)
               v3b1(ie+3,j,k) =-v3b1(ie-1,j,k)
               v1b3(ie+1,j,k) = 0.0
               v1b3(ie+2,j,k) =-v1b3(ie  ,j,k)
               v1b3(ie+3,j,k) =-v1b3(ie-1,j,k)
             endif
             if (noib3(j,k) .eq. 2) then
               v3b1(ie+2,j,k) = v3b1(ie+1,j,k)
               v3b1(ie+3,j,k) = v3b1(ie+1,j,k)
c*if alias PROBLEM.eq.advect
c               v1b3(ie+2,j,k) = 3.0 * ( v1b3(ie+1,j,k) - v1b3(ie  ,j,k) )
c       1                      + v1b3(ie-1,j,k)
c               v1b3(ie+3,j,k) = 3.0 * ( v1b3(ie+2,j,k) - v1b3(ie+1,j,k) )
c       1                      + v1b3(ie  ,j,k)
c*else
               v1b3(ie+2,j,k) = v1b3(ie+1,j,k)
               v1b3(ie+3,j,k) = v1b3(ie+1,j,k)
c*endif
             endif
             if (noib3(j,k) .eq. 3) then
               v3b1(ie+1,j,k) = 0.0
               v3b1(ie+2,j,k) = 0.0
               v3b1(ie+3,j,k) = 0.0
               v1b3(ie+1,j,k) =-emf2oib(j,k,1)
               v1b3(ie+2,j,k) =-emf2oib(j,k,2)
               v1b3(ie+3,j,k) =-emf2oib(j,k,3)
             endif
             if (noib3(j,k) .eq. 4) then
               v3b1(ie+2,j,k) = v3b1(is+1,j,k)
               v3b1(ie+3,j,k) = v3b1(is+2,j,k)
               v1b3(ie+2,j,k) = v1b3(is+1,j,k)
               v1b3(ie+3,j,k) = v1b3(is+2,j,k)
             endif
             if (noib3(j,k) .eq. 5) then
               v3b1(ie+2,j,k) = v3b1(ie  ,j,k)
               v3b1(ie+3,j,k) = v3b1(ie-1,j,k)
               v1b3(ie+1,j,k) = 0.0
               v1b3(ie+2,j,k) = v1b3(ie  ,j,k)
               v1b3(ie+3,j,k) = v1b3(ie-1,j,k)
             endif
           enddo
         enddo
#endif /* NO MPI */
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////              B V A L E M F 2              \\\\\\\\\\
c
c=======================================================================
c
c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              B V A L E M F 3              \\\\\\\\\\
c
c=======================================================================
c
       subroutine bvalemf3 ( v1b2, v2b1 )
c
c    dac:zeus3d.bvalemf3 <-------------- boundary values for 3-emf terms
c                                                         february, 1990
c
c    written by: David Clarke, February 1990
c    modified 1: September, 1990 by David Clarke; moved magnetic fields
c                to face-centres.
c    modified 2: minimal rewrite for ZEUS-MP by M-MML 10.3.98
c    modified 3: added restrictions for symmetry about J and K axes;
c                John Hayes, 10/2005
c
c  PURPOSE: This routine sets boundary values for the two terms in the
c  3-emf (centred on the 3-edges).  The active zones for "emf3" are:
c
c    i = is to ie+1;  j = js to je+1;  k = ks to ke
c
c  In order to update both the active and ghost zones of the 1- and
c  2-magnetic field components, all edges in the boundary regions are
c  required.  This gives a complete grid of values for "emf3".  Thus,
c  the ranges for the boundary values are:
c
c    k-boundaries:   i = is  , ie+1   j = js  , je+1
c    i-boundaries:                    j = js  , je+1   k = ks-2, ke+2
c    j-boundaries:   i = is-2, ie+3                    k = ks-2, ke+2
c
c  Note that the boundary values must be set even if ks > ksmn, etc.
c  because the emfs are stored in worker arrays and it is likely that
c  the boundary values have been overwritten.
c
c  See comments in BVALD.
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
      integer  :: i, j, k, ii, jj, kk, kone, jone
c
      real(rl) :: v1b2(in,jn,kn), v2b1(in,jn,kn)
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c------------------------  K - B O U N D A R Y  ------------------------
c-----------------------------------------------------------------------
c
      if(ldimen .eq. 1) then
       jone = 0
       kone = 0
      else
       jone = 1
       if(ldimen .eq. 3) then
        kone = 1
       else
        kone = 0
       endif
      endif
c
      if(ldimen .ne. 3) go to 111
c
c      Inner k boundary.
c
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1

       if (niks(1).eq.0 .or. niks(1).eq.4) then
         do kk = 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v1b2(1, 1, ks  +kk), 1, k_slice, n3m, 
     &                   9600+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v2b1(1, 1, ks  +kk), 1, k_slice, n3m, 
     &                   9700+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v1b2(1, 1, ks-2+kk), 1, k_slice, n3m, 
     &                   9400+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v2b1(1, 1, ks-2+kk), 1, k_slice, n3m, 
     &                   9500+nsub+kk, comm3d , req(nreq), ierr)
         enddo
       else
         do j=js,je+1
cdir$ ivdep
           do i=is,ie+1
             if ( abs(nikb12(i,j)) .eq. 1) then
               v1b2(i,j,ks-1) = v1b2(i,j,ks  )
               v1b2(i,j,ks-2) = v1b2(i,j,ks+1)
               v2b1(i,j,ks-1) = v2b1(i,j,ks  )
               v2b1(i,j,ks-2) = v2b1(i,j,ks+1)
             endif
             if (nikb12(i,j) .eq. 2) then
               v1b2(i,j,ks-1) = v1b2(i,j,ks  )
               v1b2(i,j,ks-2) = v1b2(i,j,ks  )
               v2b1(i,j,ks-1) = v2b1(i,j,ks  )
               v2b1(i,j,ks-2) = v2b1(i,j,ks  )
             endif
             if (nikb12(i,j) .eq. 3) then
               v1b2(i,j,ks-1) = emf3ikb(i,j,1)
               v1b2(i,j,ks-2) = emf3ikb(i,j,2)
               v2b1(i,j,ks-1) = 0.0
               v2b1(i,j,ks-2) = 0.0
             endif
             if (nikb12(i,j) .eq. 5) then
               v1b2(i,j,ks-1) =-v1b2(i,j,ks  )
               v1b2(i,j,ks-2) =-v1b2(i,j,ks+1)
               v2b1(i,j,ks-1) =-v2b1(i,j,ks  )
               v2b1(i,j,ks-2) =-v2b1(i,j,ks+1)
             endif
           enddo
         enddo
       endif
c 
c      Outer k boundary.
c 
       if (noks(1).eq.0 .or. noks(1).eq.4) then
         do kk = 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v1b2(1, 1, ke-1+kk), 1, k_slice, n3p, 
     &                   9400+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v2b1(1, 1, ke-1+kk), 1, k_slice, n3p, 
     &                   9500+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v1b2(1, 1, ke+1+kk), 1, k_slice, n3p, 
     &                   9600+nsub+kk, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v2b1(1, 1, ke+1+kk), 1, k_slice, n3p, 
     &                   9700+nsub+kk, comm3d , req(nreq), ierr)
         enddo
       else
         do j=js,je+1
cdir$ ivdep
           do i=is,ie+1
             if ( abs(nokb12(i,j)) .eq. 1) then
               v1b2(i,j,ke+1) = v1b2(i,j,ke  )
               v1b2(i,j,ke+2) = v1b2(i,j,ke-1)
               v2b1(i,j,ke+1) = v2b1(i,j,ke  )
               v2b1(i,j,ke+2) = v2b1(i,j,ke-1)
             endif
             if (nokb12(i,j) .eq. 2) then
               v1b2(i,j,ke+1) = v1b2(i,j,ke  )
               v1b2(i,j,ke+2) = v1b2(i,j,ke  )
               v2b1(i,j,ke+1) = v2b1(i,j,ke  )
               v2b1(i,j,ke+2) = v2b1(i,j,ke  )
             endif
             if (nokb12(i,j) .eq. 3) then
               v1b2(i,j,ke+1) = emf3okb(i,j,1)
               v1b2(i,j,ke+2) = emf3okb(i,j,2)
               v2b1(i,j,ke+1) = 0.0
               v2b1(i,j,ke+2) = 0.0
             endif
             if (nokb12(i,j) .eq. 5) then
               v1b2(i,j,ke+1) =-v1b2(i,j,ke  )
               v1b2(i,j,ke+2) =-v1b2(i,j,ke-1)
               v2b1(i,j,ke+1) =-v2b1(i,j,ke  )
               v2b1(i,j,ke+2) =-v2b1(i,j,ke-1)
             endif
c 
           enddo
         enddo
       endif
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif /* MPI */
#ifndef MPI_USED
         do j=js,je+1
cdir$ ivdep
           do i=is,ie+1
             if ( abs(nikb12(i,j)) .eq. 1) then
               v1b2(i,j,ks-1) = v1b2(i,j,ks  )
               v1b2(i,j,ks-2) = v1b2(i,j,ks+1)
               v2b1(i,j,ks-1) = v2b1(i,j,ks  )
               v2b1(i,j,ks-2) = v2b1(i,j,ks+1)
             endif
             if (nikb12(i,j) .eq. 2) then
               v1b2(i,j,ks-1) = v1b2(i,j,ks  )
               v1b2(i,j,ks-2) = v1b2(i,j,ks  )
               v2b1(i,j,ks-1) = v2b1(i,j,ks  )
               v2b1(i,j,ks-2) = v2b1(i,j,ks  )
             endif
             if (nikb12(i,j) .eq. 3) then
               v1b2(i,j,ks-1) = emf3ikb(i,j,1)
               v1b2(i,j,ks-2) = emf3ikb(i,j,2)
               v2b1(i,j,ks-1) = 0.0
               v2b1(i,j,ks-2) = 0.0
             endif
             if (nikb12(i,j) .eq. 4) then
               v1b2(i,j,ks-1) = v1b2(i,j,ke  )
               v1b2(i,j,ks-2) = v1b2(i,j,ke-1)
               v2b1(i,j,ks-1) = v2b1(i,j,ke  )
               v2b1(i,j,ks-2) = v2b1(i,j,ke-1)
             endif
             if (nikb12(i,j) .eq. 5) then
               v1b2(i,j,ks-1) =-v1b2(i,j,ks  )
               v1b2(i,j,ks-2) =-v1b2(i,j,ks+1)
               v2b1(i,j,ks-1) =-v2b1(i,j,ks  )
               v2b1(i,j,ks-2) =-v2b1(i,j,ks+1)
             endif
           enddo
         enddo
c 
c      Outer k boundary.
c 
         do j=js,je+1
cdir$ ivdep
           do i=is,ie+1
             if ( abs(nokb12(i,j)) .eq. 1) then
               v1b2(i,j,ke+1) = v1b2(i,j,ke  )
               v1b2(i,j,ke+2) = v1b2(i,j,ke-1)
               v2b1(i,j,ke+1) = v2b1(i,j,ke  )
               v2b1(i,j,ke+2) = v2b1(i,j,ke-1)
             endif
             if (nokb12(i,j) .eq. 2) then
               v1b2(i,j,ke+1) = v1b2(i,j,ke  )
               v1b2(i,j,ke+2) = v1b2(i,j,ke  )
               v2b1(i,j,ke+1) = v2b1(i,j,ke  )
               v2b1(i,j,ke+2) = v2b1(i,j,ke  )
             endif
             if (nokb12(i,j) .eq. 3) then
               v1b2(i,j,ke+1) = emf3okb(i,j,1)
               v1b2(i,j,ke+2) = emf3okb(i,j,2)
               v2b1(i,j,ke+1) = 0.0
               v2b1(i,j,ke+2) = 0.0
             endif
             if (nokb12(i,j) .eq. 4) then
               v1b2(i,j,ke+1) = v1b2(i,j,ks  )
               v1b2(i,j,ke+2) = v1b2(i,j,ks+1)
               v2b1(i,j,ke+1) = v2b1(i,j,ks  )
               v2b1(i,j,ke+2) = v2b1(i,j,ks+1)
             endif
             if (nokb12(i,j) .eq. 5) then
               v1b2(i,j,ke+1) =-v1b2(i,j,ke  )
               v1b2(i,j,ke+2) =-v1b2(i,j,ke-1)
               v2b1(i,j,ke+1) =-v2b1(i,j,ke  )
               v2b1(i,j,ke+2) =-v2b1(i,j,ke-1)
             endif
c 
           enddo
         enddo
#endif /* NO MPI */
c
c-----------------------------------------------------------------------
c------------------------  I - B O U N D A R Y  ------------------------
c-----------------------------------------------------------------------
c
111   CONTINUE
c
c      Inner i boundary.
c
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1
       if (niis(1).eq.0 .or. niis(1).eq.4) then
         do ii = 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v1b2(is+1+ii, 1, 1), 1, i_slice, n1m, 
     &                   10000+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v2b1(is+1+ii, 1, 1), 1, i_slice, n1m, 
     &                   10100+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v1b2(is-2+ii, 1, 1), 1, i_slice, n1m, 
     &                   9800+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v2b1(is-2+ii, 1, 1), 1, i_slice, n1m, 
     &                   9900+nsub+ii, comm3d , req(nreq), ierr)
         enddo
       else
         do k=ks-2*kone,ke+2*kone
cdir$ ivdep
           do j=js,je+jone
             if (niib2(j,k) .eq. 1) then
               v1b2(is  ,j,k) = 0.0
               v1b2(is-1,j,k) =-v1b2(is+1,j,k)
               v1b2(is-2,j,k) =-v1b2(is+2,j,k)
               v2b1(is  ,j,k) = 0.0
               v2b1(is-1,j,k) =-v2b1(is+1,j,k)
               v2b1(is-2,j,k) =-v2b1(is+2,j,k)
             endif
             if (niib2(j,k) .eq.-1) then
              if(lgeom .eq. 3) then
               v1b2(is  ,j,k) = 0.0
               v1b2(is-1,j,k) = v1b2(is+1,j,k)
               v1b2(is-2,j,k) = v1b2(is+2,j,k)
               v2b1(is  ,j,k) = 0.0
               v2b1(is-1,j,k) = v2b1(is+1,j,k)
               v2b1(is-2,j,k) = v2b1(is+2,j,k)
              else
               v1b2(is  ,j,k) = 0.0
               v1b2(is-1,j,k) =-v1b2(is+1,j,k)
               v1b2(is-2,j,k) =-v1b2(is+2,j,k)
               v2b1(is  ,j,k) = 0.0
               v2b1(is-1,j,k) =-v2b1(is+1,j,k)
               v2b1(is-2,j,k) =-v2b1(is+2,j,k)
              endif
             endif
             if (niib2(j,k) .eq. 2) then
               v1b2(is-1,j,k) = v1b2(is  ,j,k)
               v1b2(is-2,j,k) = v1b2(is  ,j,k)
               v2b1(is-1,j,k) = v2b1(is  ,j,k)
               v2b1(is-2,j,k) = v2b1(is  ,j,k)
             endif
             if (niib2(j,k) .eq. 3) then
               v1b2(is  ,j,k) = emf3iib(j,k,1)
               v1b2(is-1,j,k) = emf3iib(j,k,2)
               v1b2(is-2,j,k) = emf3iib(j,k,3)
               v2b1(is  ,j,k) = 0.0
               v2b1(is-1,j,k) = 0.0
               v2b1(is-2,j,k) = 0.0
             endif
            if (niib2(j,k) .eq. 5) then
               v1b2(is  ,j,k) = 0.0
               v1b2(is-1,j,k) = v1b2(is+1,j,k)
               v1b2(is-2,j,k) = v1b2(is+2,j,k)
               v2b1(is-1,j,k) = v2b1(is+1,j,k)
               v2b1(is-2,j,k) = v2b1(is+2,j,k)
             endif
           enddo
         enddo
       endif
c 
c      Outer i boundary.
c 
       if (nois(1).eq.0 .or. nois(1).eq.4) then
         do ii = 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v1b2(ie-1+ii, 1, 1), 1, i_slice, n1p, 
     &                   9800+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v2b1(ie-1+ii, 1, 1), 1, i_slice, n1p, 
     &                   9900+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v1b2(ie+2+ii, 1, 1), 1, i_slice, n1p, 
     &                   10000+nsub+ii, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v2b1(ie+2+ii, 1, 1), 1, i_slice, n1p, 
     &                   10100+nsub+ii, comm3d , req(nreq), ierr)
         enddo
       else
         do k=ks-2*kone,ke+2*kone
cdir$ ivdep
           do j=js,je+jone
             if ( abs(noib2(j,k)) .eq. 1) then
               v1b2(ie+1,j,k) = 0.0
               v1b2(ie+2,j,k) =-v1b2(ie  ,j,k)
               v1b2(ie+3,j,k) =-v1b2(ie-1,j,k)
               v2b1(ie+1,j,k) = 0.0
               v2b1(ie+2,j,k) =-v2b1(ie  ,j,k)
               v2b1(ie+3,j,k) =-v2b1(ie-1,j,k)
             endif
             if (noib2(j,k) .eq. 2) then
c*if alias PROBLEM.eq.advect
c               v1b2(ie+2,j,k) = 3.0 * ( v1b2(ie+1,j,k) - v1b2(ie  ,j,k) )
c       1                      + v1b2(ie-1,j,k)
c               v1b2(ie+3,j,k) = 3.0 * ( v1b2(ie+2,j,k) - v1b2(ie+1,j,k) )
c       1                      + v1b2(ie  ,j,k)
c*else
               v1b2(ie+2,j,k) = v1b2(ie+1,j,k)
               v1b2(ie+3,j,k) = v1b2(ie+1,j,k)
c*endif
               v2b1(ie+2,j,k) = v2b1(ie+1,j,k)
               v2b1(ie+3,j,k) = v2b1(ie+1,j,k)
             endif
             if (noib2(j,k) .eq. 3) then
               v1b2(ie+1,j,k) = emf3oib(j,k,1)
               v1b2(ie+2,j,k) = emf3oib(j,k,2)
               v1b2(ie+3,j,k) = emf3oib(j,k,3)
               v2b1(ie+1,j,k) = 0.0
               v2b1(ie+2,j,k) = 0.0
               v2b1(ie+3,j,k) = 0.0
             endif
             if (noib2(j,k) .eq. 5) then
               v1b2(ie+1,j,k) = 0.0
               v1b2(ie+2,j,k) = v1b2(ie  ,j,k)
               v1b2(ie+3,j,k) = v1b2(ie-1,j,k)
               v2b1(ie+2,j,k) = v2b1(ie  ,j,k)
               v2b1(ie+3,j,k) = v2b1(ie-1,j,k)
             endif
c c 
           enddo
         enddo
       endif
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif /* MPI */
#ifndef MPI_USED
         do k=ks-2*kone,ke+2*kone
cdir$ ivdep
           do j=js,je+jone
             if (niib2(j,k) .eq. 1) then
               v1b2(is  ,j,k) = 0.0
               v1b2(is-1,j,k) =-v1b2(is+1,j,k)
               v1b2(is-2,j,k) =-v1b2(is+2,j,k)
               v2b1(is  ,j,k) = 0.0
               v2b1(is-1,j,k) =-v2b1(is+1,j,k)
               v2b1(is-2,j,k) =-v2b1(is+2,j,k)
             endif
             if (niib2(j,k) .eq.-1) then
              if(lgeom .eq. 3) then
               v1b2(is  ,j,k) = 0.0
               v1b2(is-1,j,k) = v1b2(is+1,j,k)
               v1b2(is-2,j,k) = v1b2(is+2,j,k)
               v2b1(is  ,j,k) = 0.0
               v2b1(is-1,j,k) = v2b1(is+1,j,k)
               v2b1(is-2,j,k) = v2b1(is+2,j,k)
              else
               v1b2(is  ,j,k) = 0.0
               v1b2(is-1,j,k) =-v1b2(is+1,j,k)
               v1b2(is-2,j,k) =-v1b2(is+2,j,k)
               v2b1(is  ,j,k) = 0.0
               v2b1(is-1,j,k) =-v2b1(is+1,j,k)
               v2b1(is-2,j,k) =-v2b1(is+2,j,k)
              endif
             endif
             if (niib2(j,k) .eq. 2) then
               v1b2(is-1,j,k) = v1b2(is  ,j,k)
               v1b2(is-2,j,k) = v1b2(is  ,j,k)
               v2b1(is-1,j,k) = v2b1(is  ,j,k)
               v2b1(is-2,j,k) = v2b1(is  ,j,k)
             endif
             if (niib2(j,k) .eq. 3) then
               v1b2(is  ,j,k) = emf3iib(j,k,1)
               v1b2(is-1,j,k) = emf3iib(j,k,2)
               v1b2(is-2,j,k) = emf3iib(j,k,3)
               v2b1(is  ,j,k) = 0.0
               v2b1(is-1,j,k) = 0.0
               v2b1(is-2,j,k) = 0.0
             endif
             if (niib2(j,k) .eq. 4) then
               v1b2(is-1,j,k) = v1b2(ie  ,j,k)
               v1b2(is-2,j,k) = v1b2(ie-1,j,k)
               v2b1(is-1,j,k) = v2b1(ie  ,j,k)
               v2b1(is-2,j,k) = v2b1(ie-1,j,k)
             endif
            if (niib2(j,k) .eq. 5) then
               v1b2(is  ,j,k) = 0.0
               v1b2(is-1,j,k) = v1b2(is+1,j,k)
               v1b2(is-2,j,k) = v1b2(is+2,j,k)
               v2b1(is-1,j,k) = v2b1(is+1,j,k)
               v2b1(is-2,j,k) = v2b1(is+2,j,k)
             endif
           enddo
         enddo
c 
c      Outer i boundary.
c 
         do k=ks-2*kone,ke+2*kone
cdir$ ivdep
           do j=js,je+jone
             if ( abs(noib2(j,k)) .eq. 1) then
               v1b2(ie+1,j,k) = 0.0
               v1b2(ie+2,j,k) =-v1b2(ie  ,j,k)
               v1b2(ie+3,j,k) =-v1b2(ie-1,j,k)
               v2b1(ie+1,j,k) = 0.0
               v2b1(ie+2,j,k) =-v2b1(ie  ,j,k)
               v2b1(ie+3,j,k) =-v2b1(ie-1,j,k)
             endif
             if (noib2(j,k) .eq. 2) then
c*if alias PROBLEM.eq.advect
c               v1b2(ie+2,j,k) = 3.0 * ( v1b2(ie+1,j,k) - v1b2(ie  ,j,k) )
c       1                      + v1b2(ie-1,j,k)
c               v1b2(ie+3,j,k) = 3.0 * ( v1b2(ie+2,j,k) - v1b2(ie+1,j,k) )
c       1                      + v1b2(ie  ,j,k)
c*else
               v1b2(ie+2,j,k) = v1b2(ie+1,j,k)
               v1b2(ie+3,j,k) = v1b2(ie+1,j,k)
c*endif
               v2b1(ie+2,j,k) = v2b1(ie+1,j,k)
               v2b1(ie+3,j,k) = v2b1(ie+1,j,k)
             endif
             if (noib2(j,k) .eq. 3) then
               v1b2(ie+1,j,k) = emf3oib(j,k,1)
               v1b2(ie+2,j,k) = emf3oib(j,k,2)
               v1b2(ie+3,j,k) = emf3oib(j,k,3)
               v2b1(ie+1,j,k) = 0.0
               v2b1(ie+2,j,k) = 0.0
               v2b1(ie+3,j,k) = 0.0
             endif
             if (noib2(j,k) .eq. 4) then
               v1b2(ie+2,j,k) = v1b2(is+1,j,k)
               v1b2(ie+3,j,k) = v1b2(is+2,j,k)
               v2b1(ie+2,j,k) = v2b1(is+1,j,k)
               v2b1(ie+3,j,k) = v2b1(is+2,j,k)
             endif
             if (noib2(j,k) .eq. 5) then
               v1b2(ie+1,j,k) = 0.0
               v1b2(ie+2,j,k) = v1b2(ie  ,j,k)
               v1b2(ie+3,j,k) = v1b2(ie-1,j,k)
               v2b1(ie+2,j,k) = v2b1(ie  ,j,k)
               v2b1(ie+3,j,k) = v2b1(ie-1,j,k)
             endif
c c 
           enddo
         enddo
#endif /* NO MPI */
c-----------------------------------------------------------------------
c------------------------  J - B O U N D A R Y  ------------------------
c-----------------------------------------------------------------------
c
      if(ldimen .eq. 1) return
c
c      Inner j boundary.
c
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1
c
       if (nijs(1).eq.0 .or. nijs(1).eq.4) then
         do jj = 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v1b2(1, js+1+jj, 1), 1, j_slice, n2m, 
     &                   10400+nsub+jj, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v2b1(1, js+1+jj, 1), 1, j_slice, n2m, 
     &                   10500+nsub+jj, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v1b2(1, js-2+jj, 1), 1, j_slice, n2m, 
     &                   10200+nsub+jj, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v2b1(1, js-2+jj, 1), 1, j_slice, n2m, 
     &                   10300+nsub+jj, comm3d , req(nreq), ierr)
         enddo
       else
         do k=ks-2,ke+2
cdir$ ivdep
           do i=is-2,ie+3
             if ( abs(nijb1(i,k)) .eq. 1) then
               v1b2(i,js  ,k) = 0.0
               v1b2(i,js-1,k) =-v1b2(i,js+1,k)
               v1b2(i,js-2,k) =-v1b2(i,js+2,k)
               v2b1(i,js  ,k) = 0.0
               v2b1(i,js-1,k) =-v2b1(i,js+1,k)
               v2b1(i,js-2,k) =-v2b1(i,js+2,k)
             endif
             if (nijb1(i,k) .eq. 2) then
               v1b2(i,js-1,k) = v1b2(i,js  ,k)
               v1b2(i,js-2,k) = v1b2(i,js  ,k)
               v2b1(i,js-1,k) = v2b1(i,js  ,k)
               v2b1(i,js-2,k) = v2b1(i,js  ,k)
             endif
             if (nijb1(i,k) .eq. 3) then
               v1b2(i,js  ,k) = 0.0
               v1b2(i,js-1,k) = 0.0
               v1b2(i,js-2,k) = 0.0
               v2b1(i,js  ,k) =-emf3ijb(i,k,1)
               v2b1(i,js-1,k) =-emf3ijb(i,k,2)
               v2b1(i,js-2,k) =-emf3ijb(i,k,3)
             endif
             if (nijb1(i,k) .eq. 5) then
               v1b2(i,js-1,k) = v1b2(i,js+1,k)
               v1b2(i,js-2,k) = v1b2(i,js+2,k)
               v2b1(i,js  ,k) = 0.0
               v2b1(i,js-1,k) = v2b1(i,js+1,k)
               v2b1(i,js-2,k) = v2b1(i,js+2,k)
             endif
           enddo
         enddo
       endif
c 
c      Outer j boundary.
c 
c
       if (nojs(1).eq.0 .or. nojs(1).eq.4) then
         do jj = 0, 1
           nreq = nreq + 1
           call MPI_ISEND(v1b2(1, je-1+jj, 1), 1, j_slice, n2p, 
     &                   10200+nsub+jj, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_ISEND(v2b1(1, je-1+jj, 1), 1, j_slice, n2p, 
     &                   10300+nsub+jj, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v1b2(1, je+2+jj, 1), 1, j_slice, n2p, 
     &                   10400+nsub+jj, comm3d , req(nreq), ierr)
           nreq = nreq + 1
           call MPI_IRECV(v2b1(1, je+2+jj, 1), 1, j_slice, n2p, 
     &                   10500+nsub+jj, comm3d , req(nreq), ierr)
         enddo
       else
         do k=ks-2,ke+2
cdir$ ivdep
           do i=is-2,ie+3
             if ( abs(nojb1(i,k)) .eq. 1) then
               v1b2(i,je+1,k) = 0.0
               v1b2(i,je+2,k) =-v1b2(i,je  ,k)
               v1b2(i,je+3,k) =-v1b2(i,je-1,k)
               v2b1(i,je+1,k) = 0.0
               v2b1(i,je+2,k) =-v2b1(i,je  ,k)
               v2b1(i,je+3,k) =-v2b1(i,je-1,k)
             endif
             if (nojb1(i,k) .eq. 2) then
               v1b2(i,je+2,k) = v1b2(i,je+1,k)
               v1b2(i,je+3,k) = v1b2(i,je+1,k)
c*if alias PROBLEM.eq.advect
c               v2b1(i,je+2,k) = 2.0 * v2b1(i,je+1,k) - v2b1(i,je  ,k)
c               v2b1(i,je+3,k) = 2.0 * v2b1(i,je+2,k) - v2b1(i,je+1,k)
c*else
               v2b1(i,je+2,k) = v2b1(i,je+1,k)
               v2b1(i,je+3,k) = v2b1(i,je+1,k)
c*endif
             endif
             if (nojb1(i,k) .eq. 3) then
               v1b2(i,je+1,k) = 0.0
               v1b2(i,je+2,k) = 0.0
               v1b2(i,je+3,k) = 0.0
               v2b1(i,je+1,k) =-emf3ojb(i,k,1)
               v2b1(i,je+2,k) =-emf3ojb(i,k,2)
               v2b1(i,je+3,k) =-emf3ojb(i,k,3)
             endif
             if (nojb1(i,k) .eq. 5) then
               v1b2(i,je+2,k) = v1b2(i,je  ,k)
               v1b2(i,je+3,k) = v1b2(i,je-1,k)
               v2b1(i,je+1,k) = 0.0
               v2b1(i,je+2,k) = v2b1(i,je  ,k)
               v2b1(i,je+3,k) = v2b1(i,je-1,k)
             endif
           enddo
         enddo
       endif
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif /* MPI */
#ifndef MPI_USED
         do k=ks-2,ke+2
cdir$ ivdep
           do i=is-2,ie+3
             if ( abs(nijb1(i,k)) .eq. 1) then
               v1b2(i,js  ,k) = 0.0
               v1b2(i,js-1,k) =-v1b2(i,js+1,k)
               v1b2(i,js-2,k) =-v1b2(i,js+2,k)
               v2b1(i,js  ,k) = 0.0
               v2b1(i,js-1,k) =-v2b1(i,js+1,k)
               v2b1(i,js-2,k) =-v2b1(i,js+2,k)
             endif
             if (nijb1(i,k) .eq. 2) then
               v1b2(i,js-1,k) = v1b2(i,js  ,k)
               v1b2(i,js-2,k) = v1b2(i,js  ,k)
               v2b1(i,js-1,k) = v2b1(i,js  ,k)
               v2b1(i,js-2,k) = v2b1(i,js  ,k)
             endif
             if (nijb1(i,k) .eq. 3) then
               v1b2(i,js  ,k) = 0.0
               v1b2(i,js-1,k) = 0.0
               v1b2(i,js-2,k) = 0.0
               v2b1(i,js  ,k) =-emf3ijb(i,k,1)
               v2b1(i,js-1,k) =-emf3ijb(i,k,2)
               v2b1(i,js-2,k) =-emf3ijb(i,k,3)
             endif
             if (nijb1(i,k) .eq. 4) then
               v1b2(i,js-1,k) = v1b2(i,je  ,k)
               v1b2(i,js-2,k) = v1b2(i,je-1,k)
               v2b1(i,js-1,k) = v2b1(i,je  ,k)
               v2b1(i,js-2,k) = v2b1(i,je-1,k)
             endif
             if (nijb1(i,k) .eq. 5) then
               v1b2(i,js-1,k) = v1b2(i,js+1,k)
               v1b2(i,js-2,k) = v1b2(i,js+2,k)
               v2b1(i,js  ,k) = 0.0
               v2b1(i,js-1,k) = v2b1(i,js+1,k)
               v2b1(i,js-2,k) = v2b1(i,js+2,k)
             endif
           enddo
         enddo
c 
c      Outer j boundary.
c 
         do k=ks-2,ke+2
cdir$ ivdep
           do i=is-2,ie+3
             if ( abs(nojb1(i,k)) .eq. 1) then
               v1b2(i,je+1,k) = 0.0
               v1b2(i,je+2,k) =-v1b2(i,je  ,k)
               v1b2(i,je+3,k) =-v1b2(i,je-1,k)
               v2b1(i,je+1,k) = 0.0
               v2b1(i,je+2,k) =-v2b1(i,je  ,k)
               v2b1(i,je+3,k) =-v2b1(i,je-1,k)
             endif
             if (nojb1(i,k) .eq. 2) then
               v1b2(i,je+2,k) = v1b2(i,je+1,k)
               v1b2(i,je+3,k) = v1b2(i,je+1,k)
c*if alias PROBLEM.eq.advect
c               v2b1(i,je+2,k) = 2.0 * v2b1(i,je+1,k) - v2b1(i,je  ,k)
c               v2b1(i,je+3,k) = 2.0 * v2b1(i,je+2,k) - v2b1(i,je+1,k)
c*else
               v2b1(i,je+2,k) = v2b1(i,je+1,k)
               v2b1(i,je+3,k) = v2b1(i,je+1,k)
c*endif
             endif
             if (nojb1(i,k) .eq. 3) then
               v1b2(i,je+1,k) = 0.0
               v1b2(i,je+2,k) = 0.0
               v1b2(i,je+3,k) = 0.0
               v2b1(i,je+1,k) =-emf3ojb(i,k,1)
               v2b1(i,je+2,k) =-emf3ojb(i,k,2)
               v2b1(i,je+3,k) =-emf3ojb(i,k,3)
             endif
             if (nojb1(i,k) .eq. 4) then
               v1b2(i,je+2,k) = v1b2(i,js+1,k)
               v1b2(i,je+3,k) = v1b2(i,js+2,k)
               v2b1(i,je+2,k) = v2b1(i,js+1,k)
               v2b1(i,je+3,k) = v2b1(i,js+2,k)
             endif
             if (nojb1(i,k) .eq. 5) then
               v1b2(i,je+2,k) = v1b2(i,je  ,k)
               v1b2(i,je+3,k) = v1b2(i,je-1,k)
               v2b1(i,je+1,k) = 0.0
               v2b1(i,je+2,k) = v2b1(i,je  ,k)
               v2b1(i,je+3,k) = v2b1(i,je-1,k)
             endif
           enddo
         enddo
#endif /* NO MPI */
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////              B V A L E M F 3              \\\\\\\\\\
c
c=======================================================================
c
