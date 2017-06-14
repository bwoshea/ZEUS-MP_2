c======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////               O P A C _ D                 \\\\\\\\\\
c======================================================================
c
c
      subroutine opac_d
c
c     Driver routine for computing opacities when radiation is included
c
c     Written by: John Hayes
c     modified: 4-13-99; added entry for "kem" to OPACITY calls
c     modified: 5-26-2003; rewritten for F90 code
c
c
      use real_prec
      use config
      use param
      use root
      use field
      use grid
      use bndry
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use radiation
      use opac
      use cons
c
      implicit NONE
c
      integer  :: i,j,k
c
      real(rl) :: temp
c
c----------------------------------------------------------------------
c
c     update ghost zones along 1-coordinate
c
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1
#endif /* MPI_USED */
       call bvald  (3,3,0,0,0,0,d)
       call bvale  (3,3,0,0,0,0,e)
c
c     compute opacities in all interior zones
c
       call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt,
     .              is, ie, js, je, ks, ke)
#ifdef MPI_USED
c
c     wait for communications to complete
c
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif /* MPI_USED */
      if(ldimen .eq. 1) then
c
c     compute opacities in 1-coordinate ghost zones and return
c
       call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt,
     .              is-2, is-1, js, je, ks, ke)
       call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt,
     .              ie+1, ie+2, js, je, ks, ke)
c
       go to 999
      endif ! ldimen = 1
      if(ldimen .gt. 1) then
c
c     update ghost zones along 2-coordinate
c
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1
#endif /* MPI_USED */
       call bvald  (0,0,3,3,0,0,d)
       call bvale  (0,0,3,3,0,0,e)
c
c     compute opacities in 1-coordinate ghost zones
c
       call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt,
     .              is-2, is-1, js, je, ks, ke)
       call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt,
     .              ie+1, ie+2, js, je, ks, ke)
#ifdef MPI_USED
c
c     wait for communications to complete
c
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif /* MPI_USED */
       if(ldimen .eq. 2) then
c
c     compute opacities in 2-coordinate ghost zones and return
c
        call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt,
     .               is, ie, js-2, js-1, ks, ke)
        call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt,
     .               is, ie, je+1, je+2, ks, ke)
c
        go to 999
       else ! ldimen > 2
#ifdef MPI_USED
        nreq = 0
        nsub = nsub + 1
#endif /* MPI_USED */
        call bvald  (0,0,0,0,3,3,d)
        call bvale  (0,0,0,0,3,3,e)
        call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt,
     .               is, ie, js-2, js-1, ks, ke)
        call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt,
     .               is, ie, je+1, je+2, ks, ke)
#ifdef MPI_USED
        if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif /* MPI_USED */
        call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt,
     .               is, ie, js, je, ks-2, ks-1)
        call opacity(e, d, gamma, kapr, kap, sig, dkapdt, kem, dkemdt,
     .               is, ie, js, je, ke+1, ke+2)
c
        go to 999
       endif ! ldimen > 2
      endif ! ldimen > 1
c
999   return
      end
