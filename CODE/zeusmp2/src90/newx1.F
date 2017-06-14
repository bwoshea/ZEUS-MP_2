c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 N E W X 1                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
      subroutine newx1
c
c  PURPOSE: Computes "new" x1 grid variables (grid variables at advanced
c  timestep) to be used in TRANSPRT.  Grid values are calculated for
c  i=is-2 to ie+2, except for dvl1a (i=is,ie+2) and dvl1b (i=is+1,ie+2).
c
c     implemented in ZEUS-MP by John Hayes
c
c-----------------------------------------------------------------------
c
      use real_prec
      use config
      use param
      use grid
      use root
      use scratch
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
      integer  :: i, ibeg
      real(rl) :: vol1an(in) , vol1bn(in), qa, qb, qc, qd, vfac
      real(rl) :: ibuf(2)
c
c=======================================================================
c
      if(lgeom .eq. 1) vfac = 1.0D0
      if(lgeom .eq. 2) vfac = 1.0D0
      if(lgeom .eq. 3) vfac = 1.0D0/3.0D0
#ifdef MPI_USED
c
c-----------------------------------------------------------------------
c     post receive requests for updated grid values near tile boundaries
c-----------------------------------------------------------------------
c
      nreq = 0
      if(nois(1) .eq. 0) then
       nreq = nreq + 1
       call mpi_irecv(ibuf(1), 1 , MPI_FLOAT,
     .                n1p, 666, comm3d, req(nreq), ierr)
      endif
#endif /* MPI_USED */
c
c-----------------------------------------------------------------------
c     update x1an straddling inner boundary...
c-----------------------------------------------------------------------
c
      x1an(is-2) = x1a(is-2) + vg1(is-2)*dt
      do 10 i=is-1,is+2
         x1an(i  ) = x1a (i) + vg1(i)*dt
        dx1an(i-1) = x1an(i) - x1an(i-1)
10    continue
#ifdef MPI_USED
c
c-----------------------------------------------------------------------
c     ...and if this boundary is a tile boundary, send the updated
c     dx1an to the neighboring processor.
c-----------------------------------------------------------------------
c
      if(niis(1) .eq. 0) then
       nreq = nreq + 1
       call mpi_isend(dx1an(is+1), 1, MPI_FLOAT,
     .                n1m, 666, comm3d, req(nreq), ierr)
      endif
#endif /* MPI_USED */
c
c-----------------------------------------------------------------------
c     finish updating x1an, dx1an away from inner tile boundary...
c-----------------------------------------------------------------------
c
      do 11 i=is+3,ie+2
         x1an(i  ) = x1a (i) + vg1(i)*dt
        dx1an(i-1) = x1an(i) - x1an(i-1)
11    continue
#ifdef MPI_USED
c
c-----------------------------------------------------------------------
c     ...and give all tiles a chance to collect the boundary values
c     and finish up.
c-----------------------------------------------------------------------
c
      if(nreq .ne. 0) then
       call mpi_waitall(nreq, req, stat, ierr)
       nreq = 0
      endif
#endif /* MPI_USED */
      if(nois(1) .ne. 0) then
       dx1an(ie+2) = (dx1an(ie+1)/dx1an(ie)) * dx1an(ie+1)
      else
       dx1an(ie+2) = ibuf(1)
      endif
c
c-----------------------------------------------------------------------
c     repeat process for "b" grid
c-----------------------------------------------------------------------
c
      if(niis(1) .eq. 0) then
       ibeg = is-2 
      else
       ibeg = is-1
      endif
#ifdef MPI_USED
      if(niis(1) .eq. 0) then
       nreq = nreq + 1
       call mpi_irecv(ibuf(1), 2, MPI_FLOAT,
     .                n1m, 667, comm3d, req(nreq), ierr)
      endif
#endif /* MPI_USED */
      do 20 i=ie-2, ie+2
         x1bn(i) = x1an(i) + 0.5*dx1an(i)
20    continue
#ifdef MPI_USED
      if(nois(1) .eq. 0) then
       nreq = nreq + 1
       call mpi_isend(x1bn(ie-2), 2, MPI_FLOAT,
     .                n1p, 667, comm3d, req(nreq), ierr)
      endif
#endif /* MPI_USED */
      do 22 i=ibeg, ie-3
         x1bn(i) = x1an(i) + 0.5*dx1an(i)
22    continue
#ifdef MPI_USED
      if(nreq .ne. 0) then
       call mpi_waitall(nreq, req, stat, ierr)
       nreq = 0
      endif
#endif /* MPI_USED */
c
      if(niis(1) .eq. 0) then
       x1bn(is-2) = ibuf(2)
       do i=is-1, ie+2
         dx1bn(i) = x1bn(i) - x1bn(i-1)
       enddo
       dx1bn(is-2) = x1bn(is-2) - ibuf(1)
      else
       dx1bn(is-2) = dx1an(is-2)
        x1bn(is-2) =  x1an(is-1) - 0.5*dx1an(is-2)
       do i=is-1, ie+2
         dx1bn(i) = x1bn(i) - x1bn(i-1)
       enddo
      endif
c
c-----------------------------------------------------------------------
c     compute inverse coordinate values at advanced time...
c-----------------------------------------------------------------------
c
      do i = is-2,ie+2
       x1ani(i) = 1.0/(x1an(i)+tiny)
       x1bni(i) = 1.0/(x1bn(i)+tiny)
      enddo
c
c-----------------------------------------------------------------------
c     ...and update the time-centered metric arrays AND their inverses
c-----------------------------------------------------------------------
c
      do 30 i=is-2,ie+2
       if(lgeom .eq. 1) then
        g2 a h (i) = 1.0
        g2 b h (i) = 1.0
        g31a h (i) = 1.0
        g31b h (i) = 1.0
        g2 a n (i) = 1.0
        g2 b n (i) = 1.0
        g31a n (i) = 1.0
        g31b n (i) = 1.0
       endif ! XYZ
       if(lgeom .eq. 2) then
        g2 a h (i) = 1.0
        g2 b h (i) = 1.0
        g31a h (i) = 1.0
        g31b h (i) = 1.0
        g2 a n (i) = 1.0
        g2 b n (i) = 1.0
        g31a n (i) = 1.0
        g31b n (i) = 1.0
       endif ! ZRP
       if(lgeom .eq. 3) then
        g2 a h (i) = 0.5*(x1a(i) + x1an(i))
        g2 b h (i) = 0.5*(x1b(i) + x1bn(i))
        g31a h (i) = 0.5*(x1a(i) + x1an(i))
        g31b h (i) = 0.5*(x1b(i) + x1bn(i))
        g2 a n (i) = x1a n (i)
        g2 b n (i) = x1b n (i)
        g31a n (i) = x1a n (i)
        g31b n (i) = x1b n (i)
       endif ! RTP
30    continue
      do i = is-2, ie+2
       g2 ani(i) = 1.0/(g2 an(i)+tiny)
       g2 bni(i) = 1.0/(g2 bn(i)+tiny)
       g31ani(i) = 1.0/(g31an(i)+tiny)
       g31bni(i) = 1.0/(g31bn(i)+tiny)
      enddo
c
c-----------------------------------------------------------------------
c     New volume factors
c-----------------------------------------------------------------------
c
      vol1an(is-2) = g2an(is-2)*g31an(is-2)*x1an(is-2)
      do 40 i=is-2,ie+1
        vol1an(i+1) = g2an(i+1)*g31an(i+1)*x1an(i+1)
        dvl1an(i  ) = vfac*(vol1an(i+1) - vol1an(i))
40    continue
      do i = is-2, ie+1
       dvl1ani(i) = 1.0D0/dvl1an(i)
      enddo
c
      vol1bn(is-2) = g2bn(is-2)*g31bn(is-2)*x1bn(is-2)
      do 50 i=is-2,ie+1
        vol1bn(i+1) = g2bn(i+1)*g31bn(i+1)*x1bn(i+1)
        dvl1bn(i+1) = vfac*(vol1bn(i+1) - vol1bn(i))
50    continue
c
      do i = is-1, ie+2
       dvl1bni(i) = 1.0D0/dvl1bn(i)
      enddo
c
      return
      end
