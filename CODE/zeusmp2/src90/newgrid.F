c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////               N E W G R I D               \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
      subroutine newgrid
c
c  PURPOSE: Controls grid motion by computing grid velocities and "new"
c  grid variables (variables at advanced time, to be used in TRANSPRT).
c
c  EXTERNALS: NEWVG
c             NEWX1
c             NEWX2
c             NEWX3
c
c  modified for F90 version of ZEUS-MP by J. Hayes, 5-2003
c-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use grid
      use bndry
      use field
      use root
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif /* MPI_USED */
      use mpipar
c
      implicit NONE
c
      integer  :: i, j, k
c
c=======================================================================
c  return if there is no grid motion in any direction
c
      if ((x1fac .eq. 0.0) .and. (x2fac .eq. 0.0) .and.
     .    (x3fac .eq. 0.0)) return
c
c  update "X1" grid
c
      if (x1fac .ne. 0.0) then
       do i = 1, 6
        bvstat(1,i) = 0
        bvstat(3,i) = 0
       enddo
       nreq = 0
       nsub = nsub + 1
       call bval d(3,3,0,0,0,0, d)
       call bvalv1(3,3,0,0,0,0,v1)
c
        call scopy (in,  x1a n ,1,  x1a ,1)
        call scopy (in,  x1b n ,1,  x1b ,1)
        call scopy (in, dx1a n ,1, dx1a ,1)
        call scopy (in, dx1b n ,1, dx1b ,1)
        call scopy (in, g2 a n ,1, g2 a ,1)
        call scopy (in, g2 b n ,1, g2 b ,1)
        call scopy (in, g31a n ,1, g31a ,1)
        call scopy (in, g31b n ,1, g31b ,1)
        call scopy (in,dvl1a n ,1,dvl1a ,1)
        call scopy (in,dvl1b n ,1,dvl1b ,1)
c
        do 10 i=is-2,ie+2
         x1ai  (i) = 1.0/(x1a  (i)+tiny)
         x1bi  (i) = 1.0/(x1b  (i)+tiny)
         dx1ai (i) = 1.0/(dx1a (i)+tiny)
         dx1bi (i) = 1.0/(dx1b (i)+tiny)
         g2ai  (i) = 1.0/(g2a  (i)+tiny)
         g2bi  (i) = 1.0/(g2b  (i)+tiny)
         g31ai (i) = 1.0/(g31a (i)+tiny)
         g31bi (i) = 1.0/(g31b (i)+tiny)
         dvl1ai(i) = 1.0/(dvl1a(i)+tiny)
         dvl1bi(i) = 1.0/(dvl1b(i)+tiny)
10      continue
c
c  update grid velocities
c
#ifdef MPI_USED
       if(nreq .ne. 0) then
        call mpi_waitall(nreq, req, stat, ierr)
        nreq = 0
       endif
#endif /* MPI_USED */
       call newvg
c
       call newx1
c
       do i = is-2, ie+2
        dx1ani(i) = 1.0D0/dx1an(i)
        dx1bni(i) = 1.0D0/dx1bn(i)
       enddo
c
       do k = 1, kn
        do j = 1, jn
         v3oib(j,k,1) = 0.3*x1b(ie+1)*sin(x2b(j))
         v3oib(j,k,2) = 0.3*x1b(ie+2)*sin(x2b(j))
        enddo
       enddo
      endif ! x1fac
c
c  update "X2" grid 
c
      if (x2fac .ne. 0.0) then
        call scopy (jn,  x2a n ,1,  x2a ,1)
        call scopy (jn,  x2b n ,1,  x2b ,1)
        call scopy (jn, dx2a n ,1, dx2a ,1)
        call scopy (jn, dx2b n ,1, dx2b ,1)
        call scopy (jn, g32a n ,1, g32a ,1)
        call scopy (jn, g32b n ,1, g32b ,1)
        call scopy (jn, g4 a n ,1, g4 a ,1)
        call scopy (jn, g4 b n ,1, g4 b ,1)
        call scopy (jn,dvl2a n ,1,dvl2a ,1)
        call scopy (jn,dvl2b n ,1,dvl2b ,1)
        call scopy (jn,dvl2a ni,1,dvl2ai,1)
        call scopy (jn,dvl2b ni,1,dvl2bi,1)
        do 20 j=js-2,je+2
         dx2ai(j) = 1.0/dx2a(j)
         dx2bi(j) = 1.0/dx2b(j)
         if(lgeom .eq. 3) then
          dg32ad2(j) = cos(x2a(j))
          dg32bd2(j) = cos(x2b(j))
         endif
20      continue
        call newx2
      endif ! x2fac
c
      if(x3fac .ne. 0)then
        call scopy (kn,  x3a n ,1,  x3a ,1)
        call scopy (kn,  x3b n ,1,  x3b ,1)
        call scopy (kn, dx3a n ,1, dx3a ,1)
        call scopy (kn, dx3b n ,1, dx3b ,1)
        call scopy (kn,dvl3a n ,1,dvl3a ,1)
        call scopy (kn,dvl3b n ,1,dvl3b ,1)
        call scopy (kn,dvl3a ni,1,dvl3ai,1)
        call scopy (kn,dvl3b ni,1,dvl3bi,1)
        do 30 k=ks-2,ke+2
         dx3ai(k) = 1.0/dx3a(k)
         dx3bi(k) = 1.0/dx3b(k)
30      continue
        call newx3
      endif ! x3fac
c
      return
      end
c
      subroutine  scopy(n,sx,incx,sy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c
      use real_prec
      use param
c
      implicit NONE
c
      real(rl) :: sx(*),sy(*)
      integer  :: i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return
      end
