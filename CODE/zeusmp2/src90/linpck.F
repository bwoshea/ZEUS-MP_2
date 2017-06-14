c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
      integer function isamax(n,sx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      implicit NONE
      real*8 sx(*),smax
      integer i,incx,ix,n
c
      isamax = 0
      if( n .lt. 1 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = abs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = abs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = abs(sx(1))
      do 30 i = 2,n
         if(abs(sx(i)).le.smax) go to 30
         isamax = i
         smax = abs(sx(i))
   30 continue
      return
      end
c
      integer function ismax(n,sx,incx)
c
c     finds the index of element having max. value.
c     jack dongarra, linpack, 3/11/78. (added to this file by j. stone)
c
      implicit NONE
      real*8 sx(*),smax
      integer i,incx,ix,n
c
      ismax = 0
      if( n .lt. 1 ) return
      ismax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = sx(1)
      ix = ix + incx
      do 10 i = 2,n
         if(sx(ix).le.smax) go to 5
         ismax = i
         smax = sx(ix)
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = sx(1)
      do 30 i = 2,n
         if(sx(i).le.smax) go to 30
         ismax = i
         smax = sx(i)
   30 continue
      return
      end
      integer function ismin(n,sx,incx)
c
c     finds the index of element having min. value.
c     jack dongarra, linpack, 3/11/78.  (added to this file by j. stone)
c
      implicit NONE
      real*8 sx(*),smin
      integer i,incx,ix,n
c
      ismin = 0
      if( n .lt. 1 ) return
      ismin = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smin = sx(1)
      ix = ix + incx
      do 10 i = 2,n
         if(sx(ix).ge.smin) go to 5
         ismin = i
         smin = sx(ix)
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smin = sx(1)
      do 30 i = 2,n
         if(sx(i).ge.smin) go to 30
         ismin = i
         smin = sx(i)
   30 continue
      return
      end
      real*8 function sasum(n,sx,incx)
c
c     takes the sum of the absolute values.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit NONE
      real*8 sx(*),stemp
      integer i,incx,m,mp1,n,nincx
c
      sasum = 0.0e0
      stemp = 0.0e0
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        stemp = stemp + abs(sx(i))
   10 continue
      sasum = stemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + abs(sx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        stemp = stemp + abs(sx(i)) + abs(sx(i + 1)) + abs(sx(i + 2))
     *  + abs(sx(i + 3)) + abs(sx(i + 4)) + abs(sx(i + 5))
   50 continue
   60 sasum = stemp
      return
      end
