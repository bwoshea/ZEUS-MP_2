c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                  M N M X                  \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine mnmx ( qty, i1, j1, k1, i2, j2, k2
     1                 , qmin, imin, jmin, kmin
     2                 , qmax, imax, jmax, kmax )
c
c    dac:zeus3d.mnmx <------------- finds extrema of a regular 3-d array
c    from mln:zeus04.minmax; jms:zeus2d.minmax            february, 1990
c
c    written by: David Clarke
c    modified 1: converted to ZEUS-MP by Mordecai-Mark Mac Low M
c
c  PURPOSE: This subroutine returns the maximum and minimum value of a
c  rectangular 3-D array, along with the coordinates of the extrema.
c
c  INPUT VARIABLES:
c    qty             the 3-D array to be searched for extrema.
c    i1              inner i index.
c    j1              inner j index.
c    k1              inner k index.
c    i2              outer i index.
c    j2              outer j index.
c    k2              outer k index.
c
c  OUTPUT VARIABLES:
c    qmin            minimum value
c    imin,jmin,kmin  coordinates of minimum value
c    qmax            maximum value
c    imax,jmax,kmax  coordinates of maximum value
c
c  LOCAL VARIABLES:
c
c  j-sweep
c    qmaxj (qminj)   vector of maximum (minimum) values of "qty" from
c                    each i-sweep.  This vector is filled during a
c                    j-sweep.
c    imaxj           i-index of each "qmaxj"
c    iminj           i-index of each "qminj"
c
c  k-sweep
c    qmaxk (qmink)   vector of maximum (minimum) values of "qmaxj"
c                    ("qminj") from each j-sweep.
c    imaxk           i-index of each "qmaxk"
c    jmaxk           j-index of each "qmaxk"
c    imink           i-index of each "qmink"
c    jmink           j-index of each "qmink"
c
c  grand maximum
c    qmax (qmin)     maximum (minimum) value of "qmaxk" ("qmink").
c    imax            i-index of "qmax"
c    jmax            j-index of "qmax"
c    kmax            k-index of "qmax"
c    imin            i-index of "qmin"
c    jmin            j-index of "qmin"
c    kmin            k-index of "qmin"
c
c----------------------------------------------------------------------
c
      use real_prec
      use config
      use param
      use scratch
      use grid
c
      integer  :: j, k, i1, j1, k1, i2, j2, k2, imin, jmin,
     .            kmin, imax, jmax, kmax
      real(rl) :: qmin, qmax
c
      integer  :: imaxj(ijkn), iminj(ijkn),
     .            imaxk(ijkn), imink(ijkn),
     .            jmaxk(ijkn), jmink(ijkn)
c
      real(rl) :: qmaxj(ijkn), qminj(ijkn),
     1            qmaxk(ijkn), qmink(ijkn)
c
      real(rl) :: qty(in,jn,kn)
c
      integer  :: ijkl, ijkx, ijks, i
c
      real(rl) :: q(in*jn*kn)
c
c WARNING: Scratch array wg3d is used by this routine (mnmx)!
c
c
c      External statements
c
       integer       ISMIN   , ISMAX
       external      ISMIN   , ISMAX
c
c-----------------------------------------------------------------------
c
       k = k1
       j = j1
       ijkl = 0
       do 20 k=k1,k2
         do 10 j=j1,j2
           do 5 i=i1,i2
             ijkl = ijkl + 1
             q(ijkl) = qty(i,j,k)
5          continue
10       continue
20     continue

       ijkl = (i2 - i1 + 1) * (j2 - j1 + 1) * (k2 - k1 + 1)
       ijkx = ISMAX(ijkl,q,1)
       qmax = q(ijkx)
       kmax = ijkx/((i2-i1+1)*(j2-j1+1)) + k1
       jmax = (ijkx-(kmax-k1)*(i2-i1+1)*(j2-j1+1))/(i2-i1+1) + j1
       imax = ijkx-(kmax-k1)*(i2-i1+1)*(j2-j1+1)
     1      - (jmax-j1)*(i2-i1+1) + i1
       ijks = ISMIN(ijkl,q,1)
       qmin = q(ijks)
       kmin = ijks/((i2-i1+1)*(j2-j1+1)) + k1
       jmin = (ijks-(kmin-k1)*(i2-i1+1)*(j2-j1+1))/(i2-i1+1) + j1
       imin = ijks-(kmin-k1)*(i2-i1+1)*(j2-j1+1)
     1      - (jmin-j1)*(i2-i1+1) + i1
c
       return
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                  M N M X                  \\\\\\\\\\
c
c=======================================================================
