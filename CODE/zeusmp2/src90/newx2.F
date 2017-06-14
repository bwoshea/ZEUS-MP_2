c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 N E W X 2                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
      subroutine newx2
c
c  PURPOSE: Computes "new" x2 grid variables (grid variables at advanced
c  time) to be used in TRANSPRT.  Grid values are calculated for 
c  j=js-2 to je+2, except for dvl2a (j=js,je+2) and dvl2b (j=js+1,je+2).
c  Note similarity of the expressions used to those in the grid
c  generator routine GRIDI.
c
c  EXTERNALS: [none]
c
c  LOCALS:
c-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use grid
      use root
c
      implicit NONE
c
      integer  :: j
      real(rl) :: qa,qb,qc,qd
c
c=======================================================================
      x2an(js-2) = x1an(js-2) + vg2(js-2)*dt
      do 10 j=js-1,je+2
         x2an(j  ) = x2a (j) + vg2(j)*dt
        dx2an(j-1) = x2an(j) - x2an(j-1)
10    continue
      dx2an(je+2) = (dx2an(je+1)/dx2an(je))*dx2an(je+1)
c
      dx2bn(js-2) = dx2an(js-2)
       x2bn(js-2) =  x2an(js-1) - 0.5*dx2an(js-2)
      do 20 j=js-1,je+2
         x2bn(j) = x2an(j) + 0.5*dx2an(j)
        dx2bn(j) = x2bn(j) - x2bn(j-1)
20    continue
c
      do 30 j=js-2,je+2
       if(lgeom .eq. 1) then
        g32a h (j) = 1.0
        g32b h (j) = 1.0
        g4 a h (j) = 0.5*(x2a(j) + x2an(j))
        g4 b h (j) = 0.5*(x2b(j) + x2bn(j))
        g32a n (j) = 1.0
        g32b n (j) = 1.0
        g4 a n (j) = x2an(j)
        g4 b n (j) = x2bn(j)
       endif ! XYZ
       if(lgeom .eq. 2) then
        g32a h (j) =   0.5*(x2a(j) + x2an(j))
        g32b h (j) =   0.5*(x2b(j) + x2bn(j))
        g4 a h (j) = 0.125*(x2a(j) + x2an(j))**2
        g4 b h (j) = 0.125*(x2b(j) + x2bn(j))**2
        g32a n (j) =     x2an(j)
        g32b n (j) =     x2bn(j)
        g4 a n (j) = 0.5*x2an(j)**2
        g4 b n (j) = 0.5*x2bn(j)**2
       endif ! ZRP
       if(lgeom .eq. 3) then
        g32a h (j) =  sin(0.5*(x2a(j) + x2an(j)))
        g32b h (j) =  sin(0.5*(x2b(j) + x2bn(j)))
        g4 a h (j) = -cos(0.5*(x2a(j) + x2an(j)))
        g4 b h (j) = -cos(0.5*(x2b(j) + x2bn(j)))
        g32a n (j) =  sin(x2a n (j))
        g32b n (j) =  sin(x2b n (j))
        g4 a n (j) = -cos(x2a n (j))
        g4 b n (j) = -cos(x2b n (j))
       endif ! RTP
30    continue
c
c  New volume factors
c
      do 40 j=js-2,je+1 
        dvl2ah(j  ) = g4ah(j+1) - g4ah(j)
        dvl2an(j  ) = g4an(j+1) - g4an(j)
40    continue
      do 50 j=js-2,je+1
        dvl2bh(j+1) = g4bh(j+1) - g4bh(j)
        dvl2bn(j+1) = g4bn(j+1) - g4bn(j)
50    continue
      do j = js-2, je+1
       dvl2ahi(j) = 1.0D0/dvl2ah(j)
       dvl2ani(j) = 1.0D0/dvl2an(j)
       dvl2bhi(j) = 1.0D0/dvl2bh(j)
       dvl2bni(j) = 1.0D0/dvl2bn(j)
      enddo
c
      return
      end
