c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                D I V E R G                \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine diverg ( c1, c2, c3, inorm, isum, div, sumd )
c
c    dac:zeus3d.diverg <------------ computes divergence of vector field
c    from jms:zeus2d.divb                                september, 1990
c
c    written by: David Clarke
c    modified 1: 11-17-99 by PSLi, for ZeusMP.
c
c  PURPOSE:  Computes divergence of vector field (c1, c2, c3), where
c  each of c1, c2, and c3 are face-centred quantities.  Therefore, the
c  divergence will be a zone-centred quantity, and is computed in the
c  range: ismn to iemx, jsmn to jemx, and ksmn to kemx.
c
c  INPUT VARIABLES:
c    c1       1-component of vector field
c    c2       2-component of vector field
c    c3       3-component of vector field
c    inorm    =1 => normalise divergence
c             =0 => do not normalise divergence
c    isum     =1 => perform reduction sum on "div"
c             =0 => no reduction
c
c  OUTPUT VARIABLES:
c    div      divergence of vector field
c    sumd     sum of "div"
c
c  LOCAL VARIABLES:
c
c  EXTERNALS: [NONE]
c
c-----------------------------------------------------------------------
c
c
      use real_prec
      use config
      use param
      use grid
c
      implicit NONE
c
      integer  :: i, j, k, l, ip1, jp1, kp1, inorm, isum
      integer  :: kone, km1  !asif
c
      real(rl) :: sumd, sumn, nmax, nmaxi
c
      real(rl) :: norm(ijkn), sumdk(ijkn), sumnk(ijkn), nmaxk(ijkn)
c
      real(rl) :: c1 (in,jn,kn), c2(in,jn,kn), c3(in,jn,kn),
     .            div(in,jn,kn)
c
c      Careful!  "wa1d" - "wi1d" are used by BININT.
c
C       equivalence   ( norm    , ww1d     ), ( sumdk   , wx1d     )
C     1             , ( sumnk   , wy1d     ), ( nmaxk   , wz1d     )
c
c-----------------------------------------------------------------------
c	asif
        if (ldimen .eq. 3 )then
        kone=1
        else
        kone=0
        endif
c
c      Compute divergence of vector field ("div").
c
       do 30 k=ks-2*kone,ke+2*kone
         kp1 = k + kone
         do 20 j=js-2,je+2
           jp1 = j + 1
           do 10 i=is-2,ie+2
             ip1        = i + 1
             div(i,j,k) = ( g2a(ip1) * g31a(ip1) * c1(ip1,j,k)
     1                    - g2a(i  ) * g31a(i  ) * c1(i  ,j,k) )
     2                  *                         dvl1ai(i)
     3                  + ( g32a(jp1) * c2(i,jp1,k)
     4                    - g32a(j  ) * c2(i,j  ,k) )
     5                  *   g2bi(i)             * dvl2ai(j)
     6                  + ( c3(i,j,kp1) - c3(i,j,k) )
     7                  *   g31bi(i) * g32bi(j) * dvl3ai(k)
10         continue
20       continue
30     continue
c
c      Perform reduction on "div" if desired.  The reduction is done
c  first over each k-sweep, then for all k.  This is done to aid the
c  EDITOR autotasking process.
c
       sumd = 0.0
       if (isum .eq. 1) then
         do 60 k=ks,ke
           do 50 j=js,je
             do 40 i=is,ie
               sumd = sumd + div(i,j,k)
40           continue
50         continue
60       continue
       endif
c
c      Normalise divergence field if desired.
c
       if (inorm .eq. 1) then
c
c      Evaluate two normalising constants:
c
c  1.  sumn = sum over all grid zones the ratio: (average absolute
c      magnetic field) / (sum of grid zone dimensions), and
c
c  2.  nmax = maximum over the grid of the above ratios.
c
         sumn = 0.0
         nmax = 0.0
         do 90 k=ks,ke
           kp1      = k + kone
           do 80 j=js,je
             jp1 = j + 1
             do 70 i=is,ie
               ip1      = i + 1
               norm (1) = 0.5 * ( abs ( c1(ip1,j  ,k  ) + c1(i,j,k) )
     1                          + abs ( c2(i  ,jp1,k  ) + c2(i,j,k) )
     2                          + abs ( c3(i  ,j  ,kp1) + c3(i,j,k) ) )
     3                  / (dx1a(i) + g2b(i) * dx2a(j)   +
     4                     g31b(i) * g32b(j) * dx3a(k) )
               sumn     = sumn     + norm(1)
               nmax     =   max ( nmax    , norm(1) )
70           continue
80         continue
90       continue
c
c      Apply normalising constant "sumn" to the scalar "sumd" and
c  "nmax" to the array "div".
c
         if (sumn .eq. 0.0) sumd  = 0.0
         if (sumn .ne. 0.0) sumd  = sumd / sumn
         if (nmax .eq. 0.0) nmaxi = 0.0
         if (nmax .ne. 0.0) nmaxi = 1.0 / nmax
         do 120 k=ks,ke
           do 110 j=js,je
             do 100 i=is,ie
               div(i,j,k) = div(i,j,k) * nmaxi
100          continue
110        continue
120      continue
c
       endif
c
       return
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                D I V E R G                \\\\\\\\\\
c
c=======================================================================
c
