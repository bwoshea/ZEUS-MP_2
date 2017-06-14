c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 N E W D T                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine newdt (ibeg,iend,jbeg,jend,kbeg,kend
     &                ,  imin,jmin,kmin
     &                ,  dtcsm,dtv1m,dtv2m,dtv3m,dtalm,dttoi2m)
c
c    mln:zeus3d.nudt <------------------------ mhd time step calculation
c                                                              may, 1986
c
c    written by: Mike Norman
c    modified 1: May, 1986 by David Clarke; adapted for mhd.
c    modified 2: October, 1987 by Mike Norman; reworked for covariant
c                formulation.
c    modified 3: June, 1988 by Jim Stone; incorporated into ZEUS2D.
c    modified 4: February, 1990 by David Clarke; incorporated into
c                ZEUS3D.
c    modified 5: September, 1990 by David Clarke; moved magnetic fields
c                to face-centres.
c    modified 6: Feb. 15, 1996 by Robert Fiedler; completely rewritten
c                for ZEUS-MP.
c    modified 7: Dec. 23, 1996 by Robert Fiedler; added radiation.
c
c  PURPOSE:  Computes the new timestep for explicit calculations from
c  the values of the field variables updated by the source and transport
c  steps.
c
c  In explicit calculations, the timestep is given by:
c
c     dt = courno * sqrt [ min ( dtcs**2 + dtv1**2 + dtv2**2 + dtv3**2
c                              + dtal**2 + dtqq**2 ) ]
c
c  where the variable names are described below.  The timestep can be
c  reduced in size by any amount, but can be larger than the old timstep
c  by no more than a factor of 1.26.
c
c BOUNDARY VALUES USED:
c
c    Macro defined  var   ii    oi    ij    oj    ik    ok
c    -------------  ---  ----  ----  ----  ----  ----  ----
c                    d   is-1        js-1        ks-1
c    TOTAL_ENERGY    d   is-1  ie+1  js-1  je+1  ks-1  ke+1
c    TOTAL_ENERGY    s1        ie+1
c    TOTAL_ENERGY    s2                    je+1
c    TOTAL_ENERGY    s3                                ke+1
c
c    Note that s1,s2,s3 are stored in w3da,w3db,w3dc.
c
c  LOCAL VARIABLES:
c
c  i-sweep
c  dr*i      inverse of distance across zone in 1-, 2-, or 3-direction
c  drimax    maximum of dr1i, dr2i, and dr3i
c  dt**i2i   square of the inverse time step of the ** physical process
c            gathered during the i-sweep.  Possible values of ** are:
c                cs = sound speed
c                v1 = fluid motion in x1 direction
c                v2 = fluid motion in x2 direction
c                v3 = fluid motion in x2 direction
c                al = Alfven speed
c                qq = artificial viscosity
c            The first five are vectors, while the last is a scalar
c            which has been computed in ARTIFICIALVISC (passed in root).
c  dttoi2i   sum of dt**i2i (without artificial viscosity contribution)
c
c  j-sweep
c  dttoi2j   vector of maximum values of dttoi2i from each i-sweep.
c            This vector is filled during a j-sweep.
c  imaxj     i-index of each dttoi2j
c  dt**i2j   values of above at zone i=imaxj
c
c  k-sweep
c  dttoi2k   vector of maximum values of dttoi2j from each j-sweep.
c  imaxk     i-index of each dttoi2k.
c  jmaxk     j-index of each dttoi2k.
c  dt**i2j   values of above at zone i=imaxk, j=jmaxk
c
c  grand maximum inverse time
c  imax      i-index where dttoi2k is a maximum
c  jmax      j-index where dttoi2k is a maximum
c  kmax      k-index where dttoi2k is a maximum
c  dt**      time step of the ** physical process at i=imax, j=jmax,
c            and k=kmax
c  dtnew     new timestep, which is limited to be no greater than 1.26
c            times the previous timestep.
c
c  EXTERNALS:
c
c-----------------------------------------------------------------------
c
      use real_prec
      use config
      use param
      use field
      use root
      use grid
      use scratch
c
       implicit NONE
c
      integer  :: ibeg,iend,jbeg,jend,kbeg,kend
      integer  :: i   , j   , k   , n, jone, kone,
     .            imin, jmin, kmin
c
      real(rl) :: se
      real(rl) :: dr1is, dr2is, dr3is, drimaxs
      real(rl) :: dtcsm, dtv1m, dtv2m, dtv3m, dtalm,
     .            dttoi2, dttoi2m
c
      real(rl) :: tota, totai
      integer  :: nbeg
c-----------------------------------------------------------------------
c
c Find the minimum time step required by the Courant condition for
c this tile.  We will first compute 1/dt**2 required by each of the
c various physical processes in the calculation, and save the maximum
c value of their sum at each zone.  In the process, we will compute
c the updated velocity.
c
       if(ldimen .eq. 1) then
        jone = 0
       else
        jone = 1
       endif
       if(ldimen .eq. 3) then
        kone = 1
       else
        kone = 0
       endif
c
       do 30 k=kbeg,kend
        do 20 j=jbeg,jend
         do 10 i=ibeg,iend
          dr1is      = dx1ai(i)
          dr2is      = g2bi (i) * dx2ai(j)
          dr3is      = g31bi(i) * g32bi(j) * dx3ai(k)
          drimaxs    =   max ( dr1is  , dr2is  , dr3is   )
          v1(i,j,k)  = w3da(i,j,k) * 2.0 /(d(i-1,j     ,k)+d(i,j,k))
          if(xvgrid) then
           v2(i,j,k)  = w3db(i,j,k)*2.0/(d(i  ,j-jone,k  )+d(i,j,k))
     1                * g2 bni(i)
           v3(i,j,k)  = w3dc(i,j,k)*2.0/(d(i  ,j  ,k-kone)+d(i,j,k))
     1                * g31bni(i) * g32bi(j)
           else
            v2(i,j,k)  = w3db(i,j,k)*2.0/(d(i  ,j-jone,k  )+d(i,j,k))
     1                 * g2 b i(i)
            v3(i,j,k)  = w3dc(i,j,k)*2.0/(d(i  ,j  ,k-kone)+d(i,j,k))
     1                 * g31b i(i) * g32bi(j)
           endif
           if(xiso) then
            dtcs       = ( ciso * drimaxs   )**2 
           else ! xiso
            if(xtotnrg) then
             se         = e(i,j,k) / d(i,j,k)
     1                  - ( ( v1(i,j,k) + w3da(i+1,j  ,k  ) * 2.0
     &                        / ( d(i  ,j  ,k  ) + d(i+1,j  ,k  ) ) )**2
     2                    + ( v2(i,j,k) + w3db(i  ,j+1,k  ) * 2.0
     &                        / ( d(i  ,j  ,k  ) + d(i  ,j+1,k  ) ) )**2
     3                    + ( v3(i,j,k) + w3dc(i  ,j  ,k+1) * 2.0
     &                        / ( d(i  ,j  ,k  ) + d(i  ,j  ,k+1) ) )**2
     &                    ) * 0.125
            else ! xtotnrg
             se         = e(i,j,k) / d(i,j,k)
            endif ! xtotnrg
            dtcs       = gamma * gamm1 * se        * drimaxs**2
           endif ! xiso
           dtv1       = ( (v1(i,j,k) - vg1(i)) * dr1is   )**2
           dtv2       = ( (v2(i,j,k) - vg2(j)) * dr2is   )**2
           if(ldimen .eq. 3) then
            dtv3       = ( (v3(i,j,k) - vg3(k)) * dr3is   )**2
           else
            dtv3       = 0.0
           endif
c
c Copy updated er/d from w3dh to er.
c
           if(lrad .ne. 0) then
            er(i,j,k)  = w3dh(i,j,k) * d(i,j,k)
           endif
c
c     convert partial densities back into mass fractions
c
           if(nspec .gt. 1) then
            do n = 1, nspec
             abun(i,j,k,n) = w4da(i,j,k,n)
            enddo
           endif
c
           if(xmhd) then
            dtal       = ( ( b1(i,j,k) + b1(i+1,j     ,k  ) )**2
     1                   + ( b2(i,j,k) + b2(i  ,j+jone,k  ) )**2
     2                   + ( b3(i,j,k) + b3(i  ,j     ,k+kone) )**2 )
     3                 * 0.25 * drimaxs**2 / d(i,j,k)
           else ! xmhd
            dtal       = 0.0
           endif ! xmhd
           dttoi2     = dtcs          + dtv1
     1                + dtv2          + dtv3
     2                + dtal
           if (dttoi2 .gt. dttoi2m) then
             dttoi2m = dttoi2
             dtcsm   = dtcs
             dtv1m   = dtv1
             dtv2m   = dtv2
             dtv3m   = dtv3
             dtalm   = dtal
             imin    = i
             jmin    = j
             kmin    = k
           endif
10       continue
20      continue
30     continue
c
       return
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                 N E W D T                 \\\\\\\\\\
c
c=======================================================================
c
c
