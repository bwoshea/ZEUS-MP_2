c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 H S M O C                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                  University of California at San Diego
c
c=======================================================================
c
       subroutine hsmoc_2D( emf1, emf2, emf3 )
c
c    dac:zeus3d.mocemfs <-------------------------- MoC estimate of emfs
c                                                          october, 1992
c
c    written by: David Clarke
c    modified 1: Byung-Il Jun, July 1994
c                implemented John Hawley and Jim Stone's scheme to
c                fix pt. explosion of magnetic field in passive field.
c                Basically, this scheme mixes emfs computed with simple
c                upwinding(Evans and Hawley) and MoC.
c                The upwinded values are used to compute the wave
c                speeds for the characteristic cones for the MoC part.
c    modified 2: Robert Fiedler, February 1995
c                upgraded to ZEUS-3D version 3.4 -- improved cache
c                utilization and added parallelization directives for 
c                SGI multiprocessors.
c    modified 3: Mordecai-Mark Mac Low, December 1997 - March 1998
c                rewritten for ZEUS-MP without overlapping.  Calls to 
c                interpolation routines have been inlined.
c    modified 4: PSLi, December 1999
c                minor modications to prevent scratch arrays overwritten.
c
c  PURPOSE:  Uses the Method of Characteristics (MoC, invented by Jim
c  Stone, John Hawley, Chuck Evans, and Michael Norman; see Stone and
c  Norman, ApJS, v80, p791) to evaluate the velocity and magnetic field
c  needed to estimate emfs that are properly upwinded in the character-
c  istic velocities for the set of equations describing transverse
c  Alfven waves.  This is *not* the full MHD characteristic problem, but
c  a subset which has been found (reference above) to yield good results
c  for the propagation of Alfven waves.  This routine differs from the
c  previous routines MOC1, MOC2, and MOC3 in version 3.1 in that the
c  Lorentz forces are computed *before* the emfs are estimated.  Thus,
c  the emfs now use the velocities which have been updated with all the
c  source terms, including the transverse Lorenz forces.
c
c  The characteristic equations governing the propagation of Alfven
c  waves in the 1-direction are (see ZEUS3D notes "Method of Character-
c  istics"):
c
c  "plus" characteristic (C+):
c
c    ( db1/dt + (v2 - a2) * db1/dx2 + (v3 - a3) * db1/dx3 ) / sqrt(d)
c  + ( dv1/dt + (v2 - a2) * dv1/dx2 + (v3 - a3) * dv1/dx3 )  =  S    (1)
c
c  "minus" characteristic (C-):
c
c    ( db1/dt + (v2 + a2) * db1/dx2 + (v3 + a3) * db1/dx3 ) / sqrt(d)
c  - ( dv1/dt + (v2 + a2) * dv1/dx2 + (v3 + a3) * dv1/dx3 )  = -S    (2)
c
c  where   a2 = b2/sqrt(d) is the Alfven velocity in the 2-direction
c          a3 = b3/sqrt(d) is the Alfven velocity in the 3-direction
c          g1, g2, g3 are the metric factors
c          S = b1 * ( b2/g2 * dg1/dx2 + b3/g3 * dg1/dx3 )
c
c  Equations (1) and (2) can be written in Lagrangian form:
c
c      1    D+/-           D+/-
c   ------- ----(b1)  +/-  ----(v1)  =  +/- S                        (3)
c   sqrt(d)  Dt             Dt
c
c  where the Lagrangian derivatives are given by
c
c  D+/-     d                   d                   d
c  ----  =  --  +  (v2 -/+ a2) ---  +  (v3 -/+ a3) ---               (4)
c   Dt      dt                 dx2                 dx3
c
c  Differencing equations (3) [e.g. D+(b1) = b* - b+; D-(b1) = b* - b-],
c  and then solving for the advanced time values of b* and v*, one gets:
c                             _                                _
c           sqrt (d+ * d-)   |     b+          b-               |
c  b* =  ------------------- |  -------- + --------- + v+ - v-  |    (5)
c        sqrt(d+) + sqrt(d-) |_ sqrt(d+)    sqrt(d-)           _|
c
c                 1
c  v* =  ------------------- [ v+*sqrt(d+) + v-*sqrt(d-) + b+ - b- ]
c        sqrt(d+) + sqrt(d-)                                         (6)
c
c     + S Dt
c
c  where b+(-), and v+(-) are the upwinded values of the magnetic field
c  and velocity interpolated to the time-averaged bases of C+ (C-), and
c  d+(-) are estimates of the density along each characteristic path
c  during the time interval Dt.
c
c  Equations (1) and (2) would suggest that when estimating "emf2" for
c  example, that the interpolated values for "v1" and "b1" be upwinded
c  in both the 2- and 3- components of the characteristic velocity.  It
c  turns out that this is impractical numerically, and so only the
c  "partial" characteristics are tracked.  While estimating "emf2", "v1"
c  and "b1" are upwinded only in the 3-component of the characteristic
c  velocities.  Conversely, while estimating "emf3", "v1" and "b1" are
c  upwinded only in the 2-component of the characteristic velocities.
c  Since updating "b1" requires both "emf2" and "emf3", the evolution of
c  "b1" will ultimately depend upon the full characteristics.  This
c  amounts to nothing more than directionally splitting the full MoC
c  algorithm.  The effects of such a directionally split implementation
c  are not fully known.  What is known is:
c
c  1) A non-directionally split implementation of the MoC algorithm is
c     not possible without either relocating the emfs to the zone
c     corners or the magnetic field components to the zone centres.
c     The former has been tried (change deck mocemf) and was found to
c     generate intolerable diffusion of magnetic field.  In addition,
c     the algorithm is not unconditionally stable.  The latter has not
c     been tried, but is dismissed on the grounds that div(B) will be
c     determined by truncation errors rather than machine round-off
c     errors.
c
c  2) A directionally split algorithm that is not also operator split
c     (ie, performing the Lorentz updates of the velocities separately
c     from the MoC estimation of the emfs) does not allow stable Alfven
c     wave propagation in 2-D.  Operator splitting the MoC algorithm so
c     that the transverse Lorentz forces are computed *before* the
c     magnetic field update does allow Alfven waves to propagate stably
c     in multi-dimensions but appears to introduce more diffusion for
c     sub-Alfvenic flow.  On the other hand, super-Alfvenic flow does
c     *not* appear to be more diffusive in the operator split scheme.
c
c  INPUT VARIABLES:
c
c  OUTPUT VARIABLES:
c    emf1      emf along 1-edge computed using MoC estimates of v2, b2,
c              v3, and b3.
c    emf2      emf along 2-edge computed using MoC estimates of v3, b3,
c              v1, and b1.
c    emf3      emf along 3-edge computed using MoC estimates of v1, b1,
c              v2, and b2.
c
c  LOCAL VARIABLES:
c
c    1-D variables
c    bave      spatially averaged magnetic field at edge
c    srdp      sqrt(density) along the plus  characteristic (C+)
c    srdm      sqrt(density) along the minus characteristic (C-)
c    vchp      characteristic velocity along C+ (v - va)
c    vchm      characteristic velocity along C- (v + va)
c    vpch      velocity interpolated to the time-centred footpoint of C+
c    vmch      velocity interpolated to the time-centred footpoint of C-
c    bpch      B-field  interpolated to the time-centred footpoint of C+
c    bmch      B-field  interpolated to the time-centred footpoint of C-
c    vsnm1     MoC estimate of v[n-1] used to evaluate emf[n], n=1,2,3
c    bsnm1     MoC estimate of b[n-1] used to evaluate emf[n], n=1,2,3
c
c    2-D variables
c    v3intj   upwinded v3 in 2-direction
c    b3intj   upwinded b3 in 2-direction
c    etc..
c
c    3-D variables
c    srd[n]    sqrt of spatially averaged density at [n]-face n=1,2,3
c    vsnp1     MoC estimate of v[n+1] used to evaluate emf[n], n=1,2,3
c    bsnp1     MoC estimate of b[n+1] used to evaluate emf[n], n=1,2,3
c    vsnp1     is reused as the vsnp1*bsnm1 term in emf[n]
c    bsnp1     is reused as the vsnm1*bsnp1 term in emf[n]
c
c  EXTERNALS:
c    BVALEMF1, BVALEMF2, BVALEMF3
c
c I have inlined these M-MML 8.3.98
c    X1ZC1D  , X2ZC1D  , X3ZC1D
c    X1INT1D , X2INT1D , X3INT1D
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
      use lor_scr
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer  :: i, j, k    
      integer  :: kone, km1,kp1  !asif
c
      real(rl) :: absb, sgnp, sgnm, q1, q2, src,
     .            qv1, qv2, qb1, qb2, q3 ,dqm ,dqp ,xi ,fact,
     .            dv(ijkn), db(ijkn)
c
      real(rl) :: bave (ijkn), srdp (ijkn), srdm (ijkn),
     .            srdpi(ijkn), srdmi(ijkn), vchp (ijkn),
     .            vchm (ijkn), vtmp (ijkn), btmp (ijkn),
     .            vpch (ijkn), vmch (ijkn), bpch (ijkn),
     .            bmch (ijkn), vsnm1(ijkn), bsnm1(ijkn),
     .            vave (ijkn), aave (ijkn)
c
      real(rl) :: vfl  (ijkn), vt   (ijkn), bt   (ijkn),
     .            vint (ijkn), bint (ijkn)
c
      real(rl) :: v3intj(jn,kn), b3intj(jn,kn),
     .            v2intk(jn,kn), b2intk(jn,kn),
     .            v1intk(kn,in), b1intk(kn,in),
     .            v3inti(kn,in), b3inti(kn,in),
     .            v2inti(in,jn), b2inti(in,jn),
     .            v1intj(in,jn), b1intj(in,jn)
c
      real(rl) :: emf1 (in,jn,kn), emf2 (in,jn,kn), emf3(in,jn,kn)
      real(rl) :: vsnp1(in,jn,kn), bsnp1(in,jn,kn)
c
C       equivalence   
C     .               ( bave    , w1da     ), 
C     .               ( srdp    , w1db     ),
C     1               ( srdm    , w1dc     ), 
C     .               ( srdpi   , w1dd     ),
C     1               ( srdmi   , w1de     ), 
C     .               ( vchp    , w1df     ),
C     1               ( vchm    , w1dg     ), 
C     .               ( vpch    , w1dh     ),
C     1               ( vmch    , w1di     ), 
C     .               ( bpch    , w1dj     ),
C     1               ( bmch    , w1dk     ),
C     1               ( vtmp    , w1dl     ),
C     1               ( btmp    , w1dm     ),
C     1               ( vave    , w1dn     ),
C     1               ( aave    , w1do     ),
C     1               ( vsnm1   , w1dp     ),
C     1               ( bsnm1   , w1dq     ),
C     1               ( vt      , w1dq     ),
C     1               ( bt      , w1dr     ),
C     1               ( vint    , w1ds     ),
C     1               ( bint    , w1dt     ),
C     1               ( vfl     , w1du     )
c
c there are no 2D scratch arrays like in ZEUS-3D, but these can all still 
c be equivalenced to each other. M-MML 4 Mar 98
C       equivalence   ( v3intj, v2intk, v1intk, v3inti, v2inti
C     2               , v1intj                          )
C     2             , ( b3intj, b2intk, b1intk, b3inti, b2inti
C     2               , b1intj                          )
c
c      Careful!  "wa3d" through "wc3d" are equivalenced in CT.
c                wa3d -> emf1    wb3d -> emf2    wc3d -> emf3
c
c      The worker arrays "we3d" and "wf3d" should still contain "srd2"
c  and "srd3" from LORENTZ. The worker array  "wd3d" will contain
c  "srd1", but "wd3d" is needed for "bsnp1" and "term2".  Thus "srd1"
c  will be recomputed once "srd3" is no longer needed.
c
cRAF We have more array space with wg3d, so don't recompute srd1.
cRAF SGIMP does not like equivalenced variables in the same loop nest,
cRAF so eliminate term1 and term2.
c
CPS
C       equivalence   
C                  ( srd1              , w3di     )
C     1             , ( srd2              , w3dj     )
C     1             , ( srd3              , w3df     ),
C     1               ( bsnp1             , w3dg     )
C
c
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c---- 1.  emf1 ---------------------------------------------------------
c-----------------------------------------------------------------------
c
c    BEGINNING OF I LOOPS
c
c
c      By following the characteristics of the flow in the 3-direction,
c  determine values for "v2" and "b2" ("vsnp1" and "bsnp1") to be used
c  in evaluating "emf1".
c
c     Compute upwinded b3 and v3 in the 2-direction and use them to
c     compute the wave speeds for the characteristic cones for the MoC.
c
CPS   Initialize vsnp1
       do k=1,kn
          do j=1,jn
             do i=1,in
                vsnp1(i,j,k)=0.
             enddo
          enddo
       enddo
C
       kone = 0
       km1  = ks
       kp1  = ks
       k    = ks
         do 9 i=is,ie
           do 7 j=js-2,je+2
             vfl   (j) = 0.5 * (v2(i,j,k) + v2(i,j,km1)) - vg2(j)
             vt    (j) = v3(i,j,k) - vg3(k)
             bt    (j) = b3(i,j,k)
7          continue
c
c           call x2int1d ( bt, vfl, iordb3, istpb3, k, i
c     1                  , g2b, g2bi, bint )
c           call x2int1d ( vt, vfl, iords3, istps3, k, i
c     1                  , g2b, g2bi, vint )

c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
           do 807 j=js-1,je+1
             dqm      = ( vt(j  ) - vt(j-1) ) * dx2bi(j  )
             dqp      = ( vt(j+1) - vt(j  ) ) * dx2bi(j+1)
             dv (j  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )

             dqm     = ( bt(j  ) - bt(j-1) ) * dx2bi(j  )
             dqp     = ( bt(j+1) - bt(j  ) ) * dx2bi(j+1)
             db (j  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )
 807       continue
c
c  2.  Perform an upwinded interpolation of "q" to the time-centred
c      bases of the characteristics.
c
           fact = dt * g2bi(i)
           do 808 j=js,je+1
             qv1    = vt(j-1) + dx2a(j-1) * dv (j-1)
             qv2    = vt(j  ) - dx2a(j  ) * dv (j  )
             qb1    = bt(j-1) + dx2a(j-1) * db (j-1)
             qb2    = bt(j  ) - dx2a(j  ) * db (j  )
c
             xi     = vfl(j) * fact
             q3     = sign ( haf, xi )
             vint(j)= ( 0.5 + q3 ) * ( qv1 - xi * dv (j-1) )
     1              + ( 0.5 - q3 ) * ( qv2 - xi * dv (j  ) )
             bint(j)= ( 0.5 + q3 ) * ( qb1 - xi * db (j-1) )
     1              + ( 0.5 - q3 ) * ( qb2 - xi * db (j  ) )
 808       continue

           do 8 j=js,je+1
             emf1(i,j,k) = vint(j)
             emf2(i,j,k) = bint(j)
8          continue
9      continue
c
          do 6 j=js,je+1
          do 6 i=is,ie
           vsnp1(i,j,k) = (v2(i,j,k) - vg2(j)) * emf2(i,j,k)
           bsnp1(i,j,k) = b2(i,j,k) * emf1(i,j,k)
6        continue
c
c-----------------------------------------------------------------------
c
c      By following the characteristics of the flow in the 2-direction,
c  determine values for "v3" and "b3" ("vsnm1" and "bsnm1") to be used
c  in evaluating "emf1".
c
       src = 0.0
c
c       Compute upwinded b2 and v2 in the 3-direction and use them to
c       compute the wave speeds for the chracteristic cones for the MoC
c
       do 59 j=js,je+1
         do 59 i=is,ie
             vfl   (k) = 0.5 * (v3(i,j,k) + v3(i,j-1,k)) - vg3(k)
             vt    (k) = v2(i,j,k) - vg2(j)
             bt    (k) = b2(i,j,k)
57         continue
c           call x3int1d ( bt, vfl, iordb2, istpb2, i, j
c     1                 , g31b, g31bi, g32a, g32ai, bint )
c           call x3int1d ( vt, vfl, iords2, istps2, i, j
c     1                 , g31b, g31bi, g32a, g32ai, vint )

c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
             dqm      = ( vt(k  ) - vt(km1) ) * dx3bi(k  )
             dqp      = ( vt(kp1) - vt(k  ) ) * dx3bi(kp1)
             dv (k  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )

             dqm      = ( bt(k  ) - bt(km1) ) * dx3bi(k  )
             dqp      = ( bt(kp1) - bt(k  ) ) * dx3bi(kp1)
             db (k  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )
 857       continue
c
c  2.  Perform an upwinded interpolation of "q" to the time-centred
c      bases of the characteristics.
c
           fact = dt * g31bi(i) * g32ai(j)
             qv1    = vt(km1) + dx3a(km1) * dv (km1)
             qv2    = vt(k  ) - dx3a(k  ) * dv (k  )
             qb1    = bt(km1) + dx3a(km1) * db (km1)
             qb2    = bt(k  ) - dx3a(k  ) * db (k  )
c  
             xi     = vfl(k) * fact
             q3     = sign ( haf, xi )
             vint(k)= ( 0.5 + q3 ) * ( qv1 - xi * dv (km1) )
     1              + ( 0.5 - q3 ) * ( qv2 - xi * dv (k  ) )
             bint(k)= ( 0.5 + q3 ) * ( qb1 - xi * db (km1) )
     1              + ( 0.5 - q3 ) * ( qb2 - xi * db (k  ) )
 858       continue
c
             emf1(i,j,k) = vint(k)
             emf2(i,j,k) = bint(k)
58         continue
59     continue
c
c      Select an effective density and determine the characteristic
c  velocity using upwinded values for each characteristic in the
c  2-direction.
c
         do 90 i=is,ie
           do 60 j=js,je+1
             vave (j) = emf1(i,j,k)
             bave (j) = emf2(i,j,k)
c            vave (j) = 0.5 * ( v2(i,j,k) + v2(i,j,km1) ) - vg2(j)
c            bave (j) = 0.5 * ( b2(i,j,k) + b2(i,j,km1) )
             absb     = abs ( bave(j) )
             aave (j) = 0.5 * absb * ( srd3(i,j,k) + srd3(i,j-1,k) )
     1                             / ( srd3(i,j,k) * srd3(i,j-1,k) )
             sgnp     = sign ( haf, vave(j) - aave(j) )
             sgnm     = sign ( haf, vave(j) + aave(j) )
             srdp (j) = ( 0.5 + sgnp ) * srd3(i,j-1,k)
     1                + ( 0.5 - sgnp ) * srd3(i,j,  k)
             srdm (j) = ( 0.5 + sgnm ) * srd3(i,j-1,k)
     1                + ( 0.5 - sgnm ) * srd3(i,j,  k)
             srdpi(j) = 1.0 / srdp(j)
             srdmi(j) = 1.0 / srdm(j)
             vchp (j) = vave(j) - absb * srdpi(j)
             vchm (j) = vave(j) + absb * srdmi(j)
60         continue
c
c      Interpolate 1-D vectors of "v3" and "b3" in the 2-direction to
c  the footpoints of both characteristics.
c
           do 70 j=js-2,je+2
             vtmp(j) = v3(i,j,k) - vg3(k)
             btmp(j) = b3(i,j,k)
70         continue
c           call x2zc1d ( vtmp, vchp, vchm, iords3, istps3, k, i
c     1                 , g2b, g2bi, vpch, vmch )
c           call x2zc1d ( btmp, vchp, vchm, iordb3, istpb3, k, i
c     1                 , g2b, g2bi, bpch, bmch )
c
           do 870 j=js-1,je+1
             dqm      = ( vtmp(j  ) - vtmp(j-1) ) * dx2bi(j  )
             dqp      = ( vtmp(j+1) - vtmp(j  ) ) * dx2bi(j+1)
             dv (j  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )

             dqm     = ( btmp(j  ) - btmp(j-1) ) * dx2bi(j  )
             dqp     = ( btmp(j+1) - btmp(j  ) ) * dx2bi(j+1)
             db (j  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )
 870       continue
c
c  2.  Perform an upwinded interpolation of "q" to the time-centred
c      bases of the characteristics.
c
           fact = dt * g2bi(i)
           do 880 j=js,je+1
             qv1    = vtmp(j-1) + dx2a(j-1) * dv (j-1)
             qv2    = vtmp(j  ) - dx2a(j  ) * dv (j  )
             qb1    = btmp(j-1) + dx2a(j-1) * db (j-1)
             qb2    = btmp(j  ) - dx2a(j  ) * db (j  )
c
             xi     = vchp(j) * fact
             q3     = sign ( haf, xi )
             vpch(j)= ( 0.5 + q3 ) * ( qv1 - xi * dv (j-1) )
     1              + ( 0.5 - q3 ) * ( qv2 - xi * dv (j  ) )
             bpch(j)= ( 0.5 + q3 ) * ( qb1 - xi * db (j-1) )
     1              + ( 0.5 - q3 ) * ( qb2 - xi * db (j  ) )
c
             xi     = vchm(j) * fact
             q3     = sign ( haf, xi )
             vmch(j)= ( 0.5 + q3 ) * ( qv1 - xi * dv (j-1) )
     1              + ( 0.5 - q3 ) * ( qv2 - xi * dv (j  ) )
             bmch(j)= ( 0.5 + q3 ) * ( qb1 - xi * db (j-1) )
     1              + ( 0.5 - q3 ) * ( qb2 - xi * db (j  ) )

 880       continue
c
c      Evaluate "vsnm1" and "bsnm1" by solving the characteristic
c  equations.  The source term is non-zero for RTP and ZRP coordinates
c  since dg32/dx2 .ne. 0.  Compute the two terms in "emf1".
c
           q1 = dt * g2bi(i)
           do 80 j=js,je+1
             q2           = sign ( one, bave(j) )
             if(lgeom .ne. 1) then
              src          = q1 * dg32bd2(j) * g32ai(j) * bave(j)
     1                     * ( b3  (i,j,k)    + b3  (i,j-1,k)    )
     1                     / ( srd3(i,j,k)**2 + srd3(i,j-1,k)**2 )
             endif ! lgeom
             vsnm1(  j  ) = ( vpch (j) * srdp (j) + vmch (j) * srdm (j)
     1                      + q2 * ( bpch(j) - bmch(j) ) )
     3                    / ( srdp (j) + srdm (j) ) + src
             bsnm1(  j  ) = ( bpch (j) * srdpi(j) + bmch (j) * srdmi(j)
     1                      + q2 * ( vpch(j) - vmch(j) ) )
     3                    / ( srdpi(j) + srdmi(j) )
             vsnm1(  j  ) = vsnm1(j) * bave(j)
             bsnm1(  j  ) = bsnm1(j) * vave(j)
c
             vsnp1(i,j,k) = 0.5*(  vsnp1(i,j,k) + bsnm1( j ) )
             bsnp1(i,j,k) = 0.5*(  vsnm1( j ) + bsnp1(i,j,k) )
80         continue
90     continue
c
c     END OF I LOOP
c
100    continue
c
c-----------------------------------------------------------------------
c
c      Set boundary values for "term1" and "term2".
c
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1
#endif

       call bvalemf1 ( vsnp1, bsnp1 )
c
c  Wait for communications to complete.
c
#ifdef MPI_USED
       if(nreq .ne. 0)
     .     call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c      Compute "emf1" for all 1-edges, including the ghost zones.
c
         do 120 j=js-2,je+3
           do 110 i=is-2,ie+2
            if(xvgrid) then
             emf1(i,j,k) = ( vsnp1(i,j,k) - bsnp1(i,j,k) )
     1                   * dx1ah(i)
            else
             emf1(i,j,k) = ( vsnp1(i,j,k) - bsnp1(i,j,k) )
     1                   * dx1a (i)
            endif
110        continue
120      continue
130    continue
c
c-----------------------------------------------------------------------
c---- 2.  emf2 ---------------------------------------------------------
c-----------------------------------------------------------------------
c
c     BEGINNING OF FIRST J LOOP
c
       do 180 j=js,je
c
c      By following the characteristics of the flow in the 1-direction,
c  determine values for "v3" and "b3" ("vsnp1" and "bsnp1") to be used
c  in evaluating "emf2".
c
         src = 0.0
c
c      Compute upwinded b1 and v1 in the 3-direction and use them to
c      compute the wave speeds for the chracteristic cones for the MoC.
c
         do 139 i=is,ie+1
             vfl   (k) = 0.5*(v3(i,j,k) + v3(i-1,j,k)) - vg3(k)
             vt    (k) = v1(i,j,k) - vg1(i)
             bt    (k) = b1(i,j,k)
137        continue
c
c
c           call x3int1d ( bt, vfl, iordb1, istpb1, i, j
c     1                  , g31b, g31bi, g32a, g32ai, bint )
c           call x3int1d ( vt, vfl, iordb1, istpb1, i, j
c     1                  , g31b, g31bi, g32a, g32ai, vint )

c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
             dqm      = ( vt(k  ) - vt(km1) ) * dx3bi(k  )
             dqp      = ( vt(kp1) - vt(k  ) ) * dx3bi(kp1)
             dv (k  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )

             dqm      = ( bt(k  ) - bt(km1) ) * dx3bi(k  )
             dqp      = ( bt(kp1) - bt(k  ) ) * dx3bi(kp1)
             db (k  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )
 937       continue
c
c  2.  Perform an upwinded interpolation of "q" to the time-centred
c      bases of the characteristics.
c
           fact = dt * g31bi(i) * g32ai(j)
             qv1    = vt(km1) + dx3a(km1) * dv (km1)
             qv2    = vt(k  ) - dx3a(k  ) * dv (k  )
             qb1    = bt(km1) + dx3a(km1) * db (km1)
             qb2    = bt(k  ) - dx3a(k  ) * db (k  )
c  
             xi     = vfl(k) * fact
             q3     = sign ( haf, xi )
             vint(k)= ( 0.5 + q3 ) * ( qv1 - xi * dv (km1) )
     1              + ( 0.5 - q3 ) * ( qv2 - xi * dv (k  ) )
             bint(k)= ( 0.5 + q3 ) * ( qb1 - xi * db (km1) )
     1              + ( 0.5 - q3 ) * ( qb2 - xi * db (k  ) )
 938       continue
c
             v1intk(k,i) = vint(k)
             b1intk(k,i) = bint(k)
138        continue
139      continue
c
c      Select an effective density and determine the characteristic
c  velocity using upwinded values for each characteristic in the
c  1-direction.
c
           do 140 i=is,ie+1
             vave (i) = v1intk(k,i)
             bave (i) = b1intk(k,i)
c            vave (i) = 0.5 * ( v1(i,j,k) + v1(i,j,km1) ) - vg1(i)
c            bave (i) = 0.5 * ( b1(i,j,k) + b1(i,j,km1) )
             absb     = abs ( bave(i) )
             aave (i) = 0.5 * absb * ( srd3(i,j,k) + srd3(i-1,j,k) )
     1                             / ( srd3(i,j,k) * srd3(i-1,j,k) )
             sgnp     = sign ( haf, vave(i) - aave(i) )
             sgnm     = sign ( haf, vave(i) + aave(i) )
             srdp (i) = ( 0.5 + sgnp ) * srd3(i-1,j,k)
     1                + ( 0.5 - sgnp ) * srd3(i  ,j,k)
             srdm (i) = ( 0.5 + sgnm ) * srd3(i-1,j,k)
     1                + ( 0.5 - sgnm ) * srd3(i  ,j,k)
             srdpi(i) = 1.0 / srdp(i)
             srdmi(i) = 1.0 / srdm(i)
             vchp (i) = vave(i) - absb * srdpi(i)
             vchm (i) = vave(i) + absb * srdmi(i)
140        continue
c
c      Interpolate 1-D vectors of "v3" and "b3" in the 1-direction to
c  the footpoints of both characteristics.
c
           do 150 i=is-2,ie+2
             vtmp(i) = v3(i,j,k) - vg3(k)
             btmp(i) = b3(i,j,k)
150        continue
c           call x1zc1d ( vtmp, vchp, vchm, iords3, istps3, j, k
c     1                 , vpch, vmch )
c           call x1zc1d ( btmp, vchp, vchm, iordb3, istpb3, j, k
c     1                 , bpch, bmch )
c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
           do 950 i=is-1,ie+1
             dqm      = ( vtmp(i  ) - vtmp(i-1) ) * dx1bi(i  )
             dqp      = ( vtmp(i+1) - vtmp(i  ) ) * dx1bi(i+1)
             dv (i  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )

             dqm      = ( btmp(i  ) - btmp(i-1) ) * dx1bi(i  )
             dqp      = ( btmp(i+1) - btmp(i  ) ) * dx1bi(i+1)
             db (i  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )
 950       continue
c
c  2.  Perform an upwinded interpolation of "q" to the time-centred
c      bases of the characteristics.
c
           do 960 i=is,ie+1
             qv1    = vtmp(i-1) + dx1a(i-1) * dv (i-1)
             qv2    = vtmp(i  ) - dx1a(i  ) * dv (i  )
             qb1    = btmp(i-1) + dx1a(i-1) * db (i-1)
             qb2    = btmp(i  ) - dx1a(i  ) * db (i  )
c
             xi     = vchp(i) * dt
             q3     = sign ( haf, xi )
             vpch(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) )
     1              + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
             bpch(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) )
     1              + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
c
             xi     = vchm(i) * dt
             q3     = sign ( haf, xi )
             vmch(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) )
     1              + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
             bmch(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) )
     1              + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
 960       continue
c
c      Evaluate "vsnp1" and "bsnp1" by solving the characteristic
c  equations.  The source term is non-zero for RTP coordinates since
c  dg31/dx1 = 1.0.
c
           do 160 i=is,ie+1
             q2           = sign ( one, bave(i) )
             if(lgeom .eq. 3) then
              src          = dt * dg31bd1(i) * g31ai(i) * bave(i)
     2                     * ( b3  (i,j,k)    + b3  (i-1,j,k)    )
     1                     / ( srd3(i,j,k)**2 + srd3(i-1,j,k)**2 )
             endif ! lgeom
             vsnp1(i,j,k) = ( vpch (i) * srdp (i) + vmch (i) * srdm (i)
     1                       + q2 * ( bpch(i) - bmch(i) ) )
     3                      / ( srdp (i) + srdm (i) ) + src
             bsnp1(i,j,k) = ( bpch (i) * srdpi(i) + bmch (i) * srdmi(i)
     1                       + q2 * ( vpch (i) - vmch (i) ) )
     3                      / ( srdpi(i) + srdmi(i) )
             vsnp1(i,j,k) = vsnp1(i,j,k) * bave(i)
             bsnp1(i,j,k) = bsnp1(i,j,k) * vave(i)
160        continue
c
170      continue
c
c        END OF FIRST J LOOP
c
180    continue
c
c-----------------------------------------------------------------------
c
c
c       BEGINNING OF SECOND J LOOP
c
c
       do 230 j=js,je
c
c      By following the characteristics of the flow in the 3-direction,
c  determine values for "v1" and "b1" ("vsnm1" and "bsnm1") to be used
c  in evaluating "emf2".
c
c       Compute upwinded b3 and v3 in the 1-direction and use them to
c       compute the wave speeds for the chracteristic cones for the MoC
c
           do 188 i=is-2,ie+2
             vfl   (i) = 0.5*(v1(i,j,k) + v1(i,j,km1)) - vg1(i)
             vt    (i) = v3(i,j,k) - vg3(k)
             bt    (i) = b3(i,j,k)
188        continue
c
c           call x1int1d ( bt, vfl, iordb1, istpb1, j, k
c     1                  , bint )
c           call x1int1d ( vt, vfl, iords1, istps1, j, k
c     1                  , vint )
c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
           do 988 i=is-1,ie+1
             dqm      = ( vt(i  ) - vt(i-1) ) * dx1bi(i  )
             dqp      = ( vt(i+1) - vt(i  ) ) * dx1bi(i+1)
             dv (i  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )

             dqm      = ( bt(i  ) - bt(i-1) ) * dx1bi(i  )
             dqp      = ( bt(i+1) - bt(i  ) ) * dx1bi(i+1)
             db (i  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )
 988       continue
c
c  2.  Perform an upwinded interpolation of "q" to the time-centred
c      bases of the characteristics.
c
           do 987 i=is,ie+1
             qv1    = vt(i-1) + dx1a(i-1) * dv (i-1)
             qv2    = vt(i  ) - dx1a(i  ) * dv (i  )
             qb1    = bt(i-1) + dx1a(i-1) * db (i-1)
             qb2    = bt(i  ) - dx1a(i  ) * db (i  )
c
             xi     = vfl(i) * dt
             q3     = sign ( haf, xi )
             vint(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) )
     1              + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
             bint(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) )
     1              + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
 987       continue
c
           do 187 i=is,ie+1
             v3inti(k,i) = vint(i)
             b3inti(k,i) = bint(i)
187        continue
189      continue
c
       do 186 i=is,ie+1
           vsnp1(i,j,k) = 0.5*( vsnp1(i,j,k)
     1                       + b1(i,j,k)*v3inti(k,i) )
           bsnp1(i,j,k) = 0.5*( bsnp1(i,j,k)
     1                       + ( v1(i,j,k) - vg1(i) )*b3inti(k,i))
186     continue
c
c       END OF SECOND J LOOP
c
230    continue
c
c-----------------------------------------------------------------------
c
c      Set boundary values for "term1" and "term2".
c
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1
#endif
       call bvalemf2 ( vsnp1, bsnp1 )
c
c  Wait for communications to complete.
c
#ifdef MPI_USED
       if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c      Compute "emf2" for all 2-edges, including the ghost zones.
c
         do 250 j=js-2,je+2
           do 240 i=is-2,ie+3
            if(xvgrid) then
             emf2(i,j,k) = ( vsnp1(i,j,k) - bsnp1(i,j,k) )
     1                   * dx2ah(j) * g2ah(i)
            else
             emf2(i,j,k) = ( vsnp1(i,j,k) - bsnp1(i,j,k) )
     1                   * dx2a (j) * g2a (i)
            endif
240        continue
250      continue
260    continue
c
c-----------------------------------------------------------------------
c---- 3.  emf3 ---------------------------------------------------------
c-----------------------------------------------------------------------
c
c      BEGINNING OF K LOOP
c
c
c      By following the characteristics of the flow in the 2-direction,
c  determine values for "v1" and "b1" ("vsnp1" and "bsnp1") to be used
c  in evaluating "emf3".
c
c      Compute upwinded b2 and v2 in the 1-direction and use them to
c      compute the wave speeds for the chracteristic cones for MoC.
c
         do 269 j=js,je+1
           do 267 i=is-2,ie+2
             vfl   (i) = 0.5 * (v1(i,j,k) + v1(i,j-1,k)) - vg1(i)
             vt    (i) = v2(i,j,k) - vg2(j)
             bt    (i) = b2(i,j,k)
267        continue
c           call x1int1d ( bt, vfl, iordb2, istpb2, j, k
c     1                  , bint )
c           call x1int1d ( vt, vfl, iords2, istps2, j, k
c     1                  , vint )
c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
           do 1067 i=is-1,ie+1
             dqm      = ( vt(i  ) - vt(i-1) ) * dx1bi(i  )
             dqp      = ( vt(i+1) - vt(i  ) ) * dx1bi(i+1)
             dv (i  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )

             dqm      = ( bt(i  ) - bt(i-1) ) * dx1bi(i  )
             dqp      = ( bt(i+1) - bt(i  ) ) * dx1bi(i+1)
             db (i  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )
 1067      continue
c
c  2.  Perform an upwinded interpolation of "q" to the time-centred
c      bases of the characteristics.
c
           do 1068 i=is,ie+1
             qv1    = vt(i-1) + dx1a(i-1) * dv (i-1)
             qv2    = vt(i  ) - dx1a(i  ) * dv (i  )
             qb1    = bt(i-1) + dx1a(i-1) * db (i-1)
             qb2    = bt(i  ) - dx1a(i  ) * db (i  )
c
             xi     = vfl(i) * dt
             q3     = sign ( haf, xi )
             vint(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) )
     1              + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
             bint(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) )
     1              + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
 1068      continue
c
           do 268 i=is,ie+1
             v2inti(i,j) = vint(i)
             b2inti(i,j) = bint(i)
268        continue
269      continue
c
c      Select an effective density and determine the characteristic
c  velocity using upwinded values for each characteristic in the
c  2-direction.
c
         do 300 i=is,ie+1
           do 270 j=js,je+1
             vave (j) = v2inti(i,j)
             bave (j) = b2inti(i,j)
c            vave (j) = 0.5 * ( v2(i,j,k) + v2(im1,j,k) ) - vg2(j)
c            bave (j) = 0.5 * ( b2(i,j,k) + b2(im1,j,k) )
             absb     = abs ( bave(j) )
             aave (j) = 0.5 * absb * ( srd1(i,j,k) + srd1(i,j-1,k) )
     1                             / ( srd1(i,j,k) * srd1(i,j-1,k) )
             sgnp     = sign ( haf, vave(j) - aave(j) )
             sgnm     = sign ( haf, vave(j) + aave(j) )
             srdp (j) = ( 0.5 + sgnp ) * srd1(i,j-1,k)
     1                + ( 0.5 - sgnp ) * srd1(i,j  ,k)
             srdm (j) = ( 0.5 + sgnm ) * srd1(i,j-1,k)
     1                + ( 0.5 - sgnm ) * srd1(i,j  ,k)
             srdpi(j) = 1.0 / srdp(j)
             srdmi(j) = 1.0 / srdm(j)
             vchp (j) = vave(j) - absb * srdpi(j)
             vchm (j) = vave(j) + absb * srdmi(j)
270        continue
c
c      Interpolate 1-D vectors of "v1" and "b1" in the 2-direction to
c  the footpoints of both characteristics.
c
           do 280 j=js-2,je+2
             vtmp(j) = v1(i,j,k) - vg1(i)
             btmp(j) = b1(i,j,k)
280        continue
c           call x2zc1d ( vtmp, vchp, vchm, iords1, istps1, k, i
c     1                 , g2a, g2ai, vpch, vmch )
c           call x2zc1d ( btmp, vchp, vchm, iordb1, istpb1, k, i
c     1                 , g2a, g2ai, bpch, bmch )
c
           do 1080 j=js-1,je+1
             dqm      = ( vtmp(j  ) - vtmp(j-1) ) * dx2bi(j  )
             dqp      = ( vtmp(j+1) - vtmp(j  ) ) * dx2bi(j+1)
             dv (j  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )

             dqm     = ( btmp(j  ) - btmp(j-1) ) * dx2bi(j  )
             dqp     = ( btmp(j+1) - btmp(j  ) ) * dx2bi(j+1)
             db (j  ) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )
 1080      continue
c
c  2.  Perform an upwinded interpolation of "q" to the time-centred
c      bases of the characteristics.
c
           fact = dt * g2ai(i)
           do 1090 j=js,je+1
             qv1    = vtmp(j-1) + dx2a(j-1) * dv (j-1)
             qv2    = vtmp(j  ) - dx2a(j  ) * dv (j  )
             qb1    = btmp(j-1) + dx2a(j-1) * db (j-1)
             qb2    = btmp(j  ) - dx2a(j  ) * db (j  )
c
             xi     = vchp(j) * fact
             q3     = sign ( haf, xi )
             vpch(j)= ( 0.5 + q3 ) * ( qv1 - xi * dv (j-1) )
     1              + ( 0.5 - q3 ) * ( qv2 - xi * dv (j  ) )
             bpch(j)= ( 0.5 + q3 ) * ( qb1 - xi * db (j-1) )
     1              + ( 0.5 - q3 ) * ( qb2 - xi * db (j  ) )
c
             xi     = vchm(j) * fact
             q3     = sign ( haf, xi )
             vmch(j)= ( 0.5 + q3 ) * ( qv1 - xi * dv (j-1) )
     1              + ( 0.5 - q3 ) * ( qv2 - xi * dv (j  ) )
             bmch(j)= ( 0.5 + q3 ) * ( qb1 - xi * db (j-1) )
     1              + ( 0.5 - q3 ) * ( qb2 - xi * db (j  ) )

 1090      continue
c
c      Evaluate "vsnp1" and "bsnp1" by solving the characteristic
c  equations.  There is no source term since the metric factor "g1" has
c  no explicit dependence on x2.
c
           do 290 j=js,je+1
             q2           = sign ( one, bave(j) )
             vsnp1(i,j,k) = ( vpch (j) * srdp (j) + vmch (j) * srdm (j)
     1                      + q2 * ( bpch(j) - bmch(j) ) )
     2                    / ( srdp (j) + srdm (j) )
             bsnp1(i,j,k) = ( bpch (j) * srdpi(j) + bmch (j) * srdmi(j)
     1                      + q2 * ( vpch(j) - vmch(j) ) )
     2                    / ( srdpi(j) + srdmi(j) )
             vsnp1(i,j,k) = vsnp1(i,j,k) * bave(j)
             bsnp1(i,j,k) = bsnp1(i,j,k) * vave(j)
290        continue
c
300    continue
c
c-----------------------------------------------------------------------
c
c      By following the characteristics of the flow in the 1-direction,
c  determine values for "v2" and "b2" ("vsnm1" and "bsnm1") to be used
c  in evaluating "emf3".
c
       src = 0.0
c
c      Compute upwinded b1 and v1 in the 2-direction and use them to
c      compute the wave speeds for the chracteristic cones for Moc
c
       do 319 i=is,ie+1
         do 317 j=js-2,je+2
           vfl   (j) = 0.5 * (v2(i,j,k) + v2(i-1,j,k)) - vg2(j)
           vt    (j) = v1(i,j,k) - vg1(i)
           bt    (j) = b1(i,j,k)
317      continue
c         call x2int1d ( bt, vfl, iordb1, istpb1, k, i
c     1                , g2b, g2bi, bint )
c         call x2int1d ( vt, vfl, iords1, istps1, k, i
c     1                , g2b, g2bi, vint )
c
         do 1117 j=js-1,je+1
           dqm      = ( vt(j  ) - vt(j-1) ) * dx2bi(j  )
           dqp      = ( vt(j+1) - vt(j  ) ) * dx2bi(j+1)
           dv (j  ) = max ( dqm * dqp, zro )
     1              * sign ( one, dqm + dqp )
     2              / max ( abs ( dqm + dqp ), tiny )

           dqm     = ( bt(j  ) - bt(j-1) ) * dx2bi(j  )
           dqp     = ( bt(j+1) - bt(j  ) ) * dx2bi(j+1)
           db (j  ) = max ( dqm * dqp, zro )
     1              * sign ( one, dqm + dqp )
     2              / max ( abs ( dqm + dqp ), tiny )
 1117    continue
c
c  2.  Perform an upwinded interpolation of "q" to the time-centred
c      bases of the characteristics.
c
         fact = dt * g2bi(i)
         do 1118 j=js,je+1
           qv1    = vt(j-1) + dx2a(j-1) * dv (j-1)
           qv2    = vt(j  ) - dx2a(j  ) * dv (j  )
           qb1    = bt(j-1) + dx2a(j-1) * db (j-1)
           qb2    = bt(j  ) - dx2a(j  ) * db (j  )
c
           xi     = vfl(j) * fact
           q3     = sign ( haf, xi )
           vint(j)= ( 0.5 + q3 ) * ( qv1 - xi * dv (j-1) )
     1            + ( 0.5 - q3 ) * ( qv2 - xi * dv (j  ) )
           bint(j)= ( 0.5 + q3 ) * ( qb1 - xi * db (j-1) )
     1            + ( 0.5 - q3 ) * ( qb2 - xi * db (j  ) )
 1118    continue
c
         do 318 j=js,je+1
           b1intj(i,j) = bint(j)
           v1intj(i,j) = vint(j)
318      continue
319    continue
c
c      Select an effective density and determine the characteristic
c  velocity using upwinded values for each characteristic in the
c  1-direction.
c
       do 350 j=js,je+1
         do 320 i=is,ie+1
           vave (i) = v1intj(i,j)
           bave (i) = b1intj(i,j)
c          vave (i) = 0.5 * ( v1(i,j,k) + v1(i,jm1,k) ) - vg1(i)
c          bave (i) = 0.5 * ( b1(i,j,k) + b1(i,jm1,k) )
           absb     = abs ( bave(i) )
           aave (i) = 0.5 * absb * ( srd2(i,j,k) + srd2(i-1,j,k) )
     1                           / ( srd2(i,j,k) * srd2(i-1,j,k) )
           sgnp     = sign ( haf, vave(i) - aave(i) )
           sgnm     = sign ( haf, vave(i) + aave(i) )
           srdp (i) = ( 0.5 + sgnp ) * srd2(i-1,j,k)
     1              + ( 0.5 - sgnp ) * srd2(i  ,j,k)
           srdm (i) = ( 0.5 + sgnm ) * srd2(i-1,j,k)
     1              + ( 0.5 - sgnm ) * srd2(i  ,j,k)
           srdpi(i) = 1.0 / srdp(i)
           srdmi(i) = 1.0 / srdm(i)
           vchp (i) = vave(i) - absb * srdpi(i)
           vchm (i) = vave(i) + absb * srdmi(i)
320      continue
c
c      Interpolate 1-D vectors of "v2" and "b2" in the 1-direction to
c  the footpoints of both characteristics.
c
         do 330 i=is-2,ie+2
           vtmp(i) = v2(i,j,k) - vg2(j)
           btmp(i) = b2(i,j,k)
330      continue
c         call x1zc1d ( vtmp, vchp, vchm, iords2, istps2, j, k
c     1               , vpch, vmch )
c         call x1zc1d ( btmp, vchp, vchm, iordb2, istpb2, j, k
c     1               , bpch, bmch )
c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
         do 1130 i=is-1,ie+1
           dqm      = ( vtmp(i  ) - vtmp(i-1) ) * dx1bi(i  )
           dqp      = ( vtmp(i+1) - vtmp(i  ) ) * dx1bi(i+1)
           dv (i  ) = max ( dqm * dqp, zro )
     1              * sign ( one, dqm + dqp )
     2              / max ( abs ( dqm + dqp ), tiny )

           dqm      = ( btmp(i  ) - btmp(i-1) ) * dx1bi(i  )
           dqp      = ( btmp(i+1) - btmp(i  ) ) * dx1bi(i+1)
           db (i  ) = max ( dqm * dqp, zro )
     1              * sign ( one, dqm + dqp )
     2              / max ( abs ( dqm + dqp ), tiny )
 1130    continue
c
c  2.  Perform an upwinded interpolation of "q" to the time-centred
c      bases of the characteristics.
c
         do 1140 i=is,ie+1
           qv1    = vtmp(i-1) + dx1a(i-1) * dv (i-1)
           qv2    = vtmp(i  ) - dx1a(i  ) * dv (i  )
           qb1    = btmp(i-1) + dx1a(i-1) * db (i-1)
           qb2    = btmp(i  ) - dx1a(i  ) * db (i  )
c          
           xi     = vchp(i) * dt
           q3     = sign ( haf, xi )
           vpch(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) )
     1            + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
           bpch(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) )
     1            + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
c
           xi     = vchm(i) * dt
           q3     = sign ( haf, xi )
           vmch(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) )
     1            + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
           bmch(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) )
     1            + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
1140    continue
c
c      Evaluate "vsnm1" and "bsnm1" by solving the characteristic
c  equations.  The source term is non-zero for RTP coordinates since
c  dg2/dx1 = 1.0.  Compute the two terms in "emf3".
c
         do 340 i=is,ie+1
           q2           = sign ( one, bave(i) )
           if(lgeom .eq. 3) then
            src          = dt * dg2bd1(i) * g2ai(i) * bave(i)
     1                   * ( b2  (i,j,k)    + b2  (i-1,j,k)    )
     1                   / ( srd2(i,j,k)**2 + srd2(i-1,j,k)**2 )
           endif ! lgeom
           vsnm1(i    ) = ( vpch (i) * srdp (i) + vmch (i) * srdm (i)
     1                    + q2 * ( bpch(i) - bmch(i) ) )
     3                  / ( srdp (i) + srdm (i) ) + src
           bsnm1(i    ) = ( bpch (i)* srdpi(i) + bmch (i)* srdmi(i)
     1                    + q2 * ( vpch(i) - vmch(i) ) )
     3                  / ( srdpi(i) + srdmi(i) )
c
           vsnm1(i    ) = vsnm1(i) * bave(i)
           bsnm1(i    ) = bsnm1(i) * vave(i)
c
           vsnp1(i,j,k) = 0.5*(  vsnp1(i,j,k) + bsnm1(i) )
           bsnp1(i,j,k) = 0.5*(  vsnm1(i)     + bsnp1(i,j,k) )
340      continue
350    continue
c
c       END OF K LOOP
c
360    continue
c
c-----------------------------------------------------------------------
c
c      Set boundary values for "term1" and "term2".
c
#ifdef MPI_USED
       nreq = 0
       nsub = nsub + 1
#endif
       call bvalemf3 ( vsnp1, bsnp1 )
c
c  Wait for communications to complete.
c
#ifdef MPI_USED
      if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c      Compute "emf3" for all 3-edges, including the ghost zones.
c
         do 380 j=js-2,je+3
           do 370 i=is-2,ie+3
            if(xvgrid) then
             emf3(i,j,k) = ( vsnp1(i,j,k) - bsnp1(i,j,k) )
     1                   * dx3ah(k) * g31ah(i) * g32ah(j)
            else
             emf3(i,j,k) = ( vsnp1(i,j,k) - bsnp1(i,j,k) )
     1                   * dx3a (k) * g31a (i) * g32a (j)
            endif
370        continue
380      continue
390    continue
c
       return
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                 H S M O C                 \\\\\\\\\\
c
c=======================================================================
