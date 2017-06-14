c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////               L O R E N T Z               \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine lorentz_1D(ibeg,iend,jbeg,jend,kbeg,kend
     &                   ,u1,u2,u3,w1,w2,w3)
c
c    dac:zeus3d.lorentz <----- MoC estimate of transverse lorentz forces
c                                                          october, 1992
c
c    written by: David Clarke
c    modified 1: Mordecai-Mark Mac Low, December 1997 - March 1998
c                rewritten for ZEUS-MP 
c    modified 2: John Hayes, October 2005; created lorentz_1D clone
c                which assumes symmetry about the J and K axes; after
c                ZEUS3D.
c
c  PURPOSE:  Uses the Method of Characteristics (MoC) to compute the
c  transverse components of the Lorentz force and then accelerate the
c  velocities accordingly.  After a suggestion by John Hawley, these
c  source terms are computed separately from the MoC estimate of the
c  emfs (MOCEMFS) and before the magnetic field update.  This improves
c  stability of multi-dimensional Alfven waves.
c
c  MOCEMFS solves the induction equation without the use of operator
c  splitting.  Thus, the relevant characteristic velocities are v+/-va
c  and includes the flow velocity v (from the advection term v grad(B))
c  and the alfven velocity va (from the induction term va grad(B)).
c  Conversely, Euler's equation is operator split, and the advection
c  term (v grad(v)) is handled elsewhere (MOMX*).  Thus, the effective
c  characteristic velocities when estimating the transverse Lorentz
c  terms is simply +/-va (from the induction term va grad(v)).
c
c  See comments in MOCEMFS for further ideas regarding the Method of
c  Characteristics.
c
c  INPUT VARIABLES:
c
c    u[1,2,3] = old velocities
c
c  OUTPUT VARIABLES:
c
c    w[1,2,3] = velocities from ibeg to iend, jbeg to jend, kbeg to kend
c
c  LOCAL VARIABLES:
c    srd[n]    sqrt of spatially averaged density at [n]-face n=1,2,3
c    bave      spatially averaged density at edge
c    srdpi     1.0/srd[n] along the plus  alfven characteristic (A+)
c    srdmi     1.0/srd[n] along the minus alfven characteristic (A-)
c    valp      characteristic velocity along A+ (-va)
c    valm      characteristic velocity along A- (+va)
c    vpal      velocity interpolated to the time-centred footpoint of A+
c    vmal      velocity interpolated to the time-centred footpoint of A-
c    bpal      B-field  interpolated to the time-centred footpoint of A+
c    bmal      B-field  interpolated to the time-centred footpoint of A-
c    bstar     MoC estimate of b[n] used to evaluate st[n], n=1,2,3
c    st[n]     Transverse Lorentz force on the [n]-velocity, n=1,2,3
c
c  EXTERNALS:
c    
c   -----> these are now inlined (M-MML): X1ZC1D  , X2ZC1D  , X3ZC1D
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
      integer  :: i       , j       , k      
      integer  :: ibeg,iend,jbeg,jend,kbeg,kend
      integer  :: kone,km1,kp1,jone,jm1  !asif/jch
c
      real(rl) :: u1(in,jn,kn), u2(in,jn,kn), u3(in,jn,kn),
     .            w1(in,jn,kn), w2(in,jn,kn), w3(in,jn,kn)
c
      real(rl) :: absb ,q1, qv1, qv2, qb1, qb2, q3, dqm, dqp,
     .            xi, fact,
     .            dv(ijkn), db(ijkn)
c
      real(rl) :: bave(ijkn), srdpi(ijkn), srdmi(ijkn),
     .            valp(ijkn), valm (ijkn), vtmp (ijkn),
     .            btmp(ijkn), vpal (ijkn), vmal (ijkn),
     .            bpal(ijkn), bmal (ijkn), bstar(ijkn)
c
      real(rl) :: st1(in,jn,kn), st2(in,jn,kn),
     .            st3(in,jn,kn)
c
C       equivalence   
C     .               ( bave    , w1da     ),
C     .               ( srdpi   , w1db     ),
C     1               ( srdmi   , w1dc     ), ( vtmp    , w1dd     ),
C     2               ( vpal    , w1de     ), ( vmal    , w1df     )
C     2               ( valp    , bpal    , w1dg     ),
C     1               ( valm    , bmal    , w1dh     ),
C     8               ( btmp    , bstar   , w1di     ),
C     9               ( db      , w1dj     ),  ( dv     , w1dk     )
c
c      Routine MOCEMFS is expecting to find "srd1", "srd2", and "srd3" 
c  intact in worker arrays "wd3d", "we3d", and "wf3d" respectively.
c  [but MOCEMFS is in ZEUS-3D. does HSMOC also need these? M-MML 9.12.97 
c       so far, yes 4.3.98]
c
CPS: prevent overwritten of momenta.
C       equivalence   ( st1     , w3da     ), ( st2     , w3db     )
C     1             , ( st3     , w3dc     )
CPS -- end
c
C       equivalence   ( srd1    , w3di     ), 
C     .               ( srd2    , w3dj     ),
C     1               ( srd3    , w3df     )
c
c      External statements
c
c
c
c-----------------------------------------------------------------------
c
c      Compute face-centred averages of density, initialise source terms
c  to zero.
c
      kone=0
      jone=0
c
      j   = js
      k   = ks
      km1 = k-kone
      jm1 = j-jone
c
      do 1 i=ibeg-1,iend+1
       srd1(i,j,k) = sqrt ( 0.5 * ( d (i,j,k) + d (i-1,j  ,k) ) )
       srd2(i,j,k) = sqrt ( 0.5 * ( d (i,j,k) + d (i  ,jm1,k) ) )
       srd3(i,j,k) = sqrt ( 0.5 * ( d (i,j,k) + d (i  ,j  ,km1) ) )
       st1 (i,j,k) = 0.0
       st2 (i,j,k) = 0.0
       st3 (i,j,k) = 0.0
1     continue
c
c-----------------------------------------------------------------------
c---- 1-force ----------------------------------------------------------
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c---- 2-force ----------------------------------------------------------
c-----------------------------------------------------------------------
c
c      By following the Alfven velocity in the 1-direction, evaluate
c  "bstar" from "b2" to estimate the 1-transverse Lorentz force.
c
c
c      Select an effective density and determine the Alfven speed for
c  each Alfven characteristic in the 1-direction.
c
           do 170 i=ibeg,iend+1
             bave (i) = 0.5 * ( b1(i,j,k) + b1(i,jm1,k) )
             absb     = abs ( bave(i) )
             srdpi(i) = 1.0 / srd2(i  ,j,k)
             srdmi(i) = 1.0 / srd2(i-1,j,k)
             valp (i) =-absb * srdpi(i)
             valm (i) = absb * srdmi(i)
170        continue
c
c      Interpolate 1-D vectors of "u2" and "b2" in the 1-direction to
c  the footpoints of both Alfven characteristics.
c
           do 180 i=ibeg-2,iend+2
             vtmp(i) = u2(i,j,k) - vg2(j)
             btmp(i) = b2(i,j,k)
180        continue
c          call x1zc1d ( vtmp, valp, valm, iords2, istps2, j, k
c    1                 , vpal, vmal )
c          call x1zc1d ( btmp, valp, valm, iordb2, istpb2, j, k
c    1                 , bpal, bmal )

c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
           do 182 i=ibeg-1,iend+1
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
 182       continue
c
c  2.  Perform an upwinded interpolation of "q" to the time-centred
c      bases of the characteristics.
c
           do 184 i=ibeg,iend+1
             qv1    = vtmp(i-1) + dx1a(i-1) * dv (i-1)
             qv2    = vtmp(i  ) - dx1a(i  ) * dv (i  )
             qb1    = btmp(i-1) + dx1a(i-1) * db (i-1)
             qb2    = btmp(i  ) - dx1a(i  ) * db (i  )
c
             xi     = valp(i) * dt
             q3     = sign ( haf, xi )
             vpal(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) )
     1              + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
             bpal(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) )
     1              + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
c
             xi     = valm(i) * dt
             q3     = sign ( haf, xi )
             vmal(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) )
     1              + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
             bmal(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) )
     1              + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
184        continue

c
c      Evaluate "bstar" by solving the characteristic equations.
c
           do 190 i=ibeg,iend+1
             q1       = sign ( one, bave(i) )
             bstar(i) = ( bpal(i) * srdpi(i) + bmal(i) * srdmi(i)
     1                  + q1 * ( vpal(i) - vmal(i) ) ) * g2a(i)
     2                / ( srdpi(i) + srdmi(i) )
190        continue
c
c      Evaluate transverse Lorentz force.  Note that the metric
c  factor "g2a" has been absorbed into "bstar".
c
           do 200 i=ibeg,iend
             st2(i,j,k) = st2(i,j,k)
     2                  + ( bave (i+1) + bave (i) ) * g2bi(i)
     4                  * ( bstar(i+1) - bstar(i) ) * dx1ai(i)
200        continue
205      continue
210    continue
c
220    continue
c
c-----------------------------------------------------------------------
c---- 3-force ----------------------------------------------------------
c-----------------------------------------------------------------------
c
c
c      By following the Alfven velocity in the 1-direction, evaluate
c  "bstar" from "b3" to estimate the 1-transverse Lorentz force.
c
c      Select an effective density and determine the Alfven speed for
c  each Alfven characteristic in the 1-direction.
c
           do 230 i=ibeg,iend+1
             bave (i) = 0.5 * ( b1(i,j,k) + b1(i,j,km1) )
             absb     = abs ( bave(i) )
             srdpi(i) = 1.0 / srd3(i  ,j,k)
             srdmi(i) = 1.0 / srd3(i-1,j,k)
             valp (i) =-absb * srdpi(i)
             valm (i) = absb * srdmi(i)
230        continue
c
c      Interpolate 1-D vectors of "u3" and "b3" in the 1-direction to
c  the footpoints of both Alfven characteristics.
c
           do 240 i=ibeg-2,iend+2
             vtmp(i) = u3(i,j,k) - vg3(k)
             btmp(i) = b3(i,j,k)
240        continue
c           call x1zc1d ( vtmp, valp, valm, iords3, istps3, j, k
c     1                 , vpal, vmal )
c           call x1zc1d ( btmp, valp, valm, iordb3, istpb3, j, k
c     1                 , bpal, bmal )
c
c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
           do 242 i=ibeg-1,iend+1
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
 242       continue
c
c  2.  Perform an upwinded interpolation of "q" to the time-centred
c      bases of the characteristics.
c
           do 244 i=ibeg,iend+1
             qv1    = vtmp(i-1) + dx1a(i-1) * dv (i-1)
             qv2    = vtmp(i  ) - dx1a(i  ) * dv (i  )
             qb1    = btmp(i-1) + dx1a(i-1) * db (i-1)
             qb2    = btmp(i  ) - dx1a(i  ) * db (i  )
c
             xi     = valp(i) * dt
             q3     = sign ( haf, xi )
             vpal(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) )
     1              + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
             bpal(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) )
     1              + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
c
             xi     = valm(i) * dt
             q3     = sign ( haf, xi )
             vmal(i)= ( 0.5 + q3 ) * ( qv1 - xi * dv (i-1) )
     1              + ( 0.5 - q3 ) * ( qv2 - xi * dv (i  ) )
             bmal(i)= ( 0.5 + q3 ) * ( qb1 - xi * db (i-1) )
     1              + ( 0.5 - q3 ) * ( qb2 - xi * db (i  ) )
244        continue
c
c      Evaluate "bstar" by solving the characteristic equations.
c
           do 250 i=ibeg,iend+1
             q1       = sign ( one, bave(i) )
             bstar(i) = ( bpal(i) * srdpi(i) + bmal(i) * srdmi(i)
     1                  + q1 * ( vpal(i) - vmal(i) ) ) * g31a(i)
     2                / ( srdpi(i) + srdmi(i) )
250        continue
c
c      Evaluate transverse Lorentz force.  Note that the metric
c  factor "g31a" has been absorbed into "bstar".
c
           do 260 i=ibeg,iend
CJH         
            if(lgeom .eq. 3) then
             if(x1a(i) .eq. 0.0) then
              st3(i,j,k) = 0.0D0
              go to 260
             endif
            endif
             st3(i,j,k) = ( bave (i+1) + bave (i) ) * g31bi(i)
     4                  * ( bstar(i+1) - bstar(i) ) * dx1ai(i)
260        continue
270      continue
c
330    continue
c
c-----------------------------------------------------------------------
c
c      Accelerate velocities and set boundaries.
c
      q1 = 0.5 * dt
      do 340 i=ibeg,iend
       w1(i,j,k) = u1(i,j,k) + q1 * st1(i,j,k) / srd1(i,j,k)**2
       w2(i,j,k) = u2(i,j,k) + q1 * st2(i,j,k) / srd2(i,j,k)**2
       w3(i,j,k) = u3(i,j,k) + q1 * st3(i,j,k) / srd3(i,j,k)**2
340   continue
c
       return
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////               L O R E N T Z               \\\\\\\\\\
c
c=======================================================================
c
