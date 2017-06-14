c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 M O M X 3                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine momx3 ( ibeg,iend,jbeg,jend,kbeg,kend
     &                  , s1,s2,s3,mflx,
     .                    atwid1,atwid2,atwid3,atwidj1,atwidj2,atwidj3,
     .                    sflx,dq)
c
c    dac:zeus3d.momx3 <--------------- transports momenta in 3-direction
c                                                         february, 1990
c
c    written by: David Clarke
c    modified 1: November, 1992 by David Clarke; momenta are now updated
c                between and including i=is,ie, j=js,je, and k=ks,ke to
c                allow for proper treatment of periodic boundaries.
c    modified 2: Feb. 22, 1996 by Robert Fiedler; completely rewritten
c                for ZEUS-MP.
c    modified 3: July 6, 1996 by Robert Fiedler; blocked the k-loop
c                to improve performance on systems with a small TLB.
c    modified 4: Aug. 2, 1996 by Robert Fiedler; unrolled the i-loop
c                to improve performance on systems with a small cache.
c
c  PURPOSE:  Transports the three components of the momentum density in
c  the 3-direction using the consistent transport algorithm, including
c  the effects of grid compression.  The transported fluxes are thus
c  given by the mass fluxes times the time centred area of the control
c  volume faces times the interpolated velocities.  Interpolations are
c  performed in-line.
c
c  INPUT VARIABLES:
c    mflx    mass flux in 3-direction
c    s1      momentum density in 1-direction
c    s2      momentum density in 2-direction
c    s3      momentum density in 3-direction
c
c BOUNDARY VALUES USED:
c
c    Macro defined  var   ii    oi    ij    oj    ik    ok
c    -------------  ---  ----  ----  ----  ----  ----  ----
c                  mflx  is-1        js-1        ks-1  ke+1
c                    u1                          ks-2  ke+2
c                    u2                          ks-2  ke+2
c                    u3  is-1        js-1        ks-2  ke+2
c
c  OUTPUT VARIABLES:
c    s1      momentum density in 1-direction updated in the 3-direction
c    s2      momentum density in 2-direction updated in the 3-direction
c    s3      momentum density in 3-direction updated in the 3-direction
c
c  LOCAL VARIABLES:
c
c  EXTERNALS:
c
c-----------------------------------------------------------------------
c
      use real_prec
      use config
      use param
      use root
      use field
      use grid
      use scratch
c
      implicit NONE
c
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend
c
      real(rl) :: xi, q1, q2, vel, dqm,dqp
c
      real(rl) :: atwid1 (ijkn), atwid2 (ijkn), atwid3 (ijkn),
     .            atwidj1(ijkn), atwidj2(ijkn), atwidj3(ijkn)
c
      real(rl) :: sflx(ijkn,1), dq(ijkn,1)
c
      real(rl) :: mflx(in,jn,kn), s1(in,jn,kn),
     .            s2  (in,jn,kn), s3(in,jn,kn)
c.......................................................................
c
c      Compute time-centred area factors.
c
       do 10 i=ibeg,iend
         atwid1(i) = 0.5 * g2a (i) * dx1b(i) * dvl1bi(i)
         atwid2(i) = 0.5 * g2b (i) * g2b (i) * dx1a  (i) * dvl1ai(i)
         atwid3(i) = 0.5 * g31b(i) * g2b (i) * dx1a  (i) * dvl1ai(i)
10     continue
       do 20 j=jbeg,jend
         atwidj1(j) = dx2a(j) * dvl2ai(j)
         atwidj2(j) = dx2b(j) * dvl2bi(j)
         atwidj3(j) = g32b(j) * dx2a  (j) * dvl2ai(j)
20     continue
c
c.......................................................................
c
       do 3000 j=jbeg,jend
         do 2000 i=ibeg,iend
c
c------------------------------ TRANSPORT S1 ---------------------------
c
c      Interpolate "v1" at the 3-interfaces.
c
c       call x3zc3d ( v1, vel3, is, js, ie, je, iords1, istps1
c     1             , g31a, g31ai, g32b, g32bi, sflx 3, p      )
c
c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
           do 1070 k=kbeg-1,kend+1
             dqm        = (v1 (i  ,j,k  ) - v1 (i  ,j,k-1)) * dx3bi(k  )
             dqp        = (v1 (i  ,j,k+1) - v1 (i  ,j,k  )) * dx3bi(k+1)
             dq(k,1)    = max ( dqm * dqp, zro )
     1                  * sign ( one, dqm + dqp )
     2                  / max ( abs ( dqm + dqp ), tiny )
1070       continue
c
c  2.  Choose time averaged, upwinded interface value.
c
c      Construct an i-average of "v3" to be used for interpolation.
c
c      Construct the 1-momentum flux at the 3-interfaces and perform
c  1-momentum advection.  Note that the timestep "dt" is hidden in the
c  mass flux.
c
           do 1100 k=kbeg,kend+1
             vel         = 0.5 * ( v3(i-1,j,k  ) + v3(i  ,j,k  ) )
              q2         = dt * g31ai(i  ) * g32bi(j)
              xi         = ( vel           - vg3(k  ) ) * q2
              q1         = sign ( haf, xi )
              sflx (k,1) = ( 0.5 + q1 ) * ( v1 (i  ,j,k-1)
     1                   + ( dx3a(k-1) - xi ) * dq (k-1,1) )
     2                   + ( 0.5 - q1 ) * ( v1 (i  ,j,k  )
     3                   - ( dx3a(k  ) + xi ) * dq (k  ,1) )
              sflx (k,1) = ( mflx (i-1,j,k  ) + mflx (i  ,j,k  ) )
     1                   * sflx (k,1)    * atwid1(i  ) * atwidj1(j)
1100       continue
c
           do 1170 k=kbeg,kend
            if(xvgrid) then
             s1(i  ,j,k)  = ( s1(i  ,j,k) * dvl3a(k)
     1                    - sflx(k+1,1) + sflx(k,1) ) * dvl3ani(k)
            else
             s1(i  ,j,k)  = ( s1(i  ,j,k) * dvl3a(k)
     1                    - sflx(k+1,1) + sflx(k,1) ) * dvl3a i(k)
            endif
1170       continue
c
c------------------------------ TRANSPORT S2 ---------------------------
c
c      Interpolate "v2" at the 3-interfaces.
c
c       call x3zc3d ( v2, vel3, is, js, ie, je, iords2, istps2
c     1             , g31b, g31bi, g32a, g32ai, sflx 3, p      )
c
c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
           do 1270 k=kbeg-1,kend+1
             dqm       = (v2 (i  ,j,k  ) - v2 (i  ,j,k-1)) * dx3bi(k  )
             dqp       = (v2 (i  ,j,k+1) - v2 (i  ,j,k  )) * dx3bi(k+1)
             dq(k,1)   = max ( dqm * dqp, zro )
     1                 * sign ( one, dqm + dqp )
     2                 / max ( abs ( dqm + dqp ), tiny )
1270       continue
c
c  2.  Choose time averaged, upwinded interface value.
c
c      Construct a j-average of "v3" to be used for interpolation.
c
c      Construct the 2-momentum flux at the 3-interfaces and perform
c  2-momentum advection.  Note that the timestep "dt" is hidden in the
c  mass flux.
c
           do 1300 k=kbeg,kend+1
             vel         = 0.5 * ( v3(i  ,j-1,k  ) + v3(i  ,j,k  ) )
              q2         = dt * g31bi(i  ) * g32ai(j)
              xi         = ( vel           - vg3(k  ) ) * q2
              q1         = sign ( haf, xi )
              sflx (k,1) = ( 0.5 + q1 ) * ( v2 (i  ,j,k-1)
     1                   + ( dx3a(k-1) - xi ) * dq (k-1,1) )
     2                   + ( 0.5 - q1 ) * ( v2 (i  ,j,k  )
     3                   - ( dx3a(k  ) + xi ) * dq (k  ,1) )
              sflx (k,1) = ( mflx (i  ,j-1,k  ) + mflx (i  ,j  ,k  ) )
     1                   * sflx (k,1) * atwid2(i  ) * atwidj2(j)
1300       continue
c
           do 1370 k=kbeg,kend
            if(xvgrid) then
             s2(i  ,j,k) = ( s2(i  ,j,k) * dvl3a(k)
     1                   - sflx(k+1,1) + sflx(k  ,1) ) * dvl3ani(k)
            else
             s2(i  ,j,k) = ( s2(i  ,j,k) * dvl3a(k)
     1                   - sflx(k+1,1) + sflx(k  ,1) ) * dvl3a i(k)
            endif
1370       continue
c
c------------------------------ TRANSPORT S3 ---------------------------
c
c
c      Interpolate "v3" at the zone centers.
c
c       call x3fc3d ( v3, vel3, is, js, ie, je, iords3, sflx 3 )
c
c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
           do 1470 k=kbeg-1,kend+1
             dqm       = (v3 (i  ,j,k  ) - v3 (i  ,j,k-1)) * dx3ai(k-1)
             dqp       = (v3 (i  ,j,k+1) - v3 (i  ,j,k  )) * dx3ai(k  )
             dq(k,1)   = max ( dqm * dqp, zro )
     1                 * sign ( one, dqm + dqp )
     2                 / max ( abs ( dqm + dqp ), tiny )
1470       continue
c
c  2.  Choose time averaged, upwinded interface value.
c
c      Construct a k-average of "v3-vg3" to be used for interpolation.
c
c      Construct the 3-momentum flux at the 3-interfaces and perform
c  3-momentum advection.  Note that the timestep "dt" is hidden in the
c  mass flux.
c
           do 1500 k=kbeg-1,kend
             vel         = 0.5 * ( v3(i  ,j,k  ) - vg3(k  )
     1                           + v3(i  ,j,k+1) - vg3(k+1) )
              xi         = vel         * dt * g31bi(i  ) * g32bi(j)
              q1         = sign ( haf, xi )
              sflx (k,1) = ( 0.5 + q1 ) * ( v3(i  ,j,k  )
     1                   + ( dx3b(k  ) - xi ) * dq(k  ,1) )
     2                   + ( 0.5 - q1 ) * ( v3(i  ,j,k+1)
     3                   - ( dx3b(k+1) + xi ) * dq(k+1,1) )
              sflx (k,1) = ( mflx (i  ,j,k  ) + mflx (i  ,j,k+1) )
     1                   * sflx (k,1) * atwid3(i  ) * atwidj3(j)
1500       continue
c
           do 1570 k=kbeg,kend
            if(xvgrid) then
             s3(i  ,j,k) = ( s3(i  ,j,k) * dvl3b(k)
     1                   - sflx(k  ,1) + sflx(k-1,1) ) * dvl3bni(k)
            else
             s3(i  ,j,k) = ( s3(i  ,j,k) * dvl3b(k)
     1                   - sflx(k  ,1) + sflx(k-1,1) ) * dvl3b i(k)
            endif
1570       continue
2000     continue
3000   continue
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                 M O M X 3                 \\\\\\\\\\
c
c=======================================================================
c
c
