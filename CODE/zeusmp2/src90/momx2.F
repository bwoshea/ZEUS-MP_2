c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 M O M X 2                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine momx2 (ibeg,iend,jbeg,jend,kbeg,kend
     &                  ,s1,s2,s3,mflx,atwid1,atwid2,atwid3,atwidj)
c
c    dac:zeus3d.momx2 <--------------- transports momenta in 2-direction
c    from jms:zeus2d.momx2, mln:zeus04.momr                    may, 1990
c
c    written by: David Clarke
c    modified 1: November, 1992 by David Clarke; momenta are now updated
c                between and including i=is,ie, j=js,je, and k=ks,ke to
c                allow for proper treatment of periodic boundaries.
c    modified 2: Feb. 20, 1996 by Robert Fiedler; completely rewritten
c                for ZEUS-MP.
c    modified 3: by John Hayes; rewritten for F90
c    modified 4: 03/29/2006 by John Hayes; corrected typos in "atwid2"
c                indexing in s2 update
c    modified 5: 06/02/2006 by John Hayes; rerolled I loops and corrected
c                bug in "sflx" formula which manifested itself on reflecting
c                2-boundaries in spherical geometry.
c                
c
c  PURPOSE:  Transports the three components of the momentum density in
c  the 2-direction using the consistent transport algorithm, including
c  the effects of grid compression.  The transported fluxes are thus
c  given by the mass fluxes times the time centred area of the control
c  volume faces times the interpolated velocities.  Interpolations are
c  performed in-line.
c
c  INPUT VARIABLES:
c    mflx    mass flux in 2-direction; computed by TRANX2
c    s1      momentum density in 1-direction
c    s2      momentum density in 2-direction
c    s3      momentum density in 3-direction
c
c BOUNDARY VALUES USED:
c
c    Macro defined  var   ii    oi    ij    oj    ik    ok
c    -------------  ---  ----  ----  ----  ----  ----  ----
c                  mflx  is-1        js-1  je+1  ks-1
c                    u1              js-2  je+2
c                    u2  is-1        js-2  je+2  ks-1
c                    u3              js-2  je+2
c
c  OUTPUT VARIABLES:
c    s1      momentum density in 1-direction updated in the 2-direction
c    s2      momentum density in 2-direction updated in the 2-direction
c    s3      momentum density in 3-direction updated in the 2-direction
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
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend, km1
c
      real(rl) :: dqm, dqp, vel, xi, q1
      real(rl) :: atwid1(ijkn), atwid2(ijkn), atwid3(ijkn),
     .            atwidj(ijkn)
      real(rl) :: sflx(ijkn), dq(ijkn)
      real(rl) :: mflx(in,jn,kn), s1(in,jn,kn),
     .            s2  (in,jn,kn), s3(in,jn,kn)
      real(rl) :: qty1, qty2, qty3, qty4
c
      integer  iriter
c.......................................................................
c
c      Compute time-centred area factors.
c
       do 10 i=ibeg,iend
         atwid1(i) = 0.5 * g31a(i)           * dx1b(i) * dvl1bi(i)
         atwid2(i) = 0.5 * g2b (i) * g31b(i) * dx1a(i) * dvl1ai(i)
         atwid3(i) = 0.5 * g31b(i) * g31b(i) * dx1a(i) * dvl1ai(i)
10     continue
c
       do 20 j=jbeg,jend+1
        if(xvgrid) then
         atwidj(j) = g32a(j) * g32ah(j)
        else
         atwidj(j) = g32a(j) * g32a (j)
        endif
20     continue
c
c.......................................................................
c
       do 3000 k=kbeg,kend

c Do the iterations which will not be done in unrolled loops first.
c
         iriter = mod(iend-ibeg+1,4)
         do 2000 i=ibeg,iend
c
c---------------------------- TRANSPORT S1 -----------------------------
c
           do 1070 j=jbeg-1,jend+1
c
c      Interpolate "v1" at the 2-interfaces.
c
c       call x2zc3d ( v1, vel2, ks, is, ke, ie, iords1, istps1
c     1             , g2a, g2ai, vtwid2, p      )
c
c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
             dqm      = (v1 (i  ,j  ,k) - v1 (i  ,j-1,k)) * dx2bi  (j  )
             dqp      = (v1 (i  ,j+1,k) - v1 (i  ,j  ,k)) * dx2bi  (j+1)
             dq(j)  = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )
c
1070       continue
c
c  2.  Choose time averaged, upwinded interface value.
c
           do 1100 j=jbeg,jend+1
c
c      Construct an i-average of "v2" to be used for interpolation.
c
c      Construct the 1-momentum flux at the 2-interfaces and perform
c  1-momentum advection.  Note that the timestep "dt" is hidden in the
c  mass flux.
c
              vel        = 0.5 * ( v2(i-1,j  ,k) + v2(i  ,j  ,k) )
              xi         = ( vel         - vg2(j  ) ) * dt * g2ai(i  )
              q1         = sign ( haf, xi )
              sflx (j) = ( 0.5 + q1 ) * ( v1 (i  ,j-1,k)
     1                   + ( dx2a(j-1) - xi ) * dq (j-1) )
     2                   + ( 0.5 - q1 ) * ( v1 (i  ,j  ,k)
     3                   - ( dx2a(j  ) + xi ) * dq (j  ) )
             if(xvgrid) then
              sflx (j) = ( mflx (i-1,j  ,k) + mflx (i  ,j  ,k) )
     1                   * sflx (j)   * g32ah(j  ) * atwid1(i  )
             else
              sflx (j) = ( mflx (i-1,j  ,k) + mflx (i  ,j  ,k) )
     1                   * sflx (j)   * g32a (j  ) * atwid1(i  )
             endif
c
1100       continue
c
           do 1170 j=jbeg,jend
            if(xvgrid) then
             s1(i  ,j,k)  = ( s1(i  ,j,k) * dvl2a(j)
     1                    - sflx (j+1)   + sflx (j)   ) * dvl2ani(j)
            else
             s1(i  ,j,k)  = ( s1(i  ,j,k) * dvl2a(j)
     1                    - sflx (j+1)   + sflx (j)   ) * dvl2a i(j)
            endif
c
1170       continue
          if(ldimen .gt. 1) then
c
c---------------------------- TRANSPORT S2 -----------------------------
c
c      Interpolate "v2" to the zone centers.
c
c       call x2fc3d ( v2, vel2, ks, is, ke, ie, iords2, sflx2 )
c
c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
           do 1270 j=jbeg-1,jend+1
             dqm        = (v2 (i  ,j  ,k) - v2 (i  ,j-1,k)) * dx2ai(j-1)
             dqp        = (v2 (i  ,j+1,k) - v2 (i  ,j  ,k)) * dx2ai(j  )
             dq(j)    = max ( dqm * dqp, zro )
     1                  * sign ( one, dqm + dqp )
     2                  / max ( abs ( dqm + dqp ), tiny )
c
1270       continue
c
c  2.  Choose time averaged, upwinded interface value.
c
c      Construct a j-average of "v2-vg2" to be used for interpolation.
c
c      Construct the 2-momentum flux at the zone centres and perform
c  2-momentum advection.  Note that the timestep "dt" is hidden in the
c  mass flux.
c
           do 1300 j=jbeg-1,jend
              vel        = 0.5 * ( v2(i  ,j  ,k) - vg2(j  )
     1                           + v2(i  ,j+1,k) - vg2(j+1) )
              xi         = vel         * dt * g2bi(i  )
              q1         = sign ( haf, xi )
              sflx (j) = ( 0.5 + q1 ) * ( v2 (i  ,j  ,k)
     1                   + ( dx2b(j  ) - xi ) * dq (j  ) )
     2                   + ( 0.5 - q1 ) * ( v2 (i  ,j+1,k) 
     3                   - ( dx2b(j+1) + xi ) * dq (j+1) )
             if(i .eq. 11 .and. j .eq. 2) then
              qty1 = sflx(j)
             endif
CJH
CJH  taking abs value of g32b (or g32bh) to account for 
CJH  coordinate reflection in RTP geometry
CJH
             if(xvgrid) then
              sflx (j) = ( mflx (i  ,j  ,k) + mflx (i  ,j+1,k) )
     1                   * sflx (j) * abs(g32bh(j))* atwid2(i)
             else
              sflx (j) = ( mflx (i  ,j  ,k) + mflx (i  ,j+1,k) )
     1                   * sflx (j) * abs(g32b (j))* atwid2(i)
             endif
c
1300       continue
c
           do 1370 j=jbeg,jend
            if(xvgrid) then
             s2(i  ,j,k)  = ( s2(i  ,j,k) * dvl2b(j)
     1                    - sflx (j)   + sflx (j-1) )
     2                    * dvl2bni(j)
            else
             s2(i  ,j,k)  = ( s2(i  ,j,k) * dvl2b(j)
     1                    - sflx (j)   + sflx (j-1) )
     2                    * dvl2b i(j)
            endif
c
1370       continue
         endif ! ldimen
c
c---------------------------- TRANSPORT S3 -----------------------------
C         if(ldimen .eq. 3) then
c
c      Interpolate "v3" at the 2-interfaces.
c
c       call x2zc3d ( v3, vel2, ks, is, ke, ie, iords3, istps3
c     1             , g2b, g2bi, sflx2, p      )
c
c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
           do 1470 j=jbeg-1,jend+1
            dqm         = (v3 (i  ,j  ,k) - v3 (i  ,j-1,k)) * dx2bi(j  )
            dqp         = (v3 (i  ,j+1,k) - v3 (i  ,j  ,k)) * dx2bi(j+1)
            dq(j)     = max ( dqm * dqp, zro )
     1                  * sign ( one, dqm + dqp )
     2                  / max ( abs ( dqm + dqp ), tiny )
c
1470       continue
c
c  2.  Choose time averaged, upwinded interface value.
c
           do 1500 j=jbeg,jend+1
c
c      Construct a  k-average of "v2" to be used for interpolation.
c
c      Construct the 3-momentum flux at the 2-interfaces and perform
c  3-momentum advection.  Note that the timestep "dt" is hidden in the
c  mass flux.
c
             if(ldimen .eq. 3) then
              vel        = 0.5 * ( v2(i  ,j  ,k-1) + v2(i  ,j  ,k) )
             else
              vel        = v2(i  ,j  ,k)
             endif
              xi         = ( vel         - vg2(j  ) ) * dt * g2bi(i  )
              q1         = sign ( haf, xi )
              sflx (j) = ( 0.5 + q1 ) * ( v3 (i  ,j-1,k)
     1                   + ( dx2a(j-1) - xi ) * dq (j-1) )
     2                   + ( 0.5 - q1 ) * ( v3 (i  ,j  ,k)
     3                   - ( dx2a(j  ) + xi ) * dq (j  ) )
             if(ldimen .eq. 3) then
              sflx (j) = ( mflx (i  ,j  ,k-1) + mflx (i  ,j  ,k  ) )
     1                   * sflx (j  ) * atwid3(i  ) * atwidj(j  )
             else
              sflx (j) = 2.0*( mflx (i  ,j  ,k  ) )
     1                   * sflx (j  ) * atwid3(i  ) * atwidj(j  )
             endif
c
1500       continue
c
           do 1570 j=jbeg,jend
            if(xvgrid) then
             s3(i  ,j,k)    = ( s3(i  ,j,k) * dvl2a(j)
     1                      - sflx (j+1)     + sflx (j  )   )
     2                    *  dvl2ani(j)
            else
             s3(i  ,j,k)    = ( s3(i  ,j,k) * dvl2a(j)
     1                      - sflx (j+1)     + sflx (j  )   )
     2                    *  dvl2a i(j)
            endif
c
1570       continue
2000     continue
c
3000   continue
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                 M O M X 2                 \\\\\\\\\\
c
c=======================================================================
c
c
