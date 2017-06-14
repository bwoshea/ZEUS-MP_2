c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 M O M X 1                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine momx1 (ibeg,iend,jbeg,jend,kbeg,kend
     &                  ,s1,s2,s3,mflx,atwid1,vtwid,sflx,dq)
c
c    dac:zeus3d.momx1 <--------------- transports momenta in 1-direction
c    from jms:zeus2d.momx1, mln:zeus04.momz                    may, 1990
c
c    written by: David Clarke
c    modified 1: November, 1992 by David Clarke; momenta are now updated
c                between and including i=is,ie, j=js,je, and k=ks,ke to
c                allow for proper treatment of periodic boundaries.
c    modified 2: Feb. 20, 1996 by Robert Fiedler; completely rewritten
c                for ZEUS-MP
c
c  PURPOSE:  Transports the three components of the momentum density in
c  the 1-direction using the consistent transport algorithm, including
c  the effects of grid compression.  The transported fluxes are thus
c  given by the mass fluxes times the time centred area of the control
c  volume faces times the interpolated velocities.  Interpolations are
c  performed in-line.
c
c  INPUT VARIABLES:
c    mflx    mass flux in 1-direction (computed in TRANX1)
c    s1      momentum density in 1-direction
c    s2      momentum density in 2-direction
c    s3      momentum density in 3-direction
c
c BOUNDARY VALUES USED:
c
c    Macro defined  var   ii    oi    ij    oj    ik    ok
c    -------------  ---  ----  ----  ----  ----  ----  ----
c                  mflx  is-1  ie+1  js-1        ks-1
c                    u1  is-2  ie+2  js-1        ks-1
c                    u2  is-2  ie+2
c                    u3  is-2  ie+2
c
c  OUTPUT VARIABLES:
c    s1      momentum density in 1-direction updated in the 1-direction
c    s2      momentum density in 2-direction updated in the 1-direction
c    s3      momentum density in 3-direction updated in the 1-direction
c
c  LOCAL VARIABLES:
c    vel     velocity used for upwinding in interpolation routine
c    vtwid   interpolated velocity
c    sflx    momentum fluxes
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
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend, jm1
c
      real(rl) :: dqm, dqp, xi, q1, vel
c
      real(rl) :: atwid1(ijkn), atwid2(ijkn), atwid3(ijkn)
      real(rl) :: vtwid (ijkn), sflx  (ijkn),
     .            dq    (ijkn)
c
      real(rl) :: mflx(in,jn,kn), s1(in,jn,kn),
     .            s2  (in,jn,kn), s3(in,jn,kn)
c
c---------------------------- TRANSPORT S1 -----------------------------
c
c      Compute time-centred area factors.
c
       do 10 i=ibeg-1,iend+1
        if(xvgrid) then
         atwid1(i) = 0.5 * g2bh(i) * g31bh(i)
         atwid2(i) = 0.5 * g2a (i) * g2ah (i) * g31ah(i)
         atwid3(i) = 0.5 * g31a(i) * g2ah (i) * g31ah(i)
        else
         atwid1(i) = 0.5 * g2b (i) * g31b (i)
         atwid2(i) = 0.5 * g2a (i) * g2a  (i) * g31a (i)
         atwid3(i) = 0.5 * g31a(i) * g2a  (i) * g31a (i)
        endif
10     continue
      do 590 k=kbeg,kend
       do 580 j=jbeg,jend
        if(ldimen .eq. 1) then
         jm1 = j
        else
         jm1 = j - 1
        endif
c
c      Interpolate "v1" at the zone centres.
c
c       call x1fc3d ( v1, vel1, js, ks, je, ke, iords1, vtwid1 )
c
c  1.  Evaluate monotonised, van Leer difference in "q" across the zone.
c
           do 70 i=ibeg-1,iend+1
             dqm        = (v1 (i  ,j,k) - v1 (i-1,j,k)) * dx1ai(i-1)
             dqp        = (v1 (i+1,j,k) - v1 (i  ,j,k)) * dx1ai(i  )
             dq(i)      = max ( dqm * dqp, zro )
     1                  * sign ( one, dqm + dqp )
     2                  / max ( abs ( dqm + dqp ), tiny )
70         continue
c
c  2.  Choose time averaged, upwinded interface value.
c
           do 100 i=ibeg-1,iend
c
c      Construct an i-average of "v1-vg1" to be used for interpolation.
c
             vel        = 0.5 * ( v1(i  ,j,k) - vg1(i  )
     1                          + v1(i+1,j,k) - vg1(i+1) )
             xi        = vel         * dt
             q1        = sign ( haf, xi )
             vtwid (i) = ( 0.5 + q1 ) * ( v1 (i  ,j,k)
     1                 + ( dx1b(i  ) - xi ) * dq  (i      ) )
     2                 + ( 0.5 - q1 ) * ( v1 (i+1,j,k)
     3                 - ( dx1b(i+1) + xi ) * dq  (i+1    ) )
c
c      Construct the 1-momentum flux at the zone centres and perform
c  1-momentum advection.  Note that the timestep "dt" is hidden in the
c  mass flux.
c
             sflx (i    ) = ( mflx (i,j,k) + mflx (i+1,j,k) )
     1                    * vtwid (i    ) * atwid1(i)
100        continue
           do 170 i=ibeg,iend
            if(xvgrid) then
             s1   (i,j,k) = ( s1(i,j,k) * dvl1b(i)
     1                    - sflx (i    ) + sflx (i-1    ) ) 
     2                    * dvl1bni(i)
            else
             s1   (i,j,k) = ( s1(i,j,k) * dvl1b(i)
     1                    - sflx (i    ) + sflx (i-1    ) ) 
     2                    * dvl1b i(i)
            endif
170        continue
C         if(ldimen .gt. 1) then
c
c---------------------------- TRANSPORT S2 -----------------------------
c
c      Interpolate "v2" at the 1-interfaces.
c
c       call x1zc3d ( v2, vel1, js, ks, je, ke, iords2, istps2
c     1             , vtwid1, p      )
c
c     1.  Evaluate monotonised, van Leer differences across the zone.
c
           do 270 i=ibeg-1,iend+1
             dqm       = (v2 (i  ,j,k) - v2 (i-1,j,k)) * dx1bi(i  )
             dqp       = (v2 (i+1,j,k) - v2 (i  ,j,k)) * dx1bi(i+1)
             dq(i)     = max ( dqm * dqp, zro )
     1                 * sign ( one, dqm + dqp )
     2                 / max ( abs ( dqm + dqp ), tiny )
270        continue
c
c     2.  Choose time averaged, upwinded interface value.
c
           do 300 i=ibeg,iend+1
             vel       = 0.5 * ( v1(i,jm1,k) + v1(i,j,k) )
             xi        = ( vel         - vg1(i) ) * dt
             q1        = sign ( haf, xi )
             vtwid (i) = ( 0.5 + q1 ) * ( v2 (i-1,j,k)
     1                 + ( dx1a(i-1) - xi ) * dq   (i-1) )
     2                 + ( 0.5 - q1 ) * ( v2 (i  ,j,k)
     3                 - ( dx1a(i  ) + xi ) * dq   (i  ) )
c
c      Construct the 2-momentum flux at the 1-interfaces and perform
c  2-momentum advection.  Note that the timestep "dt" is hidden in the
c  mass flux.
c
             sflx (i    ) = ( mflx (i,jm1,k) + mflx (i,j,k) )
     1                    * vtwid (i    ) * atwid2(i)
300        continue
           do 370 i=ibeg,iend
            if(xvgrid) then
             s2(i,j,k) = ( s2(i,j,k) * dvl1a(i)
     1                   - sflx (i+1    ) + sflx (i    ) ) * dvl1ani(i)
            else
             s2(i,j,k) = ( s2(i,j,k) * dvl1a(i)
     1                   - sflx (i+1    ) + sflx (i    ) ) * dvl1a i(i)
            endif
370        continue
C         endif ! ldimen
c
CJH         if(ldimen .eq. 3) then
c
c---------------------------- TRANSPORT S3 -----------------------------
c
c      Interpolate "v3" at the 1-interfaces.
c
c       call x1zc3d ( v3, vel1, js, ks, je, ke, iords3, istps3
c     1             , vtwid1, p      )
c
c     1.  Evaluate monotonised, van Leer differences across the zone.
c
           do 470 i=ibeg-1,iend+1
             dqm       = (v3 (i  ,j,k) - v3 (i-1,j,k)) * dx1bi(i  )
             dqp       = (v3 (i+1,j,k) - v3 (i  ,j,k)) * dx1bi(i+1)
             dq(i)     = max ( dqm * dqp, zro )
     1                 * sign ( one, dqm + dqp )
     2                 / max ( abs ( dqm + dqp ), tiny )
470        continue
c
c     2.  Choose time averaged, upwinded interface value.
c
           do 500 i=ibeg,iend+1
c
c      Construct a k-average of "v1" to be used for interpolation.
c
            if(ldimen .eq. 3) then
             vel       = 0.5 * ( v1(i,j,k-1) + v1(i,j,k) )
            else
             vel       = v1(i,j,k)
            endif
             xi        = ( vel         - vg1(i) ) * dt
             q1        = sign ( haf, xi )
             vtwid (i) = ( 0.5 + q1 ) * ( v3 (i-1,j,k)
     1                 + ( dx1a(i-1) - xi ) * dq   (i-1) )
     2                 + ( 0.5 - q1 ) * ( v3 (i  ,j,k)
     3                 - ( dx1a(i  ) + xi ) * dq   (i  ) )
c
c      Construct the 3-momentum flux at the 1-interfaces and perform
c  3-momentum advection.  Note that the timestep "dt" is hidden in the
c  mass flux.
c
            if(ldimen .eq. 3) then
             sflx (i    ) = ( mflx (i,j,k-1) + mflx (i,j,k) )
     1                    * vtwid (i    ) * g32b(j) * atwid3(i)
            else
             sflx (i    ) = 2.0*( mflx (i,j,k) )
     1                    * vtwid (i    ) * g32b(j) * atwid3(i)
            endif
500        continue
           do 570 i=ibeg,iend
            if(xvgrid) then
             s3(i,j,k) = ( s3(i,j,k) * dvl1a(i)
     1                   - sflx (i+1    ) + sflx (i    ) ) * dvl1ani(i)
            else
             s3(i,j,k) = ( s3(i,j,k) * dvl1a(i)
     1                   - sflx (i+1    ) + sflx (i    ) ) * dvl1a i(i)
            endif
570        continue
C         else ! ldimen
C          do i = ibeg, iend
C           s3(i,j,k) = 0.0D0
C          enddo
C         endif ! ldimen
580    continue
590   continue
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                 M O M X 1                 \\\\\\\\\\
c
c=======================================================================
c
c
