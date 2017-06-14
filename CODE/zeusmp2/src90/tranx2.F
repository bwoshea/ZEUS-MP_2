c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                T R A N X 2                \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine tranx2 (ibeg,iend,jbeg,jend,kbeg,kend
     &                   ,dlo,den
     &                   ,eod,edn
     &                   ,ero,ern
     &                   ,abo,abn
     &                   ,mflx,atwid)
c
c    dac:zeus3d.tranx2 <----- transports zone-centred variables along x2
c    from jms:zeus2d.tranx2, mln:zeus04.tranz                  may, 1990
c
c    written by: David Clarke
c    modified 1: June 1992, by David Clarke; added the total energy
c                option originally designed by Byung-IL Jun.
c    modified 2: Feb. 20, 1996 by Robert Fiedler; completely rewritten
c                for ZEUS-MP.
c    modified 3: Aug. 8, 1996 by Robert Fiedler; unrolled i-loops.
c    modified 4: Dec. 23, 1996 by Robert Fiedler; added radiation.
c    modified 5: Sep. 8, 2003 by John Hayes; re-rolled i-loops.
c
c  PURPOSE:  Transports all zone centred variables in the 2-direction
c  only.  Currently transported are:
c
c                      mass   density
c                      energy density
c
c  The consistent advection algorithm, in which mass fluxes are used to
c  construct the fluxes of all variables across the interfaces, is used
c  for all hydrodynamical variables.  Thus, the mass fluxes are passed
c  to MOMX2 on order to transport the momenta as well.  The magnetic
c  field components are treated separately from the mass fluxes in CT.
c  Interpolations are done in-line.
c
c  INPUT VARIABLES: 
c    ibeg,iend,jbeg,jend,kbeg,kend  index ranges to cover.
c    dlo      mass            density at previous substep.
c    eod      specific energy density at previous substep; equals
c             (e+p)/d  if TOTAL_ENERGY is defined.
c    ero      specific radiation energy density at previous substep.
c
c BOUNDARY VALUES USED:
c
c    Macro defined  var   ii    oi    ij    oj    ik    ok
c    -------------  ---  ----  ----  ----  ----  ----  ----
c                    d   is-1        js-3  je+2  ks-1
c                   e/d              js-2  je+2
c    TOTAL_ENERGY    u1        ie+1                         
c                    u2  is-1        js-1  je+1  ks-1
c    TOTAL_ENERGY    u3                                ke+1
c
c  OUTPUT VARIABLES:
c    den      updated mass            density.
c    edn      updated specific energy density.
c    ern      updated specific radiation energy density.
c    mflx     mass flux (in the 2-direction)
c
c  LOCAL VARIABLES:
c    atwid    effective cross sectional area of the 1-interfaces
c    etwid    interpolated specific energy densities (e/d) at all
c             2-interfaces
c    rtwid    interpolated specific radiation energy densities at all
c             2-interfaces
c    eflx     energy density flux across all 2-interfaces  (reuse etwid)
c    dtwid    interpolated mass densities at all 2-interfaces
c    dflx     mass density flux across all 2-interfaces    (reuse dtwid)
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
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend, kp1, l,
     .            kstart, n, m
c
      real(rl) :: dqm,dqp, xi,q1, atwid(ijkn)
c
      real(rl) :: mflux(ijkn),
     .            dtwid(ijkn), dd  (ijkn),
     .            etwid(ijkn), deod(ijkn)
c
      real(rl) :: dlo(in,jn,kn), den(in,jn,kn), mflx(in,jn,kn),
     .            eod(in,jn,kn), edn(in,jn,kn), 
     .            ero(in,jn,kn),
     .            ern(in,jn,kn),
     .            abo(in,jn,kn,nspec), abn(in,jn,kn,nspec)
c
      real(rl) :: rtwid(ijkn), dero(ijkn), 
     .            dxo(ijkn,nspec),
     .            xtwid(ijkn,nspec)
c-----------------------------------------------------------------------
c
c Compute time-centered area factors.
c
      do 10 i=ibeg-1,iend+1
        atwid (i)         = g31b(i) * dx1a(i) * dvl1ai(i)
10    continue
c
c Transport all zone-centered quantities in the 2 direction.
c
      if(ldimen .ne. 3) then
       kstart = ks
      else
       kstart = kbeg-1
      endif
      do 2100 k=kstart,kend
       if(ldimen .eq. 3) then
        kp1 = k+1
       else
        kp1 = ks
       endif
       if ( (k.ge.kbeg) .or. (k.eq.ks-1) ) then
       do 1090 i=ibeg-1,iend

c
c   Interpolate to obtain zone-centered quantities at zone faces.
c
c     1.  Evaluate monotonised, van Leer differences across the zone.
c
        if (jbeg .eq. js) then  !  Need d(js-3) from neighbor.
         j         = js - 2
         dqm       = (dlo(i,j  ,k) - dijb (i,k,3)) * dx2bi(j  )
         dqp       = (dlo(i,j+1,k) - dlo(i,j  ,k)) * dx2bi(j+1)
         dd  (j) = max ( dqm * dqp, zro )
     1             * sign ( one, dqm + dqp )
     2             / max ( abs ( dqm + dqp ), tiny )
         if(xiso .eqv. .false.) then
          deod(j) = zro  ! Not valid, but we don't use it.
         endif
         if(lrad .ne. 0) then
           dero(j) = zro
         endif
        endif
        do 1030 j=max(jbeg-2,js-1),jend+1
         dqm       = (dlo(i,j  ,k) - dlo(i,j-1,k)) * dx2bi(j  )
         dqp       = (dlo(i,j+1,k) - dlo(i,j  ,k)) * dx2bi(j+1)
         dd  (j) = max ( dqm * dqp, zro )
     1             * sign ( one, dqm + dqp )
     2             / max ( abs ( dqm + dqp ), tiny )
         if(nspec .gt. 1) then
          do n = 1, nspec
           dqm        = (abo(i,j  ,k,n)-abo(i,j-1,k,n))*dx2bi(j  )
           dqp        = (abo(i,j+1,k,n)-abo(i,j  ,k,n))*dx2bi(j+1)
           dxo(j,n) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )
          enddo ! n
         endif
         if(xiso .eqv. .false.) then
          dqm       = (eod(i,j  ,k) - eod(i,j-1,k)) * dx2bi(j  )
          dqp       = (eod(i,j+1,k) - eod(i,j  ,k)) * dx2bi(j+1)
          deod(j) = max ( dqm * dqp, zro )
     1              * sign ( one, dqm + dqp )
     2              / max ( abs ( dqm + dqp ), tiny )
         endif ! xiso
         if(lrad .ne. 0) then
          dqm       = (ero(i,j  ,k)-ero(i,j-1,k))*dx2bi(j  )
          dqp       = (ero(i,j+1,k)-ero(i,j  ,k))*dx2bi(j+1)
          dero(j) = max ( dqm * dqp, zro )
     1              * sign ( one, dqm + dqp )
     2              / max ( abs ( dqm + dqp ), tiny )
         endif
1030    continue
c
c     2.  Choose time averaged, upwinded interface value.
c
c  For the purposes of consistent advection, construct the mass
c  flux across each 1-interface.  The mass flux will be used to create
c  the fluxes of all variables, including the momenta which are updated
c  in MOMX1.
c
        do 1040 j=jbeg-1,jend+1
         xi          = ( v2  (i,j  ,k) - vg2(j  ) ) * dt * g2bi(i)
         q1          = sign ( haf, xi )
         dtwid (j) = ( 0.5 + q1 ) * ( dlo(i,j-1,k)
     1               + ( dx2a(j-1) - xi ) * dd   (j-1) )
     2               + ( 0.5 - q1 ) * ( dlo(i,j  ,k)
     3               - ( dx2a(j  ) + xi ) * dd   (j  ) )
c
         if(nspec .gt. 1) then
          do n = 1, nspec
           xtwid(j,n) = (0.5 + q1)*(abo(i,j-1,k,n)
     1                  + (dx2a(j-1) - xi ) * dxo(j-1,n) )
     2                  + (0.5 - q1)*(abo(i,j  ,k,n)
     3                  - (dx2a(j  ) + xi ) * dxo(j,n  ) )
          enddo
         endif ! nspec
c
         mflux (j) = dtwid (j) * ( v2(i,j  ,k) - vg2(j  ) ) * dt
c
         if(xvgrid) then
          dtwid (j) = mflux (j) * atwid (i) * g32ah(j  )
         else
          dtwid (j) = mflux (j) * atwid (i) * g32a (j  )
         endif
c
         if(nspec .gt. 1) then
          do n = 1, nspec
           xtwid(j,n) = xtwid(j,n)*dtwid(j)
          enddo
         endif
c
         if(xiso .eqv. .false.) then
          etwid (j)= ( 0.5 + q1 ) * ( eod(i,j-1,k)
     1               + ( dx2a(j-1) - xi ) * deod (j-1) )
     2               + ( 0.5 - q1 ) * ( eod(i,j  ,k)
     3               - ( dx2a(j  ) + xi ) * deod (j  ) )
c
          etwid (j)= dtwid (j) * etwid (j)
         endif ! xiso
         if(lrad .ne. 0) then
          rtwid (j)= ( 0.5 + q1 ) * ( ero(i,j-1,k)
     1               + ( dx2a(j-1) - xi ) * dero (j-1) )
     2               + ( 0.5 - q1 ) * ( ero(i,j,k)
     3               - ( dx2a(j  ) + xi ) * dero (j) )
          rtwid (j)= dtwid(j)*rtwid(j)
         endif
1040    continue
c
c  Save the mass flux outside (ibeg:iend,jbeg:jend,kbeg:kend)
c  only for zones next to the inner orders.
c
        if( (i.eq.is-1 .and. k.ge.kbeg) .or.
     &      (i.ge.ibeg .and. k.eq.ks-1) .or.
     &      (i.eq.is-1 .and. k.eq.ks-1) ) then
         if (jbeg .eq. js) mflx(i,js-1,k) = mflux (js-1)
         do 1060 j=jbeg,jend
          mflx(i,j,k) = mflux(j)
1060     continue
         if (jend .eq. je) mflx(i,je+1,k) = mflux (je+1)
        endif
c
c  Perform mass density and energy density advection.  Note that
c  the timestep "dt" is hidden the fluxes "dflx" and "eflx".
c
        if (k.ge.kbeg .and. i.ge.ibeg) then
         if (jbeg .eq. js) mflx(i,js-1,k) = mflux (js-1)
         if (jend .eq. je) mflx(i,je+1,k) = mflux (je+1)
         do 1080 j=jbeg,jend
          mflx(i,j,k) = mflux (j)
c
          if(xvgrid) then
c
            den(i,j,k) = ( dlo(i,j,k) * dvl2a(j)
     1                 - dtwid (j+1)  + dtwid (j)  ) * dvl2ani(j)
c
          else
c
            den(i,j,k) = ( dlo(i,j,k) * dvl2a(j)
     1                 - dtwid (j+1)  + dtwid (j)  ) * dvl2a i(j)
c
          endif ! xvgrid
C
          if(nspec .gt. 1) then
           do n = 1, nspec
            if(xvgrid) then
c
              abn(i,j,k,n) = (abo(i,j,k,n)*dlo(i,j,k)*dvl2a(j)
     .                     - xtwid(j+1,n)+xtwid(j,n))*dvl2ani(j)
     .                     / den(i,j,k)
c
            else ! xvgrid
c
              abn(i,j,k,n) = (abo(i,j,k,n)*dlo(i,j,k)*dvl2a(j)
     .                     - xtwid(j+1,n)+xtwid(j,n))*dvl2a i(j)
     .                     / den(i,j,k)
            endif ! xvgrid
           enddo ! n
          endif ! nspec
C
          if(xiso .eqv. .false.) then
           if(xvgrid) then
c
             e  (i,j,k) = ( e  (i,j,k) * dvl2a(j)
     .                  - etwid (j+1)  + etwid (j)  ) * dvl2ani(j)
c
           else
c
             e  (i,j,k) = ( e  (i,j,k) * dvl2a(j)
     .                  - etwid (j+1)  + etwid (j)  ) * dvl2a i(j)
           endif
c
c Compute e/d for the next substep.
c
           if(xtotnrg .eqv. .false.) then
            edn(i,j,k) = e(i,j,k) / den(i,j,k)
           else ! xtotnrg
            edn(i,j,k) = gamma * e(i,j,k) / den(i,j,k)
     .                 - gamm1 * ( ( v1(i,j,k) + v1(i+1,j  ,k  ) )**2
     .                           + ( v2(i,j,k) + v2(i  ,j+1,k  ) )**2
     .                           + ( v3(i,j,k) + v3(i  ,j  ,kp1) )**2
     .                           ) * 0.125
            if(xgrav .or. xgrvfft) then
             edn(i,j,k) = edn(i,j,k)
     .                  + gamm1 * gp(i,j,k)
            endif
           endif ! xtotnrg
          endif ! xiso
C
          if(lrad .ne. 0) then
           if(xvgrid) then
            ern(i  ,j,k) = ( ero(i  ,j,k) * dlo(i  ,j,k) * dvl2a(j)
     1                   - rtwid (j+1) + rtwid (j) ) * dvl2ani(j)
           else
            ern(i  ,j,k) = ( ero(i  ,j,k) * dlo(i  ,j,k) * dvl2a(j)
     1                   - rtwid (j+1) + rtwid (j) ) * dvl2a i(j)
           endif
c
c Work with er/d for the next sweep.
c
           ern(i  ,j,k) = ern(i  ,j,k) / den(i  ,j,k)
          endif
1080     continue
        endif
1090   continue
c
       endif  !  k=ks-1 || k>=kbeg
2100   continue
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                T R A N X 2                \\\\\\\\\\
c
c=======================================================================
c
c
