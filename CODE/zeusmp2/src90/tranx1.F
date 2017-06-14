c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                T R A N X 1                \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine tranx1 (ibeg,iend,jbeg,jend,kbeg,kend
     &                   ,dlo,den
     &                   ,eod,edn
     &                   ,ero,ern
     &                   ,abo,abn
     &                   ,mflx,atwid,dtwid,etwid,mflux,dd,deod)
c
c    dac:zeus3d.tranx1 <----- transports zone-centred variables along x1
c    from jms:zeus2d.tranx1, mln:zeus04.tranz                  may, 1990
c
c    written by: David Clarke
c    modified 1: June 1992, by David Clarke; added the total energy
c                option originally designed by Byung-IL Jun.
c    modified 2: Feb. 20, 1996 by Robert Fiedler; completely rewritten
c                for ZEUS-MP.
c    modified 3: Dec. 19, 1996 by Robert Fiedler; added radiation
c    modified 4: October 2005 by John Hayes; corrected transposed subscripts
c                in "dxo" array
c
c  PURPOSE:  Transports all zone centred variables in the 1-direction
c  only.  Currently transported are:
c
c                      mass   density
c                      energy density
c
c  The consistent advection algorithm, in which mass fluxes are used to
c  construct the fluxes of all variables across the interfaces, is used
c  for all hydrodynamical variables.  Thus, the mass fluxes are passed
c  to MOMX1 on order to transport the momenta as well.  The magnetic
c  field components are treated separately from the mass fluxes in CT.
c  Interpolations are done in-line.
c
c  INPUT VARIABLES: 
c    ibeg,iend,jbeg,jend,kbeg,kend  index ranges to cover.
c    dlo      mass            density at previous substep.
c    eod      specific energy density at previous substep; equals
c             (e+p)/d  if TOTAL_ENERGY is defined.
c
c BOUNDARY VALUES USED:
c
c    Macro defined  var   ii    oi    ij    oj    ik    ok
c    -------------  ---  ----  ----  ----  ----  ----  ----
c                    d   is-3  ie+2  js-1        ks-1
c                   e/d  is-2  ie+2
c                    u1  is-1  ie+1  js-1        ks-1
c    TOTAL_ENERGY    u2                    je+1
c    TOTAL_ENERGY    u3                                ke+1
c
c  OUTPUT VARIABLES:
c    den      updated mass            density.
c    edn      updated specific energy density.
c    mflx     mass flux (in the 1-direction)
c
c  LOCAL VARIABLES:
c    atwid    effective cross sectional area of the 1-interfaces
c    etwid    interpolated specific energy densities (e/d) at all
c             1-interfaces
c    eflx     energy density flux across all 1-interfaces  (reuse etwid)
c    dtwid    interpolated mass densities at all 1-interfaces
c    dflx     mass density flux across all 1-interfaces    (reuse dtwid)
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
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend, n, kp1
c
      real(rl) :: dqm, dqp, xi,q1
      real(rl) :: atwid(ijkn), mflux(ijkn),
     .            dtwid(ijkn), dd   (ijkn),
     .            etwid(ijkn), deod (ijkn),
     .            rtwid(ijkn), dero (ijkn),
     .            xtwid(ijkn,nspec), dxo (ijkn,nspec)
c
      real(rl) :: mflx(in,jn,kn)
      real(rl) :: dlo(in,jn,kn), den(in,jn,kn),
     .            eod(in,jn,kn), edn(in,jn,kn),
     .            ero(in,jn,kn), ern(in,jn,kn),
     .            abo(in,jn,kn,nspec), abn(in,jn,kn,nspec)
c-----------------------------------------------------------------------
c
c Compute time-centered area factors.
c
      do 10 i=ibeg-1,iend+1
       if(xvgrid) then
        atwid (i)  =       g2ah(i) * g31ah(i)
       else
        atwid (i)  =       g2a (i) * g31a (i)
       endif
10    continue
c
c Transport all zone-centered quantities in the 1 direction.
c Note that transporting v1 in MOMX1 will require the mass flux at 
c x1a(is-1).  To get it from the field variables, we need d at is-3.
c We also need mflx at js-1 and ks-1 for i=is,ie+1 for v2 and v3.
c Extend loop indices to compute mass fluxes beyond inner borders.
c Be careful to assign values to mflx, (and den, edn) only
c within the range (ibeg:iend,jbeg:jend,kbeg:kend) when these
c indices are not on the borders, so that they can't get overwritten
c when this routine is called with various index ranges.
c
      do 100 k=kbeg-1,kend
       do 90 j=jbeg-1,jend
c
c   Interpolate to obtain zone-centered quantities at zone faces.
c
c     1.  Evaluate monotonised, van Leer differences across the zone.
c
        if (ibeg .eq. is) then  !  Need d(is-3) from neighbor.
         i         = is - 2
         dqm       = (dlo(i  ,j,k) - diib (j,k,3)) * dx1bi(i  )
         dqp       = (dlo(i+1,j,k) - dlo(i  ,j,k)) * dx1bi(i+1)
         dd(i)     = max ( dqm * dqp, zro )
     1             * sign ( one, dqm + dqp )
     2             / max ( abs ( dqm + dqp ), tiny )
         if(xiso .eqv. .false.) then
          deod  (i) = zro  ! Not valid, but we don't use it.
         endif
         if(lrad .ne.0) dero  (i) = zro
        endif ! ibeg
        do 30 i=max(ibeg-2,is-1),iend+1
         dqm       = (dlo(i  ,j,k) - dlo(i-1,j,k)) * dx1bi(i  )
         dqp       = (dlo(i+1,j,k) - dlo(i  ,j,k)) * dx1bi(i+1)
         dd(i)     = max ( dqm * dqp, zro )
     1             * sign ( one, dqm + dqp )
     2             / max ( abs ( dqm + dqp ), tiny )
         if(nspec .gt. 1) then
          do n = 1, nspec
           dqm       = (abo(i  ,j,k,n)-abo(i-1,j,k,n))*dx1bi(i  )
           dqp       = (abo(i+1,j,k,n)-abo(i  ,j,k,n))*dx1bi(i+1)
           dxo(i,n) = max ( dqm * dqp, zro )
     1               * sign ( one, dqm + dqp )
     2               / max ( abs ( dqm + dqp ), tiny )
          enddo
         endif ! nspec
         if(xiso .eqv. .false.) then
          dqm       = (eod(i  ,j,k) - eod(i-1,j,k)) * dx1bi(i  )
          dqp       = (eod(i+1,j,k) - eod(i  ,j,k)) * dx1bi(i+1)
          deod(i)   = max ( dqm * dqp, zro )
     1              * sign ( one, dqm + dqp )
     2              / max ( abs ( dqm + dqp ), tiny )
         endif ! xiso
         if(lrad .ne. 0) then
          dqm       = (ero(i  ,j,k) - ero(i-1,j,k)) * dx1bi(i  )
          dqp       = (ero(i+1,j,k) - ero(i  ,j,k)) * dx1bi(i+1)
          dero(i)   = max ( dqm * dqp, zro )
     1              * sign ( one, dqm + dqp )
     2              / max ( abs ( dqm + dqp ), tiny )
         endif ! lrad
30      continue
c
c     2.  Choose time averaged, upwinded interface values.
c
        do 40 i=ibeg-1,iend+1
         xi        = ( v1  (i,j,k) - vg1(i) ) * dt
         q1        = sign ( haf, xi )
         dtwid (i) = ( 0.5 + q1 ) * ( dlo(i-1,j,k)
     1              + ( dx1a(i-1) - xi ) * dd   (i-1) )
     2              + ( 0.5 - q1 ) * ( dlo(i  ,j,k)
     3              - ( dx1a(i  ) + xi ) * dd   (i  ) )
         if(nspec .gt. 1) then
          do n = 1, nspec
           xtwid(i,n) = ( 0.5 + q1 )*(abo(i-1,j,k,n)
     1                + ( dx1a(i-1) - xi )*dxo(i-1,n) )
     2                + ( 0.5 - q1 )*(abo(i  ,j,k,n)
     3                - ( dx1a(i  ) + xi )*dxo(i  ,n) )
          enddo
         endif
         if(xiso .eqv. .false.) then
          etwid (i) = ( 0.5 + q1 ) * ( eod(i-1,j,k)
     1              + ( dx1a(i-1) - xi ) * deod (i-1) )
     2              + ( 0.5 - q1 ) * ( eod(i  ,j,k)
     3              - ( dx1a(i  ) + xi ) * deod (i  ) )
         endif ! xiso
         if(lrad .ne. 0) then
          rtwid (i) = ( 0.5 + q1 ) * ( ero(i-1,j,k)
     1              + ( dx1a(i-1) - xi ) * dero (i-1) )
     2              + ( 0.5 - q1 ) * ( ero(i  ,j,k)
     3              - ( dx1a(i  ) + xi ) * dero (i  ) )
         endif ! lrad
40      continue
c
c  For the purposes of consistent advection, construct the mass
c  flux across each 1-interface.  The mass flux will be used to create
c  the fluxes of all variables, including the momenta which are updated
c  in MOMX1.
c
        do 50 i=ibeg-1,iend+1
         mflux (i    ) = dtwid (i    ) * ( v1(i,j,k) - vg1(i) ) * dt
         dtwid (i    ) = mflux (i    ) * atwid (i    )
         if(nspec .gt. 1)then
          do n = 1, nspec
           xtwid(i,n) = xtwid(i,n) * dtwid(i)
          enddo
         endif
         if(xiso .eqv. .false.) then
          etwid (i    ) = dtwid (i    ) * etwid (i    )
         endif
         if(lrad .ne. 0)rtwid (i    ) = dtwid (i    ) * rtwid (i    )
50      continue
c
c  Save the mass flux outside (ibeg:iend,jbeg:jend,kbeg:kend)
c  only for zones next to the inner borders.
c
        if ( (j.eq.js-1 .and. k.ge.kbeg)   .or.
     &       (j.ge.jbeg .and. k.eq.ks-1)   .or.
     &       (j.eq.js-1 .and. k.eq.ks-1) ) then
         if (ibeg.eq.is) mflx(is-1,j,k) = mflux (is-1)
         do 60 i=ibeg,iend
           mflx(i,j,k) = mflux (i      )
60       continue
         if (iend.eq.ie) mflx(ie+1,j,k) = mflux (ie+1)
        endif
c
c  Perform mass density and energy density advection.  Note that
c  the timestep "dt" is hidden the fluxes "dflx" and "eflx".
c  Do only zones inside (ibeg:iend,jbeg:jend,kbeg:kend).
c
        if (j.ge.jbeg .and. k.ge.kbeg) then
         if (ibeg.eq.is) mflx(is-1,j,k) = mflux (is-1)
          do 80 i=ibeg,iend
           mflx(i,j,k)= mflux (i      )
           if(xvgrid) then
            den(i,j,k) = ( dlo(i,j,k) * dvl1a(i)
     1                 - dtwid(i+1    )+dtwid (i    ) )*dvl1ani(i)
           else
            den(i,j,k) = ( dlo(i,j,k) * dvl1a(i)
     1                 - dtwid(i+1    )+dtwid (i    ) )*dvl1a i(i)
           endif
           if(nspec .gt. 1) then
            do n = 1, nspec
             if(xvgrid) then
              abn(i,j,k,n) = (dlo(i,j,k)*abo(i,j,k,n)*dvl1a(i)
     .                     -  xtwid(i+1,n) + xtwid(i,n))*dvl1ani(i)
     .                     / den(i,j,k)
             else
              abn(i,j,k,n) = (dlo(i,j,k)*abo(i,j,k,n)*dvl1a(i)
     .                     -  xtwid(i+1,n) + xtwid(i,n))*dvl1a i(i)
     .                     / den(i,j,k)
             endif
            enddo ! n
           endif ! nspec
           if(xiso .eqv. .false.) then
            if(xvgrid) then
             e(i,j,k) = ( e  (i,j,k) * dvl1a(i)
     1                  - etwid(i+1) + etwid (i) ) * dvl1ani(i)
            else
             e(i,j,k) = ( e  (i,j,k) * dvl1a(i)
     1                -   etwid(i+1) + etwid (i) ) * dvl1a i(i)
            endif
c
c Compute e/d for the next substep.
c
            if(ldimen .eq. 3) then
             kp1 = k+1
            else
             kp1 = ks
            endif
            if(xtotnrg .eqv. .false.) then
             edn(i,j,k) = e(i,j,k) / den(i,j,k)
            else ! xtotnrg
             edn(i,j,k) = gamma * e(i,j,k) / den(i,j,k)
     1                  - gamm1 * ( (v1(i,j,k)+v1(i+1,j  ,k  ))**2
     2                            + (v2(i,j,k)+v2(i  ,j+1,k  ))**2
     3                            + (v3(i,j,k)+v3(i  ,j  ,kp1))**2 )
     4                          * 0.125
            endif ! xtotnrg
           endif ! xiso
           if(lrad .ne. 0) then
            if(xvgrid) then
             ern(i,j,k) = ( ero(i,j,k) * dlo(i,j,k) * dvl1a(i)
     1                  -   rtwid(i+1) + rtwid(i) ) * dvl1ani(i)
            else
             ern(i,j,k) = ( ero(i,j,k) * dlo(i,j,k) * dvl1a(i)
     1                  -   rtwid(i+1) + rtwid(i) ) * dvl1a i(i)
            endif ! xvgrid
c
c Work with er/d for the next sweep.
c
            ern(i,j,k) = ern(i,j,k) / den(i,j,k)
           endif ! lrad
80        continue
          if (iend.eq.ie) mflx(ie+1,j,k) = mflux (ie+1)
         endif
90      continue
100    continue
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                T R A N X 1                \\\\\\\\\\
c
c=======================================================================
c
c
