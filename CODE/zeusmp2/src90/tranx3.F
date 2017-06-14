c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                T R A N X 3                \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine tranx3 (ibeg,iend,jbeg,jend,kbeg,kend
     &                   ,dlo,den
     &                   ,eod,edn
     &                   ,ero,ern
     &                   ,abo,abn
     &                   ,mflx,atwid,atwid2,dtwid,dd,mflux,etwid,deod)
c
c    dac:zeus3d.tranx3 <----- transports zone-centred variables along x3
c    from jms:zeus2d.tranx2, mln:zeus04.tranz                  may, 1990
c
c    written by: David Clarke
c    modified 1: June 1992, by David Clarke; added the total energy
c                option originally designed by Byung-IL Jun.
c    modified 2: Feb. 22, 1996 by Robert Fiedler; completely rewritten
c                for ZEUS-MP.
c    modified 3: Aug. 6,  1996 by Robert Fiedler; rearranged for
c                maximum efficiency, and unrolled i-loops.
c    modified 3: Dec. 18,  1996 by Robert Fiedler; added radiation.
c    modified 4: by John Hayes; rewrote for F90 and removed CPP logic
c    modified 5: by John Hayes; restored metric factors multiplying "xi"
c    modified 6: by John Hayes; corrected transposed array subscripts in
c                "dxo" array
c
c  PURPOSE:  Transports all zone centred variables in the 3-direction
c  only.  Currently transported are:
c
c                      mass   density
c                      energy density
c
c  The consistent advection algorithm, in which mass fluxes are used to
c  construct the fluxes of all variables across the interfaces, is used
c  for all hydrodynamical variables.  Thus, the mass fluxes are passed
c  to MOMX3 in order to transport the momenta as well.  The magnetic
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
c                    d   is-1        js-1        ks-3  ke+2
c                   e/d                          ks-2  ke+2
c    TOTAL_ENERGY    u1        ie+1
c    TOTAL_ENERGY    u2                    je+1
c                    u3  is-1        js-1        ks-1  ke+1
c
c  OUTPUT VARIABLES:
c    den      updated mass            density.
c    edn      updated specific energy density.
c    mflx     mass flux (in the 3-direction)
c
c  LOCAL VARIABLES:
c    atwid    effective cross sectional area of the 1-interfaces
c    atwid2   effective cross sectional area of the 2-interfaces
c    etwid    interpolated specific energy densities (e/d) at all
c             3-interfaces
c    eflx     energy density flux across all 3-interfaces (reuse etwid)
c    dtwid    interpolated mass densities at all 3-interfaces
c    dflx     mass density flux across all 3-interfaces   (reuse dtwid)
c
c  EXTERNALS:
c
      use real_prec
      use config
      use param
      use root
      use field
      use bndry
      use grid
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend, n
c
      real(rl) :: dqm, dqp, xi, q1
c
      real(rl) :: atwid2(ijkn),  atwid(ijkn)
      real(rl) :: mflux (ijkn,1),
     .            dtwid (ijkn,1), dd  (ijkn,1),
     .            etwid (ijkn,1), deod(ijkn,1),
     .            rtwid (ijkn,1), dero(ijkn,1),
     .            xtwid(ijkn,nspec),dxo(ijkn,nspec)
c
      real(rl) :: mflx(in,jn,kn),
     .            dlo (in,jn,kn), den(in,jn,kn),
     .            eod (in,jn,kn), edn(in,jn,kn),
     .            ero (in,jn,kn), ern(in,jn,kn),
     .            abo (in,jn,kn,nspec), abn(in,jn,kn,nspec)
c-----------------------------------------------------------------------
c
c Compute time-centered area factors.
c
      do 10 i=ibeg-1,iend+1
       atwid (i)           = g2b(i) * dx1a(i) * dvl1ai(i)
10    continue
      do 20 j=jbeg-1,jend+1
       atwid2(j)           = dx2a(j) * dvl2ai(j)
20    continue
c
c Transport all zone-centered quantities in the 3 direction.
c
      do 2100 j=jbeg-1,jend
       if ( (j.ge.jbeg) .or. (j.eq.js-1) ) then
c.......................................................................
c
c Split off the i=ibeg-1 iteration to ease unrolling.
c
        i = ibeg - 1
c
c   Interpolate to obtain zone-centered quantities at zone faces.
c
c     1.  Evaluate monotonised, van Leer differences across the zone.
c
        if (kbeg .eq. ks) then  !  Need d(ks-3) from neighbor.
         k        = ks - 2
         dqm      = (dlo(i  ,j,k  ) - dikb (i  ,j,3)) * dx3bi(k  )
         dqp      = (dlo(i  ,j,k+1) - dlo(i  ,j,k  )) * dx3bi(k+1)
         dd(k,1)  = max ( dqm * dqp, zro )
     1            * sign ( one, dqm + dqp )
     2            / max ( abs ( dqm + dqp ), tiny )
        endif
        do 30 k=max(kbeg-2,ks-1),kend+1
         dqm      = (dlo(i  ,j,k  ) - dlo(i  ,j,k-1)) * dx3bi(k  )
         dqp      = (dlo(i  ,j,k+1) - dlo(i  ,j,k  )) * dx3bi(k+1)
         dd(k,1)  = max ( dqm * dqp, zro )
     1            * sign ( one, dqm + dqp )
     2            / max ( abs ( dqm + dqp ), tiny )
         if(nspec .gt. 1) then
          do n = 1, nspec
           dqm         = (abo(i,j,k  ,n)-abo(i,j,k-1,n))*dx3bi(k  )
           dqp         = (abo(i,j,k+1,n)-abo(i,j,k  ,n))*dx3bi(k+1)
           dxo(k,n) = max ( dqm * dqp, zro )
     1                * sign ( one, dqm + dqp )
     2                / max ( abs ( dqm + dqp ), tiny )
          enddo
         endif
30      continue
c
c     2.  Choose time averaged, upwinded interface value.
c
c  For the purposes of consistent advection, construct the mass
c  flux across each 1-interface.  The mass flux will be used to create
c  the fluxes of all variables, including the momenta which are updated
c  in MOMX1.
c
        do 40 k=kbeg-1,kend+1
         xi          = ( v3  (i  ,j,k  ) - vg3(k  ) ) * dt
     .               * g31bi(i) * g32bi(j)
         q1          = sign ( haf, xi )
         dtwid (k,1) = ( 0.5 + q1 ) * ( dlo(i  ,j,k-1)
     1               + ( dx3a(k-1) - xi ) * dd   (k-1,1) )
     2               + ( 0.5 - q1 ) * ( dlo(i  ,j,k  )
     3               - ( dx3a(k  ) + xi ) * dd   (k  ,1) )
c
         mflux (k,1) = dtwid (k,1) * ( v3(i  ,j,k) - vg3(k) ) * dt
40      continue
c
c  Save the mass flux outside (ibeg:iend,jbeg:jend,kbeg:kend)
c  only for zones along the inner borders.
c
        if (kbeg .eq. ks) mflx(i  ,j,ks-1) = mflux (ks-1,1)
         do 60 k=kbeg,kend
          mflx(i  ,j,k) = mflux (k   ,1)
60       continue
         if (kend .eq. ke) mflx(i  ,j,ke+1) = mflux (ke+1,1)
c
c.......................................................................
c
         do 1090 i=ibeg,iend
c
c   Interpolate to obtain zone-centered quantities at zone faces.
c
c     1.  Evaluate monotonised, van Leer differences across the zone.
c
          if (kbeg .eq. ks) then  !  Need d(ks-3) from neighbor.
           k        = ks - 2
           dqm      = (dlo(i  ,j,k  ) - dikb (i  ,j,3)) * dx3bi(k  )
           dqp      = (dlo(i  ,j,k+1) - dlo(i  ,j,k  )) * dx3bi(k+1)
           dd(k,1)  = max ( dqm * dqp, zro )
     1              * sign ( one, dqm + dqp )
     2              / max ( abs ( dqm + dqp ), tiny )
           if(xiso .eqv. .false.) then
            deod(k,1) = zro  ! Not valid, but we don't use it.
           endif
           if(lrad .ne. 0) dero(k,1) = zro
          endif
          do 1030 k=max(kbeg-2,ks-1),kend+1
           dqm      = (dlo(i  ,j,k  ) - dlo(i  ,j,k-1)) * dx3bi(k  )
           dqp      = (dlo(i  ,j,k+1) - dlo(i  ,j,k  )) * dx3bi(k+1)
           dd(k,1)  = max ( dqm * dqp, zro )
     1              * sign ( one, dqm + dqp )
     2              / max ( abs ( dqm + dqp ), tiny )
           if(nspec .gt. 1) then
            do n = 1, nspec
             dqm      = (abo(i,j,k  ,n)-abo(i,j,k-1,n))*dx3bi(k  )
             dqp      = (abo(i,j,k+1,n)-abo(i,j,k  ,n))*dx3bi(k+1)
             dxo(k,n) = max ( dqm * dqp, zro )
     1                  * sign ( one, dqm + dqp )
     2                  / max ( abs ( dqm + dqp ), tiny )
            enddo
           endif
           if(xiso .eqv. .false.) then
            dqm      = (eod(i  ,j,k  ) - eod(i  ,j,k-1)) * dx3bi(k  )
            dqp      = (eod(i  ,j,k+1) - eod(i  ,j,k  )) * dx3bi(k+1)
            deod(k,1)= max ( dqm * dqp, zro )
     1               * sign ( one, dqm + dqp )
     2               / max ( abs ( dqm + dqp ), tiny )
           endif
           if(lrad .ne. 0) then
            dqm      = (ero(i  ,j,k  ) - ero(i  ,j,k-1)) * dx3bi(k  )
            dqp      = (ero(i  ,j,k+1) - ero(i  ,j,k  )) * dx3bi(k+1)
            dero(k,1)= max ( dqm * dqp, zro )
     1               * sign ( one, dqm + dqp )
     2               / max ( abs ( dqm + dqp ), tiny )
           endif
1030      continue
c
c     2.  Choose time averaged, upwinded interface value.
c
c  For the purposes of consistent advection, construct the mass
c  flux across each 1-interface.  The mass flux will be used to create
c  the fluxes of all variables, including the momenta which are updated
c  in MOMX1.
c
          do 1040 k=kbeg-1,kend+1
           xi          = ( v3  (i  ,j,k  ) - vg3(k  ) ) * dt
     .                 * g31bi(i) * g32bi(j)
           q1          = sign ( haf, xi )
           dtwid (k,1) = ( 0.5 + q1 ) * ( dlo(i  ,j,k-1)
     1                 + ( dx3a(k-1) - xi ) * dd   (k-1,1) )
     2                 + ( 0.5 - q1 ) * ( dlo(i  ,j,k  )
     3                 - ( dx3a(k  ) + xi ) * dd   (k  ,1) )
c
           mflux (k,1) = dtwid (k,1) * ( v3(i  ,j,k) - vg3(k) ) * dt
           dtwid (k,1) = mflux (k,1) * atwid (i  ) * atwid2(j)
           if(nspec .gt. 1) then
            do n = 1, nspec
             xtwid(k,n) = ( 0.5 + q1 ) * (abo(i,j,k-1,n)
     1                    + ( dx3a(k-1) - xi ) * dxo(k-1,n) )
     2                    + ( 0.5 - q1 ) * (abo(i,j,k  ,n)
     3                    - ( dx3a(k  ) + xi ) * dxo(k  ,n) )
            enddo
            do n = 1, nspec
             xtwid(k,n) = xtwid(k,n)*dtwid(k,1)
            enddo
           endif ! nspec
           if(xiso .eqv. .false.) then
            etwid (k,1) = ( 0.5 + q1 ) * ( eod(i  ,j,k-1)
     1                  + ( dx3a(k-1) - xi ) * deod (k-1,1) )
     2                  + ( 0.5 - q1 ) * ( eod(i  ,j,k  )
     3                  - ( dx3a(k  ) + xi ) * deod (k  ,1) )
c
            etwid (k,1) = dtwid (k,1) * etwid (k,1)
           endif
           if(lrad .ne. 0) then
            rtwid (k,1) = ( 0.5 + q1 ) * ( ero(i  ,j,k-1)
     1                  + ( dx3a(k-1) - xi ) * dero (k-1,1) )
     2                  + ( 0.5 - q1 ) * ( ero(i  ,j,k  )
     3                  - ( dx3a(k  ) + xi ) * dero (k  ,1) )
c
            rtwid (k,1) = dtwid (k,1) * rtwid (k,1)
           endif
1040      continue
c
c  Save the mass flux outside (ibeg:iend,jbeg:jend,kbeg:kend)
c  only for zones along the inner borders.
c
          if (j.eq.js-1) then
           if (kbeg .eq. ks) mflx(i  ,j,ks-1) = mflux (ks-1,1)
           do 1060 k=kbeg,kend
            mflx(i  ,j,k) = mflux (k   ,1)
1060       continue
           if (kend .eq. ke) mflx(i  ,j,ke+1) = mflux (ke+1,1)
          endif
c
c  Perform mass density and energy density advection.  Note that
c  the timestep "dt" is hidden the fluxes "dtwid" and "etwid".
c
          if (j.ge.jbeg) then
           if (kbeg .eq. ks) mflx(i  ,j,ks-1) = mflux (ks-1,1)
           do 1080 k=kbeg,kend
            mflx(i  ,j,k)= mflux (k   ,1)
            if(xvgrid) then
             den(i,j,k) = ( dlo(i  ,j,k) * dvl3a(k)
     .                  - dtwid(k+1,1)+dtwid (k,1)) * dvl3ani(k)
            else
             den(i,j,k) = ( dlo(i  ,j,k) * dvl3a(k)
     .                  - dtwid(k+1,1)+dtwid (k,1)) * dvl3a i(k)
            endif
            if(nspec .gt. 1) then
             do n = 1, nspec
              if(xvgrid) then
               abn(i,j,k,n) = (abo(i,j,k,n)*dlo(i,j,k)*dvl3a(k)
     .                      -  xtwid(k+1,n)+xtwid(k,n)) *
     .                         dvl3ani(k)/den(i,j,k)
              else
               abn(i,j,k,n) = (abo(i,j,k,n)*dlo(i,j,k)*dvl3a(k)
     .                      -  xtwid(k+1,n)+xtwid(k,n)) *
     .                         dvl3a i(k)/den(i,j,k)
              endif ! xvgrid
             enddo
            endif ! nspec
            if(xiso .eqv. .false.) then
             if(xvgrid) then
              e(i,j,k) = ( e  (i  ,j,k) * dvl3a(k)
     .                 - etwid (k+1,1) + etwid (k,1) ) * dvl3ani(k)
             else
              e(i,j,k) = ( e  (i  ,j,k) * dvl3a(k)
     .                 - etwid (k+1,1) + etwid (k,1) ) * dvl3a i(k)
             endif
c
c Compute e/d for the next substep.
c
             if(xtotnrg .eqv. .false.) then
              edn(i,j,k) =         e(i,j,k) / den(i,j,k)
             else !xtotnrg
              edn(i,j,k) = gamma * e(i,j,k) / den(i,j,k) - gamm1
     .                   * (
     .                     ( v1(i  ,j,k) + v1(i+1,j  ,k  ) )**2
     .                   + ( v2(i  ,j,k) + v2(i  ,j+1,k  ) )**2
     .                   + ( v3(i  ,j,k) + v3(i  ,j  ,k+1) )**2
     .                      ) * 0.125
             endif ! xtotnrg
            endif ! xiso
            if(lrad .ne. 0) then
             if(xvgrid) then
              ern(i,j,k) = (ero(i,j,k) * dlo(i  ,j,k) * dvl3a(k)
     1                     - rtwid(k+1,1)+rtwid (k,1) ) * dvl3ani(k)
             else
              ern(i,j,k) = (ero(i,j,k) * dlo(i  ,j,k) * dvl3a(k)
     1                   -   rtwid(k+1,1)+rtwid (k,1) ) * dvl3a i(k)
            endif
c
c Work with er/d for the next sweep.
c
            ern(i  ,j,k) = ern(i  ,j,k) / den(i  ,j,k)
           endif ! lrad
1080      continue
          if (kend .eq. ke) mflx(i  ,j,ke+1) = mflux (ke+1,1)
         endif
1090    continue
       endif  !  j>=jbeg || j=js-1
2100  continue
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                T R A N X 3                \\\\\\\\\\
c
c=======================================================================
c
c
