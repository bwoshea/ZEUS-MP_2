c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 A V I S C                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine avisc (ibeg,iend,jbeg,jend,kbeg,kend,dvdxmn
     &                  ,w1,w2,w3,u1,u2,u3,s1,s2,s3)
c
c    mln:zeus3d.viscous <------------- artificial viscosity source terms
c                                                        ?????????, 19??
c
c    written by: Mike Norman
c    modified 1: June, 1988 by Jim Stone; incorporated into ZEUS2D
c    modified 2: February, 1990 by David Clarke; incorporated into
c                ZEUS3D
c    modified 3: June, 1992 by David Clarke; expunged "ISMIN", thereby
c                decreasing execution time by 30%.
c    modified 4: Oct., 1994 by Robert Fiedler to run in parallel on SGIs
c    modified 5: Totally rewritten 2/27/96 by RAF for ZEUS-MP.
c    modified 6: October 2005 by John Hayes; returned linear viscosity
c                terms, after ZEUS-3D.
c
c  PURPOSE: Computes the artificial viscosity source terms in the
c  momentum and energy equations.  i.e., it computes
c
c             dv / dt = -DIV(Q) / rho             for v1, v2, and v3
c      and    de / dt = -Q * delta(v) / delta(x)
c
c  This routine uses the von Neumann-Richtmyer form of the artificial
c  viscosity.  This means that geometric terms are not included in
c  DIV(Q).
c
c
c  LOCAL VARIABLES:
c    qqs        diagonal elements of viscous tensor.  Thus, this is a
c               linear treatment of the artificial viscosity.
c    dvelb      v1(i+1,j,k)-v1(i,j,k) for i-sweep
c               v2(i,j+1,k)-v2(i,j,k) for j-sweep
c               v3(i,j,k+1)-v3(i,j,k) for k-sweep
c    dvela      min ( zro, dvelb ) - ensures that only compressional
c               waves (shocks) are affected.
c    dvdxmn     min ( ( delta(v1) / delta(x1) ),
c                     ( delta(v2) / delta(x2) ),
c                     ( delta(v3) / delta(x3) ) )
c    w1,w2,w3   velocity values prior to viscosity update.
c    u1,u2,u3   velocity values after    viscosity update.
c    s1,s2,s3   updated momentum densities for transport step.
c
c BOUNDARY VALUES USED:
c
c  Macro defined  var   ii    oi    ij    oj    ik    ok
c  -------------  ---  ----  ----  ----  ----  ----  ----
c                  d   is-1        js-1        ks-1
c                  u1  is-1  ie+1
c                  u2              js-1  je+1
c                  u3                          ks-1  ke+1
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
      implicit none
c
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend, km1,
     .            jm1
      real(rl) :: q1, q3, g3i, q2, q11p, q11m, q22p, q22m, q33p, q33m,
     .            qt
      real(rl) :: qqs,qqsm, qqsp, dvelas,dvelasm, dvelasp, dvdxmn,
     .            cs000, csm00, cs0m0, cs00m, dvelbs, dvelbsm
      real(rl) :: w1(in,jn,kn),w2(in,jn,kn),w3(in,jn,kn),
     .            u1(in,jn,kn),u2(in,jn,kn),u3(in,jn,kn),
     .            s1(in,jn,kn),s2(in,jn,kn),s3(in,jn,kn)
c
c
c-----------------------------------------------------------------------
c      Start artificial viscosity update.
c-----------------------------------------------------------------------
c
      if(xsubav) then
       q1 = avisc_dt * qcon
       q2 = avisc_dt * qlin
      else
       q1 = dt * qcon
       q2 = dt * qlin
      endif
c
      do 160 k=kbeg,kend
       if(ldimen .eq. 3) then
        km1 = k - 1
       else
        km1 = k
       endif
       do 150 j=jbeg,jend
        if(ldimen .gt. 1) then
         jm1 = j - 1
        else
         jm1 = j
        endif
        do 140 i=ibeg,iend
         if(qlin .ne. 0.0) then
          cs000 = sqrt(gamm1*gamma*e(i  ,j  ,k  )/d(i  ,j  ,k  ))
          csm00 = sqrt(gamm1*gamma*e(i-1,j  ,k  )/d(i-1,j  ,k  ))
          cs0m0 = sqrt(gamm1*gamma*e(i  ,jm1,k  )/d(i  ,jm1,k  ))
          cs00m = sqrt(gamm1*gamma*e(i  ,j  ,km1)/d(i  ,j  ,km1))
         endif
c
c      Do v1.
c
         dvelasm = min ( zro, w1(i  ,j,k) - w1  (i-1,j,k) )
         dvelas  = min ( zro, w1(i+1,j,k) - w1  (i  ,j,k) )
         if(qlin .eq. 0.0) then
          qqsm    = q1 * d(i-1,j,k) * dvelasm * dvelasm
          qqs     = q1 * d(i  ,j,k) * dvelas  * dvelas
          q3      = dvelas  * dx1ai(i)
         else ! qlin
          dvelbsm = w1(i  ,j,k) - w1  (i-1,j,k)
          dvelbs  = w1(i+1,j,k) - w1  (i  ,j,k)
          qqsm    = d(i-1,j,k)*dvelbsm
     .            * (q1 * dvelasm - q2 * csm00)
          qqs     = d(i  ,j,k)*dvelbs
     .            * (q1 * dvelas  - q2 * cs000)
          q3      = dvelbs  * dx1ai(i)
         endif ! qlin
c
         dvdxmn  = min ( dvdxmn, q3 )
         if(xiso .eqv. .false.) then
          if(xtotnrg .eqv. .false.) then
           e  (i,j,k) = e(i,j,k) - q3 * qqs
          else
           dvelasp = min ( zro, w1(i+2,j,k) - w1  (i+1,j,k) )
           qqsp    = q1 * d(i+1,j,k) * dvelasp * dvelasp
           q11p    = 0.5D0*(qqs + qqsp)
           q11m    = 0.5D0*(qqs + qqsm)
           qt      = (g2a(i+1)*g31a(i+1) * q11p * w1(i+1,j,k)
     1             -  g2a(i  )*g31a(i  ) * q11m * w1(i  ,j,k) )
     2             *                                dvl1ai(i)
           e(i,j,k) = e(i,j,k) - qt
          endif
         endif
         u1 (i,j,k) = w1(i,j,k)
     1              - ( qqs  - qqsm   ) * dx1bi(i)
     2              * 2.0 / ( d(i-1,j,k) + d(i,j,k) )
c
c      Do v2.
c
         if(ldimen .gt. 1) then
          dvelasm = min ( zro, w2(i,j  ,k) - w2  (i,j-1,k) )
          dvelas  = min ( zro, w2(i,j+1,k) - w2  (i,j  ,k) )
          if(qlin .eq. 0.0) then
           qqsm    = q1 * d(i,j-1,k) * dvelasm  * dvelasm
           qqs     = q1 * d(i,j  ,k) * dvelas   * dvelas
           q3      = dvelas   * dx2ai(j) * g2bi(i)
          else ! qlin
           dvelbsm = w2(i,j  ,k) - w2  (i,j-1,k)
           dvelbs  = w2(i,j+1,k) - w2  (i,j  ,k)
           qqsm    = d(i,j-1,k) * dvelbsm
     .             * (q1 * dvelasm - q2 * cs0m0)
           qqs     = d(i,j  ,k) * dvelbs
     .             * (q1 * dvelas  - q2 * cs000)
           q3      = dvelbs   * dx2ai(j) * g2bi(i)
          endif ! qlin
          dvdxmn  = min ( dvdxmn, q3 )
          if(xiso .eqv. .false.) then
           if(xtotnrg .eqv. .false.) then
            e(i,j,k) = e(i,j,k) - q3 * qqs
           else
            dvelasp = min ( zro, w2(i,j+2,k) - w2  (i,j+1,k) )
            qqsp    = q1 * d(i,j+1,k) * dvelasp * dvelasp
            q22p    = 0.5D0*(qqs + qqsp)
            q22m    = 0.5D0*(qqs + qqsm)
            qt      = ( g32a(j+1) * q22p * w2(i,j+1,k)
     4                - g32a(j  ) * q22m * w2(i,j  ,k) )
     5              *   g2bi(i)             * dvl2ai(j)
            e(i,j,k) = e(i,j,k) - qt
           endif ! xtotnrg
          endif
          u2(i,j,k) = w2(i,j,k)
     1              - ( qqs   - qqsm    ) * dx2bi(j) * g2bi(i)
     2              * 2.0 / ( d(i,j-1,k) + d(i,j,k) )
         else ! ldimen
          u2(i,j,k)  = w2(i,j,k)
         endif ! ldimen
c
c      Do v3.
c
         IF(LDIMEN .EQ. 3) THEN
         g3i     = g31bi(i) * g32bi(j)
         dvelasm = min ( zro, w3(i,j,k  ) - w3  (i,j,k-1) )
         dvelas  = min ( zro, w3(i,j,k+1) - w3  (i,j,k  ) )
         if(qlin .eq. 0.0) then
          qqsm    = q1 * d(i,j,k-1) * dvelasm  * dvelasm
          qqs     = q1 * d(i,j,k  ) * dvelas   * dvelas
          q3      = dvelas   * dx3ai(k) * g3i
         else ! qlin
          dvelbsm = w3(i,j,k  ) - w3(i,j,k-1)
          dvelbs  = w3(i,j,k+1) - w3(i,j,k  )
          qqsm    = d(i,j,k-1) * dvelbsm
     .            * (q1 * dvelasm - q2 * cs00m)
          qqs     = d(i,j,k  ) * dvelbs
     .            * (q1 * dvelas  - q2 * cs000)
          q3      = dvelbs   * dx3ai(k) * g3i
         endif ! qlin
         dvdxmn  = min ( dvdxmn, q3 )
         if(xtotnrg) q33(i,j,k) = qcon*dvelas*dvelas
         if((xiso .eqv. .false.) .and. 
     .      (xtotnrg .eqv. .false.)) then
C             if(x1dflag(i) .eqv. .false.) e(i,j,k) = e(i,j,k) - q3*qqs
              e(i,j,k) = e(i,j,k) - q3*qqs
         endif
         u3(i,j,k) = w3(i,j,k)
     1             - ( qqs   - qqsm    ) * dx3bi(k) * g3i
     2             *  2.0 / ( d(i,j,k-1) + d(i,j,k) )
         ELSE
         u3(i,j,k) = w3(i,j,k)
         ENDIF
c
c Save the updated momentum densities for the transport step, since 
c everything should be in the cache already.
c
         s1(i,j,k) = u1(i,j,k) * 0.5 * (d(i-1,j  ,k  ) + d(i,j,k))
         s2(i,j,k) = u2(i,j,k) * 0.5 * (d(i  ,jm1,k  ) + d(i,j,k))
     1             * g2b(i)
         s3(i,j,k) = u3(i,j,k) * 0.5 * (d(i  ,j  ,km1) + d(i,j,k))
     1             * g31b(i) * g32b(j)
140     continue
150    continue
160   continue
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                 A V I S C                 \\\\\\\\\\
c
c=======================================================================
c
c
