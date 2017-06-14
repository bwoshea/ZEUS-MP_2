c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                F O R C E S                \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine forces (ibeg,iend,jbeg,jend,kbeg,kend
     &                   ,u1,u2,u3,w1,w2,w3)
c
c
c Computes Pressure, Gravity, and Pseudo-Rotational Forces.
c   and includes magnetic pressure forces (longitudinal Lorentz forces)
c
c Arrays u1 , u2 , u3  hold the old velocity values, while
c        w1 , w2 , w3  receive the updated values.
c
c BOUNDARY VALUES USED:
c
c  Macro defined  var   ii    oi    ij    oj    ik    ok
c  -------------  ---  ----  ----  ----  ----  ----  ----
c                  d   is-1        js-1        ks-1
c  ZRP             d   is-1        js-1        ks-1  ke+1
c  RTP             d   is-1        js-1  je+1  ks-1  ke+1
c                  e   is-1        js-1        ks-1
c  TOTAL_ENERGY    u1  is-1  ie+1  js-1        ks-1
c  RTP             u2  is-1              je+1
c  TOTAL_ENERGY    u2  is-1        js-1  je+1  ks-1
c  TOTAL_ENERGY    u3  is-1        js-1        ks-1  ke+1
c  ZRP             u3              js-1              ke+1
c  RTP             u3  is-1        js-1              ke+1
c  RAD             er  is-1        js-1        ks-1
c  MHD             b1  is-1  ie+1  js-1  je+1  ks-1  ke+1
c  MHD             b2  is-1  ie+1  js-1  je+1  ks-1  ke+1
c  MHD             b3  is-1  ie+1  js-1  je+1  ks-1  ke+1
c Written by RAF; modified 3/13/97 (JCH); 
c     modified 2 Mar 1998 to add MHD (M-MML)
c
c Rewritten for F90 by John Hayes, God knows when, exactly.
c
c Modified: May 2006, by John Hayes: corrected k indexing error
c in expression for derdx3 in v2 update.
c Modified: 09/01/2006, by John Hayes; implemented Sean Matt's
c corrections to components of point-mass potential in RTP geometry
c Modified: 12/20/2006, by John Hayes; corrected typo in formula
c for corrections to v2 due to rotational pseudoforces.  This bug
c was inherited from ZEUS-3D!!
c 
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
      use gravmod
      use cons
      use opac
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar

      implicit NONE
c
      integer  :: ibeg,iend,jbeg,jend,kbeg,kend,km1,kp1
      integer  :: jp1, jm1
      integer  :: i,j,k
c
      real(rl) :: u1 (in,jn,kn),u2 (in,jn,kn),u3 (in,jn,kn),
     .            w1 (in,jn,kn),w2 (in,jn,kn),w3 (in,jn,kn)
c
      real(rl) :: poo(in),pmo(in),pom(in)
      real(rl) :: st,rhoi
c
      real(rl) :: s2oo(in),s2po(in)
c
      real(rl) :: s3oo(in),s3mo(in),s3op(in),s3mp(in)
c
      real(rl) :: r2,r2i,prpq,pzpq,prpt,pzpt
      real(rl) :: d1b2oo(in), d1b3oo(in), d2b1oo(in), d2b3oo(in),
     &            d3b1oo(in), d3b2oo(in)
      real(rl) :: d1b2po(in), d1b3op(in), d2b1po(in), d2b3op(in),
     &            d3b1po(in), d3b2op(in)
c
      real(rl) :: rdmcp, ros_mfp, flx_lim, big_r, derdx1, derdx2,
     .            derdx3, rdmcm, grad_er
c
c-----------------------------------------------------------------------
c
c Main loop for pressure, magnetic pressure, self-gravity, and
c     rotational pseudo-forces.  
c
c Active zones for all velocity components are (is:ie, js:je, ks:ke).
c
c In its most general configuration, this loop needs 1 layer of 
c boundary data for each of d, e, u1, u2, u3, and gp.
c
c If GRAV or GRAV_FFT is not defined, we don't need gp.
c If ISO  is     defined, we don't need e .
c
c For CARTESIAN grids and TOTAL_ENERGY not defined, we don't need u1,u2,u3
c boundary values.  All we do need in this case is the inner boundary 
c data for d (and e).
c
c Of course, TOTAL_ENERGY and ISO are mutually exclusive. 
c
c **> at the moment, TOTAL_ENERGY and MHD also can't be used together
c       (M-MML 4 Mar 98)
c
c
       do 40 k=kbeg,kend
        if(ldimen .eq. 3) then
         km1 = k-1
         kp1 = k+1
        else
         km1 = ks
         kp1 = ks
        endif
        do 30 j=jbeg,jend
         if(ldimen .eq. 1) then
          jm1 = js
          jp1 = js
         else
          jm1 = j-1
          jp1 = j+1
         endif
         do 10 i=ibeg-1,iend
c
c Compute thermal pressure
c
           poo(i) = p(i  ,j  ,k  )
           pmo(i) = p(i  ,jm1,k  )
           pom(i) = p(i  ,j  ,km1)
c.......................................................................
c   
          if(lgeom .ne. 1) then
c
c      Construct momentum densities from velocities.
c
           s3oo(i)         = u3(i  ,j  ,k  )*0.5*g31b(i  )*g32b(j  )
     1                     * ( d(i  ,j  ,km1) + d(i  ,j  ,k  ) )
           s3op(i)         = u3(i  ,j  ,kp1)*0.5*g31b(i  )*g32b(j  )
     1                     * ( d(i  ,j  ,k  ) + d(i  ,j  ,kp1) )
           s3mo(i)         = u3(i  ,jm1,k  )*0.5*g31b(i  )*g32b(jm1)
     1                     * ( d(i  ,jm1,km1) + d(i  ,jm1,k  ) )
           s3mp(i)         = u3(i  ,jm1,kp1)*0.5*g31b(i  )*g32b(jm1)
     1                     * ( d(i  ,jm1,k  ) + d(i  ,jm1,kp1) )
          endif ! lgeom 
          if(lgeom .eq. 3) then
           s2oo(i)         = u2(i  ,j  ,k  )*0.5*g2b (i  )
     1                     * ( d(i  ,jm1,k  ) + d(i  ,j  ,k  ) )
           s2po(i)         = u2(i  ,jp1,k  )*0.5*g2b (i  )
     1                     * ( d(i  ,j  ,k  ) + d(i  ,jp1,k  ) )
          endif ! lgeom
 10      continue
c.......................................................................
         if(xmhd) then
c
c        Compute differences in the squares of the magnetic field
c        components for the longitudinal Lorentz forces.
c
          do 11 i=ibeg-1,iend
           d1b2oo(i) = ( ( g2b (i  ) * b2(i  ,j  ,k  ) )**2
     1                 - ( g2b (i-1) * b2(i-1,j  ,k  ) )**2 )
     2                  * g2ai (i) * g2ai (i)
           d1b2po(i) = ( ( g2b (i  ) * b2(i  ,jp1,k  ) )**2
     1                 - ( g2b (i-1) * b2(i-1,jp1,k  ) )**2 )
     2                  * g2ai (i) * g2ai (i)
           d1b3oo(i) = ( ( g31b(i  ) * b3(i  ,j  ,k  ) )**2
     1                 - ( g31b(i-1) * b3(i-1,j  ,k  ) )**2 )
     2                  * g31ai(i) * g31ai(i)
           d1b3op(i) = ( ( g31b(i  ) * b3(i  ,j  ,kp1) )**2
     1                 - ( g31b(i-1) * b3(i-1,j  ,kp1) )**2 )
     2                  * g31ai(i) * g31ai(i)
           d2b3oo(i) = ( ( g32b(j  ) * b3(i  ,j  ,k  ) )**2
     1                 - ( g32b(jm1) * b3(i  ,jm1,k  ) )**2 )
     2                  * g32ai(j) * g32ai(j)
           d2b3op(i) = ( ( g32b(j  ) * b3(i  ,j  ,kp1) )**2
     1                 - ( g32b(jm1) * b3(i  ,jm1,kp1) )**2 )
     2                  * g32ai(j) * g32ai(j)
           d2b1oo(i) = b1(i  ,j  ,k  ) * b1(i  ,j  ,k  )
     1               - b1(i  ,jm1,k  ) * b1(i  ,jm1,k  )
           d2b1po(i) = b1(i+1,j  ,k  ) * b1(i+1,j  ,k  )
     1               - b1(i+1,jm1,k  ) * b1(i+1,jm1,k  )
           d3b1oo(i) = b1(i  ,j  ,k  ) * b1(i  ,j  ,k  )
     1               - b1(i  ,j  ,km1) * b1(i  ,j  ,km1)
           d3b1po(i) = b1(i+1,j  ,k  ) * b1(i+1,j  ,k  )
     1               - b1(i+1,j  ,km1) * b1(i+1,j  ,km1)
           d3b2oo(i) = b2(i  ,j  ,k  ) * b2(i  ,j  ,k  )
     1               - b2(i  ,j  ,km1) * b2(i  ,j  ,km1)
           d3b2op(i) = b2(i  ,jp1,k  ) * b2(i  ,jp1,k  )
     1               - b2(i  ,jp1,km1) * b2(i  ,jp1,km1)
11        continue
         endif ! xmhd

         do 20 i=ibeg,iend
c.......................................................................
c
c Perform an explicit update for v1
c
#ifndef ROTOR
          rhoi    = 2.0 / ( d(i-1,j,k) + d(i,j,k) )
c
c  1.  pressure gradient
c
          st      = - rhoi
     1            * ( poo(i)    - poo(i-1)    ) * dx1bi(i)
c
c  2.  gravitational potential gradient
c
          if(xgrav .or. xgrvfft) then
           if(xsphgrv) then
            if(x1a(i) .ne. 0.0) then
             st      = st - guniv*intm(i)*x1ai(i)*x1ai(i)
            endif
           else
            st      = st 
     .              + ( gp(i,j,k) - gp(i-1,j,k) ) * dx1bi(i)
           endif
          endif
c
c  3.  rotational pseudo-forces
c
          if(lgeom .eq. 3) then
           st      = st 
     1             + ( s2oo(i  )     + s2po(i  )
     2               + s2oo(i-1)     + s2po(i-1)   )
     3             * ( u2(i  ,j,k) + u2(i  ,jp1,k  )
     4               + u2(i-1,j,k) + u2(i-1,jp1,k  ) ) / 16.0
     5             * rhoi * dg2bd1 (i) * g2a i(i)    * g2a i(i)
           st      = st 
     1             + ( s3oo(i  )     + s3op(i  )
     2               + s3oo(i-1)     + s3op(i-1)   )
     3             * ( u3(i  ,j,k) + u3(i  ,j  ,kp1)
     4               + u3(i-1,j,k) + u3(i-1,j  ,kp1) ) / 16.0
     5             * rhoi * dg31bd1(i) * g31ai(i)**2 * g32bi(j)
          endif ! lgeom
          if(xptmass) then
c
c  4.  gravitational point mass 
c
           if(lgeom .eq. 1) then
             r2        = ( x1a(i) - x1ptm )**2
     1                 + ( x2b(j) - x2ptm )**2
     2                 + ( x3b(k) - x3ptm )**2
             r2i       = ( x1a(i) - x1ptm )
     1                 / ( r2    * sqrt(r2   ) + tiny )
           endif ! CARTESIAN
           if(lgeom .eq. 2) then
             r2        = ( x1a(i) - x1ptm )**2
     1                 + ( x2b(j) - x2ptm )**2
     2                 + 2.0 * x2b(j) * x2ptm
     3                 * ( 1.0 - cos(x3b(k) - x3ptm) )
             r2i       = ( x1a(i) - x1ptm )
     1                 / ( r2    * sqrt(r2   ) + tiny )
           endif ! CYLINDRICAL
           if(lgeom .eq. 3) then
CJH
CJH  correction due to Sean Matt -- enforce condition that radial
CJH  force due to point mass depend only on radius.
CJH
            if(x1ptm .eq. 0.0) then
             r2i = x1ai(i)**2
            else !x1ptm
             prpq      = x1a(i) * sin(x2b(j))
             pzpq      = x1a(i) * cos(x2b(j))
             prpt      = x1ptm  * sin(x2ptm)
             pzpt      = x1ptm  * cos(x2ptm)
c
             r2        = ( pzpq - pzpt )**2
     1                 + ( prpq - prpt )**2
     2                 + 2.0 * prpq * prpt
     3                 * ( 1.0 - cos(x3b(k) - x3ptm) )
             r2i       = ( ( pzpq - pzpt ) * cos(x2b(j))
     1                 + ( prpq - prpt * cos(x3b(k) - x3ptm))
     2                 * sin(x2b(j))) / ( r2    * sqrt(r2   ) + tiny )
            endif ! x1ptm
           endif ! SPHERICAL
            st        = st - guniv * ptmass * r2i
          endif ! POINT MASS GRAVITY
c
c  5.  Radiation force
c
          if(lrad .eq. 1) then
           derdx1 = ( er(i  ,j  ,k  ) - er(i-1,j  ,k  ) )
     1            * dx1bi(i)
           if(ldimen .eq. 1) then
             derdx2 = 0.0
           else
            derdx2 = ( er(i-1,j  ,k  ) - er(i-1,j-1,k  ) )
     1             * dx2bi(j)
     2             + ( er(i  ,j  ,k  ) - er(i  ,j-1,k  ) )
     3             * dx2bi(j)
            derdx2 = derdx2 * 0.5 * g2ai(i)
           endif ! ldimen = 1
           if(ldimen .lt. 3) then
            derdx3 = 0.0
           else
            derdx3 = ( er(i-1,j  ,k  ) - er(i-1,j  ,km1) )
     1             * dx3bi(k)
     2             + ( er(i  ,j  ,k  ) - er(i  ,j  ,km1) )
     3             * dx3bi(k)
            derdx3 = derdx3 * 0.5 * g31ai(i) * g32bi(j)
           endif ! ldimen < 3
           grad_er    = sqrt( derdx1**2 + derdx2**2 + derdx3**2 )
           rdmcm   = kapr(i-1,j,k)
           rdmcp   = kapr(i  ,j,k)
           ros_mfp = 2.0 / (rdmcm + rdmcp)
           big_R   = 2.0 * grad_er * ros_mfp
     1             / max(er(i  ,j,k)+er(i-1,j,k), tiny)
c
#ifdef MINERBO
c
c          Minerbo flux limiter
c
           if (big_R .le. 1.5) then
            flx_lim = 2.0/(3.0 + sqrt(9.0 + 12.0*big_R**2) )
           else
            flx_lim = 1.0/(1.0 + big_R + sqrt(1.0 + 2.0*big_R) )
           endif
#else /* MINERBO */
c
c          LP flux limiter
c
           flx_lim = (2.0 + big_R)/(6.0 + 3.0*big_R + big_R**2)
#endif /* MINERBO */
           st         = st - rhoi * flx_lim * derdx1
          endif ! lrad
c
c  6. Magnetic pressure (longitudinal Lorentz force)
c
          if(xmhd) then
             st = st - (d1b2oo(i) + d1b2po(i) + d1b3oo(i) + d1b3op(i))
     3                 * rhoi * dx1bi(i) * 0.25
          endif ! xmhd
#else /* ROTOR */
          st = 0.0 ! no compressive waves in this test
#endif /* ROTOR */
c
          w1(i,j,k) = u1(i,j,k) + dt * st
c.......................................................................
c
c Perform an explicit update for v2 
c
         if(ldimen .ne. 1) then
          rhoi    = 2.0 / ( d(i,j-1,k) + d(i,j,k) )
c
c  1.  pressure gradient
c
          if(poo(i  ) .lt. 0.0D0)
     .     write(2,"('Pressure = ',1pd12.4,' at i,j,k = ',3i4,
     .               ' on PID ',i3,'; nhy = ',i6)")poo(i),i,j,k,myid,nhy
          if((pmo(i  ) .lt. 0.0D0) .and. (j .eq. js))
     .     write(2,"('Pressure = ',1pd12.4,' at i,j,k = ',3i4,
     .               ' on PID ',i3,'; nhy = ',i6)")pmo(i),i,j-1,k,myid,
     .                 nhy
          st      = - rhoi
     1            * ( poo(i)     - pmo(i)     ) * dx2bi(j)
     2            * g2bi(i)
c  2.  gravitational potential gradient
c
          if(xgrav .or. xgrvfft) then
           if(.not. xsphgrv) then
            st      = st 
     .              + ( gp(i,j,k) - gp(i,j-1,k) ) * dx2bi(j)
     .              * g2bi(i)
           endif
          endif
c
c  3.  rotational pseudo-force
c
#ifndef ROTOR
          if(lgeom .ne. 1) then
           if(x2a(j) .ne. 0.0) then
             st      = st 
     1               + ( s3oo(i)     + s3op(i)
     2                 + s3mo(i)     + s3mp(i)      )
     3               * ( u3(i,j  ,k) + u3(i,j  ,kp1)
     4                 + u3(i,j-1,k) + u3(i,j-1,kp1) ) / 16.0
     5               * rhoi       * g2b i(i)
     6               * dg32ad2(j) * g31bi(i) * g32ai(j) * g32ai(j)
CJH  -- changed dg32bd2 to dg32ad2 in line above
           endif
          endif ! lgeom
#endif /* NOT ROTOR */
         else ! ldimen
          st = 0.0
         endif ! ldimen
          if(xptmass) then
c
c  4.  gravitational point mass
c
             if(lgeom .eq. 1) then
              r2        = ( x1b(i) - x1ptm )**2
     1                  + ( x2a(j) - x2ptm )**2
     2                  + ( x3b(k) - x3ptm )**2
              r2i       = ( x2a(j) - x2ptm )
     1                  / ( r2    * sqrt(r2   ) + tiny )
             endif ! CARTESIAN
             if(lgeom .eq. 2) then
              r2        = ( x1b(i) - x1ptm )**2
     1                  + ( x2a(j) - x2ptm )**2
     2                  + 2.0 * x2a(j) * x2ptm
     3                  * ( 1.0 - cos(x3b(k) - x3ptm) )
              r2i       = ( x2a(j) - x2ptm * cos(x3b(k) - x3ptm))
     1                  / ( r2    * sqrt(r2   ) + tiny )
             endif ! CYLINDRICAL
             if(lgeom .eq. 3) then
CJH
CJH  correction due to Sean Matt -- enforces condition that theta
CJH  component of force be zero if point mass located at origin.
CJH
              if(x1ptm .eq. 0.0) then
               r2i       = 0.0
              else ! x1ptm
               prpq      = x1b(i) * sin(x2a(j))
               pzpq      = x1b(i) * cos(x2a(j))
               prpt      = x1ptm  * sin(x2ptm)
               pzpt      = x1ptm  * cos(x2ptm)
c
               r2        = ( pzpq - pzpt )**2
     1                   + ( prpq - prpt )**2
     2                   + 2.0 * prpq * prpt
     3                   * ( 1.0 - cos(x3b(k) - x3ptm) )
               r2i       = (-( pzpq - pzpt ) * sin(x2a(j))
     1                   + ( prpq - prpt * cos(x3b(k) - x3ptm))
     2                   * cos(x2a(j))) / ( r2    * sqrt(r2   ) + tiny )
              endif ! x1ptm
             endif ! SPHERICAL
             st        = st - guniv * ptmass * r2i
          endif ! xptmass
c
c  5.  Radiation force 
c
          if(ldimen .ne. 1) then
          if(lrad .eq. 1) then
           derdx2 = ( er     (i  ,j  ,k  ) - er     (i  ,j-1,k  ) )
     1            * dx2bi(j) * g2bi(i)
           derdx1 = ( er     (i  ,j-1,k  ) - er     (i-1,j-1,k  ) )
     1            * dx1bi(i)
     2            + ( er     (i  ,j  ,k  ) - er     (i-1,j  ,k  ) )
     3            * dx1bi(i)
           derdx1 = derdx1 * 0.5
           if(ldimen .eq. 2) then
            derdx3 = 0.0
           else
            derdx3 = ( er     (i  ,j-1,k  ) - er     (i  ,j-1,km1) )
     1             * dx3bi(k)
     2             + ( er     (i  ,j  ,k  ) - er     (i  ,j  ,km1) )
     3             * dx3bi(k)
            derdx3 = derdx3 * 0.5 * g31bi(i) * g32ai(j)
           endif ! ldimen
           grad_er = sqrt( derdx1**2 + derdx2**2 + derdx3**2 )
           rdmcm   = kapr(i,j-1,k)     
           rdmcp   = kapr(i,j  ,k)
           ros_mfp = 2.0 / (rdmcm + rdmcp)
           big_R   = 2.0 * grad_er * ros_mfp
     1             / max(er(i,j  ,k)+er(i,j-1,k), tiny)
#ifdef MINERBO
c
c   Minerbo flux limiter
c
           if (big_R .le. 1.5) then
            flx_lim = 2.0/(3.0 + sqrt(9.0 + 12.0*big_R**2) )
           else
            flx_lim = 1.0/(1.0 + big_R + sqrt(1.0 + 2.0*big_R) )
           endif
#else /* MINERBO */
c
c   LP flux limiter
c
           flx_lim = (2.0 + big_R)/(6.0 + 3.0*big_R + big_R**2)
#endif /* MINERBO */
           st        = st - rhoi * flx_lim * derdx2
          endif ! lrad
c
c  6. Magnetic pressure (longitudinal Lorentz force)
c
          if(xmhd) then
           st = st - (d2b3oo(i) + d2b3op(i) + d2b1oo(i) + d2b1po(i))
     3        * rhoi * dx2bi(j) * g2bi(i) * 0.25
          endif ! xmhd
         endif ! ldimen
c
         w2(i,j,k) = u2(i,j,k) + dt * st
c
c.......................................................................
c
c Perform an explicit update for v3
c
         st = 0.0
         if(ldimen .ne. 1) then
          rhoi = 2.0 / ( d(i,j,km1) + d(i,j,k) )
c
c  1.  pressure gradient
c
          st = - rhoi
     1       * ( poo(i)    - pom(i)      ) * dx3bi(k)
     2       * g31bi(i) * g32bi(j)
c
c  2.  gravitational potential gradient
c
          if(xgrav .or. xgrvfft) then
           if(.not. xsphgrv) then
            st      = st 
     .             + ( gp(i,j,k) - gp(i,j,km1) ) * dx3bi(k)
     .             * g31bi(i) * g32bi(j)
           endif
          endif
         endif ! ldimen
c
c  4.  gravitational point mass
c
          if(xptmass) then
           if(lgeom .eq. 1) then
             r2        = ( x1b(i) - x1ptm )**2
     1                 + ( x2b(j) - x2ptm )**2
     2                 + ( x3a(k) - x3ptm )**2
             r2i       = ( x3a(k) - x3ptm )
     1                 / ( r2    * sqrt(r2   ) + tiny )
           endif ! CARTESIAN
           if(lgeom .eq. 2) then
             r2        = ( x1b(i) - x1ptm )**2
     1                 + ( x2b(j) - x2ptm )**2
     2                 + 2.0 * x2b(j) * x2ptm
     3                 * ( 1.0 - cos(x3a(k) - x3ptm) )
             r2i       = x2ptm * sin(x3a(k) - x3ptm)
     1                 / ( r2    * sqrt(r2   ) + tiny )
           endif ! CYLINDRICAL
           if(lgeom .eq. 3) then
CJH
CJH  correction due to Sean Matt -- enforces condition that phi
CJH  component of force be zero if point mass located at origin.
CJH
            if(x1ptm .eq. 0.0) then
             r2i       = 0.0
            else ! x1ptm
             prpq      = x1b(i) * sin(x2b(j))
             pzpq      = x1b(i) * cos(x2b(j))
             prpt      = x1ptm  * sin(x2ptm)
             pzpt      = x1ptm  * cos(x2ptm)
c
             r2        = ( pzpq - pzpt )**2
     1                 + ( prpq - prpt )**2
     2                 + 2.0 * prpq * prpt
     3                 * ( 1.0 - cos(x3a(k) - x3ptm) )
             r2i       = prpt * sin(x3a(k) - x3ptm)
     1                 / ( r2    * sqrt(r2   ) + tiny )
            endif ! x1ptm
           endif ! SPHERICAL
             st        = st - guniv * ptmass * r2i
          endif ! xptmass
c
c  5. Radiation forces
c
         if(ldimen .ne. 1) then
          if(lrad .eq. 1) then
           derdx3 = ( er     (i  ,j  ,k  ) - er     (i  ,j  ,km1) )
     1            * dx3bi(k) * g31bi(i) * g32bi(j)
           derdx1 = ( er     (i  ,j  ,km1) - er     (i-1,j  ,km1) )
     1            * dx1bi(i)
     2            + ( er     (i  ,j  ,k  ) - er     (i-1,j  ,k  ) )
     3            * dx1bi(i)
           derdx1 = derdx1 * 0.5
           derdx2 = ( er     (i  ,j  ,km1) - er     (i  ,j-1,km1) )
     1            * dx2bi(j)
     2            + ( er     (i  ,j  ,k  ) - er     (i  ,j-1,k  ) )
     3            * dx2bi(j)
           derdx2 = derdx2 * 0.5 * g2bi(i)  
           grad_er = sqrt( derdx1**2 + derdx2**2 + derdx3**2 )
           rdmcm   = kapr(i,j,km1)
           rdmcp   = kapr(i,j,k  )       
           ros_mfp = 2.0 / (rdmcm + rdmcp)
           big_R   = 2.0 * grad_er * ros_mfp
     1            / max(er(i,j,k  )+er(i,j,k-1), tiny)
#ifdef MINERBO
           if (big_R .le. 1.5) then
            flx_lim = 2.0/(3.0 + sqrt(9.0 + 12.0*big_R**2) )
           else
            flx_lim = 1.0/(1.0 + big_R + sqrt(1.0 + 2.0*big_R) )
           endif
#else /* MINERBO */
           flx_lim = (2.0 + big_R)/(6.0 + 3.0*big_R + big_R**2)
#endif /* MINERBO */
           st        = st - rhoi * flx_lim * derdx3
          endif
c
          if(xmhd) then
c
c  6. Magnetic pressure (longitudinal Lorentz force)
c
           st = st - (d3b1oo(i) + d3b1po(i) + d3b2oo(i) + d3b2op(i))
     3        * rhoi * dx3bi(k) * g31bi(i) * g32bi(j) * 0.25
          endif ! xmhd
c
          endif ! ldimen
          w3(i,j,k) = u3(i,j,k) + dt * st
20      continue
30     continue
40    continue
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                F O R C E S                \\\\\\\\\\
c
c=======================================================================
c
c
