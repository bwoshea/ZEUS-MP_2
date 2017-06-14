c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////             G C O L L A P S E             \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine gcollapse
c
c    gcollapse <----------- initialises spherical gas cloud for
c                           self gravitational collapse, november, 1999
c
c    written by: PSLi
c
c  PURPOSE: Sets up a spherical region at a specified point on
c  the grid (x10, x20, x30) with the specified radius (r) whose density
c  differ from the rest of the initial grid.
c
c  LOCAL VARIABLES:
c    r            initial radius of overpressured region
c    x10,x20,x30  coordinates of centre of overpressured region.
c    drat         ratio of density  across blast front
c    prat         ratio of pressure across blast front
c    d0           density          in ambient medium (default = 1.0)
c    p0           pressure         in ambient medium (default = 0.6)
c    e0           internal energy  in ambient medium (default = 0.9)
c    er0          radiation energy in ambient medium (default = 1.0)
c    v10          1-velocity       in ambient medium
c    v20          2-velocity       in ambient medium
c    v30          3-velocity       in ambient medium
c    b10          1-magnetic field on entire grid
c    b20          2-magnetic field on entire grid
c    b30          3-magnetic field on entire grid
c    d1           density          in central region (default = 1.0)
c    p1           pressure         in central region (default = 0.6)
c    e1           internal energy  in central region (default = 0.9)
c    er1          radiation energy in central region (default = 30.)
c    v11          1-velocity       in central region
c    v21          2-velocity       in central region
c    v31          3-velocity       in central region
c    m,drs,drc    parameters for specifying a sphere whose surface is
c                 sinusoidally perturbed (spherical coordinates only
c                 For an unperturbed sphere, set all values to zero
c                 (default).
c
c  EXTERNALS:
c    OVERLAP      [XYZ only]
c    BNDYALL
c    BSETMAG
c
c-----------------------------------------------------------------------
c
      use real_prec
      use config
      use param
      use field
      use bndry
      use grid
      use root
      use scratch
      use cons
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      integer  :: i, j, k, ip1, jp1, kp1, m
c
      real(rl) :: r      , x10, x20    , x30   , drat,
     .            prat   , d0 , p0     , e0    , v10,
     .            v20    , v30, b10    , b20   , b30,
     .            d1     , p1 , e1     , v11   , v21,
     .            v31    , drs, cofrac , mass  ,
     .            drc    , rsq, rin    , rout  , frac,
     .            er0    , er1, ros_mfp, dx_min, flx_lim,
     .            dmc_max
c
      integer  :: iin (ijkn), iout(ijkn), jin (ijkn),
     1            jout(ijkn), kin (ijkn), kout(ijkn)
c
      real(rl) :: massk(ijkn), sasum, overlap
c
      implicit NONE
c
      namelist / pgen     /
     .             r   , x10, x20     , x30     , drat,
     .             prat, d0 , p0      , e0 ,er0,
     .             v10 , v20, v30,
     .             b10 , b20, b30,
     .             d1  , p1 , e1 , er1,
     .             v11 , v21, v31,
     .             drs , drc, m
c
c-----------------------------------------------------------------------
c
       r    = 1.0
       x10  = 0.0
       x20  = 0.0
       x30  = 0.0
       drat = 0.0
       prat = 0.0
       d0   = 1.0
       p0   = 0.6
       e0   = 0.0
       v10  = 0.0
       v20  = 0.0
       v30  = 0.0
       b10  = 0.0
       b20  = 0.0
       b30  = 0.0
       d1   = 1.0
       p1   = 0.6
       e1   = 0.0
       v11  = 0.0
       v21  = 0.0
       v31  = 0.0
       drs  = 0.0
       drc  = 0.0
       m    = 0
       er0  = 0.0
       er1  = 0.0
c
       if (myid .eq. 0) then
         read (1, pgen)
         write (2, pgen)
#ifdef MPI_USED
         buf_in( 1) = r   
         buf_in( 2) = x10 
         buf_in( 3) = x20 
         buf_in( 4) = x30 
         buf_in( 5) = drat
         buf_in( 6) = prat
         buf_in( 7) = d0  
         buf_in( 8) = p0  
         buf_in( 9) = e0  
         buf_in(10) = v10 
         buf_in(11) = v20 
         buf_in(12) = v30 
         buf_in(13) = b10 
         buf_in(14) = b20 
         buf_in(15) = b30 
         buf_in(16) = d1  
         buf_in(17) = p1  
         buf_in(18) = e1  
         buf_in(19) = v11 
         buf_in(20) = v21 
         buf_in(21) = v31 
         buf_in(22) = drs 
         buf_in(23) = drc 
         buf_in(24) = er0
         buf_in(25) = er1
         ibuf_in( 1) = m   
#endif
       endif ! myid
#ifdef MPI_USED
        call MPI_BCAST( buf_in, 25, MPI_FLOAT
     &                , 0, comm3d, ierr )
        call MPI_BCAST( ibuf_in, 1, MPI_INTEGER
     &                , 0, comm3d, ierr )
        if (myid .ne. 0) then
         r    = buf_in( 1)
         x10  = buf_in( 2)
         x20  = buf_in( 3)
         x30  = buf_in( 4)
         drat = buf_in( 5)
         prat = buf_in( 6)
         d0   = buf_in( 7)
         p0   = buf_in( 8)
         e0   = buf_in( 9)
         v10  = buf_in(10)
         v20  = buf_in(11)
         v30  = buf_in(12)
         b10  = buf_in(13)
         b20  = buf_in(14)
         b30  = buf_in(15)
         d1   = buf_in(16)
         p1   = buf_in(17)
         e1   = buf_in(18)
         v11  = buf_in(19)
         v21  = buf_in(20)
         v31  = buf_in(21)
         drs  = buf_in(22)
         drc  = buf_in(23)
         er0  = buf_in(24)
         er1  = buf_in(25)
         m    = ibuf_in( 1)
        endif ! myid
#endif
c
c      Set up atmosphere.
c
       if(e0 .ne. 0.0) p0 = e0 * gamm1
       if(e0 .eq. 0.0) e0 = p0 / gamm1
       if(lrad .ne. 0) er0 = rad_con * (gamm1*e0*mmw*mh/(d0*boltz))**4
C       print *, myid_w,x1a(is),x1a(ie+1)
       do 30 k=ks,ke
         do 20 j=js,je
           do 10 i=is,ie
             d (i,j,k) = d0
             v1(i,j,k) = v10
             v2(i,j,k) = v20
             v3(i,j,k) = v30
             if(xiso .eqv. .false.) e (i,j,k) = e0
             if(lrad .ne. 0) er(i,j,k) = er0
10         continue
20       continue
30     continue
      if(xmhd) then
       do 31 k=ks-2,ke+3
         do 21 j=js-2,je+3
           do 11 i=is-2,ie+3
             b1(i,j,k) = b10
             b2(i,j,k) = b20
             b3(i,j,k) = b30
11         continue
21       continue
31     continue
      endif ! xmhd
c
c      Set up central region.
c
       do 40 i=is,ie
         ip1 = i + 1
         if ( abs(x1a(i)-x10) .lt. abs(x1a(ip1)-x10) ) then
           iin (i) = i
           iout(i) = ip1
         else
           iin (i) = ip1
           iout(i) = i
         endif
40     continue
c
       do 50 j=js,je
         jp1 = j + 1
         if ( abs(x2a(j)-x20) .lt. abs(x2a(jp1)-x20) ) then
           jin (j) = j
           jout(j) = jp1
         else
           jin (j) = jp1
           jout(j) = j
         endif
50     continue
c
       do 60 k=ks,ke
         kp1 = k + 1
         if ( abs(x3a(k)-x30) .lt. abs(x3a(kp1)-x30) ) then
           kin (k) = k
           kout(k) = kp1
         else
           kin (k) = kp1
           kout(k) = k
         endif
         massk(k) = 0.0
60     continue
c
       if (drat .ne. 0.0) d1 = d0 * drat
       if (prat .ne. 0.0) p1 = p0 * prat
       if (e1   .ne. 0.0) then
         p1 = e1 * gamm1
       else
         e1 = p1 / gamm1
       endif
       if(lrad .ne. 0) er1 = rad_con * (gamm1*e1*mmw*mh/(d1*boltz))**4
c
       do 90 k=ks,ke
         do 80 j=js,je
           do 70 i=is,ie
            if(lgeom .eq. 1) then
             rsq  = r**2
             rin  = ( x1a(iin (i)) - x10 )**2
     1            + ( x2a(jin (j)) - x20 )**2
     2            + ( x3a(kin (k)) - x30 )**2
             rout = ( x1a(iout(i)) - x10 )**2
     1            + ( x2a(jout(j)) - x20 )**2
     2            + ( x3a(kout(k)) - x30 )**2
            endif ! CARTESIAN
            if(lgeom .eq. 2) then
             rsq  = r**2
             rin  = ( x1a(iin (i)) - x10 )**2
     1            + ( x2a(jin (j)) - x20 )**2
             rout = ( x1a(iout(i)) - x10 )**2
     1            + ( x2a(jout(j)) - x20 )**2
            endif ! CYLINDRICAL
            if(lgeom .eq. 3) then
             rsq  = r**2 * ( 1.0 + drs * sin (m * x2a(j))
     1                           + drc * cos (m * x2a(j)) )**2
             rin  = ( x1a(iin (i)) - x10 )**2
             rout = ( x1a(iout(i)) - x10 )**2
            endif ! SPHERICAL
             if ( (rin .lt. rsq) .and. (rout .le. rsq) ) then
               d (i,j,k) = d1
               v1(i,j,k) = v11
               v2(i,j,k) = v21
               v3(i,j,k) = v31
               if(xiso .eqv. .false.) e (i,j,k) = e1
               if(lrad .ne. 0) er(i,j,k) = er1
               massk(k) = massk(k) + d1 * dvl1a(i) * dvl2a(j) * dvl3a(k)
             endif
             if ( (rin .lt. rsq) .and. (rout .gt. rsq) ) then
              if(lgeom .eq. 1) then
               frac     = overlap ( 1, r, x10, x20, x30
     1                      , x1a(iin (i)), x2a(jin (j)), x3a(kin (k))
     2                      , x1a(iout(i)), x2a(jout(j)), x3a(kout(k)) )
              else ! lgeom
               frac     = ( rsq - rin ) / ( rout - rin )
              endif ! lgeom
               cofrac   = 1.0 - frac
               d(i,j,k) = d1 * frac + d0 * cofrac
               if(xiso .eqv. .false.) e(i,j,k) = e1 * frac + e0 * cofrac
               if(lrad .ne. 0) er(i,j,k) = er1 * frac + er0 * cofrac
               massk(k) = massk(k)
     1                  + d1 * frac * dvl1a(i) * dvl2a(j) * dvl3a(k)
             endif
70         continue
80       continue
90     continue
       mass = SASUM ( nx3z, massk(ks), 1 )
C       write(*,"('total mass = ',1pd12.4)")mass
c
       call gravity
c
       return
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                 B L A S T                 \\\\\\\\\\
c
c=======================================================================
c
c
c=======================================================================
c
c    \\\\\\\\\\        B E G I N   F U N C T I O N        //////////
c    //////////               O V E R L A P               \\\\\\\\\\
c
c=======================================================================
c
       real(rl) function overlap ( ishp, rad, x0, y0, z0, xin, yin
     1                       , zin, xout, yout, zout )
c
c    dac:zeus3d.overlap <--------- overlap of region over Cartesian zone
c                                                            april, 1990
c
c    written by: David Clarke
c    modified 1:
c
c  PURPOSE:  Determines the fraction of a Cartesian zone that overlaps
c  the specified geometrical region (sphere or right cylinder).  This
c  is done by dividing the zone into 20**3 "subzones", and finding the
c  fraction of subzone centres lying inside the surface of the region.
c
c  INPUT VARIABLES:
c
c    ishp            =1 => sphere
c                    =2 => right cylinder
c    rad             radius of region
c    x0,y0,z0        coordinates of centre of curvature.
c    xin,yin,zin     coordinates of zone corner known to lie inside
c                    region.
c    xout,yout,zout  coordinates of zone corner diametrically opposed to
c                    zone corner known to lie inside region.
c
c  OUTPUT VARIABLES:
c
c  LOCAL VARIABLES:
c
c  EXTERNALS: [NONE]
c
c-----------------------------------------------------------------------
c
      use real_prec
      use config
      use param
c
      integer  :: i, j, k, nx, ny, nz, ishp
c
      real(rl) :: rad   , x0  , y0  , z0  , xin ,
     .            yin   , zin , xout, yout, zout,
     .            delx  , dely, delz, r   , fact,
     .            scount
c
      real(rl) :: xsq(20), ysq(20), zsq(20), count(20)
c
      implicit NONE
c-----------------------------------------------------------------------
c
      if(lgeom .eq. 1) then
c      Number of subzones in the x-direction is "nx", etc. for "ny" and
c  "nz".  Increment between subzones in x-direction is "delx", etc. for
c  "dely" and "delz".
c
       nx   = 20
       ny   = 20
       nz   = 20
       delx = ( xout - xin ) /  real( nx )
       dely = ( yout - yin ) /  real( ny )
       delz = ( zout - zin ) /  real( nz )
c
c      Set up subgrid inside zone.
c
       do 10 i=1,nx
         xsq  (i) = ( xin + ( 0.5 +  real(i-1) ) * delx - x0 )**2
         count(i) = 0.0
10     continue
       do 20 j=1,ny
         ysq  (j) = ( yin + ( 0.5 +  real(j-1) ) * dely - y0 )**2
20     continue
       do 30 k=1,nz
         zsq  (k) = ( zin + ( 0.5 +  real(k-1) ) * delz - z0 )**2
30     continue
c
c      Count the number of subzones lying inside the surface of the
c  region which passes through the zone.
c
       fact  = 1.0
       if (ishp .eq. 2) fact = 0.0
       do 60 k=1,nz
         do 50 j=1,ny
           do 40 i=1,nx
             r = sqrt ( fact * xsq(i) + ysq(j) + zsq(k) )
             if (r .le. rad) count(i) = count(i) + 1.0
40         continue
50       continue
60     continue
       scount = 0.0
       do 70 i=1,nx
         scount = scount + count(i)
70     continue
       scount =   max ( one, scount )
c
c      Set the fraction of the zone which overlaps the region.
c
       overlap = scount /  real( nx * ny * nz )
      else ! lgeom
       overlap = 1.0
      endif ! lgeom
c
       return
       end
c
c=======================================================================
c
c    \\\\\\\\\\          E N D   F U N C T I O N          //////////
c    //////////               O V E R L A P               \\\\\\\\\\
c
c=======================================================================
c
c
