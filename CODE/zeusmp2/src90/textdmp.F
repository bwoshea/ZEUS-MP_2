c=======================================================================
c/////////////////////////  EXAMPLE USERDUMP  \\\\\\\\\\\\\\\\\\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
      subroutine textdmp
c
c  written by:   Robert Fiedler
c
c  PURPOSE:  Example USERDUMP routine to generate ASCII output
c            files during a ZEUS-MP run.
c
c  USAGE:  To use this routine, include "-u textdmp" on the
c          Make_zeusmp command line.
c
c  DESCRIPTION:  This example USERDUMP routine prints out the
c                coordinates and the density, energy density,
c                and velocity components for a 3-D
c                hydrodynamics simulation.  It must be modified
c                to list magnetic field components, 
c                gravitational potential, radiation energy density, 
c                etc.  
c
c-----------------------------------------------------------------------
      use real_prec
      use param
      use config
      use field
      use grid
      use bndry
      use root
      use cons
      use scratch
      use opac_law
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer :: i, j, k
      real(rl) :: etotal, ekin, eint, etot_glb, dvb, totdvb,
     .            glb_j, jtot, glb_ek(3),
     .            ek(3), v1scr(in), v2scr(in), b2scr(in), b3scr(in),
     .            v3scr(in), b1scr(in), dscr(in), pscr(in)
c
      real(rl) :: so_eps , tau    , eriibn, timarg
c
      common /soui/ so_eps, tau, eriibn
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c Write out the coordinates and the density for each zone in the
c computational domain to unit 12 (usrfile).
c
12    format(1p5e16.8)
13    format(1p6e16.8)
14    format(1p7e16.8)
15    format(1p8e16.8)
      open(12,file=usrfile)
c
c-----------------------------------------------------------------------
c     if MHD, then check to verify that div(B) = 0 everywhere
c-----------------------------------------------------------------------
c
      if(xmhd) then
       call diverg ( b1, b2, b3, 1, 1, w3da, dvb )
#ifdef MPI_USED
       call MPI_Reduce(dvb, totdvb, 1,
     &                 MPI_FLOAT, MPI_SUM, 0, comm3d, ierr)
#else
       totdvb = dvb
#endif /* MPI_USED */
      endif ! xmhd
c
      if(coords(1) .eq. 0) then
       write(12,"('ZEUS-MP field variables at time = ',1pe15.8)") time
       write(12,"(a)")
      endif
c
      if(lrad .eq. 0) then
       if(xmhd) then
        write(12,"('Total div(B) = ',1pd12.4/)")totdvb
       endif ! xmhd
#ifndef ROTOR
#ifndef MHDSHKTUBE
#ifndef VORTEX
       write(12,"('         x1b             x2b             x3b',
     1            '         density         e               v1 ')")
#else /* VORTEX */
       write(12,"('      x1b             x2b             x3b',
     1            '           density           Pgas            Pmag',
     .            '           div(B)     ')")
#endif /* VORTEX */
#endif /* NOT MHDSHKTUBE */
       if(xtotnrg .eqv. .false.) then
#ifndef MHDSHKTUBE
#ifndef VORTEX
        write(12,13) (((x1b(i), x2b(j), x3b(k),
     .                  d(i,j,k), 
     .                  e(i,j,k),
C     .                  gamm1*mmw*mh*e(i,j,k)/(boltz*d(i,j,k)),
     .                  v1(i,j,k),
     .                  i=is,ie), j=js,je ), k=ks,ke)
#else /* VORTEX */
        write(12,14) (((x1b(i), x2b(j), x3b(k),
     .                  d(i,j,k), 
     .                  gamm1*e(i,j,k),
     .                  0.125*(  (b1(i,j,k) + b1(i+1,j  ,k  ))**2
     .                         + (b2(i,j,k) + b2(i  ,j+1,k  ))**2
     .                         + (b3(i,j,k) + b3(i  ,j  ,k  ))**2  ),
     .                  w3da(i,j,k),
     .                  i=is,ie), j=js,je ), k=ks,ke)
#endif /* VORTEX */
#else /* MHDSHKTUBE */
       write(12,"('         x1b           density            P ',
     1            '            V1           V2              V3 ',
     .            '            B2           B3')")
       do i = is, ie
        d scr(i) = d (i,js,ks)
        p scr(i) = p (i,js,ks)
        v1scr(i) = v1(i,js,ks)
        v2scr(i) = v2(i,js,ks)
        v3scr(i) = v3(i,js,ks)
        b1scr(i) = b1(i,js,ks)
        b2scr(i) = b2(i,js,ks)
        b3scr(i) = b3(i,js,ks)
        if(abs(v1(i,js,ks)) .lt. 1.0d-10) v1scr(i) = 0.0
        if(abs(v2(i,js,ks)) .lt. 1.0d-10) v2scr(i) = 0.0
        if(abs(v3(i,js,ks)) .lt. 1.0d-10) v3scr(i) = 0.0
        if(abs(b1(i,js,ks)) .lt. 1.0d-10) b1scr(i) = 0.0
        if(abs(b2(i,js,ks)) .lt. 1.0d-10) b2scr(i) = 0.0
        if(abs(b3(i,js,ks)) .lt. 1.0d-10) b3scr(i) = 0.0
       enddo
        write(12,15) (x1b(i), 
     .                  dscr(i), pscr(i),
     .                  v1scr(i), 
     .                  v2scr(i),
     .                  v3scr(i),
     .                  b2scr(i),
     .                  b3scr(i),
     .                  i=is,ie)
#endif /* MHDSHKTUBE */
       else
        write(12,13) (((x1b(i), x2b(j), x3b(k),
     .                  d(i,j,k), 
     .                  e(i,j,k) - 0.125D0*d(i,j,k)*(
     .                      (v1(i,j,k)+v1(i+1,j,k))**2 
     .                    + (v2(i,j,k)+v2(i,j+1,k))**2 
     .                    + (v3(i,j,k)+v3(i,j,k+1))**2        
     .                                               ),
     .                  v1(i,j,k),
     .                  i=is,ie), j=js,je), k=ks,ke)
       endif
#endif /* NOT ROTOR */
      else ! lrad
#ifdef RSHOCK
       write(12,"(' x1b (comoving)          x2b             x3b',
     1            '         density         Tgas            Trad ')")
       write(12,13) ((( 
     .                 (x1b(i)-fois(3)*time)/1.D10,
     .                 x2b(j), x3b(k), d (i,j,k), 
     .                 gamm1*mmw*mh*e(i,j,k)/(boltz*d(i,j,k)),
     .                 sqrt(sqrt(er(i,j,k)/rad_con)),
     .                 i=is,ie), j=js,je ), k=ks,ke)
#else
#ifndef MARSHAK
       write(12,"('         x1b             x2b             x3b',
     1            '         density         Tgas            Trad ')")
       write(12,13) (((x1b(i), x2b(j), x3b(k), d(i,j,k),
     .                 gamm1*mmw*mh*e(i,j,k)/(everg*d(i,j,k)),
     .                 boltz*sqrt(sqrt(er(i,j,k)/rad_con))/everg,
     .                 i=is,ie), j=js,je ), k=ks,ke)
#else /* MARSHAK */
       write(12,"('         x1b             x2b             x3b',
     1            '         density         egas            erad ')")
       write(12,13) (((x1b(i), 
     .                 x1b(i)/rmfp0,
     .                 sqrt(3.0D0)*x1b(i)/rmfp0,
     .                 d(i,j,k),
     .          so_eps*e (i,j,k),
     .                 er(i,j,k),
     .                 i=is,ie), j=js,je ), k=ks,ke)
#endif /* MARSHAK */
#endif /* RSHOCK */
      endif ! lrad
c
22    format(1p4e16.8)
23    format(1p5e16.8)
#ifndef MHDSHKTUBE
#ifndef VORTEX
#ifndef ROTOR
      write(12,"(a)")
      write(12,"('                     x1b             x2a        ',
     1 '     x3b          v2')")
      write(12,32) ( ( (i,j,k,x1b(i), x2a(j), x3b(k), v2(i,j,k),
     1                 i=is,ie), j=js,je ), k=ks,ke )
      write(12,"(a)")
      write(12,"('         x1b             x2b             x3a',
     1 '          v3')")
      write(12,42) ( ( ( x1b(i), x2b(j), x3a(k), 
     .                   v3(i,j,k),
     1                 i=is,ie), j=js,je ), k=ks,ke )
#else
      write(12,"('         x1b             x2b             v3 ',
     1         '          Omega             B3')")
      write(12,52) ( ( ( x1b(i), x2b(j), v3(i,j,k),
     .                   v3(i,j,k)*x2bi(j),b3(i,j,k),
     .                   i=is,ie), j=je,je ), k=ks,ke )
#endif /* ROTOR */
#endif /* VORTEX */
#endif /* MHDSHKTUBE */
c
32    format(3i4,1p4e16.8)
42    format(1p4e16.8)
52    format(1p5e16.8)
c
      close(12)
c
      return
      end
