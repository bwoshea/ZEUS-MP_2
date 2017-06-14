c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                    C T                    \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine ct_2D
c
c    dac:zeus3d.ct <-------- updates B-field using constrained transport
c    from jms:zeus2d.ct                                    october, 1989
c
c    written by: David Clarke
c    modified 1: May, 1990 by David Clarke; reworked to call the new
c                interpolation routines which need only be called once.
c    modified 2: August, 1990 by David Clarke; moved magnetic fields to
c                the face-centres (cospatial with the velocities).
c                Implemented a method of characteristics for evaluating
c                the emf's.  The transverse Lorentz accelerations are
c                now applied to the velocities during this step.
c    modified 3: November, 1990 by David Clarke; resurrected non-MoC
c                algorithm for evaluating emfs.  Added EDITOR alias MOC
c                which needs to be defined if MoC is to be used.
c    modified 4: June, 1992 by David Clarke; reworked singularity
c                formalism.
c    modified 5: December, 1992 by David Clarke (as suggested by John
c                Hawley); split velocity update from emf computation.
c                Velocities are now Lorentz-accelerated with *old*
c                magnetic field values. emf's are then estimated with
c                Lorentz-updated velocities.
c    modified 6: 3 March 1998, by Mordecai-Mark Mac Low; translated 
c                into ZEUS-MP form.
c    modified 7: John Hayes, October 2005; created ct_2D clone which
c                assumes symmetry about the K axis; after ZEUS3D.
c
c    modified 8: John Hayes, May 15, 2006; changed the call to LORENTZ_D
c                so that old velocities are passed in via a locally
c                defined scratch array rather than through the v[123] field
c                arrays; the latter now receive updated values only. This
c                change, combined with additional edits in LORENTZ_D
c                itself, correct an old error in which two variables
c                were assigned the same space in memory.
c
c    modified 9: extended ranges of [i,j,k] DO loops copying v[1,2,3]
c                values into u[1,2,3] arrays.  Avoids passing uninitialized
c                u[1,2,3] values into bvalv[1,2,3] -- JHayes, 09/20/2006
c
c  PURPOSE:  This routine transports the three components of the
c  magnetic field using a variation of the non-relativistic Constrained
c  Transport scheme (CT), developed by Chuck Evans and John Hawley (Ap.
c  J., 342, 700).  In this implementation, the magnetic field components
c  are face-centred, cospatial with the velocities and are updated using
c  edge-centred emf's which are cospatial with the current densities.
c
c  The emf's are evaluated by HSMOC in which the velocities and
c  magnetic field components required to compute the emf's are estimated
c  using the Method of Characteristics (MoC) algorithm developed by Jim
c  Stone et al. for 2-D.  For self-consistency, the transverse Lorentz
c  accelerations have been removed from STV* and are applied to the
c  velocities in LORENTZ using the pre-updated magnetic fields.  By
c  experimentation, it has been determined that performing the Lorentz
c  update after the magnetic field update is unstable.
c
c  LOCAL VARIABLES:
c    emf1      emf along the 2-3 edges of the grid (= v2*b3 - v3*b2)
c    emf2      emf along the 3-1 edges of the grid (= v3*b1 - v1*b3)
c    emf3      emf along the 1-2 edges of the grid (= v1*b2 - v2*b1)
c
c BOUNDARY VALUES USED:
c
c    var    ii    oi    ij    oj    ik    ok
c    ----  ----  ----  ----  ----  ----  ----
c    emf1  is-2  ie+2  js-2  je+3  ks-2  ke+3
c    emf2  is-2  ie+3  js-2  je+2  ks-2  ke+3
c    emf3  is-2  ie+3  js-2  je+3  ks-2  ke+2
c
c  EXTERNALS:
c    LORENTZ_D, HSMOC 
c
c-----------------------------------------------------------------------
c
      use config
      use param
      use grid
      use field
      use root
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
      integer  :: i, ip1, j, jp1, k, kp1
      integer  :: kone, km1   !asif
c
      real(rl) :: qty1(ijkn), qty1ni(ijkn), qty2(ijkn),
     .            qty2ni  (ijkn)
c
      real(rl) :: emf1(in,jn,kn), emf2(in,jn,kn),
     .            emf3(in,jn,kn)
c
      real(rl) :: u1(in,jn,kn), u2(in,jn,kn), u3(in,jn,kn)
c
c-----------------------------------------------------------------------
c	asif 
       if (ldimen .eq. 3) then
       kone=1
       else
       kone=0
       endif

      if(xforce) then
c
c-----------------------------------------------------------------------
c-------------------------> Update velocities <-------------------------
c-----------------------------------------------------------------------
c
c      Compute the transverse Lorentz forces and accelerate the
c  velocities.
c
       do k = 1, kn
        do j = 1, jn
         do i = 1, in
          u1(i,j,k) = v1(i,j,k)
          u2(i,j,k) = v2(i,j,k)
          u3(i,j,k) = v3(i,j,k)
         enddo
        enddo
       enddo
       call lorentz_d (u1, u2, u3, v1, v2, v3)
      endif ! xforce
c
c-----------------------------------------------------------------------
c-------------------------------> emfs <--------------------------------
c-----------------------------------------------------------------------
c
c  update all the boundary values of v before going into hsmoc to
c  compute the emfs.  Boy could this ever benefit from some
c  overlapping! think about this... (M-MML/MLN 25.3.98)

       nreq=0
       nsub = nsub + 1
       call bvalv1(3,3,0,0,0,0,v1)
       call bvalv2(3,3,0,0,0,0,v2)
       call bvalv3(3,3,0,0,0,0,v3)

#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
       nreq=0
       nsub = nsub + 1
       call bvalv1(0,0,3,3,0,0,v1)
       call bvalv2(0,0,3,3,0,0,v2)
       call bvalv3(0,0,3,3,0,0,v3)

#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
       if(ldimen.eq.3) then   !asif
       nreq=0
       nsub = nsub + 1
       call bvalv1(0,0,0,0,3,3,v1)
       call bvalv2(0,0,0,0,3,3,v2)
       call bvalv3(0,0,0,0,3,3,v3)

#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
       endif   !asif, ldimen
c
       call hsmoc_2D( emf1, emf2, emf3 )
c
c-----------------------------------------------------------------------
c
c      The emf's are finished.  Update b1, b2, and b3 using the emf's.
c  Since the same emf's are used throughout the grid, div(b) will be
c  conserved numerically to within truncation error.
c
c      Coordinate-imposed boundary conditions (e.g., reflecting at x2a=0
c  in ZRP, periodic at x3a=2*pi in RTP and ZRP) are communicated to the
c  magnetic field by the emfs.  Both the old and new zone face areas are
c  used to account for grid compression.
c
c-----------------------------------------------------------------------
c-----------------------------> Update b1 <-----------------------------
c-----------------------------------------------------------------------
c
       do 10 i=is,ie+1
         qty1  (i) = g2a   (i) * g31a  (i)
        if(xvgrid) then
         qty1ni(i) = g2ani (i) * g31ani(i)
        else
         qty1ni(i) = g2ai (i) * g31ai(i)
        endif
10     continue
       do 20 j=js-2,je+2
         qty2  (j) = g32b  (j) * dx2a  (j)
        if(xvgrid) then
         qty2ni(j) = g32bni(j) * dx2ani(j)
        else
         qty2ni(j) = g32bi(j) * dx2ai(j)
        endif
20     continue

       k   = ks
       km1 = k-kone
       kp1 = k+kone
         do 40 j=js-2,je+2
           do 30 i=is,ie+1
            if(xvgrid) then
             b1(i,j,k) = ( b1(i,j,k) * qty1(i) * qty2(j) * dx3a(k)
     1                   + dt * ( emf3(i,j+1,k  ) - emf3(i,j,k)
     2                          - emf2(i,j  ,kp1) + emf2(i,j,k) ) )
     3                 * qty1ni(i) * qty2ni(j) * dx3ani(k)
            else
             b1(i,j,k) = ( b1(i,j,k) * qty1(i) * qty2(j) * dx3a(k)
     1                   + dt * ( emf3(i,j+1,k  ) - emf3(i,j,k)
     2                          - emf2(i,j  ,kp1) + emf2(i,j,k) ) )
     3                 * qty1ni(i) * qty2ni(j) * dx3ai(k)
            endif
30         continue
40       continue
50     continue
c
c-----------------------------------------------------------------------
c-----------------------------> Update b2 <-----------------------------
c-----------------------------------------------------------------------
c
       do 60 i=is-2,ie+2
         qty1  (i) = g31b  (i) * dx1a  (i)
        if(xvgrid) then
         qty1ni(i) = g31bni(i) * dx1ani(i)
        else
         qty1ni(i) = g31bi(i) * dx1ai(i)
        endif
60     continue


         do 90 j=js,je+1
           do 80 i=is-2,ie+2
            if(xvgrid) then
             b2(i,j,k) = ( b2(i,j,k) * qty1(i) * g32a(j) * dx3a(k)
     1                   + dt * ( emf1(i  ,j,kp1) - emf1(i,j,k)
     2                          - emf3(i+1,j,k  ) + emf3(i,j,k) ) )
     3                 * qty1ni(i) * g32ani(j) * dx3ani(k)
            else
             b2(i,j,k) = ( b2(i,j,k) * qty1(i) * g32a(j) * dx3a(k)
     1                   + dt * ( emf1(i  ,j,kp1) - emf1(i,j,k)
     2                          - emf3(i+1,j,k  ) + emf3(i,j,k) ) )
     3                 * qty1ni(i) * g32ai(j) * dx3ai(k)
            endif
80         continue
90       continue
100    continue
c
c-----------------------------------------------------------------------
c-----------------------------> Update b3 <-----------------------------
c-----------------------------------------------------------------------
c
       do 110 i=is-2,ie+2
         qty1  (i) = g2b   (i) * dx1a  (i)
        if(xvgrid) then
         qty1ni(i) = g2bni (i) * dx1ani(i)
        else
         qty1ni(i) = g2bi (i) * dx1ai(i)
        endif
110    continue

         do 140 j=js-2,je+2
           do 130 i=is-2,ie+2
            if(xvgrid) then
             b3(i,j,k) = ( b3(i,j,k) * qty1(i) * dx2a(j)
     1                   + dt * ( emf2(i+1,j  ,k) - emf2(i,j,k)
     2                          - emf1(i  ,j+1,k) + emf1(i,j,k) ) )
     3                 * qty1ni(i) * dx2ani(j)
            else
             b3(i,j,k) = ( b3(i,j,k) * qty1(i) * dx2a(j)
     1                   + dt * ( emf2(i+1,j  ,k) - emf2(i,j,k)
     2                          - emf1(i  ,j+1,k) + emf1(i,j,k) ) )
     3                 * qty1ni(i) * dx2ai(j)
            endif
130        continue
140      continue
150    continue
c
       return
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                    C T                    \\\\\\\\\\
c
c=======================================================================
c
