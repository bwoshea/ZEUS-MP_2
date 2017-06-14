c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              T R A N S P R T              \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine transprt
c
c    mln:zeus3d.transprt <------------------ controls the transport step
c                                                          october, 1987
c
c    written by: Mike Norman
c    modified 1: June, 1988 by Jim Stone; incorporated into ZEUS2D
c    modified 2: February, 1990 by David Clarke; incorporated into
c                ZEUS3D
c    modified 3: Feb. 15, 1996 by Robert Fiedler; completely
c                rewritten for ZEUS-MP.
c    modified 4: Dec. 20, 1996 by Robert Fiedler; added radiation.
c    modified 5: Jan. 21, 1997 by Robert Fiedler; NO_TRANSPORT switch
c    modified 6: Dec. 30, 1999 by PSLi; added update of momenta.
c
c  PURPOSE: This subroutine transports the field variables through the
c  mesh in a directionally split manner.  In each succesive call to
c  TRANSPRT, the order of the directions is permuted (resulting in
c  XYZ...YXZ...YZX...ZYX...ZXY...XZY...XYZ...etc.).  This MAY be better
c  than leaving the order the same each time (XYZ...XYZ...etc), and
c  seems to be better than unsplit schemes (Hawley).  Momenta are
c  computed from velocities in "avisc" and then transported.  Velocities
c  are not updated until the end of the transport step.  
c
c  The magnetic field components are updated by CT which is a merger (as
c  implemented by Jim Stone) of the method of characteristics and a
c  variant of the constrained transport algorithm of Evans and Hawley.
c
c  Note that the order in which variables are transported is important
c  (especially d).  
c
c  LOCAL VARIABLES:
c
c  EXTERNALS:
c    CT
c    ADVX1   , ADVX2   , ADVX3
c
c-----------------------------------------------------------------------
c
      use config
      use param
      use root
      use grid
      use field
      use bndry
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
      integer i,j,k
c
c-----------------------------------------------------------------------
      if(xhydro .eqv. .false.) go to 666
c
      if(xmhd) then
c
c      Transport the three components of B using Constrained Transport.
c
       call ct
      endif ! xmhd
c
c      Momentum densities were computed from velocities in the
c      artificial viscosity substep (which must therefore not be
c      skipped, even if qcon = 0.)  Momentum density boundary
c      values are not needed.
c
CPS 
      DO k=ks,ke
         DO j=js,je
           DO i=is,ie
             w3da(i,j,k) = v1(i,j,k) * 0.5 * (d(i-1,j  ,k  ) + d(i,j,k))
             w3db(i,j,k) = v2(i,j,k) * 0.5 * (d(i  ,j-1,k  ) + d(i,j,k))
     1                     * g2b(i)
             w3dc(i,j,k) = v3(i,j,k) * 0.5 * (d(i  ,j  ,k-1) + d(i,j,k))
     1                     * g31b(i) * g32b(j)
           ENDDO
         ENDDO
       ENDDO
CX       if(nhy .eq. 21) then
CX        write(*,"('TRANSPRT: w3db = ',1p2d16.8)")w3db(4,4,3),w3db(4,4,4)
CX       endif
C
c---------------- directional split in X1-X2-X3 fashion ----------------
c
       nseq = 0  ! in /root/
       if (ix1x2x3 .eq. 1) then
c
c       subroutine advx1 (dlo,den
c     &                  ,eod,edn
c     &                  ,ero,ern
c     &                  ,abo,abn
c     &                  ,mflx,s1,s2,s3)
c        
         call advx1 (w3dd,d   
     &              ,w3de,w3dg
     &              ,er  ,w3dh
     .              ,abun,w4da
     &              ,w3df,w3da,w3db,w3dc)
         call advx2 (d   ,w3dd
     &              ,w3dg,w3de
     &              ,w3dh,er
     .              ,w4da,abun
     &              ,w3df,w3da,w3db,w3dc)
         call advx3 (w3dd,d
     &              ,w3de,w3dg
     &              ,er  ,w3dh
     .              ,abun,w4da
     &              ,w3df,w3da,w3db,w3dc)
c
         ix1x2x3 = 2
         goto 10
c
c---------------- directional split in X2-X1-X3 fashion ----------------
c
       else if (ix1x2x3 .eq. 2) then
c
         call advx2 (w3dd,d   
     &              ,w3de,w3dg
     &              ,er  ,w3dh
     .              ,abun,w4da
     &              ,w3df,w3da,w3db,w3dc)
         call advx1 (d   ,w3dd
     &              ,w3dg,w3de
     &              ,w3dh,er
     .              ,w4da,abun
     &              ,w3df,w3da,w3db,w3dc)
         call advx3 (w3dd,d
     &              ,w3de,w3dg
     &              ,er  ,w3dh
     .              ,abun,w4da
     &              ,w3df,w3da,w3db,w3dc)
c
         ix1x2x3 = 3
         goto 10
c
c---------------- directional split in X2-X3-X1 fashion ----------------
c
       else if (ix1x2x3 .eq. 3) then
c
c       subroutine advx1 (dlo,den
c     &                  ,eod,edn
c     &                  ,mflx,s1,s2,s3)
c
         call advx2 (w3dd,d   
     &              ,w3de,w3dg
     &              ,er  ,w3dh
     .              ,abun,w4da
     &              ,w3df,w3da,w3db,w3dc)
         call advx3 (d   ,w3dd
     &              ,w3dg,w3de
     &              ,w3dh,er
     .              ,w4da,abun
     &              ,w3df,w3da,w3db,w3dc)
         call advx1 (w3dd,d
     &              ,w3de,w3dg
     &              ,er  ,w3dh
     .              ,abun,w4da
     &              ,w3df,w3da,w3db,w3dc)
c
         ix1x2x3 = 4
         goto 10
c
c---------------- directional split in X3-X2-X1 fashion ----------------
c
       else if (ix1x2x3 .eq. 4) then
c
c       subroutine advx1 (dlo,den
c     &                  ,eod,edn
c     &                  ,mflx,s1,s2,s3)
c
         call advx3 (w3dd,d   
     &              ,w3de,w3dg
     &              ,er  ,w3dh
     .              ,abun,w4da
     &              ,w3df,w3da,w3db,w3dc)
CX         if(nhy .eq. 21)
CX     .    write(*,"('ADVX3: w3db = ',1p2d16.8)")w3db(4,4,3),w3db(4,4,4)
         call advx2 (d   ,w3dd
     &              ,w3dg,w3de
     &              ,w3dh,er
     .              ,w4da,abun
     &              ,w3df,w3da,w3db,w3dc)
CX         if(nhy .eq. 21)
CX     .    write(*,"('ADVX2: w3db = ',12pd16.8)")w3db(4,4,3),w3db(4,4,4)
         call advx1 (w3dd,d
     &              ,w3de,w3dg
     &              ,er  ,w3dh
     .              ,abun,w4da
     &              ,w3df,w3da,w3db,w3dc)
CX         if(nhy .eq. 21)
CX     .    write(*,"('ADVX1: w3db = ',1p2d16.8)")w3db(4,4,3),w3db(4,4,4)
c
         ix1x2x3 = 5
         goto 10
c
c---------------- directional split in X3-X1-X2 fashion ----------------
c
       else if (ix1x2x3 .eq. 5) then
c
c       subroutine advx1 (dlo,den
c     &                  ,eod,edn
c     &                  ,mflx,s1,s2,s3)
c
         call advx3 (w3dd,d   
     &              ,w3de,w3dg
     &              ,er  ,w3dh
     .              ,abun,w4da
     &              ,w3df,w3da,w3db,w3dc)
         call advx1 (d   ,w3dd
     &              ,w3dg,w3de
     &              ,w3dh,er
     .              ,w4da,abun
     &              ,w3df,w3da,w3db,w3dc)
         call advx2 (w3dd,d
     &              ,w3de,w3dg
     &              ,er  ,w3dh
     .              ,abun,w4da
     &              ,w3df,w3da,w3db,w3dc)
c
         ix1x2x3 = 6
         goto 10
c
c---------------- directional split in X1-X3-X2 fashion ----------------
c
       else ! if (ix1x2x3 .eq. 6) then
c
c       subroutine advx1 (dlo,den
c     &                  ,eod,edn
c     &                  ,mflx,s1,s2,s3)
c
         call advx1 (w3dd,d   
     &              ,w3de,w3dg
     &              ,er  ,w3dh
     .              ,abun,w4da
     &              ,w3df,w3da,w3db,w3dc)
         call advx3 (d   ,w3dd
     &              ,w3dg,w3de
     &              ,w3dh,er
     .              ,w4da,abun
     &              ,w3df,w3da,w3db,w3dc)
         call advx2 (w3dd,d
     &              ,w3de,w3dg
     &              ,er  ,w3dh
     .              ,abun,w4da
     &              ,w3df,w3da,w3db,w3dc)
c
         ix1x2x3 = 1
         goto 10
       endif
c
c-----------------------------------------------------------------------
c
10     continue
c
c Mark momentum density (velocity) boundary values out of date.  
c The d and e boundary values were maked out of date in advx*.
c
       do 20 i=1,6
         bvstat(i,3) = 0  !  v1
         bvstat(i,4) = 0  !  v2
         bvstat(i,5) = 0  !  v3
20     continue
c
c The velocities need to be computed from momentum densities.
c This will be done in nudt/newdt.
c
      go to 777
666   continue
      if(lrad .ne. 0) then
c
c Skipping the transport step.  Emulate it by copying er to w3dh and 
c marking the er boundary data out of date.
c
       do k=ks,ke
         do j=js,je
           do i=is,ie
             w3dh(i,j,k) = er(i,j,k)  !  Really er/d; fixed in newdt
           enddo ! i
         enddo ! j
       enddo ! k
c
       do i=1,6
         bvstat(i,6) = 0  !  er
       enddo
c
      endif ! lrad
c
777   return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////              T R A N S P R T              \\\\\\\\\\
c
c=======================================================================
c
