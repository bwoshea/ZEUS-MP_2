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
       subroutine transprt_1D
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
      integer i,j,k,l
c
c-----------------------------------------------------------------------
      if(xhydro .eqv. .false.) go to 666
c
      if(xmhd) then
c
c      Transport the three components of B using Constrained Transport.
c
       call ct_1D
      endif ! xmhd
c
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
             w3db(i,j,k) = v2(i,j,k) * 0.5 * (d(i  ,j  ,k  ) + d(i,j,k))
     1                     * g2b(i)
             w3dc(i,j,k) = v3(i,j,k) * 0.5 * (d(i  ,j  ,k  ) + d(i,j,k))
     1                     * g31b(i) * g32b(j)
           ENDDO
         ENDDO
       ENDDO
C
c---------------- directional split in X1-X2-X3 fashion ----------------
c
       nseq = 0  ! in /root/
c
c       subroutine advx1 (dlo,den
c     &                  ,eod,edn
c     &                  ,ero,ern
c     &                  ,mflx,s1,s2,s3)
c        
         call advx1 (w3dd,d   
     &              ,w3de,w3dg
     &              ,er  ,w3dh
     .              ,abun,w4da
     &              ,w3df,w3da,w3db,w3dc)
c
10     continue
c
c Mark momentum density (velocity) boundary values out of date.  
c The d and e boundary values were maked out of date in advx*.
c
       do 20 i=1,6
         bvstat(i,3) = 0  !  v1
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
      endif ! lrad
c
c      Also need an array copy for ABUN and W4DA...
c
      if(nspec .gt. 1) then
       do l =  1, nspec
        do k = ks, ke
         do j = js, je
          do i = is, ie
           w4da(i,j,k,l) = abun(i,j,k,l)
          enddo
         enddo
        enddo
       enddo
      endif
c
      do i=1,6
       bvstat(i,6) = 0  !  er
       bvstat(i,7) = 0  !  abun
      enddo
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
