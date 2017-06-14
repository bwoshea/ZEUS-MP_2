c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 M A P I N                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c     PURPOSE: maps restart file arrays into ZEUS-MP arrays
c
c     Written by: John Hayes
c
c=======================================================================
      subroutine mapin
c
      use real_prec
      use config
      use param
      use grid
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use field
      use bndry
      use root
      use gravmod
      use restart_arrays
c
      implicit NONE
c
      character*15 :: filename
c
      integer  :: i, j, k, index, l
c
c=======================================================================
c     grid.h arrays
c=======================================================================
c
c-----------------------------------------------------------------------
c     integers
c-----------------------------------------------------------------------
c
      is    = ntgrdvr(1)
      ie    = ntgrdvr(2)
      js    = ntgrdvr(3)
      je    = ntgrdvr(4)
      ks    = ntgrdvr(5)
      ke    = ntgrdvr(6)
      iga   = ntgrdvr(7)
      jga   = ntgrdvr(8)
      kga   = ntgrdvr(9)
      igcon = ntgrdvr(10)
      nx1z  = ntgrdvr(11)
      nx2z  = ntgrdvr(12)
      nx3z  = ntgrdvr(13)
c
c-----------------------------------------------------------------------
c     Fixed coordinate arrays
c-----------------------------------------------------------------------
c
      index = 0
c
      call get_grid_array(x1a,in,index)
      call get_grid_array(x2a,jn,index)
      call get_grid_array(x3a,kn,index)
c
      call get_grid_array(x1ai,in,index)
      call get_grid_array(x2ai,jn,index)
      call get_grid_array(x3ai,kn,index)
c
      call get_grid_array(dx1a,in,index)
      call get_grid_array(dx2a,jn,index)
      call get_grid_array(dx3a,kn,index)
c
      call get_grid_array(dx1ai,in,index)
      call get_grid_array(dx2ai,jn,index)
      call get_grid_array(dx3ai,kn,index)
c
      call get_grid_array(vol1a,in,index)
      call get_grid_array(vol2a,jn,index)
      call get_grid_array(vol3a,kn,index)
c
      call get_grid_array(dvl1a,in,index)
      call get_grid_array(dvl2a,jn,index)
      call get_grid_array(dvl3a,kn,index)
c
      call get_grid_array(dvl1ai,in,index)
      call get_grid_array(dvl2ai,jn,index)
      call get_grid_array(dvl3ai,kn,index)
c
      call get_grid_array(g2a,in,index)
      call get_grid_array(g31a,in,index)
      call get_grid_array(dg2ad1,in,index)
      call get_grid_array(g2ai,in,index)
      call get_grid_array(g31ai,in,index)
      call get_grid_array(dg31ad1,in,index)
c
      call get_grid_array(g32a,jn,index)
      call get_grid_array(g32ai,jn,index)
      call get_grid_array(dg32ad2,jn,index)
      call get_grid_array(g4a,jn,index)
c
      call get_grid_array(x1b,in,index)
      call get_grid_array(x2b,jn,index)
      call get_grid_array(x3b,kn,index)
c
c
      call get_grid_array(x1bi,in,index)
      call get_grid_array(x2bi,jn,index)
      call get_grid_array(x3bi,kn,index)
c
      call get_grid_array(dx1b,in,index)
      call get_grid_array(dx2b,jn,index)
      call get_grid_array(dx3b,kn,index)
c
      call get_grid_array(dx1bi,in,index)
      call get_grid_array(dx2bi,jn,index)
      call get_grid_array(dx3bi,kn,index)
c
      call get_grid_array(vol1b,in,index)
      call get_grid_array(vol2b,jn,index)
      call get_grid_array(vol3b,kn,index)
c
      call get_grid_array(dvl1b,in,index)
      call get_grid_array(dvl2b,jn,index)
      call get_grid_array(dvl3b,kn,index)
c
      call get_grid_array(dvl1bi,in,index)
      call get_grid_array(dvl2bi,jn,index)
      call get_grid_array(dvl3bi,kn,index)
c
c
      call get_grid_array(g2b,in,index)
      call get_grid_array(g31b,in,index)
      call get_grid_array(dg2bd1,in,index)
      call get_grid_array(g2bi,in,index)
      call get_grid_array(g31bi,in,index)
      call get_grid_array(dg31bd1,in,index)
c
      call get_grid_array(g32b,jn,index)
      call get_grid_array(g32bi,jn,index)
      call get_grid_array(dg32bd2,jn,index)
      call get_grid_array(g4b,jn,index)
c
c-----------------------------------------------------------------------
c     grid velocities
c-----------------------------------------------------------------------
c
      call get_grid_array(vg1,in,index)
      call get_grid_array(vg2,jn,index)
      call get_grid_array(vg3,kn,index)
c
      x1fac = rlgrdvr(index+1)
      x2fac = rlgrdvr(index+2)
      x3fac = rlgrdvr(index+3)
      index = index + 3
c
c-----------------------------------------------------------------------
c     moving grid arrays -- 1 coordinate
c-----------------------------------------------------------------------
c
      call get_grid_array(x1ah,in,index)
      call get_grid_array(x1an,in,index)
      call get_grid_array(dx1ah,in,index)
      call get_grid_array(dx1an,in,index)
      call get_grid_array(dvl1ah,in,index)
      call get_grid_array(dvl1an,in,index)
      call get_grid_array(g2ah,in,index)
      call get_grid_array(g2an,in,index)
      call get_grid_array(g31ah,in,index)
      call get_grid_array(g31an,in,index)
      call get_grid_array(x1ahi,in,index)
      call get_grid_array(x1ani,in,index)
      call get_grid_array(dx1ahi,in,index)
      call get_grid_array(dx1ani,in,index)
      call get_grid_array(dvl1ahi,in,index)
      call get_grid_array(dvl1ani,in,index)
      call get_grid_array(g2ahi,in,index)
      call get_grid_array(g2ani,in,index)
      call get_grid_array(g31ahi,in,index)
      call get_grid_array(g31ani,in,index)
c
      call get_grid_array(x1bh,in,index)
      call get_grid_array(x1bn,in,index)
      call get_grid_array(dx1bh,in,index)
      call get_grid_array(dx1bn,in,index)
      call get_grid_array(dvl1bh,in,index)
      call get_grid_array(dvl1bn,in,index)
      call get_grid_array(g2bh,in,index)
      call get_grid_array(g2bn,in,index)
      call get_grid_array(g31bh,in,index)
      call get_grid_array(g31bn,in,index)
      call get_grid_array(x1bhi,in,index)
      call get_grid_array(x1bni,in,index)
      call get_grid_array(dx1bhi,in,index)
      call get_grid_array(dx1bni,in,index)
      call get_grid_array(dvl1bhi,in,index)
      call get_grid_array(dvl1bni,in,index)
      call get_grid_array(g2bhi,in,index)
      call get_grid_array(g2bni,in,index)
      call get_grid_array(g31bhi,in,index)
      call get_grid_array(g31bni,in,index)
c
c-----------------------------------------------------------------------
c     moving grid arrays -- 2 coordinate
c-----------------------------------------------------------------------
c
      call get_grid_array(x2ah,jn,index)
      call get_grid_array(x2an,jn,index)
      call get_grid_array(dx2ah,jn,index)
      call get_grid_array(dx2an,jn,index)
      call get_grid_array(dvl2ah,jn,index)
      call get_grid_array(dvl2an,jn,index)
      call get_grid_array(g32ah,jn,index)
      call get_grid_array(g32an,jn,index)
      call get_grid_array(g4ah,jn,index)
      call get_grid_array(g4an,jn,index)
      call get_grid_array(x2ahi,jn,index)
      call get_grid_array(x2ani,jn,index)
      call get_grid_array(dx2ahi,jn,index)
      call get_grid_array(dx2ani,jn,index)
      call get_grid_array(dvl2ahi,jn,index)
      call get_grid_array(dvl2ani,jn,index)
      call get_grid_array(g32ahi,jn,index)
      call get_grid_array(g32ani,jn,index)
c
      call get_grid_array(x2bh,jn,index)
      call get_grid_array(x2bn,jn,index)
      call get_grid_array(dx2bh,jn,index)
      call get_grid_array(dx2bn,jn,index)
      call get_grid_array(dvl2bh,jn,index)
      call get_grid_array(dvl2bn,jn,index)
      call get_grid_array(g32bh,jn,index)
      call get_grid_array(g32bn,jn,index)
      call get_grid_array(g4bh,jn,index)
      call get_grid_array(g4bn,jn,index)
      call get_grid_array(x2bhi,jn,index)
      call get_grid_array(x2bni,jn,index)
      call get_grid_array(dx2bhi,jn,index)
      call get_grid_array(dx2bni,jn,index)
      call get_grid_array(dvl2bhi,jn,index)
      call get_grid_array(dvl2bni,jn,index)
      call get_grid_array(g32bhi,jn,index)
      call get_grid_array(g32bni,jn,index)
c
c-----------------------------------------------------------------------
c     moving grid arrays -- 3 coordinate
c-----------------------------------------------------------------------
c
      call get_grid_array(x3ah,kn,index)
      call get_grid_array(x3an,kn,index)
      call get_grid_array(dx3ah,kn,index)
      call get_grid_array(dx3an,kn,index)
      call get_grid_array(dvl3ah,kn,index)
      call get_grid_array(dvl3an,kn,index)
      call get_grid_array(x3ahi,kn,index)
      call get_grid_array(x3ani,kn,index)
      call get_grid_array(dx3ahi,kn,index)
      call get_grid_array(dx3ani,kn,index)
      call get_grid_array(dvl3ahi,kn,index)
      call get_grid_array(dvl3ani,kn,index)
c
      call get_grid_array(x3bh,kn,index)
      call get_grid_array(x3bn,kn,index)
      call get_grid_array(dx3bh,kn,index)
      call get_grid_array(dx3bn,kn,index)
      call get_grid_array(dvl3bh,kn,index)
      call get_grid_array(dvl3bn,kn,index)
      call get_grid_array(x3bhi,kn,index)
      call get_grid_array(x3bni,kn,index)
      call get_grid_array(dx3bhi,kn,index)
      call get_grid_array(dx3bni,kn,index)
      call get_grid_array(dvl3bhi,kn,index)
      call get_grid_array(dvl3bni,kn,index)
c
      if(.false.) then
       write(*,"('GRID arrays on PE ',i2,' read ',i8,' expected ',i8)")
     .            myid,index,ngridr+mgridr
      endif
c
c=======================================================================
c     field.h arrays
c=======================================================================
c
      index = 1
      call get_field_reals ( d,index)
      call get_field_reals ( e,index)
      call get_field_reals (v1,index)
      call get_field_reals (v2,index)
      call get_field_reals (v3,index)
c
      if(xmhd) then
       call get_  field_reals (b1,index)
       call get_  field_reals (b2,index)
       call get_  field_reals (b3,index)
      endif ! xmhd
c
      if(lrad .ne. 0) then
       do k = 1, kn
        do j = 1, jn
         do i = 1, in
          er(i,j,k) = rlfldvr(index)
          index = index + 1
         enddo
        enddo
       enddo
      endif ! lrad
c
      if(XGRAV) then
       call get_  field_reals (gp,index)
      endif
c
      if(nspec .gt. 1) then
       do l = 1, nspec
        do k = 1, kn
         do j = 1, jn
          do i = 1, in
           abun(i,j,k,l) = rlfldvr(index)
           index = index + 1
          enddo
         enddo
        enddo
       enddo
      endif ! NSPEC
c
      if(.false.) then
       write(*,"('FIELD REALS on PE ',i2,' read ',i8,' expected ',i8)")
     .            myid,index-1,nfieldr
      endif
c
c=======================================================================
c     bndry.h arrays
c=======================================================================
c
c-----------------------------------------------------------------------
c     real arrays
c-----------------------------------------------------------------------
c
      index = 1
      do i = 1, nbvar
       fiis(i) = rlbdryvr(index)
       index           = index + 1
      enddo
      do i = 1, nbvar
       fois(i) = rlbdryvr(index)
       index           = index + 1
      enddo
      do j = 1, nbvar
       fijs(j) = rlbdryvr(index)
       index           = index + 1
      enddo
      do j = 1, nbvar
       fojs(j) = rlbdryvr(index)
       index           = index + 1
      enddo
      do k = 1, nbvar
       fiks(k) = rlbdryvr(index)
       index           = index + 1
      enddo
      do k = 1, nbvar
       foks(k) = rlbdryvr(index)
       index           = index + 1
      enddo
c
      call get_L_MN_layers(diib,3,jn,kn,index)
      call get_L_MN_layers(doib,3,jn,kn,index)
      call get_L_MN_layers(dijb,3,in,kn,index)
      call get_L_MN_layers(dojb,3,in,kn,index)
      call get_L_MN_layers(dikb,3,in,jn,index)
      call get_L_MN_layers(dokb,3,in,jn,index)
c
      call get_L_MN_layers(eiib,2,jn,kn,index)
      call get_L_MN_layers(eoib,2,jn,kn,index)
      call get_L_MN_layers(eijb,2,in,kn,index)
      call get_L_MN_layers(eojb,2,in,kn,index)
      call get_L_MN_layers(eikb,2,in,jn,index)
      call get_L_MN_layers(eokb,2,in,jn,index)
c
      call get_L_MN_layers(v1iib,2,jn,kn,index)
      call get_L_MN_layers(v1oib,2,jn,kn,index)
      call get_L_MN_layers(v1ijb,2,in,kn,index)
      call get_L_MN_layers(v1ojb,2,in,kn,index)
      call get_L_MN_layers(v1ikb,2,in,jn,index)
      call get_L_MN_layers(v1okb,2,in,jn,index)
c
      call get_L_MN_layers(v2iib,2,jn,kn,index)
      call get_L_MN_layers(v2oib,2,jn,kn,index)
      call get_L_MN_layers(v2ijb,2,in,kn,index)
      call get_L_MN_layers(v2ojb,2,in,kn,index)
      call get_L_MN_layers(v2ikb,2,in,jn,index)
      call get_L_MN_layers(v2okb,2,in,jn,index)
c
      call get_L_MN_layers(v3iib,2,jn,kn,index)
      call get_L_MN_layers(v3oib,2,jn,kn,index)
      call get_L_MN_layers(v3ijb,2,in,kn,index)
      call get_L_MN_layers(v3ojb,2,in,kn,index)
      call get_L_MN_layers(v3ikb,2,in,jn,index)
      call get_L_MN_layers(v3okb,2,in,jn,index)
c
      if(xmhd) then
       call get_L_MN_layers(b1iib,2,jn,kn,index)
       call get_L_MN_layers(b1oib,2,jn,kn,index)
       call get_L_MN_layers(b1ijb,2,in,kn,index)
       call get_L_MN_layers(b1ojb,2,in,kn,index)
       call get_L_MN_layers(b1ikb,2,in,jn,index)
       call get_L_MN_layers(b1okb,2,in,jn,index)
c
       call get_L_MN_layers(b2iib,2,jn,kn,index)
       call get_L_MN_layers(b2oib,2,jn,kn,index)
       call get_L_MN_layers(b2ijb,2,in,kn,index)
       call get_L_MN_layers(b2ojb,2,in,kn,index)
       call get_L_MN_layers(b2ikb,2,in,jn,index)
       call get_L_MN_layers(b2okb,2,in,jn,index)
c
       call get_L_MN_layers(b3iib,2,jn,kn,index)
       call get_L_MN_layers(b3oib,2,jn,kn,index)
       call get_L_MN_layers(b3ijb,2,in,kn,index)
       call get_L_MN_layers(b3ojb,2,in,kn,index)
       call get_L_MN_layers(b3ikb,2,in,jn,index)
       call get_L_MN_layers(b3okb,2,in,jn,index)
c
       call get_L_MN_layers(emf1iib,3,jn,kn,index)
       call get_L_MN_layers(emf1oib,3,jn,kn,index)
       call get_L_MN_layers(emf1ijb,3,in,kn,index)
       call get_L_MN_layers(emf1ojb,3,in,kn,index)
       call get_L_MN_layers(emf1ikb,3,in,jn,index)
       call get_L_MN_layers(emf1okb,3,in,jn,index)
c
       call get_L_MN_layers(emf2iib,3,jn,kn,index)
       call get_L_MN_layers(emf2oib,3,jn,kn,index)
       call get_L_MN_layers(emf2ijb,3,in,kn,index)
       call get_L_MN_layers(emf2ojb,3,in,kn,index)
       call get_L_MN_layers(emf2ikb,3,in,jn,index)
       call get_L_MN_layers(emf2okb,3,in,jn,index)
c
       call get_L_MN_layers(emf3iib,3,jn,kn,index)
       call get_L_MN_layers(emf3oib,3,jn,kn,index)
       call get_L_MN_layers(emf3ijb,3,in,kn,index)
       call get_L_MN_layers(emf3ojb,3,in,kn,index)
       call get_L_MN_layers(emf3ikb,3,in,jn,index)
       call get_L_MN_layers(emf3okb,3,in,jn,index)
      endif ! xmhd
c
      if(lrad .ne. 0) then
       call get_L_MN_layers(eriib,2,jn,kn,index)
       call get_L_MN_layers(eroib,2,jn,kn,index)
       call get_L_MN_layers(erijb,2,in,kn,index)
       call get_L_MN_layers(erojb,2,in,kn,index)
       call get_L_MN_layers(erikb,2,in,jn,index)
       call get_L_MN_layers(erokb,2,in,jn,index)
      endif ! lrad
c
      if(XGRAV) then
       call get_L_MN_layers(gpiib,2,jn,kn,index)
       call get_L_MN_layers(gpoib,2,jn,kn,index)
       call get_L_MN_layers(gpijb,2,in,kn,index)
       call get_L_MN_layers(gpojb,2,in,kn,index)
       call get_L_MN_layers(gpikb,2,in,jn,index)
       call get_L_MN_layers(gpokb,2,in,jn,index)
      endif ! XGRAV
c
      if(nspec .gt. 1) then
       call get_nspec_MN_layers(abiib,nspec,jn,kn,index)
       call get_nspec_MN_layers(aboib,nspec,jn,kn,index)
       call get_nspec_MN_layers(abijb,nspec,in,kn,index)
       call get_nspec_MN_layers(abojb,nspec,in,kn,index)
       call get_nspec_MN_layers(abikb,nspec,in,jn,index)
       call get_nspec_MN_layers(abokb,nspec,in,jn,index)
c
      endif ! NSPEC
c
      if(.false.) then
       write(*,"('BNDRY REALS on PE ',i2,' read ',i8,' expected ',i8)")
     .            myid,index-1,nbdryr
      endif
c
c-----------------------------------------------------------------------
c     integer arrays
c-----------------------------------------------------------------------
c
      index = 1
      do i = 1, 3
       niis(i) = ntbdryvr(index)
       index           = index + 1
      enddo
      do i = 1, 3
       nois(i) = ntbdryvr(index)
       index           = index + 1
      enddo
      do i = 1, 3
       nijs(i) = ntbdryvr(index)
       index           = index + 1
      enddo
      do i = 1, 3
       nojs(i) = ntbdryvr(index)
       index           = index + 1
      enddo
      do i = 1, 3
       niks(i) = ntbdryvr(index)
       index           = index + 1
      enddo
      do i = 1, 3
       noks(i) = ntbdryvr(index)
       index           = index + 1
      enddo
      do i = 1, 6
       do j = 1, nbvar
        bvstat(i,j) = ntbdryvr(index)
        index          = index + 1
       enddo
      enddo
c
      call get_MN_integers(niib  ,jn,kn,index)
      call get_MN_integers(niib2 ,jn,kn,index)
      call get_MN_integers(niib3 ,jn,kn,index)
      call get_MN_integers(niib23,jn,kn,index)
      call get_MN_integers(noib  ,jn,kn,index)
      call get_MN_integers(noib2 ,jn,kn,index)
      call get_MN_integers(noib3 ,jn,kn,index)
      call get_MN_integers(noib23,jn,kn,index)
c
      call get_MN_integers(nijb  ,in,kn,index)
      call get_MN_integers(nijb3 ,in,kn,index)
      call get_MN_integers(nijb1 ,in,kn,index)
      call get_MN_integers(nijb31,in,kn,index)
      call get_MN_integers(nojb  ,in,kn,index)
      call get_MN_integers(nojb3 ,in,kn,index)
      call get_MN_integers(nojb1 ,in,kn,index)
      call get_MN_integers(nojb31,in,kn,index)
c
      call get_MN_integers(nikb  ,in,jn,index)
      call get_MN_integers(nikb1 ,in,jn,index)
      call get_MN_integers(nikb2 ,in,jn,index)
      call get_MN_integers(nikb12,in,jn,index)
      call get_MN_integers(nokb  ,in,jn,index)
      call get_MN_integers(nokb1 ,in,jn,index)
      call get_MN_integers(nokb2 ,in,jn,index)
      call get_MN_integers(nokb12,in,jn,index)
c
      if(lrad .ne. 0) then
       call get_MN_integers(liib ,jn,kn,index)
       call get_MN_integers(loib ,jn,kn,index)
       call get_MN_integers(lijb ,in,kn,index)
       call get_MN_integers(lojb ,in,kn,index)
       call get_MN_integers(likb ,in,jn,index)
       call get_MN_integers(lokb ,in,jn,index)
      endif ! LRAD
c
      if(.false.) then
       write(*,"('BNDRY INTS on PE ',i2,' read ',i8,' expected ',i8)")
     .            myid,index-1,nbdryi
      endif
c
c=======================================================================
c     root.h arrays
c=======================================================================
c
c-----------------------------------------------------------------------
c     reals
c-----------------------------------------------------------------------
c
      b1floor = rlrtvr( 1)
      b2floor = rlrtvr( 2)
      b3floor = rlrtvr( 3)
      ciso    = rlrtvr( 4)
      courno  = rlrtvr( 5)
      dfloor  = rlrtvr( 6)
      dtal    = rlrtvr( 7)
      dtcs    = rlrtvr( 8)
      dtv1    = rlrtvr( 9)
      dtv2    = rlrtvr(10)
      dtv3    = rlrtvr(11)
      dtqq    = rlrtvr(12)
      dtnew   = rlrtvr(13)
      avisc_dt = rlrtvr(14)
      dtrd    = rlrtvr(15)
      dt      = rlrtvr(16)
      dtdump  = rlrtvr(17)
      dthdf   = rlrtvr(18)
      dthist  = rlrtvr(19)
      dtmin   = rlrtvr(20)
      dttsl   = rlrtvr(21)
      dtqqi2  = rlrtvr(22)
      dtusr   = rlrtvr(23)
      dtnri2  = rlrtvr(24)
      dtrdi2  = rlrtvr(25)
      dtimrdi2= rlrtvr(26)
      efloor  = rlrtvr(27)
      erfloor = rlrtvr(28)
      gamma   = rlrtvr(29)
      gamm1   = rlrtvr(30)
      qcon    = rlrtvr(31)
      qlin    = rlrtvr(32)
      tdump   = rlrtvr(33)
      thdf    = rlrtvr(34)
      thist   = rlrtvr(35)
      time    = rlrtvr(36)
      tlim    = rlrtvr(37)
      cpulim  = rlrtvr(38)
      trem    = rlrtvr(39)
      tsave   = rlrtvr(40)
CJH   we don't store tused to ensure that restart files can be subjected
CJH   to a repeatability (diff) test.
      tused   = 0.0 ! rlrtvr(41)
      tusr    = rlrtvr(42)
      ttsl    = rlrtvr(43)
      v1floor = rlrtvr(44)
      v2floor = rlrtvr(45)
      v3floor = rlrtvr(46)
      emf1floor = rlrtvr(47)
      emf2floor = rlrtvr(48)
      emf3floor = rlrtvr(49)
      gpfloor = rlrtvr(50)
c
c-----------------------------------------------------------------------
c     integers
c-----------------------------------------------------------------------
c
      ifsen(1) = ntrtvr( 1)
      ifsen(2) = ntrtvr( 2)
      ifsen(3) = ntrtvr( 3)
      ifsen(4) = ntrtvr( 4)
      ifsen(5) = ntrtvr( 5)
      ifsen(6) = ntrtvr( 6)
      idebug   = ntrtvr( 7)
      iordb1   = ntrtvr( 8)
      iordb2   = ntrtvr( 9)
      iordb3   = ntrtvr(10)
      iordd    = ntrtvr(11)
      iorde    = ntrtvr(12)
      iorder   = ntrtvr(13)
      iords1   = ntrtvr(14)
      iords2   = ntrtvr(15)
      iords3   = ntrtvr(16)
      istpb1   = ntrtvr(17)
      istpb2   = ntrtvr(18)
      istpb3   = ntrtvr(19)
      istpd    = ntrtvr(20)
      istpe    = ntrtvr(21)
      istper   = ntrtvr(22)
      istps1   = ntrtvr(23)
      istps2   = ntrtvr(24)
      istps3   = ntrtvr(25)
      ix1x2x3  = ntrtvr(26)
      jx1x2x3  = ntrtvr(27)
      nhy      = ntrtvr(28)
      nlim     = ntrtvr(29)
      nred     = ntrtvr(30)
      mbatch   = ntrtvr(31)
      nwarn    = ntrtvr(32)
      nseq     = ntrtvr(33)
      flstat   = ntrtvr(34)
      ioinp    = ntrtvr(35)
      iotsl    = ntrtvr(36)
      iolog    = ntrtvr(37)
      iohst    = ntrtvr(38)
      iomov    = ntrtvr(39)
      iores    = ntrtvr(40)
      ioshl    = ntrtvr(41)
c
c-----------------------------------------------------------------------
c     characters
c-----------------------------------------------------------------------
c
      hdffile = chrtvr(1:15)
      hstfile = chrtvr(16:30)
      resfile = chrtvr(31:45)
      usrfile = chrtvr(46:60)
      tslfile = chrtvr(61:68)
c
c=======================================================================
c     gravity.h arrays
c=======================================================================
c
      tgrav      = rlgrvvr(1)
      ptmass     = rlgrvvr(2)
      x1ptm      = rlgrvvr(3)
      x2ptm      = rlgrvvr(4)
      x3ptm      = rlgrvvr(5)
c
      gsup       = ntgrvvr(1)
c
      return
      end
c
c=======================================================================
c=======================================================================
c
      subroutine get_MN_integers(target,jn,kn,index)
c
      use real_prec
      use restart_arrays
c
      implicit NONE
c
      integer :: jn, kn
      integer :: target(jn,kn), index, i, j, k
c
      do k = 1, kn
       do j = 1, jn
        target(j,k) = ntbdryvr(index)
        index       = index + 1
       enddo
      enddo
c
      return
      end
c
      subroutine put_MN_integers(source,jn,kn,index)
c
      use real_prec
      use restart_arrays
c
      implicit NONE
c
      integer :: jn, kn
      integer :: source(jn,kn), index, i, j, k
c
      do k = 1, kn
       do j = 1, jn
        ntbdryvr(index) = source(j,k)
        index           = index + 1
       enddo
      enddo
c
      return
      end
c
      subroutine get_nspec_MN_layers(target,nspec,jn,kn,index)
c
      use real_prec
      use restart_arrays
c
      implicit NONE
c
      integer  :: index, i, j, k, l, jn, kn, nspec
      real(rl) :: target(jn,kn,2,nspec)
c
      do l = 1, nspec
       do k = 1, kn
        do j = 1, jn
         target(j,k,1,l) = rlbdryvr(index)
         target(j,k,2,l) = rlbdryvr(index+1)
         index = index + 2
        enddo
       enddo
      enddo
c
      return
      end
c
      subroutine put_nspec_MN_layers(source,nspec,jn,kn,index)
c
      use real_prec
      use restart_arrays
c
      implicit NONE
c
      integer  :: index, i, j, k, l, jn, kn, nspec
      real(rl) :: source(jn,kn,2,nspec)
c
      do l = 1, nspec
       do k = 1, kn
        do j = 1, jn
         rlbdryvr(index)   = source(j,k,1,l)
         rlbdryvr(index+1) = source(j,k,2,l)
         index = index + 2
        enddo
       enddo
      enddo
c
      return
      end
c
      subroutine get_L_MN_layers(target,l,m,n,index)
c
      use real_prec
      use param
      use restart_arrays
c
      implicit NONE
c
      integer  :: index, i, j, k, l, m, n
      real(rl) :: target(m,n,l)
c
      do i = 1, l
       do k = 1, n
        do j = 1, m
         target(j,k,i) = rlbdryvr(index)
         index = index + 1
        enddo
       enddo
      enddo
c
      return
      end
c
      subroutine put_L_MN_layers(source,l,m,n,index)
c
      use real_prec
      use param
      use restart_arrays
c
      implicit NONE
c
      integer  :: index, i, j, k, l, m, n
      real(rl) :: source(m,n,l)
c
      do i = 1, l
       do k = 1, n
        do j = 1, m
         rlbdryvr(index) = source(j,k,i)
         index = index + 1
        enddo
       enddo
      enddo
c
      return
      end
c
      subroutine get_grid_array(target,in,index)
c
      use real_prec
      use restart_arrays
c
      implicit NONE
c
      integer  :: index, i, in
      real(rl) :: target(in)
c
      do i = 1, in
       target(i) = rlgrdvr(i+index)
      enddo
      index = index + in
c
      return
      end
c
      subroutine put_grid_array(source,in,index)
c
      use real_prec
      use restart_arrays
c
      implicit NONE
c
      integer  :: index, i, in
      real(rl) :: source(in)
c
      do i = 1, in
       rlgrdvr(i+index) = source(i)
      enddo
      index = index + in
c
      return
      end
c
      subroutine get_field_reals(target,index)
c
      use real_prec
      use param
      use restart_arrays
c
      implicit NONE
c
      integer  :: index, i, j, k
      real(rl) :: target(in,jn,kn)
c
      do k = 1, kn
       do j = 1, jn
        do i = 1, in
         target(i,j,k) = rlfldvr(index)
         index = index + 1
        enddo
       enddo
      enddo
c
      return
      end
c
      subroutine put_field_reals(source,index)
c
      use real_prec
      use param
      use restart_arrays
c
      implicit NONE
c
      integer  :: index, i, j, k
      real(rl) :: source(in,jn,kn)
c
      do k = 1, kn
       do j = 1, jn
        do i = 1, in
         rlfldvr(index) = source(i,j,k)
         index = index + 1
        enddo
       enddo
      enddo
c
      return
      end
c
      subroutine shift_bndry_values(object,jnnew)
c
      use real_prec
      use param
c
      implicit NONE
c
      real(rl) :: object(in,jn,kn)
      integer  :: jnnew, i, k
c
      do k = 1, kn
       do i = 1, in
        object(i,jnnew-2,k) = object(i,4,k)
        object(i,jnnew-1,k) = object(i,5,k)
        object(i,jnnew  ,k) = object(i,6,k)
       enddo
      enddo
c
      return
      end
c
      subroutine fill_twod_fields(object,jnnew)
c
      use real_prec
      use param
c
      implicit NONE
c
      real(rl) :: object(in,jn,kn)
      integer  :: i,j,k,jnnew
c
      do k = 1, kn
       do j = 4, jnnew-3
        do i = 1, in
         object(i,j,k) = object(i,3,k)
        enddo
       enddo
      enddo
c
      return
      end
c
