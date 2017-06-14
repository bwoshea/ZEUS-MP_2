c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                M A P O U T                \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c     PURPOSE: maps ZEUS-MP arrays into restart file arrays
c
c     Written by: John Hayes
c
c     Modified:
c
c     1) January 5, 2007: added line to write nijb1 and nojb1 to restart
c        dump file (the corresponding read statements are already present
c        in mapin.F).
c
c=======================================================================
      subroutine mapout
c
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
      integer :: i, j, k, index, l
c
c=======================================================================
c     grid.h arrays
c=======================================================================
c
c-----------------------------------------------------------------------
c     integers
c-----------------------------------------------------------------------
c
      ntgrdvr(1) = is
      ntgrdvr(2) = ie
      ntgrdvr(3) = js
      ntgrdvr(4) = je
      ntgrdvr(5) = ks
      ntgrdvr(6) = ke
      ntgrdvr(7) = iga
      ntgrdvr(8) = jga
      ntgrdvr(9) = kga
      ntgrdvr(10) = igcon
      ntgrdvr(11) = nx1z
      ntgrdvr(12) = nx2z
      ntgrdvr(13) = nx3z
c
c-----------------------------------------------------------------------
c     Fixed coordinate arrays
c-----------------------------------------------------------------------
c
      index = 0
      call put_grid_array(x1a,in,index)
      call put_grid_array(x2a,jn,index)
      call put_grid_array(x3a,kn,index)
      call put_grid_array(x1ai,in,index)
      call put_grid_array(x2ai,jn,index)
      call put_grid_array(x3ai,kn,index)
      call put_grid_array(dx1a,in,index)
      call put_grid_array(dx2a,jn,index)
      call put_grid_array(dx3a,kn,index)
      call put_grid_array(dx1ai,in,index)
      call put_grid_array(dx2ai,jn,index)
      call put_grid_array(dx3ai,kn,index)
      call put_grid_array(vol1a,in,index)
      call put_grid_array(vol2a,jn,index)
      call put_grid_array(vol3a,kn,index)
      call put_grid_array(dvl1a,in,index)
      call put_grid_array(dvl2a,jn,index)
      call put_grid_array(dvl3a,kn,index)
      call put_grid_array(dvl1ai,in,index)
      call put_grid_array(dvl2ai,jn,index)
      call put_grid_array(dvl3ai,kn,index)
      call put_grid_array(g2a,in,index)
      call put_grid_array(g31a,in,index)
      call put_grid_array(dg2ad1,in,index)
      call put_grid_array(g2ai,in,index)
      call put_grid_array(g31ai,in,index)
      call put_grid_array(dg31ad1,in,index)
      call put_grid_array(g32a,jn,index)
      call put_grid_array(g32ai,jn,index)
      call put_grid_array(dg32ad2,jn,index)
      call put_grid_array(g4a,jn,index)
c
      call put_grid_array(x1b,in,index)
      call put_grid_array(x2b,jn,index)
      call put_grid_array(x3b,kn,index)
      call put_grid_array(x1bi,in,index)
      call put_grid_array(x2bi,jn,index)
      call put_grid_array(x3bi,kn,index)
      call put_grid_array(dx1b,in,index)
      call put_grid_array(dx2b,jn,index)
      call put_grid_array(dx3b,kn,index)
      call put_grid_array(dx1bi,in,index)
      call put_grid_array(dx2bi,jn,index)
      call put_grid_array(dx3bi,kn,index)
      call put_grid_array(vol1b,in,index)
      call put_grid_array(vol2b,jn,index)
      call put_grid_array(vol3b,kn,index)
      call put_grid_array(dvl1b,in,index)
      call put_grid_array(dvl2b,jn,index)
      call put_grid_array(dvl3b,kn,index)
      call put_grid_array(dvl1bi,in,index)
      call put_grid_array(dvl2bi,jn,index)
      call put_grid_array(dvl3bi,kn,index)
      call put_grid_array(g2b,in,index)
      call put_grid_array(g31b,in,index)
      call put_grid_array(dg2bd1,in,index)
      call put_grid_array(g2bi,in,index)
      call put_grid_array(g31bi,in,index)
      call put_grid_array(dg31bd1,in,index)
      call put_grid_array(g32b,jn,index)
      call put_grid_array(g32bi,jn,index)
      call put_grid_array(dg32bd2,jn,index)
      call put_grid_array(g4b,jn,index)
c-----------------------------------------------------------------------
c     grid velocities
c-----------------------------------------------------------------------
      call put_grid_array(vg1,in,index)
      call put_grid_array(vg2,jn,index)
      call put_grid_array(vg3,kn,index)
c
      rlgrdvr(index+1) = x1fac
      rlgrdvr(index+2) = x2fac
      rlgrdvr(index+3) = x3fac
      index = index + 3
c
c-----------------------------------------------------------------------
c     moving grid arrays -- 1 coordinate
c-----------------------------------------------------------------------
c
      call put_grid_array(x1ah,in,index)
      call put_grid_array(x1an,in,index)
      call put_grid_array(dx1ah,in,index)
      call put_grid_array(dx1an,in,index)
      call put_grid_array(dvl1ah,in,index)
      call put_grid_array(dvl1an,in,index)
      call put_grid_array(g2ah,in,index)
      call put_grid_array(g2an,in,index)
      call put_grid_array(g31ah,in,index)
      call put_grid_array(g31an,in,index)
      call put_grid_array(x1ahi,in,index)
      call put_grid_array(x1ani,in,index)
      call put_grid_array(dx1ahi,in,index)
      call put_grid_array(dx1ani,in,index)
      call put_grid_array(dvl1ahi,in,index)
      call put_grid_array(dvl1ani,in,index)
      call put_grid_array(g2ahi,in,index)
      call put_grid_array(g2ani,in,index)
      call put_grid_array(g31ahi,in,index)
      call put_grid_array(g31ani,in,index)
c
      call put_grid_array(x1bh,in,index)
      call put_grid_array(x1bn,in,index)
      call put_grid_array(dx1bh,in,index)
      call put_grid_array(dx1bn,in,index)
      call put_grid_array(dvl1bh,in,index)
      call put_grid_array(dvl1bn,in,index)
      call put_grid_array(g2bh,in,index)
      call put_grid_array(g2bn,in,index)
      call put_grid_array(g31bh,in,index)
      call put_grid_array(g31bn,in,index)
      call put_grid_array(x1bhi,in,index)
      call put_grid_array(x1bni,in,index)
      call put_grid_array(dx1bhi,in,index)
      call put_grid_array(dx1bni,in,index)
      call put_grid_array(dvl1bhi,in,index)
      call put_grid_array(dvl1bni,in,index)
      call put_grid_array(g2bhi,in,index)
      call put_grid_array(g2bni,in,index)
      call put_grid_array(g31bhi,in,index)
      call put_grid_array(g31bni,in,index)
c
c-----------------------------------------------------------------------
c     moving grid arrays -- 2 coordinate
c-----------------------------------------------------------------------
c
      call put_grid_array(x2ah,jn,index)
      call put_grid_array(x2an,jn,index)
      call put_grid_array(dx2ah,jn,index)
      call put_grid_array(dx2an,jn,index)
      call put_grid_array(dvl2ah,jn,index)
      call put_grid_array(dvl2an,jn,index)
      call put_grid_array(g32ah,jn,index)
      call put_grid_array(g32an,jn,index)
      call put_grid_array(g4ah,jn,index)
      call put_grid_array(g4an,jn,index)
      call put_grid_array(x2ahi,jn,index)
      call put_grid_array(x2ani,jn,index)
      call put_grid_array(dx2ahi,jn,index)
      call put_grid_array(dx2ani,jn,index)
      call put_grid_array(dvl2ahi,jn,index)
      call put_grid_array(dvl2ani,jn,index)
      call put_grid_array(g32ahi,jn,index)
      call put_grid_array(g32ani,jn,index)
c
      call put_grid_array(x2bh,jn,index)
      call put_grid_array(x2bn,jn,index)
      call put_grid_array(dx2bh,jn,index)
      call put_grid_array(dx2bn,jn,index)
      call put_grid_array(dvl2bh,jn,index)
      call put_grid_array(dvl2bn,jn,index)
      call put_grid_array(g32bh,jn,index)
      call put_grid_array(g32bn,jn,index)
      call put_grid_array(g4bh,jn,index)
      call put_grid_array(g4bn,jn,index)
      call put_grid_array(x2bhi,jn,index)
      call put_grid_array(x2bni,jn,index)
      call put_grid_array(dx2bhi,jn,index)
      call put_grid_array(dx2bni,jn,index)
      call put_grid_array(dvl2bhi,jn,index)
      call put_grid_array(dvl2bni,jn,index)
      call put_grid_array(g32bhi,jn,index)
      call put_grid_array(g32bni,jn,index)
c
c-----------------------------------------------------------------------
c     moving grid arrays -- 3 coordinate
c-----------------------------------------------------------------------
c
      call put_grid_array(x3ah,kn,index)
      call put_grid_array(x3an,kn,index)
      call put_grid_array(dx3ah,kn,index)
      call put_grid_array(dx3an,kn,index)
      call put_grid_array(dvl3ah,kn,index)
      call put_grid_array(dvl3an,kn,index)
      call put_grid_array(x3ahi,kn,index)
      call put_grid_array(x3ani,kn,index)
      call put_grid_array(dx3ahi,kn,index)
      call put_grid_array(dx3ani,kn,index)
      call put_grid_array(dvl3ahi,kn,index)
      call put_grid_array(dvl3ani,kn,index)
c
      call put_grid_array(x3bh,kn,index)
      call put_grid_array(x3bn,kn,index)
      call put_grid_array(dx3bh,kn,index)
      call put_grid_array(dx3bn,kn,index)
      call put_grid_array(dvl3bh,kn,index)
      call put_grid_array(dvl3bn,kn,index)
      call put_grid_array(x3bhi,kn,index)
      call put_grid_array(x3bni,kn,index)
      call put_grid_array(dx3bhi,kn,index)
      call put_grid_array(dx3bni,kn,index)
      call put_grid_array(dvl3bhi,kn,index)
      call put_grid_array(dvl3bni,kn,index)
c
      if(.false.) then
       write(*,"('GRID arrays on PE ',i2,' wrote ',i8,' expected ',i8)")
     .            myid,index,ngridr+mgridr
      endif
c
c=======================================================================
c     field.h arrays
c=======================================================================
c
      index = 1
      call put_field_reals( d, index)
      call put_field_reals( e, index)
      call put_field_reals(v1, index)
      call put_field_reals(v2, index)
      call put_field_reals(v3, index)
      if(xmhd) then
       call put_field_reals(b1, index)
       call put_field_reals(b2, index)
       call put_field_reals(b3, index)
      endif
      if(lrad .ne. 0) call put_field_reals(er, index)
      if(xgrav) call put_field_reals(gp, index)
      if(nspec .gt. 1) then
       do l = 1, nspec
        do k = 1, kn
         do j = 1, jn
          do i = 1, in
           rlfldvr(index) = abun(i,j,k,l)
           index = index + 1
          enddo
         enddo
        enddo
       enddo
      endif ! NSPEC
      if(.false.) then
       write(*,"('FIELD REALS on PE ',i2,' wrote ',i8,' expected ',i8)")
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
       rlbdryvr(index) = fiis(i)
       index           = index + 1
      enddo
      do i = 1, nbvar
       rlbdryvr(index) = fois(i)
       index           = index + 1
      enddo
      do j = 1, nbvar
       rlbdryvr(index) = fijs(j)
       index           = index + 1
      enddo
      do j = 1, nbvar
       rlbdryvr(index) = fojs(j)
       index           = index + 1
      enddo
      do k = 1, nbvar
       rlbdryvr(index) = fiks(k)
       index           = index + 1
      enddo
      do k = 1, nbvar
       rlbdryvr(index) = foks(k)
       index           = index + 1
      enddo
c
      call put_L_MN_layers(diib,3,jn,kn,index)
      call put_L_MN_layers(doib,3,jn,kn,index)
      call put_L_MN_layers(dijb,3,in,kn,index)
      call put_L_MN_layers(dojb,3,in,kn,index)
      call put_L_MN_layers(dikb,3,in,jn,index)
      call put_L_MN_layers(dokb,3,in,jn,index)
c
      call put_L_MN_layers(eiib,2,jn,kn,index)
      call put_L_MN_layers(eoib,2,jn,kn,index)
      call put_L_MN_layers(eijb,2,in,kn,index)
      call put_L_MN_layers(eojb,2,in,kn,index)
      call put_L_MN_layers(eikb,2,in,jn,index)
      call put_L_MN_layers(eokb,2,in,jn,index)
c
      call put_L_MN_layers(v1iib,2,jn,kn,index)
      call put_L_MN_layers(v1oib,2,jn,kn,index)
      call put_L_MN_layers(v1ijb,2,in,kn,index)
      call put_L_MN_layers(v1ojb,2,in,kn,index)
      call put_L_MN_layers(v1ikb,2,in,jn,index)
      call put_L_MN_layers(v1okb,2,in,jn,index)
c
      call put_L_MN_layers(v2iib,2,jn,kn,index)
      call put_L_MN_layers(v2oib,2,jn,kn,index)
      call put_L_MN_layers(v2ijb,2,in,kn,index)
      call put_L_MN_layers(v2ojb,2,in,kn,index)
      call put_L_MN_layers(v2ikb,2,in,jn,index)
      call put_L_MN_layers(v2okb,2,in,jn,index)
c
      call put_L_MN_layers(v3iib,2,jn,kn,index)
      call put_L_MN_layers(v3oib,2,jn,kn,index)
      call put_L_MN_layers(v3ijb,2,in,kn,index)
      call put_L_MN_layers(v3ojb,2,in,kn,index)
      call put_L_MN_layers(v3ikb,2,in,jn,index)
      call put_L_MN_layers(v3okb,2,in,jn,index)
c
      if(xmhd) then
       call put_L_MN_layers(b1iib,2,jn,kn,index)
       call put_L_MN_layers(b1oib,2,jn,kn,index)
       call put_L_MN_layers(b1ijb,2,in,kn,index)
       call put_L_MN_layers(b1ojb,2,in,kn,index)
       call put_L_MN_layers(b1ikb,2,in,jn,index)
       call put_L_MN_layers(b1okb,2,in,jn,index)
c
       call put_L_MN_layers(b2iib,2,jn,kn,index)
       call put_L_MN_layers(b2oib,2,jn,kn,index)
       call put_L_MN_layers(b2ijb,2,in,kn,index)
       call put_L_MN_layers(b2ojb,2,in,kn,index)
       call put_L_MN_layers(b2ikb,2,in,jn,index)
       call put_L_MN_layers(b2okb,2,in,jn,index)
c
       call put_L_MN_layers(b3iib,2,jn,kn,index)
       call put_L_MN_layers(b3oib,2,jn,kn,index)
       call put_L_MN_layers(b3ijb,2,in,kn,index)
       call put_L_MN_layers(b3ojb,2,in,kn,index)
       call put_L_MN_layers(b3ikb,2,in,jn,index)
       call put_L_MN_layers(b3okb,2,in,jn,index)
c
       call put_L_MN_layers(emf1iib,3,jn,kn,index)
       call put_L_MN_layers(emf1oib,3,jn,kn,index)
       call put_L_MN_layers(emf1ijb,3,in,kn,index)
       call put_L_MN_layers(emf1ojb,3,in,kn,index)
       call put_L_MN_layers(emf1ikb,3,in,jn,index)
       call put_L_MN_layers(emf1okb,3,in,jn,index)
c
       call put_L_MN_layers(emf2iib,3,jn,kn,index)
       call put_L_MN_layers(emf2oib,3,jn,kn,index)
       call put_L_MN_layers(emf2ijb,3,in,kn,index)
       call put_L_MN_layers(emf2ojb,3,in,kn,index)
       call put_L_MN_layers(emf2ikb,3,in,jn,index)
       call put_L_MN_layers(emf2okb,3,in,jn,index)
c
       call put_L_MN_layers(emf3iib,3,jn,kn,index)
       call put_L_MN_layers(emf3oib,3,jn,kn,index)
       call put_L_MN_layers(emf3ijb,3,in,kn,index)
       call put_L_MN_layers(emf3ojb,3,in,kn,index)
       call put_L_MN_layers(emf3ikb,3,in,jn,index)
       call put_L_MN_layers(emf3okb,3,in,jn,index)
      endif ! xmhd
      if(lrad .ne. 0) then
       call put_L_MN_layers(eriib,2,jn,kn,index)
       call put_L_MN_layers(eroib,2,jn,kn,index)
       call put_L_MN_layers(erijb,2,in,kn,index)
       call put_L_MN_layers(erojb,2,in,kn,index)
       call put_L_MN_layers(erikb,2,in,jn,index)
       call put_L_MN_layers(erokb,2,in,jn,index)
      endif
      if(xgrav) then
       call put_L_MN_layers(gpiib,2,jn,kn,index)
       call put_L_MN_layers(gpoib,2,jn,kn,index)
       call put_L_MN_layers(gpijb,2,in,kn,index)
       call put_L_MN_layers(gpojb,2,in,kn,index)
       call put_L_MN_layers(gpikb,2,in,jn,index)
       call put_L_MN_layers(gpokb,2,in,jn,index)
      endif
      if(nspec .gt. 1) then
       call put_nspec_MN_layers(abiib,nspec,jn,kn,index)
       call put_nspec_MN_layers(aboib,nspec,jn,kn,index)
       call put_nspec_MN_layers(abijb,nspec,in,kn,index)
       call put_nspec_MN_layers(abojb,nspec,in,kn,index)
       call put_nspec_MN_layers(abikb,nspec,in,jn,index)
       call put_nspec_MN_layers(abokb,nspec,in,jn,index)
      endif ! NSPEC
c
      if(.false.) then
       write(*,"('BNDRY REALS on PE ',i2,' wrote ',i8,' expected ',i8)")
     .            myid,index-1,nbdryr
      endif
c
c-----------------------------------------------------------------------
c     integer arrays
c-----------------------------------------------------------------------
c
      index = 1
      do i = 1, 3
       ntbdryvr(index) = niis(i)
       index           = index + 1
      enddo
      do i = 1, 3
       ntbdryvr(index) = nois(i)
       index           = index + 1
      enddo
      do i = 1, 3
       ntbdryvr(index) = nijs(i)
       index           = index + 1
      enddo
      do i = 1, 3
       ntbdryvr(index) = nojs(i)
       index           = index + 1
      enddo
      do i = 1, 3
       ntbdryvr(index) = niks(i)
       index           = index + 1
      enddo
      do i = 1, 3
       ntbdryvr(index) = noks(i)
       index           = index + 1
      enddo
      do i = 1, 6
       do j = 1, nbvar
        ntbdryvr(index) = bvstat(i,j)
        index           = index + 1
       enddo
      enddo
c
      call put_MN_integers(niib  ,jn,kn,index)
      call put_MN_integers(niib2 ,jn,kn,index)
      call put_MN_integers(niib3 ,jn,kn,index)
      call put_MN_integers(niib23,jn,kn,index)
      call put_MN_integers(noib  ,jn,kn,index)
      call put_MN_integers(noib2 ,jn,kn,index)
      call put_MN_integers(noib3 ,jn,kn,index)
      call put_MN_integers(noib23,jn,kn,index)
c
      call put_MN_integers(nijb  ,in,kn,index)
      call put_MN_integers(nijb3 ,in,kn,index)
      call put_MN_integers(nijb1 ,in,kn,index)
      call put_MN_integers(nijb31,in,kn,index)
      call put_MN_integers(nojb  ,in,kn,index)
      call put_MN_integers(nojb3 ,in,kn,index)
      call put_MN_integers(nojb1 ,in,kn,index)
      call put_MN_integers(nojb31,in,kn,index)
c
      call put_MN_integers(nikb  ,in,jn,index)
      call put_MN_integers(nikb1 ,in,jn,index)
      call put_MN_integers(nikb2 ,in,jn,index)
      call put_MN_integers(nikb12,in,jn,index)
      call put_MN_integers(nokb  ,in,jn,index)
      call put_MN_integers(nokb1 ,in,jn,index)
      call put_MN_integers(nokb2 ,in,jn,index)
      call put_MN_integers(nokb12,in,jn,index)
c
      if(lrad .ne. 0) then
       call put_MN_integers(liib ,jn,kn,index)
       call put_MN_integers(loib ,jn,kn,index)
       call put_MN_integers(lijb ,in,kn,index)
       call put_MN_integers(lojb ,in,kn,index)
       call put_MN_integers(likb ,in,jn,index)
       call put_MN_integers(lokb ,in,jn,index)
      endif ! lrad
      if(.false.) then
       write(*,"('BNDRY INTS on PE ',i2,' wrote ',i8,' expected ',i8)")
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
      rlrtvr( 1) = b1floor
      rlrtvr( 2) = b2floor
      rlrtvr( 3) = b3floor
      rlrtvr( 4) = ciso
      rlrtvr( 5) = courno
      rlrtvr( 6) = dfloor
      rlrtvr( 7) = dtal
      rlrtvr( 8) = dtcs
      rlrtvr( 9) = dtv1
      rlrtvr(10) = dtv2
      rlrtvr(11) = dtv3
      rlrtvr(12) = dtqq
      rlrtvr(13) = dtnew
      rlrtvr(14) = avisc_dt
      rlrtvr(15) = dtrd
      rlrtvr(16) = dt
      rlrtvr(17) = dtdump
      rlrtvr(18) = dthdf
      rlrtvr(19) = dthist
      rlrtvr(20) = dtmin
      rlrtvr(21) = dttsl
      rlrtvr(22) = dtqqi2
      rlrtvr(23) = dtusr
      rlrtvr(24) = dtnri2
      rlrtvr(25) = dtrdi2
      rlrtvr(26) = dtimrdi2
      rlrtvr(27) = efloor
      rlrtvr(28) = erfloor
      rlrtvr(29) = gamma
      rlrtvr(30) = gamm1
      rlrtvr(31) = qcon
      rlrtvr(32) = qlin
      rlrtvr(33) = tdump
      rlrtvr(34) = thdf
      rlrtvr(35) = thist
      rlrtvr(36) = time
      rlrtvr(37) = tlim
      rlrtvr(38) = cpulim
      rlrtvr(39) = trem
      rlrtvr(40) = tsave
CJH  -- storing tused causes restart files to fail "diff" test for
CJH     repeatability check
      rlrtvr(41) = 0.0 ! tused
      rlrtvr(42) = tusr
      rlrtvr(43) = ttsl
      rlrtvr(44) = v1floor
      rlrtvr(45) = v2floor
      rlrtvr(46) = v3floor
      rlrtvr(47) = emf1floor
      rlrtvr(48) = emf2floor
      rlrtvr(49) = emf3floor
      rlrtvr(50) = gpfloor
c
c-----------------------------------------------------------------------
c     integers
c-----------------------------------------------------------------------
c
      ntrtvr( 1) = ifsen(1)
      ntrtvr( 2) = ifsen(2)
      ntrtvr( 3) = ifsen(3)
      ntrtvr( 4) = ifsen(4)
      ntrtvr( 5) = ifsen(5)
      ntrtvr( 6) = ifsen(6)
      ntrtvr( 7) = idebug
      ntrtvr( 8) = iordb1
      ntrtvr( 9) = iordb2
      ntrtvr(10) = iordb3
      ntrtvr(11) = iordd
      ntrtvr(12) = iorde
      ntrtvr(13) = iorder
      ntrtvr(14) = iords1
      ntrtvr(15) = iords2
      ntrtvr(16) = iords3
      ntrtvr(17) = istpb1
      ntrtvr(18) = istpb2
      ntrtvr(19) = istpb3
      ntrtvr(20) = istpd
      ntrtvr(21) = istpe
      ntrtvr(22) = istper
      ntrtvr(23) = istps1
      ntrtvr(24) = istps2
      ntrtvr(25) = istps3
      ntrtvr(26) = ix1x2x3
      ntrtvr(27) = jx1x2x3
      ntrtvr(28) = nhy
      ntrtvr(29) = nlim
      ntrtvr(30) = nred
      ntrtvr(31) = mbatch
      ntrtvr(32) = nwarn
      ntrtvr(33) = nseq
      ntrtvr(34) = flstat
      ntrtvr(35) = ioinp
      ntrtvr(36) = iotsl
      ntrtvr(37) = iolog
      ntrtvr(38) = iohst
      ntrtvr(39) = iomov
      ntrtvr(40) = iores
      ntrtvr(41) = ioshl
c
c-----------------------------------------------------------------------
c     characters
c-----------------------------------------------------------------------
c
      chrtvr = hdffile // hstfile // resfile // usrfile
     .      // tslfile
c
c=======================================================================
c     gravity.h arrays
c=======================================================================
c
      rlgrvvr(1) = tgrav
      rlgrvvr(2) = ptmass
      rlgrvvr(3) = x1ptm
      rlgrvvr(4) = x2ptm
      rlgrvvr(5) = x3ptm
c
      ntgrvvr(1) = gsup
c
      return
      end
