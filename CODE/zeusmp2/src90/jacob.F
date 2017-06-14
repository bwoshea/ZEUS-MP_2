c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                  J A C O B                \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c=======================================================================
c
      subroutine jacob(eold, erold, dd, rhs, fn1, fn2, 
     .                 d_fn1_d_er, d_fn1_d_eg,
     .                 d_fn2_d_er,
     .                 d_fn2_d_eg, grdvcf, deldotv,deldotf)
c
c     completes the jacobian and the RHS for the Newton-Raphson
c     solution to the implicit radiation energy equation
c
c     written by: John C. Hayes; March, 1997
c
c     modified for F90 by J. Hayes; May, 2003
c
      use real_prec
      use config
      use param
      use field
      use bndry
      use grid
      use root
      use radiation
      use opac
      use cons
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
c     loop indices and loop bounds
c
      integer  :: i, j, k
      integer  :: jsm2, jep1, jep2, ksm2, kep1, kep2
c
c     matrix arrays
c
      real(rl) :: dd(neqm,neqm,in,jn,kn)
c
c     the right-hand side
c
      real(rl) :: rhs(neqm,in,jn,kn)
c
c     scalar work variables
c
      real(rl) :: d_f1p_d_e, d_f1_d_e,
     .            d_f2p_d_e, d_f2_d_e,
     .            d_f3p_d_e, d_f3_d_e,
     .            cancel_term
c
c     work arrays
c
      real(rl) :: f1           (in,jn,kn)
      real(rl) :: f2           (in,jn,kn)
      real(rl) :: f3           (in,jn,kn)
      real(rl) :: deldotf      (in,jn,kn)
      real(rl) :: deldotv      (in,jn,kn)
      real(rl) :: d_deldotf_d_e(in,jn,kn)
      real(rl) :: d_grdvcp _d_e(in,jn,kn)
      real(rl) :: d_fn1_d_er   (in,jn,kn)
      real(rl) :: d_fn1_d_eg   (in,jn,kn)
      real(rl) :: d_fn2_d_er   (in,jn,kn)
      real(rl) :: d_fn2_d_eg   (in,jn,kn)
      real(rl) :: fn1          (in,jn,kn)
      real(rl) :: fn2          (in,jn,kn)
      real(rl) :: grdvcf       (in,jn,kn)
      real(rl) :: grdvcp       (in,jn,kn)
      real(rl) :: erold        (in,jn,kn)
      real(rl) :: e old        (in,jn,kn)
      real(rl) :: erb          (in,jn,kn)
      real(rl) :: kapb         (in,jn,kn)
      real(rl) :: kemb         (in,jn,kn)
      real(rl) :: bbb          (in,jn,kn)
      real(rl) :: pb           (in,jn,kn)
c
c     time-centering function
c
      real(rl) :: thbar, xn, xnp1
c
      thbar(xn,xnp1) = radth*xnp1 + (1.D0-radth)*xn
c
c======================================================================
c     Do all the points owned by this tile, and then alter any
c     values that lie on a boundary, if appropriate.
c======================================================================
c
c     individual terms for the main diagonal and the RHS
c
      if(ldimen .eq. 1) then
       jsm2 = js
       jep1 = je
       jep2 = je
       ksm2 = ks
       kep1 = ke
       kep2 = ke
      endif
      if(ldimen .eq. 2) then
       jsm2 = js-2
       jep1 = je+1
       jep2 = je+2
       ksm2 = ks
       kep1 = ke
       kep2 = ke
      endif
      if(ldimen .eq. 3) then
       jsm2 = js-2
       jep1 = je+1
       jep2 = je+2
       ksm2 = ks-2
       kep1 = ke+1
       kep2 = ke+2
      endif
c
      do k = ksm2, kep2
       do j = jsm2, jep2
        do i = is-2, ie+2
         erb(i,j,k) = thbar(erold(i,j,k),ern(i,j,k))
        enddo
       enddo
      enddo
c
      do k = ks, ke
       do j = js, je
        do i = is, ie
         kapb(i,j,k) = thbar(kap(i,j,k),kapn(i,j,k))
         kemb(i,j,k) = thbar(kem(i,j,k),kemn(i,j,k))
         bbb (i,j,k) = thbar(bb (i,j,k),bbn (i,j,k))
         pb  (i,j,k) = thbar(p  (i,j,k),pn  (i,j,k))
        enddo
       enddo
      enddo
c
      do k = ks, kep1
       do j = js, jep1
        do i = is, ie+1
         f1(i,j,k) = -dr1(i,j,k) * dx1bi(i) * 
     .                            (erb(i,j,k) - erb(i-1,j  ,k  ))
         if(ldimen .gt. 1) then
          f2(i,j,k) = -dr2(i,j,k) * dx2bi(j) * g2bi(i) *
     .                             (erb(i,j,k) - erb(i  ,j-1,k  ))
          if(ldimen .eq. 3) then
           f3(i,j,k) = -dr3(i,j,k) * dx3bi(k) * g31bi(i) * g32bi(j) *
     .                              (erb(i,j,k) - erb(i  ,j  ,k-1))
          else
           f3(i,j,k) = 0.0D0
          endif ! ldimen = 3
         else ! ldimen > 1
          f2(i,j,k) = 0.0D0
         endif ! ldimen > 1
        enddo
       enddo
      enddo
c
      do k = ks, ke
       do j = js, je
        do i = is, ie
         grdvcp (i,j,k) = grdvcf(i,j,k) * erb(i,j,k)
         deldotf(i,j,k) = dvl1ai(i) * 
     .                    ( g2a(i+1)*g31a(i+1)*f1(i+1,j  ,k  ) -
     .                      g2a(i  )*g31a(i  )*f1(i  ,j  ,k  )     )
        enddo
       enddo
      enddo
      if(ldimen .gt. 1) then
       do k = ks, ke
        do j = js, je
         do i = is, ie
          deldotf(i,j,k) = deldotf(i,j,k)
     .                   + g2bi(i) * dvl2ai(j) *
     .                     (          g32a(j+1)*f2(i  ,j+1,k  ) -
     .                                g32a(j  )*f2(i  ,j  ,k  )     )
         enddo
        enddo
       enddo
      endif
      if(ldimen .eq. 3) then
       do k = ks, ke
        do j = js, je
         do i = is, ie
          deldotf(i,j,k) = deldotf(i,j,k)
     .                  + g31bi(i) * g32bi(j) * dvl3ai(k) *
     .                    (                    f3(i  ,j  ,k+1) -
     .                                         f3(i  ,j  ,k  )     )
         enddo
        enddo
       enddo
      endif
c
      do k = ks, ke
       do j = js, je
        do i = is, ie
         d_deldotf_d_e(i,j,k) = 
     .                          dvl1ai(i) * 
     .     (  (g2a(i+1)*g31a(i+1))**2 * dr1(i+1,j  ,k  )*dvl1bi(i+1)
     .       +(g2a(i  )*g31a(i  ))**2 * dr1(i  ,j  ,k  )*dvl1bi(i  ) )
     .                        +
     .                          dvl2ai(j) * g2bi(i)**2 *
     .     (  (g32a(j+1)         )**2 * dr2(i  ,j+1,k  )*dvl2bi(j+1)
     .       +(g32a(j  )         )**2 * dr2(i  ,j  ,k  )*dvl2bi(j  ) )
     .                        +
     .                          dvl3ai(k) * g31bi(i)**2*g32bi(j)**2 *
     .     (                            dr3(i  ,j  ,k+1)*dvl3bi(k+1)
     .       +                          dr3(i  ,j  ,k  )*dvl3bi(k  ) )
         d_grdvcp_d_e(i,j,k)  = grdvcf(i,j,k)
        enddo
       enddo
      enddo
c
      do k = ks, ke
       do j = js, je
        do i = is, ie
c
         d_fn1_d_er(i,j,k) = 1.D0 + dt*radth*( clight*kemb(i,j,k)
     .                                        + d_deldotf_d_e(i,j,k) 
     .                                        + d_grdvcp_d_e(i,j,k))
c
         d_fn1_d_eg(i,j,k) =      - dt*radth*(
     .                        4.D0*pi*(bbb (i,j,k)*dkapde(i,j,k) +
     .                                 kapb(i,j,k)*dbbde (i,j,k) )
     .                            - clight* erb (i,j,k)*dkemde(i,j,k))
c
         d_fn2_d_er(i,j,k) =      -dt*radth*clight*kemb(i,j,k)
c
         d_fn2_d_eg(i,j,k) = 1.D0 + dt*radth*(
     .                       4.D0*pi*(bbb (i,j,k)*dkapde(i,j,k) + 
     .                                kapb(i,j,k)*dbbde (i,j,k) )
     .                           - clight* erb (i,j,k)*dkemde(i,j,k)
     .                           + deldotv(i,j,k)*dpde  (i,j,k)    )
        enddo
       enddo
      enddo
c
      do k = ks, ke
       do j = js, je
        do i = is, ie
c
         fn1(i,j,k) = ern(i,j,k) - erold(i,j,k) - dt * (
     .
     .                4.D0*pi*kapb(i,j,k)*bbb(i,j,k)
     .
     .              -      clight*kemb(i,j,k)*erb(i,j,k)    
     .
     .              -             deldotf(i,j,k) 
     .
     .              -             grdvcp (i,j,k)       )
c
         fn2(i,j,k) = en(i,j,k) - eold(i,j,k) - dt * (
     .
     .              - 4.D0*pi*kapb(i,j,k)*bbb(i,j,k)
     .
     .              +      clight*kemb(i,j,k)*erb(i,j,k)
     .
     .              -         pb(i,j,k)*deldotv(i,j,k))
        enddo
       enddo
      enddo
c
c     compute the main diagonal and the RHS
c
      do k = ks, ke
       do j = js, je
        do i = is, ie
         dd  (1,1,i,j,k) = (d_fn1_d_er(i,j,k) -
     .                      d_fn1_d_eg(i,j,k)*d_fn2_d_er(i,j,k) /
     .                                        d_fn2_d_eg(i,j,k))
     .                   *dvl1a(i)*dvl2a(j)*dvl3a(k)
         rhs (1,  i,j,k) = (d_fn1_d_eg(i,j,k)/d_fn2_d_eg(i,j,k)
     .                      * fn2     (i,j,k)  -  fn1(i,j,k))
     .                   *dvl1a(i)*dvl2a(j)*dvl3a(k)
        enddo
       enddo
      enddo
c
c======================================================================
c     Physical boundaries.  In the case of reflecting or in-flow
c     b.c.'s, the value of the main diagonal at physical boundaries
c     must be altered.  For the other b.c. types currently supported,
c     no further change is necessary as the derivatives are unaffected
c     in those cases.
c======================================================================
c
#ifdef MPI_USED
      if(coords(1) .eq. 0) then              ! inner 1-boundary
#endif /* MPI_USED */
       do k = ks, ke
        do j = js, je
         if(liib(j,k) .eq. 1 .or. liib(j,k) .eq. 2) then
          cancel_term    = dt * radth * dvl1ai(is) * g2a(is) *
     .                     g31a(is) * (-dr1(is,j,k)*dx1bi(is))
          dd(1,1,is,j,k) = dd(1,1,is,j,k) + cancel_term
     .                   * dvl1a(is)*dvl2a(j)*dvl3a(k)
         endif
        enddo
       enddo
#ifdef MPI_USED
      endif
      if(coords(1) .eq. ntiles(1)-1) then    ! outer 1-boundary
#endif /* MPI_USED */
       do k = ks, ke
        do j = js, je
         if(loib(j,k) .eq. 1 .or. loib(j,k) .eq. 2) then
          cancel_term    = -dt * radth * dvl1ai(ie) * g2a(ie+1) *
     .                      g31a(ie+1) * ( dr1(ie+1,j,k) * 
     .                      dx1bi(ie+1))
          dd(1,1,ie,j,k) = dd(1,1,ie,j,k) + cancel_term
     .                   * dvl1a(ie)*dvl2a(j)*dvl3a(k)
         endif
        enddo
       enddo
#ifdef MPI_USED
      endif
#endif /* MPI_USED */
      if(ldimen .gt. 1) then
#ifdef MPI_USED
      if(coords(2) .eq. 0) then              ! inner 2-boundary
#endif /* MPI_USED */
       do k = ks, ke
        do i = is, ie
         if(lijb(i,k) .eq. 1 .or. lijb(i,k) .eq. 2) then
          cancel_term      =  dt * radth * dvl2ai(js) * g2bi(i) *
     .                        g32a(js) * (-dr2(i,js,k) * 
     .                        dx2bi(js) * g2bi(i))
          dd(1,1,i,js,k) = dd(1,1,i,js,k) + cancel_term
     .                   * dvl1a(i)*dvl2a(js)*dvl3a(k)
         endif
        enddo
       enddo
#ifdef MPI_USED
      endif
      if(coords(2) .eq. ntiles(2)-1) then    ! outer 2-boundary
#endif /* MPI_USED */
       do k = ks, ke
        do i = is, ie
         if(lojb(i,k) .eq. 1 .or. lojb(i,k) .eq. 2) then
          cancel_term    = -dt * radth * dvl2ai(je) * g2bi(i) *
     .                      g32a(je+1) * ( dr2(i,je+1,k) * 
     .                      dx2bi(je+1) * g2bi(i))
          dd(1,1,i,je,k) = dd(1,1,i,je,k) + cancel_term
     .                   * dvl1a(i)*dvl2a(je)*dvl3a(k)
         endif
        enddo
       enddo
#ifdef MPI_USED
      endif
#endif /* MPI_USED */
      if(ldimen .eq. 3) then
#ifdef MPI_USED
      if(coords(3) .eq. 0) then              ! inner 3-boundary
#endif /* MPI_USED */
       do j = js, je
        do i = is, ie
         if(likb(i,j) .eq. 1 .or. likb(i,j) .eq. 2) then
          cancel_term     =  dt * radth * 
     .                       dvl3ai(ks) * g31bi(i) * g32bi(j) *
     .                       (-dr3(i,j,ks) * dx3bi(ks) *
     .                       g31bi(i) * g32bi(j))
         dd(1,1,i,j,ks) =  dd(1,1,i,j,ks) + cancel_term
     .                   * dvl1a(i)*dvl2a(j)*dvl3a(ks)
         endif
        enddo
       enddo
#ifdef MPI_USED
      endif
      if(coords(3) .eq. ntiles(3)-1) then    ! outer 3-boundary
#endif /* MPI_USED */
       do j = js, je
        do i = is, ie
         if(lokb(i,j) .eq. 1 .or. lokb(i,j) .eq. 2) then
          cancel_term   = -dt * radth *
     .                     dvl3ai(ke) * g31bi(i) * g32bi(j) *
     .                     ( dr3(i,j,ke+1) * dx3bi(ke+1) *
     .                     g31bi(i) * g32bi(j))
         dd(1,1,i,j,ke) =  dd(1,1,i,j,ke) + cancel_term
     .                   * dvl1a(i)*dvl2a(j)*dvl3a(ke)
         endif
        enddo
       enddo
#ifdef MPI_USED
      endif
#endif /* MPI_USED */
      endif ! ldimen = 3
      endif ! ldimen > 1
c
      return
      end
c=======================================================================c
c    \\\\\\\\\\        E N D    S U B R O U T I N E       //////////
c    //////////                  J A C O B                \\\\\\\\\\
c
c=======================================================================c
