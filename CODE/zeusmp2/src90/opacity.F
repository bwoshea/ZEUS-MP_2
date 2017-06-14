/* #define ICF */
c=======================================================================
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                  O P A C                  \\\\\\\\\\
c
c=======================================================================
c
      subroutine opacity(e, d, gam, kr, kp, sg, dkpdt, km, dkedt,
     .                   ibeg, iend, jbeg, jend, kbeg, kend)
c
      use real_prec
      use config
      use param
      use cons
      use grid
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
      integer  :: i, ibeg, iend, j, jbeg, jend, k, kbeg, kend
c
      real(rl) :: e    (in,jn,kn), d   (in,jn,kn), t    (in,jn,kn), 
     .            p    (in,jn,kn), kp   (in,jn,kn),
     .            dtde (in,jn,kn), dpde(in,jn,kn), dbbdt(in,jn,kn),
     .            dkpdt(in,jn,kn), sg  (in,jn,kn), er   (in,jn,kn),
     .            kr   (in,jn,kn), km  (in,jn,kn), dkedt(in,jn,kn)
c
      real(rl) :: gam, coef, temp, binv, gam1, stef
      real(rl) :: rmfp0i
c
      real(rl) :: kap_floor, kap_ceil, sp_xscn, hundred
c
      real(rl) :: so_eps , tau    , eriibn, timarg
c
      common /soui/ so_eps, tau, eriibn
c
      hundred = 100.0d0
      stef = clight * rad_con / (4.D0 * pi)
      gam1 = gam - 1.0D0
c
      binv = 1.D0 / boltz
      do 20 k = kbeg, kend
      do 20 j = jbeg, jend
      do 20 i = ibeg, iend
#ifndef MARSHAK
       t(i,j,k) = (gam1)*mmw*mh*binv*e(i,j,k)/d(i,j,k)
#else
       t(i,j,k) = sqrt( sqrt( so_eps*e(i,j,k)/rad_con) )
#endif
20    continue
c
      rmfp0i  = 1.D0 / rmfp0
      do 50 k = kbeg, kend
      do 50 j = jbeg, jend
      do 50 i = ibeg, iend
#ifndef MARSHAK
         kr  (i,j,k) = rmfp0i * (d(i,j,k)/rho0)**xnu *
     .                         (t(i,j,k)/ t_0)**(-powr)
#ifdef ICF
         kr  (i,j,k) = max(min_coef, kr(i,j,k))
         kr  (i,j,k) = min(max_coef, kr(i,j,k))
#endif
         sg   (i,j,k) = 0.D0
#ifndef ICF
         kp   (i,j,k) = kr(i,j,k)
#else
         kp   (i,j,k) = min(hundred,kpfrac*kr(i,j,k))
#endif /* ICF */
         dkpdt(i,j,k) = -powr * kp(i,j,k) / t(i,j,k)
#else
         kr   (i,j,k) = rmfp0i
         kp   (i,j,k) = kr(i,j,k)
         dkpdt(i,j,k) = 0.0
#endif /* MARSHAK */
50    continue
      do k = kbeg, kend
       do j = jbeg, jend
        do i = ibeg, iend
         km   (i,j,k) = kp   (i,j,k)
         dkedt(i,j,k) = dkpdt(i,j,k)
        enddo
       enddo
      enddo
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D  S U B R O U T I N E      //////////
c    //////////                  O P A C               \\\\\\\\\\
c
c=======================================================================
