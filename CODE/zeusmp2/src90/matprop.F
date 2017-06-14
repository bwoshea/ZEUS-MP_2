c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////               M A T P R O P               \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c=======================================================================
c
      subroutine matprop(e, er, d, gam, t, dtde, p, bb, dbbdt,
     .                   kr, kp, sg, dkpdt, km, dkedt, dpde,
     .                   ibeg, iend, jbeg, jend, 
     .                   kbeg, kend)
c
c     This routine is a conglomeration of the old TEMP, EOS, PLANCK,
c     ABSORP, and SCATT routines.
c
      use real_prec
      use config
      use param
      use root
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
      real(rl) :: e     (in,jn,kn), d   (in,jn,kn), t    (in,jn,kn), 
     .            p     (in,jn,kn), bb  (in,jn,kn), kp   (in,jn,kn),
     .            dtde  (in,jn,kn), dpde(in,jn,kn), dbbdt(in,jn,kn),
     .            dkpdt (in,jn,kn), sg  (in,jn,kn), er   (in,jn,kn),
     .            dkedt (in,jn,kn), km  (in,jn,kn), kr   (in,jn,kn)
c
      real(rl) :: gam, coef, temp, binv, gam1, power, stef
      real(rl) :: t_max, rdmcm1, rdmcp1, rdmcm2, rdmcp2, rdmcm3, rdmcp3,
     .            ros_mfp
      real(rl) :: rmfp0i, hundred
c
      real(rl) :: so_eps , tau    , eriibn, timarg
c
      common /soui/ so_eps, tau, eriibn
c
      hundred = 100.0d0
      gam1 = gam - 1.0
      stef = clight * rad_con / (4.D0 * pi)
      if(leos .eq. 1) then
       do 10 k = kbeg, kend
       do 10 j = jbeg, jend
       do 10 i = ibeg, iend
        p   (i,j,k) = gam1*e(i,j,k)
        dpde(i,j,k) = gam1
10     continue
      else
       if(myid .eq. 0) write(*,"('MATPROP: Illegal EOS option!!')")
#ifdef MPI_USED
       call mpi_finalize(ierr)
#endif
       stop
      endif
c
      binv = 1.D0 / boltz
      do 20 k = kbeg, kend
      do 20 j = jbeg, jend
      do 20 i = ibeg, iend
#ifndef MARSHAK
       t   (i,j,k) = gam1*mmw*mh*binv*e(i,j,k)/d(i,j,k)
       dtde(i,j,k) = gam1*mmw*mh*binv         /d(i,j,k)
#else
       t   (i,j,k) = sqrt(sqrt( so_eps*e(i,j,k)/rad_con ))
#endif
C       t   (i,j,k) = sqrt(sqrt( e(i,j,k)/rad_con ))
C       dtde(i,j,k) = 0.25D0*t(i,j,k)/e(i,j,k)
20    continue
c
      do 30 k = kbeg, kend
      do 30 j = jbeg, jend
      do 30 i = ibeg, iend
        dbbdt(i,j,k) = stef *                t(i,j,k)**3
        bb   (i,j,k) =       dbbdt(i,j,k) * t(i,j,k)
        dbbdt(i,j,k) = 4.0 * dbbdt(i,j,k)
30    continue
c
      rmfp0i = 1.D0 / rmfp0
c
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
         kp   (i,j,k) = kpfrac*kr(i,j,k)
         kp   (i,j,k) = min(hundred,kp(i,j,k))
#endif /* ICF */
         dkpdt(i,j,k) = -powr * kp(i,j,k) / t(i,j,k)
#else /* MARSHAK */
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
c    //////////               M A T P R O P            \\\\\\\\\\
c
c=======================================================================
