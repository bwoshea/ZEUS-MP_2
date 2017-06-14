c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////               D I F C O 2                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c     PURPOSE: compute diffusion coefficients for 2-component of
c              radiation flux
c
c     Written by: Robert Fiedler and John Hayes
c
c=======================================================================
c
c=======================================================================
c
      subroutine difco2(erold,dfc,ibeg,iend,jbeg,jend,kbeg,kend)
c
      use real_prec
      use config
      use param
      use grid
      use field
      use root
      use scratch
      use radiation
      use opac
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use cons
c
      implicit NONE
c
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend
c
      real(rl) :: der1, der2, der3, lmda, qa, qb, grad_er, ros_mfp
      real(rl) :: r(in,jn,kn), dfc(in,jn,kn), erold(in,jn,kn)
      real(rl) :: onethird, rdmcm, rdmcp, flux1, flux2
c
#ifndef MARSHAK
c
c=======================================================================
c  Compute radiation energy gradient vector and R
c=======================================================================
c
      if(ldimen .eq. 1) then
       do k = kbeg, kend
        do j = jbeg, jend
         do i = ibeg, iend
          dfc(i,j,k) = 0.0D0
         enddo
        enddo
       enddo
       return
      endif ! ldimen = 1
      do k = kbeg, kend
       do j = jbeg, jend
        do i = ibeg, iend
         der1     = (erold(i  ,j-1,k) - erold(i-1,j-1,k  )) * dx1bi(i)
     .            + (erold(i  ,j  ,k) - erold(i-1,j  ,k  )) * dx1bi(i)
         der1     = der1 * 0.5 
         der2     = (erold(i  ,j  ,k) - erold(i  ,j-1,k  )) * dx2bi(j)
     .            * g2bi(i)
         if(ldimen .eq. 2) then
          der3     = 0.0
         else
          der3     = (erold(i  ,j-1,k) - erold(i  ,j-1,k-1)) * dx3bi(k)
     .             + (erold(i  ,j  ,k) - erold(i  ,j  ,k-1)) * dx3bi(k)
          der3     = der3 * 0.5 * g31bi(i) * g32ai(j)
         endif
         grad_er  = sqrt(der1**2 + der2**2 + der3**2)
         rdmcm    = kapr(i,j-1,k)
         rdmcp    = kapr(i,j  ,k)
         ros_mfp  = 2.0 / (rdmcm + rdmcp)
         r(i,j,k) = 2.0 * grad_er * ros_mfp / 
     .              max(erold(i,j,k)+erold(i,j-1,k), tiny)
c
c  Compute FLD constant and components of Eddington tensor for 
c  Chapman-Enskog theory (ifld=1)
c
         if (ifld .eq. 1) then
          lmda       = (2.0+r(i,j,k)) / (6.0+3.0*r(i,j,k)+r(i,j,k)**2)
         endif
c
c  Compute FLD constant and components of Eddington tensor for
c  piecewise-linear Minerbo theory (ifld=2)
c
         if (ifld .eq. 2) then
	  qa   = 2.0/(3.0 + sqrt(9.0+12.0*r(i,j,k)**2))
	  qb   = 1.0/(1.0 + r(i,j,k) + sqrt(1.0+2.0*r(i,j,k)))
          if (r(i,j,k) .le. 1.5) then
            lmda = qa
          else
            lmda = qb
          endif
         endif
         dfc(i,j,k) = clight * lmda * ros_mfp
c
        enddo ! i
       enddo ! j
      enddo ! k
c
#else /* MARSHAK */
      onethird = 1. / 3.
      do k = kbeg, kend
       do j = jbeg, jend
        do i = ibeg, iend
         rdmcm    = kapr(i,j-1,k)
         rdmcp    = kapr(i,j  ,k)
         ros_mfp  = 2.0 / (rdmcm + rdmcp)
         dfc(i,j,k) = clight * onethird * ros_mfp
        enddo
       enddo
      enddo
#endif /* MARSHAK */
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////               D I F C O 2                 \\\\\\\\\\
c
c=======================================================================
c
