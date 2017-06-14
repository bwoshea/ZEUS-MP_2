c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////               D I F C O 1                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                  University of California at San Diego
c
c     PURPOSE: compute diffusion coefficients for 1-component of
c              radiation flux
c
c     Written by: Robert Fiedler and John Hayes
c
c=======================================================================
c=======================================================================
c
      subroutine difco1(erold,dfc,ibeg,iend,jbeg,jend,kbeg,kend)
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
      REAL(rl) :: r(in,jn,kn), dfc(in,jn,kn), erold(in,jn,kn)
      REAL(rl) :: onethird, rdmcm, rdmcp, flux1, flux2
c
#ifndef MARSHAK
c
c=======================================================================
c  Compute radiation energy gradient vector and R
c=======================================================================
c
      do k = kbeg, kend
       do j = jbeg, jend
        do i = ibeg, iend
         der1     = (erold(i  ,j  ,k) - erold(i-1,j  ,k  )) * dx1bi(i)
         if(ldimen .eq. 1) then
          der2     = 0.D0
         else
          der2     = (erold(i-1,j  ,k) - erold(i-1,j-1,k  )) * dx2bi(j)
     .             + (erold(i  ,j  ,k) - erold(i  ,j-1,k  )) * dx2bi(j)
          der2     = der2 * 0.5 * g2ai(i)
         endif ! ldimen = 1
         if(ldimen .ne. 3) then
          der3     = 0.0
         else
          der3     = (erold(i-1,j  ,k) - erold(i-1,j  ,k-1)) * dx3bi(k)
     .             + (erold(i  ,j  ,k) - erold(i  ,j  ,k-1)) * dx3bi(k)
          der3     = der3 * 0.5 * g31ai(i) * g32bi(j)
         endif ! ldimen /= 3
         grad_er  = sqrt(der1**2 + der2**2 + der3**2)
         rdmcm    = kapr(i-1,j,k)
         rdmcp    = kapr(i  ,j,k)
         ros_mfp  = 2.0 / (rdmcm + rdmcp)
         r(i,j,k) = 2.0 * grad_er * ros_mfp / 
     .              max(erold(i,j,k)+erold(i-1,j,k), tiny)
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
         rdmcm    = kapr(i-1,j,k)
         rdmcp    = kapr(i  ,j,k)
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
c    //////////               D I F C O 1                 \\\\\\\\\\
c
c=======================================================================
c
