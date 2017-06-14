c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////               L S _ D P R D               \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c     PURPOSE:  compute a vector dot product
c
c     Written by: F. Douglas SwestyF. Douglas Swesty
c
c=======================================================================
      real*8 function ls_dprd(isx,iex,isy,iey,isz,iez,v1,v2)
c
      use real_prec
      use config
      use param
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit none
c
      integer  :: isx, iex, isy, iey, isz, iez
c
      real(rl) :: v1(neqm,in,jn,kn), 
     .            v2(neqm,in,jn,kn)
      real(rl) ::  psum, dotprd
c
      integer  :: ix, iy, iz, j, ierror
c
      psum = 0.0d0
      do iz = isz,iez,1
       do iy = isy,iey,1
        do ix = isx,iex,1
         do j = 1,neqm,1
           psum = psum+v1(j,ix,iy,iz)*v2(j,ix,iy,iz)
         enddo
        enddo
       enddo
      enddo
c
#ifdef MPI_USED
      call mpi_allreduce(psum,dotprd,1,MPI_FLOAT,mpi_sum,
     &                   comm3d,ierror) 
#endif
#ifndef MPI_USED
      dotprd = psum
#endif
c
      ls_dprd = dotprd
c
 999  return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D    S U B R O U T I N E      //////////
c    //////////                L S _ D P R D             \\\\\\\\\\
c
c=======================================================================
