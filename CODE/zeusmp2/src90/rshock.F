c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                R S H O C K                \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c=======================================================================
      subroutine rshock
c
      use real_prec
      use config
      use param
      use cons
      use grid
      use field
      use radiation
      use opac
      use root
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer  :: i, j, k
c
      real(rl) :: d0, e0, temp, dtemp
c
      namelist /pgen/ d0, e0
c
c     initialize and read in parameters from PGEN namelist
c
      d0 = 1.D0
      e0 = 1.D-4
c
      if (myid_w .eq. 0) then
       read (1,pgen)
       write(2,pgen)
#ifdef MPI_USED
       buf_in(1) = d0
       buf_in(2) = e0
      endif
      call MPI_BCAST( buf_in, 2, MPI_FLOAT
     &               , 0, MPI_COMM_WORLD, ierr )
      if(myid_w .ne. 0) then
       d0 = buf_in(1)
       e0 = buf_in(2)
#endif /* MPI_USED */
      endif
c
      gamm1 = gamma - 1.0D0
c
c     set up density distribution
c
      do 1 k = 1, kn
      do 1 j = 1, jn
      do 1 i = 1, in
       d(i,j,k) = d0
       p(i,j,k) = -1.0
1     continue
c
c     initialize gas and radiation energy densities, velocities
c
      fois( 2) = e0*d0*boltz/((gamma-1.0)*mmw*mh)
      temp     = e0
      fois(12) = rad_con * temp**4
      dtemp = -75.D0/7.D10
      do 2 k=1,kn
      do 2 j=1,jn
      do 2 i=is,ie
       temp      = e0 + dtemp*(x1b(i)-7.0D10)
       e (i,j,k) = temp*d0*boltz/((gamma-1.0)*mmw*mh)
       if(lrad .eq. 1) then
        er(i,j,k) = rad_con * temp**4
       endif
       v1(i,j,k) = fois(3)
       v2(i,j,k) = 0.D0
       v3(i,j,k) = 0.D0
2     continue
c
c     set up outer 1-boundary values
c
      do k = 1, kn
      do j = 1, jn
       e oib(j,k,1) = fois( 2)
       e oib(j,k,2) = fois( 2)
       if(lrad .eq. 1) then
        eroib(j,k,1) = fois(12)
        eroib(j,k,2) = fois(12)
       endif
      enddo
      enddo
c
      return
      end
      subroutine rshockres
c
      use real_prec
      use config
      use param
      use cons
      use grid
      use field
      use radiation
      use opac
      use root
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer  :: i, j, k
c
      real(rl) :: d0, e0, temp, dtemp
c
      namelist /pgen/ d0, e0
c
c     initialize and read in parameters from PGEN namelist
c
      d0 = 1.D0
      e0 = 1.D-4
c
      if (myid_w .eq. 0) then
       read (1,pgen)
       write(2,pgen)
#ifdef MPI_USED
       buf_in(1) = d0
       buf_in(2) = e0
      endif
      call MPI_BCAST( buf_in, 2, MPI_FLOAT
     &               , 0, MPI_COMM_WORLD, ierr )
      if(myid_w .ne. 0) then
       d0 = buf_in(1)
       e0 = buf_in(2)
#endif /* MPI_USED */
      endif
c
      gamm1 = gamma - 1.0D0
c
      return
      end
