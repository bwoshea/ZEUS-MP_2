c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                   N O H                   \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c=======================================================================
c
       subroutine noh
c
      use real_prec
      use param
      use field
      use bndry
      use grid
      use root
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
      integer  :: i  , j, k, ip1, jp1,
     .            kp1, idirect
c
      real(rl) ::
     .            d0 , p0 , e0, v0
c
      namelist / pgen     /
     .               p0, d0, v0, idirect
c
c-----------------------------------------------------------------------
c
       d0 =  1.0D0
       p0 =  1.0D-6
       v0 = -1.0D0
c
       if (myid .eq. 0) then
         read (1, pgen)
         write (2, pgen)
#ifdef MPI_USED
         buf_in(1)   = d0
         buf_in(2)   = p0 
         buf_in(3)   = v0  
         ibuf_in( 1) = idirect
       endif
       call MPI_BCAST( buf_in, 3, MPI_FLOAT
     &               , 0, comm3d, ierr )
       call MPI_BCAST( ibuf_in, 1, MPI_INTEGER
     &               , 0, comm3d, ierr )
       if (myid .ne. 0) then
         d0      = buf_in(1)
         p0      = buf_in(2)
         v0      = buf_in(3)
         idirect = ibuf_in( 1)
#endif /* MPI_USED */
       endif
c
c      Set up initial state.
c
      do k=ks,ke
       do j=js,je
        do i=is,ie
         d(i,j,k) = d0
         e(i,j,k) = p0/gamm1
        enddo
       enddo
      enddo
c
      do k = 1, kn
       do j = 1, jn
        do i = 1, in
         if(idirect .eq. 1) then
          v1(i,j,k) = v0
          v2(i,j,k) = 0.0
          v3(i,j,k) = 0.0
         else if(idirect .eq. 2) then
          v2(i,j,k) = v0
          v1(i,j,k) = 0.0
          v3(i,j,k) = 0.0
         else
          v3(i,j,k) = v0
          v1(i,j,k) = 0.0
          v2(i,j,k) = 0.0
         endif
        enddo
       enddo
      enddo
c
      if(idirect .eq. 1) then
       do k = 1, kn
        do j = 1, jn
         v1oib(j,k,1) = v0
         v1oib(j,k,2) = v0
         v2oib(j,k,1) = 0.0
         v2oib(j,k,2) = 0.0
         v3oib(j,k,1) = 0.0
         v3oib(j,k,2) = 0.0
         d oib(j,k,1) = d0
         d oib(j,k,2) = d0
         e oib(j,k,1) = p0/gamm1
         e oib(j,k,2) = p0/gamm1
        enddo
       enddo
      else if(idirect .eq. 2) then
       do k = 1, kn
        do i = 1, in
         v2ojb(i,k,1) = v0
         v2ojb(i,k,2) = v0
         v1ojb(i,k,1) = 0.0
         v1ojb(i,k,2) = 0.0
         v3ojb(i,k,1) = 0.0
         v3ojb(i,k,2) = 0.0
         d ojb(i,k,1) = d0
         d ojb(i,k,2) = d0
         e ojb(i,k,1) = p0/gamm1
         e ojb(i,k,2) = p0/gamm1
        enddo
       enddo
      else
       do j = 1, jn
        do i = 1, jn
         v3okb(i,j,1) = v0
         v3okb(i,j,2) = v0
         v1okb(i,j,1) = 0.0
         v1okb(i,j,2) = 0.0
         v2okb(i,j,1) = 0.0
         v2okb(i,j,2) = 0.0
         d okb(i,j,1) = d0
         d okb(i,j,2) = d0
         e okb(i,j,1) = p0/gamm1
         e okb(i,j,2) = p0/gamm1
        enddo
       enddo
      endif
c
      return
      end
