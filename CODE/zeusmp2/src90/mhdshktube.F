c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////            M H D S H K T U B E            \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c      PURPOSE: initializes Brio & Wu MHD shock tube problem.
c
c      Written by: John Hayes
c
c=======================================================================
c
       subroutine mhdshktube
c
      use real_prec
      use param
      use config
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
      real(rl) :: x10, x20, x30,
     .            d0 , p0 , e0 , b1_0, b2_0,
     .            d1 , p1 , e1 , b1_1, b2_1, b3_0, b3_1, line, bt_0,
     .            bt_2
c
      namelist / pgen     /
     .               x10, x20, p0, d0, p1, d1, idirect, x30,
     .               b1_0, b2_0, b1_1, b2_1, b3_0, b3_1
c
c----------------------------------------------------------------------- 
c 
       x10  = 400.0D0
       x20  = 0.5D0
       x30  = 0.5D0
       d0   = 1.0D0
       p0   = 1.0D0
       b1_0 = 0.75D0
       b2_0 = 1.00D0
       b3_0 = 0.00D0
       b1_1 = 0.75D0
       b2_1 =-1.00D0
       b3_1 = 0.00D0
       d1   = 0.125D0
       p1   = 0.1D0
c
       if (myid .eq. 0) then
         read (1, pgen)
         write (2, pgen)
#ifdef MPI_USED
         buf_in(1) = x10 
         buf_in(2) = x20 
         buf_in(3) = d0  
         buf_in(4) = p0  
         buf_in(5) = d1  
         buf_in(6) = p1  
         buf_in(7) = b1_0
         buf_in(8) = b1_1
         buf_in(9) = b2_0
         buf_in(10) = b2_1
         buf_in(11) = b3_0
         buf_in(12) = b3_1
         ibuf_in( 1) = idirect
       endif
       call MPI_BCAST( buf_in, 12, MPI_FLOAT
     &               , 0, comm3d, ierr )
       call MPI_BCAST( ibuf_in, 1, MPI_INTEGER
     &               , 0, comm3d, ierr )
       if (myid .ne. 0) then
         x10  = buf_in(1)
         x20  = buf_in(2)
         d0   = buf_in(3)
         p0   = buf_in(4)
         d1   = buf_in(5)
         p1   = buf_in(6)
         b1_0 = buf_in(7)
         b1_1 = buf_in(8)
         b2_0 = buf_in(9)
         b2_1 = buf_in(10)
         b3_0 = buf_in(11)
         b3_1 = buf_in(12)
         idirect = ibuf_in( 1)
#endif /* MPI_USED */
       endif
c
c      Set up tube.
c
      e0 = p0 / gamm1
      e1 = p1 / gamm1
c
      do k=1,kn
       do j=1,jn
        do i=1,in
         v1(i,j,k) = 0.0D0
         v2(i,j,k) = 0.0D0
         v3(i,j,k) = 0.0D0
        enddo
       enddo
      enddo
c
      do k = 1, kn
       do j = 1, jn
        do i = 1, in
         if(idirect .eq. 1) then
          if(x1a(i) .le. x10) then
           e (i,j,k) = e0
           d (i,j,k) = d0
           b1(i,j,k) = b1_0
           b2(i,j,k) = b2_0
           b3(i,j,k) = b3_0
          else
           e (i,j,k) = e1
           d (i,j,k) = d1
           b1(i,j,k) = b1_1
           b2(i,j,k) = b2_1
           b3(i,j,k) = b3_1
          endif
         endif ! idirect
         if(idirect .eq. 2) then
          if(x2a(j) .le. x20) then
           e (i,j,k) = e0
           d (i,j,k) = d0
           b1(i,j,k) = b1_0
           b2(i,j,k) = b2_0
           b3(i,j,k) = b3_0
          else
           e (i,j,k) = e1
           d (i,j,k) = d1
           b1(i,j,k) = b1_1
           b2(i,j,k) = b2_1
           b3(i,j,k) = b3_1
          endif
         endif ! idirect
         if(idirect .eq. 3) then
          if(x3a(k) .le. x30) then
           e (i,j,k) = e0
           d (i,j,k) = d0
           b1(i,j,k) = b1_0
           b2(i,j,k) = b2_0
           b3(i,j,k) = b3_0
          else
           e (i,j,k) = e1
           d (i,j,k) = d1
           b1(i,j,k) = b1_1
           b2(i,j,k) = b2_1
           b3(i,j,k) = b3_1
          endif
         endif ! idirect
        enddo
       enddo
      enddo
c
      return
      end
