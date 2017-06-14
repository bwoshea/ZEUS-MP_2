c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                   S O D                   \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c=======================================================================
c
       subroutine sod
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
      real(rl) :: x10, x20, x30,
     .            d0 , p0 , e0 ,
     .            d1 , p1 , e1
c
      integer  :: iin (ijkn), iout(ijkn), jin (ijkn),
     .            jout(ijkn), kin (ijkn), kout(ijkn)
c
      real(rl) :: massk(ijkn)
c
      namelist / pgen     /
     .               x10, x20, p0, d0, p1, d1, idirect, x30
c
c-----------------------------------------------------------------------
c
#ifdef IGNORE
       do k = 1, kn
       do j = 1, jn
       do i = 1, in
        v1(i,j,k) = 1.0D10
       enddo
       enddo
       enddo
#endif
       x10  = 0.5D0
       x20  = 0.5D0
       x30  = 0.5D0
       d0   = 1.0D0
       p0   = 1.0D0
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
         ibuf_in( 1) = idirect
       endif
       call MPI_BCAST( buf_in, 6, MPI_FLOAT
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
         idirect = ibuf_in( 1)
#endif /* MPI_USED */
       endif
c
c      Set up tube.
c
      e0 = p0 / gamm1
      e1 = p1 / gamm1
      do k=1,kn
       do j=1,jn
        do i=1,in
         v1(i,j,k) = 0.0D0
         v2(i,j,k) = 0.0D0
         v3(i,j,k) = 0.0D0
        enddo
       enddo
      enddo
      if(.false.) then
      do k = 1, kn
       do j = 1, jn
        v1oib(j,k,1) = -1.0D0
        v1oib(j,k,2) = -1.0D0
        d oib(j,k,1) = d1
        d oib(j,k,2) = d1
        e oib(j,k,1) = e1
        e oib(j,k,2) = e1
       enddo
      enddo
      endif
      if(.false.) then
      do k = 1, kn
       do i = 1, in
        v2ojb(i,k,1) = -1.0D0
        v2ojb(i,k,2) = -1.0D0
        d ojb(i,k,1) = d1
        d ojb(i,k,2) = d1
        e ojb(i,k,1) = e1
        e ojb(i,k,2) = e1
       enddo
      enddo
      endif
      if(.false.) then
      do j = 1, jn
       do i = 1, in
        v3okb(i,j,1) = -1.0D0
        v3okb(i,j,2) = -1.0D0
        d okb(i,j,1) = d1
        d okb(i,j,2) = d1
        e okb(i,j,1) = e1
        e okb(i,j,2) = e1
       enddo
      enddo
      endif
c
      if(idirect .eq. 1) then
       do k = 1, kn
        do j = 1, jn
         do i = 1, in
          if(x1a(i) .le. x10) then
           e(i,j,k) = e0
           d(i,j,k) = d0
          else
           e(i,j,k) = e1
           d(i,j,k) = d1
          endif
         enddo
        enddo
       enddo
      endif ! idirect
      if(idirect .eq. 2) then
       do k = 1, kn
        do j = 1, jn
         do i = 1, in
          if(x2a(j) .le. x20) then
           e(i,j,k) = e0
           d(i,j,k) = d0
          else
           e(i,j,k) = e1
           d(i,j,k) = d1
          endif
         enddo
        enddo
       enddo
      endif ! idirect
      if(idirect .eq. 3) then
       do k = 1, kn
        do j = 1, jn
         do i = 1, in
          if(x3a(k) .le. x30) then
           e(i,j,k) = e0
           d(i,j,k) = d0
          else
           e(i,j,k) = e1
           d(i,j,k) = d1
          endif
         enddo
        enddo
       enddo
      endif ! idirect
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                   S O D                   \\\\\\\\\\
c
c=======================================================================
c
c
      subroutine sodres
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
      integer  :: i, j, k, ip1, jp1, kp1, idirect
      real(rl) :: x10, x20, x30,
     .            d0 , p0 , e0 ,
     .            d1 , p1 , e1
c
      integer  :: iin (ijkn), iout(ijkn), jin (ijkn),
     .            jout(ijkn), kin (ijkn), kout(ijkn)
      real(rl) :: massk(ijkn)
c
      namelist / pgen     /
     .               x10, x20, p0, d0, p1, d1, idirect
c
c-----------------------------------------------------------------------
c
       x10  = 0.5D0
       x20  = 0.5D0
       x30  = 0.5D0
       d0   = 1.0D0
       p0   = 1.0D0
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
         ibuf_in( 1) = idirect
       endif
       call MPI_BCAST( buf_in, 6, MPI_FLOAT
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
         idirect = ibuf_in( 1)
#endif /* MPI_USED */
       endif
c
       return
       end
