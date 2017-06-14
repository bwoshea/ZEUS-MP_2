c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////            L S _ M U L _ B N D            \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c     PURPOSE:  compute a symmetric MATVEC
c
c     Written by: F. Douglas Swesty and John Hayes
c
c=======================================================================
      subroutine sym_mul_bnd(isx,iex,isy,iey,isz,iez,
     .                      dd, ddp1,
     .                          ddp2,
     .                          ddp3,
     .                      x,rhs)
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
      integer  :: isxm1, isxp1, iexm1, iexp1, isym1, isyp1, ieym1, 
     .            ieyp1, iszm1, iszp1, iezm1, iezp1
c
      real(rl) ::   dd(neqm,neqm,in,jn,kn)
      real(rl) :: ddp1(neqm,neqm,in,jn,kn)
      real(rl) :: ddp2(neqm,neqm,in,jn,kn)
      real(rl) :: ddp3(neqm,neqm,in,jn,kn)
      real(rl) :: x   (neqm,       in,jn,kn), 
     .            rhs (neqm,       in,jn,kn)
c
c                            local variables
c
      real(rl) :: sum
c
c                            loop indices
c
      integer  :: i, jx, jy, jz, k
c
c     x-faces
c
      do jz = isz,iez,1
       do jx = isx,iex,1
        do i = 1,neqm,1
         sum = 0.0d0
         do k = 1,neqm,1
          sum = sum + 
     .                ddp1(i,k,jx-1,isy,jz) * x(k,jx-1,isy  ,jz  )+
     .                  dd(i,k,jx  ,isy,jz) * x(k,jx  ,isy  ,jz  )+
     .                ddp1(i,k,jx  ,isy,jz) * x(k,jx+1,isy  ,jz  )
     .               +ddp2(i,k,jx,isy-1,jz) * x(k,jx  ,isy-1,jz  )
     .               +ddp2(i,k,jx,isy  ,jz) * x(k,jx  ,isy+1,jz  )
     .               +ddp3(i,k,jx,isy,jz-1) * x(k,jx  ,isy  ,jz-1)
     .               +ddp3(i,k,jx,isy,jz  ) * x(k,jx  ,isy  ,jz+1)
         enddo
         rhs(i,jx,isy,jz) = sum
        enddo
        do i = 1,neqm,1
         sum = 0.0d0
         do k = 1,neqm,1
          sum = sum + 
     .                ddp1(i,k,jx-1,iey,jz) * x(k,jx-1,iey  ,jz  )+
     .                  dd(i,k,jx  ,iey,jz) * x(k,jx  ,iey  ,jz  )+
     .                ddp1(i,k,jx  ,iey,jz) * x(k,jx+1,iey  ,jz  )
     .               +ddp2(i,k,jx,iey-1,jz) * x(k,jx  ,iey-1,jz  )
     .               +ddp2(i,k,jx,iey  ,jz) * x(k,jx  ,iey+1,jz  )
     .               +ddp3(i,k,jx,iey,jz-1) * x(k,jx  ,iey  ,jz-1)
     .               +ddp3(i,k,jx,iey,jz  ) * x(k,jx  ,iey  ,jz+1)
         enddo
         rhs(i,jx,iey,jz) = sum
        enddo
       enddo
      enddo
      if(ldimen .gt. 1) then
c
c     y-faces
c
       do jz = isz,iez,1
        do jy = isy+1,iey-1,1
         do i = 1,neqm,1
          sum = 0.0d0
          do k = 1,neqm,1
           sum = sum + 
     .                 ddp1(i,k,isx-1,jy,jz) * x(k,isx-1,jy  ,jz  )+
     .                   dd(i,k,isx  ,jy,jz) * x(k,isx  ,jy  ,jz  )+
     .                 ddp1(i,k,isx  ,jy,jz) * x(k,isx+1,jy  ,jz  )+
     .                 ddp2(i,k,isx,jy-1,jz) * x(k,isx  ,jy-1,jz  )+
     .                 ddp2(i,k,isx,jy  ,jz) * x(k,isx  ,jy+1,jz  )
     .                +ddp3(i,k,isx,jy,jz-1) * x(k,isx  ,jy  ,jz-1)
     .                +ddp3(i,k,isx,jy,jz  ) * x(k,isx  ,jy  ,jz+1)
          enddo
          rhs(i,isx,jy,jz) = sum
         enddo
         do i = 1,neqm,1
          sum = 0.0d0
          do k = 1,neqm,1
           sum = sum + 
     .                 ddp1(i,k,iex-1,jy,jz) * x(k,iex-1,jy  ,jz  )+
     .                   dd(i,k,iex  ,jy,jz) * x(k,iex  ,jy  ,jz  )+
     .                 ddp1(i,k,iex  ,jy,jz) * x(k,iex+1,jy  ,jz  )+
     .                 ddp2(i,k,iex,jy-1,jz) * x(k,iex  ,jy-1,jz  )+
     .                 ddp2(i,k,iex,jy  ,jz) * x(k,iex  ,jy+1,jz  )
     .                +ddp3(i,k,iex,jy,jz-1) * x(k,iex  ,jy  ,jz-1)
     .                +ddp3(i,k,iex,jy,jz  ) * x(k,iex  ,jy  ,jz+1)
          enddo
          rhs(i,iex,jy,jz) = sum
         enddo
        enddo
       enddo
       if(ldimen .gt. 2) then
c
c     z-faces
c
        do jx = isx+1,iex-1,1
         do jy = isy+1,iey-1,1
          do i = 1,neqm,1
           sum = 0.0d0
           do k = 1,neqm,1
            sum = sum + ddp3(i,k,jx,jy,isz-1) * x(k,jx  ,jy  ,isz-1)+
     .                  ddp2(i,k,jx,jy-1,isz) * x(k,jx  ,jy-1,isz  )+
     &                  ddp1(i,k,jx-1,jy,isz) * x(k,jx-1,jy  ,isz  )+
     &                    dd(i,k,jx,jy,isz) * x(k,jx  ,jy  ,isz  )+
     &                  ddp1(i,k,jx,jy,isz) * x(k,jx+1,jy  ,isz  )+
     &                  ddp2(i,k,jx,jy,isz) * x(k,jx  ,jy+1,isz  )+
     .                  ddp3(i,k,jx,jy,isz) * x(k,jx  ,jy  ,isz+1)
           enddo
           rhs(i,jx,jy,isz) = sum
          enddo
          do i = 1,neqm,1
           sum = 0.0d0
           do k = 1,neqm,1
            sum = sum + ddp3(i,k,jx,jy,iez-1) * x(k,jx  ,jy  ,iez-1)+
     .                  ddp2(i,k,jx,jy-1,iez) * x(k,jx  ,jy-1,iez  )+
     &                  ddp1(i,k,jx-1,jy,iez) * x(k,jx-1,jy  ,iez  )+
     &                    dd(i,k,jx,jy,iez) * x(k,jx  ,jy  ,iez  )+
     &                  ddp1(i,k,jx,jy,iez) * x(k,jx+1,jy  ,iez  )+
     &                  ddp2(i,k,jx,jy,iez) * x(k,jx  ,jy+1,iez  )+
     .                  ddp3(i,k,jx,jy,iez) * x(k,jx  ,jy  ,iez+1)
           enddo
           rhs(i,jx,jy,iez) = sum
          enddo
         enddo
        enddo
       endif ! ldimen > 2
      endif ! ldimen > 1
c
c                        zero the boundary elements
      if(ldimen .eq. 1) then
       isym1 = isy
       isyp1 = isy
       ieym1 = isy
       ieyp1 = isy
       iszm1 = isz
       iszp1 = isz
       iezm1 = isz
       iezp1 = isz
      endif
      if(ldimen .eq. 2) then
       isym1 = isy-1
       isyp1 = isy+1
       ieym1 = iey-1
       ieyp1 = iey+1
       iszm1 = isz
       iszp1 = isz
       iezm1 = isz
       iezp1 = isz
      endif
      if(ldimen .eq. 3) then
       isym1 = isy-1
       isyp1 = isy+1
       ieym1 = iey-1
       ieyp1 = iey+1
       iszm1 = isz-1
       iszp1 = isz+1
       iezm1 = iez-1
       iezp1 = iez+1
      endif
      do jz = iszm1, iezp1
       do jy = isy, iey
        do k = 1,neqm,1
           rhs(k,isx-1,jy,jz) = 0.0D0
           rhs(k,iex+1,jy,jz) = 0.0D0
        enddo
       enddo
      enddo
      if(ldimen .gt. 1) then
c
       do jz = iszm1, iezp1
        do jx = isx-1, iex+1
         do k = 1,neqm,1
           rhs(k,jx,isym1,jz) = 0.0d0
           rhs(k,jx,ieyp1,jz) = 0.0d0
         enddo
        enddo
       enddo
       if(ldimen .gt. 2) then
c
        do jy = isy, iey
         do jx = isx, iex
          do k = 1,neqm,1
             rhs(k,jx,jy,iszm1) = 0.0D0
             rhs(k,jx,jy,iezp1) = 0.0D0
          enddo
         enddo
        enddo
       endif ! ldimen > 2
      endif ! ldimen > 1
c
 999  return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////            L S _ M U L _ B N D            \\\\\\\\\\
c
c=======================================================================
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////            L S _ M U L _ I N T            \\\\\\\\\\
c
c=======================================================================
      subroutine sym_mul_int(isx,iex,isy,iey,isz,iez,
     .                      dd, ddp1,
     .                          ddp2,
     .                          ddp3,
     .                      x,rhs)
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
      real(rl) ::   dd(neqm,neqm,in,jn,kn)
      real(rl) :: ddp1(neqm,neqm,in,jn,kn)
      real(rl) :: ddp2(neqm,neqm,in,jn,kn)
      real(rl) :: ddp3(neqm,neqm,in,jn,kn)
      real(rl) :: x   (neqm,       in,jn,kn), 
     .            rhs (neqm,       in,jn,kn)
c
c                            local variables
c
      real(rl) :: sum
c
c                            loop indices
c
      integer  :: i, jx, jy, jz, k
c
      integer  :: isxm1, isxp1, iexm1, iexp1, isym1, isyp1, ieym1, 
     .            ieyp1, iszm1, iszp1, iezm1, iezp1
c
      if(ldimen .eq. 1) then
       isym1 = isy
       isyp1 = isy
       ieym1 = isy
       ieyp1 = isy
       iszm1 = isz
       iszp1 = isz
       iezm1 = isz
       iezp1 = isz
      endif
      if(ldimen .eq. 2) then
       isym1 = isy-1
       isyp1 = isy+1
       ieym1 = iey-1
       ieyp1 = iey+1
       iszm1 = isz
       iszp1 = isz
       iezm1 = isz
       iezp1 = isz
      endif
      if(ldimen .eq. 3) then
       isym1 = isy-1
       isyp1 = isy+1
       ieym1 = iey-1
       ieyp1 = iey+1
       iszm1 = isz-1
       iszp1 = isz+1
       iezm1 = iez-1
       iezp1 = iez+1
      endif 
c
      do jz = iszp1,iezm1,1
       do jy = isyp1,ieym1,1
        do jx = isx+1,iex-1,1
         do i = 1,neqm,1
          sum = 0.0d0
          do k = 1,neqm,1
           sum = sum + 
     .                 ddp1(i,k,jx-1,jy,jz) * x(k,jx-1,jy  ,jz  )+
     .                   dd(i,k,jx  ,jy,jz) * x(k,jx  ,jy  ,jz  )+
     .                 ddp1(i,k,jx  ,jy,jz) * x(k,jx+1,jy  ,jz  )
     .                +ddp2(i,k,jx,jy-1,jz) * x(k,jx  ,jy-1,jz  )
     .                +ddp2(i,k,jx,jy  ,jz) * x(k,jx  ,jy+1,jz  )
     .                +ddp3(i,k,jx,jy,jz-1) * x(k,jx  ,jy  ,jz-1)
     .                +ddp3(i,k,jx,jy,jz  ) * x(k,jx  ,jy  ,jz+1)
          enddo
          rhs(i,jx,jy,jz) = sum
         enddo
        enddo
       enddo
      enddo
c
 999  return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////            L S _ M U L _ I N T            \\\\\\\\\\
c
c=======================================================================
