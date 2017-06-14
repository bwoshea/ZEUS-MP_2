c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////           L S _ D I A G _ B N D           \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c     PURPOSE:  compute diagonal preconditioning
c
c     Written by: F. Douglas Swesty and John Hayes
c
c=======================================================================
      subroutine sym_diag_bnd(isx,iex,isy,iey,isz,iez,
     .                        ipcflag,itrans,
     .                        dd, ddp1, ddp2,ddp3,x,rhs)
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
      integer  :: isx, iex, isy, iey, isz, iez, 
     .            itrans, ipcflag
c
      real(rl) ::   dd(neqm,neqm,in,jn,kn)
      real(rl) :: ddp1(neqm,neqm,in,jn,kn)
      real(rl) :: ddp2(neqm,neqm,in,jn,kn)
      real(rl) :: ddp3(neqm,neqm,in,jn,kn)
      real(rl) :: x   (neqm      ,in,jn,kn)
      real(rl) :: rhs (neqm      ,in,jn,kn)
c
c                            loop indices
c
      integer  :: i, jx, jy, jz, k
c
      integer  :: isxm1, isxp1, iexm1, iexp1, isym1, isyp1, ieym1, 
     .            ieyp1,
     .            iszm1, iszp1, iezm1, iezp1
c
c                            diagonal scale factor
      real(rl) :: atil
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
      do jz = isz,iez,1
       do jx = isx,iex,1
        do i = 1,neqm,1
         if(ipcflag .eq. 0) then
          atil = 1.0D0
         endif
         if(ipcflag .eq. 2) then
         if( dabs(dd(i,i,jx,isy,jz)) .gt. 1.0d-60 ) then
          atil = dd(i,i,jx,isy,jz)
         else
          atil = sign( 1.0d0 , dd(i,i,jx,isy,jz) )
         endif
         x(i,jx,isy,jz) = rhs(i,jx,isy,jz)/atil
         if( dabs(dd(i,i,jx,iey,jz)) .gt. 1.0d-60 ) then
          atil = dd(i,i,jx,iey,jz)
         else
          atil = sign( 1.0d0 , dd(i,i,jx,iey,jz) )
         endif
         x(i,jx,iey,jz) = rhs(i,jx,iey,jz)/atil
         endif ! ipcflag
        enddo
       enddo
       if(ldimen .gt. 1) then
        do jy = isy+1,iey-1,1
         do i = 1,neqm,1
          if(ipcflag .eq. 0) then
           atil = 1.0D0
          endif
          if(ipcflag .eq. 2) then
          if( dabs(dd(i,i,isx,jy,jz)) .gt. 1.0d-60 ) then
           atil = dd(i,i,isx,jy,jz)
          else
           atil = sign( 1.0d0 , dd(i,i,isx,jy,jz) )
          endif
          x(i,isx,jy,jz) = rhs(i,isx,jy,jz)/atil
          if( dabs(dd(i,i,iex,jy,jz)) .gt. 1.0d-60 ) then
           atil = dd(i,i,iex,jy,jz)
          else
           atil = sign( 1.0d0 , dd(i,i,iex,jy,jz) )
          endif
          endif ! ipcflag
          x(i,iex,jy,jz) = rhs(i,iex,jy,jz)/atil
         enddo
        enddo
       endif ! ldimen > 1
      enddo
      if(ldimen .eq. 3) then
       do jx = isx,iex,1
        do jy = isy,iey,1
         do i = 1,neqm,1
          if(ipcflag .eq. 0) then
           atil = 1.0d0
          endif
          if(ipcflag .eq. 2) then
          if( dabs(dd(i,i,jx,jy,isz)) .gt. 1.0d-60 ) then
           atil = dd(i,i,jx,jy,isz)
          else
           atil = sign( 1.0d0 , dd(i,i,jx,jy,isz) )
          endif
          x(i,jx,jy,isz) = rhs(i,jx,jy,isz)/atil
          if( dabs(dd(i,i,jx,jy,iez)) .gt. 1.0d-60 ) then
           atil = dd(i,i,jx,jy,iez)
          else
           atil = sign( 1.0d0 , dd(i,i,jx,jy,iez) )
          endif
          endif ! ipcflag
          x(i,jx,jy,iez) = rhs(i,jx,jy,iez)/atil
         enddo
        enddo
       enddo
      endif ! ldimen = 3
c                         zero the boundary elements
      do jz = iszm1, iezp1
       do jy = isy, iey
        do k = 1,neqm,1
         x(k,isx-1,jy,jz) = 0.0d0
         x(k,iex+1,jy,jz) = 0.0d0
        enddo
       enddo
      enddo
c
      if(ldimen .gt. 1) then
       do jz = iszm1, iezp1
        do jx = isx-1, iex+1
         do k = 1,neqm,1
          x(k,jx,isy-1,jz) = 0.0d0
          x(k,jx,iey+1,jz) = 0.0d0
         enddo
        enddo
       enddo
       if(ldimen .gt. 2) then
        do jx = isx, iex
         do jy = isy, iey
          do k = 1,neqm,1
           x(k,jx,jy,isz-1) = 0.0D0
           x(k,jx,jy,iez+1) = 0.0D0
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
c    //////////           L S _ D I A G _ B N D           \\\\\\\\\\
c
c=======================================================================
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////           L S _ D I A G _ I N T           \\\\\\\\\\
c
c=======================================================================
      subroutine sym_diag_int(isx,iex,isy,iey,isz,iez,ipcflag,itrans,
     .                        dd, ddp1, ddp2, ddp3, x, rhs)
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
      integer  :: isx, iex, isy, iey, isz, iez, 
     .            itrans, ipcflag
c
      real(rl) ::   dd(neqm,neqm,in,jn,kn)
      real(rl) :: ddp1(neqm,neqm,in,jn,kn)
      real(rl) :: ddp2(neqm,neqm,in,jn,kn)
      real(rl) :: ddp3(neqm,neqm,in,jn,kn)
      real(rl) :: x   (neqm,     in,jn,kn)
      real(rl) :: rhs (neqm,     in,jn,kn)
c
c                            loop indices
c
      integer  :: i, jx, jy, jz, k
c
      integer  :: isxm1, isxp1, iexm1, iexp1, isym1, isyp1, ieym1, 
     .            ieyp1,
     .            iszm1, iszp1, iezm1, iezp1
c
c                            diagonal scale factor
      real(rl) :: atil
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
          if(ipcflag .eq. 0) then
           atil = 1.0d0
          endif
          if(ipcflag .eq. 2) then
          if( dabs(dd(i,i,jx,jy,jz)) .gt. 1.0d-60 ) then
           atil = dd(i,i,jx,jy,jz)
          else
           atil = sign( 1.0d0 , dd(i,i,jx,jy,jz) )
          endif
          endif ! ipcflag
          x(i,jx,jy,jz) = rhs(i,jx,jy,jz)/atil
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
c    //////////           L S _ D I A G _ I N T           \\\\\\\\\\
c
c=======================================================================
