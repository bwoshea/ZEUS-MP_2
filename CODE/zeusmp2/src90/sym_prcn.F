c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////           L S _ P R C N _ B N D           \\\\\\\\\\
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
      subroutine sym_prcn_bnd(isx,iex,isy,iey,isz,iez,
     .                       ipcflag,itrn,
     .                       d, dp1,
     .                          dp2,
     .                          dp3,
     .                       x,      rhs)
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
      integer :: isx, iex, isy, iey, isz, iez,
     .           ipcflag, itrn
c
      real(rl) ::   d(neqm,neqm,in,jn,kn)
      real(rl) :: dp1(neqm,neqm,in,jn,kn)
      real(rl) :: dp2(neqm,neqm,in,jn,kn)
      real(rl) :: dp3(neqm,neqm,in,jn,kn)
      real(rl) :: x  (neqm,     in,jn,kn)
      real(rl) :: rhs(neqm,     in,jn,kn)
c
c                            loop indices
c
      integer :: i, j, k
c
c                            diagonal scale factor
c
      real(rl) :: atil
c
      if(ipcflag.eq.2 .or. ipcflag .eq. 0) then
        call sym_diag_bnd(isx,iex,isy,iey,isz,iez,ipcflag,itrn,
     .                    d, dp1,dp2,dp3,x,rhs)
      else
        write(*,*)  ' ls_prcn: precon. flag ',ipcflag, 'not supported'
        write(*,*)  ' ls_prcn: i am stopping '
        stop
      endif
c
 999  return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////           L S _ P R C N _ B N D           \\\\\\\\\\
c
c=======================================================================
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////           L S _ P R C N _ I N T           \\\\\\\\\\
c
c=======================================================================
      subroutine sym_prcn_int(isx,iex,isy,iey,isz,iez,
     .                       ipcflag,itrn,
     .                       d, dp1,
     .                          dp2,
     .                          dp3,
     .                       x,      rhs)
c
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
     .            ipcflag, itrn
c
      real(rl) ::   d(neqm,neqm,in,jn,kn)
      real(rl) :: dp1(neqm,neqm,in,jn,kn)
      real(rl) :: dp2(neqm,neqm,in,jn,kn)
      real(rl) :: dp3(neqm,neqm,in,jn,kn)
      real(rl) :: x  (neqm,     in,jn,kn)
      real(rl) :: rhs(neqm,     in,jn,kn)
c
c                            loop indices
c
      integer  :: i, j, k
c
c                            diagonal scale factor
c
      real(rl) :: atil
c
      if(ipcflag .eq. 0 .or. ipcflag.eq.2) then
        call sym_diag_int(isx,iex,isy,iey,isz,iez,ipcflag,itrn,
     .                   d, dp1,dp2,dp3,x,rhs)
      else
        write(*,*)  ' ls_prcn: precon. flag ',ipcflag, 'not supported'
        write(*,*)  ' ls_prcn: i am stopping '
        stop
      endif
c
 999  return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////           L S _ P R C N _ I N T           \\\\\\\\\\
c
c=======================================================================
