c=======================================================================
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
      subroutine intchk(iswres,iswhdf,iswtsl,iswhst,iswusr)
c
c  PURPOSE:  Reads the buffer for valid interrupt messages and takes
c    appropriate action. Also checks stopping criteria, and sets ifsen=1
c    if stop condition is detected.
c
c  INPUT ARGUMENTS: none
c
c  OUTPUT ARGUMENTS: iswres,iswhdf,iswtsl,iswhst,iswusr=switches for restart,
c  hdf, history, and USER dumps; set to 1 if dump to be made in DATAIO
c
c  LAST MODIFIED: RAF 9/09/96 for ZEUS-MP
c  LAST MODIFIED: efh 04/15/99 including tslice
c-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use root
      use clockmod
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c IBM
#ifdef ARCH_IBM
#  define CHECKIN checkin_
#  define BCDFLT  bcdflt_
#else
#  define CHECKIN checkin
#  define BCDFLT  bcdflt
#endif
c
      implicit NONE
c
      real(rl4) :: cputime, wclock
c
      integer   :: iswres,iswhdf,iswtsl,iswhst,iswusr
c
      integer   :: i,nchar,istrt,iend
      real(rl)  :: valnew
      character*80 :: msg
      integer  :: maxwarn
      integer  :: icheckn, CHECKIN
c
c  List of valid interrupt messages
c
      character*3 :: intmsg(17)
      data intmsg /  'sto','?  ','pau','abo','tli','cpu','nli','dum'
     &  ,'dtd','hdf','dtf','tsl','dtt','hst','dth','usr','dtu' /
c
cRAF: Number of warnings to issue before giving up
c
      data maxwarn /20/
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c  Check stopping criteria
c
c We check to see if the CPU time limit is approaching only in
c interactive mode (mbatch = 0), because this requires
c broadcasting the master's value of the CPU time used.
c
c Do we have to make the following system call every time step?
c
      if (myid_w .eq. 0) then
        call clocks (cputime, wclock)
#ifndef ARCH_CRAY
        tused = real(cputime)
#else
        tused = wclock
#endif
      endif
c
      if (tlim .gt. 0.0 .and. time .ge. tlim) then
        if (myid_w .eq. 0) 
     &    write(6,"(/1x,'terminating on physical time limit',
     &    /1x,'tlim=',1pe12.5,'   nlim=',i7,'  cpulim=',1pe12.5,
     &    /1x,'time=',1pe12.5,'  cycle=',i7,'   tused=',1pe12.5)") 
     &    tlim,nlim,cpulim,time,nhy,tused
        ifsen(1) = 1
      endif
      if (nlim .gt. 0    .and. nhy .ge. nlim) then
        if (myid_w .eq. 0) 
     &    write(6,"(/1x,'terminating on cycle limit',
     &    /1x,'tlim=',1pe12.5,'   nlim=',i7,'  cpulim=',1pe12.5,
     &    /1x,'time=',1pe12.5,'  cycle=',i7,'   tused=',1pe12.5)") 
     &    tlim,nlim,cpulim,time,nhy,tused
        ifsen(1) = 1
      endif
c
cRAF: Stop if more than maxwarn warnings have been issued.
c
      if (nwarn .gt. maxwarn) then
        if (myid_w .eq. 0) 
     &    write(6,"(/1x,'terminating after ',i4,' warnings',
     &    /1x,'tlim=',1pe12.5,'   nlim=',i7,'  cpulim=',1pe12.5,
     &    /1x,'time=',1pe12.5,'  cycle=',i7,'   tused=',1pe12.5)") 
     &    maxwarn,tlim,nlim,cpulim,time,nhy,tused
        ifsen(1) = 1
      endif
c
      if (ifsen(1) .eq. 1) return
      if (mbatch .ne. 0) return
#ifndef ARCH_CRAY
c
c.......................................................................
c
c FOR INTERACTIVE MODE ONLY -- REQUIRES COMMUNICATION
c
c Stop if the CPU time is approaching the limit.  See if there's 
c enough time remaining to take an average time step.
c
#ifdef MPI_USED
      call MPI_BCAST(tused, 1, MPI_FLOAT, 0, comm3d, ierr)
c
      if (cpulim - tsave - tused .lt. tused/real(nhy)) then
        if (myid_w .eq. 0) 
     &    write(6,"(/1x,'terminating on single-process CPU time limit',
     &    /1x,'tlim=',1pe12.5,'   nlim=',i7,'  cpulim=',1pe12.5,
     &    /1x,'time=',1pe12.5,'  cycle=',i7,'   tused=',1pe12.5)") 
     &    tlim,nlim,cpulim,time,nhy,tused
        ifsen(1) = 1
        return
      endif
#endif
c
c  Check for interrupt messages if not in batch mode.  If none or 
c  illegal message found then return
c
      icheckn = 1
      if (myid .eq. 0) nchar = CHECKIN(msg,icheckn)
c
#ifdef MPI_USED
       call MPI_BCAST (nchar, 1, MPI_INTEGER, 0, comm3d, ierr)
#endif
c
      if (nchar .eq. 0) return
c
#ifdef MPI_USED
       call MPI_BCAST (msg, 3, MPI_CHARACTER, 0, comm3d, ierr)
#endif
c
      do i=1,15
        if (msg(1:3) .eq. intmsg(i)) goto 20
      enddo
      if (myid .eq. 0)
     &write(6,"(1x,a3,' is not an interrupt message.  Legal messages ',
     &'are:',/1x,'sto ? pau abo tli nli cpu dum dtd hdf dtf hst dth',
     &' usr dtu')") msg(1:3)
      return
20    continue
c
c  Legal interrupt message found, process it
c
c  stop command
c
      if (msg(1:3) .eq. 'sto') then
        if (myid .eq. 0)
     &  write(6,"(1x,a3,': execution stopped with',
     &  /1x,'tlim=',1pe12.5,'   nlim=',i7,'  cpulim=',1pe12.5,
     &  /1x,'time=',1pe12.5,'  cycle=',i7,'   tused=',1pe12.5)") 
     &  msg,tlim,nlim,cpulim,time,nhy,tused
        ifsen(1) = 1
        return
      endif
c
c  status command
c
      if (msg(1:3) .eq. '?  ') then
        if (myid .eq. 0) 
     &  write(6,"(1x,a3,': execution continuing with',
     &  /1x,'tlim=',1pe12.5,'   nlim=',i7,'  cpulim=',1pe12.5,
     &  /1x,'time=',1pe12.5,'  cycle=',i7,'   tused=',1pe12.5,
     .  /1x,'dt  =',1pe12.5)") 
     &  msg,tlim,nlim,cpulim,time,nhy,tused,dt
        return
      endif
c
c  pause command
c
      if (msg(1:3) .eq. 'pau') then
        if (myid .eq. 0)
     &  write(6,"(1x,a3,': execution halted with',
     &  /1x,'tlim=',1pe12.5,'   nlim=',i7,'  cpulim=',1pe12.5,
     &  /1x,'time=',1pe12.5,'  cycle=',i7,'   tused=',1pe12.5,
     &  /1x,'Hit any key to restart execution')") 
     &  msg,tlim,nlim,cpulim,time,nhy,tused
      icheckn = 0
        if (myid .eq. 0) nchar = CHECKIN(msg,icheckn)
c
#ifdef MPI_USED
        call MPI_BARRIER (comm3d, ierr)
#endif
c
        return
      endif
c
c  abort command
c
      if (msg(1:3) .eq. 'abo') then
        if (myid .eq. 0)
     &  write(6,"(a3,': ABORT! do you really want to abort execution?',
     &   ' (type yes or no)')") msg
      icheckn = 0
        if (myid .eq. 0) nchar = CHECKIN(msg,icheckn)
c
#ifdef MPI_USED
        call MPI_BCAST (nchar, 1, MPI_INTEGER, 0, comm3d, ierr)
#endif
c
        if (nchar .eq. 0) return
c
#ifdef MPI_USED
        call MPI_BCAST (msg, 3, MPI_CHARACTER, 0, comm3d, ierr)
#endif
c
        if (msg(1:3) .ne. 'yes') then
          if (myid .eq. 0)
     &    write(6,"('Abort cancelled, continuing execution ...')")
          return
        else
          if (myid .eq. 0)
     &    write(6,"('ABORT.................')")
c
#ifdef MPI_USED
          call MPI_FINALIZE (ierr)
#endif
c
          stop
        endif
      endif
c
c  reset physical time limit (tlim) command
c
      if (msg(1:3) .eq. 'tli') then
        call findno(msg,istrt,iend)
        if (istrt .lt. 0) then
          if (myid .eq. 0) write(6,2000) msg
          return
        endif
        call BCDFLT(msg,valnew,(istrt-1),(iend-istrt+1))
        if (valnew .lt. 0.0 .or. valnew .ge. huge) then
          if (myid .eq. 0) write(6,2000) msg
2000      format(1x,a3,': could not read reset number; execution ',
     &    'continuing')
          return
        endif
        tlim = valnew
        if (myid .eq. 0)
     &  write(6,"(a3,': tlim reset to ',1pe12.5)") msg,tlim
        return
      endif
c
c  reset CPU time limit (cpulim) command
c
      if (msg(1:3) .eq. 'cpu') then
        call findno(msg,istrt,iend)
        if (istrt .lt. 0) then
          if (myid .eq. 0) write(6,2000) msg
          return
        endif
        call BCDFLT(msg,valnew,(istrt-1),(iend-istrt+1))
        if (valnew .lt. 0.0 .or. valnew .ge. huge) then
          if (myid .eq. 0) write(6,2000) msg
          return
        endif
        cpulim = valnew
        if (myid .eq. 0)
     &  write(6,"(a3,': cpulim reset to ',1pe12.5)") msg,cpulim
        return
      endif
c
c  reset cycle limit (nlim) command
c
      if (msg(1:3) .eq. 'nli') then
        call findno(msg,istrt,iend)
        if (istrt .lt. 0) then
          if (myid .eq. 0) write(6,2000) msg
          return
        endif
        call BCDFLT(msg,valnew,(istrt-1),(iend-istrt+1))
        if (valnew .lt. 0.0 .or. valnew .ge. huge) then
          if (myid .eq. 0) write(6,2000) msg
          return
        endif
        nlim = nint(valnew)
        if (myid .eq. 0) write(6,"(a3,': nlim reset to ',i12)") msg,nlim
      endif
c
c  turn restart dump switch on
c
      if (msg(1:3) .eq. 'dum') then
        if (myid .eq. 0) write(6,"(a3,': restart dump switch on')") msg
        iswres = 1
        return
      endif
c
c  reset dump frequency (dtdump) command
c
      if (msg(1:3) .eq. 'dtd') then
        call findno(msg,istrt,iend)
        if (istrt .lt. 0) then
          if (myid .eq. 0) write(6,2000) msg
          return
        endif
        call BCDFLT(msg,valnew,(istrt-1),(iend-istrt+1))
        if (valnew .lt. 0.0 .or. valnew .ge. huge) then
          if (myid .eq. 0) write(6,2000) msg
          return
        endif
        dtdump = valnew
        if (myid .eq. 0)
     &  write(6,"(a3,': dtdump reset to ',1pe12.5)") msg,dtdump
      endif
c
c  turn hdf dumps on
c
      if (msg(1:3) .eq. 'hdf') then
        if (myid .eq. 0) write(6,"(a3,': hdf dump switch on')") msg
        iswhdf = 1
        return
      endif
c
c  reset hdf dump frequency (dthdf) command
c
      if (msg(1:3) .eq. 'dtf') then
        call findno(msg,istrt,iend)
        if (istrt .lt. 0) then
          if (myid .eq. 0) write(6,2000) msg
          return
        endif
        call BCDFLT(msg,valnew,(istrt-1),(iend-istrt+1))
        if (valnew .lt. 0.0 .or. valnew .ge. huge) then
          if (myid .eq. 0) write(6,2000) msg
          return
        endif
        dthdf = valnew
        if (myid .eq. 0)
     &  write(6,"(a3,': dthdf reset to ',1pe12.5)") msg,dthdf
      endif
c
c  turn tsl dumps on
c
      if (msg(1:3) .eq. 'tsl') then
        if (myid .eq. 0) write(6,"(a3,': tsl dump switch on')") msg
        iswtsl = 1
        return
      endif
c
c  reset tsl dump frequency (dttsl) command
c
      if (msg(1:3) .eq. 'dtt') then
        call findno(msg,istrt,iend)
        if (istrt .lt. 0) then
          if (myid .eq. 0) write(6,2000) msg
          return
        endif
        call BCDFLT(msg,valnew,(istrt-1),(iend-istrt+1))
        if (valnew .lt. 0.0 .or. valnew .ge. huge) then
          if (myid .eq. 0) write(6,2000) msg
          return
        endif
        dttsl = valnew
        if (myid .eq. 0)
     &  write(6,"(a3,': dttsl reset to ',1pe12.5)") msg,dttsl
      endif
c
c  turn history dumps on
c
      if (msg(1:3) .eq. 'hst') then
        if (myid .eq. 0) write(6,"(a3,': history dump switch on')") msg
        iswhst = 1
        return
      endif
c
c  reset history dump frequency (dth) command
c
      if (msg(1:3) .eq. 'dth') then
        call findno(msg,istrt,iend)
        if (istrt .lt. 0) then  
          if (myid .eq. 0) write(6,2000) msg  
          return 
        endif  
        call BCDFLT(msg,valnew,(istrt-1),(iend-istrt+1))
        if (valnew .lt. 0.0 .or. valnew .ge. huge) then
          if (myid .eq. 0) write(6,2000) msg
          return 
        endif  
        dthist = valnew   
        if (myid .eq. 0)
     &  write(6,"(a3,': dthist reset to ',1pe12.5)") msg,dthist
      endif
c
c  turn USER dumps on
c
      if (msg(1:3) .eq. 'usr') then
        if (myid .eq. 0) write(6,"(a3,': USER dump switch on')") msg
        iswusr = 1
        return
      endif
c
c  reset USER dump frequency (dtu) command
c
      if (msg(1:3) .eq. 'dtu') then
        call findno(msg,istrt,iend)
        if (istrt .lt. 0) then
          if (myid .eq. 0) write(6,2000) msg
          return
        endif
        call BCDFLT(msg,valnew,(istrt-1),(iend-istrt+1))
        if (valnew .lt. 0.0 .or. valnew .ge. huge) then
          if (myid .eq. 0) write(6,2000) msg
          return
        endif
        dtusr = valnew   
        if (myid .eq. 0)
     &  write(6,"(a3,': dtusr reset to ',1pe12.5)") msg,dtusr
      endif
c
#endif /* NOT ARCH_CRAY */
      return
      end
