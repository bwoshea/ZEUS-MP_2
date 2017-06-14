c=======================================================================
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
      subroutine dataio(iswres,iswhdf,iswtsl,iswhst,iswusr)
c
c  PURPOSE:  Controls data I/O for restart, HDF, history, and ASCII
c  dumps.
c
c  INPUT ARGUMENTS: iswres,iswhdf,iswhst=switches for restart,hdf, and
c    history dumps.  Values of 1 ensure dumps will be made.
c
c  OUTPUT ARGUMENTS: none
c
c
c  LAST MODIFIED:  by John Hayes; 9/9/97
c  LAST MODIFIED:  by efh 04/15/99 including call to tslice
c  LAST MODIFIED:  by John Hayes, 05/26/2005; removed call to tslice
c-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use root
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c-----------------------------------------------------------------------
      implicit NONE
c
      integer  :: iswres,iswhdf,iswhst,iswusr,iswshl,iswtsl
      integer  :: incr,strtoi
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
cFH   if (dtdump .gt. 0.0 .and. tused .ge. (tdump+dtdump)) then
      if (dtdump .gt. 0.0 .and. time .ge. (tdump+dtdump)) then
        tdump = tdump + dtdump
        iswres = 1
      endif
      if (dthdf  .gt. 0.0 .and. time  .ge. (thdf+dthdf)) then
        thdf = thdf + dthdf
        iswhdf = 1
      endif
      if (dttsl  .gt. 0.0 .and. time  .ge. (ttsl+dttsl)) then
        ttsl = ttsl + dttsl
        iswtsl = 1
      endif
      if (dthist .gt. 0.0 .and. time  .ge. (thist+dthist)) then
        thist = thist + dthist
        iswhst = 1
      endif
      if (dtusr .gt. 0.0 .and. time  .ge. (tusr+dtusr)) then
        tusr = tusr + dtusr
        iswusr = 1
      endif
c
c  restart dump
c
      if(xrestart) then
      if (iswres .eq. 1) then
        call msave(resfile)
        incr = strtoi(resfile,13,15) + 1
        write(resfile,"(a3,a2,3i2.2,'.',i3.3)") 'res',id,coords(1)
     &                                      ,coords(2),coords(3),incr
        iswres=0
      endif
      endif ! xrestart
c
c  HDF dump
c
      if(xhdf) then
       if (iswhdf .eq. 1) then
        call hdfall(hdffile)
        incr = strtoi(hdffile,13,15) + 1
        write(hdffile,"(a3,a2,3i2.2,'.',i3.3)") 'hdf',id,coords(1)
     &                                      ,coords(2),coords(3),incr
c
c If an array of output times t_out has been defined, set the value
c of (thdf + dthdf) to the next desired output time.
c
        if (incr .le. nbuff-8) then
          if (t_out(incr).ne.0.0) then
            thdf = t_out(incr) - dthdf
          endif
        endif
        iswhdf=0
       endif
      endif ! xhdf
c
c  ascii dump
c
      if(xascii) then
      if (iswusr .eq. 1) then
        call textdmp
        incr = strtoi(usrfile,13,15) + 1
        write(usrfile,"(a3,a2,3i2.2,'.',i3.3)") 'usr',id,coords(1)
     &                                      ,coords(2),coords(3),incr
c
c If an array of output times t_out has been defined, set the value
c of (tusr + dtusr) to the next desired output time.
c
        if (incr .le. nbuff-8) then
          if (t_out(incr).ne.0.0) then
            tusr = t_out(incr) - dtusr
          endif
        endif
        iswusr=0
      endif
      endif ! xascii
c
      return
      end
