c-----------------------------------------------------------------------
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
       program zmp_pp
c
c ZEUS-MP Post-processor: Main Program.
c
c PURPOSE:  This program asks the user what post-processing taskes
c           to perform, and calls the appropriate routines.
c
c LAST MODIFIED: RAF, 9/12/96.
c
c.......................................................................
c
      use real_prec
      use param
c
      implicit NONE
c
      integer :: i, ltask
      character :: task*12
c.......................................................................
c
c Write greeting.
c
       write(*,"(/'ZMP_PP: Welcome to the ZEUS-MP Post-Processor!')")
       ltask = len( task )
c
c Loop back here for next task.
c
   10  continue
c
c Determine task.
c
       do 20 i=1,ltask
         task(i:i) = ' '
   20  continue
       write(*,"(/'What task would you like performed?')")
       write(*,"( 'Type one of the following:'/)")
       write(*,"( '   h4splice   (concatenate HDF4 files one by one)')")
       write(*,"( '   h5splice   (concatenate HDF5 files one by one)')")
       write(*,"( '   auto_h4    (automated h4splice)')")
       write(*,"( '   auto_h5    (automated h5splice)')")
       write(*,"( '   quit       (exit ZMP_PP)'/)")
       read(*,"(a)") task
c
c QUIT
c
       if (task(1:4) .eq. 'quit') then
         stop
       else if (task(1:7) .eq. 'h4splice') then
c
c SPLICE HDFDUMPs -- HDF Scientific Data Sets "hdf<id><coords>.<ndump>"
c
         call h4splice
         goto 10
       else if (task(1:8) .eq. 'h5splice') then
c
c SPLICE HDF5 DUMPs -- HDF5 Scientific Data Sets "hdf<id><coords>.<ndump>"
c
         call h5splice
         goto 10
       else if (task(1:12) .eq. 'auto_h4') then
c
c AUTO-SPLICE HDF4 DUMPs -- HDF Scientific Data Sets "hdf<id><coords>.<ndump>"
c
         call auto_h4
         goto 10
       else if (task(1:12) .eq. 'auto_h5') then
c
c AUTO-SPLICE HDF5 DUMPs -- HDF5 Scientific Data Sets "hdf<id><coords>.<ndump>"
c
         call auto_h5
         goto 10
       else
c
c WHAT?
c
         write(*,"(/'ZMP_PP: Please retype your task.')")
         goto 10
       endif
c
       end
