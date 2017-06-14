c=======================================================================
c
       subroutine clocks (cputime, wclock)
c
c This routine obtains the CPU and wall-clock times in seconds for the 
c calling process since the times (cputime0, wclock0) at which the 
c clocks were initalized (passed through common block /clocks/). 
c
c With MPI, MPI_WALL time is used except on the EXEMPLAR.
c
c For systems without a CPU timer, cputime is set to zero.
c
c Written by RAF, last modified 3/25/96.
c......................................................................
c
      use real_prec
      use param
      use config
      use clockmod
c
      implicit NONE
c
#ifndef ARCH_CRAY
      real(rl4) :: cputime, wclock
#else
      real(rl ) :: cputime, wclock
#endif
c
#ifdef ARCH_IBM
#define ETIME   etime_
#else
#define ETIME   etime
#endif
c
#ifndef ARCH_TERAGRID
      real(rl) :: ETIME
#endif
c
c
#ifdef MPI_USED
      real(rl) :: MPI_WTIME, wall
#endif
c
#ifdef MPI_USED
c
c Wall clock time is easy to get with MPI:
c
      wall   = MPI_WTIME()
#ifndef ARCH_CRAY
      wclock = real(wall) - wclock0
#else
      wclock = wall - wclock0
#endif
#endif /* MPI_USED */
c
c Get the CPU time for this process/thread.
c
c For systems with the standard UNIX etime/itime routines.  Note that
c the output from itime is an integer array with (hours,minutes,seconds)
c and therefor is accurate only to the nearest second.  The output from
c function etime itself is the sum of user plus system times, which can
c be significantly longer than the user time alone (stored in 
c tarray(1)). 
c
#ifndef MPI_USED
      cputime = 0.0
#endif
c
#ifndef ARCH_TERAGRID
#ifndef ARCH_CRAY
      cputime = ETIME ( tarray ) - cputime0
#else
      cputime = wclock
#endif
#endif
c
      end
