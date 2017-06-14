c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              M A R S H A K                \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c     PURPOSE: initializes Marshak wave test problem.
c
c     Written by: Robert Fiedler and John Hayes
c
c=======================================================================
#define ABSORPTION        0.57735
      subroutine marshak
c
      use real_prec
      use config
      use param
      use cons
      use grid
      use field
      use radiation
      use opac
      use root
      use bndry
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer  :: i, j, k
      real(rl) :: d1, t0
      real(rl) :: third
c
      real(rl) :: d0     , e0     , er0    , soeps
      real(rl) :: ros_mfp, flx_lim, dmc_max, dx_min
c
      real(rl) :: so_eps , tau    , eriibn, timarg
c
      common /soui/ so_eps, tau, eriibn
c
      REAL souis
c
      namelist /pgen/ d0, e0, er0, soeps
c
c=======================================================================
c
c
c     initialize and read in parameters from PGEN namelist
c
       d0    =  1.0      !  mass density in gm/cc
       e0    =  0.0      !  gas energy
       er0   =  0.0      !  rad energy
       soeps =  0.1      !  epsilon defined in Su and Olson eq. (13)
c
      if (myid_w .eq. 0) then
        read (1,pgen)
        write(2,pgen)
#ifdef MPI_USED
       buf_in(1) = d0
       buf_in(2) = e0
       buf_in(3) = er0
       buf_in(4) = soeps
      endif
       call MPI_BCAST( buf_in, 4, MPI_FLOAT
     &               , 0, MPI_COMM_WORLD, ierr )
      if(myid_w .ne. 0) then
       d0    = buf_in(1)
       e0    = buf_in(2)
       er0   = buf_in(3)
       soeps = buf_in(4)
#endif /* MPI_USED */
      endif
      so_eps = soeps
c
c Compute the time-dependent BC (eriib) at time = 0.0.
c
       timarg = 0.0D0
       eriibn = souis (timarg)
c
c Copy the BC into the boundary value array.
c
       do k=ks-1,ke+1
         do j=js-1,je+1
           eriib(j,k,1) = eriibn
         enddo ! j
       enddo ! k
c
c     initialize field arrays
c
      do 1 k = 1, kn
      do 1 j = 1, jn
      do 1 i = 1, in
       d (i,j,k) = d0
       v1(i,j,k) = 0.0
       v2(i,j,k) = 0.0
       v3(i,j,k) = 0.0
       e (i,j,k) = e0
       er(i,j,k) = er0
1     continue
c
      return
      end
c
      subroutine source
c
      use real_prec
      use config
      use param
      use root
      use bndry
      use grid
c
      implicit NONE
c
      integer  :: i,j,k
c
      real(rl) :: souis, eriibn
c
      eriibn = souis(time)
      do k=ks-1,ke+1
        do j=js-1,je+1
          eriib(j,k,1) = eriibn
          eriib(j,k,2) = eriibn
        enddo ! j
      enddo ! k
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   F U N C T I O N          //////////
c    //////////                 S O U I S                 \\\\\\\\\\
c
c=======================================================================
c
       real*8 function souis (time)
c
c PURPOSE
c   For non-equilibrium Marshak wave test problem (Su and Olson).
c   Returns the value of the radiation energy density (called u(x,tau))
c   at the surface x=0 at the current problem time for a given value
c   of the parameter so_eps (passed through common block soui).
c
c   Assumes that 4 * F_inc / c has been set to unity so that u(0,tau)
c   equals E(0,t).
c
c   This output is used to specify the time-dependent BC.
c
c AUTHOR
c   Robert A. Fiedler, 1/16/97
c
c LAST MODIFIED
c   1/16/97
c
c USES
c   qromo (Numerical Recipes), so_u_i1, so_u_i2 -- the integrands.
c
c TEST PROGRAM
c       program sou
cc
cc Test main program for routine souis.
cc
c       implicit NONE
c#include "param.h"
c#include "root.h"
c#include "cons.h"
c       integer i
c       REAL kap, so_eps, tau, eriibn
c       common /soui/ so_eps, tau, eriibn
c       REAL souis
c       external souis
c       REAL taus(11)
c       data taus /0.001, 0.003, 0.01, 0.03, 0.1, 0.3
c     1          ,1.0  , 3.0  , 10.0, 30.0, 100./
c       write(*,"('Enter epsilon')")
c       read (*,*) so_eps
c       kap = ABSORPTION
c       do i=1,11
c         time = taus(i) / epsilon / clight / kap
c         eriibn = souis(time)
c         write(*,"(' u(',1pe10.3,') = ',e15.5)") taus(i), eriibn
c       enddo
c       end
c.......................................................................
c
      use real_prec
      use config
      use param
c      use root
      use grid
      use cons
c
c Do not include bndry.h when using with test program sou.
c
      use bndry
c
       implicit NONE
c
       real(rl) :: so_eps, tau, eriibn
       common /soui/ so_eps, tau, eriibn
       real(rl) :: time
       real(rl) :: kap, qa, qb
       integer  :: j,k
c
       real(rl) :: so_u_i1, so_u_i2
c
       external so_u_i1, so_u_i2
c
       kap = ABSORPTION
       tau = so_eps * clight * kap * time
c
c Numerically evaluate the integrals of equation (36) in Su and Olson.
c They are functions of so_eps and tau(time).
c
       call qromo ( so_u_i1, zro, one, qa )
       call qromo ( so_u_i2, zro, one, qb )
c
c Find u(0,tau).
c
       souis = one - (sqrt(3.0)/pi) * ( two * qa + exp(-tau) * qb )
c
       return
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                S O U I S                  \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   F U N C T I O N          //////////
c    //////////               S O U I 1                   \\\\\\\\\\
c
c=======================================================================
c
       real*8 function so_u_i1 ( eta )
c
c First integrand in Su and Olson, eq. (36).
c
c Written by Robert Fiedler, 1/16/97.
c
      use real_prec
      use param
c
      implicit NONE
c
      real(rl) :: eta, g1, th
c
      real(rl) :: so_eps, tau, eriibn
      common /soui/ so_eps, tau, eriibn
c
      g1 = eta * sqrt ( so_eps + one / ( max ( one - eta**2, tiny ) ) )
      th = acos ( sqrt ( 3.0 / ( 3.0 + 4.0 * g1**2 ) ) )
      so_u_i1 = exp ( -tau * eta**2 ) * sin ( th )
     1         / ( eta * sqrt ( 3.0 + 4.0 * g1**2 ) )
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   F U N C T I O N            //////////
c    //////////                S O U I 1                  \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   F U N C T I O N          //////////
c    //////////               S O U I 2                   \\\\\\\\\\
c
c=======================================================================
c
       real*8 function so_u_i2 ( eta )
c
c Second integrand in Su and Olson, eq. (36).
c
c Written by Robert Fiedler, 1/16/97.
c
      use real_prec
      use param
      implicit NONE
c
      real(rl) :: eta, g2, th
c
      real(rl) :: so_eps, tau, eriibn
      common /soui/ so_eps, tau, eriibn
c
       g2 = sqrt ( max ( one - eta, tiny ) * ( so_eps + one 
     1    / max ( eta, tiny ) ) )
       th = acos ( sqrt ( 3.0 / ( 3.0 + 4.0 * g2**2 ) ) )
       so_u_i2 = exp ( -tau / ( so_eps * max ( eta, tiny ) ) ) 
     1         * sin ( th )
     2         / ( max ( eta, tiny ) * ( one + so_eps * eta ) 
     3         * sqrt ( 3.0 + 4.0 * g2**2 ) )
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   F U N C T I O N            //////////
c    //////////                S O U I 2                  \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 Q R O M O                 \\\\\\\\\\
c
c=======================================================================
c
       subroutine qromo ( func, x1, x2, ss )
c
c    dac:zeus3d.qromb <------------------- integrates func from x1 to x2
c    from whp:numerical recipes                           december, 1992
c
c    written by: David Clarke
c    modified 1: Robert Fiedler 1/7/97, for improper integrals; 
c                see Numerical Recipes in Fortran, 2nd ed., p. 137.
c    modified 2: Robert Fiedler 1/16/97, quit for small errors dss;
c                tuned for Su and Olson non-equilibrium Marshak wave.
c
c  PURPOSE:  Returns the definite integral of the function "func"
c  between specified limits "x1" and "x2" using Romberg's method of
c  order 2k, where k=2 is Simpson's rule (see Numerical Recipes, 1st
c  edition for FORTRAN, page 114).
c
c  INPUT VARIABLES:
c    func     name of external function describing the univariate
c             function to be integrated.
c    x1, x2   integration limits
c
c  OUTPUT VARIABLES:
c    ss       value of definite integral
c
c  EXTERNALS:
c    FUNC
c    MIDPNT
c    POLINT
c
c-----------------------------------------------------------------------
      use real_prec
      use param
c
      implicit NONE
c
      integer, parameter :: isig=8
      integer, parameter :: jmx=20
      integer, parameter :: k=10
      integer, parameter :: km=k-1
c
      real(rl),parameter :: eps=0.1**isig
c
       integer  :: j
       real(rl) :: x1      , x2      , ss      , dss
c
       real(rl) :: s (jmx+1), h (jmx+1)
c
c      External statements
c
       real(rl) :: func
       external      func    , midpnt  , polint
c
c-----------------------------------------------------------------------
c
       h(1) = 1.0
       do 10 j=1,jmx
         call midpnt ( func, x1, x2, s(j), j )
         if (j .ge. k) then
           call polint ( h(j-km), s(j-km), k, zro, ss, dss )
c           if (abs(dss).lt.eps*abs(ss)) return
c
c Change from Numerical Recipes -- give up if the error is very small.
c
           if (abs(dss).lt.eps*abs(ss) .or. abs(dss).lt.tiny) return
         endif
         s(j+1) = s(j)
         h(j+1) = h(j) / 9.0  !  Step tripling and even error series.
10     continue
       write (6, 2000) jmx, dss, ss
       return
c
c-----------------------------------------------------------------------
c----------------------- Write format statements -----------------------
c-----------------------------------------------------------------------
c
2000   format('QROMO   : *** WARNING *** Romberg Integration failed '
     1       ,'to converge in ',i2,' steps.',/
     1       ,'QROMO   : error =',1pg12.5,', definite integral =',g12.5)
c
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                 Q R O M O                 \\\\\\\\\\
c
c=======================================================================
c
c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                M I D P N T                \\\\\\\\\\
c
c=======================================================================
c
       subroutine midpnt ( func, x1, x2, s, n )
c
c    dac:zeus3d.trapzd <-------------- n'th refinement of trapezoid rule
c    from whp:numerical recipes                           december, 1992
c
c    written by: David Clarke
c    modified 1: Robert Fiedler, 1/8/97, changed to midpnt for
c                improper integrals; see Numerical Recipes in Fortran,
c                2nd. ed., p. 136.
c
c  PURPOSE:  This routine returns in "s" the n'th stage of refinement of
c  an extended trapezoidal rule.  "func" is input as the name of the
c  function to be integrated between limits "x1" and "x2", also input.
c  When called with n=1, the crudest estimate of the definite integral
c  is returned.  Subsequent calls with n=2,3,... (in that order) will
c  improve the accuracy of "s" by adding 2**(n-2) additional interior
c  points.  The value of "s" should not be modified between successive
c  calls (see Numerical Recipes, 1st edition for FORTRAN, page 111).
c
c  INPUT VARIABLES:
c    func     name of external function describing the univariate
c             function to be integrated.
c    x1, x2   integration limits
c    n        stage of refinement
c
c  OUTPUT VARIABLES:
c    s        value of definite integral after n'th stage of refinement.
c
c  EXTERNALS: [NONE]
c
c-----------------------------------------------------------------------
c
      use real_prec
      use param
c
      implicit NONE
c
      integer  :: n       , it      , j
      real(rl) :: x1      , x2      , s      , tnm     , del
     1             , x       , sum   , ddel
c
c      External statements
c
       real(rl) ::   func
       external      func
c
c-----------------------------------------------------------------------
c
       if (n .eq. 1) then
         s  = ( x2 - x1 ) * func ( haf * ( x1 + x2 ) )
       else
         it   = 3**( n - 2 )
         tnm  = real ( it )
         del  = ( x2 - x1 ) / ( 3.0 * tnm )
         ddel = del + del  !  Added pts alternate in spacing: del, ddel.
         x    = x1 + 0.5 * del
         sum  = 0.0
         do 10 j=1,it
           sum = sum + func ( x )
           x   = x + ddel
           sum = sum + func ( x )
           x   = x +  del
10       continue
         s  = ( s + ( x2 - x1 ) * sum / tnm ) / 3.0
       endif
c
       return
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                M I D P N T                \\\\\\\\\\
c
c=======================================================================
c
c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                P O L I N T                \\\\\\\\\\
c
c=======================================================================
c
       subroutine polint ( xa, ya, n, x, y, dy )
c
c    dac:zeus3d.polint <---------------------- interpolates gridded data
c    from whp:numerical recipes                           december, 1992
c
c    written by: David Clarke
c    modified 1: Robert Fiedler, 1/7/97, for ZEUS-MP.
c
c  PURPOSE:  Given arrays "xa" and "ya", each of length "n", and a given
c  value "x", this routine returns a value "y" and an error estimate
c  "dy" (see Numerical Recipes, 1st edition for FORTRAN, page 82).
c
c  INPUT VARIABLES:
c    xa       input independent variable array
c    ya       input dependent   variable array
c    n        length of "xa" and "ya"
c    x        value of independent variable at which dependent variable
c             is to be interpolated
c
c  OUTPUT VARIABLES:
c    y        interpolated value of dependent variable
c    dy       estimate of error of interpolation.
c
c  EXTERNALS: [NONE]
c
c-----------------------------------------------------------------------
c
      use real_prec
      use param
c
      implicit NONE
c
      integer, parameter :: nmax=20
c
      integer  :: n       , ns      , i       , m
      real(rl) :: x       , y       , dy      , dif     , dift
     1             , ho      , hp      , w       , den
c
      real(rl) :: xa      (   n), ya      (   n)
     1             , c       (nmax), d       (nmax)
c
c-----------------------------------------------------------------------
c
       ns  = 1
       dif = abs ( x - xa(1) )
       do 10 i=1,n
         dift = abs ( x - xa(i) )
         if (dift .lt. dif) then
           ns  = i
           dif = dift
         endif
         c(i) = ya(i)
         d(i) = ya(i)
10     continue
       y  = ya(ns)
       ns = ns - 1
       do 30 m=1,n-1
         do 20 i=1,n-m
           ho   = xa(i  ) - x
           hp   = xa(i+m) - x
           w    = c (i+1) - d(i)
           den  = w / ( ( ho - hp ) + tiny )
           c(i) = ho * den
           d(i) = hp * den
20       continue
         if (2*ns .lt. n-m) then
           dy = c(ns+1)
         else
           dy = d(ns  )
           ns = ns - 1
         endif
         y = y + dy
30     continue
c
       return
       end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                P O L I N T                \\\\\\\\\\
c
c=======================================================================
c
