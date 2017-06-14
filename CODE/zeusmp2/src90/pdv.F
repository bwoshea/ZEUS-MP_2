c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                   P D V                   \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine pdv 
     &                (ibeg,iend,jbeg,jend,kbeg,kend,dlo,eod)
c
c Updates the energy with the PdV compressional work term.  Also
c computes the specific internal energy density eod for the transport
c step.
c
c BOUNDARY VALUES USED:
c
c  Macro defined  var   ii    oi    ij    oj    ik    ok
c  -------------  ---  ----  ----  ----  ----  ----  ----
c                  u1        ie+1
c                  u2                    je+1
c                  u3                                ke+1
c
c OUTPUT
c    dlo    Mass density (copy)                      for transport step.
c    eod    Specific energy density e/d [or (e+p)/d] for transport step.
c
c Written by RAF 2/12/96.
c
c Last modified: 09/01/2006 by John Hayes; added "implicit NONE" (duh!)
c                and added declaration for div_v
c
c......................................................................
      use real_prec
      use config
      use param
      use root
      use grid
      use field
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend, kp1,
     .            jp1
c
      real(rl) :: dlo(in,jn,kn), eod(in,jn,kn)
      real(rl) :: q1, q2, div_v
c
      if(xtotnrg .eqv. .false.) then
       q1 = 0.5 * dt * gamm1
      endif ! xtotnrg
c
      do 30 k=kbeg,kend
       if(ldimen .eq. 3) then
        kp1 = k+1
       else
        kp1 = ks
       endif
       do 20 j=jbeg,jend
        if(ldimen .ge. 2) then
         jp1 = j+1
        else
         jp1 = js
        endif
        do 10 i=ibeg,iend
         if(xtotnrg .eqv. .false.) then
c
c----------------------------------------------------------------------
c      Compute divergence of velocity field and add pdv source term 
c      to energy density.
c----------------------------------------------------------------------
c
          if(ldimen .eq. 3) then           
           q2         = ( g2a(i+1) * g31a(i+1) * v1(i+1,j,k)
     1                - g2a(i  ) * g31a(i  ) * v1(i  ,j,k) )
     2                  *                         dvl1ai(i)
     3                  + ( g32a(j+1) * v2(i,j+1,k)
     4                    - g32a(j  ) * v2(i,j  ,k) )
     5                  *   g2bi(i)             * dvl2ai(j)
     6                  + ( v3(i,j,kp1) - v3(i,j,k) )
     7                  *   g31bi(i) * g32bi(j) * dvl3ai(k)
          endif 
          if(ldimen .eq. 2) then
           q2         = ( g2a(i+1) * g31a(i+1) * v1(i+1,j,k)
     1                    - g2a(i  ) * g31a(i  ) * v1(i  ,j,k) )
     2                  *                         dvl1ai(i)
     3                  + ( g32a(j+1) * v2(i,j+1,k)
     4                    - g32a(j  ) * v2(i,j  ,k) )
     5                  *   g2bi(i)             * dvl2ai(j)
          endif ! ldimen
          if(ldimen .eq. 1) then
             q2         = ( g2a(i+1) * g31a(i+1) * v1(i+1,j,k)
     1                    - g2a(i  ) * g31a(i  ) * v1(i  ,j,k) )
     2                  *                         dvl1ai(i)
          endif ! ldimen
          div_v = q2
c
c----------------------------------------------------------------------
c     Case of simple gamma-law EOS
c----------------------------------------------------------------------
c
          if(leos .eq. 1) then
           q2         = q2 * q1
           e(i,j,k)   = ( 1.0 - q2 ) / ( 1.0 + q2 ) * e(i,j,k)
          endif ! ideal gas EOS
          eod(i,j,k) = e(i,j,k) / d(i,j,k)
         else ! xtotnrg
          eod(i,j,k) = gamma * e(i,j,k) / d(i,j,k)
     1               - gamm1 * ( ( v1(i,j,k) + v1(i+1,j  ,k  ) )**2
     2                         + ( v2(i,j,k) + v2(i  ,jp1,k  ) )**2
     3                         + ( v3(i,j,k) + v3(i  ,j  ,kp1) )**2 )
     4                       * 0.125
         endif ! xtotnrg
         dlo(i,j,k) = d(i,j,k)
10      continue
20     continue
30    continue
c
      return
      end
