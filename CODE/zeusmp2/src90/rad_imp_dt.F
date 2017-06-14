c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////           R A D _ I M P _ D T             \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c=======================================================================
c
       subroutine rad_imp_dt(erold,ernew) 
c
c  INPUT
c    erold   radiation energy density at old time step
c    ernew   updated radiation energy density
c
c  OUTPUT
c    dtimrdi  proposed new timestep
c
c  EXTERNALS:
c
c  AUTHOR
c    Robert Fiedler and John Hayes
c
      use real_prec
      use config
      use param
      use field
      use grid
      use root
      use radiation
      use scratch
      use cons
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer  :: i    , j    , k,
     1            imax , jmax , kmax,
     &            iamax, jamax, kamax
c
      real(rl) :: erold(in,jn,kn) , ernew(in,jn,kn)
      real(rl) :: egold(in,jn,kn) , egnew(in,jn,kn)
      real(rl) :: dtec, er_max, ader, ader_max, fac_dt,
     .            eg_max, adeg, adeg_max, ade_max
c
      real(rl) :: delta_t_max
      real(rl) :: inp(2), out(2)
      real(rl) :: trnew, trold
c
      data fac_dt /1.0592537/
c
c-----------------------------------------------------------------------
c
c Find the minimum time step required by the Courant condition for
c this tile.
c
c
c-----------------------------------------------------------------------
c
c Find max new rad energy and max abs rad energy change for this tile.
c
      iamax = is
      jamax = js
      kamax = ks
      ader_max = abs((ernew(is,js,ks)-erold(is,js,ks))/
     .                (erold(is,js,ks)+tiny))
      do k=ks,ke
       do j=js,je
        do i=is,ie
         ader = abs((ernew(i,j,k) - erold(i,j,k))/
     .              (erold(i,j,k)+tiny))
         if (ader .gt. ader_max) then
          ader_max = ader
          iamax = i
          jamax = j
          kamax = k
         endif
        enddo ! i
       enddo ! j
      enddo ! k
#ifdef MPI_USED
c
c Now find the largest e_max and ade_max among all tiles, and send the 
c results to all in buf_out.
c
      buf_in(1) = ader_max
      buf_in(2) = float(myid)
      call MPI_ALLREDUCE(buf_in, buf_out, 1, MPI_2DOUBLE_PRECISION,
     .                    MPI_MAXLOC, comm3d, ierr)
      if(.false.) then
       if(mod(nhy,10) .eq. 0) then
        if(myid .eq. int(buf_out(2))) then
         write(*,"('RAD_IMP_DT: max dE/E = ',1pd12.4,' at i,j,k =',
     .             3i3,' on node ',i2)")
     .         buf_out(1),iamax,jamax,kamax,int(buf_out(2))
         write(*,"('erold, ernew = ',1p2d12.4)")
     .       erold(iamax,jamax,kamax),
     .       ernew(iamax,jamax,kamax)
        endif
       endif
      endif
      ader_max = buf_out(1)
#else /* MPI_USED */
      if(.false.) then
       if(mod(nhy,10) .eq. 0) then
        write(*,"('RAD_IMP_DT: max dE/E = ',1pd12.4,' at i,j,k =',
     .               3i3)")
     .           ader_max,iamax,jamax,kamax
        write(*,"('erold, ernew = ',1p2d12.4)")
     .         erold(iamax,jamax,kamax),
     .         ernew(iamax,jamax,kamax)
       endif
      endif
#endif /* MPI_USED */
c
c Convert dtec from input dt
c
      dtec = 1.0/max(sqrt(dtimrdi2),tiny)
c
c Error test
c
      ade_max = ader_max
      if (ade_max .le. epsmaxd) then
       dtec = dtec * fac_dt
      else
       dtec = dtec / fac_dt
      endif
      dtimrdi2 = 1.0/max(dtec,tiny)**2
c
      return
      end
