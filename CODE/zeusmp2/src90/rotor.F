c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 SOLAR WIND                \\\\\\\\\\
c
c                            Developed by
c                            Asif-ud-Doula
c=======================================================================
c
       subroutine rotor
c
c    magbrake <----------- initialises magnetic braking
c                                                       August 5th, 2004 
c
c    written by: Asif-ud-Doula
c
c  PURPOSE: Sets up an aligned magnetic rotator problem: Mouschovias and Paleologou 1980
c           Repeat of Stone (1992) MHD test for Zeus2D
c	    The problem is set up for cartesian coordinates: XYZ
c	    but it is equvalent to cylindrical coordinate solution.
c
c  Modified 1: relocated "implicit NONE" statements in subroutines
c              rotor and rotorres (11/21/05)
c
c-----------------------------------------------------------------------
c
      use real_prec
      use config
      use param
      use field
      use bndry
      use grid
      use root
      use scratch
      use cons
      use gravmod
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer  :: i, j, k, ip1, jp1, kp1, m,km1,jm1,kone,i_init
      real(rl) :: rhox,alpha,b0,zdisk,omega0 
c
       namelist / pgen     / rhox,alpha,b0,zdisk,omega0,i_init
c
c=======================================================================
c
      if(ldimen.eq.3) then 
       kone=1
      else
       kone=0
      endif
c
      i_init = 1    ! initial conditions flag; 1 = discontinous,
c                   !                          2 = continuous
      rhox   = 1.0  ! density of the ambient medium
      alpha  = 10.0 ! density of the disk
      b0     = 1.0  ! magnetic field along 1-directio, perpendicular to the disk
      zdisk  = 1.0  ! height/width of the disk
      omega0 = 1.0  ! rotaion rate of the disk, uniform
c
      if (myid .eq. 0) then
       read (1,pgen)
       write(2,pgen)
#ifdef MPI_USED
       buf_in( 1) = rhox   
       buf_in( 2) = alpha 
       buf_in( 3) = b0 
       buf_in( 4) = zdisk 
       buf_in( 5) = omega0
       ibuf_in(1) = i_init
#endif
      endif
#ifdef MPI_USED
      call MPI_BCAST( buf_in, 5, MPI_FLOAT
     &                , 0, comm3d, ierr )
      call MPI_BCAST(ibuf_in, 1, MPI_INTEGER,
     .                  0, comm3d, ierr)
#endif
#ifdef MPI_USED
      if (myid .ne. 0) then
       print *, 'I received the message', buf_in(1), myid
       rhox  = buf_in( 1)
       alpha  = buf_in( 2)
       b0  = buf_in( 3)
       zdisk  = buf_in( 4)
       omega0  = buf_in( 5)
       i_init = ibuf_in(1)
      endif ! myid
#endif
      do k=ks-2*kone,ke+3*kone
       do j=js-2,je+3
        do i=is-2,ie+3
         v1(i,j,k)=0.0
         v2(i,j,k)=0.0
         v3(i,j,k)=0.0
         d(i,j,k) =rhox 
         if(xmhd) then
          b1(i,j,k)=b0
          b2(i,j,k)=0.0
          b3(i,j,k)=0.0
         endif
c
c-----------------------------------------------------------------------
c     put in the disk, lies along 2-dir
c-----------------------------------------------------------------------
c
         if(abs(x1a(i)).le.zdisk) then
          d(i,j,k) = alpha
          if(i_init .eq. 1) then
           v3(i,j,k) = omega0*x2b(j) ! DIC
          else
           v3(i,j,k) = omega0*x2b(j)*0.5*(1.0 + cos(pi*x1b(i)/zdisk)) ! CIC
          endif
         endif
        enddo
       enddo
      enddo
c
      if(i_init .eq. 1) call dic(tlim,alpha,omega0,x2b(js))
      if(i_init .eq. 2) call cic(tlim,alpha,omega0,x2b(js))
c
      return
      end
c
c===================================================================
        subroutine rotorres
c
      use real_prec
      use config
      use param
      use field
      use bndry
      use grid
      use root
      use scratch
      use cons
      use gravmod
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer  :: i, j, k, ip1, jp1, kp1, m,km1,jm1
       real(rl) :: rhox,alpha,b0,zdisk,omega0
c

       namelist / pgen     / rhox,alpha,b0,zdisk,omega0
c
c-----------------------------------------------------------------------
c
          rhox=1.
          alpha=10.
          b0=1.0
          zdisk=1.0
          omega0=1.0
        if(myid.eq.0) then
#ifdef MPI_USED
         buf_in( 1) = rhox
         buf_in( 2) = alpha
         buf_in( 3) = b0
         buf_in( 4) = zdisk
         buf_in( 5) = omega0
c        ibuf_in( 1) = m
        call MPI_BCAST( buf_in, 5, MPI_FLOAT
     &                , 0, comm3d, ierr )
#endif
       endif
#ifdef MPI_USED
        if (myid .ne. 0) then
         print *, 'I received the message', buf_in(1), myid
         rhox  = buf_in( 1)
         alpha  = buf_in( 2)
         b0  = buf_in( 3)
         zdisk  = buf_in( 4)
         omega0  = buf_in( 5)
        endif ! myid
#endif
            return
            end
c===================================================================
      subroutine dic(t,rho,o0,r)
c
      use real_prec
      use config
      use param, ONLY : in
      use root,  ONLY : tlim
      use grid,  ONLY : is, ie, x1b
c
      implicit NONE
c
      integer :: i, n, m
c
      real*8  :: t, rho, o0, r, omega(in), bphi(in),
     .           sr, z, l1, l2, l3, dz
c
      sr  = sqrt(rho)
c
c-----------------------------------------------------------------------
c     do the case rho > 1.0
c-----------------------------------------------------------------------
c
      if(rho .gt. 1.0) then
c
       do i = is, ie
c
        z = x1b(i)
        if(z .le. 1.0) then
         if(t .lt. sr*(1.0-z)) then
          omega(i) = o0
          bphi (i) = 0.0
         else
          n = 0
10        continue
          l1 = sr*(2.0*float(n) + 1.0 - z)
          l2 = sr*(2.0*float(n) + 1.0 + z)
          l3 = sr*(2.0*float(n) + 3.0 - z)
          if(t .ge. l1 .and. t .lt. l2) then
           omega(i) = o0*sr/(sr+1.0)*((sr-1.0)/(sr+1.0))**n
           bphi (i) = -r*o0*sr/(sr+1.0)*((sr-1.0)/(sr+1.0))**n
           go to 20
          endif
          if(t .ge. l2 .and. t .lt. l3) then
           omega(i) = o0*((sr-1.0)/(sr+1.0))**(n+1)
           bphi (i) = 0.0
           go to 20
          endif
          n = n + 1
          go to 10
20       continue
        endif
       else ! Z
        if(t .lt. z-1.0) then
         omega(i) = 0.0
        else ! t
         m = 0
30       continue
         l1 = 2.0*float(m  )*sr + z - 1.0
         l2 = 2.0*float(m+1)*sr + z - 1.0
         if(t .ge. l1 .and. t .lt. l2) then
          omega(i) = o0*sr/(sr+1.0)*((sr-1.0)/(sr+1.0))**m
          go to 40
         endif
         m = m + 1
         go to 30
40       continue
        endif ! t
        bphi(i) = -r*omega(i)
       endif ! Z
c
       enddo ! i
c
      endif ! rho
c
c-----------------------------------------------------------------------
c     do the case rho = 1.0
c-----------------------------------------------------------------------
c
      if(rho .eq. 1.0) then
c
       do i = is, ie
c
        z = x1b(i)
        if(z .le. 1.0) then
         if(t .lt. 1.0-z) then
          omega(i) = o0
          bphi (i) = 0.0
         endif
         if(t .ge. 1.0-z .and. t .lt. 1.0+z) then
          omega(i) =  0.5*o0
          bphi (i) = -0.5*r*o0
         endif
         if(t .ge. 1.0+z) then
          omega(i) = 0.0
          bphi (i) = 0.0
         endif
        else ! z
         if(t .lt. z-1.0) then
          omega(i) = 0.0
          bphi (i) = 0.0
         endif
         if(t .ge. z-1.0 .and. t .lt. z+1.0) then
          omega(i) =  0.5*o0
          bphi (i) = -0.5*r*o0
         endif
         if(t .ge. z+1.0) then
          omega(i) = 0.0
          bphi (i) = 0.0
         endif
        endif ! z
       enddo ! i
      endif ! rho
c
      open(666,file='dic.out')
      write(666,"('     x1b        Omega        Bphi')")
      write(666,"(1p3d12.4)")(x1b(i),omega(i),bphi(i),i=is,ie)
      close(666)
c
      return
      end

c===================================================================
      subroutine cic(t,rho,o0,r)
c
      use real_prec
      use config
      use param, ONLY : in, pi
      use root,  ONLY : tlim
      use grid,  ONLY : is, ie, x1b
c
      implicit NONE
c
      integer :: i, n, m
c
      real*8  :: t, rho, o0, r, beta, omega(in), bphi(in),
     .           sr, z, i1, i2, i3, i4, i5, j1, j2, j3, j4,
     .           k1, k2, k3, l5
c
      sr   = sqrt(rho)
      beta = 1.0
c
c-----------------------------------------------------------------------
c     do the case rho > 1.0
c-----------------------------------------------------------------------
c
      if(rho .gt. 1.0) then
c
       do i = is, ie
        z = x1b(i)
        if(z .le. 1.0) then
         call geti(t,rho,z,beta,i1,i2,i3,i4,i5)
         l5       = -sin(pi*z/beta)*sin(pi*t/(sr*beta))
         omega(i) = 0.5*o0*(i1+i2+i3+i4+i5)
         bphi (i) = r*sr*0.5*o0*(i1-i2+i3-i4+l5)
        endif ! z
        if(z .gt. 1.0 .and. z .le. beta) then
         call getj(t,rho,z,beta,j1,j2,j3,j4)
         l5       = -sin(pi*z/beta)*sin(pi*t/beta)
         omega(i) = 0.5*o0*(j1+j2+j3+j4)
         bphi (i) = r*0.5*o0*(j1-j2-j3+l5)
        endif ! z
        if(z .gt. beta) then
         call getk(t,rho,z,beta,k1,k2,k3)
         omega(i) = 0.5*o0*(k1+k2+k3)
         bphi (i) = -r*omega(i)
        endif ! z
       enddo ! i
c
      endif ! rho
c
      open(666,file='cic.out')
      write(666,"(1p3d12.4)")(x1b(i),omega(i),bphi(i),i=is,ie)
      close(666)
c
      return
      end
c
      subroutine geti(t,rho,z,beta,i1,i2,i3,i4,i5)
c
      use real_prec
      use param, ONLY : pi
c
      implicit NONE
c
      real(rl) :: t, rho, z, beta, i1, i2, i3, i4, i5
c
      integer  :: ii1, ii2, ii3, ii4, ii
      real(rl) :: sr, bd, br, x1, x2, x3, x4, l1, l2, l3,
     .            dd, pp, rr, xx, eqn42b, eqn45b
c
      eqn42b(dd,pp,rr,xx,ii) =
     .   cos(pp*pi/beta)/dd*(cos(pi*(xx+1.0)/beta)
     .  -     rr**ii*cos(pi*(xx+1.0-2.0*pp*float(ii))/beta))
     .  +     pp*sin(pp*pi/beta)/dd*(sin(pi*(xx+1.0)/beta)
     .  -     rr**ii*sin(pi*(xx+1.0-2.0*pp*float(ii))/beta))
     .  - 0.5*cos(pi*xx/(pp*beta))
     .  + 0.5*rr**ii*cos(pi*(xx-2.0*pp*float(ii))/(pp*beta))
c
      eqn45b(dd,pp,rr,xx,ii) =
     .   cos(pp*pi/beta)/dd*(cos(pi*xx/beta)
     .  -    rr**ii*cos(pi*(xx-2.0*pp*float(ii))/beta))
     .  +    pp*sin(pp*pi/beta)/dd*(sin(pi*xx/beta)
     .  -    rr**ii*sin(pi*(xx-2.0*pp*float(ii))/beta))
     .  + 0.5*(rr**ii-1.0)
c
      sr = sqrt(rho)
      bd = 2.0*((cos(sr*pi/beta))**2 + rho*(sin(sr*pi/beta))**2)
      br = (sr - 1.0)/(sr + 1.0)
c
      if(t .lt. sr*(1.0-z)) then
       i1 = 0.0
      else
       x1 = t + sr*z
       ii1 = 1
10     continue
       l1 = t - sr*(1.0-z)
       l2 = 2.0*sr*float(ii1)
       l3 = t + sr*(1.0+z)
       if(l2 .gt. l1 .and. l2 .le. l3) go to 20
       ii1 = ii1 + 1
       go to 10
20     continue
       i1 = eqn42b(bd,sr,br,x1,ii1)
      endif
c
      if(t .lt. sr*(1.0+z)) then
       i2 = 0.0
      else
       x2 = t - sr*z
       ii2 = 1
 30    continue
       l1 = t - sr*(1.0+z)
       l2 = 2.0*sr*float(ii2)
       l3 = t+sr*(1.0-z)
       if(l2 .gt. l1 .and. l2 .le. l3) go to 40
       ii2 = ii2 + 1
       go to 30
 40    continue
       i2 = eqn42b(bd,sr,br,x2,ii2)
      endif
c
      if(t .lt. sr*(1.0-z)+(beta-1.0)) then
       i3 = 0.0
      else
       x3 = t + sr*z - beta+1.0
       ii3 = 1
 50    continue
       l1 = t - sr*(1.0-z) - (beta-1.0)
       l2 = 2.0*sr*float(ii3)
       l3 = t + sr*(1.0+z) - (beta-1.0)
       if(l2 .gt. l1 .and. l2 .le. l3) go to 60
       ii3 = ii3 + 1
       go to 50
 60    continue
       i3 = eqn45b(bd,sr,br,x3,ii3)
      endif
c
      if(t. lt. sr*(1.0+z)+(beta-1.0)) then
       i4 =  0.0
      else
       x4 = t - sr*z - beta + 1.0
       ii4 = 1
 70    continue
       l1 = t -sr*(1.0+z) - (beta-1.0)
       l2 = 2.0*sr*float(ii4)
       l3 = t + sr*(1.0-z) - (beta-1.0)
       if(l2 .gt. l1 .and. l2 .le. l3) go to 80
       ii4 = ii4 + 1
       go to 70
 80    continue
       i4 = eqn45b(bd,sr,br,x4,ii4)
      endif
c
      i5 = 1.0 + cos(pi*z/beta)*cos(pi*t/(sr*beta))
c
      return
      end
c
      subroutine getj(t,rho,z,beta,j1,j2,j3,j4)
c
      use real_prec
      use param, ONLY : pi
c
      real(rl) :: t, rho, z, beta, j1 j2, j3, j4,
     .            sr, dd, rr, yy, ff, xx, l1, l2, l3
c
      integer  :: jj
c
      sr = sqrt(rho)
      dd = 2.0*((cos(sr*pi/beta))**2 + rho*(sin(sr*pi/beta))**2)
      rr = (sr-1.0)/(sr+1.0)
      ff = 2.0*sr/(sr+1.0)
c
      if(t .lt. beta-z) then
       j1 = 0.0
      else
       j1 = -1.0*(cos(pi*(t+z)/(2.0*beta)))**2
      endif
c
      if(t .lt. z-1.0) then
       j2 = 0.0
      else
       yy = t - z + 1.0
       jj = 1
10     continue
       l1 = t-z+1.0
       l2 = 2.0*sr*float(jj)
       l3 = t-z+1.0+2.0*sr
       if(l2 .gt. l1 .and. l2 .le. l3) go to 20
       jj = jj + 1
       go to 10
20     continue
       xx = pi/beta
       j2 = -2.0/dd*(sin(xx)*cos(sr*xx) - sr*cos(xx)*sin(sr*xx))
     .             *(cos(sr*xx)*sin(yy*xx) - sr*sin(sr*xx)*cos(yy*xx))
     .   -   rr**(jj-1)*ff*(
     .         (sr*sin(xx)*sin(sr*xx)
     .          + cos(xx)*cos(sr*xx))*cos(xx*(yy-2.0*sr*float(jj)+sr))/dd
     .       + (sr*cos(xx)*sin(sr*xx)
     .          - sin(xx)*cos(sr*xx))*sin(xx*(yy-2.0*sr*float(jj)+sr))/dd
     .   - 0.5*cos(xx*(yy/sr-2.0*float(jj)+1.0)))
      endif
c
      if(t .lt. z+beta-2.0) then
       j3 = 0.0
      else
       yy = t-z-beta+2.0
       jj = 1
30     continue
       l1 = t-z-beta+2.0
       l2 = 2.0*sr*float(jj)
       l3 = t-z-beta+2.0+2.0*sr
       if(l2 .gt. l1 .and. l2 .le. l3) go to 40
       jj = jj + 1
       go to 30
40     continue
       xx = pi/beta
       j3 = (((cos(sr*xx))**2-rho*(sin(sr*xx))**2)*cos(xx*yy)
     .       + 2.0*sr*sin(sr*xx)*cos(sr*xx)*sin(xx*yy))/dd
     .  - 0.5 + rr**(jj-1)*ff*
     .    (0.5 - sr/dd*sin(sr*xx)*sin(xx*(yy-2.0*sr*float(jj)+sr))
     .         - cos(sr*xx)/dd*cos(xx*(yy-2.0*sr*float(jj)+sr)))
      endif
c
      j4 = 1.0 + cos(pi*z/beta)*cos(pi*t/beta)
c
      return
      end
c
      subroutine getk(t,rho,z,beta,k1,k2,k3)
c
      use real_prec
      use param, ONLY : pi
c
      implicit NONE
c
      integer  :: jj
      real(rl) :: t, rho, z, beta, k1, k2, k3,
     .            sr, dd, rr, yy, ff, xx, l1, l2, l3
c
      sr = sqrt(rho)
      dd = 2.0*((cos(sr*pi/beta))**2 + rho*(sin(sr*pi/beta))**2)
      rr = (sr-1.0)/(sr+1.0)
      ff = 2.0*sr/(sr+1.0)
c
      if(t .lt. z-1.0) then
       k1 = 0.0
      else
       yy = t - z + 1.0
       jj = 1
10     continue
       l1 = t-z+1.0
       l2 = 2.0*sr*float(jj)
       l3 = t-z+1.0+2.0*sr
       if(l2 .gt. l1 .and. l2 .le. l3) go to 20
       jj = jj + 1
       go to 10
20     continue
       xx = pi/beta
       k1 = -2.0/dd*(sin(xx)*cos(sr*xx) - sr*cos(xx)*sin(sr*xx))
     .             *(cos(sr*xx)*sin(yy*xx) - sr*sin(sr*xx)*cos(yy*xx))
     .   -   rr**(jj-1)*ff*(
     .       (sr*sin(xx)*sin(sr*xx)
     .   + cos(xx)*cos(sr*xx))*cos(xx*(yy-2.0*sr*float(jj)+sr))/dd
     .       + (sr*cos(xx)*sin(sr*xx)
     .   - sin(xx)*cos(sr*xx))*sin(xx*(yy-2.0*sr*float(jj)+sr))/dd
     .   - 0.5*cos(xx*(yy/sr-2.0*float(jj)+1.0)))
      endif
c
      if(t .lt. z+beta-2.0) then
       k2 = 0.0
      else
       yy = t-z-beta+2.0
       jj = 1
30     continue
       l1 = t-z-beta+2.0
       l2 = 2.0*sr*float(jj)
       l3 = t-z-beta+2.0+2.0*sr
       if(l2 .gt. l1 .and. l2 .le. l3) go to 40
       jj = jj + 1
       go to 30
40     continue
       xx = pi/beta
       k2 = (((cos(sr*xx))**2-rho*(sin(sr*xx))**2)*cos(xx*yy)
     .       + 2.0*sr*sin(sr*xx)*cos(sr*xx)*sin(xx*yy))/dd
     .  - 0.5 + rr**(jj-1)*ff*
     .    (0.5 - sr/dd*sin(sr*xx)*sin(xx*(yy-2.0*sr*float(jj)+sr))
     .         - cos(sr*xx)/dd*cos(xx*(yy-2.0*sr*float(jj)+sr)))
      endif
c
      if(t .lt. z-beta) then
       k3 = 0.0
      else
       k3 = (sin(pi*(t-z+beta)/(2.0*beta)))**2
      endif
c
      return
      end
