c=======================================================================
c
c                            Developed by
c                Laboratory for Computational Astrophysics
c                  University of California at San Diego
c
      subroutine gravity
c
c-----------------------------------------------------------------------
c     driver routine which selects an algorithm for computing the 
c     gravitational potential, gp(i,j,k), based upon user-specified
c     values of XGRAV, XSPHGRV, XGRVFFT, LDIMEN, and LGEOM:
c
c     General potentials:
c
c     LGEOM = 1 (XYZ)
c      (1) Dirichlet/Neuman BC's --> use grav3D_MG
c      (2) Triply-periodic  BC's --> use fftwgrav
c
c     LGEOM = 2 (ZRP):
c      (1) LDIMEN = 2 --> use grav2D_CG
c      (2) LDIMEN = 3 --> use grav3D_CG
c
c     LGEOM = 3 (RTP):
c      (1) LDIMEN = 2 --> use grav2D_CG
c      (2) LDIMEN = 3 --> use grav3D_CG
c
c     Spherically symmetric potential (GM/r):
c
c     LGEOM = 3 (RTP)
c      (1) LDIMEN = 1 --> use spherical_gravity
c      (2) LDIMEN = 2 --> use spherical_gravity
c-----------------------------------------------------------------------
c
      use real_prec
      use config
      use param
      use cons
      use domain
      use root
      use field
      use grid
      use bndry
      use impsoln
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
c-----------------------------------------------------------------------
c     spherical potential (GM/r); supported if LGEOM = 3 and 
c     LDIMEN = 1 or 2
c-----------------------------------------------------------------------
c
      if(xsphgrv) then
c
c     check for 1D or 2D...
c
       if(ldimen .eq. 3) then
        if(myid .eq. 0) then
         write(*,"('XSPHGRV: LDIMEN must equal 1 or 2!')")
         write(*,"('aborting run...')")
        endif
        go to 1
       endif
c
c     check for RTP geometry...
c
       if(lgeom .ne. 3) then
        if(myid .eq. 0) then
         write(*,"('XSPHGRV: LGEOM must equal 3!')")
         write(*,"('aborting run...')")
        endif
        go to 1
       endif
c
       call spherical_gravity
       return
c
1      continue
#ifdef MPI_USED
       call mpi_finalize(ierr)
#endif
       stop
      endif
c
c-----------------------------------------------------------------------
c     potential for triply-periodic 3D problems (LGEOM =1)
c-----------------------------------------------------------------------
c
#ifdef FFT
      if(xgrvfft) then
c
c     check for 3D...
c
       if(ldimen .ne. 3) then
        if(myid .eq. 0) then
         write(*,"('XGRVFFT: LDIMEN must equal 3!')")
         write(*,"('aborting run...')")
        endif
        go to 2
       endif
       if(lgeom .ne. 1) then
c
c     check for XYZ geometry...
c
        if(myid .eq. 0) then
         write(*,"('XGRVFFT: LGEOM must equal 1!')")
         write(*,"('aborting run...')")
        endif
        go to 2
       endif
c
c     check for periodic BC's on all axes...
c
       if(periodic(1) .eqv. .false. .or. periodic(2) .eqv. .false. .or.
     .    periodic(3) .eqv. .false.) then
        if(myid .eq. 0) then
         write(*,"('XGRVFFT: problem must be triply periodic!')")
         write(*,"('aborting run...')")
        endif
        go to 2
       endif
c
       call fftwgrav
       return
c
2      continue
#ifdef MPI_USED
       call mpi_finalize(ierr)
#endif
       stop
      endif
#endif /* FFT */
c
c-----------------------------------------------------------------------
c     Solve Poisson's equation for 2D/3D problems which are not
c     triply-periodic
c-----------------------------------------------------------------------
c
      if(ldimen .eq. 2) then
c
c     2D, so must be spherical or cylindrical
c
       if(lgeom .eq. 1) then
        if(myid .eq. 0) then
         write(*,"('XGRAV: Can not do cartesian gravity in 2D!')")
         write(*,"('aborting run...')")
        endif
        go to 3
       endif
c
       call grav2D_CG
       return
c
3      continue
#ifdef MPI_USED
       call mpi_finalize(ierr)
#endif
       stop
      endif
c
c-----------------------------------------------------------------------
c     In theory, any geometry is legal here
c-----------------------------------------------------------------------
c
      if(ldimen .eq. 3) then
       if(lgeom .eq. 1) then
#ifdef USE_MGMPI
        call grav3D_MG
#else
        call grav3D_CG
#endif
       else
        call grav3D_CG
       endif
      endif
c
      return
      end
c
c=======================================================================
      subroutine spherical_gravity
c
c     computes M(r) for use in forces.F in the case of 1D or 2D 
c     spherical problems, in which spherical gravity is assumed
c     in place of computing a potential
c
c     Written by: John Hayes
c
c     Modified 1: 12/12/05, fixed typo in "dml" expression
c     for 2D case when xwedge=false.  (deleted a "dcos2" reference).
c
      use real_prec
      use config
      use param
      use cons
      use field
      use grid
      use bndry
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
      integer  :: i, j, k, index
c
      real(rl), dimension(:,:), allocatable :: dml, dmu
      real(rl), dimension(:  ), allocatable :: intm_half
c
c
      real(rl) :: local_mass, theta1, theta2, alpha1, alpha2, 
     .            dcos1, dcos2
c
      real(rl) :: third, halfpi, onep5pi
c
      allocate(dml (in,jn))
      allocate(dmu (in,jn))
c
      do i = 1, in
       intm(i) = 0.0D0
      enddo
c
      third   = 1.0D0/3.0D0
      halfpi  = 0.5D0*pi
      onep5pi = 1.5D0*pi
c
      allocate(intm_half(2*in))
c
c-----------------------------------------------------------------------
c     differential mass elements (lower half cell and upper half cell)
c-----------------------------------------------------------------------
c
      if(ldimen .eq. 2) then
       if(xwedge) then
        do j = js, je
         theta1 = x2a(j)
         theta2 = x2a(j+1)
         if(x2a(j+1) .le. halfpi) then
          alpha1 = halfpi - x2a(j)
          alpha2 = halfpi - x2a(j+1)
          dcos1    = cos(theta1) - cos(theta2)
          dcos2    = cos(alpha2) - cos(alpha1)
         else
          alpha1 = onep5pi - x2a(j)
          alpha2 = onep5pi - x2a(j+1)
          dcos1    = cos(theta1) - cos(theta2)
          dcos2    = cos(alpha2) - cos(alpha1)
         endif
         do i = is, ie
          dmu(i,j) = 2.0D0*third*pi*d(i,j,ks)*(x1a(i+1)**3 - x1b(i)**3)
     .             * (dcos1 + dcos2)
         enddo
         if(coords(1) .eq. 0) then
          do i = is+1, ie
           dml(i,j) = 2.0D0*third*pi*d(i,j,ks)*(x1b(i)**3 - x1a(i)**3)
     .              * (dcos1 + dcos2)
          enddo
          dml(is,j) = 2.0D0*third*pi*d(is,j,ks)*x1b(is)**3
     .             * (dcos1 + dcos2)
         else ! coords(1)
          do i = is, ie
           dml(i,j) = 2.0D0*third*pi*d(i,j,ks)*(x1b(i)**3 - x1a(i)**3)
     .              * (dcos1 + dcos2)
          enddo
         endif ! coords(1)
        enddo
       else ! xwedge
        do j = js, je
         theta1 = x2a(j)
         theta2 = x2a(j+1)
         dcos1    = cos(theta1) - cos(theta2)
         do i = is, ie
          dmu(i,j) = 2.0D0*third*pi*d(i,j,ks)*(x1a(i+1)**3 - x1b(i)**3)
     .             * dcos1
         enddo
         if(coords(1) .eq. 0) then
          do i = is+1, ie
           dml(i,j) = 2.0D0*third*pi*d(i,j,ks)*(x1b(i)**3 - x1a(i)**3)
     .              * dcos1
          enddo
          dml(is,j) = 2.0D0*third*pi*d(is,j,ks)*x1b(is)**3
     .              * dcos1
         else ! coords(1)
          do i = is, ie
           dml(i,j) = 2.0D0*third*pi*d(i,j,ks)*(x1b(i)**3 - x1a(i)**3)
     .              * dcos1
          enddo
         endif ! coords(1)
        enddo
       endif ! xwedge
      else ! ldimen
       do i = is, ie
        dmu(i,js) = 4.0D0*third*pi*d(i,js,ks)*(x1a(i+1)**3 - x1b(i)**3)
       enddo
       if(coords(1) .eq. 0) then
        do i = is+1, ie
         dml(i,js) = 4.0D0*third*pi*d(i,js,ks)*(x1b(i)**3 - x1a(i)**3)
        enddo
        dml(is,js) = 4.0D0*third*pi*x1b(is)**3*d(is,js,ks)
       else
        do i = is, ie
         dml(i,js) = 4.0D0*third*pi*d(i,js,ks)*(x1b(i)**3 - x1a(i)**3)
        enddo
       endif ! coords(1)
      endif ! ldimen
c
c-----------------------------------------------------------------------
c     total interior mass (local processor) contained between x1a(is)
c     and x1a(ie+1), computed at every half cell interface
c-----------------------------------------------------------------------
c
      do i = 1, 2*in
       intm_half(i) = 0.0D0
      enddo
      local_mass = 0.0D0
c
      index = 0
      do i = is, ie
       index = index + 2
       do j = js, je
        local_mass = local_mass + dml(i,j)
       enddo
       intm_half(index  ) = local_mass
       do j = js, je
        local_mass = local_mass + dmu(i,j)
       enddo
       intm_half(index+1) = local_mass
      enddo
c
#ifdef MPI_USED
c
c-----------------------------------------------------------------------
c     sum mass arrays across rows of processors (2-coord)
c-----------------------------------------------------------------------
c
      call MPI_BARRIER(comm3d, ierr)
c
      if(ntiles(2) .gt. 1) call sum_over_two_coord(intm_half,local_mass)
c
c-----------------------------------------------------------------------
c     sum mass arrays up 1-coordinate
c-----------------------------------------------------------------------
c
      call MPI_BARRIER(comm3d, ierr)
c
      if(ntiles(1) .gt. 1) call sum_over_one_coord(intm_half,local_mass)
c
#endif /* MPI_USED */
c
c-----------------------------------------------------------------------
c     pick out interior mass at cell FACES
c-----------------------------------------------------------------------
c
      index = -1
      do i = is, ie+1
       index = index + 2
       intm(i) = intm_half(index)
      enddo
c
      deallocate(intm_half)
      deallocate(dml)
      deallocate(dmu)
c
      return
      end
c
c=======================================================================
      subroutine sum_over_two_coord(intm_half,local_mass)
c
#ifdef MPI_USED
      use real_prec
      use param
      use grid
      use mpiyes
      use mpipar
c
      implicit NONE
c
      real(rl) :: intm_half(2*in)
      real(rl), dimension(:), allocatable :: intm_buf
      real(rl) :: local_mass, loc_buf
c
      integer source, target, tag, i, l, nhalfcells, count
c
      count = 2*in
      allocate(intm_buf(count))
c
      if(coords(2) .ne. 0) then
c
c-----------------------------------------------------------------------
c     process with (coords(1),coords(2)) = (N,M) sends data to 
c     process with (coords(1),coords(2)) = (N,0) 
c-----------------------------------------------------------------------
c
       target = coords(1)*ntiles(2)
       tag    = 1000+coords(1)
       call MPI_SEND(intm_half(1), count,
     .               MPI_FLOAT,
     .               target, tag, comm3d, ierr)
      else  ! coords(2)
c
c-----------------------------------------------------------------------
c     process with coords(2) = 0 computes theta-summed values of local
c     interior mass
c-----------------------------------------------------------------------
c
       tag    = 1000+coords(1)
       do l = 1, ntiles(2)-1
       target = myid+l
        call MPI_RECV(intm_buf(1), count,
     .               MPI_FLOAT,
     .               target, tag, comm3d, stat, ierr)
        do i = 1, 2*in
         intm_half(i) = intm_half(i) + intm_buf(i)
        enddo ! i
       enddo ! l
      endif ! coords(2)
c
      call MPI_BARRIER(comm3d, ierr)
c
      if(coords(2) .eq. 0) then
c
c-----------------------------------------------------------------------
c     process   with (coords(1),coords(2)) = (N,0) sends interior mass to
c     processes with (coords(1),coords(2)) = (N,M)
c-----------------------------------------------------------------------
c
       do l = 1, ntiles(2)-1
        target = coords(1)*ntiles(2) + l
        tag    = 2000 + l
        call MPI_SEND(intm_half(1), count,
     .               MPI_FLOAT,
     .               target, tag, comm3d, ierr)
       enddo ! l
      else  ! coords(2)
c
c-----------------------------------------------------------------------
c     processes with coords(2) > 0 receive incoming messages and update
c     local values of intm_half
c-----------------------------------------------------------------------
c
       source = coords(1)*ntiles(2)
       tag    = 2000 + coords(2)
       call MPI_RECV(intm_buf(1), count,
     .              MPI_FLOAT,
     .              source, tag, comm3d, stat, ierr)
       do i = 1, 2*in
        intm_half(i) = intm_buf(i)
       enddo
      endif ! coords(2)
c
c-----------------------------------------------------------------------
c     using theta-summed interior masses, each processor updates its
c     value of local_mass, which is now the total mass contained by all
c     processors having a common value of coords(1)
c-----------------------------------------------------------------------
c
      nhalfcells = 2*(ie-is+1)
      local_mass = intm_half(nhalfcells+1)
c
      deallocate(intm_buf)
c
#endif /* MPI_USED */
      return
      end
c=======================================================================
      subroutine sum_over_one_coord(intm_half,local_mass)
c
#ifdef MPI_USED
      use real_prec
      use param
      use grid
      use mpiyes
      use mpipar
c
      implicit NONE
c
      real(rl) :: intm_half(2*in)
      real(rl) :: local_mass, loc_buf
c
      integer source, target, tag, i, l, nhalfcells
c
c-----------------------------------------------------------------------
c     all processes with coords(1) = 0 begin upward cascading sum of
c     interior mass, passing data to all processes with same value of 
c     coords(2)
c-----------------------------------------------------------------------
c
      nhalfcells = 2*(ie-is+1)
c
      if(coords(1) .eq. 0) then
c
       tag = 1000*n1p
       call MPI_SEND(local_mass, 1,
     .               MPI_FLOAT,
     .               n1p, tag, comm3d, ierr) 
c
      else ! coords(1)
c
       tag = 1000*myid
       call MPI_RECV(loc_buf, 1,
     .               MPI_FLOAT,
     .               n1m, tag, comm3d, stat, ierr)
c
       do i = 1, nhalfcells+1
        intm_half(i) = intm_half(i) + loc_buf
       enddo
c
       local_mass = intm_half(nhalfcells+1)
c
       tag = 1000*n1p
c
       if(coords(1) .ne. ntiles(1)-1) then
        call MPI_SEND(local_mass, 1,
     .                MPI_FLOAT,
     .                n1p, tag, comm3d, ierr) 
       endif ! coords(1)
c
      endif ! coords(1)
c
#endif /* MPI_USED */
      return
      end
