c=======================================================================
c
c                            Developed by
c                Laboratory for Computational Astrophysics
c                  University of California at San Diego
c
      subroutine grav2D_CG
c
c     2-D gravitational potential solver using John Hayes's
c     CGSOLVE package for the covariant form of Poisson's equation
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
      real(rl), dimension(:,:,:,:,:), allocatable, save ::
     .          ad, a0, a1, a2
      real(rl), dimension(:,:,:,:), allocatable, save ::
     .          x, b
c
      real(rl) :: toler, error, glbtoler
c
      integer :: i, j, k, prev_tot
c
c-----------------------------------------------------------------------
c     First timestep only: initialize everything required for CGSOLVE
c-----------------------------------------------------------------------
c
      if(igcall .eq. 0) then
c
c-----------------------------------------------------------------------
c      allocate matrix, solution, and RHS array storage
c-----------------------------------------------------------------------
c
       allocate(ad(neqm,neqm,IN,JN,KN))
       allocate(a0(neqm,neqm,IN,JN,KN))
       allocate(a1(neqm,neqm,IN,JN,KN))
       allocate(a2(neqm,neqm,IN,JN,KN))
       allocate(x (     neqm,IN,JN,KN))
       allocate(b (     neqm,IN,JN,KN))
c
      endif ! igcall
c
c-----------------------------------------------------------------------
c     compute BC arrays
c
c     Uncomment the "call gpbv" line if external BC's are to be computed
c     from a multipole approximation.  NOTE: gpbv uses a low-order
c     approximation which is not accurate if the external boundary is
c     not far from the mass distribution.  If the external boundary
c     potential can be computed analytically (as in the case of spherical
c     collapse), then those values should be adopted as fois(13) and
c     input in zmp_inp in the BC namelists as appropriate.
c-----------------------------------------------------------------------
c
C      call gpbv
      do k = 1, kn
       do j = 1, jn
        gpoib(j,k,1) = fois(13)
       enddo
      enddo
c
c-----------------------------------------------------------------------
c     Load arrays for matrix
c-----------------------------------------------------------------------
c
      do k = 1, kn
      do j = 1, jn
      do i = 1, in
       ad(1,1,i,j,k) = 0.0D0
       a0(1,1,i,j,k) = 0.0D0
       a1(1,1,i,j,k) = 0.0D0
       a2(1,1,i,j,k) = 0.0D0
       b (  1,i,j,k) = 0.0D0
      enddo
      enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute diagonals along i faces
c-----------------------------------------------------------------------
c
      call compute_diagonals(ad,a0,a1,is,is,js,je,ks,ke)
      call compute_diagonals(ad,a0,a1,ie,ie,js,je,ks,ke)
c
c-----------------------------------------------------------------------
c     exchange i-face values with neighboring processors
c-----------------------------------------------------------------------
c
      nreq = 0
      nsub = 0
      call updt_mtrx_bnd_1(is, ie, js, je, ks, ke, a0)
      nsub = nsub + 1
      call updt_mtrx_bnd_1(is, ie, js, je, ks, ke, a1)
      nsub = nsub + 1
      call updt_mtrx_bnd_1(is, ie, js, je, ks, ke, ad)
c
c-----------------------------------------------------------------------
c     compute diagonals at interior zones
c-----------------------------------------------------------------------
c
      call compute_diagonals(ad,a0,a1,is+1,ie-1,js+1,je-1,ks,ke)
c
c-----------------------------------------------------------------------
c     wait for i-face MPI exchanges to finish
c-----------------------------------------------------------------------
c
#ifdef MPI_USED
      if(nreq .ne. 0) then 
       call mpi_waitall(nreq, req, stat, ierr)
       nreq = 0
       nsub = 0
      endif
#endif /* MPI_USED */
c
c-----------------------------------------------------------------------
c     compute diagonals along j faces
c-----------------------------------------------------------------------
c
      call compute_diagonals(ad,a0,a1,is+1,ie-1,js,js,ks,ke)
      call compute_diagonals(ad,a0,a1,is+1,ie-1,je,je,ks,ke)
c
c-----------------------------------------------------------------------
c     send j-face values to j-neighbor processes
c-----------------------------------------------------------------------
c
      call updt_mtrx_bnd_2(is, ie, js, je, ks, ke, a0)
      nsub = nsub + 1
      call updt_mtrx_bnd_2(is, ie, js, je, ks, ke, a1)
      nsub = nsub + 1
      call updt_mtrx_bnd_2(is, ie, js, je, ks, ke, ad)
c
c-----------------------------------------------------------------------
c     compute RHS
c-----------------------------------------------------------------------
c
      DO k=ks,ke
         DO j=js,je
            DO i=is,ie
               b(1,i,j,k) = 4.0D0*d(i,j,k)*
     .                       PI*guniv*dvl1a(i)*dvl2a(j)
            END DO
         END DO
      END DO
c
c-----------------------------------------------------------------------
c     wait for j-face exchanges to finish
c-----------------------------------------------------------------------
c
#ifdef MPI_USED
      if(nreq .ne. 0) then 
       call mpi_waitall(nreq, req, stat, ierr)
       nreq = 0
      endif
#endif /* MPI_USED */
c
c-----------------------------------------------------------------------
c     modify matrix elements for boundary conditions as needed
c-----------------------------------------------------------------------
c
      if(coords(1) .eq. 0) then
       do k = ks, ke
        do j = js, je
         if((niis(3) .eq. 1) .or. (niis(3) .eq. 2)) then
           ad(1,1,is,j,k) = ad(1,1,is,j,k) + 
     .                      dvl2a(j)*g2a(is  )*
     .                      g31a(is  )*dx1bi(is  )
         endif
         if(niis(3) .eq. 3) then 
           b(1,is,j,k) = b(1,is,j,k) + guniv*
     .     dvl2a(j)*g2a(is  )*g31a(is)*dx1bi(is)*gpiib(j,k,1)
         endif
        enddo ! j
       enddo ! k
      endif ! coords(1)
      if(coords(1) .eq. ntiles(1)-1) then
       do k = ks, ke
        do j = js, je
         a0(1,1,ie,j,k) = 0.0D0
         if((nois(3) .eq. 1) .or. (nois(3) .eq. 2)) then
           ad(1,1,ie,j,k) = ad(1,1,ie,j,k) + 
     .                      dvl2a(j)*g2a(ie+1)*
     .                      g31a(ie+1)*dx1bi(ie+1)
         endif
         if(nois(3) .eq. 3) then
           b(1,ie,j,k) = b(1,ie,j,k) + guniv*
     .     dvl2a(j)*g2a(ie+1)*g31a(ie+1)*dx1bi(ie+1)*
     .     gpoib(j,k,1)
         endif
        enddo ! j
       enddo ! k
      endif ! coords(1)
      if(coords(2) .eq. 0) then
       do k = ks, ke
        do i = is, ie
         if((nijs(3) .eq. 1) .or. (nijs(3) .eq. 2)) then
           ad(1,1,i,js,k) = ad(1,1,i,js,k) + 
     .                      dvl1a(i)*
     .                      g2bi(i)**2*g32a(js  )*dx2bi(js  )
         endif ! nijs
         if(nijs(3) .eq. 3) then
           b(1,i,js,k) = b(1,i,js,k) + guniv*
     .     dvl1a(i)*g2bi(i)**2*g32a(js  )*dx2bi(js  )*
     .     gpijb(i,k,1)
         endif ! nijs
        enddo ! i
       enddo ! k
      endif ! coords(2)
      if(coords(2) .eq. ntiles(2)-1) then
       do k = ks, ke
        do i = is, ie
         if(nojs(3) .ne. 4) a1(1,1,i,je,k) = 0.0D0
         if((nojs(3) .eq. 1) .or. (nojs(3) .eq. 2)) then
           ad(1,1,i,je,k) = ad(1,1,i,je,k) + 
     .                      dvl1a(i)*
     .                      g2bi(i)**2*g32a(je+1)*dx2bi(je+1)
         endif ! nojs
         if(nojs(3) .eq. 3) then
           b(1,i,je,k) = b(1,i,je,k) + guniv*
     .     dvl1a(i)*g2bi(i)**2*g32a(je+1)*dx2bi(je+1)*
     .     gpojb(i,k,1)
         endif ! nojs
        enddo ! i
       enddo ! k
      endif ! coords(2)
c
c-----------------------------------------------------------------------
c     initialize potential
c-----------------------------------------------------------------------
c
      if(igcall .eq. 0) then
       do k = 1, kn
        do j = 1, jn
         do i = 1, in
           x(1,i,j,k) = 0.0D0
         enddo
        enddo
       enddo
      else
       do k = ks, ke
        do j = js, je
         x(1,is,j,k) = -gp(is,j,k)
         x(1,ie,j,k) = -gp(ie,j,k)
        enddo
       enddo
       nreq = 0
       nsub = 0
       call updt_vec_bnd_1(is, ie, js, je, ks, ke, x)
       do k = ks, ke
        do i = is+1, ie-1
         x(1,i,js,k) = -gp(i,js,k)
         x(1,i,je,k) = -gp(i,je,k)
        enddo
       enddo
#ifdef MPI_USED
       if(nreq .ne. 0) then
        call mpi_waitall(nreq, req, stat, ierr)
        nreq = 0
       endif
#endif
       call updt_vec_bnd_2(is, ie, js, je, ks, ke, x)
       do j = js+1, je-1
        do i = is+1, ie-1
         x(1,i,j,ks) = -gp(i,j,ks)
        enddo
       enddo
      endif ! igcall
#ifdef MPI_USED
      if(nreq .ne. 0) then
       call mpi_waitall(nreq, req, stat, ierr)
       nreq = 0
      endif
#endif
c
c     -------------------------------------------------------
c     solve the system
c     -------------------------------------------------------
c
      prev_tot = totcgit
      maxitr   = (ie-is+1)  *  (je-js+1)  *  (ke-ks+1)  *
     .           ntiles(1)  *  ntiles(2)  *  ntiles(3)
      toler    = 1.0D-8 ! *
C     .           (ie-is+1)  *  (je-js+1)  *  (ke-ks+1)  *
C     .           ntiles(1)  *  ntiles(2)  *  ntiles(3)  
C#ifdef MPI_USED
C      call MPI_ALLREDUCE(toler,glbtoler,1,MPI_FLOAT,
C     .                   mpi_min,comm3d,ierr)
C#else
      glbtoler = toler
C#endif /* MPI_USED */
c
      call cgsolve(ad, a0, a1, a2, x, b,
     .             glbtoler, error)
      if(nits .ge. maxitr .and. error .gt. glbtoler
     .       .and. myid_w .eq. 0) then
       write(*,"('**********  CGSOLVE did not converge with dt=',
     .  1pe12.5,'  **********',/1x,'(nhy,error,toler,time,nits) =',
     .  i7,3e14.5,i5)")dt,nhy,error,glbtoler,time,nits
#ifdef MPI_USED
       call mpi_finalize(ierr)
#endif
       stop
      else
       totcgit = totcgit + nits
       ncgcall = ncgcall + 1
      endif
c
c     -------------------------------------------------------
c     Finish up
c     -------------------------------------------------------
c
c     positive potential, consistent with ZEUS convention
c
      do k = ks, ke
       do j = js, je
        do i = is, ie
         gp(i,j,k) = -x(1,i,j,k)
        enddo
       enddo
      enddo
c
      nreq = 0
      nsub = nsub + 1
      call updt_gp_bnd_1(is,ie,js,je,ks,ke,gp)
c
c     scale BC arrays to CGS
c
      do k = ks, ke
       do j = js, je
        gpiib(j,k,1) = guniv*gpiib(j,k,1)
        gpiib(j,k,2) =       gpiib(j,k,1)
        gpoib(j,k,1) = guniv*gpoib(j,k,1)
        gpoib(j,k,2) =       gpoib(j,k,1)
       enddo
      enddo
#ifdef MPI_USED
      if(nreq .ne. 0) then
       call mpi_waitall(nreq, req, stat, ierr )
       nreq = 0
      endif
#endif
      nsub = nsub + 1
      call updt_gp_bnd_2(is,ie,js,je,ks,ke,gp)
c
c     scale BC arrays to CGS
c
      do k = ks, ke
       do i = is, ie
        gpijb(i,k,1) = guniv*gpijb(i,k,1)
        gpijb(i,k,2) =       gpijb(i,k,1)
        gpojb(i,k,1) = guniv*gpojb(i,k,1)
        gpojb(i,k,2) =       gpojb(i,k,1)
       enddo
      enddo
#ifdef MPI_USED
      if(nreq .ne. 0) then
       call mpi_waitall(nreq, req, stat, ierr )
       nreq = 0
      endif
#endif
      if(coords(1) .eq. 0) then
       do k = ks, ke
        do j = js, je
         if(niis(3) .eq. 1) gp(is-1,j,k) = gp(is,j,k)
         if(niis(3) .eq. 2) gp(is-1,j,k) = gp(is,j,k)
         if(niis(3) .eq. 3) gp(is-1,j,k) = gpiib(j,k,1)
        enddo
       enddo
      endif ! coords(1)
      if(coords(1) .eq. ntiles(1)-1) then
       do k = ks, ke
        do j = js, je
         if(nois(3) .eq. 1) gp(ie+1,j,k) = gp(ie,j,k)
         if(nois(3) .eq. 2) gp(ie+1,j,k) = gp(ie,j,k)
         if(nois(3) .eq. 3) gp(ie+1,j,k) = gpoib(j,k,1)
        enddo
       enddo
      endif ! coords(1)
      if(coords(2) .eq. 0) then
       do k = ks, ke
        do i = is, ie
         if(nijs(3) .eq. 1) gp(i,js-1,k) = gp(i,js,k)
         if(nijs(3) .eq. 2) gp(i,js-1,k) = gp(i,js,k)
         if(nijs(3) .eq. 3) gp(i,js-1,k) = gpijb(i,k,1)
        enddo
       enddo
      endif
      if(coords(2) .eq. ntiles(2)-1) then
       do k = ks, ke
        do i = is, ie
         if(nojs(3) .eq. 1) gp(i,je+1,k) = gp(i,je,k)
         if(nojs(3) .eq. 2) gp(i,je+1,k) = gp(i,je,k)
         if(nojs(3) .eq. 3) gp(i,je+1,k) = gpojb(i,k,1)
        enddo
       enddo
      endif
c
      igcall = 1
c
      return
      end
c
      subroutine compute_diagonals(ad,a0,a1,is,ie,js,je,ks,ke)
c
      use real_prec
      use param, ONLY: in,jn,kn
      use grid,  ONLY: g2a, g31a, g2bi, g32a,
     .                 dx1bi, dx2bi, dvl1a, dvl2a
c
      implicit NONE
c
      integer  :: i, j, k, is, ie, js, je, ks, ke
c
      real(rl) :: ad(1,1,in,jn,kn),
     .            a0(1,1,in,jn,kn),
     .            a1(1,1,in,jn,kn)
c
      do k = ks, ke
       do j = js, je
        do i = is, ie
         a0(1,1,i,j,k) = dvl2a(j)*g2a(i+1)*g31a(i+1)
     .                 * dx1bi(i+1)
         a1(1,1,i,j,k) = dvl1a(i)*g2bi(i)**2*g32a(j+1)
     .                 * dx2bi(j+1)
         ad(1,1,i,j,k) =
     .           -dvl2a(j)*(g2a(i+1)*g31a(i+1)*dx1bi(i+1)+
     .                               g2a(i  )*g31a(i  )*dx1bi(i  ))
     .           -dvl1a(i)*g2bi(i)**2*
     .                      (g32a(j)*dx2bi(j)+g32a(j+1)*dx2bi(j+1))
        enddo
       enddo
      enddo
c
      return
      end
