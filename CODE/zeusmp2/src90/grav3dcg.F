c=======================================================================
c
c                            Developed by
c                Laboratory for Computational Astrophysics
c                  University of California at San Diego
c
      subroutine grav3D_CG
c
c     3-D gravitational potential solver using John Hayes's
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
      real(rl), dimension(:,:), allocatable ::
     .          gpiibc, gpoibc, gpijbc, gpojbc, gpikbc, gpokbc
c
      real(rl) :: toler, error, glbtoler, eps, glbeps, epstol, epsmax,
     .            epsmin, glbemin, phip(in,jn,kn)
c
      integer :: i, j, k, prev_tot, imax, jmax, kmax
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


      endif ! igcall = 0
c
c-----------------------------------------------------------------------
c     compute BC arrays
c
c     NOTE: see comments in grav2dcg.F regarding use of subroutine
c           gpbv vs BC's specified through f**s(13) arrays.
c-----------------------------------------------------------------------
c
      call gpbv
c
c-----------------------------------------------------------------------
c     Scale Gravitational Constant out of BC arrays
c
c     NOTE: use only if gpbv is NOT called.  fiis(13), fois(13), etc.
c           must be in CGS units in zmp_inp.
c-----------------------------------------------------------------------
c
C     do k = ks, ke
C      do j = js, je
C       gpiib(j,k,1) = fiis(13) / guniv
C       gpiib(j,k,2) = gpiib(j,k,1)
C       gpoib(j,k,1) = fois(13) / guniv
C       gpoib(j,k,2) = gpoib(j,k,1)
C      enddo
C     enddo
C     do k = ks, ke
C      do i = is, ie
C       gpijb(i,k,1) = fijs(13) / guniv
C       gpijb(i,k,2) = gpijb(i,k,1)
C       gpojb(i,k,1) = fojs(13) / guniv
C       gpojb(i,k,2) = gpojb(i,k,1)
C      enddo
C     enddo
C     do j = js, je
C      do i = is, ie
C       gpikb(i,j,1) = fiks(13) / guniv
C       gpikb(i,j,2) = gpikb(i,j,1)
C       gpokb(i,j,1) = foks(13) / guniv
C       gpokb(i,j,2) = gpokb(i,j,1)
C      enddo
C     enddo
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
c     Compute diagonals along 1-faces
c-----------------------------------------------------------------------
c
      call compute_3d_diagonals(is, is, js, je, ks, ke,
     .                          ad, a0, a1, a2)
      call compute_3d_diagonals(ie, ie, js, je, ks, ke,
     .                          ad, a0, a1, a2)
#ifdef MPI_USED
c
c-----------------------------------------------------------------------
c     initiate MPI send/recv's for 1-face data
c-----------------------------------------------------------------------
c
      nreq = 0
      nsub = 0
      call updt_mtrx_bnd_1(is, ie, js, je, ks, ke, a0)
      nsub = nsub + 1
      call updt_mtrx_bnd_1(is, ie, js, je, ks, ke, a1)
      nsub = nsub + 1
      call updt_mtrx_bnd_1(is, ie, js, je, ks, ke, a2)
      nsub = nsub + 1
      call updt_mtrx_bnd_1(is, ie, js, je, ks, ke, ad)
#endif /* MPI_USED */
c
c-----------------------------------------------------------------------
c     in the meantime, compute 2-face data
c-----------------------------------------------------------------------
c
      call compute_3d_diagonals(is+1, ie-1, js, js, ks, ke,
     .                          ad, a0, a1, a2)
      call compute_3d_diagonals(is+1, ie-1, je, je, ks, ke,
     .                          ad, a0, a1, a2)
#ifdef MPI_USED
c
c-----------------------------------------------------------------------
c     verify that 1-face data send/recv's have finished...
c-----------------------------------------------------------------------
c
      if(nreq .ne. 0) then 
       call mpi_waitall(nreq, req, stat, ierr)
       nreq = 0
       nsub = 0
      endif
c
c-----------------------------------------------------------------------
c     and fire off the 2-face data
c-----------------------------------------------------------------------
c
      call updt_mtrx_bnd_2(is, ie, js, je, ks, ke, a0)
      nsub = nsub + 1
      call updt_mtrx_bnd_2(is, ie, js, je, ks, ke, a1)
      nsub = nsub + 1
      call updt_mtrx_bnd_2(is, ie, js, je, ks, ke, a2)
      nsub = nsub + 1
      call updt_mtrx_bnd_2(is, ie, js, je, ks, ke, ad)
#endif /* MPI_USED */
c
c-----------------------------------------------------------------------
c     compute 3-face data
c-----------------------------------------------------------------------
c
      call compute_3d_diagonals(is+1, ie-1, js+1, je-1, ks, ks,
     .                          ad, a0, a1, a2)
      call compute_3d_diagonals(is+1, ie-1, js+1, je-1, ke, ke,
     .                          ad, a0, a1, a2)
#ifdef MPI_USED
c
c-----------------------------------------------------------------------
c     verify that 2-face data has been exchanged...
c-----------------------------------------------------------------------
c
      if(nreq .ne. 0) then 
       call mpi_waitall(nreq, req, stat, ierr)
       nreq = 0
       nsub = 0
      endif
c
c-----------------------------------------------------------------------
c     ...and fire off the 3-face data
c-----------------------------------------------------------------------
c
      call updt_mtrx_bnd_3(is, ie, js, je, ks, ke, a0)
      nsub = nsub + 1
      call updt_mtrx_bnd_3(is, ie, js, je, ks, ke, a1)
      nsub = nsub + 1
      call updt_mtrx_bnd_3(is, ie, js, je, ks, ke, a2)
      nsub = nsub + 1
      call updt_mtrx_bnd_3(is, ie, js, je, ks, ke, ad)
#endif /* MPI_USED */
c
c-----------------------------------------------------------------------
c     compute the matrix elements in all zones interior to the domain
c     faces
c-----------------------------------------------------------------------
c
      call compute_3d_diagonals(is+1, ie-1, js+1, je-1, ks+1, ke-1,
     .                          ad, a0, a1, a2)
#ifdef MPI_USED
c
c-----------------------------------------------------------------------
c     verify that 3-face data exchanges are finished
c-----------------------------------------------------------------------
c
      if(nreq .ne. 0) then 
       call mpi_waitall(nreq, req, stat, ierr)
       nreq = 0
      endif
#endif /* MPI_USED */
c
c-----------------------------------------------------------------------
c     compute RHS
c-----------------------------------------------------------------------
c
      DO k=ks,ke
         DO j=js,je
            DO i=is,ie
               b(1,i,j,k) = 4.0D0*d(i,j,k)*
     .                       PI*guniv*dvl1a(i)*dvl2a(j)*
     .                                dvl3a(k)
            END DO
         END DO
      END DO
c
      epstol = 1.0D-12 ! * sqrt(float((ie-is+1)*(je-js+1)*(ke-ks+1)*
C     .                             ntiles(1)*ntiles(2)*ntiles(3) ) )
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
     .                      dvl2a(j)*dvl3a(k)*g2a(is  )*
     .                      g31a(is  )*dx1bi(is  )
         endif
         if(niis(3) .eq. 3) then 
           b(1,is,j,k) = b(1,is,j,k) + guniv*
     .     dvl2a(j)*dvl3a(k)*g2a(is  )*g31a(is)*dx1bi(is)*gpiib(j,k,1)
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
     .                      dvl2a(j)*dvl3a(k)*g2a(ie+1)*
     .                      g31a(ie+1)*dx1bi(ie+1)
         endif
         if(nois(3) .eq. 3) then
           b(1,ie,j,k) = b(1,ie,j,k) + guniv*
     .     dvl2a(j)*dvl3a(k)*g2a(ie+1)*g31a(ie+1)*dx1bi(ie+1)*
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
     .                      dvl1a(i)*dvl3a(k)*
     .                      g2bi(i)**2*g32a(js  )*dx2bi(js  )
         endif ! nijs
         if(nijs(3) .eq. 3) then
           b(1,i,js,k) = b(1,i,js,k) + guniv*
     .     dvl1a(i)*dvl3a(k)*g2bi(i)**2*g32a(js  )*dx2bi(js  )*
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
     .                      dvl1a(i)*dvl3a(k)*
     .                      g2bi(i)**2*g32a(je+1)*dx2bi(je+1)
         endif ! nojs
         if(nojs(3) .eq. 3) then
           b(1,i,je,k) = b(1,i,je,k) + guniv*
     .     dvl1a(i)*dvl3a(k)*g2bi(i)**2*g32a(je+1)*dx2bi(je+1)*
     .     gpojb(i,k,1)
         endif ! nojs
        enddo ! i
       enddo ! k
      endif ! coords(2)
      if(coords(3) .eq. 0) then
       do j = js, je
        do i = is, ie
        if((niks(3) .eq. 1) .or. (niks(3) .eq. 2)) then
          ad(1,1,i,j,ks) = ad(1,1,i,j,ks) + 
     .                     dvl1a(i)*dvl2a(j)*
     .                     g31bi(i)**2*g32bi(j)**2*dx3bi(ks  )
         endif
         if(niks(3) .eq. 3) b(1,i,j,ks) = b(1,i,j,ks) + guniv*
     .    dvl1a(i)*dvl2a(j)*g31bi(i)**2*g32bi(j)**2*dx3bi(ks  )*
     .    gpikb(i,j,1)
        enddo ! i
       enddo ! j
      endif ! coords(3)
      if(coords(3) .eq. ntiles(3)-1) then
       do j = js, je
        do i = is, ie
         if(noks(3) .ne. 4) a2(1,1,i,j,ke) = 0.0D0
         if((noks(3) .eq. 1) .or. (noks(3) .eq. 2)) then
          ad(1,1,i,j,ke) = ad(1,1,i,j,ke) + 
     .                     dvl1a(i)*dvl2a(j)*
     .                     g31bi(i)**2*g32bi(j)**2*dx3bi(ke+1)
         endif
         if(noks(3) .eq. 3) b(1,i,j,ke) = b(1,i,j,ke) + guniv*
     .    dvl1a(i)*dvl2a(j)*g31bi(i)**2*g32bi(j)**2*dx3bi(ke+1)*
     .    gpokb(i,j,1)
        enddo ! i
       enddo ! k
      endif ! coords(3)
c
      call phicheck(ad, a0, a1, a2, b, gp, phip, epsmax, epsmin,
     .              imax, jmax, kmax)
c
#ifdef MPI_USED
      call mpi_allreduce(epsmax, glbeps, 1, MPI_FLOAT,
     .                   mpi_max, comm3d, ierr)
      call mpi_allreduce(epsmin, glbemin, 1, MPI_FLOAT,
     .                   mpi_min, comm3d, ierr)
c
      epsmax = glbeps
      epsmin = glbemin
#endif
c
C      if(myid .eq. 0) then
C       write(*,"('PHICHECK: max errphi = ',1pd12.4)")epsmax
C       write(*,"('imax, jmax, kmax = ',3i4)")imax,jmax,kmax
C       write(*,"('phip, gp = ',1p2d16.8)")phip(imax,jmax,kmax),
C     .                                    gp  (imax,jmax,kmax)
C      endif
c
      if(epsmax .lt. epstol) then
       do k = ks, ke
        do j = js, je
         gpiib(j,k,1) = guniv*gpiib(j,k,1)
         gpiib(j,k,2) =       gpiib(j,k,1)
         gpoib(j,k,1) = guniv*gpoib(j,k,1)
         gpoib(j,k,2) =       gpoib(j,k,1)
        enddo
       enddo
       do k = ks, ke
        do j = js, je
         gpiib(j,k,1) = guniv*gpiib(j,k,1)
         gpiib(j,k,2) =       gpiib(j,k,1)
         gpoib(j,k,1) = guniv*gpoib(j,k,1)
         gpoib(j,k,2) =       gpoib(j,k,1)
        enddo
       enddo
       do k = ks, ke
        do j = js, je
         gpiib(j,k,1) = guniv*gpiib(j,k,1)
         gpiib(j,k,2) =       gpiib(j,k,1)
         gpoib(j,k,1) = guniv*gpoib(j,k,1)
         gpoib(j,k,2) =       gpoib(j,k,1)
        enddo
       enddo
       return
      else
       ncgcall = ncgcall + 1
      endif
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
#ifdef MPI_USED
       nreq = 0
       nsub = 0
       call updt_vec_bnd_1(is, ie, js, je, ks, ke, x)
#endif
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
       call updt_vec_bnd_2(is, ie, js, je, ks, ke, x)
#endif
       do j = js+1, je-1
        do i = is+1, ie-1
         x(1,i,j,ks) = -gp(i,j,ks)
         x(1,i,j,ke) = -gp(i,j,ke)
        enddo
       enddo
#ifdef MPI_USED
       if(nreq .ne. 0) then
        call mpi_waitall(nreq, req, stat, ierr)
        nreq = 0
       endif
       call updt_vec_bnd_3(is, ie, js, je, ks, ke, x)
#endif
       do k = ks+1, ke-1
        do j = js+1, je-1
         do i = is+1, ie-1
          x(1,i,j,k) = -gp(i,j,k)
         enddo
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
      igcall   = 1
c
      prev_tot = totcgit
      maxitr   = (ie-is+1)  *  (je-js+1)  *  (ke-ks+1)  *
     .           ntiles(1)  *  ntiles(2)  *  ntiles(3)
      toler    = 1.0D-6
CJH      toler    = 1.0D-8 * sqrt(float(
CJH     .           (ie-is+1)*(je-js+1)*(ke-ks+1)*
CJH     .           ntiles(1)*ntiles(2)*ntiles(3) ) )
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
C       if(myid .eq. 0) write(*,"('GRAVITY: cgsolve converged',
C     .    ' in ',i5,' iterations')")nits
       totcgit = totcgit + nits
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
      nsub = nsub + 1
      call updt_gp_bnd_3(is,ie,js,je,ks,ke,gp)
c
c     scale BC arrays to CGS
c
      do j = js, je
       do i = is, ie
        gpikb(i,j,1) = guniv*gpikb(i,j,1)
        gpikb(i,j,2) =       gpikb(i,j,1)
        gpokb(i,j,1) = guniv*gpokb(i,j,1)
        gpokb(i,j,2) =       gpokb(i,j,1)
       enddo
      enddo
c
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
      if(coords(3) .eq. 0) then
       do j = js, je
        do i = is, ie
         if(niks(3) .eq. 1) gp(i,j,ks-1) = gp(i,j,ks)
         if(niks(3) .eq. 2) gp(i,j,ks-1) = gp(i,j,ks)
         if(niks(3) .eq. 3) gp(i,j,ks-1) = gpikb(i,j,1)
        enddo
       enddo
      endif
      if(coords(3) .eq. ntiles(3)-1) then
       do j = js, je
        do i = is, ie
         if(noks(3) .eq. 1) gp(i,j,ke+1) = gp(i,j,ke)
         if(noks(3) .eq. 2) gp(i,j,ke+1) = gp(i,j,ke)
         if(noks(3) .eq. 3) gp(i,j,ke+1) = gpokb(i,j,1)
        enddo
       enddo
      endif
c
C      if(myid .eq. 0)
C     .write(*,"('Phi(i): ',i3,1pd16.8)")(i,gp(i,je,ks),i=is,ie)
c
      return
      end
c
      subroutine compute_3d_diagonals(ibeg,iend,jbeg,jend,kbeg,kend,
     .                                ad, a0, a1, a2)
c
      use real_prec
      use param
      use grid
c
      implicit NONE
c
      integer  :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend
c
      real(rl) :: ad(1,1,in,jn,kn), a0(1,1,in,jn,kn),
     .            a1(1,1,in,jn,kn), a2(1,1,in,jn,kn)
c
      do k = kbeg, kend
       do j = jbeg, jend
        do i = ibeg, iend
CJH         a0(1,1,i,j,k) = dvl2a(j)*dvl3a(k)*g2a(i+1)*g31a(i+1)
CJH     .                 * dx1bi(i+1)
CJH         a1(1,1,i,j,k) = dvl1a(i)*dvl3a(k)*g2bi(i)**2*g32a(j+1)
CJH     .                 * dx2bi(j+1)
         a0(1,1,i,j,k) = dvl2a(j)*dvl3a(k)*(g2a(i+1)*g31a(i+1))**2
     .                 * dvl1bi(i+1)
         a1(1,1,i,j,k) = dvl1a(i)*dvl3a(k)*(g2bi(i)*g32a(j+1))**2
     .                 * dvl2bi(j+1)
         a2(1,1,i,j,k) = dvl1a(i)*dvl2a(j)*g31bi(i)**2*g32bi(j)**2
     .                 * dx3bi(k+1)
         ad(1,1,i,j,k) =
C     .           -dvl2a(j)*dvl3a(k)*( g2a(i+1)*g31a(i+1)*
C     .                                 dx1bi(i+1)+
C     .                                g2a(i  )*g31a(i  )*
C     .                                 dx1bi(i  ) )
C     .           -dvl1a(i)*dvl3a(k)*g2bi(i)**2*
C     .              (g32a(j)*dx2bi(j) + g32a(j+1)*dx2bi(j+1))
     .           -dvl2a(j)*dvl3a(k)*( (g2a(i+1)*g31a(i+1))**2*
     .                                 dvl1bi(i+1)+
     .                                (g2a(i  )*g31a(i  ))**2*
     .                                 dvl1bi(i  ) )
     .           -dvl1a(i)*dvl3a(k)*g2bi(i)**2*
     .            (g32a(j)**2*dvl2bi(j) + g32a(j+1)**2*dvl2bi(j+1))
     .           -dvl1a(i)*dvl2a(j)*g31bi(i)**2*g32bi(j)**2*
     .                              (dx3bi(k) + dx3bi(k+1))
        enddo
       enddo
      enddo
c
      return
      end
c
      subroutine phicheck(ad, a0, a1, a2, b, gp, phip, epsmax,
     .                    epsmin, imax, jmax, kmax)
c
      use real_prec
      use param
      use grid
c
      implicit NONE
c
      integer  :: i, j, k, imax, jmax, kmax
c
      real(rl) :: ad(1,1,in,jn,kn), a0(1,1,in,jn,kn),
     .            a1(1,1,in,jn,kn), a2(1,1,in,jn,kn),
     .            b (1,  in,jn,kn)
c
      real(rl) :: gp(in,jn,kn), phip(in,jn,kn)
c
      real(rl) :: eps, epsmin, epsmax, glbeps, glbemin
c
      do k = ks, ke
       do j = js, je
        do i = is, ie
         phip(i,j,k) = -(1.0/ad(1,1,i,j,k))*
     .                 ( b (1,  i  ,j  ,k  )                 +
     .                   a2(1,1,i  ,j  ,k-1)*gp(i  ,j  ,k-1) +
     .                   a2(1,1,i  ,j  ,k  )*gp(i  ,j  ,k+1) +
     .                   a1(1,1,i  ,j-1,k  )*gp(i  ,j-1,k  ) +
     .                   a1(1,1,i  ,j  ,k  )*gp(i  ,j+1,k  ) +
     .                   a0(1,1,i-1,j  ,k  )*gp(i-1,j  ,k  ) +
     .                   a0(1,1,i  ,j  ,k  )*gp(i+1,j  ,k  )
     .                  )
        enddo
       enddo
      enddo
c
      epsmax = -1.0d30
      epsmin =  1.0d30
      do k = ks, ke
       do j = js, je
        do i = is, ie
         eps = abs(phip(i,j,k)-gp(i,j,k))/(gp(i,j,k)+tiny)
         if(eps .gt. epsmax) then
          epsmax = eps
          imax   = i
          jmax   = j
          kmax   = k
         endif
         if(eps .lt. epsmin) then
          epsmin = eps
         endif
        enddo
       enddo
      enddo
c
      return
      end
