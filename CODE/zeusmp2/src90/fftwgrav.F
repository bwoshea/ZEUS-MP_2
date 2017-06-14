c=======================================================================
c
c                            Developed by
c                Laboratory for Computational Astrophysics
c                  University of California at San Diego
c
c     author : P.S. Li
c
      subroutine fftwgrav
#ifdef FFT
c
c     3-D gravitational potential solver using FFTW package
c     (http://www.fftw.org) to solve Poisson's equation.
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
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer :: i, j, k, nx, ny, nz, iindex, jindex, kindex,
     .           snz, n, kt, ii, jj, source, target, tag, npt,
     .           target1, source1, sindex, remainder, ntotal
C#ifdef MPI_USED
C      integer :: status(MPI_STATUS_SIZE)
C#endif /* MPI_USED */
      real(rl), allocatable :: temp(:,:,:), dslice(:,:,:),
     .                         temp1(:,:,:)
      real(rl) :: h2
c
#ifdef ARCH_IBM
#define POISSON_SOLVER poisson_solver_
#else
#define POISSON_SOLVER poisson_solver
#endif /* ARCH_IBM */
C
C Initialize data
C
      nx=ie-is+1
      ny=je-js+1
      nz=ke-ks+1
      allocate (temp(nx*ntiles(1),ny*ntiles(2),nz*ntiles(3)/nprocs_w))
      allocate (temp1(in,jn,kn))
      allocate (dslice(in,jn,kn))
      ntotal=in*jn*kn
      h2=dx1a(is)*dx1a(is)*2.D0*pi*guniv
CPS
C  Block domain -> slab domain
C
#ifdef MPI_USED
      kt=ntiles(3)*nz/nprocs_w
C
C  Send self density information to others and receive others density 
C  information on the same tile layer in 3-direction.  The naming
C  convention of myid and coords sets k index runs first and i index last.
C
      nsub = nsub + 1
      nreq = 0
      do i=0,ntiles(1)-1
        do j=0,ntiles(2)-1
          if(coords(1).eq.i.AND.coords(2).eq.j) then
            do ii=0,ntiles(1)-1
              do jj=0,ntiles(2)-1
                if(coords(1).ne.ii.OR.coords(2).ne.jj) then
                  source=ii*ntiles(2)*ntiles(3)+jj*ntiles(3)+coords(3)
                  kindex=(source/ntiles(3))*kt+ks-1
                  do k=1,kt
                    nreq = nreq + 1
                    call MPI_IRECV(dslice(1,1,k+kindex),1,k_slice,
     &                source,41000+source*kt+k+nsub,comm3d,
     &                req(nreq),ierr)
                  enddo
                endif
              enddo
            enddo
          else
            target=i*ntiles(2)*ntiles(3)+j*ntiles(3)+coords(3)
            kindex=(target/ntiles(3))*kt+ks-1
            do k=1,kt
              nreq = nreq + 1
              call MPI_ISEND(d(1,1,k+kindex),1,k_slice,target,
     &          41000+myid_w*kt+k+nsub,comm3d,req(nreq),ierr)
            enddo
          endif
        enddo
      enddo
C
C  Wait for communications to complete.
C
      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
C
C  Copy self density to self temporary array.
C
      do k=1,kt
        kindex=(myid_w/ntiles(3))*kt+ks-1
        dslice(:,:,k+kindex)=d(:,:,k+kindex)
      enddo
C
C  Resequence processors in 3-direction for FFTW
C
      nreq=0
      npt=nprocs_w/ntiles(3)
      sindex=ntiles(1)*ntiles(2)-1
      remainder=mod(myid_w,ntiles(3))
      source1=(myid_w-(myid_w/npt)*sindex)*ntiles(3)+
     &        (myid_w/npt)*(1-ntiles(3))
      target1=myid_w+(myid_w-remainder)*(1-ntiles(3))/ntiles(3)+
     &        sindex*remainder
      if(source1.NE.myid_w) then
        nreq = nreq + 1
        call MPI_IRECV(temp1,ntotal,MPI_FLOAT,source1,
     &                 46000+source1+nsub,comm3d,req(nreq),ierr)
        nreq = nreq + 1
        call MPI_ISEND(dslice,ntotal,MPI_FLOAT,target1,
     &                 46000+myid_w+nsub,comm3d,req(nreq),ierr)
      else
        temp1=dslice
      endif
C
C  Wait for communications to complete.
C
      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
C
C  Block to slab and corp the boundary regions for FFTW.
C
      do i=0,ntiles(1)-1
        do j=0,ntiles(2)-1
          kindex=(i*ntiles(2)+j)*kt+ks-1
          do k=1,kt
            do jj=1,ny
              do ii=1,nx
                temp(i*nx+ii,j*ny+jj,k)=temp1(ii+is-1,jj+js-1,
     &            k+kindex)*h2
              enddo
            enddo
          enddo
        enddo
      enddo
#endif /* MPI_USED */
C
C   Call the Poisson Solver subroutine
C
      call POISSON_SOLVER(myid_w,nprocs_w,ntiles(1)*nx,
     &                    ntiles(2)*ny,kt,temp)

#ifdef MPI_USED
CPS
C  Slab domain -> Block domain
C
C  Slab to block and attach boundary regions back.
C
      do i=0,ntiles(1)-1
        do j=0,ntiles(2)-1
          kindex=(i*ntiles(2)+j)*kt+ks-1
          do k=1,kt
            do jj=1,ny
              do ii=1,nx
                dslice(ii+is-1,jj+js-1,k+kindex)=
     &            temp(i*nx+ii,j*ny+jj,k)
              enddo
            enddo
          enddo
        enddo
      enddo
C
C  Resequence processors in 3-direction for block domain restoration.
C
      nreq=0
      if(source1.NE.myid_w) then
        nreq = nreq + 1
        call MPI_IRECV(temp1,ntotal,MPI_FLOAT,target1,
     &                 26000+target1+nsub,comm3d,req(nreq),ierr)
        nreq = nreq + 1
        call MPI_ISEND(dslice,ntotal,MPI_FLOAT,source1,
     &                 26000+myid_w+nsub,comm3d,req(nreq),ierr)
      else
        temp1=dslice
      endif
C
C  Wait for communications to complete.
C
      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
C
C  Final restoration of potential in the same tile layer in 3-direction.
C
      do k=1,kt
        kindex=(myid_w/ntiles(3))*kt+ks-1
        gp(is:ie,js:je,k+kindex)=temp1(is:ie,js:je,k+kindex)
      enddo

      nreq = 0
      do i=0,ntiles(1)-1
        do j=0,ntiles(2)-1
          if(coords(1).eq.i.AND.coords(2).eq.j) then
            do ii=0,ntiles(1)-1
              do jj=0,ntiles(2)-1
                if(coords(1).ne.ii.OR.coords(2).ne.jj) then
                  source=ii*ntiles(2)*ntiles(3)+jj*ntiles(3)+coords(3)
                  kindex=(source/ntiles(3))*kt+ks-1
                  do k=1,kt
                    nreq = nreq + 1
                    call MPI_IRECV(gp(1,1,k+kindex),1,k_slice,
     &                source,21000+source*kt+k+nsub,comm3d,
     &                req(nreq),ierr)
                  enddo
                endif
              enddo
            enddo
          else
            target=i*ntiles(2)*ntiles(3)+j*ntiles(3)+coords(3)
            kindex=(target/ntiles(3))*kt+ks-1
            do k=1,kt
              nreq = nreq + 1
              call MPI_ISEND(temp1(1,1,k+kindex),1,k_slice,
     &          target,21000+myid_w*kt+k+nsub,comm3d,
     &          req(nreq),ierr)
            enddo
          endif
        enddo
      enddo
C
C  Wait for communications to complete.
C
      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )

#endif /* MPI_USED */
C
C Update boundary of gravitational potential.
C
#ifdef MPI_USED
C      nsub = 1
      nreq = 0
      nreq = nreq + 1
      call MPI_IRECV(gp(is-1,1,1),1,i_slice,n1m
     &              ,31000+nsub, comm3d, req(nreq), ierr)
      nreq = nreq + 1
      call MPI_IRECV(gp(ie+1,1,1),1,i_slice,n1p
     &              ,31100+nsub, comm3d, req(nreq), ierr)
      nreq = nreq + 1
      call MPI_ISEND(gp(is,1,1),1,i_slice,n1m
     &              ,31100+nsub, comm3d, req(nreq), ierr)
      nreq = nreq + 1
      call MPI_ISEND(gp(ie,1,1),1,i_slice,n1p
     &              ,31000+nsub, comm3d, req(nreq), ierr)

      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )

      nreq = 0
      nreq = nreq + 1
      call MPI_IRECV(gp(1,js-1,1),1,j_slice,n2m
     &              ,31200+nsub, comm3d, req(nreq), ierr)
      nreq = nreq + 1
      call MPI_IRECV(gp(1,je+1,1),1,j_slice,n2p
     &              ,31300+nsub, comm3d, req(nreq), ierr)
      nreq = nreq + 1
      call MPI_ISEND(gp(1,js,1),1,j_slice,n2m
     &              ,31300+nsub, comm3d, req(nreq), ierr)
      nreq = nreq + 1
      call MPI_ISEND(gp(1,je,1),1,j_slice,n2p
     &              ,31200+nsub, comm3d, req(nreq), ierr)

      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
      
      nreq = 0
      nreq = nreq + 1
      call MPI_IRECV(gp(1,1,ks-1),1,k_slice,n3m
     &              ,31400+nsub, comm3d, req(nreq), ierr)
      nreq = nreq + 1
      call MPI_IRECV(gp(1,1,ke+1),1,k_slice,n3p
     &              ,31500+nsub, comm3d, req(nreq), ierr)
      nreq = nreq + 1
      call MPI_ISEND(gp(1,1,ks),1,k_slice,n3m
     &              ,31500+nsub, comm3d, req(nreq), ierr)
      nreq = nreq + 1
      call MPI_ISEND(gp(1,1,ke),1,k_slice,n3p
     &              ,31400+nsub, comm3d, req(nreq), ierr)

      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
C
C Cleaning
C
      deallocate(temp)
      deallocate(temp1)
      deallocate(dslice)
C
#endif /* MPI_USED */
C
#endif /* FFT */
      return
      end
