c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              G R E Y _ F L D              \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c=======================================================================
c
      subroutine grey_fld(dcopy, eold, eod , erold, erod)
c
c     call       grey_fld(w3dd , e   , w3de, w3dh , er)
c
c
c     driver routine for computing an implicit solution to the
c     radiation diffusion equation.  The structure of
c     the N-R iteration and convergence control is taken from ZEUS-2D.
c
c     written by: John C. Hayes; June, 1997
c     re-written for F90 by J. Hayes; May, 2003
c
      use real_prec
      use config
      use param
      use field
      use bndry
      use grid
      use root
      use radiation
      use cons
      use opac
      use opac_law
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif /* MPI_USED */
      use mpipar
      use clockmod
      use impsoln
c
      implicit NONE
c
c     matrix arrays
c
      real(rl) :: dd  (neqm,neqm,in,jn,kn)
      real(rl) :: ddp1(neqm,neqm,in,jn,kn)
      real(rl) :: ddp2(neqm,neqm,in,jn,kn)
      real(rl) :: ddp3(neqm,neqm,in,jn,kn)
      real(rl) :: rhs (neqm,     in,jn,kn)
      real(rl) :: x   (neqm,     in,jn,kn)
c
c     work arrays
c
      real(rl) ::   fn1     (in,jn,kn)
      real(rl) ::   fn2     (in,jn,kn)
      real(rl) :: d_fn1_d_er(in,jn,kn)
      real(rl) :: d_fn1_d_eg(in,jn,kn)
      real(rl) :: d_fn2_d_er(in,jn,kn)
      real(rl) :: d_fn2_d_eg(in,jn,kn)
      real(rl) :: grdvcf    (in,jn,kn)
c
c     matrix solver control
c
      real(rl) ::    toler, error
      real(rl) ::    lssave
c
c     N-R iteration control
c
      real(rl) :: dthydro, trad, epsemx, epsermx, q1, demx, dermx,
     .            enorm, ernorm, glbepsemx, glbepsermx,
     .            glbenorm, glbernorm, glbdemx, glbdermx, alpha
c
      real(rl) :: iterhst(21,7), inp(2), erout(2), egout(2)
c
      logical boost
c
c     thermodynamic scratch variables
c
      real(rl) :: dertot (in,jn,kn), dersave(in,jn,kn),
     .     erold  (in,jn,kn), ernew  (in,jn,kn), deldotf(in,jn,kn),
     .     deldotv(in,jn,kn), erod   (in,jn,kn),
     .     t      (in,jn,kn), dtde   (in,jn,kn), dbbdt  (in,jn,kn), 
     .     dcopy  (in,jn,kn), eold   (in,jn,kn), eod    (in,jn,kn)
c
      real(rl) :: dersq, ersq, glbdersq, glbersq
c
c     DO loop control
c
      integer  :: i, j, k, iter, lswp
c
      integer  :: igmax, jgmax, kgmax, irmax, jrmax, krmax
c
      lssave = totlsit
c
      do k = ks, ke
       do j = js, je
        do i = is, ie
         erold(i,j,k) = erold(i,j,k) * d(i,j,k)
        enddo
       enddo
      enddo
c
      do i = 1, 6
       bvstat(i,6) = 0    ! er
      enddo
c
      do k = ks, ke
       do j = js, je
        do i = is, ie
         dersave(i,j,k) = der(i,j,k)
        enddo
       enddo
      enddo
c
      nreq = 0
      nsub = nsub + 1
      call bvald  (1,1,0,0,0,0,    d)
      call bvale  (3,3,0,0,0,0, eold)
      call bvalers(3,3,0,0,0,0,erold)
      if(xhydro) then
       call bvalv1 (1,1,0,0,0,0,   v1) !CJH
       call bvalv2 (1,1,0,0,0,0,   v2)
       call bvalv3 (1,1,0,0,0,0,   v3)
      endif ! xhydro
c
      call matprop(eold, erold, d, gamma, t, dtde, p, bb, dbbdt,
     .             kapr, kap, sig, dkapdt, kem, dkemdt, dpde,
     .             is, ie, js, je, ks, ke)
c
      do k = ks, ke
       do j = js, je
        do i = is, ie
         pn    (i,j,k) = p     (i,j,k)
         bbn   (i,j,k) = bb    (i,j,k)
         kapn  (i,j,k) = kap   (i,j,k)
         kemn  (i,j,k) = kem   (i,j,k)
         dbb de(i,j,k) = dbb dt(i,j,k) * dtde(i,j,k)
         dkapde(i,j,k) = dkapdt(i,j,k) * dtde(i,j,k)
         dkemde(i,j,k) = dkemdt(i,j,k) * dtde(i,j,k)
        enddo
       enddo
      enddo
c
c-----------------------------------------------------------------------
c     wait for communications to complete
c-----------------------------------------------------------------------
c
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif /* MPI_USED */
c
c-----------------------------------------------------------------------
c     if not 1D, then post sends and receives along j boundaries
c-----------------------------------------------------------------------
c
      if(ldimen .gt. 1) then
       nreq = 0
       nsub = nsub + 1
       call bvald  (0,0,1,1,0,0,    d)
       call bvale  (0,0,1,1,0,0, eold)
       call bvalers(0,0,3,3,0,0,erold)
       if(xhydro) then
        call bvalv1 (0,0,1,1,0,0,v1)
        call bvalv2 (0,0,1,1,0,0,v2) !CJH
        call bvalv3 (0,0,1,1,0,0,v3)
       endif ! xhydro
      endif ! ldimen
c
c-----------------------------------------------------------------------
c     finish up calculations at i boundaries
c-----------------------------------------------------------------------
c
      call matprop(eold, erold, d, gamma, t, dtde, p, bb, dbbdt,
     .             kapr, kap, sig, dkapdt, kem, dkemdt, dpde,
     .             is-1, is-1, js, je, ks, ke)
      call matprop(eold, erold, d, gamma, t, dtde, p, bb, dbbdt,
     .             kapr, kap, sig, dkapdt, kem, dkemdt, dpde,
     .             ie+1, ie+1, js, je, ks, ke)
c
      i = is-1
      do k = ks, ke
       do j = js, je
         pn    (i,j,k) = p     (i,j,k)
         bbn   (i,j,k) = bb    (i,j,k)
         kapn  (i,j,k) = kap   (i,j,k)
         kemn  (i,j,k) = kem   (i,j,k)
         dbb de(i,j,k) = dbb dt(i,j,k) * dtde(i,j,k)
         dkapde(i,j,k) = dkapdt(i,j,k) * dtde(i,j,k)
         dkemde(i,j,k) = dkemdt(i,j,k) * dtde(i,j,k)
       enddo
      enddo
c
      i = ie+1
      do k = ks, ke
       do j = js, je
         pn    (i,j,k) = p     (i,j,k)
         bbn   (i,j,k) = bb    (i,j,k)
         kapn  (i,j,k) = kap   (i,j,k)
         kemn  (i,j,k) = kem   (i,j,k)
         dbb de(i,j,k) = dbb dt(i,j,k) * dtde(i,j,k)
         dkapde(i,j,k) = dkapdt(i,j,k) * dtde(i,j,k)
         dkemde(i,j,k) = dkemdt(i,j,k) * dtde(i,j,k)
       enddo
      enddo
c
      if(ldimen .gt. 1) then
c
c-----------------------------------------------------------------------
c     if j-boundary communications pending, wait for them to finish
c-----------------------------------------------------------------------
c
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif /* MPI_USED */
c
c-----------------------------------------------------------------------
c     if 3D, then post sends and receives for the k boundaries
c-----------------------------------------------------------------------
c
       if(ldimen .eq. 3) then
        nreq = 0
        nsub = nsub + 1
        call bvald  (0,0,0,0,1,1,    d)
        call bvale  (0,0,0,0,1,1, eold)
        call bvalers(0,0,0,0,3,3,erold)
        if(xhydro) then
         call bvalv1 (0,0,0,0,1,1,v1)
         call bvalv2 (0,0,0,0,1,1,v2)
         call bvalv3 (0,0,0,0,1,1,v3) !CJH
        endif ! xhydro
       endif ! ldimen = 3
c
c-----------------------------------------------------------------------
c     finish up calculations on j boundaries
c-----------------------------------------------------------------------
c
       call matprop(eold, erold, d, gamma, t, dtde, p, bb, dbbdt,
     .              kapr, kap, sig, dkapdt, kem, dkemdt, dpde,
     .              is, ie, js-1, js-1, ks, ke)
       call matprop(eold, erold, d, gamma, t, dtde, p, bb, dbbdt,
     .              kapr, kap, sig, dkapdt, kem, dkemdt, dpde,
     .              is, ie, je+1, je+1, ks, ke)
c
       j = js-1
       do k = ks, ke
        do i = is, ie
          pn    (i,j,k) = p     (i,j,k)
          bbn   (i,j,k) = bb    (i,j,k)
          kapn  (i,j,k) = kap   (i,j,k)
          kemn  (i,j,k) = kem   (i,j,k)
          dbb de(i,j,k) = dbb dt(i,j,k) * dtde(i,j,k)
          dkapde(i,j,k) = dkapdt(i,j,k) * dtde(i,j,k)
          dkemde(i,j,k) = dkemdt(i,j,k) * dtde(i,j,k)
        enddo
       enddo
c
       j = je+1
       do k = ks, ke
        do i = is, ie
          pn    (i,j,k) = p     (i,j,k)
          bbn   (i,j,k) = bb    (i,j,k)
          kapn  (i,j,k) = kap   (i,j,k)
          kemn  (i,j,k) = kem   (i,j,k)
          dbb de(i,j,k) = dbb dt(i,j,k) * dtde(i,j,k)
          dkapde(i,j,k) = dkapdt(i,j,k) * dtde(i,j,k)
          dkemde(i,j,k) = dkemdt(i,j,k) * dtde(i,j,k)
        enddo
       enddo
c
c-----------------------------------------------------------------------
c     if 3D, then wait for pending communications on k boundaries
c-----------------------------------------------------------------------
c
       if(ldimen .eq. 3) then
#ifdef MPI_USED
        if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif /* MPI_USED */
c
c-----------------------------------------------------------------------
c     finish up calculations on k boundaries
c-----------------------------------------------------------------------
c
        call matprop(eold, erold, d, gamma, t, dtde, p, bb, dbbdt,
     .               kapr, kap, sig, dkapdt, kem, dkemdt, dpde,
     .               is, ie, js, je, ks-1, ks-1)
        call matprop(eold, erold, d, gamma, t, dtde, p, bb, dbbdt,
     .               kapr, kap, sig, dkapdt, kem, dkemdt, dpde,
     .               is, ie, js, je, ke+1, ke+1)
c
        k = ks-1
        do j = js, je
         do i = is, ie
           pn    (i,j,k) = p     (i,j,k)
           bbn   (i,j,k) = bb    (i,j,k)
           kapn  (i,j,k) = kap   (i,j,k)
           kemn  (i,j,k) = kem   (i,j,k)
           dbb de(i,j,k) = dbb dt(i,j,k) * dtde(i,j,k)
           dkapde(i,j,k) = dkapdt(i,j,k) * dtde(i,j,k)
           dkemde(i,j,k) = dkemdt(i,j,k) * dtde(i,j,k)
         enddo
        enddo
c
        k = ke+1
        do j = js, je
         do i = is, ie
           pn    (i,j,k) = p     (i,j,k)
           bbn   (i,j,k) = bb    (i,j,k)
           kapn  (i,j,k) = kap   (i,j,k)
           kemn  (i,j,k) = kem   (i,j,k)
           dbb de(i,j,k) = dbb dt(i,j,k) * dtde(i,j,k)
           dkapde(i,j,k) = dkapdt(i,j,k) * dtde(i,j,k)
           dkemde(i,j,k) = dkemdt(i,j,k) * dtde(i,j,k)
         enddo
        enddo
c
       endif ! ldimen = 3
      endif ! ldimen > 1
c
c======================================================================
c     evaluate grad(v):P term
c======================================================================
c
      call grdv(is, ie, js, je, ks, ke, erold, grdvcf,
     .          deldotv)
c
c-----------------------------------------------------------------------
c     initialize timestep and boost iteration control
c-----------------------------------------------------------------------
c
      nred    = 0
      dthydro = dt
 1    continue
c
      trad    = 0.0
c
 2    continue
c
      boost   = .false.
c
c======================================================================
c     establish initial guess for new radiation and gas energy
c     densities; initialize running sums for cumulative changes
c     in both variables
c======================================================================
c
      do k = ks, ke
       do j = js, je
        do i = is, ie
         ern(i,j,k) = erold(i,j,k) + der(i,j,k)
         e n(i,j,k) = e old(i,j,k) + de (i,j,k)
        enddo
       enddo
      enddo
c
c=======================================================================
c     update old values of radiation energy density in ghost zones;
c     overlap with initialization of guess for radiation energy density
c     at new timestep
c=======================================================================
c
#ifdef MPI_USED
      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif /* MPI_USED */
c
c=======================================================================
c     calculate diffusion coefficients from old values of radiation
c     energy density; overlap with update of new radiation energy
c     density in ghost zones
c=======================================================================
c
      do i = 1, 6
       bvstat(i,6) = 0    ! er
      enddo
      nreq = 0
      nsub = nsub + 1
      call bvalers(3,3,0,0,0,0,ern)
      if(ldimen .gt. 1) then
       nsub = nsub + 1
       call bvalers(0,0,3,3,0,0,ern)
       if(ldimen .eq. 3) then
        nsub = nsub + 1
        call bvalers(0,0,0,0,3,3,ern)
       endif ! ldimen = 3
      endif ! ldimen > 1
c
      call difco1(erold,dr1,is,ie+1,js,je  ,ks,ke  )
      call difco2(erold,dr2,is,ie  ,js,je+1,ks,ke  )
      call difco3(erold,dr3,is,ie  ,js,je  ,ks,ke+1)
c
      do k = 1, kn
       do j = 1, jn
        do i = 1, in
         dertot(i,j,k) = 0.D0
         x   (1,i,j,k) = 0.D0 ! der(i,j,k)
        enddo
       enddo
      enddo
c
#ifdef MPI_USED
      if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
      nreq = 0
#endif /* MPI_USED */
c
c
c======================================================================
c     evaluate sub- and super-diagonal elements of the jacobian;
c     all of which are constant during the N-R iteration. 
c======================================================================
c
      call offdiag(ddp1, ddp2, ddp3)
#ifdef MPI_USED
      nreq = 0
      nsub = nsub + 1
      call msendrec_bnd(is,ie,js,je,ks,ke,ddp1)
      if(ldimen .gt. 1) then
       nsub = nsub + 1
       call msendrec_bnd(is,ie,js,je,ks,ke,ddp2)
       if(ldimen .eq. 3) then
        nsub = nsub + 1
        call msendrec_bnd(is,ie,js,je,ks,ke,ddp3)
       endif ! ldimen = 3
      endif ! > ldimen > 1
#endif /* MPI_USED */
c
c======================================================================
c     Start the N-R iteration
c======================================================================
c
      do 1000 iter = 1, nmeiter
c
c======================================================================
c     complete the jacobian equation with the main diagonal elements
c     and the right-hand side
c======================================================================
c
       call jacob(eold, erold, dd, rhs, fn1, fn2,
     .                  d_fn1_d_er, d_fn1_d_eg,
     .                  d_fn2_d_er,
     .                  d_fn2_d_eg, grdvcf, deldotv, deldotf)
c
#ifdef MPI_USED
       nsub = nsub + 1
       call msendrec_bnd(is,ie,js,je,ks,ke,dd)
       nsub = nsub + 1
       call sendrec_bnd(is,ie,js,je,ks,ke,x)
#endif /* MPI_USED */
C      if(nhy .eq. nlim-1) then
C       if(myid .eq. 0) open(666,file='look.0')
C       if(myid .eq. 1) open(666,file='look.1')
C       if(myid .eq. 2) open(666,file='look.2')
C       if(myid .eq. 3) open(666,file='look.3')
C       if(myid .eq. 4) open(666,file='look.4')
C       if(myid .eq. 5) open(666,file='look.5')
C       if(myid .eq. 6) open(666,file='look.6')
C       if(myid .eq. 7) open(666,file='look.7')
C       j = js
C       k = ks
C       write(666,"(1p4d12.4)")
C     .           (x1b(i),grdvcf(i,j,k),deldotv(i,j,k),rhs(1,i,j,k),
C     .           i=is,ie)
C       call mpi_finalize(ierr)
C       stop
C      endif
c
c======================================================================
c     solve the matrix for the changes to Erad
c======================================================================
c
       maxitr  = (ie-is+1)  *  (je-js+1)  *  (ke-ks+1)  *
     .           ntiles(1)  *  ntiles(2)  *  ntiles(3)
       toler   = epsrad *
     .           (ie-is+1)  *  (je-js+1)  *  (ke-ks+1)  *
     .           ntiles(1)  *  ntiles(2)  *  ntiles(3)
c
       call cgsolve(dd, ddp1,
     .                  ddp2,
     .                  ddp3,
     .              x, rhs, toler, error)
c
       totlsit = totlsit + float(nits)
       if(nits .ge. maxitr .and. error .gt. toler
     .        .and. myid_w .eq. 0) then
        write(*,"('**********  CGSOLVE did not converge with dt=',
     .   1pe12.5,'  **********',/1x,'(nhy,error,toler,time,nits) =',
     .   i7,3e14.5,i5)")dt,nhy,error,toler,time,nits
        write(*,"('i,dd,ddp1,x,rhs = ',i3,1p4d12.4)")(i,
     .   dd(1,1,i,3,3),ddp1(1,1,i,3,3),x(1,i,3,3),
     .   rhs (1,  i,3,3), i=is,ie)
        stop
        goto 1001
       else
       endif
       iterhst(iter,6) = float(nits)
       iterhst(iter,7) = error
c
c======================================================================
c     Map the solution vector from LSSOLVE into the correction
c     array
c======================================================================
c
       do k = ks, ke
        do j = js, je
         do i = is, ie
          der(i,j,k) = x(1,i,j,k)
          de (i,j,k) = 
     .                 -(fn2(i,j,k)+d_fn2_d_er(i,j,k)*der(i,j,k)) /
     .                              d_fn2_d_eg(i,j,k)
         enddo
        enddo
       enddo
c
c======================================================================
c     restrict maximum change in a variable by q1
c     sum changes in er into running totals
c======================================================================
c
       epsermx = 0.D0
       epse mx = 0.D0
       dersq   = 0.D0
       ersq    = 0.D0
c
       do k = ks, ke
        do j = js, je
         do i = is, ie
          dertot(i,j,k) = dertot(i,j,k) + der(i,j,k)
c
c     gas energy convergence criteria
c
          if(cnvcrit .le. 2) then
           if(abs(de(i,j,k)/en(i,j,k)) .gt. epsemx) then
            igmax = i
            jgmax = j
            kgmax = k
           endif
           epsemx  = max(abs(de(i,j,k)/en(i,j,k)),epsemx)
          else
           if(abs(de(i,j,k)/ennom) .gt. epsemx) then
            igmax = i
            jgmax = j
            kgmax = k
           endif
           epsemx  = max(abs(de(i,j,k)/ennom),epsemx)
          endif
c
c     radiation energy convergence criteria
c
          irmax = is
          jrmax = js
          krmax = ks
          if(cnvcrit .eq. 1) then
           if(abs(der(i,j,k)/ern(i,j,k)) .gt. epsermx) then
              irmax      = i
              jrmax      = j
              krmax      = k
           endif
           epsermx       = max(abs(der(i,j,k)/ern(i,j,k)),epsermx)
          else if(cnvcrit .eq. 2) then
           if(abs(der(i,j,k)/ernom) .gt. epsermx) then
            irmax      = i
            jrmax      = j
            krmax      = k
           endif
           epsermx       = max(abs(der(i,j,k)/ernom),epsermx)
          else if(cnvcrit .eq. 3) then
           ersq  =  ersq + ern(i,j,k)**2
           dersq = dersq + der(i,j,k)**2
          endif
         enddo
        enddo
       enddo
C       if(nhy .eq. nlim-1) then
C        write(666,"('er,erm1,kr,krm1,dr1 = ',1p5d12.4)")
C     .   (erold(i,js,ks),erold(i-1,js,ks),kapr(i,js,ks),
C     .    kapr(i-1,js,ks),dr1(i,js,ks),i=is,ie)
C        call mpi_finalize(ierr)
C        stop
C       endif
#ifdef MPI_USED
       call mpi_allreduce(epsemx, glbepsemx, 1,
     .          MPI_FLOAT,mpi_max, comm3d, ierr)
       call mpi_allreduce(epsermx, glbepsermx, 1,
     .          MPI_FLOAT,mpi_max, comm3d, ierr)
       epse mx = glbepse mx
       epsermx = glbepsermx
#endif /* MPI_USED */
       q1 = min(one,min(demax/epsemx,dermax/epsermx))
c
c======================================================================
c     apply corrections to variables
c======================================================================
c
       do k = ks, ke
        do j = js, je
         do i = is, ie
          e n(i,j,k) = e n(i,j,k) + q1*de (i,j,k)
          ern(i,j,k) = ern(i,j,k) + q1*der(i,j,k)
          x(1,i,j,k) = 0.D0
         enddo
        enddo
       enddo
       totnrit = totnrit + 1.D0
c
c======================================================================
c     mark boundaries out of date
c======================================================================
c
       do i = 1, 6
        bvstat(i,2) = 0    ! e
        bvstat(i,6) = 0    ! er
       enddo
c
c======================================================================
c     update boundaries
c======================================================================
c
       nreq = 0
       nsub = nsub + 1
       call bvalers(3,3,0,0,0,0,ern)
       call bvale  (1,1,0,0,0,0,e n)
       if(ldimen .gt. 1) then
        nsub = nsub + 1
        call bvalers(0,0,3,3,0,0,ern)
        call bvale  (0,0,1,1,0,0,e n)
        if(ldimen .eq. 3) then
         nsub = nsub + 1
         call bvalers(0,0,0,0,3,3,ern)
         call bvale  (0,0,0,0,1,1,e n)
        endif ! ldimen = 3
       endif ! ldimen > 1
c
c======================================================================
c     if converged, exit
c======================================================================
c
       if(epsermx .lt. epsme
     .    .and. epsemx .le. epsme
     .                         ) then
#ifdef MPI_USED
        if(nreq .ne. 0) call MPI_WAITALL (nreq,req,stat,ierr)
        nreq = 0
#endif /* MPI_USED */
        if ((boost .eqv. .true.) .or. (epsermx .lt. epsme)
     .                           .or. (epse mx .lt. epsme) ) 
     .         goto 1002
        boost = .true.
       endif
c
c======================================================================
c     otherwise, set up for next iteration
c======================================================================
c
c
       call matprop(en, ern, d, gamma, t, dtde, pn, bbn, dbbdt,
     .              kapr, kapn, sig, dkapdt, kemn, dkemdt, dpde,
     .              is, ie, js, je, ks, ke)
c
       do k = ks, ke
        do j = js, je
         do i = is, ie
          dbb de(i,j,k) = dbb dt(i,j,k) * dtde(i,j,k)
          dkapde(i,j,k) = dkapdt(i,j,k) * dtde(i,j,k)
          dkemde(i,j,k) = dkemdt(i,j,k) * dtde(i,j,k)
         enddo
        enddo
       enddo
#ifdef MPI_USED
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif /* MPI_USED */
       call matprop(en, ern, d, gamma, t, dtde, pn, bbn, dbbdt,
     .              kapr, kapn, sig, dkapdt, kemn, dkemdt, dpde,
     .              is-1, is-1, js, je, ks, ke)
       call matprop(en, ern, d, gamma, t, dtde, pn, bbn, dbbdt,
     .              kapr, kapn, sig, dkapdt, kemn, dkemdt, dpde,
     .              ie+1, ie+1, js, je, ks, 
     .              ke)
       if(ldimen .gt. 1) then
        call matprop(en, ern, d, gamma, t, dtde, pn, bbn, dbbdt,
     .               kapr, kapn, sig, dkapdt, kemn, dkemdt, dpde,
     .               is, ie, js-1, js-1, ks, 
     .               ke)
c
        j = js-1
        do k = ks, ke
         do i = is, ie
          dbb de(i,j,k) = dbb dt(i,j,k) * dtde(i,j,k)
          dkapde(i,j,k) = dkapdt(i,j,k) * dtde(i,j,k)
          dkemde(i,j,k) = dkemdt(i,j,k) * dtde(i,j,k)
         enddo
        enddo
c
        call matprop(en, ern, d, gamma, t, dtde, pn, bbn, dbbdt,
     .               kapr, kapn, sig, dkapdt, kemn, dkemdt, dpde,
     .               is, ie, je+1, je+1, ks, 
     .               ke)
c
        j = je+1
        do k = ks, ke
         do i = is, ie
          dbb de(i,j,k) = dbb dt(i,j,k) * dtde(i,j,k)
          dkapde(i,j,k) = dkapdt(i,j,k) * dtde(i,j,k)
          dkemde(i,j,k) = dkemdt(i,j,k) * dtde(i,j,k)
         enddo
        enddo
        if(ldimen .eq. 3) then
         call matprop(en, ern, d, gamma, t, dtde, pn, bbn, dbbdt,
     .                kapr, kapn, sig, dkapdt, kemn, dkemdt, dpde,
     .                is, ie, js, je, ks-1, 
     .                ks-1)
c
         k = ks-1
         do j = js, je
          do i = is, ie
           dbb de(i,j,k) = dbb dt(i,j,k) * dtde(i,j,k)
           dkapde(i,j,k) = dkapdt(i,j,k) * dtde(i,j,k)
           dkemde(i,j,k) = dkemdt(i,j,k) * dtde(i,j,k)
          enddo
         enddo
c
         call matprop(en, ern, d, gamma, t, dtde, pn, bbn, dbbdt,
     .                kapr, kapn, sig, dkapdt, kemn, dkemdt, dpde,
     .                is, ie, js, js, ke+1, 
     .                ke+1)
c
         k = ke+1
         do j = js, je
          do i = is, ie
           dbb de(i,j,k) = dbb dt(i,j,k) * dtde(i,j,k)
           dkapde(i,j,k) = dkapdt(i,j,k) * dtde(i,j,k)
           dkemde(i,j,k) = dkemdt(i,j,k) * dtde(i,j,k)
          enddo
         enddo
        endif ! ldimen = 3
       endif ! ldimen > 1
 1000 continue
c
c=======================================================================
c------------------  End Newton-Raphson iteration loop  ----------------
c
c     If this point is reached, then Newton-Raphson did not converge 
c     within the maximum allowed number of iterations.  Try reducing
c     timestep.
c
c     1001 = target if total correction exceeds limit or ICCG failed
c            to converge => reduce timestep
c=======================================================================
c
      if(myid_w .eq. 0) then
       write(*,"('NR failed to converge with dt=',1pe12.5,
     &  '  ',/1x,'(nhy,iter,der,de,time) =',i7,i3,3e14.5)")dt,
     &  nhy,iter,epsermx,epsemx,time
      endif
C      call mpi_finalize(ierr)
C      stop
1001  continue
      nred = nred + 1
      if (nred .gt. 6) then
         write(2,"('**********  TIMESTEP REDUCTION FAILED:  ABORTING  ',
     .             '**********')")
        ifsen(1) = 1
        return
      endif
      dt = 0.25D0 * dt
      goto 1
c
c=======================================================================
c     NR apparently converged.  Check that total changes in variables
c     are less than max allowed
c=======================================================================
c
1002  continue
c
c=======================================================================
c     final energy density update
c=======================================================================
c
      do k = 1, kn
       do j = 1, jn
        do i = 1, in
         dcopy(i,j,k) = d(i,j,k)
        enddo
       enddo
      enddo
c
      do k = ks, ke
       do j = js, je
        do i = is, ie
         der  (i,j,k) = (ern  (i,j,k) - erold(i,j,k))
         de   (i,j,k) = (e n  (i,j,k) - e old(i,j,k))
         erod (i,j,k) = ern(i,j,k) / d(i,j,k)
         e od (i,j,k) = e n(i,j,k) / d(i,j,k)
         erold(i,j,k) = erold(i,j,k) / d(i,j,k)
         e old(i,j,k) = e n(i,j,k)
        enddo
       enddo
      enddo
c
c=======================================================================
c    mark the boundary values as out of date
c=======================================================================
c
      do i = 1, 6
       bvstat(i,2) = 0    ! e
       bvstat(i,6) = 0    ! er
      enddo
c
c=======================================================================
c     check that entire timestep is finished
c=======================================================================
c
      trad = trad + dt
      if (trad .lt. dthydro) then
         dt = min(dt,dthydro-trad)
         goto 2
      endif
      dt = dthydro
      if(myid_w .eq. 0) then
       if(mod(nhy,100) .eq. 0) then
        write(*,"('nhy, time, dt = ',i6,1p2d12.4)")nhy,time,dt
        write(*,"('NR / CG iter cts: ',i2,' / ',0pf7.2)")iter,
     .        (totlsit-lssave)/float(iter)
       endif
      endif
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////              I M P _ C O U P              \\\\\\\\\\
c
c=======================================================================
