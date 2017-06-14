c=======================================================================
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
      subroutine restart
c
c  PURPOSE:  Sets up restarted run
c
c  EXTERNALS:
c     PROBRES -- macroname which is defined to be the user-supplied
c                subroutine name which re-initializes user-defined
c                variables for the problem to be restarted.
c                PROBRES is undefined by default for backward 
c                compatiability with earlier versions of ZEUS-2D.
c     nudt    -- computes timestep
c
c  LOCALS:
c
c  Modified: 
c
c  Jan. 7, 1999 by JCH for ZEUS-MP
c  Mar     1999 by efh, including debug
c  January 3, 2007 by JCH: commented out unnecessary "call nudt"
c  statement at end of routine.
c-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use root
      use grid
      use field
      use bndry
      use scratch
      use gravmod
      use radiation
      use impsoln
      use opac
      use opac_law
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
      integer  :: i   , j, k, n, jone, kone
      integer  :: iord, istp
      real(rl) :: dtrat
c
      integer  :: nbl, igrid
      real(rl) :: x1min, x1max, x1rat, dx1min,
     .            x2min, x2max, x2rat, dx2min,
     .            x3min, x3max, x3rat, dx3min
c
      logical  lgrid
c
      namelist /pcon/ nlim,tlim,cpulim,tsave,mbatch
      namelist /hycon/
     . qcon,qlin,courno,dtrat,iord,istp,
     . iordd,iorde,iords1,iords2,iords3,iordb1,iordb2,iordb3,iorder,
     . istpd,istpe,istps1,istps2,istps3,istpb1,istpb2,istpb3,istper,
     . dfloor,efloor,v1floor,v2floor,v3floor,b1floor,b2floor,b3floor,
     . emf1floor,emf2floor,emf3floor,erfloor,gpfloor
c
      namelist /ggen1/ nbl,x1min,x1max,igrid,x1rat,dx1min,lgrid
      namelist /ggen2/ nbl,x2min,x2max,igrid,x2rat,dx2min,lgrid
      namelist /ggen3/ nbl,x3min,x3max,igrid,x3rat,lgrid
c
      namelist /iib/     niis, fiis
      namelist /oib/     nois, fois
      namelist /ijb/     nijs, fijs
      namelist /ojb/     nojs, fojs
      namelist /ikb/     niks, fiks
      namelist /okb/     noks, foks
c
      namelist /mpitop/ ntiles, periodic
      namelist /grvcon/tgrav,ptmass,x1ptm,x2ptm,x3ptm,gsup,guniv
      namelist /radcon/ ifld,epsme,demax,dermax,nmeiter,radth,epsrad,
     .     cnvcrit,ernom,ennom,epsmaxd,cgerrcrit,
     .     ipcflag,rmfp0,xnu,powr,rho0,t_0,kpfrac
      namelist /eqos/ gamma, ciso, mmw
      namelist /gcon/ x1fac,x2fac,x3fac,iga,jga,kga,igcon
c
c=======================================================================
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c
c------------------------  PROB CONTROL  -------------------------------
c
c   nlim   = cycles to run
c   tlim   = physical (problem) time to stop calculation
c cpulim   = CPU time in seconds to stop the calculation (default  3.6M)
c  tsave   = CPU time to reserve for terminating the run (default 30.0s)
c mbatch   = 0 interactive mode                          (default     0)
c          = 1 batch mode (does not scan for keyboard input)
c
      nlim   = 1 000 000
      tlim   = 0.0
      cpulim = 3 600 000.0
      tsave  = 30.0
      mbatch = 0
c
      nred   = 0
      if (myid_w .eq. 0) then
       read (1,pcon)
       write(2,pcon)
      endif
c
c------------------------  HYDRO CONTROL  ------------------------------
c
c  qcon   = quadratic artificial viscosity (q) constant
c  qlin   = linear    artificial viscosity (q) constant
c  courno = courant number
c  dtrat  = ratio of initial dt to dtmin (used to compute dtmin below)
c  iord   = default order of advection scheme for all variables
c  istp   = default steepening flag for all variables
c  iord** = order of advection scheme to be used for variable **
c  iostp**= steepening flag for 3rd order advection.  When istp**=1,
c           use the discontinuity detection to steepen shocks during
c           interpolation for variable ** in X1INT,X1INTFC,X2INT,X2INTFC
c  **floor = smallest value desired for variable ** on grid
c            Can also be used to set a default value for initialization.
c            Note that no attempt is made to ensure that actual values
c            stay above the floor values.
c
      if (myid_w .eq. 0) then
        qcon   = 2.0
        qlin   = 0.0
        courno = 0.5
        dtrat  = 1.0e-3
        iord   = 2
        iordd  = 0
        iorde  = 0
        iords1 = 0
        iords2 = 0
        iords3 = 0
        iordb1 = 0
        iordb2 = 0
        iordb3 = 0
        iorder = 0
        istp   = 0
        istpd  = 2
        istpe  = 2
        istps1 = 2
        istps2 = 2
        istps3 = 2
        istpb1 = 2
        istpb2 = 2
        istpb3 = 2
        istper = 2
        d floor = tiny
        e floor = tiny
        v1floor = 0.0
        v2floor = 0.0
        v3floor = 0.0
        b1floor = 0.0
        b2floor = 0.0
        b3floor = 0.0
        emf1floor= 0.0
        emf2floor= 0.0
        emf3floor= 0.0
        erfloor = tiny
        gpfloor = 0.0
c
      read (1,hycon)
      write(2,hycon)
c
c Set flags to default values unless they were set in the input deck.
c
        if(iordd  .eq. 0) iordd  = iord
        if(iorde  .eq. 0) iorde  = iord
        if(iords1 .eq. 0) iords1 = iord
        if(iords2 .eq. 0) iords2 = iord
        if(iords3 .eq. 0) iords3 = iord
        if(iordb1 .eq. 0) iordb1 = iord
        if(iordb2 .eq. 0) iordb2 = iord
        if(iordb3 .eq. 0) iordb3 = iord
        if(iorder .eq. 0) iorder = iord
c
        if(istpd  .eq. 2) istpd  = istp
        if(istpe  .eq. 2) istpe  = istp
        if(istps1 .eq. 2) istps1 = istp
        if(istps2 .eq. 2) istps2 = istp
        if(istps3 .eq. 2) istps3 = istp
        if(istps1 .eq. 2) istps1 = istp
        if(istps2 .eq. 2) istps2 = istp
        if(istps3 .eq. 2) istps3 = istp
        if(istper .eq. 2) istper = istp
c
c copy input flags to a buffer for later use and broadcasting.
c
#ifdef MPI_USED
         buf_in( 1) = qcon
         buf_in( 2) = qlin
         buf_in( 3) = courno
         buf_in( 4) = dtrat
         buf_in( 5) = dfloor
         buf_in( 6) = efloor
         buf_in( 7) = v1floor
         buf_in( 8) = v2floor
         buf_in( 9) = v3floor
         buf_in(10) = b1floor
         buf_in(11) = b2floor
         buf_in(12) = b3floor
         buf_in(13) = emf1floor
         buf_in(14) = emf2floor
         buf_in(15) = emf3floor
         buf_in(16) = erfloor
         buf_in(17) = gpfloor
         buf_in(18) = tlim
         buf_in(19) = cpulim
         buf_in(20) = tsave
c
        ibuf_in( 1) = iordd
        ibuf_in( 2) = iorde
        ibuf_in( 3) = iords1
        ibuf_in( 4) = iords2
        ibuf_in( 5) = iords3
        ibuf_in( 6) = iordb1
        ibuf_in( 7) = iordb2
        ibuf_in( 8) = iordb3
        ibuf_in( 9) = iorder
        ibuf_in(10) = istpd
        ibuf_in(11) = istpe
        ibuf_in(12) = istps1
        ibuf_in(13) = istps2
        ibuf_in(14) = istps3
        ibuf_in(15) = istpb1
        ibuf_in(16) = istpb2
        ibuf_in(17) = istpb3
        ibuf_in(18) = istper
        ibuf_in(19) = nlim
        ibuf_in(20) = mbatch
        ibuf_in(21) = idebug
#endif
      endif
c
c Broadcast pcon and hycon to the others (use arrays).
c
#ifdef MPI_USED
      call MPI_BCAST( buf_in,20,MPI_FLOAT  ,0
     &                                       ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ibuf_in,21,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if (myid_w .ne. 0) then
        qcon    =  buf_in( 1)
        qlin    =  buf_in( 2)
        courno  =  buf_in( 3)
        dtrat   =  buf_in( 4)
        d floor =  buf_in( 5)
        e floor =  buf_in( 6)
        v1floor =  buf_in( 7)
        v2floor =  buf_in( 8)
        v3floor =  buf_in( 9)
        b1floor =  buf_in(10)
        b2floor =  buf_in(11)
        b3floor =  buf_in(12)
        emf1floor=  buf_in(13)
        emf2floor=  buf_in(14)
        emf3floor=  buf_in(15)
        erfloor =  buf_in(16)
        gpfloor =  buf_in(17)
        tlim    =  buf_in(18)
        cpulim  =  buf_in(19)
        tsave   =  buf_in(20)
c
        iordd   = ibuf_in( 1)
        iorde   = ibuf_in( 2)
        iords1  = ibuf_in( 3)
        iords2  = ibuf_in( 4)
        iords3  = ibuf_in( 5)
        iordb1  = ibuf_in( 6)
        iordb2  = ibuf_in( 7)
        iordb3  = ibuf_in( 8)
        iorder  = ibuf_in( 9)
        istpd   = ibuf_in(10)
        istpe   = ibuf_in(11)
        istps1  = ibuf_in(12)
        istps2  = ibuf_in(13)
        istps3  = ibuf_in(14)
        istpb1  = ibuf_in(15)
        istpb2  = ibuf_in(16)
        istpb3  = ibuf_in(17)
        istper  = ibuf_in(18)
        nlim    = ibuf_in(19)
        mbatch  = ibuf_in(20)
        idebug  = ibuf_in(21)
      endif
#endif
c
c-------------------------  SET MPI DATATYPES  -------------------------c
c Define MPI Derived Datatypes for passing 2-D slices of 3-D arrays.
c
#ifdef MPI_USED
      call MPI_TYPE_VECTOR (  jn*kn,1 ,in
     &, MPI_FLOAT, i_slice,ierr)
      call MPI_TYPE_COMMIT (i_slice,ierr)
c
      call MPI_TYPE_VECTOR (  jn*kn,1*neqm ,in*neqm
     &, MPI_FLOAT, ils_slice,ierr)
      call MPI_TYPE_COMMIT (ils_slice,ierr)
c
      call MPI_TYPE_VECTOR (  jn*kn*nspec,1 ,in
     &, MPI_FLOAT, iab_slice,ierr)
      call MPI_TYPE_COMMIT (iab_slice,ierr)
c
      call MPI_TYPE_VECTOR (  jn*kn,1*neqm*neqm ,in*neqm*neqm
     &, MPI_FLOAT, ilsm_slice,ierr)
      call MPI_TYPE_COMMIT (ilsm_slice,ierr)
      call MPI_TYPE_VECTOR (    kn  ,in,in*jn,MPI_FLOAT
     &                                 ,j_slice,ierr)
      call MPI_TYPE_COMMIT (j_slice,ierr)
      call MPI_TYPE_VECTOR (    kn  ,in*neqm,in*jn*neqm,
     .                      MPI_FLOAT
     &                                 ,jls_slice,ierr)
      call MPI_TYPE_COMMIT (jls_slice,ierr)
      call MPI_TYPE_VECTOR (    kn  ,in*neqm*neqm,in*jn*neqm*neqm,
     .                      MPI_FLOAT
     &                                 ,jlsm_slice,ierr)
      call MPI_TYPE_COMMIT (jlsm_slice,ierr)
c
      call MPI_TYPE_VECTOR (    kn*nspec  ,in,in*jn,
     .                      MPI_FLOAT
     &                                 ,jab_slice,ierr)
      call MPI_TYPE_COMMIT (jab_slice,ierr)
c
      call MPI_TYPE_VECTOR (1,in*jn ,1    ,MPI_FLOAT
     &                                 ,k_slice,ierr)
      call MPI_TYPE_COMMIT (k_slice,ierr)
      call MPI_TYPE_VECTOR (1,in*jn*neqm ,1    ,
     .                      MPI_FLOAT
     &                                 ,kls_slice,ierr)
      call MPI_TYPE_COMMIT (kls_slice,ierr)
      call MPI_TYPE_VECTOR (1,in*jn*neqm*neqm ,1    ,
     .                      MPI_FLOAT
     &                                 ,klsm_slice,ierr)
      call MPI_TYPE_COMMIT (klsm_slice,ierr)
c
      call MPI_TYPE_VECTOR (nspec,in*jn ,in*jn*kn,
     .                      MPI_FLOAT
     &                                 ,kab_slice,ierr)
      call MPI_TYPE_COMMIT (kab_slice,ierr)
#endif /* MPI */
c
c-------------------------  READ GRID NAMELISTS  -----------------------
c
CJH   --- GRID data handled in RESTART ---
c
c------------------------  GRAVITY CONTROL  ----------------------------
c
c  Gravitational self-potentials can be included in both 1-D and 2-D
c  problems by solving the Poisson equation in the GRAVITY module.
c  Point mass potentials are included directly in the momentum eqn
c  by using a non-zero value for the variable ptmass.  Point mass
c  potentials do not require defining GRAV, do not call the GRAVITY
c  module, and are not included in the array phi but are explicitely
c  added to the momentum eqn terms in the routines STV1 and STV2.
c     tgrav  = time when gravitation is switched on
c     ptmass = fixed central point mass object
c     izero  = i index of Z=0 (x1=0) for odd symmetry case
c              (cylindrical geometry only)
c     igrijb = ijb flag (0 for     symmetric (Dirichlet) boundary      )
c                       (1 for non-symmetric             boundary whose
c                        value is calculated using multipole expansion )
c     igrojb = ojb flag ("  "      "          "         "              )
c     epsgrv = error limit              for ICCGAF
c     maxgrv = maximum iteration count  for ICCGAF
c     ks0grv = level of cyclic reduction in ICCGAF
c     phibverr = error criteria for multipole moments in PHIBV
c     phibvnmx = max number of moments to be taken in PHIBV
c     graverr = error criteria for computing a soln in GRAVITY
c
      if(myid_w .eq. 0) then
        tgrav     = 0.0        ! efh 99/05/14
        ptmass    = 0.0
        x1ptm     = x1a(is)
        x2ptm     = x2a(js)
        x3ptm     = x3a(ks)
c
       read (1,grvcon)
       write(2,grvcon)
      endif
c
c------------------------  RADIATION CONTROL  --------------------------
c
      if(myid_w .eq. 0) then
       ifld  = 1
       epsme = 1.0D-8
       demax = 0.2D0
       dermax = 0.2D0
       nmeiter = 20
       radth = 1.0D0
       epsrad = 1.0D-8
       cnvcrit = 1
       ernom   = 1.0D0
       ennom   = 1.0D0
       epsmaxd = 0.05D0
       cgerrcrit = 666
       ipcflag   = 666
       rmfp0     = huge
       xnu       = huge
       powr      = huge
       rho0      = huge
       t_0       = huge
       read (1,radcon)
       write(2,radcon)
      endif ! myid_w
c
c------------------------  EQUATION OF STATE  --------------------------
c
c      gamma = ratio of specific heats
c      ciso  = isothermal sound speed
c
      if(myid_w .eq. 0) then
        gamma = 0.0
        ciso  = 0.0
c
        read (1,eqos)
        write(2,eqos)
      endif
#ifdef MPI_USED
c
c Broadcast grvcon, radcon, and eqos to the others.
c
      if (myid_w .eq. 0) then
         buf_in( 1) = guniv
         buf_in( 2) = tgrav
         buf_in( 3) = ptmass
         buf_in( 4) = x1ptm
         buf_in( 5) = x2ptm
         buf_in( 6) = x3ptm
         buf_in( 7) = epsme
         buf_in( 8) = demax
         buf_in( 9) = dermax
         buf_in(10) = radth
         buf_in(11) = epsrad
         buf_in(12) = epsmaxd
         buf_in(13) = epsmaxd ! there used to be an obsolete parameter here
         buf_in(14) = rmfp0
         buf_in(15) = xnu
         buf_in(16) = powr
         buf_in(17) = rho0
         buf_in(18) = t_0
         buf_in(19) =  gamma
         buf_in(20) =  ciso
         buf_in(21) =  mmw
         buf_in(22) =  ernom
         buf_in(23) =  ennom
        ibuf_in( 1) = gsup
        ibuf_in( 2) = nmeiter
        ibuf_in( 3) = cnvcrit
        ibuf_in( 4) = cgerrcrit
        ibuf_in( 5) = ipcflag
        ibuf_in( 6) = ifld
      endif
      call MPI_BCAST( buf_in,23,MPI_FLOAT  ,0,comm3d,ierr)
      call MPI_BCAST(ibuf_in,6,MPI_INTEGER,0,comm3d,ierr)
      call MPI_BCAST(xwedge ,1,MPI_LOGICAL,0,comm3d,ierr)
      if (myid_w .ne. 0) then
        guniv      =  buf_in( 1)
        tgrav      =  buf_in( 2)
        ptmass     =  buf_in( 3)
        x1ptm      =  buf_in( 4)
        x2ptm      =  buf_in( 5)
        x3ptm      =  buf_in( 6)
        epsme      =  buf_in( 7)
        demax      =  buf_in( 8)
        dermax     =  buf_in( 9)
        radth      =  buf_in(10)
        epsrad     =  buf_in(11)
        epsmaxd    =  buf_in(12)
        epsmaxd    =  buf_in(13) ! there used to be an obsolete parameter here
        rmfp0      =  buf_in(14)
        xnu        =  buf_in(15)
        powr       =  buf_in(16)
        rho0       =  buf_in(17)
        t_0        =  buf_in(18)
        gamma      =  buf_in(19)
        ciso       =  buf_in(20)
        mmw        =  buf_in(21)
        ernom      =  buf_in(22)
        ennom      =  buf_in(23)
        gsup       = ibuf_in( 1)
        nmeiter    = ibuf_in( 2)
        cnvcrit    = ibuf_in( 3)
        cgerrcrit  = ibuf_in( 4)
        ipcflag    = ibuf_in( 5)
        ifld       = ibuf_in( 6)
      endif
#endif /* MPI */
      gamm1 = gamma - 1.0
c
c--------------RADIATION TILE BOUNDARY ARRAYS -------------------------
c
c  For the loops below, set the boundary flags to zero for internal
c  boundaries.  Copy the constant values into the 2-D arrays.
c
c----------------------  Constant values for IIB  ---------------------
c
      if(lrad .ne. 0) then
       do 40 k=ks-2,ke+3
         do 30 j=js-2,je+3
           if (coords(1) .eq. 0) then
             liib    (j,k) = niis(2)
           else
             liib    (j,k) = 0
           endif
30       continue
40     continue
c
c-----------------------  Constant values for OIB  ---------------------
c
       do 60 k=ks-2,ke+3
         do 50 j=js-2,je+3
           if (coords(1) .eq. ntiles(1) - 1) then
             loib    (j,k) = nois(2)
           else
             loib    (j,k) = 0
           endif
50       continue
60     continue
c
c-----------------------  Constant values for IJB  ---------------------
c
       do 80 k=ks-2,ke+3
         do 70 i=is-2,ie+3
           if (coords(2) .eq. 0) then
             lijb    (i,k) = nijs(2)
           else
             lijb    (i,k) = 0
           endif
70       continue
80     continue
c
c-----------------------  Constant values for OJB  ---------------------
c
       do 100 k=ks-2,ke+3
         do 90 i=is-2,ie+3
           if (coords(2) .eq. ntiles(2) - 1) then
             lojb    (i,k) = nojs(2)
           else
             lojb    (i,k) = 0
           endif
90       continue
100    continue
c
c-----------------------  Constant values for IKB  ---------------------
c
       do 120 j=js-2,je+3
         do 110 i=is-2,ie+3
           if (coords(3) .eq. 0) then
             likb    (i,j) = niks(2)
           else
             likb    (i,j) = 0
           endif
110      continue
120    continue
c
c-----------------------  Constant values for OKB  ---------------------
c
       do 140 j=js-2,je+3
         do 130 i=is-2,ie+3
           if (coords(3) .eq. ntiles(3) - 1) then
             lokb    (i,j) = noks(2)
           else
             lokb    (i,j) = 0
           endif
130      continue
140    continue
      endif ! lrad
c
c-------------------------  PROBLEM RESTART  ---------------------------
c
c  PROBRES is a user defined cpp macroname representing a subroutine
c  which re-intializes all user-defined variables for the particular 
c  problem to be restarted.  PROBRES should re-initialize only the
c  variables defined by the user in PROBLEM which are not in any
c  ZEUS-2D common blocks and are therefore not saved in the restart
c  dump.  PROBRES must read the same namelist(s) as PROBLEM.
c
      call PROBRES
c
c-------------------------  GRID MOTION CONTROL  -----------------------
c
c        igcon = 1 and x1fac < 0 gives "lagrangian" tracking in x1 lines
c        igcon = 1 and x2fac < 0 gives "lagrangian" tracking in x2 lines
c        igcon = 2 for input grid boundary speeds
c                vg1(io) = x1fac * central soundspeed
c                vg2(jo) = x2fac * central soundspeed
c        igcon = 3 for constant motion at x1[2]fac
c
      if(myid_w .eq. 0) then
        read (1,gcon)
        write(2,gcon)
      endif
c
c----------------------  COMPUTE TIMESTEP  -----------------------------
      if(ldimen .eq. 3) then
       kone = 1
      else
       kone = 0
      endif 
      if(ldimen .eq. 1) then
       jone = 0 
      else
       jone = 1
      endif
c
      do 370 k=ks,ke
       do 360 j=js,je
        do 350 i=is,ie
         w3da(i,j,k) = v1(i,j,k)*0.5*(d(i-1,j     ,k     )+d(i,j,k))
         w3db(i,j,k) = v2(i,j,k)*0.5*(d(i  ,j-jone,k     )+d(i,j,k))
     1                  * g2b(i)
         w3dc(i,j,k) = v3(i,j,k)*0.5*(d(i  ,j     ,k-kone)+d(i,j,k))
     1                  * g31b(i) * g32b(j)
         if(lrad .ne. 0) then
          w3dh(i,j,k) = er(i,j,k) / d(i,j,k)
         endif ! lrad
350     continue
360    continue
370   continue
c
      if(nspec .gt. 1) then
       do n = 1, nspec
        do k = ks, ke
         do j = js, je
          do i = is, ie
           w4da(i,j,k,n) = abun(i,j,k,n)
          enddo
         enddo
        enddo
       enddo
      endif
c
      if(.not. xvgrid) then
       do i = 1, in
        vg1(i) = 0.0D0
       enddo
       do j = 1, jn
        vg2(j) = 0.0D0
       enddo
       do k = 1, kn
        vg3(k) = 0.0D0
       enddo
      endif
c
C -- JCH: this shouldn't be here.  The correct timestep is contained in the
C         restart dump!
C
C      call nudt
c
      return
      end
