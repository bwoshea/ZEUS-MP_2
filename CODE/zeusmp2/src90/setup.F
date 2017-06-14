c=======================================================================
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
      subroutine setup
c
c  PURPOSE: Sets up execution of a new run by initializing all variables
c  according to the flags and values in the input deck.  Calls PROBLEM
c  a CPP macro) to initialize field variables for the particular
c  problem to be studied, otherwise field variables are set to "tiny"
c  (parameter defined to be smallest number possible on machine).  Note
c  user must define PROBLEM to be the appropriate subroutine name in the
c  file zeusmp.def or in the Make_zeusmp command line.
c
c  Order in which namelists are read has been changed from ZEUS-2D so 
c  that the boundary conditions are known before the MPI virtual 
c  topology is defined and before the grid is computed.  MPI needs to 
c  know if the grid is periodic, while the grid should be symmetric 
c  across reflecting boundaries but constant across flow in/out and 
c  periodic for periodic ones.
c
c  EXTERNALS:
c     ggen    -- initializes grid according to input deck
c     bval*   -- boundary value routines
c     nudt    -- computes initial timestep
c     newgrid -- computes new grid position for moving grid
c     scopy
c
c  LOCALS:
c
c  LAST MODIFIED: JCH 3/13/97.
c-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use root
      use grid
      use field
      use bndry
      use radiation
      use opac
      use opac_law
      use scratch
      use gravmod
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use cons
      use impsoln
c
      implicit NONE
c
      integer  :: i, j, k, n, jone, kone
      integer  :: iord,istp
      real(rl) :: dtrat
c
      integer imax,ISMAX
      real(rl) :: dmax
c
      namelist /pcon/ nlim,tlim,cpulim,tsave,mbatch
      namelist /hycon/
     . qcon,qlin,courno,dtrat,iord,istp,
     . iordd,iorde,iords1,iords2,iords3,iordb1,iordb2,iordb3,iorder,
     . istpd,istpe,istps1,istps2,istps3,istpb1,istpb2,istpb3,istper,
     . dfloor,efloor,v1floor,v2floor,v3floor,b1floor,b2floor,b3floor,
     . emf1floor,emf2floor,emf3floor,erfloor,gpfloor
c
      namelist /iib/     niis, fiis
      namelist /oib/     nois, fois
      namelist /ijb/     nijs, fijs
      namelist /ojb/     nojs, fojs
      namelist /ikb/     niks, fiks
      namelist /okb/     noks, foks
c
      namelist /mpitop/ ntiles, periodic
      namelist /grvcon/ guniv,tgrav,ptmass,x1ptm,x2ptm,x3ptm,xwedge,
     .                  gsup
      namelist /radcon/ ifld,epsme,demax,dermax,nmeiter,radth,epsrad,
     .     cnvcrit,ernom,ennom,epsmaxd,cgerrcrit,
     .     ipcflag,rmfp0,xnu,powr,rho0,t_0,kpfrac
      namelist /eqos/ gamma, ciso, mmw
      namelist /gcon/ x1fac,x2fac,x3fac,iga,jga,kga,igcon
c=======================================================================
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c
c------------------------  PROBLEM CONTROL  ---------------------------
c
c   nlim   = cycles to run                               (default    1M)
c   tlim   = physical (problem) time to stop calculation (default  0.0 )
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
c  qcon   = quadratic artificial viscosity (q) constant (default 2.0)
c  qlin   = linear    artificial viscosity (q) constant (default 0.0)
c  courno = courant number                              (default 0.5)
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
      endif
#ifdef MPI_USED
c
c Broadcast pcon and hycon to the others (use arrays).
c
      call MPI_BCAST( buf_in,20,MPI_FLOAT  ,0
     &                                       ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ibuf_in,20,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
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
      endif
#endif /* MPI */
c
c------------------------  BOUNDARY CONTROL ----------------------------
c
c  The following points describecoords(1), how boundaries are handled:
c
c  1)  Any of 6 mhd boundary conditions may be specified independently
c  at every zone on the physical problem boundary.  The boundary type is
c  specified by nflo, where
c
c      nflo = 0  =>  interior boundary (get data from neighboring tile)
c           = 1  =>  reflecting (v(normal) = b(normal) = 0)
c           =-1  =>  reflecting (XYZ: same as 1; ZRP: same as 1 with
c                    inversion of 3-components at ijb; RTP: same as 1
c                    with inversion of 2- and 3-components at iib and
c                    inversion of 3-components at ijb and ojb.)
c           = 2  =>  flow out
c           = 3  =>  flow in
c           = 4  =>  periodic
c           = 5  =>  reflecting (v(normal) = 0, b(tangential) = 0)
c
c  Note that in ZRP and RTP, some boundary conditions are implied by
c  the choice of limits.  e.g., if 0 .le. x3a .le. 2*pi in either ZRP
c  or RTP, then periodic boundary conditions should be imposed.
c  Set "niib" to -1 (reflecting with inversion of 2- and 3-components) 
c  if the inner i boundary is at the   origin (RTP).
c  Set "nijb" to -1 (reflecting with inversion of 3-components) if
c     the inner j boundary is on the "Z" axis (ZRP or RTP).
c  Set "nojb" to -1 (reflecting with inversion of 3-components) if
c     the outer j boundary is on the "Z" axis (RTP).
c
c  Since the grid is staggered, the boundary conditions apply over
c  slightly different regions, depending on the centering of the 
c  variable c  and the type of boundary condition used.  Thus, "niib" 
c  is applied to zone or 1-face centered quantities (d, v1, e, gp, b1),
c  "niib2" is applied to 2-face centred quantities (v2, b2), 
c  "niib3" is applied to 3-face centred quantities (v3, b3), and 
c  "niib23" is appled to corner centred quantities (emf1, emf2, emf3).
c  Note that the secondary boundary integer flags are determined from 
c  "niib" automatically.
c
c  Since constant values are often used and are easy to input via 
c  namelists, a single scalar variable is read in for nflo on each 
c  boundary.  These are read in as:
c
c      niis(1),nois(1),nijs(1),nojs(1),niks(1),noks(1)
c
c  For more complicated boundary conditions, values of nflo at each 
c  boundary zone are stored in the 2-D arrays:
c
c      niib,noib,nijb,nojb,nikb,nokb
c
c  These arrays are automatically set to the constant input values.  For
c  more complicated problems, these arrays may be altered in the Problem
c  Generator.
c
c  2) Since the radiation boundary types may differ from the fluid 
c  boundary types, we may specify the former independently of the 
c  latter.  Thus, the radiation boundary types are specified by lflo,
c  where
c
c      lflo = 0  =>  interior boundary (get data from neighboring tile)
c           = 1  =>  reflecting
c           = 2  =>  flow out
c           = 3  =>  flow in
c           = 4  =>  periodic
c
c  Constant values for lflo are read in as:
c
c      niis(2),nois(2),nijs(2),nojs(2),niks(2),noks(2)
c
c  while the 2-D arrays are liib,loib,lijb,lojb,likb,lokb.
c
c  3) Boundary conditions on the gravitational potential are
c     specified by igr, where
c
c       igr = 0  =>  interior boundary (get data from neighboring tile)
c           = 1  =>  reflecting (dgp/d(normal) = 0 "von Neumann")
c           = 2  =>  outflow (equivalent to reflecting)
c           = 3  =>  gp specified (Dirichlet)
c           = 4  =>  periodic
c
c     The flags are defined by analogy with the hydro nflo flag.  This
c     is quite different from ZEUS-2D.  The flags igr are read in as:
c
c      niis(3),nois(3),nijs(3),nojs(3),niks(3),noks(3)
c
c     No 2-D integer arrays for gravity BCs, I guess.
c
c  4) For flow-in boundaries, boundary values of d,e,v1,v2,[v3],...
c  must be input.  Since constant values are often used and are easy to 
c  input via namelists, a single scalar variable is read in for each 
c  field variable on each boundary, for example, the quantities
c
c    fiis(1),fois(1),fijs(1),fojs(1),fiks(1),foks(1) 
c
c  give the constant boundary values for d.  There is a set of f's
c  for e, v, b, er, and gp.
c
c  For more complicated inflow boundary conditions, arrays are used to 
c  store specified values at each boundary zone for each function; for 
c  example:
c
c    diib(j,k,2) is inner i boundary density for sweep j,k at ism2
c    diib(j,k,1) is inner i boundary density for sweep j,k at ism1
c    doib(j,k,1) is outer i boundary density for sweep j,k at iep1
c    doib(j,k,2) is outer i boundary density for sweep j,k at iep2
c
c    dijb(k,i,2) is inner j boundary density for sweep k,i at jsm2
c    dijb(k,i,1) is inner j boundary density for sweep k,i at jsm1
c    dojb(k,i,1) is outer j boundary density for sweep k,i at jep1
c    dojb(k,i,2) is outer j boundary density for sweep k,i at jep2
c
c    dikb(i,j,2) is inner k boundary density for sweep i,j at ksm2
c    dikb(i,j,1) is inner k boundary density for sweep i,j at ksm1
c    dokb(i,j,1) is outer k boundary density for sweep i,j at kep1
c    dokb(i,j,2) is outer k boundary density for sweep i,j at kep2
c
c    These arrays are initialized automatically to constant values;
c    the user's Problem Generator may alter these arrays, e.g. to
c    specify an inlet for a jet or for certain advection tests.
c
c 5) For corner zones [(ii-1,ji-1),(io+1,ji-1),etc] there is a
c    "pecking order" of precedence which attempts to pick the BC
c    that provides the most stable solution.
c    THE USER MAY WANT TO OVERIDE THIS CHOICE in the Problem Generator.
c
c  Set defaults and read in boundary values.  buf_in still contains
c  the "floor" or default values for each field variable.
c
      if (myid_w .eq. 0) then
        do 10 i=1,3
          niis(i) = 0
          nijs(i) = 0
          niks(i) = 0
          nois(i) = 0
          nojs(i) = 0
          noks(i) = 0
10      continue
        do 20 i=1,nbvar
          fiis(i) = buf_in(i+4)
          fois(i) = buf_in(i+4)
          fijs(i) = buf_in(i+4)
          fojs(i) = buf_in(i+4)
          fiks(i) = buf_in(i+4)
          foks(i) = buf_in(i+4)
20      continue
c
        read (1,iib)
        write(2,iib)
c
        read (1,oib)
        write(2,oib)
c
        read (1,ijb)
        write(2,ijb)
c
        read (1,ojb)
        write(2,ojb)
c
        read (1,ikb)
        write(2,ikb)
c
        read (1,okb)
        write(2,okb)
c
      endif
c
#ifdef MPI_USED
c
c Tell the others what the master has read.
c
      call MPI_BCAST(niis, 3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nois, 3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nijs, 3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nojs, 3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(niks, 3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(noks, 3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
      call MPI_BCAST(fiis,nbvar,MPI_FLOAT  ,0
     &                                    ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fois,nbvar,MPI_FLOAT  ,0
     &                                    ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fijs,nbvar,MPI_FLOAT  ,0
     &                                    ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fojs,nbvar,MPI_FLOAT  ,0
     &                                    ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fiks,nbvar,MPI_FLOAT  ,0
     &                                    ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(foks,nbvar,MPI_FLOAT  ,0
     &                                    ,MPI_COMM_WORLD,ierr)
#endif /* MPI */
c
c-------------------------  SET BOUNDARY TYPES  ------------------------
c
      if (coords(1) .gt. 0) then
        niis(1) = 0
        niis(2) = 0
        niis(3) = 0
      endif
      if (coords(2) .gt. 0) then
        nijs(1) = 0
        nijs(2) = 0
        nijs(3) = 0
      endif
      if (coords(3) .gt. 0) then
        niks(1) = 0
        niks(2) = 0
        niks(3) = 0
      endif
      if (coords(1) .lt. ntiles(1) - 1) then
        nois(1) = 0
        nois(2) = 0
        nois(3) = 0
      endif
      if (coords(2) .lt. ntiles(2) - 1) then
        nojs(1) = 0
        nojs(2) = 0
        nojs(3) = 0
      endif
      if (coords(3) .lt. ntiles(3) - 1) then
        noks(1) = 0
        noks(2) = 0
        noks(3) = 0
      endif
c
c-------------------------  INITIALIZE GRID  ---------------------------
c
c  Routine ggen reads ggen1, ggen2, ggen3; computes the grid for each 
c  tile.
c
#ifdef MPI_USED
      call MPI_BARRIER(comm3d, ierr)
#endif
      call ggen
#ifdef MPI_USED
c
c-------------------------  SET MPI DATATYPES  -------------------------
c
c Define MPI Derived Datatypes for passing 2-D slices of 3-D arrays.
c
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
c
      call MPI_TYPE_VECTOR (    kn  ,in,in*jn,MPI_FLOAT
     &                                 ,j_slice,ierr)
      call MPI_TYPE_COMMIT (j_slice,ierr)
c
      call MPI_TYPE_VECTOR (    kn  ,in*neqm,in*jn*neqm,
     .                      MPI_FLOAT
     &                                 ,jls_slice,ierr)
      call MPI_TYPE_COMMIT (jls_slice,ierr)
c
      call MPI_TYPE_VECTOR (    kn*nspec  ,in,in*jn,
     .                      MPI_FLOAT
     &                                 ,jab_slice,ierr)
      call MPI_TYPE_COMMIT (jab_slice,ierr)
c
      call MPI_TYPE_VECTOR (    kn  ,in*neqm*neqm,in*jn*neqm*neqm,
     .                      MPI_FLOAT
     &                                 ,jlsm_slice,ierr)
      call MPI_TYPE_COMMIT (jlsm_slice,ierr)
c
      call MPI_TYPE_VECTOR (1,in*jn ,1    ,MPI_FLOAT
     &                                 ,k_slice,ierr)
      call MPI_TYPE_COMMIT (k_slice,ierr)
c
      call MPI_TYPE_VECTOR (1,in*jn*neqm ,1    ,
     .                      MPI_FLOAT
     &                                 ,kls_slice,ierr)
      call MPI_TYPE_COMMIT (kls_slice,ierr)
c
      call MPI_TYPE_VECTOR (nspec,in*jn ,in*jn*kn,
     .                      MPI_FLOAT
     &                                 ,kab_slice,ierr)
      call MPI_TYPE_COMMIT (kab_slice,ierr)
c
      call MPI_TYPE_VECTOR (1,in*jn*neqm*neqm ,1    ,
     .                      MPI_FLOAT
     &                                 ,klsm_slice,ierr)
      call MPI_TYPE_COMMIT (klsm_slice,ierr)
#endif /* MPI */
c
c------------------------  GRAVITY CONTROL  ----------------------------
c
c  Self-gravity can be included
c  by solving the Poisson equation in the GRAVITY module.
c  Point mass potentials are included directly in the momentum eqn
c  by using a non-zero value for the variable ptmass.  Point mass
c  potentials do not require defining GRAV, do not call the GRAVITY
c  module, and are not included in the array phi but are explicitly
c  added to the momentum eqn terms in FORCES.
c     g              gravitational constant
c     tgrav          time when gravitation is switched on
c     ptmass         mass of fixed point mass 
c     x1ptm          x1 coordinate of the point mass
c     x2ptm
c     x3ptm
c
c  Boundary flags more general than the ones below are input in the
c  BC namelists described above.  Keep these for now.  Strange, though.
c
c     igrijb = ijb flag (0 for     symmetric (Dirichlet) boundary      )
c                       (1 for non-symmetric             boundary whose
c                        value is calculated using multipole expansion )
c     igrojb = ojb flag ("  "      "          "         "              )
c
c  These are for the multigrid solver.  R: real, I: integer
c
c     gsup       I   Cycle number to calculate gravitational potential
      if (myid_w .eq. 0) then
        guniv     = 6.667d-8
        tgrav     = 0.0
        gsup      = 1
        xwedge    = .false.
c
        read (1,grvcon)
        write(2,grvcon)
      endif
      igcall = 0
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
      if (myid_w .eq. 0) then
        gamma = 0.0
        ciso  = 0.0
        mmw   = 1.0
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
c-----------------------  SET CONSTANT VALUES IN BOUNDARY ARRAYS -------
c
c  For the loops below, set the boundary flags to zero for internal
c  boundaries.  Copy the constant values into the 2-D arrays.
c
c-----------------------  Constant values for IIB  ---------------------
c
       do 40 k=ks-2,ke+3
         do 30 j=js-2,je+3
           if (coords(1) .eq. 0) then
             niib    (j,k) = niis(1)
            if(lrad .ne. 0) then
             liib    (j,k) = niis(2)
            endif
C             if(xmhd) then
              niib2   (j,k) = niis(1)
              niib3   (j,k) = niis(1)
              niib23  (j,k) = niis(1)
C             endif
           else
             niib    (j,k) = 0
            if(lrad .ne. 0) then
             liib    (j,k) = 0
            endif
C             if(xmhd) then
              niib2   (j,k) = 0
              niib3   (j,k) = 0
              niib23  (j,k) = 0
C             endif
           endif
           diib    (j,k,1) = fiis(1)
           diib    (j,k,2) = huge
           diib    (j,k,3) = huge
          if(xiso .eqv. .false.) then
           eiib    (j,k,1) = fiis(2)
           eiib    (j,k,2) = huge
          endif
           v1iib   (j,k,1) = fiis(3)
           v1iib   (j,k,2) = huge
           v2iib   (j,k,1) = fiis(4)
           v2iib   (j,k,2) = huge
           v3iib   (j,k,1) = fiis(5)
           v3iib   (j,k,2) = huge
          if(xmhd) then
           b1iib   (j,k,1) = fiis(6)
           b1iib   (j,k,2) = huge
           b2iib   (j,k,1) = fiis(7)
           b2iib   (j,k,2) = huge
           b3iib   (j,k,1) = fiis(8)
           b3iib   (j,k,2) = huge
           emf1iib (j,k,1) = fiis(9)
           emf1iib (j,k,2) = huge
           emf1iib (j,k,3) = huge
           emf2iib (j,k,1) = fiis(10)
           emf2iib (j,k,2) = huge
           emf2iib (j,k,3) = huge
           emf3iib (j,k,1) = fiis(11)
           emf3iib (j,k,2) = huge
           emf3iib (j,k,3) = huge
          endif ! xmhd
          if(lrad .ne. 0) then
           eriib   (j,k,1) = fiis(12)
           eriib   (j,k,2) = huge
          endif
          if(xgrav .or. xgrvfft)
     .     gpiib   (j,k,1) = fiis(13)
30       continue
40     continue
c
c-----------------------  Constant values for OIB  ---------------------
c
       do 60 k=ks-2,ke+3
         do 50 j=js-2,je+3
           if (coords(1) .eq. ntiles(1) - 1) then
             noib    (j,k) = nois(1)
            if(lrad .ne. 0) then
             loib    (j,k) = nois(2)
            endif
C             if(xmhd) then
              noib2   (j,k) = nois(1)
              noib3   (j,k) = nois(1)
              noib23  (j,k) = nois(1)
C             endif
           else
             noib    (j,k) = 0
            if(lrad .ne. 0) then
             loib    (j,k) = 0
            endif
C             if(xmhd) then
              noib2   (j,k) = 0
              noib3   (j,k) = 0
              noib23  (j,k) = 0
C             endif
           endif
           doib    (j,k,1) = fois(1)
           doib    (j,k,2) = huge
           doib    (j,k,3) = huge
          if(xiso .eqv. .false.) then
           eoib    (j,k,1) = fois(2)
           eoib    (j,k,2) = huge
          endif ! xiso
           v1oib   (j,k,1) = fois(3)
           v1oib   (j,k,2) = huge
           v2oib   (j,k,1) = fois(4)
           v2oib   (j,k,2) = huge
           v3oib   (j,k,1) = fois(5)
           v3oib   (j,k,2) = huge
          if(xmhd) then
           b1oib   (j,k,1) = fois(6)
           b1oib   (j,k,2) = huge
           b2oib   (j,k,1) = fois(7)
           b2oib   (j,k,2) = huge
           b3oib   (j,k,1) = fois(8)
           b3oib   (j,k,2) = huge
           emf1oib (j,k,1) = fois(9)
           emf1oib (j,k,2) = huge
           emf1oib (j,k,3) = huge
           emf2oib (j,k,1) = fois(10)
           emf2oib (j,k,2) = huge
           emf2oib (j,k,3) = huge
           emf3oib (j,k,1) = fois(11)
           emf3oib (j,k,2) = huge
           emf3oib (j,k,3) = huge
          endif ! xmhd
          if(lrad .ne. 0) then
           eroib   (j,k,1) = fois(12)
           eroib   (j,k,2) = huge
          endif
          if(xgrav .or. xgrvfft)
     .     gpoib   (j,k,1) = fois(13)
50       continue
60     continue
c
c-----------------------  Constant values for IJB  ---------------------
c
       do 80 k=ks-2,ke+3
         do 70 i=is-2,ie+3
           if (coords(2) .eq. 0) then
             nijb    (i,k) = nijs(1)
            if(lrad .ne. 0) then
             lijb    (i,k) = nijs(2)
            endif
C             if(xmhd) then
              nijb3   (i,k) = nijs(1)
              nijb1   (i,k) = nijs(1)
              nijb31  (i,k) = nijs(1)
C             endif
           else
             nijb    (i,k) = 0
            if(lrad .ne. 0) then
             lijb    (i,k) = 0
            endif
C             if(xmhd) then
              nijb3   (i,k) = 0
              nijb1   (i,k) = 0
              nijb31  (i,k) = 0
C             endif
           endif
           dijb    (i,k,1) = fijs(1)
           dijb    (i,k,2) = huge
           dijb    (i,k,3) = huge
          if(xiso .eqv. .false.) then
           eijb    (i,k,1) = fijs(2)
           eijb    (i,k,2) = huge
          endif ! xiso
           v1ijb   (i,k,1) = fijs(3)
           v1ijb   (i,k,2) = huge
           v2ijb   (i,k,1) = fijs(4)
           v2ijb   (i,k,2) = huge
           v3ijb   (i,k,1) = fijs(5)
           v3ijb   (i,k,2) = huge
          if(xmhd) then
           b1ijb   (i,k,1) = fijs(6)
           b1ijb   (i,k,2) = huge
           b2ijb   (i,k,1) = fijs(7)
           b2ijb   (i,k,2) = huge
           b3ijb   (i,k,1) = fijs(8)
           b3ijb   (i,k,2) = huge
           emf1ijb (i,k,1) = fijs(9)
           emf1ijb (i,k,2) = huge
           emf1ijb (i,k,3) = huge
           emf2ijb (i,k,1) = fijs(10)
           emf2ijb (i,k,2) = huge
           emf2ijb (i,k,3) = huge
           emf3ijb (i,k,1) = fijs(11)
           emf3ijb (i,k,2) = huge
           emf3ijb (i,k,3) = huge
          endif ! xmhd
          if(lrad .ne. 0) then
           erijb   (i,k,1) = fijs(12)
           erijb   (i,k,2) = huge
          endif
          if(xgrav .or. xgrvfft)
     .     gpijb   (i,k,1) = fijs(13)
70       continue
80     continue
c
c-----------------------  Constant values for OJB  ---------------------
c
       do 100 k=ks-2,ke+3
         do 90 i=is-2,ie+3
           if (coords(2) .eq. ntiles(2) - 1) then
             nojb    (i,k) = nojs(1)
            if(lrad .ne. 0) then
             lojb    (i,k) = nojs(2)
            endif
C             if(xmhd) then
              nojb3   (i,k) = nojs(1)
              nojb1   (i,k) = nojs(1)
              nojb31  (i,k) = nojs(1)
C             endif
           else
             nojb    (i,k) = 0
            if(lrad .ne. 0) then
             lojb    (i,k) = 0
            endif
C             if(xmhd) then
              nojb3   (i,k) = 0
              nojb1   (i,k) = 0
              nojb31  (i,k) = 0
C             endif
           endif
           dojb    (i,k,1) = fojs(1)
           dojb    (i,k,2) = huge
           dojb    (i,k,3) = huge
          if(xiso .eqv. .false.) then
           eojb    (i,k,1) = fojs(2)
           eojb    (i,k,2) = huge
          endif ! xiso
           v1ojb   (i,k,1) = fojs(3)
           v1ojb   (i,k,2) = huge
           v2ojb   (i,k,1) = fojs(4)
           v2ojb   (i,k,2) = huge
           v3ojb   (i,k,1) = fojs(5)
           v3ojb   (i,k,2) = huge
          if(xmhd) then
           b1ojb   (i,k,1) = fojs(6)
           b1ojb   (i,k,2) = huge
           b2ojb   (i,k,1) = fojs(7)
           b2ojb   (i,k,2) = huge
           b3ojb   (i,k,1) = fojs(8)
           b3ojb   (i,k,2) = huge
           emf1ojb (i,k,1) = fojs(9)
           emf1ojb (i,k,2) = huge
           emf1ojb (i,k,3) = huge
           emf2ojb (i,k,1) = fojs(10)
           emf2ojb (i,k,2) = huge
           emf2ojb (i,k,3) = huge
           emf3ojb (i,k,1) = fojs(11)
           emf3ojb (i,k,2) = huge
           emf3ojb (i,k,3) = huge
          endif
          if(lrad .ne. 0) then
           erojb   (i,k,1) = fojs(12)
           erojb   (i,k,2) = huge
          endif
          if(xgrav .or. xgrvfft)
     .     gpojb   (i,k,1) = fojs(13)
90       continue
100    continue
c
c-----------------------  Constant values for IKB  ---------------------
c
       do 120 j=js-2,je+3
         do 110 i=is-2,ie+3
           if (coords(3) .eq. 0) then
             nikb    (i,j) = niks(1)
            if(lrad .ne. 0) then
             likb    (i,j) = niks(2)
            endif
C             if(xmhd) then
              nikb1   (i,j) = niks(1)
              nikb2   (i,j) = niks(1)
              nikb12  (i,j) = niks(1)
C             endif
           else
             nikb    (i,j) = 0
            if(lrad .ne. 0) then
             likb    (i,j) = 0
            endif
C            if(xmhd) then
             nikb1   (i,j) = 0
             nikb2   (i,j) = 0
             nikb12  (i,j) = 0
C            endif
           endif
           dikb    (i,j,1) = fiks(1)
           dikb    (i,j,2) = huge
           dikb    (i,j,3) = huge
          if(xiso .eqv. .false.) then
           eikb    (i,j,1) = fiks(2)
           eikb    (i,j,2) = huge
          endif
           v1ikb   (i,j,1) = fiks(3)
           v1ikb   (i,j,2) = huge
           v2ikb   (i,j,1) = fiks(4)
           v2ikb   (i,j,2) = huge
           v3ikb   (i,j,1) = fiks(5)
           v3ikb   (i,j,2) = huge
          if(xmhd) then
           b1ikb   (i,j,1) = fiks(6)
           b1ikb   (i,j,2) = huge
           b2ikb   (i,j,1) = fiks(7)
           b2ikb   (i,j,2) = huge
           b3ikb   (i,j,1) = fiks(8)
           b3ikb   (i,j,2) = huge
           emf1ikb (i,j,1) = fiks(9)
           emf1ikb (i,j,2) = huge
           emf1ikb (i,j,3) = huge
           emf2ikb (i,j,1) = fiks(10)
           emf2ikb (i,j,2) = huge
           emf2ikb (i,j,3) = huge
           emf3ikb (i,j,1) = fiks(11)
           emf3ikb (i,j,2) = huge
           emf3ikb (i,j,3) = huge
          endif ! xmhd
          if(lrad .ne. 0) then
           erikb   (i,j,1) = fiks(12)
           erikb   (i,j,2) = huge
          endif
          if(xgrav .or. xgrvfft)
     .     gpikb   (i,j,1) = fiks(13)
110      continue
120    continue
c
c-----------------------  Constant values for OKB  ---------------------
c
       do 140 j=js-2,je+3
         do 130 i=is-2,ie+3
           if (coords(3) .eq. ntiles(3) - 1) then
            nokb    (i,j) = noks(1)
            if(lrad .ne. 0) then
             lokb    (i,j) = noks(2)
            endif
C            if(xmhd) then
             nokb1   (i,j) = noks(1)
             nokb2   (i,j) = noks(1)
             nokb12  (i,j) = noks(1)
C            endif
           else
            nokb    (i,j) = 0
            if(lrad .ne. 0) then
             lokb    (i,j) = 0
            endif
C            if(xmhd) then
             nokb1   (i,j) = 0
             nokb2   (i,j) = 0
             nokb12  (i,j) = 0
C            endif
           endif
           dokb    (i,j,1) = foks(1)
           dokb    (i,j,2) = huge
           dokb    (i,j,3) = huge
          if(xiso .eqv. .false.) then
           eokb    (i,j,1) = foks(2)
           eokb    (i,j,2) = huge
          endif
           v1okb   (i,j,1) = foks(3)
           v1okb   (i,j,2) = huge
           v2okb   (i,j,1) = foks(4)
           v2okb   (i,j,2) = huge
           v3okb   (i,j,1) = foks(5)
           v3okb   (i,j,2) = huge
          if(xmhd) then
           b1okb   (i,j,1) = foks(6)
           b1okb   (i,j,2) = huge
           b2okb   (i,j,1) = foks(7)
           b2okb   (i,j,2) = huge
           b3okb   (i,j,1) = foks(8)
           b3okb   (i,j,2) = huge
           emf1okb (i,j,1) = foks(9)
           emf1okb (i,j,2) = huge
           emf1okb (i,j,3) = huge
           emf2okb (i,j,1) = foks(10)
           emf2okb (i,j,2) = huge
           emf2okb (i,j,3) = huge
           emf3okb (i,j,1) = foks(11)
           emf3okb (i,j,2) = huge
           emf3okb (i,j,3) = huge
          endif ! xmhd
          if(lrad .ne. 0) then
           erokb   (i,j,1) = foks(12)
           erokb   (i,j,2) = huge
          endif
          if(xgrav .or. xgrvfft)
     .     gpokb   (i,j,1) = foks(13)
130      continue
140    continue
c
c-------------------------  PROBLEM GENERATOR  -------------------------
c
c  PROBLEM is a user-defined cpp macro name representing a subroutine
c  which intializes all field variables for the particular problem to
c  be studied.  PROBLEM should initialize the field variable arrays 
c  for both active zones and at least the first layer of boundary zones,
c  unless the default or input constant values already specify the 
c  desired problem.
c
c  For non-uniform initial magnetic field configurations, to satisfy 
c  the constraint DIV(B)=0 PROBLEM should initialize
c  b1, b2, b3 by differencing a vector potential.
c
c  First initialize all field variables to default (input) values.
c
      do 170 k=ks-2,ke+3
        do 160 j=js-2,je+3
          do 150 i=is-2,ie+3
            d  (i,j,k) = dfloor
            if(xiso .eqv. .false.) e  (i,j,k) = efloor
            v1 (i,j,k) = v1floor
            v2 (i,j,k) = v2floor
            v3 (i,j,k) = v3floor
           if(xmhd) then
            b1 (i,j,k) = b1floor * g2bi(i) * g31bi(i)
            b2 (i,j,k) = b2floor * g32bi(j)
            b3 (i,j,k) = b3floor * g2bi(i)
           endif ! xmhd
           if(lrad .ne. 0) then
            er (i,j,k) = erfloor
           endif
150       continue
160     continue
170   continue
c
c Set the bvstat array to 0, indicating that the boundary values
c will need to be updated.
c
      do 175 j=1,nbvar
        do 175 i=1,6
          bvstat(i,j) = 0
175   continue
      dt     = huge
c
      call PROBLEM
c
c Note that the 2-D boundary arrays were set to constant values above,
c which should be overwritten with any desired spatially-dependent 
c functions in the Problem Generator.
c
c Set the scalar boundary value flags niis, etc. to nonzero values if 
c the corresponding boundaries are physical ones.
c
       if (niib(js,ks).ne.0 .and. coords(1).eq.0        ) 
     &                              niis(1) = niib(js,ks)
       if (noib(js,ks).ne.0 .and. coords(1).eq.ntiles(1) - 1) 
     &                              nois(1) = noib(js,ks)
       if (nijb(is,ks).ne.0 .and. coords(2).eq.0        ) 
     &                              nijs(1) = nijb(is,ks)
       if (nojb(is,ks).ne.0 .and. coords(2).eq.ntiles(2) - 1) 
     &                              nojs(1) = nojb(is,ks)
       if (nikb(is,js).ne.0 .and. coords(3).eq.0        ) 
     &                              niks(1) = nikb(is,js)
       if (nokb(is,js).ne.0 .and. coords(3).eq.ntiles(3) - 1) 
     &                              noks(1) = nokb(is,js)
c
c-----------------------  INITIALIZE ADDITIONAL BOUNDARY ARRAYS --------
c
c      Set the secondary boundary flags ("niib2", "niib3", "niib23",
c  etc.) from the primary boundary flags ("niib", etc.).
c
       call bndyflgs
c
c      If the second (and third) layer boundaries were not set by the
c  user, set them equal to the first layer boundary values.
c
       do 190 k=ks-2,ke+3
         do 180 j=js-2,je+3
           if (   diib(j,k,2) .eq. huge)    diib(j,k,2) =    diib(j,k,1)
           if (   diib(j,k,3) .eq. huge)    diib(j,k,3) =    diib(j,k,2)
          if(xiso .eqv. .false.) then
           if (   eiib(j,k,2) .eq. huge)    eiib(j,k,2) =    eiib(j,k,1)
          endif ! xiso
           if (  v1iib(j,k,2) .eq. huge)   v1iib(j,k,2) =   v1iib(j,k,1)
           if (  v2iib(j,k,2) .eq. huge)   v2iib(j,k,2) =   v2iib(j,k,1)
           if (  v3iib(j,k,2) .eq. huge)   v3iib(j,k,2) =   v3iib(j,k,1)
          if(xmhd) then
           if (  b1iib(j,k,2) .eq. huge)   b1iib(j,k,2) =   b1iib(j,k,1)
           if (  b2iib(j,k,2) .eq. huge)   b2iib(j,k,2) =   b2iib(j,k,1)
           if (  b3iib(j,k,2) .eq. huge)   b3iib(j,k,2) =   b3iib(j,k,1)
           if (emf1iib(j,k,2) .eq. huge) emf1iib(j,k,2) = emf1iib(j,k,1)
           if (emf1iib(j,k,3) .eq. huge) emf1iib(j,k,3) = emf1iib(j,k,2)
           if (emf2iib(j,k,2) .eq. huge) emf2iib(j,k,2) = emf2iib(j,k,1)
           if (emf2iib(j,k,3) .eq. huge) emf2iib(j,k,3) = emf2iib(j,k,2)
           if (emf3iib(j,k,2) .eq. huge) emf3iib(j,k,2) = emf3iib(j,k,1)
           if (emf3iib(j,k,3) .eq. huge) emf3iib(j,k,3) = emf3iib(j,k,2)
          endif ! xmhd
          if(lrad .ne. 0) then
           if (  eriib(j,k,2) .eq. huge)   eriib(j,k,2) =   eriib(j,k,1)
          endif
           if (   doib(j,k,2) .eq. huge)    doib(j,k,2) =    doib(j,k,1)
           if (   doib(j,k,3) .eq. huge)    doib(j,k,3) =    doib(j,k,2)
          if(xiso .eqv. .false.) then
           if (   eoib(j,k,2) .eq. huge)    eoib(j,k,2) =    eoib(j,k,1)
          endif
           if (  v1oib(j,k,2) .eq. huge)   v1oib(j,k,2) =   v1oib(j,k,1)
           if (  v2oib(j,k,2) .eq. huge)   v2oib(j,k,2) =   v2oib(j,k,1)
           if (  v3oib(j,k,2) .eq. huge)   v3oib(j,k,2) =   v3oib(j,k,1)
          if(xmhd) then
           if (  b1oib(j,k,2) .eq. huge)   b1oib(j,k,2) =   b1oib(j,k,1)
           if (  b2oib(j,k,2) .eq. huge)   b2oib(j,k,2) =   b2oib(j,k,1)
           if (  b3oib(j,k,2) .eq. huge)   b3oib(j,k,2) =   b3oib(j,k,1)
           if (emf1oib(j,k,2) .eq. huge) emf1oib(j,k,2) = emf1oib(j,k,1)
           if (emf1oib(j,k,3) .eq. huge) emf1oib(j,k,3) = emf1oib(j,k,2)
           if (emf2oib(j,k,2) .eq. huge) emf2oib(j,k,2) = emf2oib(j,k,1)
           if (emf2oib(j,k,3) .eq. huge) emf2oib(j,k,3) = emf2oib(j,k,2)
           if (emf3oib(j,k,2) .eq. huge) emf3oib(j,k,2) = emf3oib(j,k,1)
           if (emf3oib(j,k,3) .eq. huge) emf3oib(j,k,3) = emf3oib(j,k,2)
          endif ! xmhd
          if(lrad .ne. 0) then
           if (  eroib(j,k,2) .eq. huge)   eroib(j,k,2) =   eroib(j,k,1)
          endif
180      continue
190    continue
c
       do 210 k=ks-2,ke+3
         do 200 i=is-2,ie+3
           if (   dijb(i,k,2) .eq. huge)    dijb(i,k,2) =    dijb(i,k,1)
           if (   dijb(i,k,3) .eq. huge)    dijb(i,k,3) =    dijb(i,k,2)
          if(xiso .eqv. .false.) then
           if (   eijb(i,k,2) .eq. huge)    eijb(i,k,2) =    eijb(i,k,1)
          endif
           if (  v1ijb(i,k,2) .eq. huge)   v1ijb(i,k,2) =   v1ijb(i,k,1)
           if (  v2ijb(i,k,2) .eq. huge)   v2ijb(i,k,2) =   v2ijb(i,k,1)
           if (  v3ijb(i,k,2) .eq. huge)   v3ijb(i,k,2) =   v3ijb(i,k,1)
          if(xmhd) then
           if (  b1ijb(i,k,2) .eq. huge)   b1ijb(i,k,2) =   b1ijb(i,k,1)
           if (  b2ijb(i,k,2) .eq. huge)   b2ijb(i,k,2) =   b2ijb(i,k,1)
           if (  b3ijb(i,k,2) .eq. huge)   b3ijb(i,k,2) =   b3ijb(i,k,1)
           if (emf1ijb(i,k,2) .eq. huge) emf1ijb(i,k,2) = emf1ijb(i,k,1)
           if (emf1ijb(i,k,3) .eq. huge) emf1ijb(i,k,3) = emf1ijb(i,k,2)
           if (emf2ijb(i,k,2) .eq. huge) emf2ijb(i,k,2) = emf2ijb(i,k,1)
           if (emf2ijb(i,k,3) .eq. huge) emf2ijb(i,k,3) = emf2ijb(i,k,2)
           if (emf3ijb(i,k,2) .eq. huge) emf3ijb(i,k,2) = emf3ijb(i,k,1)
           if (emf3ijb(i,k,3) .eq. huge) emf3ijb(i,k,3) = emf3ijb(i,k,2)
          endif ! xmhd
          if(lrad .ne. 0) then
           if (  erijb(i,k,2) .eq. huge)   erijb(i,k,2) =   erijb(i,k,1)
          endif
           if (   dojb(i,k,2) .eq. huge)    dojb(i,k,2) =    dojb(i,k,1)
           if (   dojb(i,k,3) .eq. huge)    dojb(i,k,3) =    dojb(i,k,2)
          if(xiso .eqv. .false.) then
           if (   eojb(i,k,2) .eq. huge)    eojb(i,k,2) =    eojb(i,k,1)
          endif
           if (  v1ojb(i,k,2) .eq. huge)   v1ojb(i,k,2) =   v1ojb(i,k,1)
           if (  v2ojb(i,k,2) .eq. huge)   v2ojb(i,k,2) =   v2ojb(i,k,1)
           if (  v3ojb(i,k,2) .eq. huge)   v3ojb(i,k,2) =   v3ojb(i,k,1)
          if(xmhd) then
           if (  b1ojb(i,k,2) .eq. huge)   b1ojb(i,k,2) =   b1ojb(i,k,1)
           if (  b2ojb(i,k,2) .eq. huge)   b2ojb(i,k,2) =   b2ojb(i,k,1)
           if (  b3ojb(i,k,2) .eq. huge)   b3ojb(i,k,2) =   b3ojb(i,k,1)
           if (emf1ojb(i,k,2) .eq. huge) emf1ojb(i,k,2) = emf1ojb(i,k,1)
           if (emf1ojb(i,k,3) .eq. huge) emf1ojb(i,k,3) = emf1ojb(i,k,2)
           if (emf2ojb(i,k,2) .eq. huge) emf2ojb(i,k,2) = emf2ojb(i,k,1)
           if (emf2ojb(i,k,3) .eq. huge) emf2ojb(i,k,3) = emf2ojb(i,k,2)
           if (emf3ojb(i,k,2) .eq. huge) emf3ojb(i,k,2) = emf3ojb(i,k,1)
           if (emf3ojb(i,k,3) .eq. huge) emf3ojb(i,k,3) = emf3ojb(i,k,2)
          endif ! xmhd
          if(lrad .ne. 0) then
           if (  erojb(i,k,2) .eq. huge)   erojb(i,k,2) =   erojb(i,k,1)
          endif
200      continue
210    continue
c
       do 230 j=js-2,je+3
         do 220 i=is-2,ie+3
           if (   dikb(i,j,2) .eq. huge)    dikb(i,j,2) =    dikb(i,j,1)
           if (   dikb(i,j,3) .eq. huge)    dikb(i,j,3) =    dikb(i,j,2)
          if(xiso .eqv. .false.) then
           if (   eikb(i,j,2) .eq. huge)    eikb(i,j,2) =    eikb(i,j,1)
          endif
           if (  v1ikb(i,j,2) .eq. huge)   v1ikb(i,j,2) =   v1ikb(i,j,1)
           if (  v2ikb(i,j,2) .eq. huge)   v2ikb(i,j,2) =   v2ikb(i,j,1)
           if (  v3ikb(i,j,2) .eq. huge)   v3ikb(i,j,2) =   v3ikb(i,j,1)
          if(xmhd) then
           if (  b1ikb(i,j,2) .eq. huge)   b1ikb(i,j,2) =   b1ikb(i,j,1)
           if (  b2ikb(i,j,2) .eq. huge)   b2ikb(i,j,2) =   b2ikb(i,j,1)
           if (  b3ikb(i,j,2) .eq. huge)   b3ikb(i,j,2) =   b3ikb(i,j,1)
           if (emf1ikb(i,j,2) .eq. huge) emf1ikb(i,j,2) = emf1ikb(i,j,1)
           if (emf1ikb(i,j,3) .eq. huge) emf1ikb(i,j,3) = emf1ikb(i,j,2)
           if (emf2ikb(i,j,2) .eq. huge) emf2ikb(i,j,2) = emf2ikb(i,j,1)
           if (emf2ikb(i,j,3) .eq. huge) emf2ikb(i,j,3) = emf2ikb(i,j,2)
           if (emf3ikb(i,j,2) .eq. huge) emf3ikb(i,j,2) = emf3ikb(i,j,1)
           if (emf3ikb(i,j,3) .eq. huge) emf3ikb(i,j,3) = emf3ikb(i,j,2)
          endif ! xmhd
          if(lrad .ne. 0) then
           if (  erikb(i,j,2) .eq. huge)   erikb(i,j,2) =   erikb(i,j,1)
          endif
           if (   dokb(i,j,2) .eq. huge)    dokb(i,j,2) =    dokb(i,j,1)
           if (   dokb(i,j,3) .eq. huge)    dokb(i,j,3) =    dokb(i,j,2)
          if(xiso .eqv. .false.) then
           if (   eokb(i,j,2) .eq. huge)    eokb(i,j,2) =    eokb(i,j,1)
          endif
           if (  v1okb(i,j,2) .eq. huge)   v1okb(i,j,2) =   v1okb(i,j,1)
           if (  v2okb(i,j,2) .eq. huge)   v2okb(i,j,2) =   v2okb(i,j,1)
           if (  v3okb(i,j,2) .eq. huge)   v3okb(i,j,2) =   v3okb(i,j,1)
          if(xmhd) then
           if (  b1okb(i,j,2) .eq. huge)   b1okb(i,j,2) =   b1okb(i,j,1)
           if (  b2okb(i,j,2) .eq. huge)   b2okb(i,j,2) =   b2okb(i,j,1)
           if (  b3okb(i,j,2) .eq. huge)   b3okb(i,j,2) =   b3okb(i,j,1)
           if (emf1okb(i,j,2) .eq. huge) emf1okb(i,j,2) = emf1okb(i,j,1)
           if (emf1okb(i,j,3) .eq. huge) emf1okb(i,j,3) = emf1okb(i,j,2)
           if (emf2okb(i,j,2) .eq. huge) emf2okb(i,j,2) = emf2okb(i,j,1)
           if (emf2okb(i,j,3) .eq. huge) emf2okb(i,j,3) = emf2okb(i,j,2)
           if (emf3okb(i,j,2) .eq. huge) emf3okb(i,j,2) = emf3okb(i,j,1)
           if (emf3okb(i,j,3) .eq. huge) emf3okb(i,j,3) = emf3okb(i,j,2)
          endif ! xmhd
          if(lrad .ne. 0) then
           if (  erokb(i,j,2) .eq. huge)   erokb(i,j,2) =   erokb(i,j,1)
          endif
220      continue
230    continue
c
c-------------------------  GRID MOTION CONTROL  -----------------------
c
c  x1fac     x1 motion factor
c            < 0 gives "Lagrangian" tracking in x1 lines
c  x2fac     x2 motion factor
c            < 0 gives "Lagrangian" tracking in x2 lines
c  x3fac     x3 motion factor
c            < 0 gives "Lagrangian" tracking in x3 lines
c  iga       i<ia => zone ratio is preserved in x1 lines
c  jga       j<ja => zone ratio is preserved in x2 lines
c  kga       k<ka => zone ratio is preserved in x3 lines
c  igcon     selects grid treatment:
c            =0 => separate motion
c            =1 => averaged motion
c            =2 => tracking x1, x2, and x3 boundaries
c            =3 => averaged boundary tracking
c            =4 => input grid boundary speeds
c                  vg1(io) = x1fac * central sound speed
c                  vg2(jo) = x2fac * central sound speed
c                  vg3(ko) = x3fac * central sound speed
c
      if (myid_w .eq. 0) then
        x1fac = 0.0
        x2fac = 0.0
        x3fac = 0.0
        iga    = is
        jga    = js
        kga    = ks
        igcon = 0
c
        read (1,gcon)
        write(2,gcon)
c
        iga =  min (  max ( iga, is ), ie )
        jga =  min (  max ( jga, js ), je )
        kga =  min (  max ( kga, ks ), ke )
c
      endif
#ifdef MPI_USED
      if (myid_w .eq. 0) then
         buf_in(1) = x1fac
         buf_in(2) = x2fac
         buf_in(3) = x3fac
        ibuf_in(1) = iga
        ibuf_in(2) = jga
        ibuf_in(3) = kga
        ibuf_in(4) = igcon
      endif
      call MPI_BCAST( buf_in,3,MPI_FLOAT  ,0,comm3d,ierr)
      call MPI_BCAST(ibuf_in,4,MPI_INTEGER,0,comm3d,ierr)
      if (myid_w .ne. 0) then
        x1fac =  buf_in(1)
        x2fac =  buf_in(2)
        x3fac =  buf_in(3)
        iga    = ibuf_in(1)
        jga    = ibuf_in(2)
        kga    = ibuf_in(3)
        igcon = ibuf_in(4)
      endif
#endif /* MPI */
C      if(xvgrid) then
       do 320 i=is-2,ie+3
         vg1(i) = 0.0
320    continue
       do 330 j=js-2,je+3
         vg2(j) = 0.0
330    continue
       do 340 k=ks-2,ke+3
         vg3(k) = 0.0
340    continue
C      endif
c
c----------------------  INITIAL TIMESTEP  -----------------------------
c
       nreq = 0
c
c Compute the momentum densities from the initial conditions.  We
c need some density boundary values, but we won't try to be efficient
c here.  Be sure to call bvald in sequence so that the values on edges
c get set properly.  
c
#ifdef MPI_USED
       call MPI_BARRIER(comm3d, ierr)
#endif
       nreq = 0
       nsub = nsub + 1
       call bvald(1,0,0,0,0,0,d)
#ifdef MPI_USED
       if(nreq .ne. 0) call MPI_WAITALL(nreq,req,stat,ierr)
#endif
       nreq = 0
       nsub = nsub + 1
       call bvald(0,0,1,0,0,0,d)
#ifdef MPI_USED
       if(nreq .ne. 0) call MPI_WAITALL(nreq,req,stat,ierr)
#endif
       nreq = 0
       nsub = nsub + 1
       call bvald(0,0,0,0,1,0,d)
#ifdef MPI_USED
       if(nreq .ne. 0) call MPI_WAITALL(nreq,req,stat,ierr)
#endif
       if(ldimen .eq. 3) then
        kone = 1
       else
        kone = 0
       endif
       if(ldimen .gt. 1) then
        jone = 1
       else
        jone = 0
       endif
c
       do 370 k=ks,ke
        do 360 j=js,je
         do 350 i=is,ie
          w3da(i,j,k) = v1(i,j,k)*0.5*(d(i-1,j     ,k     )+d(i,j,k))
          w3db(i,j,k) = v2(i,j,k)*0.5*(d(i  ,j-jone,k     )+d(i,j,k))
     1                * g2b(i)
          w3dc(i,j,k) = v3(i,j,k)*0.5*(d(i  ,j     ,k-kone)+d(i,j,k))
     1                * g31b(i) * g32b(j)
350      continue
360     continue
370    continue
c
c      initialize radiation scratch array
c
      if(lrad .ne. 0) then
       do k = ks, ke
        do j = js, je
         do i = is, ie
          w3dh(i,j,k) = er(i,j,k)/d(i,j,k)
         enddo
        enddo
       enddo
      endif ! lrad
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
      endif ! nspec
c
c Set the artificial viscosity time step dtqqi2 to zero.
c
       dtqqi2 = 0.0
       dtmin  = tiny
c
c
      if(lrad .ne. 0) then
       dtimrdi2 = 1.0/(dtrat * tlim)**2
       dtnri2   = 1.0 / dt**2
      endif
       call nudt
       dtmin = dtrat * dt
       dt    = 10.0 * dtmin
      if(lrad .ne. 0) then
       dtnri2   = 1.0 / dt**2
       dtrdi2   = 1.0 / dt**2
       dtimrdi2   = 1.0 / dt**2
      endif
c
c-----------------------  INITIALIZE NEW GRID  -------------------------
c
c  Compute grid velocities, and new grid positions [in routine newgrid].
c  Note newgrid will recompute "n+1/2" and "n+1" grid lines, but only
c  if x1fac or x2fac .ne. 0.  The "n+1/2" and "n+1" grid lines are
c  initialized in ggen to the old values in case the grid never moves.  Note
c  newgrid must be called after nudt since it needs the timestep.  Thus,
c  the initial timestep computed above did not account for grid motion.
c
      if(xvgrid) call newgrid
c
c------------------------  Initialize everything else  -----------------
c
      ix1x2x3 = 1
      jx1x2x3 = 1  !  For rad_pf* routines
c
      return
      end
