c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////               O P T I O N S               \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c=======================================================================
      subroutine options
c
      use param
      use config
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      character*39 :: GEOMETRY_OPTION
      character*39 :: GRID_OPTION      
      character*39 :: PT_MASS_OPTION   
      character*39 :: GRAVITY_OPTION   
      character*39 :: PS_OPTION      
      character*39 :: MHD_OPTION      
      character*39 :: RAD_OPTION   
      character*39 :: AV_OPTION   
      character*39 :: MPI_OPTION      
      character*39 :: RES_DUMPS_OPTION 
      character*39 :: HDF_DUMP_OPTION 
      character*39 :: TSL_DUMPS_OPTION
      character*39 :: TXT_DUMPS_OPTION
c
      if(lgeom .eq. 1) then
      GEOMETRY_OPTION  = "              * geometry:        XYZ   "
      endif
      if(lgeom .eq. 2) then
      GEOMETRY_OPTION  = "              * geometry:        ZRP   "
      endif
      if(lgeom .eq. 3) then
      GEOMETRY_OPTION  = "              * geometry:        RTP   "
      endif
c
      if(xvgrid) then
      GRID_OPTION      = "              * moving grid      ON    "
      else
      GRID_OPTION      = "              * moving grid      OFF   "
      endif
c
      if(xptmass) then
      PT_MASS_OPTION   = "              * point masses     ON    "
      else
      PT_MASS_OPTION   = "              * point masses     OFF   "
      endif
c
      if(xgrav .or. xgrvfft) then
      GRAVITY_OPTION   = "              * self-gravity     ON    "
      if(xgrvfft) then
       PS_OPTION       = "              * Poisson solver: FFTW   "
      else
       PS_OPTION       = "              * Poisson solver: MGMPI  "
      endif
      else
      GRAVITY_OPTION   = "              * self-gravity     OFF   "
      endif
c
      if(xmhd) then
      MHD_OPTION       = "              * magnetic fields  ON    "
      else
      MHD_OPTION       = "              * magnetic fields  OFF   "
      endif
c
      if(lrad .eq. 0) then
      RAD_OPTION       = "              * rad transport    OFF   "
      endif
      if(lrad .eq. 1) then
      RAD_OPTION       = "              * rad transport = greyFLD"
      endif
      if(lrad .eq. 2) then
      RAD_OPTION       = "              * rad transport = MGFLD  "
      endif
      if(lrad .eq. 3) then
      RAD_OPTION       = "              * rad transport = MGVTEF "
      endif
c
      if(xsubav) then
       AV_OPTION       = "              * A.V. sub-cycling ON    "
      else
       AV_OPTION       = "              * A.V. sub-cycling OFF   "
      endif
c
#ifdef MPI_USED
      MPI_OPTION       = "              * message passing  ON    "
#else
      MPI_OPTION       = "              * message passing  OFF   "
#endif
c
      if(xrestart) then
      RES_DUMPS_OPTION = "              * restart dumps    ON    "
      else
      RES_DUMPS_OPTION = "              * restart dumps    OFF   "
      endif
c
      if(xhdf) then
#ifdef USE_HDF4
      HDF_DUMP_OPTION =  "              * HDF4 VIZ dumps   ON    "
#endif
#ifdef USE_HDF5
      HDF_DUMP_OPTION =  "              * HDF5 VIZ dumps   ON    "
#endif
      else
      HDF_DUMP_OPTION =  "              * HDF VIZ dumps    OFF   "
      endif
c
      if(xtsl) then
      TSL_DUMPS_OPTION = "              * TSL dumps        ON    "
      else
      TSL_DUMPS_OPTION = "              * TSL dumps        OFF   "
      endif
c
      if(xascii) then
      TXT_DUMPS_OPTION = "              * text dumps       ON    "
      else
      TXT_DUMPS_OPTION = "              * text dumps       OFF   "
      endif
c
c  write greeting.
c
      write(6,"(///10x,'ZZZZZ EEEEE U   U  SSSS     M   M PPPP ')")
      write(6,   "(10x,'   Z  E     U   U S         MM MM P   P')")
      write(6,   "(10x,'  Z   EEEE  U   U  SSS  === M M M PPPP ')")
      write(6,   "(10x,' Z    E     U   U     S     M   M P    ')")
      write(6,   "(10x,'ZZZZZ EEEEE  UUU  SSSS      M   M P    ')")
      write(6,"()")
      write(6,"()")
      write(6,   "(10x,'    DEVELOPED BY ROBERT A. FIEDLER')")
      write(6,   "(10x,'      ZEUS-MP V2.1.2 - 1/25/07')")
      write(6,"()")
      write(6,   "(10x,'    RUNNING AS ',  i4,' PROCESS(ES)')")nprocs_w
      write(6,   "(10x,'   WITH THE FOLLOWING CONFIGURATION:')")
      write(6,"()")
      write(6,1) GEOMETRY_OPTION
      write(6,1) GRID_OPTION  
      write(6,1) PT_MASS_OPTION   
      write(6,1) GRAVITY_OPTION   
      if(xgrvfft) then
       write(6,1) PS_OPTION   
      endif
      write(6,1) MHD_OPTION       
      write(6,1) RAD_OPTION 
      write(6,1) MPI_OPTION 
      write(6,1) RES_DUMPS_OPTION 
      write(6,1) HDF_DUMP_OPTION 
      write(6,1) TSL_DUMPS_OPTION 
      write(6,1) TXT_DUMPS_OPTION 
      write(6,"()")
c
1     format(a39)
c
      return
      end
