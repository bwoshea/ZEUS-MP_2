c-----------------------------------------------------------------------
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
       subroutine auto_h5
#ifdef USE_HDF5
c
c ZEUS-MP Post-processor: AUTO_H5
c
c      Purpose:  this is like h5splice, except that the DO loop
c      over dump number is automated such that all dumps beginning
c      with Nbeg and ending with Nend are processed, where Nbeg
c      and Nend are input by the user.
c
c Written:  By John Hayes, May 2006
c
c.......................................................................
c
      use real_prec
      use config
      use param
      use mpino
      use mpipar
      use grid
      use root
c
      use hdf5
c
      implicit NONE
c
c-----------------------------------------------------------------------
c     hdf5-specific parameters, identifiers, etc.
c-----------------------------------------------------------------------
c
      integer        :: rank, error
      integer(hid_t) :: tfile_id             ! tile file IDs
      integer(hid_t) :: cfile_id             ! combined file ID
      integer(hsize_t), dimension(7) :: dims
c
      character*7    :: ax_label(0:3), scr_ch, scr_1
      character*11, dimension(:), allocatable :: dsetname 
c
c-----------------------------------------------------------------------
c     zmp_inp namelists/parameters
c-----------------------------------------------------------------------
c
      integer   :: irestart
c
      real(rl) :: x1min,x1max,x1rat,dx1min,dfndx1r,x1r,deltx1r,errx1r,
     &            x2min,x2max,x2rat,dx2min,dfndx2r,x2r,deltx2r,errx2r,
     &            x3min,x3max,x3rat,dx3min,dfndx3r,x3r,deltx3r,errx3r
c
      integer   :: nbl,igrid,imin,imax,jmin,jmax,kmin,kmax
      logical   :: lgrid
c
      namelist /geomconf/  lgeom, ldimen
      namelist /physconf/  lrad   , xhydro , xgrav, xmhd    , xgrvfft,
     .                     xptmass, xtotnrg, xiso , xvgrid  , xsubav,
     .                     xforce , xsphgrv, leos , nspec
      namelist /ioconf/    xascii , xhdf, xrestart, xtsl
      namelist /preconf/   small_no, large_no
      namelist /arrayconf/ izones, jzones, kzones, maxijk
      namelist /rescon/ irestart,tdump,dtdump,id,resfile
      namelist /mpitop/ ntiles, periodic
      namelist /ggen1/  nbl,x1min,x1max,igrid,x1rat,dx1min,lgrid
      namelist /ggen2/  nbl,x2min,x2max,igrid,x2rat,dx2min,lgrid
      namelist /ggen3/  nbl,x3min,x3max,igrid,x3rat,dx3min,lgrid
c
c-----------------------------------------------------------------------
c     file, function, array counters
c-----------------------------------------------------------------------
c
      integer   :: iskip,jskip,kskip,izone,jzone,kzone
      integer   :: n, indx, icmb, ndump, jcmb, kcmb
      integer   :: it,jt,kt,i,j,k,itil
      integer   :: nfunc,lfunc, nbeg, nend
c
c-----------------------------------------------------------------------
c     local/global data arrays
c-----------------------------------------------------------------------
c
      real(rl4)  :: tval
      real(rl4), dimension(:), allocatable :: xscale, yscale, zscale,
     .                                        scr_in, scr_jn, scr_kn
c
      real(rl4), dimension(:), allocatable :: xscmb, yscmb, zscmb
      real(rl4), dimension(:    ), allocatable :: data
      real(rl4), dimension(:,:,:), allocatable :: dcmb
c
      integer, parameter :: max_pe = 256
      integer, parameter :: max_fc = 257
c
c-----------------------------------------------------------------------
c     file names and misc.
c-----------------------------------------------------------------------
c
      character*8 phrase
      character*15 hdf(max_pe), cmb
      character*120 line
c
c=======================================================================
c     Read namelist data
c=======================================================================
c
c-----------------------------------------------------------------------
c     Read ZEUS-MP input file "zmp_inp" to determine the number of zones
c     in each direction in the physical mesh (nx1z,nx2z,nx3z), and the 
c     number of tiles in each direction "ntiles".
c-----------------------------------------------------------------------
c
      open(1,file='zmp_inp',status='old')
c
c-----------------------------------------------------------------------
c     read the code configuration namelists
c-----------------------------------------------------------------------
c
      read(1,geomconf)
      read(1,physconf)
      read(1,ioconf)
      read(1,preconf)
      read(1,arrayconf)
c
      in = izones + 5
      jn = jzones + 5
      kn = kzones + 5
c
      allocate(data (izones*
     .               jzones*
     .               kzones))
c
      allocate(xscale(izones))
      allocate(scr_in(izones))
      allocate(yscale(jzones))
      allocate(scr_jn(jzones))
      allocate(zscale(kzones))
      allocate(scr_kn(kzones))
c
c-----------------------------------------------------------------------
c   Look for and read the "mpitop" and "ggen*" namelists.
c-----------------------------------------------------------------------
c
       do 10 n=1,100
         read(1,"(a8)") phrase
         if (phrase(3:8) .eq. 'mpitop' .or. 
     .       phrase(3:8) .eq. 'MPITOP') goto 20
   10  continue
       write(*,"(/'H5SPLICE: ERROR -- Did not find mpitop namelist.')")
       close (1)
       return
c
   20  continue
       backspace (1)
       read(1,mpitop)
       write(*,"('H5SPLICE: ntiles is ',3i6)") ntiles(1),ntiles(2)
     & , ntiles(3)
c
      allocate(xscmb(izones*ntiles(1)))
      allocate(yscmb(jzones*ntiles(2)))
      allocate(zscmb(kzones*ntiles(3)))
c
      allocate(dcmb (izones*ntiles(1),
     .               jzones*ntiles(2),
     .               kzones*ntiles(3)))
c
c-----------------------------------------------------------------------
c   Read the "rescon" namelist to find the run's "id".
c-----------------------------------------------------------------------
c
       irestart = 0
        tdump = 0.0
       dtdump = 0.0
       id     = 'aa'
       resfile= 'resaa000000.000'
       read (1,rescon)
       write(*,"(/'H5SPLICE: id is ',a)") id
c
       imax  = 3
       imin = imax
       nbl   = 1
       lgrid = .false.
c
   30  continue
       read(1,ggen1)
       imax = imax + nbl
       if (.not. lgrid) go to 30
c
       nx1z = ( imax - imin ) / ntiles(1)
CPS
c
       jmax  = 3
       jmin  = jmax
       nbl   = 1
       lgrid = .false.
c
   40  continue
       read(1,ggen2)
       jmax = jmax + nbl
       if (.not. lgrid) go to 40
c
       nx2z = ( jmax - jmin ) / ntiles(2)
C
c
       kmax  = 3
       kmin = kmax
       nbl   = 1
       lgrid = .false.
c
   50  continue
       read(1,ggen3)
       kmax = kmax + nbl
       if (.not. lgrid) go to 40
c
       nx3z = ( kmax - kmin ) / ntiles(3)
CPS
C
       write(*,"('H5SPLICE: nx1z,nx2z,nx3z are ',3i6)") nx1z,nx2z,nx3z
       close(1)
c
      nfunc = 5
      if(xmhd) nfunc = nfunc + 3
#ifdef VORTEX
      nfunc = nfunc + 1 ! including dump of div(v)
#endif
      if(lrad .ne. 0) nfunc = nfunc + 1
c
      allocate(dsetname(nfunc))
c
      dsetname(1) = " i velocity"
      dsetname(2) = " j velocity"
      dsetname(3) = " k velocity"
      if(xmhd) then
       dsetname(4) = "i mag field"
       dsetname(5) = "j mag field"
       dsetname(6) = "k mag field"
       dsetname(7) = "gas density"
       dsetname(8) = " gas energy"
       if(lrad .gt. 0) dsetname(9) = " rad energy"
      else
       dsetname(4) = "gas density"
       dsetname(5) = " gas energy"
       if(lrad .gt. 0) dsetname(6) = " rad energy"
      endif
c
c
c=======================================================================
c     Initialize FORTRAN interface.
c=======================================================================
c
      CALL h5open_f (error)
c
c=======================================================================
c     Main loop over output dump number
c=======================================================================
c
c-----------------------------------------------------------------------
c     Ask for the dump number.  Loop back to this point when finished
c     with this dump number.  Negative input dump number terminates
c     execution.
c-----------------------------------------------------------------------
c
60    CONTINUE
c
      write(*,"(/'AUTO_H5: Enter beginning and ending dump',
     .            ' numbers'/)")
c
      read(*,*) nbeg, nend
c
c Generate names of hdf files to read.
c
      DO NDUMP = nbeg, nend
c
c-----------------------------------------------------------------------
c     Generate names of hdf files to read.
c-----------------------------------------------------------------------
c
       indx = 0
       do 90 kt=0,ntiles(3)-1
         do 80 jt=0,ntiles(2)-1
           do 70 it=0,ntiles(1)-1
             indx = indx + 1
             write(hdf(indx),"(a3,a2,3i2.2,'.',i3.3)") 'hdf',id
     &       , it,jt,kt,ndump
   70      continue
   80    continue
   90  continue
       nprocs = indx
c
c-----------------------------------------------------------------------
c     Read each function from each tile's dump and write it into
c     the combined file 'cmb' ("hdf<id>.<ndump>").
c-----------------------------------------------------------------------
c
      write(cmb,"(a3,a2,'.',i3.3)") 'hdf',id,ndump
c
      CALL h5fcreate_f(cmb, H5F_ACC_TRUNC_F, cfile_id, error)
c
c=======================================================================
c     Begin loop over data arrays (coords, velocities, rho, etc.)
c=======================================================================
c
      DO LFUNC = 1, NFUNC  ! NFUNC counts only the 3-D data arrays.
c
c-----------------------------------------------------------------------
c     Begin loop over tile files
c-----------------------------------------------------------------------
c
       indx = 0
       DO KT = 0,ntiles(3)-1
        DO JT = 0,ntiles(2)-1
         DO IT = 0,ntiles(1)-1
c
          indx = indx + 1
c
          call h5fopen_f(hdf(indx),h5f_acc_rdonly_f,tfile_id,error)
c
c-----------------------------------------------------------------------
c         for lfunc = 1, also read coordinate arrays and save.  when
c         lfunc > 1, read arrays but discard.
c-----------------------------------------------------------------------
c
          RANK      = 1
          dims(1  ) = 1
          dims(2:7) = 0
          if(lfunc .eq. 1 .and. indx .eq. 1) then
           call read_viz(tfile_id,rank,dims,"   time",tval)
          else
           call read_viz(tfile_id,rank,dims,"   time",scr_1)
          endif
c
          dims(1  ) = izones
          dims(2:7) = 0
          if(lfunc .eq. 1) then
           call read_viz(tfile_id,rank,dims,"i coord",xscale)
           do i = 1, nx1z
            xscmb(i+it*nx1z) = xscale(i)
           enddo
          else
           call read_viz(tfile_id,rank,dims,"i coord",scr_in)
          endif
c
          dims(1  ) = jzones
          dims(2:7) = 0
          if(lfunc .eq. 1) then
           call read_viz(tfile_id,rank,dims,"j coord",yscale)
           do j = 1, nx2z
            yscmb(j+jt*nx2z) = yscale(j)
           enddo
          else
           call read_viz(tfile_id,rank,dims,"j coord",scr_jn)
          endif
c
          dims(1  ) = kzones
          dims(2:7) = 0
          if(lfunc .eq. 1) then
           call read_viz(tfile_id,rank,dims,"k coord",zscale)
           do k = 1, nx3z
            zscmb(k+kt*nx3z) = zscale(k)
           enddo
          else
           call read_viz(tfile_id,rank,dims,"k coord",scr_kn)
          endif
c
          RANK = 3
          dims(1  ) = izones
          dims(2  ) = jzones
          dims(3  ) = kzones
c
          call read_viz(tfile_id,rank,dims,dsetname(lfunc),data)
c
          do k=1,nx3z
           do j=1,nx2z
            do i=1,nx1z
             itil = i + (j-1)*nx1z + (k-1)*nx1z*nx2z
C             icmb       = i + it*nx1z 
C     1                      + (j+jt*nx2z-1)*(ntiles(1)*nx1z)
C     2                      + (k+kt*nx3z-1)*(ntiles(1)*nx1z)
C     3                      * (ntiles(2)*nx2z)
             icmb                 = i + it*nx1z
             jcmb                 = j + jt*nx2z
             kcmb                 = k + kt*nx3z
             dcmb(icmb,jcmb,kcmb) = data(itil)
            enddo ! i
           enddo ! j
          enddo ! k
c 
c-----------------------------------------------------------------------
c     Terminate access to the tile file.
c-----------------------------------------------------------------------
c
          CALL h5fclose_f(tfile_id, error)
c
         ENDDO ! IT
        ENDDO ! JT
       ENDDO ! KT
c
       if(lfunc .eq. 1) then
        RANK = 1
        dims(1  ) = 1
        dims(2:7) = 0
        call write_cmb(cfile_id,rank,dims,"   time",tval)
c
        dims(1  ) = izones*ntiles(1)
        dims(2:7) = 0
        call write_cmb(cfile_id,rank,dims,"i coord",xscmb)
c
        dims(1  ) = jzones*ntiles(2)
        dims(2:7) = 0
        call write_cmb(cfile_id,rank,dims,"j coord",yscmb)
c
        dims(1  ) = kzones*ntiles(3)
        dims(2:7) = 0
        call write_cmb(cfile_id,rank,dims,"k coord",zscmb)
       endif ! lfunc
c
       RANK = 3
       dims(1  ) = izones*ntiles(1)
       dims(2  ) = jzones*ntiles(2)
       dims(3  ) = kzones*ntiles(3)
       dims(4:7) = 0
       call write_cmb(cfile_id,rank,dims,dsetname(lfunc),dcmb)
c
      ENDDO ! LFUNC
c
c-----------------------------------------------------------------------
c     Terminate access to the combined data file.
c-----------------------------------------------------------------------
c
      CALL h5fclose_f(cfile_id, error)
c
c=======================================================================
c     Loop back up for the next dump.
c=======================================================================
c
      ENDDO ! NDUMP
c
#endif /* USE_HDF5 */
      return
      end
