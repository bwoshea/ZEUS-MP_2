c-----------------------------------------------------------------------
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
       subroutine auto_h4
c
c ZEUS-MP Post-processor: AUTO_H4.
c
c      Purpose:  this is like h4splice, except that the DO loop
c      over dump number is automated such that all dumps beginning
c      with Nbeg and ending with Nend are processed, where Nbeg
c      and Nend are input by the user.
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
      implicit NONE
c
      integer   :: irestart
c
      integer   :: nbl,igrid,imin,imax,jmin,jmax,kmin,kmax
      logical   :: lgrid
c
      integer   :: iskip,jskip,kskip,izone,jzone,kzone
      integer   :: n, indx, icmb, ndump, nbeg, nend
      integer   :: it,jt,kt,i,j,k,itil
      integer   :: nfunc,lfunc,ret,shape(3),rank,rcmb,shcmb(3)
c
      real(rl4), dimension(:), allocatable :: xscale,yscale,zscale
c
      integer, parameter :: max_pe = 256
      integer, parameter :: max_fc = 257
c
      real(rl4), dimension(:), allocatable :: xscmb,yscmb, zscmb
      real(rl4), dimension(:), allocatable :: data ,dcmb
c
      real(rl) :: x1min,x1max,x1rat,dx1min,dfndx1r,x1r,deltx1r,errx1r,
     &            x2min,x2max,x2rat,dx2min,dfndx2r,x2r,deltx2r,errx2r,
     &            x3min,x3max,x3rat,dx3min,dfndx3r,x3r,deltx3r,errx3r
c
      character*8 phrase
      character*15 hdf(max_pe), cmb
      character*120 line
      character strng(max_fc)*32,qunit*1,qfrmt*1,qcoord*16
c
      integer  dsgdims,dsgdisc,dsgdata,dsgdast
      integer  dssdims,dssdast,dssdisc,dsadata,dspdata
      external dsgdims,dsgdisc,dsgdata,dsgdast
      external dssdims,dssdast,dssdisc,dsadata,dspdata
CPS
      integer ntotal,inz,jnz,knz
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
c.......................................................................
c
c Read ZEUS-MP input file "zmp_inp" to determine the number of zones
c in each direction in the physical mesh (nx1z,nx2z,nx3z), and the 
c number of tiles in each direction "ntiles".
c
      open(1,file='zmp_inp',status='old')
c
c     read the code configuration namelists
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
      allocate(xscmb(max_pe*in))
      allocate(yscmb(max_pe*jn))
      allocate(zscmb(max_pe*kn))
      allocate(data (in*jn*kn))
      allocate(dcmb(max_pe*in*jn*kn))
      allocate(xscale(in))
      allocate(yscale(jn))
      allocate(zscale(kn))
c
c
c   Look for and read the "mpitop" and "ggen*" namelists.
c
       do 10 n=1,100
         read(1,"(a8)") phrase
         if (phrase(3:8) .eq. 'MPITOP' .or.
     .       phrase(3:8) .eq. 'mpitop') goto 20
   10  continue
       write(*,"(/'HSPLICE: ERROR -- Did not find mpitop namelist.')")
       close (1)
       return
c
   20  continue
       backspace (1)
       read(1,mpitop)
       write(*,"('HSPLICE: ntiles is ',3i6)") ntiles(1),ntiles(2)
     & , ntiles(3)
c
c   Read the "rescon" namelist to find the run's "id".
c
       irestart = 0
        tdump = 0.0
       dtdump = 0.0
       id     = 'aa'
       resfile= 'resaa000000.000'
       read (1,rescon)
       write(*,"(/'HSPLICE: id is ',a)") id
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
       write(*,"('HSPLICE: nx1z,nx2z,nx3z are ',3i6)") nx1z,nx2z,nx3z
       close(1)
c.......................................................................
c
c Ask for the dump number.  Loop back to this point when finished
c with this dump number.  Negative input dump number terminates
c execution.
c
   60  continue
c
      write(*,"(/'AUTO_HSPLICE: Enter beginning and ending dump',
     .            ' numbers'/)")
c
      read(*,*) nbeg, nend
c
c Generate names of hdf files to read.
c
      DO NDUMP = nbeg, nend

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
c.......................................................................
c
c Read each function from each tile's dump and write it into
c the combined file 'cmb' ("hdf<id>.<ndump>").
c
       write(cmb,"(a3,a2,'.',i3.3)") 'hdf',id,ndump
c
c   Loop over the number of functions to be read -- determine number.
c
       nfunc = 5
       if(xmhd) nfunc = nfunc + 3
#ifdef VORTEX
       nfunc = nfunc + 1
#endif
       if(lrad .ne. 0) nfunc = nfunc + 1
c
       do lfunc=1,nfunc
c
c   For each tile, read one function.  On the first pass, read
c   and save the mesh points for each tile.
c
CPS
         indx = 0
         inz=0
         jnz=0
         knz=0
         do kt=0,ntiles(3)-1
           do jt=0,ntiles(2)-1
             do it=0,ntiles(1)-1
               indx = indx + 1
               if (lfunc.eq.1) then
C
                 ret = dsgdims(hdf(indx),rank,shape,3)
                 ret = dsgdisc(1,shape(1),xscale)
                 ret = dsgdisc(2,shape(2),yscale)
                 ret = dsgdisc(3,shape(3),zscale)
c
c         Copy the mesh points for the tile into the cmb mesh arrays.
c
CPS
                 do i=1,nx1z
                   xscmb(i+it*nx1z) = xscale(i)
                 enddo ! i
                 do j=1,nx2z
                   yscmb(j+jt*nx2z) = yscale(j)
                 enddo ! j
                 do k=1,nx3z
                   zscmb(k+kt*nx3z) = zscale(k)
                 enddo ! k
               endif ! lfunc = 1
C
c
c         Read the tile HDF file for the lfunc-th function.
c
               do n=1,lfunc
                 ret = dsgdata(hdf(indx),rank,shape,data)
                 ret = dsgdast(strng(lfunc),qunit,qfrmt,qcoord)
               enddo ! n
c
c         Copy data values for this tile into the cmb data array.
c
               do k=1,nx3z
                 do j=1,nx2z
                   do i=1,nx1z
                     itil       = i + (j-1)*nx1z + (k-1)*nx1z*nx2z
CPS
                     icmb       = i + it*nx1z 
     1                          + (j+jt*nx2z-1)*(ntiles(1)*nx1z)
     2                          + (k+kt*nx3z-1)*(ntiles(1)*nx1z)
     3                          * (ntiles(2)*nx2z)
C
                     dcmb(icmb) = data(itil)
                   enddo ! i
                 enddo ! j
               enddo ! k
c
CPS
             enddo ! it
           enddo ! jt
         enddo ! kt
C
c
c   For the first function only, name the cmb mesh HDF file and
c   prepare it to receive data (after reading all tiles mesh points).
c
         rcmb     = 3
CPS
         shcmb(1) = nx1z*ntiles(1)
         shcmb(2) = nx2z*ntiles(2)
         shcmb(3) = nx3z*ntiles(3)
C
         ret = dssdims(rcmb,shcmb)
         ret = dssdisc(1,shcmb(1),xscmb)
         ret = dssdisc(2,shcmb(2),yscmb)
         ret = dssdisc(3,shcmb(3),zscmb)
#ifndef SD
         ret = dssdast(strng(lfunc),' ',' ',qcoord)
#endif
         if (lfunc.eq.1) then
c
c   Write the cmb mesh data into the combined file for this function.
c
           ret = dspdata(cmb,rcmb,shcmb,dcmb)
         else ! lfunc = 1
c
c   Append the cmb mesh data into the combined file for this function.
c
           ret = dsadata(cmb,rcmb,shcmb,dcmb)
         endif ! lfunc = 1
       enddo ! lfunc
c
      ENDDO ! NDUMP
c
      return
c.......................................................................
c
c Error messages
c
  999  continue
       write(*,"(/'HSPLICE: ERROR -- could not open file ',a15
     & , ' for tile ',i4)") hdf(indx),indx
c
       return
       end
