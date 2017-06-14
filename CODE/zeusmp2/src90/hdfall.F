#ifdef USE_HDF4
c=======================================================================
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
      subroutine hdfall(filename)
c
c  PURPOSE: Makes an hdf dump of all the active field variables.  The
c  set of field variables dumped is problem specific (depends on what
c  physics is defined).  Data is written in the Scientific Data Set
c  format to the file zhzXXNNNNNN.MMM.
c  Note that data must be stored column major and contiguously in order
c  to interface correctly to the C hdf routines.  All variables are
c  dumped as zone centered quantities.
c
c  EXTERNALS: HDF library routines
c
c  LOCALS:
c
c  LAST MODIFIED: by JCH; 3/12/97.
c-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use grid
      use field
      use root
      use scratch
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
      character*15 :: filename
      character*16 :: coordsys
      character*32 :: string
c
      integer      :: i,j,k,indx,kp1
c
c real is correct on PVP Crays, but real*4 is required on T3E.
c This is due to the kind of C floats on the two machines according to
c Albert Cheng of the NCSA HDF group.  The varying integer definitions
c I have decided on by trial and error M-MML 17.5.98
c
C---------------------------------------------------------------------
C      use on non-T3E UNICOS systems
C
C      integer rank,shape(3),ret
C      real data(in*jn*kn),xscale(in),yscale(jn),zscale(kn)
C      integer  dssdims,dssdast,dssdisc,dsadata,dspdata
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C     use on T3E
C
C      integer rank,shape(3),ret
C      real*4 data(in*jn*kn),xscale(in),yscale(jn),zscale(kn)
C      integer  dssdims,dssdast,dssdisc,dsadata,dspdata
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C     use on everything else
C
      integer(kind=4) :: rank,shape(3),ret
      real(rl4)       :: data(in*jn*kn),xscale(in),yscale(jn),
     .                   zscale(kn)
      integer(kind=4) :: dssdims,dssdast,dssdisc,dsadata,dspdata
C---------------------------------------------------------------------
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c
      if(lgeom .eq. 3) coordsys = 'spherical polar' // char(0)
      if(lgeom .eq. 2) coordsys = 'cylindrical' // char(0)
      if(lgeom .eq. 1) coordsys = 'cartesian' // char(0)
c
      do 10 i=is,ie
        xscale(i-is+1) = real(x1b(i))
10    continue
      do 20 j=js,je
        yscale(j-js+1) = real(x2b(j))
20    continue
      do 30 k=ks,ke
        zscale(k-ks+1) = real(x3b(k))
30    continue
c
      rank     = 3
      shape(1) = nx1z
      shape(2) = nx2z
      shape(3) = nx3z
      ret = dssdims(rank,shape)
      ret = dssdisc(1,shape(1),xscale)
      ret = dssdisc(2,shape(2),yscale)
      ret = dssdisc(3,shape(3),zscale)
c
c  1-velocity
c
      do 120 k=ks,ke
        do 110 j=js,je
          do 100 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
	    data(indx) = real(0.5*(v1(i,j,k) + v1(i+1,j,k)))
100       continue
110     continue
120   continue
      write(string,"('1-VELOCITY AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dspdata(filename,rank,shape,data)
c
c  2-velocity
c
      do 220 k=ks,ke
        do 210 j=js,je
          do 200 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
	    data(indx)  = real(0.5*(v2(i,j,k) + v2(i,j+1,k)))
200       continue
210     continue
220   continue
      write(string,"('2-VELOCITY AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
c
c  3-velocity
c
      do 320 k=ks,ke
        if(ldimen .eq. 3) then
         kp1 = k+1
        else
         kp1 = k
        endif
        do 310 j=js,je
          do 300 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
	    data(indx)  = real(0.5*(v3(i,j,k) + v3(i,j,kp1)))
300       continue
310     continue
320   continue
      write(string,"('3-VELOCITY AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
c
      if(xmhd) then
c
c  1-magnetic field
c
      do 420 k=ks,ke
        do 410 j=js,je
          do 400 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
            data(indx)  = real(0.5*(b1(i,j,k) + b1(i+1,j,k)))
400       continue
410     continue
420   continue
      write(string,"('1-MAG FIELD AT TIME=',1pe8.2,'    ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
c
c  2-magnetic field
c
      do 520 k=ks,ke
        do 510 j=js,je
          do 500 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
	    data(indx)  = real(0.5*(b2(i,j,k) + b2(i,j+1,k)))
500       continue
510     continue
520   continue
      write(string,"('2-MAG FIELD AT TIME=',1pe8.2,'    ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
c
c  3-magnetic field
c
      do 620 k=ks,ke
        if(ldimen .eq. 3) then
         kp1 = k+1
        else
         kp1 = k
        endif
        do 610 j=js,je
          do 600 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
	    data(indx)  = real(0.5*(b3(i,j,k) + b3(i,j,kp1)))
600       continue
610     continue
620   continue
      write(string,"('3-MAG FIELD AT TIME=',1pe8.2,'    ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
      endif ! xmhd
c
c  density
c
      do 720 k=ks,ke
        do 710 j=js,je
          do 700 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
            data(indx) = real(d(i,j,k))
700       continue
710     continue
720   continue
      write(string,"('DENSITY AT TIME=',1pe8.2,'        ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
c
c  internal energy
c
      do 820 k=ks,ke
        do 810 j=js,je
          do 800 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
            data(indx) = real(e(i,j,k))
800       continue
810     continue
820   continue
      write(string,"('GAS ENERGY AT TIME=',1pe8.2,'   ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
c
c  radiation internal energy
c
      if(lrad .ne. 0) then
      do 1020 k=ks,ke
        do 1010 j=js,je
          do 1000 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
            data(indx) = real(er(i,j,k))
1000      continue
1010    continue
1020  continue
      write(string,"('RADIATION T AT TIME=',1pe8.2,'    ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
      endif ! lrad
c
      return
      end
#elif defined USE_HDF5
c=======================================================================
c=======================================================================
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c     Purpose: File writer for HDF5 viz data files
c
c     Written by: John Hayes, February 2006
c
      subroutine hdfall(filename)
c
      use real_prec
      use config
      use param
      use grid
      use field
      use root
      use scratch
      use cons
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      use hdf5
c
      implicit NONE
c
c-----------------------------------------------------------------------
c     ZEUS-MP - specific file and data descriptors, arrays, etc.
c-----------------------------------------------------------------------
c
      character*15 :: filename
      character*16 :: coordsys
      character*32 :: string
c
      integer      :: i,j,k,indx,kp1
c
      real(rl4)    :: icoord(ie-is+1)
      real(rl4)    :: jcoord(je-js+1)
      real(rl4)    :: kcoord(ke-ks+1)
c
      real(rl4)    :: tval
      real(rl4)    :: data(ie-is+1,je-js+1,ke-ks+1)
c
c-----------------------------------------------------------------------
c     hdf5-specific parameters, identifiers, etc.
c-----------------------------------------------------------------------
c
      integer        :: rank, error
      integer(hid_t) :: file_id
      integer(hsize_t), dimension(7) :: dims
c
c-----------------------------------------------------------------------
c     Initialize FORTRAN interface.
c-----------------------------------------------------------------------
c
      CALL h5open_f (error)
c
c-----------------------------------------------------------------------
c     Create a new file using default properties.
c-----------------------------------------------------------------------
c
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
c
c-----------------------------------------------------------------------
c     Create/Write datasets
c-----------------------------------------------------------------------
c
c --- evolution time
c
      RANK = 1
c
      tval      = real(time)
c
      dims(1  ) = 1
      dims(2:7) = 0
c
      call write_viz(file_id,rank,dims,"   time",tval)
c
c --- Coordinate arrays
c
      do i = is, ie
       icoord(i-is+1) = x1b(i)
      enddo
c
      dims(1  ) = ie-is+1
      dims(2:7) = 0
c
      call write_viz(file_id,rank,dims,"i coord",icoord)
c
      do j = js, je
       jcoord(j-js+1) = x2b(j)
      enddo
c
      dims(1  ) = je-js+1
      dims(2:7) = 0
c
      call write_viz(file_id,rank,dims,"j coord",jcoord)
c
      do k = ks, ke
       kcoord(k-ks+1) = x3b(k)
      enddo
c
      dims(1  ) = ke-ks+1
      dims(2:7) = 0
c
      call write_viz(file_id,rank,dims,"k coord",kcoord)
c
c --- Field arrays
c
      RANK    = 3
      dims(1  ) = ie-is+1
      dims(2  ) = je-js+1
      dims(3  ) = ke-ks+1
      dims(4:7) = 0
c
c --- 1-Velocity
c
      do k=ks,ke
       do j=js,je
        do i=is,ie
         data(i-is+1,j-js+1,k-ks+1) = real(0.5*(v1(i,j,k)+v1(i+1,j,k)))
        enddo
       enddo
      enddo
c
      call write_viz(file_id,rank,dims," i velocity",data)
c
c --- 2-Velocity
c
      do k=ks,ke
       do j=js,je
        do i=is,ie
         data(i-is+1,j-js+1,k-ks+1) = real(0.5*(v2(i,j,k)+v2(i,j+1,k)))
        enddo
       enddo
      enddo
c
      call write_viz(file_id,rank,dims," j velocity",data)
c
c --- 3-Velocity
c
      do k=ks,ke
       if(ldimen .eq. 3) then
        kp1 = k+1
       else
        kp1 = k
       endif
       do j=js,je
        do i=is,ie
         data(i-is+1,j-js+1,k-ks+1) = real(0.5*(v3(i,j,k)+v3(i,j,kp1)))
        enddo
       enddo
      enddo
c
      call write_viz(file_id,rank,dims," k velocity",data)
c
      if(XMHD) then
c
c --- 1-B Field
c
       do k=ks,ke
        do j=js,je
         do i=is,ie
          data(i-is+1,j-js+1,k-ks+1) = real(0.5*(b1(i,j,k)+b1(i+1,j,k)))
         enddo
        enddo
       enddo
c
       call write_viz(file_id,rank,dims,"i mag field",data)
c
c --- 2-B field
c
       do k=ks,ke
        do j=js,je
         do i=is,ie
          data(i-is+1,j-js+1,k-ks+1) = real(0.5*(b2(i,j,k)+b2(i,j+1,k)))
         enddo
        enddo
       enddo
c
       call write_viz(file_id,rank,dims,"j mag field",data)
c
c --- 3-B Field
c
       do k=ks,ke
        if(ldimen .eq. 3) then
         kp1 = k+1
        else
         kp1 = k
        endif
        do j=js,je
         do i=is,ie
          data(i-is+1,j-js+1,k-ks+1) = real(0.5*(b3(i,j,k)+b3(i,j,kp1)))
         enddo
        enddo
       enddo
c
       call write_viz(file_id,rank,dims,"k mag field",data)
      endif ! XMHD
c
c --- Gas density
c
      do k=ks,ke
       do j=js,je
        do i=is,ie
         data(i-is+1,j-js+1,k-ks+1) = real(d(i,j,k))
        enddo
       enddo
      enddo
c
      call write_viz(file_id,rank,dims,"gas density",data)
c
c --- Gas energy
c
      do k=ks,ke
       do j=js,je
        do i=is,ie
         data(i-is+1,j-js+1,k-ks+1) = real(e(i,j,k))
        enddo
       enddo
      enddo
c
      call write_viz(file_id,rank,dims," gas energy",data)
c
c --- radiation energy
c
      if(lrad .gt. 0) then
       do k=ks,ke
        do j=js,je
         do i=is,ie
          data(i-is+1,j-js+1,k-ks+1) = real(er(i,j,k))
         enddo
        enddo
       enddo
c
       call write_viz(file_id,rank,dims," rad energy",data)
      endif ! lrad
c
c-----------------------------------------------------------------------
c     Terminate access to the file.
c-----------------------------------------------------------------------
c
      CALL h5fclose_f(file_id, error)
c
c-----------------------------------------------------------------------
c     Close FORTRAN interface.
c-----------------------------------------------------------------------
c
      CALL h5close_f(error)
c
      return
      end
c=======================================================================
c=======================================================================
      subroutine write_viz(file_id,rank,dims,dsetname,dset)
c
      use hdf5
c
      implicit none
c
      integer(hid_t) :: file_id                ! file identifier
      integer ::   rank                        ! dataset rank
      integer(hsize_t), dimension(7) :: dims   ! dataset dimensions
      character(len=*) :: dsetname             ! dataset name
      real :: dset
c
c
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
c
      call h5screate_simple_f(rank, dims, dspace_id, error)
c
c                      ! Get dset_id for data set
      call h5dcreate_f(file_id,dsetname,h5t_native_real,dspace_id,
     .                 dset_id,error)
c
      call h5dwrite_f(dset_id, h5t_native_real, dset, dims, error)
      call h5dclose_f(dset_id, error) ! end access to the dataset
      call h5sclose_f(dspace_id, error) ! term. access to data space
c
      return
      end
#else
      subroutine hdfall(filename)
      character*15 :: filename
      return
      end
#endif /* USE_HDF4 || USE_HDF5 */
