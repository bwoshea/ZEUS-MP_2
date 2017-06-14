c=======================================================================
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
      subroutine msave(filename)
c
c  PURPOSE:  Writes [reads] all common block variables to the logical
c  unit 4.  Currently, the common blocks written [read] are:
c            common /gridvarr/ = real    grid     variables
c            common /gridvari/ = integer grid     variables
c            common /fieldr  / = real    field    variables
c            common /fieldi  / = integer field    variables
c            common /bndryr  / = real    boundary variables
c            common /bndryi  / = integer boundary variables
c            common /rootr   / = real    root     variables
c            common /rooti   / = integer root     variables
c            common /gravcomr/ = real    gravity  variables
c            common /gravcomi/ = integer gravity  variables
c
c  The following blocks are NOT written [read]:
c            common /frmlcomr/ = real    formal soln variables
c            common /frmlcomi/ = integer formal soln variables
c            common /mmntcomr/ = real    moment soln variables
c            common /mmntcomi/ = integer moment soln variables
c 
c  EXTERNALS: [none]
c
c  LOCALS:
c
c  MODIFIED: 26 Aug. 1996 by RAF for ZEUS-MP.
c  MODIFIED:  6 Jan. 1999 by JCH
c  MODIFIED:  7 Jan. 1999 by JCH (increased nrootr for root.h;
c                                 makes room for rad t-step var.)
c  MODIFIED: 18 May. 1999 by efh (increased nrootr due to tslice)
c  MODIFIED: December, 2005 by JCH; increased nfieldr for abun array
c-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use grid
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use field
      use bndry
      use root
      use gravmod
      use restart_arrays
c
#ifdef USE_HDF5
      use hdf5
#endif
c
      implicit NONE
c
      character*15 :: filename
#ifdef USE_HDF5
c
c-----------------------------------------------------------------------
c     hdf5-specific parameters, identifiers, etc.
c-----------------------------------------------------------------------
c
      integer        :: rank, error
      integer(hid_t) :: file_id
      integer(hsize_t), dimension(7) :: dims
#endif /* USE_HDF5 */
c
c=======================================================================
c
#ifndef USE_HDF5
      open(unit=4,file=filename,status='unknown',form='unformatted')
#endif
c=======================================================================
c
      mgridr = 40*in + 36*jn + 24*kn
      ngridr = 27*in + 23*jn + 15*kn + 3
      ngridi = 13
      allocate(rlgrdvr(ngridr+mgridr))
      allocate(ntgrdvr(ngridi))
c
      nfieldr                  =           5*in*jn*kn
      if(xmhd        ) nfieldr = nfieldr + 3*in*jn*kn
      if(lrad  .ne. 0) nfieldr = nfieldr +   in*jn*kn
      if(xgrav       ) nfieldr = nfieldr +   in*jn*kn
      if(nspec .gt. 1) nfieldr = nfieldr +   in*jn*kn*nspec
      allocate(rlfldvr(nfieldr))
c
      nbdryr = 60*jn*kn + 60*in*jn + 60*in*kn + 6*nbvar
      nbdryi = 10*jn*kn + 10*in*jn + 10*in*kn + 6*nbvar
     .                  + 18
      allocate(rlbdryvr(nbdryr))
      allocate(ntbdryvr(nbdryi))
c
      call mapout
c
#ifndef USE_HDF5
      write(4)  rl grd  vr , nt grd  vr
     &         ,rl fld  vr 
     &         ,rl bdry vr , nt bdry vr
     &         ,rl rt   vr , nt rt   vr , ch rt   vr
     &         ,rl grv  vr , nt grv  vr
      if (myid .eq. 0)
     &  write(2,"(/1x,'restart dump written at time=',1pe12.5,' cycle='
     & ,i6)") time,nhy
      close(unit=4)
#else /* USE_HDF5 */
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
      RANK = 1
c
      dims(1  ) = ngridr+mgridr
      dims(2:7) = 0
c
      call write_real_res(file_id,rank,dims,"grid_real",rlgrdvr)
c
      dims(1  ) = ngridi
      dims(2:7) = 0
c
      call write_int_res(file_id,rank,dims,"grid_int",ntgrdvr)
c
      dims(1  ) = nfieldr
      dims(2:7) = 0
c
      call write_real_res(file_id,rank,dims,"field_real",rlfldvr)
c
      dims(1  ) = nbdryr
      dims(2:7) = 0
c
      call write_real_res(file_id,rank,dims,"bndry_real",rlbdryvr)
c
      dims(1  ) = nbdryi
      dims(2:7) = 0
c
      call write_int_res(file_id,rank,dims,"bndry_int",ntbdryvr)
c
      dims(1  ) = nrootr
      dims(2:7) = 0
c
      call write_real_res(file_id,rank,dims,"root_real",rlrtvr)
c
      dims(1  ) = nrooti
      dims(2:7) = 0
c
      call write_int_res(file_id,rank,dims,"root_int",ntrtvr)
c
      dims(1  ) = nrootch
      dims(2:7) = 0
c
      call write_chr_res(file_id,rank,dims,"root_chr",chrtvr)
c
      dims(1  ) = ngravr
      dims(2:7) = 0
c
      call write_real_res(file_id,rank,dims,"grav_real",rlgrvvr)
c
      dims(1  ) = ngravi
      dims(2:7) = 0
c
      call write_int_res(file_id,rank,dims,"grav_int",ntgrvvr)
#endif /* USE_HDF5 */
c
      deallocate(rlgrdvr)
      deallocate(ntgrdvr)
      deallocate(rlfldvr)
      deallocate(rlbdryvr)
      deallocate(ntbdryvr)
c
      return
c
c-----------------------------  MGET  ----------------------------------
c
      entry mget(filename)
#ifndef USE_HDF5
      open(unit=4,file=filename,status='old',form='unformatted')
#endif
c
      mgridr = 40*in + 36*jn + 24*kn
      ngridr = 27*in + 23*jn + 15*kn + 3
      ngridi = 13
      allocate(rlgrdvr(ngridr+mgridr))
      allocate(ntgrdvr(ngridi))
c
      nfieldr                  =           5*in*jn*kn
      if(xmhd        ) nfieldr = nfieldr + 3*in*jn*kn
      if(lrad  .ne. 0) nfieldr = nfieldr +   in*jn*kn
      if(xgrav       ) nfieldr = nfieldr +   in*jn*kn
      if(nspec .gt. 1) nfieldr = nfieldr +   in*jn*kn*nspec
      allocate(rlfldvr(nfieldr))
c
      nbdryr = 60*jn*kn + 60*in*jn + 60*in*kn + 6*nbvar
      nbdryi = 10*jn*kn + 10*in*jn + 10*in*kn + 6*nbvar
     .                  + 18
      allocate(rlbdryvr(nbdryr))
      allocate(ntbdryvr(nbdryi))
c
#ifndef USE_HDF5
      read(4)   rl grd  vr , nt grd  vr
     &         ,rl fld  vr
     &         ,rl bdry vr , nt bdry vr
     &         ,rl rt   vr , nt rt   vr , ch rt   vr
     &         ,rl grv  vr , nt grv  vr
      close(unit=4)
#else /* USE_HDF5 */
c
c-----------------------------------------------------------------------
c     Initialize FORTRAN interface.
c-----------------------------------------------------------------------
c
      call h5open_f (error)  ! initialize f90 interface
c
c-----------------------------------------------------------------------
c     Open restart file
c-----------------------------------------------------------------------
c
      call h5fopen_f(filename,h5f_acc_rdonly_f,file_id,error)
c
c-----------------------------------------------------------------------
c     Read datasets
c-----------------------------------------------------------------------
c
      RANK = 1
c
      dims(1  ) = ngridr+mgridr
      dims(2:7) = 0
c
      call read_real_res(file_id,rank,dims,"grid_real",rlgrdvr)
c
      dims(1  ) = ngridi
      dims(2:7) = 0
c
      call read_int_res(file_id,rank,dims,"grid_int",ntgrdvr)
c
      dims(1  ) = nfieldr
      dims(2:7) = 0
c
      call read_real_res(file_id,rank,dims,"field_real",rlfldvr)
c
      dims(1  ) = nbdryr
      dims(2:7) = 0
c
      call read_real_res(file_id,rank,dims,"bndry_real",rlbdryvr)
c
      dims(1  ) = nbdryi
      dims(2:7) = 0
c
      call read_int_res(file_id,rank,dims,"bndry_int",ntbdryvr)
      dims(1  ) = nrootr
      dims(2:7) = 0
c
      call read_real_res(file_id,rank,dims,"root_real",rlrtvr)
c
      dims(1  ) = nrooti
      dims(2:7) = 0
c
      call read_int_res(file_id,rank,dims,"root_int",ntrtvr)
c
      dims(1  ) = nrootch
      dims(2:7) = 0
c     
      call read_chr_res(file_id,rank,dims,"root_chr",chrtvr)
c
      dims(1  ) = ngravr
      dims(2:7) = 0
c
      call read_real_res(file_id,rank,dims,"grav_real",rlgrvvr)
c
      dims(1  ) = ngravi
      dims(2:7) = 0
c
      call read_int_res(file_id,rank,dims,"grav_int",ntgrvvr)
#endif /* USE_HDF5 */
c
      call mapin
c
      deallocate(rlgrdvr)
      deallocate(ntgrdvr)
      deallocate(rlfldvr)
      deallocate(rlbdryvr)
      deallocate(ntbdryvr)
c
      return
      end
#ifdef USE_HDF5
c=======================================================================
c=======================================================================
      subroutine write_real_res(file_id,rank,dims,dsetname,dset)
c
      use hdf5
c
      implicit none
c
      integer(hid_t) :: file_id                ! file identifier
      integer ::   rank                        ! dataset rank
      integer(hsize_t), dimension(7) :: dims   ! dataset dimensions
      character(len=*) :: dsetname             ! dataset name
      real(kind=8) :: dset
c
c
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
c
      call h5screate_simple_f(rank, dims, dspace_id, error)
c
c                      ! Get dset_id for data set
      call h5dcreate_f(file_id,dsetname,h5t_native_double,dspace_id,
     .                 dset_id,error)
c
      call h5dwrite_f(dset_id, h5t_native_double, dset, dims, error)
      call h5dclose_f(dset_id, error) ! end access to the dataset
      call h5sclose_f(dspace_id, error) ! term. access to data space
c
      return
      end
c=======================================================================
c=======================================================================
      subroutine write_int_res(file_id,rank,dims,dsetname,dset)
c
      use hdf5
c
      implicit none
c
      integer(hid_t) :: file_id                ! file identifier
      integer ::   rank                        ! dataset rank
      integer(hsize_t), dimension(7) :: dims   ! dataset dimensions
      character(len=*) :: dsetname             ! dataset name
      integer :: dset
c
c
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
c
      call h5screate_simple_f(rank, dims, dspace_id, error)
c
c                      ! Get dset_id for data set
      call h5dcreate_f(file_id,dsetname,h5t_native_integer,dspace_id,
     .                 dset_id,error)
c
      call h5dwrite_f(dset_id, h5t_native_integer, dset, dims, error)
      call h5dclose_f(dset_id, error) ! end access to the dataset
      call h5sclose_f(dspace_id, error) ! term. access to data space
c
      return
      end
c=======================================================================
c=======================================================================
      subroutine read_real_res(file_id,rank,dims,dsetname,dset)
c
      use hdf5
c
      implicit none
c
      integer(hid_t) :: file_id                ! file identifier
      integer ::   rank                        ! dataset rank
      integer(hsize_t), dimension(7) :: dims   ! dataset dimensions
      character(len=*) :: dsetname             ! dataset name
      real(kind=8) :: dset
c
c
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
c
      call h5screate_simple_f(rank, dims, dspace_id, error)
c
c                      ! Get dset_id for data set
      call h5dopen_f(file_id,dsetname,dset_id,error)
c
      call h5dread_f(dset_id, h5t_native_double, dset, dims, error)
      call h5dclose_f(dset_id, error) ! end access to the dataset
      call h5sclose_f(dspace_id, error) ! term. access to data space
c
      return
      end
c=======================================================================
c=======================================================================
      subroutine read_int_res(file_id,rank,dims,dsetname,dset)
c
      use hdf5
c
      implicit none
c
      integer(hid_t) :: file_id                ! file identifier
      integer ::   rank                        ! dataset rank
      integer(hsize_t), dimension(7) :: dims   ! dataset dimensions
      character(len=*) :: dsetname             ! dataset name
      integer :: dset
c
c
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
c
      call h5screate_simple_f(rank, dims, dspace_id, error)
c
c                      ! Get dset_id for data set
      call h5dopen_f(file_id,dsetname,dset_id,error)
c
      call h5dread_f(dset_id, h5t_native_integer, dset, dims, error)
      call h5dclose_f(dset_id, error) ! end access to the dataset
      call h5sclose_f(dspace_id, error) ! term. access to data space
c
      return
      end
c=======================================================================
c=======================================================================
      subroutine read_chr_res(file_id,rank,dims,dsetname,dset)
c
      use hdf5
c
      implicit none
c
      integer(hid_t) :: file_id                ! file identifier
      integer ::   rank                        ! dataset rank
      integer(hsize_t), dimension(7) :: dims   ! dataset dimensions
      character(len=*) :: dsetname             ! dataset name
      character :: dset
c
c
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
c
      call h5screate_simple_f(rank, dims, dspace_id, error)
c
c                      ! Get dset_id for data set
      call h5dopen_f(file_id,dsetname,dset_id,error)
c
      call h5dread_f(dset_id, h5t_native_character, dset, dims, error)
      call h5dclose_f(dset_id, error) ! end access to the dataset
      call h5sclose_f(dspace_id, error) ! term. access to data space
c
      return
      end
c=======================================================================
c=======================================================================
      subroutine write_chr_res(file_id,rank,dims,dsetname,dset)
c
      use hdf5
c
      implicit none
c
      integer(hid_t) :: file_id                ! file identifier
      integer ::   rank                        ! dataset rank
      integer(hsize_t), dimension(7) :: dims   ! dataset dimensions
      character(len=*) :: dsetname             ! dataset name
      character :: dset
c
c
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer :: error
c
      call h5screate_simple_f(rank, dims, dspace_id, error)
c
c                      ! Get dset_id for data set
      call h5dcreate_f(file_id,dsetname,h5t_native_character,dspace_id,
     .                 dset_id,error)
c
      call h5dwrite_f(dset_id, h5t_native_character, dset, dims, error)
      call h5dclose_f(dset_id, error) ! end access to the dataset
      call h5sclose_f(dspace_id, error) ! term. access to data space
c
      return
      end
#endif /* USE_HDF5 */
