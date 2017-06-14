c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////               S R C S T E P               \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
       subroutine srcstep
c
c    jms:zeus2d.srcstep <------------------------ source step controller
c                                                          october, 1987
c
c    written by: Jim Stone
c    modified 1: June, 1988 by Jim Stone; incorporated into ZEUS2D
c    modified 2: Spring, 1989 by Jim Stone; rewritten
c    modified 3: February, 1990 by David Clarke; incorporated into
c                ZEUS3D
c    modified 4: July, 1990 by David Clarke; because current densities
c                are not needed to compute source terms in the new CT
c                algorithm (MOCCT), workers can be used to store the
c                momenta, thereby saving redundant computations in STV1
c                and STV2.
c    modified 5: June 1992, by David Clarke; added the total energy
c                option originally designed by Byung-IL Jun.
c    modified 6: RAF, 3/27/96, completely rewritten for ZEUS-MP.
c    modified 7: RAF, 1/2/97, added radiation force, radiation
c                diffusion, and radiation-matter coupling terms.
c    modified 8: RAF, 1/22/97, moved modules into driver routines.
c    modified 9: RAF, 2/18/97, added PF and N-R timestep controllers.
c    modified 10: PSLi, 12/30/99, added subcycle of artificial viscosity.
c
c
c  PURPOSE: Controls the update of velocities (v1, v2, v3) and internal
c  energy (e) from source terms in the equation of motion and energy
c  equations respectively.
c
c  LOCAL VARIABLES:      
c    w3da     scratch 1-momentum denisty
c    w3db     scratch 2-momentum denisty
c    w3dc     scratch 3-momentum denisty
c    w3dd     scratch 1-velocity ; scratch density (pdv)
c    w3de     scratch 2-velocity ; scratch e/d     (pdv)
c    w3df     scratch 3-velocity
c    j1       1-current density
c    j2       2-current density
c    j3       3-current density
c
c  EXTERNALS:
c    AVISC_D  ,  FORCES_D  , PDV_D
c
c-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use root
      use field
      use grid
      use bndry
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use cons
c
      implicit NONE
c
      real(rl) :: subdt, etot, etot_glb
c
      integer  :: i, j, k, index, n
c
c-----------------------------------------------------------------------
c     EQUATION OF STATE
c
c     In this release, gamma-law ideal gas and isothermal equations of
c     state are included and supported by the EOS parameter LEOS=1.
c     User-added EOS modules will require LEOS parameters of 2 and 
c     higher.
c-----------------------------------------------------------------------
c
      call eos_d
c
c-----------------------------------------------------------------------
c FORCES
c
c Thermal and MHD Pressure, Gravity, and Rotational Pseudo-Forces
c (including gravitational point mass)
c
c Routine "forces" updates the three velocity components.  The arrays
c v1, v2, and v3 save the old values, while w3dd, w3de, w3df get the 
c new ones.
c
      if(lrad .gt. 0) call opac_d
c
      call forces_d (v1,v2,v3,w3dd,w3de,w3df)
c
c       subroutine forces_d
c     1 (v1old, v2old, v3old, v1new, v2new, v3new)
c
c.......................................................................
c
c ARTIFICIAL VISCOSITY
c
c  Update velocity and e from the artificial viscosity source terms.
c  We need 1 layer of updated boundary data for the velocity.
c  We use just the "m" layer of d and, unless ISO is defined, e -- they
c  should be up to date already, since we used them in "forces" but did
c  not update them.  We don't need e BVs with no "linear" viscosity.
c
c  Arrays w3dd, w3de, w3df save the old velcoity values, while v1, v2, 
c  v3 get the updated ones.
c
c  The artificial viscosity routine must also compute momentum densities
c  w3da, w3db, w3dc from v1, v2, v3 for use in the transport step.
c
CPS
C  Subcycle of artificial viscosity calculation
C
      if(xsubav) then
       index=0
       subdt=dt
       avisc_dt=courno / (SQRT(dtqqi2)+tiny)
       if(dtqqi2.eq.0.0) avisc_dt=dt
#ifdef MPI_USED
          buf_in(1) = avisc_dt
          call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1
     &                      , MPI_2DOUBLE_PRECISION
     &                      , MPI_MINLOC, comm3d, ierr)
          avisc_dt  =   buf_out(1)
#endif
       DO WHILE (subdt .GT. 0.0)
          IF(subdt.GT.avisc_dt) THEN
             subdt=subdt-avisc_dt
          ELSE
             avisc_dt=subdt
             subdt=-1.0
          ENDIF
          index=index+1
          IF(mod(index,2).GT.0) THEN
             call avisc_d (w3dd,w3de,w3df,v1,v2,v3,w3da,w3db,w3dc)
          ELSE
             call avisc_d (v1,v2,v3,w3dd,w3de,w3df,w3da,w3db,w3dc)
          ENDIF
          avisc_dt=courno / (SQRT(dtqqi2)+tiny)
#ifdef MPI_USED
          buf_in(1) = avisc_dt
          call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1
     &                      , MPI_2DOUBLE_PRECISION
     &                      , MPI_MINLOC, comm3d, ierr)
          avisc_dt  =   buf_out(1)
#endif
       ENDDO
C
C  Update velocity arrays if necessary.
C
       IF(mod(index,2).EQ.0) THEN
          DO k=ks,ke
             DO j=js,je
                DO i=is,ie
                   v1(i,j,k)=w3dd(i,j,k)
                   v2(i,j,k)=w3de(i,j,k)
                   v3(i,j,k)=w3df(i,j,k)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
      else ! xsubav
       call avisc_d (w3dd,w3de,w3df,v1,v2,v3,w3da,w3db,w3dc)

      endif ! xsubav
c
c
c      subroutine avisc_d
c     1            (v1old,v2old,v3old,v1new,v2new,v3new,s1,s2,s3)
c
c......................................................................
c     IMPLICIT GREY FLUX-LIMITED DIFFUSION SOLVER
c......................................................................
c
      if(lrad .ne. 0) call rad_solve
c
c......................................................................
c No radiation, not isothermal.
c
c  COMPRESSIONAL WORK TERM (PdV).
c
c  Finally, update the energy with the PdV source term
c  only if the internal energy equation is being solved.
c
c  We need just 1 "p" layer of updated boundary data for v1, v2, and v3,
c  but none for d and e.  
c
c  Routine pdv also saves the density and e/d in (w3dd,w3de) for use
c  in the transport step.  This is why it is being called even when
c  TOTAL_ENERGY is defined.
c
c  NOTE: PDV is not called if radiation transport (lrad > 0) is used
c        because it is assumed that the pdV term is incorporated in
c        the radiation step (as it is in the diffusion solver
c        included here).  A user-supplied radiation option which
c        neglects the pdV step will require an adjustment to the
c        IF condition which follows:
c
      if((xiso .eqv. .false.) .and. (lrad .eq. 0)) then
C       if(xtotnrg) call eos_d
       call pdv_d (w3dd ,w3de )
      endif ! xiso
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////               S R C S T E P               \\\\\\\\\\
c
c=======================================================================
c
