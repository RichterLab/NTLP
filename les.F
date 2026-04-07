      program les_mpi_2d
c
      use pars
      use fields
      use particles
      use con_data
      use con_stats
      use netcdf_io
#ifdef TECIO
      use tec_io
#endif
      use profiling
      use mod_mpi
      use mod_fft
      use mod_solver
      use mod_init
      use mod_io
      include 'mpif.h'
c
c ------------- definition of internal flags
c
c
c       iDNS    = 0; call the subgrid computation of vis_m and vis_s
c               = 1; call the molecular viscosity and diffusivity
c
c       igrdr   =  3; data comes from restart file
c               =  2; data comes from initialization (random)
c               =  1; data comes from coarser grid (or otherwise)
c
c       ibcu    = 1 ; upper boundary condition set by radiation bc
c               = 0 ; fixed value = 0.
c               = -1; value defined by coarser mesh for all variables
c
c       ibcl    = 0 ; lower boundary condition set by similarity theory (sr. setup)
c               = -1; value defined by coarser mesh for all variables
c
c       ifix_dt = 0 ; variable time step with fixed cfl number in setcon
c               = 1 ; fixed time step set in sr. get_dt
c
c       ifree   = 0 ; use spatially averaged surface conditions for MO (call lower)
c               = 1 ; use point-by-point conditions for MO free convection (call lower_free)
c
c       ihst    = nn; frequency at which global variables are output in history file
c               < 0 ; no history files
c
c
c       ismlt   = 1 ; use businger formulas in MO 
c                 0 ; use large and everyone elses formulas in MO 
c
c       iupwnd  = 0;  use skew symmetric formulas for all derivatives
c                     in scalar equations
c               = 1;  use hybrid upwind scheme for all derivatives
c                     in scalar equations
c
c       ivis0   = 0; old eddy viscosity model 
c               = 1; new eddy viscosity model 
c
c       new_vis = step; the iteration step for which the new model
c                       is turned on when ivis0=1
c               < 0; new model is on at all steps for ivis0=1
c
c       nscl  .ge. 1   number of scalars to be followed set in parameter statements
c                      change entries in sr. init, and sr. suft for surface bc's
c
c -------------------------------------------------------------------------------
c
c ---------- initialize MPI, get myid, numprocs, 
c            test if on root process
c
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
      call mpi_comm_size(mpi_comm_world,numprocs,ierr)
c
      i_root = 0
      l_root = .false.
      if(myid .eq. i_root) l_root = .true.
c
      l_debug = .false.
      if(idebug .eq. 1) l_debug = .true.
c
      ts_mpi = mpi_wtime()

c ---------- initialize the profiling modules
      call initialize_profiling
      call start_phase(measurement_id_solver)

!----- Read the input file for all necessary parameters
      call start_phase(measurement_id_setup)

c
c ------------- establish association between pointers
c               and data structures
c
      call fill_cc
      call fill_cs
      call fill_ci

      call read_input_file

c
c -------- set number of x-y slab cpus
c
!      ncpu_s = 8
c
      case_inp = 'cou'
c
      call get_units
      call gridd
      call setcon

      call change_RH_bcs_to_q
c
c -------------- scratch run
c
      if (irestart.eq.0)  then
         igrdr = 2
         itn = 0
         iti = 0
         it = iti
         case = case_inp
         call init
         call setup
         if (inetcdf .eq. 1) then
         call netcdf_init
         call netcdf_init_histog
         end if
         if (iviznetcdf .eq. 1) then
         call netcdf_init_viz
         end if

         call particle_setup
         call particle_init 
         !call read_part_res !comment out usually 

#ifdef TECIO
         call init_tecio
#endif

c
c ---------- choose routine for getting initial guess
c
         if(iocean .eq. 1) then
            call randoc
         else
            if (ifields .eq. 1) then
            call get_fields
            else
            call random   
            end if
         endif
         !Call dns_vis even when doing LES since vis_m needs to be initialized with something
         call dns_vis
         call get_max
	 
	 !Call to populate the particle statistics with ICs
         call particle_xy_stats
      else

         igrdr = 3
         call restart
         call setup
         iti = it


         !Call dns_vis even when doing LES since vis_m needs to be initialized with something
         call dns_vis
         call get_max

 
         call particle_setup
         call read_part_res

#ifdef TECIO
         call init_tecio
#endif

         if (inetcdf .eq. 1) then
         call netcdf_res
         call netcdf_res_histog
         endif

         if (iviznetcdf .eq. 1) then
         call netcdf_res_viz
         endif

      endif
      call end_phase(measurement_id_setup)

c
c --------------- time loop ------------
c
      tzero = time
      call get_dt
 9000 continue
      call start_phase( measurement_id_timestepping_loop)
      call set_sav(iti)
      if (myid==0) then
      write(*,*) 'Starting time loop'
      write(*,*) 'it,time = ',it,time
      end if

      part => first_particle
      do while (associated(part))

      if (part%pidx == 1 .and. part%procidx == 0) then
      write(*,'(a7,4e15.6)') 'xp1:  ',time,part%xp(1:3)
      write(*,'(a7,4e15.6)') 'vp1:  ',time,part%vp(1:3)
      write(*,'(a7,4e15.6)') 'uf1:  ',time,part%uf(1:3)
      write(*,'(a7,2e15.6)') 'Tp1:  ',time,part%Tp
      write(*,'(a7,2e15.6)') 'Tf1:  ',time,part%Tf
      write(*,'(a7,2e15.6)') 'rad1: ',time,part%radius 
      write(*,'(a7,2e15.6)') 'qinf1:',time,part%qinf 
      write(*,'(a7,2e15.6)') 'qstr1:',time,part%qstar 
      write(*,'(a10,2e15.6)') 'kappa_s1:',time,part%kappa_s
      write(*,'(a7,2e15.6)') 'ms1:',time,part%m_s
      write(*,'(a7,2e15.6)') 'res1:',time,part%res
      write(*,'(a7,e15.6,i)') 'mult1:',time,part%mult
      write(*,'(a7,2e15.6)') 'rc1:',time,part%rc
      end if

      part => part%next
      end do

      if (myid==0) then
      write(*,*) 'time,tnumpart:',time,tnumpart
      write(*,'(a15,3e15.6)') 'radmin,radmax:',time,radmin,radmax
      write(*,'(a15,3e15.6)') 'tempmin,tempmax:',time,tempmin,tempmax
      write(*,'(a15,3e15.6)') 'qmin,qmax:',time,qmin,qmax
      write(*,'(a15,3e15.6)') 'vp1min,vp1max:',time,vp1min,vp1max
      write(*,'(a15,3e15.6)') 'vp2min,vp2max:',time,vp2min,vp2max
      write(*,'(a15,3e15.6)') 'vp2min,vp2max:',time,vp2min,vp2max
      write(*,'(a15,3e15.6)') 'time,radavg:',time,radavg
      write(*,'(a15,e15.6,3i)'),'time,100,impos:',time,tnum100,tnumimpos
      end if



      if(it .ge. new_vis .and. ivis0 .eq. 1) then
          ivis = 1
      else
          ivis = 0
      endif


c
c ---------------- 3 stage runge-kutta time stepping
c
      t_stage_s = mpi_wtime()
      do  8999 istage=1,3
c
      dtzeta = dt*zetas(istage)
      dtgama = dt*gama(istage)
c
c ---------- compute derivatives of (u,v,w)
c

      call start_phase(measurement_id_derivatives)
      call exchange
      call get_derv
      call end_phase(measurement_id_derivatives)

c
c --------- new eddy viscosity, and bcs
c
      call start_phase(measurement_id_eddy_viscosity_and_bcs)
      if(iss .eq. 0 .and. ifree .eq. 0) then
         if (iDNS .eq. 1) then
            call lower_dns
         else
            call lower
         end if
      elseif(ifree .eq. 1) then
         if (iDNS .eq. 1) then
            call lower_dns
         else
            call lower_free
         end if
      endif

      !Fill the base state based on the current surface temp
      call fill_base(tsfcc(1))

      if(ise .eq. numprocs-1) then
         if (iDNS .eq. 1) then
            call upper_dns
         else
            call upper
         end if
      endif
      call bcast_pbc
      call bcast_lbc
      call get_means(istage)
      if(ivis .eq. 1 .and. iDNS .ne. 1) then
         call iso
	 call surfvis
      endif
      if(istage .eq. 1)then
        call xy_stats
        call tke_budget
        do iscl=1,nscl
        call Tvar_budget(iscl)
        enddo
        call extra_flux_terms
        call pbltop(itop)
      endif
      call end_phase(measurement_id_eddy_viscosity_and_bcs)

c
c ------------ save history files
c
      if(mnout .and. istage .eq. 1)  then
          if(l_debug) then
             call print(nprt,izs,ize)
          endif
          if(l_root) call print(6,1,nnz)
      endif
      if(l_root) then
         if (inetcdf .eq. 1) then
           call start_phase(measurement_id_io_history)
           if((mhis .or. it.eq.1) .and. istage .eq. 1) then
              call open_his_netcdf
              call write_his_netcdf
              call close_his_netcdf
	   endif
           call end_phase(measurement_id_io_history)
         end if
      endif

c
c ------------ save histogram files
c
      if((mhis .or. it.eq.1)  .and. istage .eq. 1) then
          call start_phase(measurement_id_io_histograms)
          call open_histog_netcdf
          call write_histograms
          call close_histog_netcdf
          call end_phase(measurement_id_io_histograms)
      end if

c
c ------------ save viz files
c
      if((msave_v .or. it.eq.1) .and. istage .eq. 1) then
         call start_phase(measurement_id_io_viz)
         if (iviznetcdf) then
            call open_viz_netcdf
            call write_viz_netcdf
            call close_viz_netcdf
         end if
         call end_phase(measurement_id_io_viz)
      endif



c
c ------------ save velocity field
c
      if(msave .and. istage .eq. 1) then
         call start_phase(measurement_id_io_particles)
         call save_v
         call save_particles
         call end_phase(measurement_id_io_particles)
      endif


c
c ------------ save pressure field
c
      if(msave .and. istage .eq. 1) then
         call start_phase(measurement_id_io_pressure)
         call save_p
         call end_phase(measurement_id_io_pressure)
      endif


c
c --------- get rhs for all equations
c
      call start_phase(measurement_id_flow_solve_1)
      call comp1(istage,it)
      call end_phase(measurement_id_flow_solve_1)
      if(istage .eq. 1) then
         if(msave .and. l_root) call save_c
      endif

     
c
c --------- solve for pressure
c
      call start_phase(measurement_id_flow_solve_p)
      call comp_p
      call end_phase(measurement_id_flow_solve_p)
c
c --------- add pressure gradient and dealias
c
      call start_phase(measurement_id_flow_solve_2)
      call comp2
      if(micut) then
         call dealias
      endif
      call end_phase(measurement_id_flow_solve_2)


c      
c -------- update particles
c
         if (ipart_method .eq. 1) then
         call start_phase(measurement_id_particle_solver)
         call particle_update_rk3(istage)
         call end_phase(measurement_id_particle_solver)
         end if





c
 8999 continue

      call start_phase(measurement_id_humidity)
      call humidity_control
      call end_phase(measurement_id_humidity)

      if (ipart_method .eq. 2) then

         call start_phase(measurement_id_particle_solver)
         call particle_update_BE
         call end_phase(measurement_id_particle_solver)

      end if

      !Do the two-way coupling. Shouldn't include SFS,reintro, or coalescence
      !since these don't represent transfer between phases

      call start_phase(measurement_id_particle_coupling)
      call particle_coupling_update
      call end_phase(measurement_id_particle_coupling)

      call start_phase(measurement_id_particle_coupling_exchange)
      call particle_coupling_exchange
      call end_phase(measurement_id_particle_coupling_exchange)

      call start_phase(measurement_id_particle_coupling)
      call apply_particle_coupling
      call end_phase(measurement_id_particle_coupling)


      !Call coalescence outside of RK loop since it's not appropriate as
      !part of RK scheme
      if (icoalesce) then

         call start_phase(measurement_id_particle_coalesce)
         call particle_coalesce
         call end_phase(measurement_id_particle_coalesce)

      end if

      if (ipartdiff) then

        call start_phase(measurement_id_particle_diff)
        if (isfs == 2) then
           !Call Weil et al. (2004) Lagrangian stochastic model
           call SFS_velocity
        elseif(isfs == 1) then
           !Call stochastic model for the position based on vis_s
           call SFS_position
        end if
        call end_phase(measurement_id_particle_diff)

      end if

      if (ireintro) then
         call start_phase(measurement_id_particle_reintro)
         call particle_reintro
         call end_phase(measurement_id_particle_reintro)
      end if


      !After dust settles, calculate particle stats
      call start_phase(measurement_id_particle_stats)
      call particle_xy_stats
      call end_phase(measurement_id_particle_stats)


      !Calculate the histogram, save the particle trajectories
      if (mhis) then

         call start_phase(measurement_id_io_histograms)
         call radius_histogram
         call end_phase(measurement_id_io_histograms)

         call start_phase(measurement_id_io_traj)
         if (itrajout) then
         call particle_write_traj
         end if
         call end_phase(measurement_id_io_traj)
      end if


#ifdef TECIO
      if (mod(it,200) .eq. 0) then
      call start_phase(measurement_id_io_tecio)
      call plt_fields
      call end_phase(measurement_id_io_tecio)
      end if
#endif


      call get_max
      call get_dt

      call end_phase(measurement_id_timestepping_loop)


      if (it.ge.itmax) go to 99000

      if (time .gt. max_time) then  !Adjust itmax if max_time has been exceeded
         itmax = floor(real(it)/real(itape))*itape + itape  !Stop at next itape
      end if

      go to 9000
c
99000 continue

#ifdef TECIO
      call finalize_tecio
#endif

      call end_phase(measurement_id_solver)

      te_mpi = mpi_wtime()
      write(6,9997) (te_mpi - ts_mpi)
 9997 format(' Job Execution Time = ',e15.6)
c
 9998 continue

c --------- report the final measurements and shutdown the profiling.  take care
c           that we only report timings from a single rank and don't have its
c           output jumbled with every other ranks'.
      call mpi_barrier(mpi_comm_world,ierr)
      if(l_root) call report_profile(6)
      call shutdown_profiling

      call mpi_finalize(ierr)
c
      stop
      end



      subroutine get_fields
c
c ----------- special routine to read just 3d fields
c             as an initial guess, easy to customize
c             if missing data, etc..
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
c
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
      real, allocatable, dimension(:,:,:) :: temp
      logical there
c
      allocate(temp(nvar,nnx,iys:iye))
c
c ---------- input file to read from
c
!      path_ran = 'XXXXXXXXX/u.le.cou000'
c
c --------------------- get restart file from local directory
c                       reuse unit number
c
      close(nvel)
c
      inquire(file=path_ran,exist=there)
      if(there) then
         if(l_root) write(6,6001) path_ran
      else
         if(l_root) write(6,6005) path_ran
         stop
      endif
c
c ---- open file
c
      call mpi_file_open(mpi_comm_world, path_ran,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, nvel, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ------------ read 3d fields, make rhs*8
c
      nsize  = int(nvar,k8)*nnx*nny
      nsize2 = int(nvar,k8)*nnx*(iys-1)
      n_read = nvar*nnx*(iye+1-iys)
c
      do k=izs,ize
         offset = int((k-1),k8)*nsize + nsize2
         call mpi_file_read_at_all(nvel,offset,temp,n_read,
     +                              mpi_real8,status,ierr)
         if (ierr /= 0) goto 9992
#if defined(SWAP)
         call byteswap(temp)
#endif
         do j=iys,iye
         do i=1,nnx
            u(i,j,k) = temp(1,i,j) 
            v(i,j,k) = temp(2,i,j)
            w(i,j,k) = temp(3,i,j)
            e(i,j,k) = temp(nvar,i,j)
         enddo
         enddo
         do is = 1,nscl
            do j = iys,iye
            do i = 1,nnx
               t(i,j,is,k) = temp(3+is,i,j)
            enddo
            enddo
         enddo
c
      enddo
c
c ---- close file
c
      call mpi_file_close(nvel, ierr)
c
      deallocate(temp)
c
      do iz=izs,ize
c
         ug(iz) = 0.0
         vg(iz) = 0.0
c
c ---------------- initial guess for pressure
c
         do iy=iys,iye
         do ix=1,nnx
            p(ix,iy,iz) = 0.0
         enddo
         enddo
      enddo
c
      return
c ---------------------------- process errors
  100 continue
      write(6,9000) path_ran, nvel
      call mpi_finalize(ierr)
      stop
c
 9992 continue
      write(6,6100) nvel,iz
      call mpi_finalize(ierr)
      stop
c ---------------------
 6001 format(' SR. GET_FIELDS: FILE READ FOR INITIALIZATION = ',a80)
 6005 format(' 6005, SR. GET_FIELDS: cannot find restart file = ',a80)
 6100 format(' SR. GET_FIELDS: file read error on unit number = ',i2,/,
     +       '               at iz = ',i4)
 9000 format(' 9000, SR. GET_FIELDS: cannot open file =',a80,/,
     +       ' to unit number = ',i2)
      end
      subroutine set_sav(istart)
c
      use pars
      use fields
      use con_data
      use mod_io
      use con_stats
      use netcdf_io
c
      data ionce /0/
      save ionce
c
      if(it .ne. istart) then
c
c ------------------- increment time if not first time through
c
         time=time+dt
         viz_t_elapsed = viz_t_elapsed + dt
      endif
c
      it=it+1
c
      dt    = dt_new
      mnout = (mod(it,imean).eq.0).or. (it.eq.1)
      mtape = (mod(it,itape).eq.0)
      micut = (mod(it,itcut).eq.0)
      if(ihst .lt. 0) then
         mhis = .false.
      else
         !mhis = (mod(it,ihst).eq.0 .and. it.ge.it_his .or. it.eq.1)
	 mhis     = (floor(time/ihst) .ne. 
     +              floor((time+dt)/ihst))
      endif

      if (i_viz .lt. 0) then
         msave_v = .false.
      else
         !msave_v = (mod(it,i_viz).eq.0 .and. it .ge. it_viz)
	 msave_v = (floor(time/i_viz) .ne. 
     +              floor((time+dt)/i_viz))
      endif
c
c ---------- decide whether velocity fields are saved
c
      msave = .false.
      if(mtape) then
         itn=itn+1
         msave = .true.
         call get_output_filenames
      endif
c
c ---------- decide whether viz fields are saved
c
      if((i_viz .gt. 0) .and. (it .ge. it_viz_nxt)) then
        if (iviznetcdf) then
           !call close_viz_netcdf
           !call open_viz_netcdf
        else
           call viz_output_filename
        end if
        it_viz_nxt = it_viz_nxt + itape
      endif
c
c --------- decide whether history files are to be saved
c
      if((ihst .gt. 0) .and. (it .ge. it_his_nxt)) then
         if (inetcdf .eq. 1) then
            !call open_his_netcdf
            !call open_histog_netcdf
         else
            call open_histograms
            call open_his
         end if
         it_his_nxt = it_his_nxt + itape
      endif
c
      return
      end
      subroutine read_input_file
      use pars
      use particles
      use con_data
      implicit none

      character(48) :: label
      character(180) :: params_dir
      namelist /step_params/ itmax,imean,ihst,itape,
     +                       i_viz,max_time

      namelist /grid_params/ ncpu_s, Uo, Ttop, Tbot,
     +         qstar, tsfcc, ugcont, vgcont, dvdr, hurr_rad,
     +         zi, zl, xl, yl, zw1, dpdx, Swall,
     +         cfl,dt_new,surf_RH,surf_p,surf_rho,
     +         rad_kappa,rad_Fct,rad_Fcb

      namelist /path_names/ path_seed,path_part,path_res,
     +         path_ran

      namelist /flags/ ismlt,ifree,isfc,iradup,
     +         iupwnd,ibuoy,ifilt,itcut,isubs,ibrcl,iocean,
     +         method,idebug,iz_space,ivis0,ifix_dt,new_vis,iDNS,
     +         icouple,iTcouple,iHcouple,ievap,ifields,irestart,ilin,
     +         ineighbor,icoalesce,ipart_method,inewpart,icase,
     +         ipartdiff,isfs,iexner,ilongwave,ihurr,
     +         inetcdf,iviznetcdf,ihumiditycontrol,ireintro,
     +         itrajout, ikernel


      namelist /constants/ nuf, Cpa, Pra, Sc,
     +         tnumpart,mult_init,rhow,part_grav,rhos,
     +         Cpp,Mw,Ms,Ru,Gam,Ion,Sal,Rd,
     +         radius_init,radius_std,kappas_init,kappas_std,
     +         grav, t00,fcor,zo,zos,
     +         vp_init,Tp_init,nprime


      !params.in contains namelists to read
      !open(12, file="./params.in", status="old")
      call get_command_argument(1,params_dir)
      open(12,file=params_dir,status="old")

      read(12,nml=step_params)
      if (myid==0) print step_params

      read(12,nml=flags)
      if (myid==0) print flags

      read(12,nml=grid_params)
      if (myid==0) print grid_params

      read(12,nml=path_names)
      if (myid==0) print path_names

      read(12,nml=constants)
      if (myid==0) print constants
      CpaCpp = Cpa/Cpp



      end subroutine read_input_file
      function ran2(idum)
      integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real :: ran2, AM, EPS, RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     +     IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     +     IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER :: idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/,iv/NTAB*0/,iy/0/

      if (idum .le. 0) then
          idum=max(-idum,1)
          idum2 = idum
          do j = NTAB+8,1,-1
             k=idum/IQ1
             idum=IA1*(idum-k*IQ1)-k*IR1
             if (idum .lt. 0) idum=idum+IM1
             if (j .le. NTAB) iv(j) = idum
          end do
          iy=iv(1)
      end if
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum .lt. 0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2 .lt. 0) idum2=idum2+IM2
      j = 1+iy/NDIV
      iy = iv(j) - idum2
      iv(j) = idum
      if (iy .lt. 1) iy = iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      end function ran2
      function gasdev(idum)
       integer :: idum
       real :: gasdev
       integer :: iset
       real :: fac, gset, rsq, v1, v2, ran2
       save iset,gset
       data iset/0/

       if (idum .lt. 0) iset=0
       if (iset .eq. 0) then
 1000    v1 = 2.*ran2(idum)-1.0
         v2 = 2.*ran2(idum)-1.0
         rsq = v1**2+v2**2
         if ( (rsq .ge. 1) .or. (rsq .eq. 0)) goto 1000
         fac = sqrt(-2.0*log(rsq)/rsq)
         gset = v1*fac
         gasdev = v2*fac
         iset = 1
       else
         gasdev = gset
         iset = 0
       end if
       return
      end
      function cdf_func_single(d,CDF,M,S)
      implicit none
      real :: d, CDF, M, S, cdf_func_single

      cdf_func_single = 1.0/2.0*
     + (1.0 + erf((log(d) - M)/S/sqrt(2.0)))  - CDF

      end function cdf_func_single
      function cdf_func(d,CDF,factor,Ma,Sa,Mc,Sc)
      implicit none
      real :: d, CDF, factor, Ma, Sa, Mc, Sc, cdf_func

      cdf_func = 1.0/(1.0+factor)/2.0*
     + ( 1.0 + erf((log(d) - Ma)/Sa/sqrt(2.0)) +
     + factor*(1 + erf((log(d) - Mc)/Sc/sqrt(2.0))) ) - CDF

      end function cdf_func
