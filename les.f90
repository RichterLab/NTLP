      program les_mpi_2d
!
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
!
! ------------- definition of internal flags
!
!
!       iDNS    = 0; call the subgrid computation of vis_m and vis_s
!               = 1; call the molecular viscosity and diffusivity
!
!       igrdr   =  3; data comes from restart file
!               =  2; data comes from initialization (random)
!               =  1; data comes from coarser grid (or otherwise)
!
!       ibcu    = 1 ; upper boundary condition set by radiation bc
!               = 0 ; fixed value = 0.
!               = -1; value defined by coarser mesh for all variables
!
!       ibcl    = 0 ; lower boundary condition set by similarity theory (sr. setup)
!               = -1; value defined by coarser mesh for all variables
!
!       ifix_dt = 0 ; variable time step with fixed cfl number in setcon
!               = 1 ; fixed time step set in sr. get_dt
!
!       ifree   = 0 ; use spatially averaged surface conditions for MO (call lower)
!               = 1 ; use point-by-point conditions for MO free convection (call lower_free)
!
!       ihst    = nn; frequency at which global variables are output in history file
!               < 0 ; no history files
!
!
!       ismlt   = 1 ; use businger formulas in MO
!                 0 ; use large and everyone elses formulas in MO
!
!       iupwnd  = 0;  use skew symmetric formulas for all derivatives
!                     in scalar equations
!               = 1;  use hybrid upwind scheme for all derivatives
!                     in scalar equations
!
!       ivis0   = 0; old eddy viscosity model
!               = 1; new eddy viscosity model
!
!       new_vis = step; the iteration step for which the new model
!                       is turned on when ivis0=1
!               < 0; new model is on at all steps for ivis0=1
!
!       nscl  .ge. 1   number of scalars to be followed set in parameter statements
!                      change entries in sr. init, and sr. suft for surface bc's
!
! -------------------------------------------------------------------------------
!
! ---------- initialize MPI, get myid, numprocs,
!            test if on root process
!
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
      call mpi_comm_size(mpi_comm_world,numprocs,ierr)
!
      i_root = 0
      l_root = .false.
      if(myid .eq. i_root) l_root = .true.
!
      l_debug = .false.
      if(idebug .eq. 1) l_debug = .true.
!
      ts_mpi = mpi_wtime()

! ---------- initialize the profiling modules
      call initialize_profiling
      call start_phase(measurement_id_solver)

!----- Read the input file for all necessary parameters
      call start_phase(measurement_id_setup)

!
! ------------- establish association between pointers
!               and data structures
!
      call fill_cc
      call fill_cs
      call fill_ci

      call read_input_file

!
! -------- set number of x-y slab cpus
!
!      ncpu_s = 8
!
      case_inp = 'cou'
!
      call get_units
      call gridd
      call setcon

      call change_RH_bcs_to_q
!
! -------------- scratch run
!
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

!
! ---------- choose routine for getting initial guess
!
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

!
! --------------- time loop ------------
!
      tzero = time
      call get_dt
      do
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
        write(*,'(a15,e15.6,3i)') 'time,100,impos:',time,tnum100,tnumimpos
        end if



        if(it .ge. new_vis .and. ivis0 .eq. 1) then
            ivis = 1
        else
            ivis = 0
        endif


!
! ---------------- 3 stage runge-kutta time stepping
!
        t_stage_s = mpi_wtime()
        do istage=1,3
!
          dtzeta = dt*zetas(istage)
          dtgama = dt*gama(istage)
!
! ---------- compute derivatives of (u,v,w)
!

          call start_phase(measurement_id_derivatives)
          call exchange
          call get_derv
          call end_phase(measurement_id_derivatives)

!
! --------- new eddy viscosity, and bcs
!
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

!
! ------------ save history files
!
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

!
! ------------ save histogram files
!
          if((mhis .or. it.eq.1)  .and. istage .eq. 1) then
              call start_phase(measurement_id_io_histograms)
              call open_histog_netcdf
              call write_histograms
              call close_histog_netcdf
              call end_phase(measurement_id_io_histograms)
          end if

!
! ------------ save viz files
!
          if((msave_v .or. it.eq.1) .and. istage .eq. 1) then
             call start_phase(measurement_id_io_viz)
             if (iviznetcdf) then
                call open_viz_netcdf
                call write_viz_netcdf
                call close_viz_netcdf
             end if
             call end_phase(measurement_id_io_viz)
          endif



!
! ------------ save velocity field
!
          if(msave .and. istage .eq. 1) then
             call start_phase(measurement_id_io_flow)
             call save_v
             call end_phase(measurement_id_io_flow)
             call start_phase(measurement_id_io_particles)
             call save_particles
             call end_phase(measurement_id_io_particles)
          endif


!
! ------------ save pressure field
!
          if(msave .and. istage .eq. 1) then
             call start_phase(measurement_id_io_flow)
             call save_p
             call end_phase(measurement_id_io_flow)
          endif


!
! --------- get rhs for all equations
!
          call start_phase(measurement_id_flow_solve_1)
          call comp1(istage,it)
          call end_phase(measurement_id_flow_solve_1)
          if(istage .eq. 1) then
             if(msave .and. l_root) then
                call start_phase(measurement_id_io_flow)
                call save_c
                call end_phase(measurement_id_io_flow)
             endif
          endif


!
! --------- solve for pressure
!
          call start_phase(measurement_id_flow_solve_p)
          call comp_p
          call end_phase(measurement_id_flow_solve_p)
!
! --------- add pressure gradient and dealias
!
          call start_phase(measurement_id_flow_solve_2)
          call comp2
          if(micut) then
             call dealias
          endif
          call end_phase(measurement_id_flow_solve_2)


!
! -------- update particles
!
          if (ipart_method .eq. 1) then
          call start_phase(measurement_id_particle_solver)
          call particle_update_rk3(istage)
          call end_phase(measurement_id_particle_solver)
          end if

        end do  ! istage

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

        if (it >= itmax) exit

        if (time .gt. max_time) then  !Adjust itmax if max_time has been exceeded
           itmax = floor(real(it)/real(itape))*itape + itape  !Stop at next itape
        end if

      end do  ! time loop

#ifdef TECIO
      call start_phase(measurement_id_io_tecio)
      call finalize_tecio
      call end_phase(measurement_id_io_tecio)
#endif

      call end_phase(measurement_id_solver)

      te_mpi = mpi_wtime()
      write(6,9997) (te_mpi - ts_mpi)
 9997 format(' Job Execution Time = ',e15.6)
!
! --------- report the final measurements and shutdown the profiling.  take care
!           that we only report timings from a single rank and don't have its
!           output jumbled with every other ranks'.
      call mpi_barrier(mpi_comm_world,ierr)
      if(l_root) call report_profile(6)
      call shutdown_profiling

      call mpi_finalize(ierr)
!
      stop
      end
