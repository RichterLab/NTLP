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
      subroutine get_max
c
c --------- routine computes max velocities as sweep through
c           the velocity field 
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      real u_send(5), u_recv(5)
c

      dx_i = 1.0/dx
      dy_i = 1.0/dy
c
      u_temp = 0.0
      v_temp = 0.0
      w_temp = 0.0
      vis_temp = 0.0
      vis_temp = 0.0
      do iz=izs,ize
c
        u_xy = 0.0
        v_xy = 0.0
        w_xy = 0.0
        vis_xym = 0.0
        vis_xys = 0.0
        vis_xy = 0.0
        do iy=iys,iye
        do ix=1,nnx
           u_xy = amax1(u_xy,abs(u(ix,iy,iz)+stokes(iz)))
           v_xy = amax1(v_xy,abs(v(ix,iy,iz)))
           w_xy = amax1(w_xy,abs(w(ix,iy,iz)))
           vis_xym = amax1(vis_xym,vis_m(ix,iy,iz))
           vis_xys = amax1(vis_xys,vis_s(ix,iy,1,iz))
           vis_xy = amax1(vis_xys,vis_xym)
        enddo
        enddo
        u_xy   = u_xy*dx_i
        v_xy   = v_xy*dy_i
        wsav   = w_xy
        w_xy   = w_xy/abs(dzw(iz))
        vis_xy = vis_xy/amin1(dx,dy,dzw(iz))**2
c
        u_temp = amax1(u_xy,u_temp)
        v_temp = amax1(v_xy,v_temp)
        w_temp = amax1(w_xy,w_temp)
        vis_temp = amax1(vis_xy,vis_temp)
c
c       if(iz .le. 15) then
c         write(6,6000) iz, wmax
c6000     format(' in get_dt iz = ',i3,' wmax = ',e15.6)
c       endif
c
      enddo
      u_send(1) = u_temp
      u_send(2) = v_temp
      u_send(3) = w_temp
      u_send(4) = wsav
      u_send(5) = vis_temp 
  
c
      call mpi_allreduce(u_send,u_recv,5,mpi_real8,
     +     mpi_max,mpi_comm_world,ierror)
c
      umax = u_recv(1)
      vmax = u_recv(2)
      wmax = u_recv(3)
      wabs = u_recv(4)
      vismax = u_recv(5)


c
      return
      end
      subroutine get_dt
c
c ---------- routine computes max time step for given cfl number
c            from max's found previously
c
      use pars
      use con_data
      use con_stats
      use particles

c

      ucfl = umax
      vcfl = vmax
      wcfl = wmax
      ucflm = ucfl
      vcflm = vcfl
      wcflm = wcfl
      vel_max = wcflm
      vel_max = amax1(ucflm,vel_max)
      vel_max = amax1(vcflm,vel_max)

      if(vel_max .le. 0.0) then
          write(6,6000) ucflm, vcflm,wcflm, vel_max
 6000     format('6000, sr. get_dt bad news, umax = ',e15.6,/,
     +           ' vmax = ',e15.6,' wmax = ',e15.6,/,
     +           ' vel_max = ',e15.5,/,
     +           ' infinite time step !!!')
          call mpi_finalize(ierr)
          stop
      endif

c
c ---------------- choose fixed or variable time step
c
      if(ifix_dt .ne. 0) then
c
c ------------- if used, change to fit your problem
c
!        dt_new = 0.5  !Now in input file
      else
c
c ------------------- new estimate of best time step
c                     from cfl constraint
c
      dt_new  = cfl/vel_max
      dt_new = amin1(dt_new, 5.0)
c     dt_new = amin1(dt_new, 10.0)
      dt_cfl = dt_new
      endif

c ---------------- compare against viscous stability limit
c
      if (ifix_dt .eq. 0) then

      if(vismax*dt_new .gt. 0.5) then
         dt_new = 0.5/vismax
         if(l_root) then
            write(6,6200) dt_new, dt_cfl, vismax
 6200       format(' 6200 get_dt: cfl time step too large',/,
     +      '   viscous time step = ',e15.6,
     +      ' cfl time step = ',e15.6,' vismax = ',e15.6)
         endif
      endif

      if (icouple.eq.1 .or.
     +    iTcouple.eq.1 .or.
     +    iHcouple.eq.1) then

          call calc_tauc_min

      if (tauc_min .lt. dt_new) then
            dt_new = tauc_min
            if (l_root) then
            write(6,6210) dt_new, dt_cfl, tauc_min
 6210       format(' 6210 cfl time step too large for droplets',/,
     +      '   tauc time step = ',e15.6,
     +      ' cfl time step = ',e15.6,' tauc = ',e15.6)
            end if
      end if

      end if  !Particle coupling
      

      end if  !ifix_dt .eq. 0

c
      return
      end
      subroutine lterp(n,zary,zpt,i,ip1,ratio)
c
c ---- linear interpolation for zpt in zary, where zary is 
c      monotonic increasing or decreasing function
c
      dimension zary(*)
      nm1 = n-1
      if(n.le.1) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
      endif
      if(zary(1) .lt. zary(2)) go to 1
                               go to 101
    1 continue
c
c **** monotonic increasing array
c
        if(zpt .lt. zary(1)) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
        else if(zpt .gt. zary(n)) then
          i = n
          ip1 = n
          ratio = 1.0
          go to 999
        endif
        do j=1,nm1
              if(zpt .ge. zary(j) .and.
     $           zpt .le. zary(j+1)) then
                 i = j
                 ip1 = j+1
                 ratio = (zpt - zary(i))/(zary(ip1) - zary(i))
                 go to 999
              endif
        enddo
c
c **** decreasing array
c
  101 continue 
        if(zpt .gt. zary(1)) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
        else if(zpt .lt. zary(n)) then
          i = n
          ip1 = n
          ratio = 1.0
          go to 999
        endif
        do j=1,nm1
              if(zpt .le. zary(j) .and.
     $           zpt .ge. zary(j+1)) then
                 i = j
                 ip1 = j+1
                 ratio = (zpt - zary(i))/(zary(ip1) - zary(I))
                 go to 999
              endif
        enddo
  999 continue
      return
      end
      subroutine setcon
c
      use pars
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      external get_zi
c
c ----------------- get machine type, can erad datadir also
c
      open(99,file='mach.file')
      read(99,9000) imach
 9000 format(i1)
c
      close(99)
c
c ----------- initialize fft
c
      call rffti(nnx,trigx(1,1))
      call rffti(nny,trigx(1,2))
      call cffti(nny,trigc(1))
c
c ----------- start step for history files
c
      if (iti.eq.0) then
         it_his = 1
         it_viz = 1
      else
         it_his = iti+2
         it_viz = iti+2
      endif

      it_his_nxt = it_his
      it_viz_nxt = it_viz
c
c ---------------- set min value of e
c
      if(iocean .eq. 1) then
c
         smal_e = 0.0
         smal_e = 1.0e-12
c        smal_e = 6.0e-03
      else
         smal_e = 1.0e-09
c        smal_e = 0.0
      endif
c
c ---------------------- set constants in eddy viscosity model
c
      ck       = 0.1
      ceps     = 0.93
      csmag    = sqrt(ck*sqrt(ck/ceps))
      stab_c   = 0.76
c
c ----------------- set stability constant
c
      stabmin = 1.0e-12
c
c ---------------- minimum dsl length constant
c 
      almin_c = 0.0001
c
c -------------- initialize grid restart flag
c
      igrdr = 1
c
c -------------- create mpi operation to find max and location
c                using local gradient method
c
      call mpi_op_create(get_zi,.true.,ziloc,ierror)
c
c ------------------- define coefficients for 3-order runge-kutta
c                     time stepping scheme, borrowed from Spalart,
c                     Moser and Rogers, J. Comp. Physics 3/21/90
c                     Note this is a simplier version since all terms
c                     are lumped in the non-linear terms.
c                     cfl number is for an entire runge-kutta step
c                     in this case three stages. cfl = max(u)*dt/dx
c
      zetas(1) = 0.0
      zetas(2) = -17.0/60.0
      zetas(3) = -5.0/12.0
      gama(1)  = 8.0/15.0
      gama(2)  = 5.0/12.0
      gama(3)  = 3.0/4.0
      etas(1)  = -1.0
      etas(2)  = -1.0 + 8.0/15.0
      etas(3)  = -1.0 + 2.0/3.0
c
c ----------- a full step, at the new time
c
      etas(4)  =  0.0
c
c
      return
c
      end

      subroutine setup
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
c
c ------------ turn on new sgs model at a particular step
c
      if(it .ge. new_vis .and. ivis0 .eq. 1) then
          ivis = 1
      else
          ivis = 0
      endif
c
      if(igrdr . eq. 3) then
         if(l_root) then
            write(6,6)iti,utau,tsfcc(1) ,qstar(1)
            write(6,510)
            write(6,520) wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +            ,cdbtm,ugcont
            call print(6,1,nnz)
         endif
          if(l_debug) then
            write(nprt,6)iti,utau,tsfcc(1) ,qstar(1)
            write(nprt,510)
            write(nprt,520) wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +            ,cdbtm,ugcont
            call print(nprt,izs,ize)
          endif
      endif
c     if(ifilt.eq.1)call filter
      if(l_root) then
         write(6,1) nnx,nny,nnz,ismlt,ifilt,iti,itmax,
     +             iupwnd,ibuoy,itcut,
     +             dt,zo,tsfcc(1),isubs,ibrcl,
     +             method, iocean, ivis
      endif
      if(l_debug) then
         write(nprt,1) nnx,nny,nnz,ismlt,ifilt,iti,itmax,
     +             iupwnd,ibuoy,itcut,
     +             dt,zo,tsfcc(1),isubs,ibrcl,
     +             method, iocean, ivis
      endif
c
c -------------- boundary condition flags 
c
      ibcu = iradup
c     ibcu = 0
      ibcl = 0
c
c -------------------- wavenumbers, introduce a normalized
c                      set of wavenumbers to eliminate computation
c                      in derivatives , xderiv, yderiv
c
      do i=1,nnx
         xkn(i) = float(i-1)*pi2/xl
         if(i.gt.ncx)xkn(i) = -float(nnx-i+1)*pi2/xl
      enddo
      fn = 1.0/float(nnx)
      do i=1,nnx
         xk(i) = xkn(i)*fn
      enddo
      do i=1,nny
         ykn(i) = float(i-1)*pi2/yl
         if(i.gt.ncy)ykn(i) = -float(nny-i+1)*pi2/yl
      enddo
      fn = 1.0/float(nny)
      do i=1,nny
         yk(i) = ykn(i)*fn
      enddo
      ii = -1
      do i=1,ncx
         ii = ii + 2
         temp = xkn(i)**2
         do j=1,nny
            temp1       = temp + ykn(j)**2
            xks(ii,j)   = temp1
            xks(ii+1,j) = temp1
         enddo
      enddo
      xnn = abs(batag*dtdzf(1))
c
c ----------- choose correct sign so gravity waves
c             propagate out of the domain
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = -1.0
      if(ibcu.eq.1) then
         do iy=1,nny
         do ix=1,nnxp2
            if(xks(ix,iy) .le. 0.) then
              wavexy(ix,iy) = 0.0
            else
              wavexy(ix,iy) = sgn*sqrt(xnn/xks(ix,iy))
            endif
         enddo
         enddo
      endif
c
c -------------------- set length scale for SGS model
c
      if(iz_space .eq. 0) then
c
c ------------- uniform vertical spacing
c
      dx32 = dx*3./2.
      dy32 = dy*3./2.
      dsl  = (abs(dx32*dy32*dzw(1)))**(1./3.)
      dslg = dsl
      if(l_root)  write(6,2000) dsl
      if(l_debug) write(nprt,2000) dsl
c
c --------------------- create dsl array for easy indexing in comp1
c
      do iz=0,nnzp1
         dsl_z(iz) = dslg
      enddo
c
c ------------- variable vertical spacing
c
      else
c
c ----------- just estimate dsl for average spacing
c
         dx32 = dx*3./2.
         dy32 = dy*3./2.
c
         dsl_max = (abs(dx32*dy32*dzw(0)))**(1./3.)
         do iz=0,nnzp1
            dsl_z(iz) = (abs(dx32*dy32*dzw(iz)))**(1./3.)
            if(dsl_z(iz) .gt. dsl_max) dsl_max = dsl_z(iz)
         enddo
c        do iz=0,nnzp1
c           dsl_z(iz) = dsl_max
c        enddo
         dsl  = dsl_max
         dslg = dsl
      endif
c
      gridr = 1.0
      sml_eg = smal_e*gridr
c -------------------- set viscosity model parameters 
      if(ivis .ne. 1) then
        viscon = 0.0
        xksurf = 0.0
        nmatch = -1
        myid_newvis = 0
        do iz=1,nnz
           dfac(iz) = 1.0
        enddo
      endif
c ------------------- set stokes velocity for atmos/oceanic flow
      call stokesv
c
c --------- can add a time factor so as to skip into any part of
c           the specified geostrophic arrays. time factor in seconds
c
      t_factor = 7200.0
c
c ---------- for print out to get more digits
c
      t_ref = 300.0
c
c -------------------- specify cooling rate and initial
c                      temperature even for restarts
c
      c_rate   = 0.25/3600.0
      t_surf_i = 265.0
c
c -------------------- do not look for zi below zi_min
c
      zi_min = 30.0
      if(iocean .eq. 1) zi_min = -5.0
      iz_min = 1
      do iz=1,nnz-1
         if(zz(iz) .lt. zi_min .and.
     +      zz(iz+1) .ge. zi_min) iz_min = iz
      enddo
      if(l_root) then
         write(6,9000) zi_min, iz_min
      endif
c
 9998 continue
      return
c --------------------------- format statements
    6 format(///,' DATA FROM RESTART FILE AT STEP =',I5,
     +       ' U_* = ',e15.6,' TS = ',e15.6,' Q_* = ',e15.6,///)
  510 format(' RESTART ***** CASE WITH : ******',/)
  520 format(' WT = ',e12.4,', U_* = ',e12.4,', L = ',e12.4,
     +       ', DTDZ FREE = ',e12.4,', ZODY = ',e12.4,/,10x,
     +       '  ZO(BTM) = ',e12.4,', CDBTM = ',e12.4,
     +       ', UG = ',e12.4)
    1 format(10x,' NNX = ',i5,',  NNY = ',i5,
     + ',  NNZ = ',i5,/,10x,' SFC SMLT = ',i1,
     + ',  FILTER = ',i1,
     + ',  ITI = ',i6,',  ITMAX = ',i8,/,10x,
     + ' IUPWIND = ',i1,',  BUYNCY = ',i1,
     + ',  ITCUT = ',i1,/,10x,
     + ' DT = ',e15.6,',  ZO = ',e15.6,',  TS = ',e15.6,
     + ',  SUBSD = ',i1,/,
     + 10x,' BRCLICITY = ',i1,',  METHOD = ',i1,',  IOCEAN = ',i1,
     + ',  IVIS = ',i1)
 2000 format(10x,' DSL = ',e15.6)
 9000 format(' Search for zi above the height = ',e15.6,/,
     +       ' iz_min = ',i5)
      end
      subroutine nblnk(word)
      parameter (nmax=304)
      character wordt*304, word*(*)
      nchar = len(word)
      if(nchar .gt. nmax) then
         write(6,6000) nchar,nmax
 6000    format(' TROUBLE, IN SR. NBLNK : NCHAR = ',i6,
     +          ' EXCEEDS NMAX = ',i6)
         stop
      endif
      jj = 0
      do j=1,nchar
         if(word(j:j) .ne. ' ') then
            jj = jj + 1
            wordt(jj:jj) = word(j:j)
         endif
         word(j:j) = ' '
      enddo
      do j=1,jj
         word(j:j) = wordt(j:j)
      enddo
c
      return
      end
      subroutine blnk(word)
      character word*(*)
      nchar = len(word)
      do j=1,nchar
         word(j:j) = ' '
      enddo
c
      return
      end
      subroutine iso
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      real sfk(1:nnz)
c
c ---- get isotropy factor and scale it to match at the matching
c      height. uses boundary conditions from lower and upper. 
c
      do iz=1,nnz
c        dfac(iz) = 1.0
         dfac(iz) = 0.0
         sfk(iz)  = 0.0
      enddo
      do iz=izs,ize
         dfac(iz) = 1.0
      enddo
c
c ------ set nmatch equal to fraction of initial zi in sr. random
c
c     nmatch = izi/2
c     nmatch = 16
      nmatch = 48
      do i=0,numprocs-1,ncpu_s
         if(nmatch .ge. iz_s(i) .and.
     +      nmatch .le. iz_e(i)) myid_newvis = i
      enddo
c
      do iz=izs,min(ize,nmatch)
         izp1 = iz + 1
         izm1 = iz - 1
         weit = dzw(iz)/(dzw(iz) + dzw(izp1))
         weit1 = 1.0 - weit
c
c ---- get fluctuating strains
c
         do j=iys,iye
         do i=1,nnx
            s11 = weit1*ux(i,j,iz)**2 + weit*ux(i,j,izp1)**2
            s22 = weit1*vy(i,j,iz)**2 + weit*vy(i,j,izp1)**2
            wz  = (w(i,j,iz)-w(i,j,izm1))*dzw_i(iz)
            wzp = (w(i,j,izp1)-w(i,j,iz))*dzw_i(izp1)
            s33 = weit*wzp**2 + weit1*wz**2
            s12 = weit1*(uy(i,j,iz) + vx(i,j,iz))**2 +
     +            weit*(uy(i,j,izp1) + vx(i,j,izp1))**2
            s13 = (((u(i,j,izp1) - u(i,j,iz) +
     +            u_mn(iz) - u_mn(izp1))*dzu_i(izp1) +
     +            wx(i,j,iz)))**2
            s23 = (((v(i,j,izp1) - v(i,j,iz) +
     +          v_mn(iz) - v_mn(izp1))*dzu_i(izp1) +
     +          wy(i,j,iz)))**2
            sfk(iz) = sfk(iz) + 2.0*(s11 + s22 + s33) +
     +                       s12 + s13 + s23
         enddo
         enddo
         sfk(iz) = sfk(iz)*fnxy
      enddo
      call mpi_sum_z(sfk,i_root,myid,nnz,1)
c
      do iz=izs,min(ize,nmatch)
         izp1 = iz + 1
         izm1 = iz - 1
c
         sfk(iz) = sqrt(sfk(iz))
         smk = sqrt((u_mn(izp1)-u_mn(iz))**2 +
     +              (v_mn(izp1)-v_mn(iz))**2)*abs(dzu_i(izp1))
         if(sfk(iz) .le. 0. .and. smk .le. 0.) then
           dfac(iz) = 1.0
         else
           dfac(iz) = sfk(iz)/(sfk(iz) + smk)
         endif
c     if(l_root) write(6,6001) iz,sfk(iz),smk,dfac(iz)
 6001 format(' iz = ',i3,' sfk = ',e15.6,
     +       ' smk = ',e15.6,' dfac = ',e15.6)
      enddo
c
c
c ---- rescale ratio to give unity at match height
c      and if nested grid match value at upper boundary
c      of coarser grid
c
      if(myid .eq. myid_newvis) then
         dfacm = dfac(nmatch)
      endif
c
      call mpi_bcast(dfacm,1,mpi_real8,
     +              myid_newvis,mpi_comm_world,ierr)
c
      do iz=izs,min(ize,nmatch)
         dfac(iz) = dfac(iz)/dfacm
         dfac(iz) = amax1(dfac(iz), 0.1)
         dfac(iz) = amin1(dfac(iz), 1.0)
      enddo
c
c --------- gather dfac on all processes for printing and use in tke_vis
c           use reduce and divide by number of slab cpus
c
      call mpi_sum_z(dfac,i_root,myid,nnz,1)
      fncpu_s = 1.0/float(ncpu_s)
      do iz=1,nnz
         dfac(iz) = dfac(iz)*fncpu_s
      enddo
c
c     call mpi_gath_root(dfac(izs),dfac(1),iz_s,iz_e,izs,ize,nnz,
c    +                   myid,numprocs,ncpu_s)
c
c     if(l_root) write(6,6000) nmatch,ivis,(iz,dfac(iz),iz=1,nnz)
 6000 format(' in sr. iso, nmatch = ',i3,/,
     +       ' ivis = ',i3,'iz',5x,'dfac',/,(i3,1x,e15.6))
c     write(nprt,3001) (iz,dfac(iz),iz=1,nnz)
 3001 format(' iz ',5x,' dfac ',/,(i5,e15.6))
      return
      end
      subroutine surfvis
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      real xkvis(nnx,iys:iye), alwk(nnx,iys:iye)
c
      real send(3), buf(3)
c
      xksurf = 0.0
      viscon = 0.0
      vise   = 0.0
c
c ----------- only root process(es) compute 
c
      if(iss .eq. 0) then

c     ck = 0.1
c     csmag = 0.18
c     xkmax  = dzdz/dt/5.
      iz   = 1
      izm1 = iz - 1
      izp1 = iz + 1
c     xkmax  = dzu(izp1)*dzu(izp1)/(5.0*dt)
      dz_i = dzu_i(izp1)
      if(iocean .eq. 1) then
         call sufto
      else
         call suft
      endif
      if(qstar(1) .eq. 0.) then
         zeta = 0.0
      else
         zeta = abs(z(1))/amonin
      endif
      if(ismlt .eq. 1) then
          call busngr(zeta,phim,phis,psim,psis)
      else
          call fzol(zeta,phim,phis,psim,psis)
      endif
      viscon = vk*abs(z(1))/(utau*phim)
      vise   = utau*vk*abs(z(1))/phim
c
c ---- get special value at z1 to match with surface layer
c
      uws = 0.0
      vws = 0.0
      do iy=iys,iye
      do ix=1,nnx
         uws = uws + 0.5*(u(ix,iy,iz)-u_mn(iz) + 
     +         u(ix,iy,izp1) - u_mn(izp1))*(w(ix,iy,iz)-w_mn(iz))
         vws = vws + 0.5*(v(ix,iy,iz)-v_mn(iz) + 
     +         v(ix,iy,izp1) - v_mn(izp1))*(w(ix,iy,iz)-w_mn(iz))
      enddo
      enddo
      uws = uws*fnxy
      vws = vws*fnxy !Brian corrected from fxy 1/15/15
c
c ---- get average fluctuating eddy viscsoity
c
      do iy=iys,iye
      do ix=1,nnx
         e(ix,iy,iz)=amax1(e(ix,iy,iz),sml_eg)
      enddo
      enddo
      dslk = amin1(dsl,vk*abs(z(iz))/csmag)
c     stabmin = 1.e-12
c     almin = 0.0001*dsl
      almin = almin_c*dsl_z(iz)
      do iy=iys,iye
      do ix=1,nnx
         alwk(ix,iy)=dslk
c
c --------no stability corrected length scales when
c         new eddy viscosity is on
c
c         stab=batag*(t(ix,iy,1,izp1)-t(ix,iy,1,iz))*dz_i
c         if(stab.gt.stabmin) then
c           als = stab_c*sqrt(e(ix,iy,iz)/stab)
c           alwk(ix,iy) = amin1(dslk,als)
c         endif
c         alwk(ix,iy)=amax1(almin,alwk(ix,iy))
         xkvis(ix,iy)=ck*alwk(ix,iy)*sqrt(e(ix,iy,iz))*dfac(1)
c        xkvis(ix,iy)=amin1(xkvis(ix,iy),xkmax)
      enddo
      enddo
c
c ---- get average viscosity
c
      xkavg = 0.0
      do iy=iys,iye
      do ix=1,nnx
         xkavg = xkavg + xkvis(ix,iy)
      enddo
      enddo
      xkavg = xkavg*fnxy
c
      buf(1) = uws
      buf(2) = vws
      buf(3) = xkavg
      call mpi_sum_xy(buf,myid,iss,ise,3)
c
      uws   = buf(1)
      vws   = buf(2)
      xkavg = buf(3)
c
      xkz1 = vise - sqrt(uws**2 + vws**2)*viscon
      xksurf =  xkz1 - xkavg
      xksurf = amax1(xksurf,0.0)
      xksurf = amin1(xksurf,vise)
c     if(l_root) write(6,6000) dfac(1), xkavg, xkz1, vise, xksurf
 6000 format(' dfac = ',e12.4,' xkavg = ',e12.4,' xkz1 = ',e12.4,/,
     +       ' vise = ',e12.4,' xksurf = ',e12.4)
c
      endif
c
c ---------- broadcast values to other processes
c
      send(1) = xksurf
      send(2) = viscon
      send(3) = vise
c
      call mpi_bcast(send,3,mpi_real8,
     +              i_root,mpi_comm_world,ierr)
c
      xksurf = send(1)
      viscon = send(2)
      vise   = send(3)
c
      return
      end
      subroutine comp1(istep,it)
c
c ----- 3-order runge-kutta time stepping and monotone scalar fluxes in x,y,z.
c       designed to use mpi in x & y directions.
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      parameter(js = 6, ns = 3, nstat = js + ns*nscl)
      real stat(0:nnz,nstat)
c
c ------ temp arrays to hold rhs from step n-1 and
c        field variables from step n 
c
      real urhs(nnx,iys:iye,izs:ize), 
     +     vrhs(nnx,iys:iye,izs:ize),
     +     wrhs(nnx,iys:iye,izs:ize),
     +     erhs(nnx,iys:iye,izs:ize),
     +     trhs(nnx,iys:iye,nscl,izs:ize)
c

      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            urhs(ix,iy,iz) = u(ix,iy,iz) + dtzeta*r1(ix,iy,iz)
            vrhs(ix,iy,iz) = v(ix,iy,iz) + dtzeta*r2(ix,iy,iz)
            wrhs(ix,iy,iz) = w(ix,iy,iz) + dtzeta*r3(ix,iy,iz)
            erhs(ix,iy,iz) = e(ix,iy,iz) + dtzeta*r5(ix,iy,iz)
         enddo
         enddo
      enddo
      do iz=izs,ize
         do l=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            trhs(ix,iy,l,iz) = t(ix,iy,l,iz) + dtzeta*r4(ix,iy,l,iz)
         enddo
         enddo
         enddo
      enddo
c
c --------- get viscosity and rhs of (e,u,v,w)-equations
c           at next step
c
      if (iDNS .eq. 1) then
         call dns_vis
         call rhs_uvw_DNS(istep)
      else
         call tke_vis(istep) 
         call rhs_uvw(istep) 
      end if

c
c -------- evaluate rhs of scalar equations
c
      do l=1,nscl
         if (iDNS .eq. 1) then
            call rhs_scl_dns(istep,l) 
         else
            call rhs_scl(istep,l) 
         end if
      enddo

c
c ---------- gather stat sums on root processor
c            using mpi_reduction over all processors
c
      if(istep .eq. 1) then
c
        do j=1,nstat
        do iz=0,nnz
           stat(iz,j) = 0.0
        enddo
        enddo
        do iz=izs,ize
           stat(iz,1) = uwsb(iz)
           stat(iz,2) = vwsb(iz)
           stat(iz,3) = wwsb(iz)
!           stat(iz,4) = tr_tau(iz)
           stat(iz,5) = triz(iz)
           stat(iz,6) = shrz(iz)
           !stat(iz,7) = t_diss(iz)
        enddo
        m1 = js
        m2 = js + nscl
        m3 = js + 2*nscl
        do l=1,nscl
           do iz=izs,ize
              stat(iz,m1+l) = utsb(iz,l)
              stat(iz,m2+l) = vtsb(iz,l)
              stat(iz,m3+l) = wtsb(iz,l)
           enddo
        enddo

        !Populate the iz=0 locations if processor lies at the bottom -- only for certain quantities
        if (iss .eq. 0) then
           stat(0,1) = uwsb(0)
           stat(0,2) = vwsb(0)
           do l=1,nscl
              stat(0,m3+l) = wtsb(0,l)
           enddo
        end if

        call mpi_sum_z(stat(0,1),i_root,myid,nstat*(nnz+1),1)
        do iz=0,nnz
           uwsb(iz)   = stat(iz,1)
           vwsb(iz)   = stat(iz,2)
        enddo
        do iz=1,nnz
           wwsb(iz)   = stat(iz,3)
!           tr_tau(iz) = stat(iz,4)
           triz(iz)   = stat(iz,5)
           shrz(iz)   = stat(iz,6)
           !t_diss(iz) = stat(iz,7)
        enddo
        do l=1,nscl
           do iz=1,nnz
              utsb(iz,l) = stat(iz,m1+l)
              vtsb(iz,l) = stat(iz,m2+l)
           enddo
           do iz=0,nnz
              wtsb(iz,l) = stat(iz,m3+l)
           enddo
        enddo
        do iz=1,nnz
           buyz(iz) = batag*wtsb(iz,1)
        enddo

c
c -------- end if block
c
      endif
c
c ------- save old rhs in field variables for RK-advancement
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = urhs(ix,iy,iz)
            v(ix,iy,iz) = vrhs(ix,iy,iz)
            w(ix,iy,iz) = wrhs(ix,iy,iz)
            e(ix,iy,iz) = erhs(ix,iy,iz)
         enddo
         enddo
      enddo
      do iz=izs,ize
         do l=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,l,iz) = trhs(ix,iy,l,iz)
         enddo
         enddo
         enddo
      enddo
c
      return
      end
      subroutine rhs_uvw(istep)
c
c ---------- get right hand sides of (u,v,w) equations
c            for pencil size (nnx, iys:iye, izs:ize) 
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
c
      real fntd(nnx,iys:iye,izs:ize)
      real fnt1(nnx,iys:iye), fnt2(nnx,iys:iye) 
      real fnt3(nnx,iys:iye), fnt4(nnx,iys:iye)
      real tau13_u(nnx,iys:iye), tau23_u(nnx,iys:iye)
      real tau13_l(nnx,iys:iye), tau23_l(nnx,iys:iye)
      real r3_sum(1:nnz)
c

      do iz=izs,ize
c
      izm1 = iz - 1
      izp1 = iz + 1
      weit  = dzw(iz)/(dzw(iz) + dzw(izp1))
      weit1 = 1.0 - weit
c
c ---------- dynamics 
c
      do iy=iys,iye
      do ix=1,nnx
         uzm = (u(ix,iy,iz)-u(ix,iy,izm1))*dzu_i(iz)
         vzm = (v(ix,iy,iz)-v(ix,iy,izm1))*dzu_i(iz)
         uz  = (u(ix,iy,izp1)-u(ix,iy,iz))*dzu_i(izp1)
         vz  = (v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
c
         u_avg = u(ix,iy,iz)*weit1 + u(ix,iy,izp1)*weit
         v_avg = v(ix,iy,iz)*weit1 + v(ix,iy,izp1)*weit
c
c ------------ advection
c
         u_adv =  v(ix,iy,iz)*(vx(ix,iy,iz)-uy(ix,iy,iz))-
     +          0.5*(w(ix,iy,iz  )*(uz - wx(ix,iy,iz))+
     +           w(ix,iy,izm1)*(uzm - wx(ix,iy,izm1)))
         v_adv = -u(ix,iy,iz)*(vx(ix,iy,iz)-uy(ix,iy,iz))+
     +          0.5*(w(ix,iy,iz  )*(wy(ix,iy,iz) - vz)+
     +           w(ix,iy,izm1)*(wy(ix,iy,izm1) - vzm))
         w_adv = u_avg*(uz - wx(ix,iy,iz))
     +           - v_avg*(wy(ix,iy,iz) - vz)
c
c ------------ coriolis, vertical and horizontal components
c
         u_cor =  fcor*v(ix,iy,iz) - fcor_h*w(ix,iy,iz)
         v_cor = -fcor*(u(ix,iy,iz) + stokes(iz))
         w_cor =  fcor_h*u(ix,iy,iz)
c
c ------------ buoyancy (with hydrostatic part)
c

!MODIFY BUOYANCY TERM TO INCLUDE HUMIDITY
!         w_buy = batag*(t(ix,iy,1,iz)*weit1 +
!     +                  t(ix,iy,1,izp1)*weit)

         w_buy = bfac*grav/theta_base(iz)*
     +        ( (t(ix,iy,1,iz)*weit1 + t(ix,iy,1,izp1)*weit)*
     +        (1.0+0.61*(weit1*t(ix,iy,2,iz) + weit*t(ix,iy,2,izp1))) - 
     +        (theta_base(iz)*weit1 + theta_base(izp1)*weit) )


c
c ------------ geostrophic wind
c
        if (ihurr == 1) then
           u_geo = -fcor*vg(iz) - vg(iz)**2/hurr_rad  +
     +                 uxym(iz)**2/hurr_rad + vxym(iz)*vg(iz)/hurr_rad
           v_geo = -uxym(iz)*dvdr - uxym(iz)*vg(iz)/hurr_rad
        else
           u_geo = -fcor*vg(iz)
           v_geo =  fcor*(ug(iz)-ugal)
        end if

c
c ------------ totals
c
         r1(ix,iy,iz) = u_adv + u_cor + u_geo
         r2(ix,iy,iz) = v_adv + v_cor + v_geo
         r3(ix,iy,iz) = w_adv + w_cor + w_buy



         !Add particle momentum coupling
         if (icouple == 1) then
         r1(ix,iy,iz) = r1(ix,iy,iz) + partsrc(ix,iy,iz,1)
         r2(ix,iy,iz) = r2(ix,iy,iz) + partsrc(ix,iy,iz,2)
         !Note: partsrc(3,ix,iy,iz) is located at u,v-locations
         !Interpolate to w-location:
         r3(ix,iy,iz) = r3(ix,iy,iz) + (weit*partsrc(ix,iy,izp1,3)+
     +                  weit1*partsrc(ix,iy,iz,3))
         end if

      enddo
      enddo
c
c ---------------- stokes term in ocean cases
c
      if(iocean .eq. 1) then
        stokavg = stokes(iz)*weit1 + stokes(izp1)*weit
        do iy=iys,iye
        do ix=1,nnx
            r2(ix,iy,iz) = r2(ix,iy,iz) + stokes(iz)*
     +                    (uy(ix,iy,iz) - vx(ix,iy,iz))
            uz = (u(ix,iy,izp1) - u(ix,iy,iz))*dzu_i(izp1)
            r3(ix,iy,iz) = r3(ix,iy,iz) + stokavg* 
     +                    (uz - wx(ix,iy,iz))
        enddo
        enddo
      endif
c
c --------- get tau_13,_23 at iz-1 
c
      if (iz.ne.1 .or. ibcl.ne.0) then
         do iy=iys,iye
         do ix=1,nnx
            uzm = (u(ix,iy,iz)-u(ix,iy,izm1))*dzu_i(iz)
            vzm = (v(ix,iy,iz)-v(ix,iy,izm1))*dzu_i(iz)
            tau13_l(ix,iy) = -vis_m(ix,iy,izm1)*(uzm + wx(ix,iy,izm1)) -
     +             vis_mean(izm1)*(u_mn(iz)-u_mn(izm1))*dzu_i(iz)
            tau23_l(ix,iy) = -vis_m(ix,iy,izm1)*(vzm + wy(ix,iy,izm1)) -
     +             vis_mean(izm1)*(v_mn(iz)-v_mn(izm1))*dzu_i(iz)
         enddo
         enddo
      else
         do iy=iys,iye
         do ix=1,nnx
            tau13_l(ix,iy) = tau13m(ix,iy)
            tau23_l(ix,iy) = tau23m(ix,iy)
         enddo
         enddo
      endif
c
c ----------- x and z horizontal SGS fluxes for u, v, w
c             tau_11, tau_12, tau_13, tau_23 at iz
c
      do iy=iys,iye
      do ix=1,nnx
         fnt1(ix,iy) = -(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                    ux(ix,iy,iz)
         fnt2(ix,iy) = -.5*(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                    (uy(ix,iy,iz)+vx(ix,iy,iz))
         uz = (u(ix,iy,izp1)-u(ix,iy,iz))*dzu_i(izp1)
         vz = (v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
         tau13_u(ix,iy) = -vis_m(ix,iy,iz)*(uz+wx(ix,iy,iz)) -
     +            vis_mean(iz)*(u_mn(izp1)-u_mn(iz))*dzu_i(izp1)
         tau23_u(ix,iy) = -vis_m(ix,iy,iz)*(vz + wy(ix,iy,iz)) -
     +            vis_mean(iz)*(v_mn(izp1)-v_mn(iz))*dzu_i(izp1)
         fnt3(ix,iy) = tau13_u(ix,iy)
      enddo
      enddo
      call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      call xderivp(fnt2(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      call xderivp(fnt3(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      do iy=iys,iye
      do ix=1,nnx
         r1(ix,iy,iz) = r1(ix,iy,iz) - fnt1(ix,iy)
     +           -(tau13_u(ix,iy)-tau13_l(ix,iy))*dzw_i(iz)
         r2(ix,iy,iz) = r2(ix,iy,iz) - fnt2(ix,iy)
     +           -(tau23_u(ix,iy)-tau23_l(ix,iy))*dzw_i(iz)
         fnt4(ix,iy) = -(vis_m(ix,iy,izm1)+vis_m(ix,iy,iz))*
     +                (w(ix,iy,iz)-w(ix,iy,izm1))*dzw_i(iz)
         fnt2(ix,iy) = -(vis_m(ix,iy,izp1)+vis_m(ix,iy,iz))*
     +                (w(ix,iy,izp1)-w(ix,iy,iz))*dzw_i(izp1)
         r3(ix,iy,iz) = r3(ix,iy,iz) - fnt3(ix,iy) -
     +                   (fnt2(ix,iy)-fnt4(ix,iy))*dzu_i(izp1)
      enddo
      enddo
c
c -------- save SGS fluxes for printout
c
      if(istep .eq. 1) then
         uwsb(iz)   = 0.0
         vwsb(iz)   = 0.0
         wwsb(iz)   = 0.0
         tr_tau(iz) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            uwsb(iz) = uwsb(iz) + tau13_u(ix,iy)
            vwsb(iz) = vwsb(iz) + tau23_u(ix,iy)
            wwsb(iz) = wwsb(iz) + fnt4(ix,iy)
            ufluc    = (u(ix,iy,izp1) - uxym(izp1))*weit +
     +                 (u(ix,iy,iz) - uxym(iz))*weit1
            vfluc    = (v(ix,iy,izp1) - vxym(izp1))*weit +
     +                 (v(ix,iy,iz) - vxym(iz))*weit1
            tr_tau(iz) = tr_tau(iz) +
     +                 tau13_u(ix,iy)*ufluc + tau23_u(ix,iy)*vfluc
         enddo
         enddo
         uwsb(iz)   = uwsb(iz)*fnxy
         vwsb(iz)   = vwsb(iz)*fnxy
         wwsb(iz)   = wwsb(iz)*fnxy
         tr_tau(iz) = tr_tau(iz)*fnxy
      endif
c
c ---------- end z loop
c
      enddo
c
c ---------- SGS fluxes tau_12, tau_22, tau_23 that depend on 
c            y-derivatives 
c
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            fntd(ix,iy,iz) = -.5*(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                   (uy(ix,iy,iz)+vx(ix,iy,iz))
         enddo
         enddo
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            r1(ix,iy,iz)   = r1(ix,iy,iz) - fntd(ix,iy,iz)
            fntd(ix,iy,iz) = -(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                    vy(ix,iy,iz)
         enddo
         enddo
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do iz=izs,ize
         izp1 = iz + 1
         do iy=iys,iye
         do ix=1,nnx
            r2(ix,iy,iz)   = r2(ix,iy,iz) - fntd(ix,iy,iz)
            vz             = (v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
            fntd(ix,iy,iz) = -vis_m(ix,iy,iz)*(vz + wy(ix,iy,iz)) -
     +            vis_mean(iz)*(v_mn(izp1)-v_mn(iz))*dzu_i(izp1)
         enddo
         enddo
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do iz=1,nnz
         r3_sum(iz) = 0.0
      enddo
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            r3(ix,iy,iz) = r3(ix,iy,iz) - fntd(ix,iy,iz)
            r3_sum(iz)   = r3_sum(iz) + r3(ix,iy,iz)
         enddo
         enddo
         r3_sum(iz) = r3_sum(iz)*fnxy
      enddo
c
      call mpi_sum_z(r3_sum,i_root,myid,nnz,1)
c
c ------- make sure <r3> = 0 and set r3 = 0 at top
c
      do iz=izs,ize
         if(iz .eq. nnz) then
            do iy=iys,iye
            do ix=1,nnx
               r3(ix,iy,iz) = 0.0
            enddo
            enddo
         else
            do iy=iys,iye
            do ix=1,nnx
               r3(ix,iy,iz) = r3(ix,iy,iz) - r3_sum(iz)
            enddo
            enddo
         endif
      enddo
c
      return
      end
      subroutine rhs_uvw_DNS(istep)
c
c ---------- get right hand sides of (u,v,w) equations
c            for pencil size (nnx, iys:iye, izs:ize) 
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
c
      real fntd(nnx,iys:iye,izs:ize)
      real fnt1(nnx,iys:iye), fnt2(nnx,iys:iye) 
      real fnt3(nnx,iys:iye), fnt4(nnx,iys:iye)
      real tau13_u(nnx,iys:iye), tau23_u(nnx,iys:iye)
      real tau13_l(nnx,iys:iye), tau23_l(nnx,iys:iye)
      real r3_sum(1:nnz)
      real sfc_flx(2)
c
      do iz=izs,ize
c
      izm1 = iz - 1
      izp1 = iz + 1
      weit  = dzw(iz)/(dzw(iz) + dzw(izp1))
      weit1 = 1.0 - weit
c
c ---------- dynamics 
c
      do iy=iys,iye
      do ix=1,nnx
         uzm = (u(ix,iy,iz)-u(ix,iy,izm1))*dzu_i(iz)
         vzm = (v(ix,iy,iz)-v(ix,iy,izm1))*dzu_i(iz)
         uz  = (u(ix,iy,izp1)-u(ix,iy,iz))*dzu_i(izp1)
         vz  = (v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
c
         u_avg = u(ix,iy,iz)*weit1 + u(ix,iy,izp1)*weit
         v_avg = v(ix,iy,iz)*weit1 + v(ix,iy,izp1)*weit
c
c ------------ advection
c
         u_adv =  v(ix,iy,iz)*(vx(ix,iy,iz)-uy(ix,iy,iz))-
     +          0.5*(w(ix,iy,iz  )*(uz - wx(ix,iy,iz))+
     +           w(ix,iy,izm1)*(uzm - wx(ix,iy,izm1)))
         v_adv = -u(ix,iy,iz)*(vx(ix,iy,iz)-uy(ix,iy,iz))+
     +          0.5*(w(ix,iy,iz  )*(wy(ix,iy,iz) - vz)+
     +           w(ix,iy,izm1)*(wy(ix,iy,izm1) - vzm))
         w_adv = u_avg*(uz - wx(ix,iy,iz))
     +           - v_avg*(wy(ix,iy,iz) - vz)
c
c ------------ coriolis, vertical and horizontal components
c
         u_cor =  fcor*v(ix,iy,iz) - fcor_h*w(ix,iy,iz)
         v_cor = -fcor*(u(ix,iy,iz) + stokes(iz))
         w_cor =  fcor_h*u(ix,iy,iz)
c
c ------------ buoyancy (with hydrostatic part)
c
!         w_buy = batag*(t(ix,iy,1,iz)*weit1 +
!     +                  t(ix,iy,1,izp1)*weit)

         w_buy = bfac*grav/theta_base(iz)*
     +        ( (t(ix,iy,1,iz)*weit1 + t(ix,iy,1,izp1)*weit)*
     +        (1.0+0.61*(weit1*t(ix,iy,2,iz) + weit*t(ix,iy,2,izp1))) - 
     +        (theta_base(iz)*weit1 + theta_base(izp1)*weit) )
c
c ------------ geostrophic wind
c
         !u_geo = -fcor*vg(iz)
         !v_geo =  fcor*(ug(iz)-ugal)
         !Instead of geostrophic wind (which is a pressure gradient)
         !make u_geo and v_geo equal to my pressure gradient:

         u_geo = -dpdx
         v_geo = 0.0

c
c ------------ totals
c
         r1(ix,iy,iz) = u_adv + u_cor + u_geo
         r2(ix,iy,iz) = v_adv + v_cor + v_geo
         r3(ix,iy,iz) = w_adv + w_cor + w_buy


         !Add particle momentum coupling
         if (icouple == 1) then
         r1(ix,iy,iz) = r1(ix,iy,iz) + partsrc(ix,iy,iz,1)
         r2(ix,iy,iz) = r2(ix,iy,iz) + partsrc(ix,iy,iz,2)
         !Note: partsrc(3,ix,iy,iz) is located at u,v-locations
         !Interpolate to w-location:
         r3(ix,iy,iz) = r3(ix,iy,iz) + (weit*partsrc(ix,iy,izp1,3)+
     +                  weit1*partsrc(ix,iy,iz,3))
         end if

      enddo
      enddo
c
c ---------------- stokes term in ocean cases
c
      if(iocean .eq. 1) then
        stokavg = stokes(iz)*weit1 + stokes(izp1)*weit
        do iy=iys,iye
        do ix=1,nnx
            r2(ix,iy,iz) = r2(ix,iy,iz) + stokes(iz)*
     +                    (uy(ix,iy,iz) - vx(ix,iy,iz))
            uz = (u(ix,iy,izp1) - u(ix,iy,iz))*dzu_i(izp1)
            r3(ix,iy,iz) = r3(ix,iy,iz) + stokavg* 
     +                    (uz - wx(ix,iy,iz))
        enddo
        enddo
      endif
c
c --------- get tau_13,_23 at iz-1 
c
!      Have it compute t13,t23 like normal, even at bottom
!      REQUIRES ghost points to be set correctly for no-slip (done in lower,upper)
!      Also, get rid of the mean correction for that 2-part model
         sfc_flx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            uzm = (u(ix,iy,iz)-u(ix,iy,izm1))*dzu_i(iz)
            vzm = (v(ix,iy,iz)-v(ix,iy,izm1))*dzu_i(iz)
            tau13_l(ix,iy) = -vis_m(ix,iy,izm1)*(uzm + wx(ix,iy,izm1))
            tau23_l(ix,iy) = -vis_m(ix,iy,izm1)*(vzm + wy(ix,iy,izm1))
            if (iz == 1) then
                sfc_flx(1) = sfc_flx(1) + tau13_l(ix,iy)
                sfc_flx(2) = sfc_flx(2) + tau23_l(ix,iy)
            end if
         enddo
         enddo
!      else
!         do iy=iys,iye
!         do ix=1,nnx
!            tau13_l(ix,iy) = tau13m(ix,iy)
!            tau23_l(ix,iy) = tau23m(ix,iy)
!         enddo
!         enddo
!      endif

!      As an aside, compute the mean sfc_flx to put into uwsfc and vwsfc
       if (iz == 1) then 
          call mpi_sum_xy(sfc_flx,myid,iss,ise,2)
          uwsfc = sfc_flx(1)*fnxy
          vwsfc = sfc_flx(2)*fnxy
	  utau = sqrt(uwsfc**2 + vwsfc**2)
       end if
c
c ----------- x and z horizontal SGS fluxes for u, v, w
c             tau_11, tau_12, tau_13, tau_23 at iz
c
      do iy=iys,iye
      do ix=1,nnx
         fnt1(ix,iy) = -(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                    ux(ix,iy,iz)
         fnt2(ix,iy) = -.5*(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                    (uy(ix,iy,iz)+vx(ix,iy,iz))
         uz = (u(ix,iy,izp1)-u(ix,iy,iz))*dzu_i(izp1)
         vz = (v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
         tau13_u(ix,iy) = -vis_m(ix,iy,iz)*(uz+wx(ix,iy,iz))
         tau23_u(ix,iy) = -vis_m(ix,iy,iz)*(vz + wy(ix,iy,iz))
         fnt3(ix,iy) = tau13_u(ix,iy)
      enddo
      enddo
      call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      call xderivp(fnt2(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      call xderivp(fnt3(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      do iy=iys,iye
      do ix=1,nnx
         r1(ix,iy,iz) = r1(ix,iy,iz) - fnt1(ix,iy)
     +           -(tau13_u(ix,iy)-tau13_l(ix,iy))*dzw_i(iz)
         r2(ix,iy,iz) = r2(ix,iy,iz) - fnt2(ix,iy)
     +           -(tau23_u(ix,iy)-tau23_l(ix,iy))*dzw_i(iz)
         fnt4(ix,iy) = -(vis_m(ix,iy,izm1)+vis_m(ix,iy,iz))*
     +                (w(ix,iy,iz)-w(ix,iy,izm1))*dzw_i(iz)
         fnt2(ix,iy) = -(vis_m(ix,iy,izp1)+vis_m(ix,iy,iz))*
     +                (w(ix,iy,izp1)-w(ix,iy,iz))*dzw_i(izp1)
         r3(ix,iy,iz) = r3(ix,iy,iz) - fnt3(ix,iy) -
     +                   (fnt2(ix,iy)-fnt4(ix,iy))*dzu_i(izp1)
      enddo
      enddo
c
c -------- save SGS fluxes for printout
!          NOTE: now uwsb is t13_viscous,vwsb is t23_viscous, wwsb is t33_viscous
c
      if(istep .eq. 1) then
         uwsb(iz)   = 0.0
         vwsb(iz)   = 0.0
         wwsb(iz)   = 0.0
!         tr_tau(iz) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            uwsb(iz) = uwsb(iz) + tau13_u(ix,iy)
            vwsb(iz) = vwsb(iz) + tau23_u(ix,iy)
            wwsb(iz) = wwsb(iz) + fnt4(ix,iy)
            ufluc    = (u(ix,iy,izp1) - uxym(izp1))*weit +
     +                 (u(ix,iy,iz) - uxym(iz))*weit1
            vfluc    = (v(ix,iy,izp1) - vxym(izp1))*weit +
     +                 (v(ix,iy,iz) - vxym(iz))*weit1
!            tr_tau(iz) = tr_tau(iz) +
!     +                 tau13_u(ix,iy)*ufluc + tau23_u(ix,iy)*vfluc
         enddo
         enddo
         uwsb(iz)   = uwsb(iz)*fnxy
         vwsb(iz)   = vwsb(iz)*fnxy
         wwsb(iz)   = wwsb(iz)*fnxy
!         tr_tau(iz) = tr_tau(iz)*fnxy

         !Save the surface viscous stresses:
         if (iz==1) then
            uwsb(0) = 0.0
            vwsb(0) = 0.0
            do iy=iys,iye
            do ix=1,nnx
               uwsb(0) = uwsb(0) + tau13_l(ix,iy)
               vwsb(0) = vwsb(0) + tau23_l(ix,iy)
            end do
            end do
            uwsb(0) = uwsb(0)*fnxy
            vwsb(0) = vwsb(0)*fnxy
         end if
         

      endif
c
c ---------- end z loop
c
      enddo
c
c ---------- SGS fluxes tau_12, tau_22, tau_23 that depend on 
c            y-derivatives 
c
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            fntd(ix,iy,iz) = -.5*(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                   (uy(ix,iy,iz)+vx(ix,iy,iz))
         enddo
         enddo
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            r1(ix,iy,iz)   = r1(ix,iy,iz) - fntd(ix,iy,iz)
            fntd(ix,iy,iz) = -(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                    vy(ix,iy,iz)
         enddo
         enddo
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do iz=izs,ize
         izp1 = iz + 1
         do iy=iys,iye
         do ix=1,nnx
            r2(ix,iy,iz)   = r2(ix,iy,iz) - fntd(ix,iy,iz)
            vz             = (v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
            fntd(ix,iy,iz) = -vis_m(ix,iy,iz)*(vz + wy(ix,iy,iz))
         enddo
         enddo
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do iz=1,nnz
         r3_sum(iz) = 0.0
      enddo
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            r3(ix,iy,iz) = r3(ix,iy,iz) - fntd(ix,iy,iz)
            r3_sum(iz)   = r3_sum(iz) + r3(ix,iy,iz)
         enddo
         enddo
         r3_sum(iz) = r3_sum(iz)*fnxy
      enddo
c
      call mpi_sum_z(r3_sum,i_root,myid,nnz,1)
c
c ------- make sure <r3> = 0 and set r3 = 0 at top
c
      do iz=izs,ize
         if(iz .eq. nnz) then
            do iy=iys,iye
            do ix=1,nnx
               r3(ix,iy,iz) = 0.0
            enddo
            enddo
         else
            do iy=iys,iye
            do ix=1,nnx
               r3(ix,iy,iz) = r3(ix,iy,iz) - r3_sum(iz)
            enddo
            enddo
         endif
      enddo
c
      return
      end subroutine rhs_uvw_DNS
      subroutine rhs_scl(istep,iscl)
c
c ------ get right hand side of scalar equation (iscl)
c        monotone scalar fluxes only in z
c        for pencil size (nnx, iys:iye, izs:ize) 
c        care is taken so that if monotone is on then
c        conservative horizontal flux form is used!
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
c
c
      real fnt1(nnx,iys:iye,izs:ize)
      real tx(nnx,iys:iye), ty(nnx,iys:iye,izs:ize)
      real flux_u(nnx,iys:iye), flux_l(nnx,iys:iye)
      real taut3_u(nnx,iys:iye,nscl), taut3_l(nnx,iys:iye,nscl)
c
c --------- set sign for ocean simulations that use monotone
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = -1.0
      upwn = 2.0
      if(iupwnd .ne. 1) upwn = 1.0
c
c --------- outer loop over z
c
      do iz=izs,ize
c
      izm2 = iz - 2
      izm1 = iz - 1
      izp1 = iz + 1
      izp2 = iz + 2
      weit  = dzw(iz)/(dzw(iz) + dzw(izp1))
      weit1 = 1.0 - weit
      weit3 = dzw(izm1)/(dzw(iz) + dzw(izm1))
      weit4 = 1.0 - weit3
      dzw2_i = 1.0/(dzw(iz) + dzw(izp1))
      dzw3_i = 2.0*dzw2_i
c
      do iy=iys,iye
      do ix=1,nnx
         tx(ix,iy) = t(ix,iy,iscl,iz)
      enddo
      enddo
      call xderivp(tx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
c
c --------- compute tau_t3 at iz-1 
c
      if (iz.ne.1 .or. ibcl.ne.0) then
         do iy=iys,iye
         do ix=1,nnx
            taut3_l(ix,iy,iscl) = -vis_s(ix,iy,iscl,izm1)*
     +              (t(ix,iy,iscl,iz) - t(ix,iy,iscl,izm1))*dzu_i(iz)
         enddo
         enddo
      else
         do iy=iys,iye
         do ix=1,nnx
            taut3_l(ix,iy,iscl) = taut3m(ix,iy,iscl)
         enddo
         enddo
      endif
c
c ---------- SGS tau_t1, tau_t3 and resolved u*theta scalar fluxes
c            skew symmetric advective term 0.5(udt/dx + d/dx(ut))
c

!A modification to set the upwards flux equal to the flux at bottom for
!getting statistically steady humidity, temp fields in hurricane PBL
      if (iz.eq.nnz .and. icase.eq.4) then

      do iy=iys,iye
      do ix=1,nnx
         taut3_u(ix,iy,iscl) = qstar(iscl)

         fnt1(ix,iy,iz) = -0.5*((vis_s(ix,iy,iscl,iz)+
     +           vis_s(ix,iy,iscl,izm1))*tx(ix,iy) - 
     +           upwn*t(ix,iy,iscl,iz)*
     +                      (u(ix,iy,iz)+stokes(iz)))
      enddo
      enddo

      else

      do iy=iys,iye
      do ix=1,nnx
         taut3_u(ix,iy,iscl) = -vis_s(ix,iy,iscl,iz)*
     +      (t(ix,iy,iscl,izp1) - t(ix,iy,iscl,iz))*dzu_i(izp1)
         fnt1(ix,iy,iz) = -0.5*((vis_s(ix,iy,iscl,iz)+
     +           vis_s(ix,iy,iscl,izm1))*tx(ix,iy) - 
     +           upwn*t(ix,iy,iscl,iz)*
     +                      (u(ix,iy,iz)+stokes(iz)))
      enddo
      enddo

      end if

      call xderivp(fnt1(1,iys,iz),trigx(1,1),xk(1),nnx,iys,iye)
      do iy=iys,iye
      do ix=1,nnx
         r4(ix,iy,iscl,iz) = - fnt1(ix,iy,iz)
     +           -(taut3_u(ix,iy,iscl)-taut3_l(ix,iy,iscl))*dzw_i(iz)
      enddo
      enddo

c
      if(iupwnd .ne. 1) then
c
c --------- skew symmetric advective form for
c           vertical flux = 0.5(wdt/dz + d/dz(wt))
c
      do iy=iys,iye
      do ix=1,nnx
         theta_u = weit1*t(ix,iy,iscl,iz) +
     +                weit*t(ix,iy,iscl,izp1)
         theta_l = weit3*t(ix,iy,iscl,iz) +
     +                weit4*t(ix,iy,iscl,izm1)
         r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) 
     +     -0.5*(u(ix,iy,iz)+stokes(iz))*tx(ix,iy)
     +     -0.5*(w(ix,iy,iz)*theta_u - w(ix,iy,izm1)*theta_l)*dzw_i(iz)
c
         r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz)
     +     -0.25*(w(ix,iy,iz)*
     +       (t(ix,iy,iscl,izp1) - t(ix,iy,iscl,iz))*dzu_i(izp1) +
     +            w(ix,iy,izm1)*
     +       (t(ix,iy,iscl,iz) - t(ix,iy,iscl,izm1))*dzu_i(iz))
      enddo
      enddo
c
      else
c
c ----------- z-direction special
c
         if(iz .eq. 1) then
              do iy=iys,iye
              do ix=1,nnx
                 flux_l(ix,iy) = sgn*0.5*w(ix,iy,izm1)*
     +                        (t(ix,iy,iscl,izm1)+t(ix,iy,iscl,iz))
                 flux_u(ix,iy) =
     +           amax1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izp1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izm1))) +
     +           amin1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,izp1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izp1),
     +                t(ix,iy,iscl,izp2)))
              enddo
              enddo
         else if(iz .eq. nnz) then
              do iy=iys,iye
              do ix=1,nnx
                 flux_u(ix,iy) = sgn*0.5*w(ix,iy,iz)*
     +                        (t(ix,iy,iscl,izp1)+t(ix,iy,iscl,iz))
                 flux_l(ix,iy) =
     +           amax1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,izm1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izm1),
     +                t(ix,iy,iscl,izm2))) +
     +           amin1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izm1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izp1)))
              enddo
              enddo
         else
              do iy=iys,iye
              do ix=1,nnx
                 flux_u(ix,iy) =
     +           amax1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izp1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izm1))) +
     +           amin1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,izp1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izp1),
     +                t(ix,iy,iscl,izp2)))
                 flux_l(ix,iy) =
     +           amax1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,izm1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izm1),
     +                t(ix,iy,iscl,izm2))) +
     +           amin1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izm1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izp1)))
              enddo
              enddo
         endif
c
c ---------- sum vertical monotone flux
c
         do iy=iys,iye
         do ix=1,nnx
            r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz)
     +          - sgn*(flux_u(ix,iy) - flux_l(ix,iy))*dzw_i(iz)
         enddo
         enddo
c
c -------- end monotone if block
c
      endif
c
c -------- save SGS fluxes for printout, gather sums on exit
c
      if(istep .eq. 1) then
         utsb(iz,iscl) = 0.0
         wtsb(iz,iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            wtsb(iz,iscl) = wtsb(iz,iscl) + taut3_u(ix,iy,iscl)
            utsb(iz,iscl) = utsb(iz,iscl) -
     +            0.5*(vis_s(ix,iy,iscl,iz)+
     +                    vis_s(ix,iy,iscl,izm1))*tx(ix,iy)
         enddo
         enddo
         utsb(iz,iscl) = utsb(iz,iscl)*fnxy
         wtsb(iz,iscl) = wtsb(iz,iscl)*fnxy
      endif
c
c ---------- end z loop
c
      enddo
c
c --------- outer loop over z for y-depenence
c
      do iz=izs,ize
      do iy=iys,iye
      do ix=1,nnx
         ty(ix,iy,iz)  = t(ix,iy,iscl,iz)
      enddo
      enddo
      enddo
c
c --------- y-derivative of t for [izs:ize]
c
      call yd_mpi(ty(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
c ------------- add skew symmetric advective flux and SGS flux
c               to y-derivative computation. check for monotone
c
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = -0.5*(
     +           (vis_s(ix,iy,iscl,iz)+vis_s(ix,iy,iscl,izm1))*
     +                 ty(ix,iy,iz) - upwn*t(ix,iy,iscl,iz)*v(ix,iy,iz))
         enddo
         enddo
         if(iupwnd .ne. 1) then
           do iy=iys,iye
           do ix=1,nnx
              r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) - 
     +                            0.5*v(ix,iy,iz)*ty(ix,iy,iz)
           enddo
           enddo
         endif
      enddo
c
c --------- y-derivatives of scalar fluxes for [izs:ize]
c
      call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +            nnx,nny,ixs,ixe,ix_s,ix_e,
     +            iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) - fnt1(ix,iy,iz)
         enddo
         enddo
      enddo


!---------add on the thermal coupling from the particles:
      if (iscl == 1) then
      if (iTcouple == 1) then
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
          if (iexner) then
          r4(ix,iy,1,iz) = r4(ix,iy,1,iz) + partTsrc(ix,iy,iz)/
     +    exner(surf_p,func_p_base(surf_p,tsfcc(1),zz(iz)))
          else
             r4(ix,iy,1,iz) = r4(ix,iy,1,iz) + partTsrc(ix,iy,iz)
          end if    
         end do
         end do
      end do
      end if
      end if


      if (iscl == 2) then
      if (iHcouple == 1) then
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
          r4(ix,iy,1,iz) = r4(ix,iy,1,iz) + partTEsrc(ix,iy,iz)
          r4(ix,iy,2,iz) = r4(ix,iy,2,iz) + partHsrc(ix,iy,iz)
         end do
         end do
      end do
      end if
      end if


!---------add on simple radiation scheme:
      if (ilongwave == 1) then
      if (iscl == 1) then

         call calc_radiation

         do iz=izs,ize
           do iy=iys,iye
           do ix=1,nnx
           !calc_radiation fills "radsrc" vector
           r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) + radsrc(iz)

           end do
           end do
         end do

      end if
      end if

c
c -------- save SGS fluxes for printout
c
      if(istep .eq. 1) then
      do iz=izs,ize
         vtsb(iz,iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vtsb(iz,iscl) = vtsb(iz,iscl) -
     +            0.5*(vis_s(ix,iy,iscl,iz)+
     +              vis_s(ix,iy,iscl,izm1))*ty(ix,iy,iz)
         enddo
         enddo
         vtsb(iz,iscl) = vtsb(iz,iscl)*fnxy
      enddo
      endif
c
      return
      end
      subroutine rhs_scl_dns(istep,iscl)
c
c ------ get right hand side of scalar equation (iscl)
c        monotone scalar fluxes only in z
c        for pencil size (nnx, iys:iye, izs:ize) 
c        care is taken so that if monotone is on then
c        conservative horizontal flux form is used!
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
c
c
      real fnt1(nnx,iys:iye,izs:ize)
      real tx(nnx,iys:iye), ty(nnx,iys:iye,izs:ize)
      real flux_u(nnx,iys:iye), flux_l(nnx,iys:iye)
      real taut3_u(nnx,iys:iye,nscl), taut3_l(nnx,iys:iye,nscl)
      real :: sfc_flx
c
c --------- set sign for ocean simulations that use monotone
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = -1.0
      upwn = 2.0
      if(iupwnd .ne. 1) upwn = 1.0
c
c --------- outer loop over z
c
      do iz=izs,ize
c
      izm2 = iz - 2
      izm1 = iz - 1
      izp1 = iz + 1
      izp2 = iz + 2
      weit  = dzw(iz)/(dzw(iz) + dzw(izp1))
      weit1 = 1.0 - weit
      weit3 = dzw(izm1)/(dzw(iz) + dzw(izm1))
      weit4 = 1.0 - weit3
      dzw2_i = 1.0/(dzw(iz) + dzw(izp1))
      dzw3_i = 2.0*dzw2_i
c
      do iy=iys,iye
      do ix=1,nnx
         tx(ix,iy) = t(ix,iy,iscl,iz)
      enddo
      enddo
      call xderivp(tx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
c
c --------- compute tau_t3 at iz-1 
c

         sfc_flx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            taut3_l(ix,iy,iscl) = -vis_s(ix,iy,iscl,izm1)*
     +              (t(ix,iy,iscl,iz) - t(ix,iy,iscl,izm1))*dzu_i(iz)
         if (iz == 1) then 
             sfc_flx = sfc_flx + taut3_l(ix,iy,iscl)
         end if
         enddo
         enddo


!      Aside - get each wtsfc value:
       if (iz == 1 .and. isfc(iscl) == 1) then
          call mpi_sum_xy(sfc_flx,myid,iss,ise,1)
          wtsfc(iscl) = sfc_flx*fnxy 
       end if
c
c ---------- SGS tau_t1, tau_t3 and resolved u*theta scalar fluxes
c            skew symmetric advective term 0.5(udt/dx + d/dx(ut))
c
      do iy=iys,iye
      do ix=1,nnx
         taut3_u(ix,iy,iscl) = -vis_s(ix,iy,iscl,iz)*
     +   (t(ix,iy,iscl,izp1) - t(ix,iy,iscl,iz))*dzu_i(izp1)
c
         fnt1(ix,iy,iz) = -0.5*((vis_s(ix,iy,iscl,iz)+
     +                  vis_s(ix,iy,iscl,izm1))*
     +                    tx(ix,iy) - upwn*t(ix,iy,iscl,iz)*
     +                      (u(ix,iy,iz)+stokes(iz)))
      enddo
      enddo
      call xderivp(fnt1(1,iys,iz),trigx(1,1),xk(1),nnx,iys,iye)
      do iy=iys,iye
      do ix=1,nnx
         r4(ix,iy,iscl,iz) = - fnt1(ix,iy,iz)
     +           -(taut3_u(ix,iy,iscl)-taut3_l(ix,iy,iscl))*dzw_i(iz)
      enddo
      enddo
c
      if(iupwnd .ne. 1) then
c
c --------- skew symmetric advective form for
c           vertical flux = 0.5(wdt/dz + d/dz(wt))
c
      do iy=iys,iye
      do ix=1,nnx
         theta_u = weit1*t(ix,iy,iscl,iz) +
     +                weit*t(ix,iy,iscl,izp1)
         theta_l = weit3*t(ix,iy,iscl,iz) +
     +                weit4*t(ix,iy,iscl,izm1)
         r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) 
     +     -0.5*(u(ix,iy,iz)+stokes(iz))*tx(ix,iy)
     +     -0.5*(w(ix,iy,iz)*theta_u - w(ix,iy,izm1)*theta_l)*dzw_i(iz)
c
         r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz)
     +     -0.25*(w(ix,iy,iz)*
     +       (t(ix,iy,iscl,izp1) - t(ix,iy,iscl,iz))*dzu_i(izp1) +
     +            w(ix,iy,izm1)*
     +       (t(ix,iy,iscl,iz) - t(ix,iy,iscl,izm1))*dzu_i(iz))
      enddo
      enddo
c
      else
c
c ----------- z-direction special
c
         if(iz .eq. 1) then
              do iy=iys,iye
              do ix=1,nnx
                 flux_l(ix,iy) = sgn*0.5*w(ix,iy,izm1)*
     +                        (t(ix,iy,iscl,izm1)+t(ix,iy,iscl,iz))
                 flux_u(ix,iy) =
     +           amax1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izp1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izm1))) +
     +           amin1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,izp1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izp1),
     +                t(ix,iy,iscl,izp2)))
              enddo
              enddo
         else if(iz .eq. nnz) then
              do iy=iys,iye
              do ix=1,nnx
                 flux_u(ix,iy) = sgn*0.5*w(ix,iy,iz)*
     +                        (t(ix,iy,iscl,izp1)+t(ix,iy,iscl,iz))
                 flux_l(ix,iy) =
     +           amax1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,izm1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izm1),
     +                t(ix,iy,iscl,izm2))) +
     +           amin1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izm1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izp1)))
              enddo
              enddo
         else
              do iy=iys,iye
              do ix=1,nnx
                 flux_u(ix,iy) =
     +           amax1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izp1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izm1))) +
     +           amin1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,izp1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izp1),
     +                t(ix,iy,iscl,izp2)))
                 flux_l(ix,iy) =
     +           amax1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,izm1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izm1),
     +                t(ix,iy,iscl,izm2))) +
     +           amin1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izm1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izp1)))
              enddo
              enddo
         endif
c
c ---------- sum vertical monotone flux
c
         do iy=iys,iye
         do ix=1,nnx
            r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz)
     +          - sgn*(flux_u(ix,iy) - flux_l(ix,iy))*dzw_i(iz)
         enddo
         enddo
c
c -------- end monotone if block
c
      endif
c
c -------- save SGS fluxes for printout, gather sums on exit
c
      if(istep .eq. 1) then
         utsb(iz,iscl) = 0.0
         wtsb(iz,iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            wtsb(iz,iscl) = wtsb(iz,iscl) + taut3_u(ix,iy,iscl)
            utsb(iz,iscl) = utsb(iz,iscl) -
     +            0.5*(vis_s(ix,iy,iscl,iz)+
     +                     vis_s(ix,iy,iscl,izm1))*tx(ix,iy)
         enddo
         enddo
         utsb(iz,iscl) = utsb(iz,iscl)*fnxy
         wtsb(iz,iscl) = wtsb(iz,iscl)*fnxy

         !Do it special for wtsb the lower surface:
         if (iz==1) then
         wtsb(0,iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            wtsb(0,iscl) = wtsb(0,iscl) + taut3_l(ix,iy,iscl)
         enddo
         enddo
         wtsb(0,iscl) = wtsb(0,iscl)*fnxy
         end if

      endif
c
c ---------- end z loop
c
      enddo
c
c --------- outer loop over z for y-depenence
c
      do iz=izs,ize
      do iy=iys,iye
      do ix=1,nnx
         ty(ix,iy,iz)  = t(ix,iy,iscl,iz)
      enddo
      enddo
      enddo
c
c --------- y-derivative of t for [izs:ize]
c
      call yd_mpi(ty(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
c ------------- add skew symmetric advective flux and SGS flux
c               to y-derivative computation. check for monotone
c
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = -0.5*((vis_s(ix,iy,iscl,iz)+
     +              vis_s(ix,iy,iscl,izm1))*ty(ix,iy,iz) -
     +                  upwn*t(ix,iy,iscl,iz)*v(ix,iy,iz))
         enddo
         enddo
         if(iupwnd .ne. 1) then
           do iy=iys,iye
           do ix=1,nnx
              r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) - 
     +                            0.5*v(ix,iy,iz)*ty(ix,iy,iz)
           enddo
           enddo
         endif
      enddo

!---------add on the thermal coupling from the particles:
      if (iscl == 1) then 
      if (iTcouple == 1) then
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
          r4(ix,iy,1,iz) = r4(ix,iy,1,iz) + partTsrc(ix,iy,iz)
         end do
         end do
      end do
      end if
      end if  

      if (iscl == 2) then 
      if (iHcouple == 1) then
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
          r4(ix,iy,1,iz) = r4(ix,iy,1,iz) + partTEsrc(ix,iy,iz)
          r4(ix,iy,2,iz) = r4(ix,iy,2,iz) + partHsrc(ix,iy,iz)
         end do
         end do
      end do
      end if
      end if 

!-----------add on the humidity wall sink 

      if (iscl==2) then

      do iz=izs,ize
        do iy=iys,iye
        do ix=1,nnx
            r4(ix,iy,2,iz) = r4(ix,iy,2,iz) + Swall
        enddo
        enddo
      enddo
      endif



c
c --------- y-derivatives of scalar fluxes for [izs:ize]
c
      call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +            nnx,nny,ixs,ixe,ix_s,ix_e,
     +            iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) - fnt1(ix,iy,iz)
         enddo
         enddo
      enddo
c
c -------- save SGS fluxes for printout
c
      if(istep .eq. 1) then
      do iz=izs,ize
         vtsb(iz,iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vtsb(iz,iscl) = vtsb(iz,iscl) -
     +            0.5*(vis_s(ix,iy,iscl,iz)+
     +                 vis_s(ix,iy,iscl,izm1))*ty(ix,iy,iz)
         enddo
         enddo
         vtsb(iz,iscl) = vtsb(iz,iscl)*fnxy
      enddo
      endif
c
      return
      end subroutine rhs_scl_dns
      subroutine dns_vis
      use particles
      use pars
      use fields
      implicit none


!     In DNS mode, just set the molecular viscosity (and scalar diffusivities)
!     Also, to make the rest of code work, set the rhs of e equation to 0

      !Both for air at the moment:
      vis_m = nuf 

      !Use Prantdl number for thermal diffusivity:
      vis_s(1:nnx,iys:iye,1,izs-1:ize+1) = nuf/Pra   ! alpha=nu/Prandtl   
      vis_s(1:nnx,iys:iye,2,izs-1:ize+1) = nuf/Sc    !Dv =nu/Sc
      r5 = 0.0
      e = 0.0

      end
      subroutine tke_vis(istep)
c
c -------------- get viscosity using deardorff tke model with
c                stability correction. fixes for surface layer. 
c                 get rhs of e-equation
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      real fnt1(nnx,iys:iye), fnt2(nnx,iys:iye,izs:ize)
      real fnt3(nnx,iys:iye)
      real ex(nnx,iys:iye), ey(nnx,iys:iye,izs:ize)
      real u_avg(nnx,iys:iye), v_avg(nnx,iys:iye), dissp(nnx,iys:iye)
      real alk(nnx,iys:iye,izs-1:ize+1)
c
      do iz=izs-1,ize+1
c
      izp1 = iz + 1
      dslk  = dsl_z(iz)
      if(iz .gt. 0) dslk  = amin1(dsl_z(iz),vk*abs(z(iz))/csmag)
      almin = almin_c*dsl_z(iz)
      if(iz .eq. 0 .or. iz .eq. nnz+1) then
         dfack = 1.0
      else
         dfack = dfac(iz)
      endif
      if(ivis .eq. 1 .and. iz .le. nmatch) then
c
c --------------- no stability corrected length scales
c
         do j=iys,iye
         do i=1,nnx
            alk(i,j,iz) = dslk
         end do
         end do
      else
         do j=iys,iye
         do i=1,nnx
            alk(i,j,iz) = dslk
            stab = batag*(t(i,j,1,izp1) - t(i,j,1,iz))*dzu_i(izp1)
            if(stab.gt.stabmin) then
              als = stab_c*sqrt(e(i,j,iz)/stab)
              alk(i,j,iz) = amin1(dslk,als)
            endif
            alk(i,j,iz)  = amax1(almin,alk(i,j,iz))
         enddo
         enddo
      endif
      do j=iys,iye
      do i=1,nnx
         vis_m(i,j,iz) = ck*alk(i,j,iz)*sqrt(e(i,j,iz))*dfack
         vis_s(i,j,1:nscl,iz) = (1.+2.*alk(i,j,iz)/dslk)*vis_m(i,j,iz) 
      enddo
      enddo
c
c -------------- special case for iz = 1
c
      if(iz.eq.1 .and. ibcl .eq. 0) then
         do iy=iys,iye
         do ix=1,nnx
            vis_m(ix,iy,iz-1) = vis_m(ix,iy,iz)
            vis_s(ix,iy,1:nscl,iz-1) = vis_s(ix,iy,1:nscl,iz)
         enddo
         enddo
      endif
c
c -------------- end z loop
c
      enddo
c
c -------------- if special 2 part surface layer model is on
c                get "mean" viscosity
c
      do iz=izs-1,ize
         izm1         = iz - 1
         izp1         = iz + 1
         vis_mean(iz) = 0.0
         if(ivis .eq. 1 .and. iz .le. nmatch) then
            if(iz .le. 1) then
              vis_mean(iz) = xksurf
            else
              stravg = sqrt((u_mn(izp1)-u_mn(iz))**2 + 
     +              (v_mn(izp1)-v_mn(iz))**2)*abs(dzu_i(izp1))
              vis_mean(iz) = xksurf*viscon*stravg
            endif
         endif
      enddo
c
c --------- update rhs of sgs e from x and z pieces
c           cube of size (nnx, iys,iye, izs:ize)
c
      do iz=izs,ize
c
      izm1   = iz - 1
      izp1   = iz + 1
      weit   = dzw(iz)/(dzw(iz) + dzw(izp1))
      weit1  = 1.0 - weit
      dzw2_i = 1.0/(dzw(iz) + dzw(izp1))
      dzw3_i = 2.0*dzw2_i
      dslk   = dsl_z(iz)
c
      do iy=iys,iye
      do ix=1,nnx
         ex(ix,iy) = e(ix,iy,iz)
      enddo
      enddo
      call xderivp(ex(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
c
c ------------ include stokes contribution in advection
c              and horizontal x-diffusion
c
      do iy=iys,iye
      do ix=1,nnx
         u_avg(ix,iy)   = (stokes(iz) + u(ix,iy,iz))*weit1 +
     +                    (stokes(izp1) + u(ix,iy,izp1))*weit
         fnt1(ix,iy)    = e(ix,iy,iz)*u_avg(ix,iy) - 
     +                    4.0*vis_m(ix,iy,iz)*ex(ix,iy)
      enddo
      enddo
      call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      do iy=iys,iye
      do ix=1,nnx
         r5(ix,iy,iz) = -fnt1(ix,iy) - 
     +         (w(ix,iy,izp1)*e(ix,iy,izp1) -
     +          w(ix,iy,izm1)*e(ix,iy,izm1))*dzw2_i
c
	r5(ix,iy,iz)=0.25*((r5(ix,iy,iz) - u_avg(ix,iy)*ex(ix,iy))*2.0
     +        - w(ix,iy,iz)*(e(ix,iy,izp1)-e(ix,iy,izm1))*dzw3_i)
      enddo
      enddo
c
c ------------- 9/1989 add ihflt=1 option--mean shear does not generate sgs tke
c
      uxymm=0.
      uxymp=0.
      vxymm=0.
      vxymp=0.
      if(ivis .eq. 1 .and. iz .le. nmatch) then
         uxymm = u_mn(iz)
         uxymp = u_mn(izp1)
         vxymm = v_mn(iz)
         vxymp = v_mn(izp1)
      endif
c
      do iy=iys,iye
      do ix=1,nnx
c
c ----------------- dissipation 
c
         dissp(ix,iy) =  (0.19+0.74*alk(ix,iy,iz)/dslk)*
     +            e(ix,iy,iz)*sqrt(e(ix,iy,iz))/alk(ix,iy,iz)
         r5(ix,iy,iz)=r5(ix,iy,iz) - dissp(ix,iy)
c
c ----------------- vertical diffusion
c
         fnt3(ix,iy) = 
     +      ((vis_m(ix,iy,izp1)+vis_m(ix,iy,iz))*
     +       (e(ix,iy,izp1)-e(ix,iy,iz))*dzw_i(izp1) -
     +       (vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +       (e(ix,iy,iz  )-e(ix,iy,izm1))*dzw_i(iz))*dzu_i(izp1)
         r5(ix,iy,iz) = r5(ix,iy,iz) + fnt3(ix,iy)
c
c ----------------- shear production
c
         s11 = weit1*ux(ix,iy,iz)**2 + weit*ux(ix,iy,izp1)**2
         s22 = weit1*vy(ix,iy,iz)**2 + weit*vy(ix,iy,izp1)**2
         wz  = (w(ix,iy,iz)-w(ix,iy,izm1))*dzw_i(iz)
         wzp = (w(ix,iy,izp1)-w(ix,iy,iz))*dzw_i(izp1)
         s33 = weit*wzp**2 + weit1*wz**2
         s12 = weit1*(uy(ix,iy,iz) + vx(ix,iy,iz))**2 +
     +         weit*(uy(ix,iy,izp1) + vx(ix,iy,izp1))**2
         uzmn=(u(ix,iy,izp1)-uxymp-u(ix,iy,iz)+uxymm)*dzu_i(izp1) 
         vzmn=(v(ix,iy,izp1)-vxymp-v(ix,iy,iz)+vxymm)*dzu_i(izp1)
         s13 = (uzmn + wx(ix,iy,iz))**2
         s23 = (vzmn + wy(ix,iy,iz))**2
c
         fnt1(ix,iy) = vis_m(ix,iy,iz)*(2.0*(s11 + s22 + s33) +
     +                                   s13 + s23 + s12)
         r5(ix,iy,iz) = r5(ix,iy,iz) + fnt1(ix,iy)
c
c ----------------- buoyancy, get tau_w*theta
c
         buoy_sgs = -vis_s(ix,iy,1,iz)*(t(ix,iy,1,izp1) -
     +                      t(ix,iy,1,iz))*dzu_i(izp1)
                     !check if iscl is needed
         r5(ix,iy,iz) = r5(ix,iy,iz) + batag*buoy_sgs
c
         enddo
         enddo
c
c ---------------- compute shear, buoyancy, diffusion
c                  terms in SGS e eqn for printout
c            **** triz is only vertical diffusion ****
c
      if(istep .eq. 1) then
         shrz(iz)   = 0.0
         triz(iz)   = 0.0
         t_diss(iz) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            shrz(iz)   = shrz(iz) + fnt1(ix,iy)
            t_diss(iz) = t_diss(iz) + dissp(ix,iy)
            triz(iz)   = triz(iz) + fnt3(ix,iy)
         enddo
         enddo
         shrz(iz)   = shrz(iz)*fnxy
         t_diss(iz) = t_diss(iz)*fnxy
         triz(iz)   = triz(iz)*fnxy
      endif
c
c -------------- end z loop
c
      enddo
c
c --------- update tendency of sgs e from y contributions
c           pencil size (nnx,iys:iye,izs:ize)
c
      do iz=izs,ize
      do iy=iys,iye
      do ix=1,nnx
         ey(ix,iy,iz) = e(ix,iy,iz)
      enddo
      enddo
      enddo
c
      call yd_mpi(ey(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
c ------ skew symmetic advection [vde/dy + d/dy(ve)]/2
c        plus SGS diffusion contribution
c
      do iz=izs,ize
      izm1   = iz - 1
      izp1   = iz + 1
      weit   = dzw(iz)/(dzw(iz) + dzw(izp1))
      weit1  = 1.0 - weit
      do iy=iys,iye
      do ix=1,nnx
         v_avg(ix,iy)   = v(ix,iy,iz)*weit1 + v(ix,iy,izp1)*weit
         fnt2(ix,iy,iz) = e(ix,iy,iz)*v_avg(ix,iy) -
     +                    4.0*vis_m(ix,iy,iz)*ey(ix,iy,iz)
         r5(ix,iy,iz)   = r5(ix,iy,iz) - 0.5*(v_avg(ix,iy)*ey(ix,iy,iz)) 
      enddo
      enddo
      enddo
c
      call yd_mpi(fnt2(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
      do iz=izs,ize
      do iy=iys,iye
      do ix=1,nnx
         r5(ix,iy,iz) = r5(ix,iy,iz) - 0.5*fnt2(ix,iy,iz)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine lower
c
c ------ setup lower boundary condition for entire plane at (iz = 1)
c        using either businger or large formulas with wind.
c        index f(.,.,2)  indicates lower. threaded version
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
      real sfc_flx(2+nscl)
      real upars(2)
c

          
      !This was for the C-FOG setup
      if (icase.eq.3) then
      !! To change the lower surface temperature at some specified time
      if (time .gt. 3600.0) then
         tsfcc(1) = 282.0
      else
         tsfcc(1)  = 284.0
      end if
      end if

      ! FATIMA ADV. CASE SETUP
      if (icase.eq.6) then
      if (time .gt. 3600.0) then
         tsfcc(1) = 286.5
         surf_rho = surf_p/Rd/tsfcc(1)
      else
         tsfcc(1) = 287.5
         surf_rho = surf_p/Rd/tsfcc(1)
      end if
      end if 


      iz   = 1
      izm1 = iz - 1
      dz_i = dzu_i(1)
c
      do iy=iys,iye
      do ix=1,nnx
         ebc(ix,iy,2)  = amax1(e(ix,iy,iz),sml_eg)
         wbc(ix,iy,2)  = 0.0
         pbc(ix,iy,2)  = 0.0
         pbc2(ix,iy,2) = 0.0
      enddo
      enddo
c
      if(iocean .eq. 1) then
         call sufto
         do iy=iys,iye
         do ix=1,nnx
            tau13m(ix,iy) = -au13m
            tau23m(ix,iy) = -au23m
         enddo
         enddo
         do iscl=1,nscl
           do iy=iys,iye
           do ix=1,nnx
              taut3m(ix,iy,iscl) = aut3m(iscl)
           enddo
           enddo
         enddo
c
      else
c
         call suft
         fac = -utau**2/(windm*sqrt(u1xy**2 + v1xy**2))
         do iy=iys,iye
         do ix=1,nnx
            tau13m(ix,iy)=fac*(windm*(u(ix,iy,iz)+ugal-u1xy)+
     +                     wind(ix,iy)*u1xy)
            tau23m(ix,iy)=fac*(windm*(v(ix,iy,iz)-v1xy)+
     +                     wind(ix,iy)*v1xy)
         enddo
         enddo
         do iscl=1,nscl
            dnom3=t10xy(iscl)*windm
            if(dnom3 .ne. 0.) then
               dnom_i = 1.0/dnom3
               do iy=iys,iye
               do ix=1,nnx
                  taut3m(ix,iy,iscl)=aut3m(iscl)*
     +                 (windm*(t(ix,iy,iscl,iz)-t1xy(iscl))+
     +                  wind(ix,iy)*(t1xy(iscl)-tsfcc(iscl)))*dnom_i
               enddo
               enddo
            else
               do iy=iys,iye
               do ix=1,nnx
                  taut3m(ix,iy,iscl) = aut3m(iscl)
               enddo
               enddo
            endif
         enddo
c
      endif
c
c -------- partial sums of surface fluxes and mean scalar
c
      sfc_flx(1) = 0.0
      sfc_flx(2) = 0.0
      do iy=iys,iye
      do ix=1,nnx
         sfc_flx(1) = sfc_flx(1) + tau13m(ix,iy)
         sfc_flx(2) = sfc_flx(2) + tau23m(ix,iy)
      enddo
      enddo
      do iscl=1,nscl
         sfc_flx(2+iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            sfc_flx(2+iscl) = sfc_flx(2+iscl) + taut3m(ix,iy,iscl)
         enddo
         enddo
      enddo
c
      call mpi_sum_xy(sfc_flx,myid,iss,ise,(2+nscl))
      uwsfc = sfc_flx(1)*fnxy
      vwsfc = sfc_flx(2)*fnxy
      do iscl=1,nscl
         wtsfc(iscl) = sfc_flx(2+iscl)*fnxy
!      write(nprt,2345) uwsfc, vwsfc, wtsfc(iscl), tsfcc(iscl)
! 2345 format(' in lower 2345 uwsfc = ',e15.6,' vwsfc = ',e15.6,
!     +       ' wtsfc = ',e15.6,' tsfcc = ',e15.6)
      enddo
c
      do iy=iys,iye
      do ix=1,nnx
         dudz     = 2.*(u(ix,iy,iz) + ugal)*dz_i
         dvdz     = 2.*v(ix,iy,iz)*dz_i
         ubc(ix,iy,2) = u(ix,iy,iz) - dudz*dzu(iz)
         vbc(ix,iy,2) = v(ix,iy,iz) - dvdz*dzu(iz)
      enddo
      enddo
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            dtdz     = 2.*(t(ix,iy,iscl,iz)-tsfcc(iscl))*dz_i
            tbc(ix,iy,iscl,2) = t(ix,iy,iscl,iz) - dtdz*dzu(iz)
         enddo
         enddo
      enddo
c
c ------------ initialize u, v, w, t and derivatives at izm1
c
      do iy=iys,iye
      do ix=1,nnx
         u(ix,iy,izm1)  = ubc(ix,iy,2)
         v(ix,iy,izm1)  = vbc(ix,iy,2)
         w(ix,iy,izm1)  = wbc(ix,iy,2)
         r3(ix,iy,izm1) =  0.0
         e(ix,iy,izm1)  = ebc(ix,iy,2)
         ux(ix,iy,izm1) = 0.0
         uy(ix,iy,izm1) = 0.0
         vx(ix,iy,izm1) = 0.0
         vy(ix,iy,izm1) = 0.0
         wx(ix,iy,izm1) = wbc(ix,iy,2)
         wy(ix,iy,izm1) = wbc(ix,iy,2)
      enddo
      enddo
c
c ------------- no need to call derivatives here since
c               wbc = 0, change for more general lower bc
c
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,iscl,izm1) = tbc(ix,iy,iscl,2)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine lower_dns

      !Make lower BC by setting w = 0, u,v equal to mirror of interior points
      !Also set the scalar boundary condition based on either Neumann or Dirichlet

      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      implicit none

      integer :: iz,izm1,iscl,ix,iy
      real :: dz_i,dtdz,flux

      iz   = 1
      izm1 = iz - 1
      dz_i = dzu_i(1)

c ------------ initialize u, v, w, t and derivatives at izm1
c
      do iy=iys,iye
      do ix=1,nnx
         u(ix,iy,izm1)  = -2.0*Uo-u(ix,iy,iz)
         v(ix,iy,izm1)  = -v(ix,iy,iz) 
         w(ix,iy,izm1)  = 0.0 
         r3(ix,iy,izm1) = 0.0
         e(ix,iy,izm1)  = 0.0 
         ux(ix,iy,izm1) = 0.0
         uy(ix,iy,izm1) = 0.0
         vx(ix,iy,izm1) = 0.0
         vy(ix,iy,izm1) = 0.0
         wx(ix,iy,izm1) = 0.0 
         wy(ix,iy,izm1) = 0.0 
           wbc(ix,iy,2) = 0.0
           ebc(ix,iy,2) = 0.0
           ubc(ix,iy,2) = u(ix,iy,izm1)
           vbc(ix,iy,2) = v(ix,iy,izm1)
           pbc(ix,iy,2) = 0.0
           pbc2(ix,iy,2)= 0.0
         
      enddo
      enddo


!NOTE: sign convention is that wtsfc is flux upwards INTO domain (not just in vertical direction)
! Set the scalar boundary condition based on isfc
c            isfc = 0, specified surface heat flux (through qstar)
c                 = 1, specified surface temperature (Tbot)

      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            if (isfc(iscl)==1) then
            t(ix,iy,iscl,izm1) = 2.0*Tbot(iscl) - t(ix,iy,iscl,iz) 
            end if

            if (isfc(iscl)==0) then 
            flux = dzu(0)*wtsfc(iscl)/vis_s(ix,iy,iscl,izm1)
            t(ix,iy,iscl,izm1) = t(ix,iy,iscl,iz)     
     +       + flux  
            end if
         enddo
         enddo
      enddo

      end subroutine lower_dns
      subroutine upper_dns

      !Make upper BC by setting w = 0, u,v equal to mirror of interior points
      !Also set the scalar boundary condition based on either Neumann or Dirichlet

      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      implicit none

      integer :: iz,izp1,iscl,ix,iy
      real :: dz_i,dtdz,flux

      iz   = nnz
      izp1 = iz + 1

c
c ------------ initialize u, v, w, t and derivatives at izp1
c
      do iy=iys,iye
      do ix=1,nnx
         u(ix,iy,izp1)  = 2.0*Uo-u(ix,iy,iz)
         v(ix,iy,izp1)  = -v(ix,iy,iz) 
         w(ix,iy,izp1)  = 0.0 
         r3(ix,iy,izp1) = 0.0
         e(ix,iy,izp1)  = 0.0 
         ux(ix,iy,izp1) = 0.0
         uy(ix,iy,izp1) = 0.0
         vx(ix,iy,izp1) = 0.0
         vy(ix,iy,izp1) = 0.0
         wx(ix,iy,izp1) = 0.0 
         wy(ix,iy,izp1) = 0.0 
           wbc(ix,iy,1) = 0.0
           ebc(ix,iy,1) = 0.0
           ubc(ix,iy,1) = u(ix,iy,izp1)
           vbc(ix,iy,1) = v(ix,iy,izp1)
           pbc(ix,iy,1) = 0.0
           pbc2(ix,iy,1)= 0.0
      enddo
      enddo

! Set the scalar boundary condition based on isfc
c            isfc = 0, specified surface heat flux (through qstar)
c                 = 1, specified surface temperature (Ttop)

!NOTE: sign convention is that wtsfc is flux upwards INTO domain (not just in vertical direction)

      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            if (isfc(iscl)==1) then
            t(ix,iy,iscl,izp1) = 2.0*Ttop(iscl)-t(ix,iy,iscl,iz) 
            end if
            if (isfc(iscl)==0) then 
            flux = dzu(izp1)*wtsfc(iscl)/vis_s(ix,iy,iscl,izp1)
            t(ix,iy,iscl,izp1) = t(ix,iy,iscl,iz)      
     +       - flux 
            end if
         enddo
         enddo
      enddo

      end subroutine upper_dns
      subroutine lower_free
c
c --------------- setup lower boundary condition for free
c                 convection where each processor applies
c                 log-law at several (ix,iy) for iz = 1
c
c                 index f(.,.,2)  indicates lower
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      real u_level1(nnx,iys:iye,2+nscl), buf(2+2*nscl)
      real sbuf(2+2*nscl,mxs:mxe,iys:iye)
      real rbuf((2+2*nscl)*nnx*(iye+1-iys))
c
c -------------- broadcast level 1 data everywhere
c
      if(iss .eq. 0) then
         do iy=iys,iye
         do ix=1,nnx
            u_level1(ix,iy,1) = u(ix,iy,1)
            u_level1(ix,iy,2) = v(ix,iy,1)
         enddo
         enddo
         do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            u_level1(ix,iy,2+iscl) = t(ix,iy,iscl,1)
         enddo
         enddo
         enddo
      endif
      num = nnx*(iye + 1 - iys)*(2+nscl)
c
c ------ send all of root data to other processors
c
      call mpi_send_root(u_level1(1,iys,1),
     +             num,myid,numprocs,ncpu_s)
c
c --------- every task gets their own fluxes and surface scalars
c
      call suft2(u_level1)
c
c --------- send surface scalars and momentum fluxes
c           back to root(s)
c
      if(numprocs .eq. 1) go to 999
c
      do iy=iys,iye
      do ix=mxs,mxe
         sbuf(1,ix,iy)  = tau13m(ix,iy)
         sbuf(2,ix,iy)  = tau23m(ix,iy)
      enddo
      enddo
      do iscl=1,nscl
      do iy=iys,iye
      do ix=mxs,mxe
         sbuf(2+iscl,ix,iy)      = taut3m(ix,iy,iscl)
         sbuf(2+nscl+iscl,ix,iy) = t_grnd(ix,iy,iscl)
      enddo
      enddo
      enddo
c
      irow_r = mod(myid,ncpu_s)
      if(myid .ge. ncpu_s) then
        num = (2+2*nscl)*(mxe+1-mxs)*(iye+1-iys)
        call mpi_send(sbuf(1,mxs,iys),num,mpi_real8,irow_r,1,
     +       mpi_comm_world,ierr)
      else
        do l=irow_r+ncpu_s,numprocs-1,ncpu_s
           num = (2+2*nscl)*(mx_e(l)+1-mx_s(l))*(iye+1-iys)
           call mpi_recv(rbuf(1),num,mpi_real8,l,1,
     +          mpi_comm_world,istatus,ierr)
c          call f_suft2(rbuf,maxnx,maxny,mx_s(l),mx_e(l),iys,iye,nscl,
           call f_suft2(rbuf,nnx,mx_s(l),mx_e(l),iys,iye,nscl,
     +                  tau13m,tau23m,taut3m,t_grnd)
        enddo
      endif
c
  999 continue
c
c ------------ only for root row = 0
c              get sums of surface conditions
c              and set surface boundary conditions
c
      if(iss .eq. 0) then
c
         buf(1) = 0.0
         buf(2) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            buf(1) = buf(1) + tau13m(ix,iy)
            buf(2) = buf(2) + tau23m(ix,iy)
         enddo
         enddo
         do iscl=1,nscl
            buf(2+iscl)      = 0.
            buf(2+nscl+iscl) = 0.
            do iy=iys,iye
            do ix=1,nnx
               buf(2+iscl)      = buf(2+iscl) + taut3m(ix,iy,iscl)
               buf(2+nscl+iscl) = buf(2+nscl+iscl) + t_grnd(ix,iy,iscl)
            enddo
            enddo
         enddo
c
         call mpi_sum_xy(buf,myid,iss,ise,2+2*nscl)
         uwsfc = buf(1)*fnxy
         vwsfc = buf(2)*fnxy
         do iscl=1,nscl
            wtsfc(iscl) = buf(2+iscl)*fnxy
            tsfcc(iscl) = buf(2+nscl+iscl)*fnxy
         enddo
c
         iz   = 1
         izm1 = iz - 1
         dz_i = dzu_i(iz)
c
         do iy=iys,iye
         do ix=1,nnx
            ebc(ix,iy,2)=amax1(e(ix,iy,iz),sml_eg)
            wbc(ix,iy,2)= 0.0
            pbc(ix,iy,2) = 0.0
            pbc2(ix,iy,2) = 0.0
         enddo
         enddo
c
         do iy=iys,iye
         do ix=1,nnx
            dudz     = 2.*u(ix,iy,iz)*dz_i
            dvdz     = 2.*v(ix,iy,iz)*dz_i
            ubc(ix,iy,2) = u(ix,iy,iz) - dudz*dzu(iz)
            vbc(ix,iy,2) = v(ix,iy,iz) - dvdz*dzu(iz)
         enddo
         enddo
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               dtdz     = 2.*(t(ix,iy,iscl,iz)-tsfcc(iscl))*dz_i
               tbc(ix,iy,iscl,2) = t(ix,iy,iscl,iz) - dtdz*dzu(iz)
            enddo
            enddo
         enddo
c
c ------------ initialize u, v, w, t and derivatives at izm1
c
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,izm1)  = ubc(ix,iy,2)
            v(ix,iy,izm1)  = vbc(ix,iy,2)
            w(ix,iy,izm1)  = wbc(ix,iy,2)
            r3(ix,iy,izm1) =  0.0
            e(ix,iy,izm1)  = ebc(ix,iy,2)
            ux(ix,iy,izm1) = 0.0
            uy(ix,iy,izm1) = 0.0
            vx(ix,iy,izm1) = 0.0
            vy(ix,iy,izm1) = 0.0
            wx(ix,iy,izm1) = wbc(ix,iy,2)
            wy(ix,iy,izm1) = wbc(ix,iy,2)
         enddo
         enddo
c
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izm1) = tbc(ix,iy,iscl,2)
            enddo
            enddo
         enddo
c
c ----- end of if block for root row
c
      endif
c
 7999 continue
c
      return
      end
      subroutine f_suft2(rbuf,nnx,mxs,mxe,iys,iye,nscl,
     +                  tau13m,tau23m,taut3m,t_grnd)
c
c ------ fill surface arrays on root processors
c
      real rbuf(2+2*nscl,mxs:mxe,iys:iye)
      real tau13m(nnx,iys:iye), tau23m(nnx,iys:iye),
     +     taut3m(nnx,iys:iye,nscl), t_grnd(nnx,iys:iye,nscl)
c
      do iy=iys,iye
      do ix=mxs,mxe
         tau13m(ix,iy) = rbuf(1,ix,iy)
         tau23m(ix,iy) = rbuf(2,ix,iy)
      enddo
      enddo
      do iscl=1,nscl
         do iy=iys,iye
         do ix=mxs,mxe
            taut3m(ix,iy,iscl) = rbuf(2+iscl,ix,iy)
            t_grnd(ix,iy,iscl) = rbuf(2+nscl+iscl,ix,iy)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine upper
c
c ---- set boundary condition on upper boundary iz=nnz
c      option for special radiation boundary condition
c                 index f(.,.,1)  indicates upper. 
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      iz   = nnz
      izm1 = iz - 1
      izm2 = iz - 2
      izp1 = iz + 1
      izp2 = iz + 2
c
      if(ibcu .eq. 0) then
c
c --------- boundary conditions are gradient conditions
c
c            dudzbc = 0.0
c            dvdzbc = 0.0
c            dtdzbc = dtdzf
c            wbc    = 0.0
c            ebc    = 0.0
c
        do iy=iys,iye
        do ix=1,nnx
           wbc(ix,iy,1) = 0.0
           ebc(ix,iy,1) = 0.0
           ubc(ix,iy,1) = u(ix,iy,iz)
           vbc(ix,iy,1) = v(ix,iy,iz)
           pbc(ix,iy,1) = 0.0
           pbc2(ix,iy,1)= 0.0
        enddo
        enddo
        do iscl=1,nscl
c
c ---------- get average scalar gradient
c
           dtdzf(iscl) = 0.0
           do iy=iys,iye
           do ix=1,nnx
              dtdzf(iscl) = dtdzf(iscl) + (t(ix,iy,iscl,nnz) -
     +                      t(ix,iy,iscl,nnz-1))*dzu_i(nnz)
           enddo
           enddo
           dtdzf(iscl) = dtdzf(iscl)*fnxy
        enddo
c
        call mpi_sum_xy(dtdzf,myid,iss,ise,nscl)
c
        do iscl=1,nscl
           do iy=iys,iye
           do ix=1,nnx
              tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) + 
     +                            dtdzf(iscl)*dzu(nnzp1)
           enddo
           enddo
        enddo
      else if(ibcu .eq. 1) then
c
c ------------- special if iradup boundary condition
c               get estimate of w from continuity and 
c               linearized relation for pressure
c
      xmeanp = 0.0
      grad_ug = ug(nnz) - ug((nnz-1))
      do iy=iys,iye
      do ix=1,nnx
         wbc(ix,iy,1) = w(ix,iy,izm1)-
     +                  (ux(ix,iy,iz)+vy(ix,iy,iz))*dzw(iz)
         pbc(ix,iy,1) = .5*(w(ix,iy,izm1)+wbc(ix,iy,1))
         ebc(ix,iy,1) = 0.0
         ubc(ix,iy,1) = u(ix,iy,iz) + grad_ug
         vbc(ix,iy,1) = v(ix,iy,iz)
         pbc2(ix,iy,1)=0.5*(u(ix,iy,iz)**2 + v(ix,iy,iz)**2) +
     +              0.25*(w(ix,iy,izm1)**2 + wbc(ix,iy,1)**2)
         xmeanp = xmeanp + pbc2(ix,iy,1)
      enddo
      enddo
      call mpi_sum_xy(xmeanp,myid,iss,ise,1)

      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) + 
     +                          dtdzf(iscl)*dzu(nnzp1)
         enddo
         enddo
      enddo
      xmeanp = xmeanp*fnxy
      do iy=iys,iye
      do ix=1,nnx
         pbc2(ix,iy,1) = pbc2(ix,iy,1) - xmeanp
      enddo
      enddo
c
c ---------- end if block
c
      endif
c
      do iy=iys,iye
      do ix=1,nnx
         w(ix,iy,iz)   = wbc(ix,iy,1)
         e(ix,iy,iz)   = ebc(ix,iy,1)
         r3(ix,iy,iz)  = 0.0
         r5(ix,iy,iz)  = 0.0
         u(ix,iy,izp1) = ubc(ix,iy,1)
         v(ix,iy,izp1) = vbc(ix,iy,1)
c ------------- note w and e nnz+1 values are not needed
         w(ix,iy,izp1) = wbc(ix,iy,1)
         e(ix,iy,izp1) = ebc(ix,iy,1)
         r3(ix,iy,izp1)= 0.0
         r5(ix,iy,izp1)= 0.0
c
c ---------- set derivatives at top of box (wx,wy not needed)
c            ux,uy,vx,vy are used in e production, but neglect
c            at top of box becuase of bc
c
         wx(ix,iy,izp1) = 0.0
         wy(ix,iy,izp1) = 0.0
         ux(ix,iy,izp1) = 0.0
         uy(ix,iy,izp1) = 0.0
         vx(ix,iy,izp1) = 0.0
         vy(ix,iy,izp1) = 0.0
c        ux(ix,iy,izp1) = ubc(ix,iy,1)
c        uy(ix,iy,izp1) = ubc(ix,iy,1)
c        vx(ix,iy,izp1) = vbc(ix,iy,1)
c        vy(ix,iy,izp1) = vbc(ix,iy,1)
      enddo
      enddo
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,iscl,izp1) = tbc(ix,iy,iscl,1)
            t(ix,iy,iscl,izp2) = tbc(ix,iy,iscl,1)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine comp_p
c
c --------- setup pressure solver
c
      use pars
      use fftwk
      use fields
      use con_data
      use con_stats
      include 'mpif.h'
      real fnt1(nnx,iys:iye,izs:ize)
      real fs(nnx,iys:iye,2), fr(nnx,iys:iye,2)
      integer istatus(mpi_status_size)
c
      gami = 1.0/dtgama
c
      nb = myid - ncpu_s
      nt = myid + ncpu_s
c
c ------------ Send both r3 and updated w (from comp1)
c              to processor above the current myid.
c
      if(iss .eq. 0) then
         nb = mpi_proc_null
      endif
      if(ise .eq. numprocs-1) then
         nt = mpi_proc_null
      endif
      nsend = 2*nnx*(iye + 1 - iys)
      nrecv = nsend
      do iy=iys,iye
      do ix=1,nnx
         fs(ix,iy,1) = r3(ix,iy,ize)
         fs(ix,iy,2) = w(ix,iy,ize)
      enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nt,2,
     +     fr(1,iys,1),nrecv,mpi_real8,nb,2,
     +     mpi_comm_world,istatus,ierr)
      if(iss .ne. 0) then
         do iy=iys,iye
         do ix=1,nnx
            r3(ix,iy,izs-1) = fr(ix,iy,1)
            w(ix,iy,izs-1)  = fr(ix,iy,2)
         enddo
         enddo
      endif
c
c ----------- setup general pressure calculation
c             relies on rhs from step n-1 being included 
c             in velocity-arrays already
c
      do iz=izs,ize
         izm1 = iz -1
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = u(ix,iy,iz)*gami + r1(ix,iy,iz)
         enddo
         enddo
         call xderivp(fnt1(1,iys,iz),trigx(1,1),xk(1),nnx,iys,iye)
c
         if(iz .eq. 1) then
            do iy=iys,iye
            do ix=1,nnx
                p(ix,iy,iz) = fnt1(ix,iy,iz) +  
     +                     ((w(ix,iy,iz) -wbc(ix,iy,2))*gami +
     +                       r3(ix,iy,iz))*dzw_i(iz)
            enddo
            enddo
         else if(iz .eq. nnz) then
            do iy=iys,iye
            do ix=1,nnx
                p(ix,iy,iz) = fnt1(ix,iy,iz) + 
     +                     ((wbc(ix,iy,1) - w(ix,iy,izm1))*gami -
     +                      r3(ix,iy,izm1))*dzw_i(iz)
            enddo
            enddo
         else 
            do iy=iys,iye
            do ix=1,nnx
                p(ix,iy,iz) = fnt1(ix,iy,iz) + 
     +                    ((w(ix,iy,iz)  - w(ix,iy,izm1))*gami +
     +                      r3(ix,iy,iz) - r3(ix,iy,izm1))*dzw_i(iz)
            enddo
            enddo
         endif
c
c --------- end z loop
c
      enddo
c
c ----------- check for radiation boundary condition, all processors
c
      if(ibcu .eq. 1) then
        do iy=iys,iye
        do ix=1,nnx
           ptop(ix,iy,1) = pbc(ix,iy,1)
           ptop(ix,iy,2) = pbc2(ix,iy,1)
        enddo
        enddo
      endif
c
c --------- now y contribution
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = v(ix,iy,iz)*gami + r2(ix,iy,iz)
         enddo
         enddo
      enddo
c
      call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
             p(ix,iy,iz) = p(ix,iy,iz) + fnt1(ix,iy,iz) 
         enddo
         enddo
      enddo
c
      call pressure
c
      return
      end
      subroutine comp2
c
c ------- add p gradients to rhs. Use already defined p
c         at ize+1 to get w (see sr. pressure).
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real fnt1(nnx,iys:iye,izs:ize), fnt2(nnx,iys:iye)
      real r3_sum(1:nnz)
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      do iz=1,nnz
         r3_sum(iz) = 0.0
      enddo
c
c --------- dp/dy at all z
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = p(ix,iy,iz)
         enddo
         enddo
      enddo
c
      call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
      do iz=izs,ize
c
         izm1  = iz - 1
         izp1  = iz + 1
c
         do iy=iys,iye
         do ix=1,nnx
            fnt2(ix,iy) = p(ix,iy,iz)
         enddo
         enddo
         call xderivp(fnt2(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         do iy=iys,iye
         do ix=1,nnx
            r1(ix,iy,iz) = r1(ix,iy,iz) - fnt2(ix,iy)
            r2(ix,iy,iz) = r2(ix,iy,iz) - fnt1(ix,iy,iz)
         enddo
         enddo
         if (iz.ne.nnz) then
            do iy=iys,iye
            do ix=1,nnx
               r3(ix,iy,iz) = r3(ix,iy,iz) -
     +            (p(ix,iy,izp1)-p(ix,iy,iz))*dzu_i(izp1)
               r3_sum(iz) = r3_sum(iz) + r3(ix,iy,iz)
            enddo
            enddo
            r3_sum(iz) = r3_sum(iz)*fnxy
         endif
c
c ------------------------ time stepping with 3-order rk method
c                          first w variables
c
      if(iz .ne. nnz) then
         do iy=iys,iye
         do ix=1,nnx
c           w(ix,iy,iz)  = w(ix,iy,iz)+dtgama*r3(ix,iy,iz)
            e(ix,iy,iz)  = e(ix,iy,iz)+dtgama*r5(ix,iy,iz)
         enddo
         enddo
      else
c
c --------- update wout and eout by setting = to bc values
c
         do iy=iys,iye
         do ix=1,nnx
            w(ix,iy,iz)  = wbc(ix,iy,1)
            e(ix,iy,iz)  = ebc(ix,iy,1)
            r3(ix,iy,iz) = 0.0
            r5(ix,iy,iz) = 0.0
         enddo
         enddo
      endif
c
c -------- now all u-variables
c
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = u(ix,iy,iz)+dtgama*r1(ix,iy,iz)
            v(ix,iy,iz) = v(ix,iy,iz)+dtgama*r2(ix,iy,iz)
         enddo
         enddo
         do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,iscl,iz)  = t(ix,iy,iscl,iz)+
     +                          dtgama*r4(ix,iy,iscl,iz)
         enddo
         enddo
         enddo
c
c -------- end z loop
c
      enddo
c
c ---------- gather partial sums for w computation
c
      call mpi_sum_z(r3_sum,i_root,myid,nnz,1)
c
      do iz=izs,min(ize,nnz-1)
         do iy=iys,iye
         do ix=1,nnx
            r3(ix,iy,iz) = r3(ix,iy,iz) - r3_sum(iz)
            w(ix,iy,iz)  = w(ix,iy,iz) + dtgama*r3(ix,iy,iz)

         enddo
         enddo
      enddo
c
      return
      end
      subroutine pressure
c
c -------- solve for pressure using a matrix transpose
c          across mpi tasks and tridiagonal solver. 
c          The transposed array
c          is dimensioned (0:nnz+1). Values 
c          (0 & nnz+1) are not needed but are useful in the 
c          matrix transpose when we return (see send_ztox).
c          On exit p is defined at all [izs-1:ize+1].
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real pfft(nny,jxs:jxe,izs-1:ize+1)
      real pt(0:nnz+1,jxs:jxe,jys:jye)
      real ptopfft(nny,jxs:jxe,1:2)
      real psum(1:nnz)
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
c ------------ Fourier analyze the right hand side
c              at all iz = izs,ize. results are in pfft
c
c


      call fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)

c
c ------------ Fourier analyze the radiation bc arrays
c
      if(ibcu .eq. 1) then
        call fft2d_mpi(ptop(1,iys,1),ptopfft(1,jxs,1),
     +           trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           1,2,myid,ncpu_s,numprocs,-2)
      endif
c
c ---------- transpose first and last index of array
c            the order of pfft is (y,x,z)
c
      call xtoz_trans(pfft,pt,nny,nnz,jys,jye,jy_s,jy_e,
     +                jxs,jxe,izs,ize,iz_s,iz_e,myid,ncpu_s,
     +                numprocs)


      call solve_trid(pt, ptopfft)

c
c ------------- transpose back
c
      call ztox_trans(pt,pfft,nny,nnz,jys,jye,jy_s,jy_e,
     +                jxs,jxe,izs,ize,iz_s,iz_e,myid,ncpu_s,
     +                numprocs)
c
      iz_ee = ize+1
      if(ise .eq. numprocs-1) then
         iz_ee = ize
      endif
c
c --------- inverse fft at all iz=izs,iz_ee to get p
c           see z indices
c
      call fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,iz_ee,myid,ncpu_s,numprocs,2)
c
c -------- partial sums for pressure
c
      do iz=1,nnz
         psum(iz) = 0.0
      enddo
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            psum(iz) = psum(iz) + p(ix,iy,iz)
         enddo
         enddo
         psum(iz) = psum(iz)*fnxy
      enddo
      call mpi_sum_z(psum,i_root,myid,nnz,1)
c
      do iz=izs,iz_ee
c        psum(iz) = -psum(iz) + engz(iz) + c23*engsbz(iz)
         do iy=iys,iye
         do ix=1,nnx
            p(ix,iy,iz) = p(ix,iy,iz) - psum(iz)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine solve_trid(pt, ptop)
c 
c --------- tridiagonal solver. odd order for ptop, ptop2
c           because of 2d-fft
c
      use pars
      use con_data
      use con_stats
c
      real ptop(nny,jxs:jxe,1:2)
      real pt(0:nnz+1,jxs:jxe,jys:jye)
      real aa(nnz,jxs:jxe),bb(nnz,jxs:jxe),
     +     dd(nnz,jxs:jxe),rh(nnz,jxs:jxe)
      real fac_u(nnz), fac_l(nnz), fac_a(nnz)
c
      do iz=1,nnz
         fac_u(iz) = 1.0/(dzw(iz)*dzu(iz+1))
         fac_l(iz) = 1.0/(dzw(iz)*dzu(iz))
         fac_a(iz) = fac_l(iz) + fac_u(iz)
      enddo
c
      do kp=jys,jye    
         do lp=jxs,jxe
         do iz=2,nnz-1
            bb(iz,lp)  = fac_l(iz)
            aa(iz,lp)  = fac_u(iz)
            dd(iz,lp)  = -xks(lp,kp) - fac_a(iz)
            rh(iz,lp)  = pt(iz,lp,kp)
         enddo
         enddo
c
c --------------- lower boundary, fill exterior pressure (not used)
c
         do lp=jxs,jxe
            bb(1,lp)  = 1.0
            aa(1,lp)  = fac_u(1)
            dd(1,lp)  = -xks(lp,kp) - fac_u(1)
            rh(1,lp)  = pt(1,lp,kp)
            pt(0,lp,kp) = 0.0
         enddo
c
c --------------- upper boundary, fill exterior pressure (not used)
c
         if(ibcu .eq. 1) then
            do lp=jxs,jxe
              bb(nnz,lp) = 0.0
              aa(nnz,lp) = 0.0
              dd(nnz,lp) = 1.0
              rh(nnz,lp) = ptop(kp,lp,1)*wavexy(lp,kp) + ptop(kp,lp,2)
              pt(nnz+1,lp,kp) = 0.0
            enddo
         else
            do lp=jxs,jxe
               bb(nnz,lp) = fac_l(nnz)
               aa(nnz,lp) = 1.0
               dd(nnz,lp) = -xks(lp,kp) - fac_l(nnz)
               rh(nnz,lp) = pt(nnz,lp,kp)
               pt(nnz+1,lp,kp) = 0.0
            enddo
         endif
c
c ---------------- special situation for zeroth mode
c                  makes mean pressure = 0
c
         if(kp .eq. 1 .and. jxs .eq. 1) then
           do iz=1,nnz
              dd(iz,1) = 1.0
              rh(iz,1) = 0.0
              aa(iz,1) = 0.0
              bb(iz,1) = 0.0
              dd(iz,2) = 1.0
              rh(iz,2) = 0.0
              aa(iz,2) = 0.0
              bb(iz,2) = 0.0
           enddo
         endif
c
c --------------- solve system
c


         call tridv(bb,dd,aa,rh,nnz,jxs,jxe)
         do lp=jxs,jxe
         do iz=1,nnz
            pt(iz,lp,kp) = rh(iz,lp)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine tridv(b,d,a,r,n,j1,j2)
c
c --- tridiagonal matrix solver with multiple vectors
c     (note j and i loops are reversed from cray version)
c
c --- input:   n   size of a,b,d and r
c              b   below diagonal elements (b(1) not used)
c              d   diagonal elements
c              a   above diagonal elements (a(n) not used)
c              r   right hand side
c              j1:j2  range of input vectors
c
c --- output:  r   solution vector
c
      real b(n,j1:j2), d(n,j1:j2), a(n,j1:j2), r(n,j1:j2)
c

      if(n .le. 1 ) then
         do j=j1,j2
            r(1,j) = r(1,j)/d(1,j)
         enddo
         go to 999
      endif
      do j=j1,j2
         d(1,j) = 1.0/d(1,j)
      enddo
      do j=j1,j2
      do i=2,n
         fac = b(i,j)*d(i-1,j)
         d(i,j) = 1.0/(d(i,j) - fac*a(i-1,j))
         r(i,j) = r(i,j) - fac*r(i-1,j)
      enddo
      enddo
      do j=j1,j2
         r(n,j) = r(n,j)*d(n,j)
      enddo
      do j=j1,j2
      do i=n-1,1,-1
         r(i,j) = d(i,j)*(r(i,j) - a(i,j)*r(i+1,j))
      enddo
      enddo
  999 continue
c
      return
      end
      subroutine get_derv
c
c ------- get ux,uy,vx,vy at all z for this node
c         using parallel fft. can be improved (?)
c         by using exchange to send derivatives
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      iz_ss = izs-1
      iz_ee = ize+1
      if(iss .eq. 0) then
         iz_ss = izs 
      endif
      if(ise .eq. numprocs-1) then
         iz_ee = ize
      endif
c
c ------- make sure <w> = 0
c
      do iz=izs-1,ize+1
         w_sum = 0.0
         do iy=iys,iye
         do ix=1,nnx
            w_sum = w_sum + w(ix,iy,iz)
         enddo
         enddo
         w_sum = w_sum*fnxy
         call mpi_sum_xy(w_sum,myid,iss,ise,1)
         do iy=iys,iye
         do ix=1,nnx
            w(ix,iy,iz) = w(ix,iy,iz) - w_sum
         enddo
         enddo
      enddo
c
c     do iz=iz_ss,iz_ee
      do iz=izs-1,ize+1
c        if(iz .eq. izs-1 .or. iz .eq. ize+1) then
c        do iy=iys,iye
c        do ix=1,nnx
c           ux(ix,iy,iz) = 0.0
c           vx(ix,iy,iz) = 0.0
c           wx(ix,iy,iz) = 0.0
c           uy(ix,iy,iz) = 0.0
c           vy(ix,iy,iz) = 0.0
c           wy(ix,iy,iz) = 0.0
c        enddo
c        enddo
c        else
         do iy=iys,iye
         do ix=1,nnx
            ux(ix,iy,iz) = u(ix,iy,iz)
            vx(ix,iy,iz) = v(ix,iy,iz)
            wx(ix,iy,iz) = w(ix,iy,iz)
            uy(ix,iy,iz) = u(ix,iy,iz)
            vy(ix,iy,iz) = v(ix,iy,iz)
            wy(ix,iy,iz) = w(ix,iy,iz)
         enddo
         enddo
c        endif
         call xderivp(ux(1,iys,iz),trigx(1,1),xk(1),
     +                 nnx,iys,iye)
         call xderivp(vx(1,iys,iz),trigx(1,1),xk(1),
     +                 nnx,iys,iye)
         call xderivp(wx(1,iys,iz),trigx(1,1),xk(1),
     +                 nnx,iys,iye)
      enddo
c
c ---------- get y derivatives for (u,v,w)
c
c     call yd_mpi(uy(1,iys,iz_ss),trigx(1,2),yk(1),
c    +           nnx,nny,ixs,ixe,ix_s,ix_e,
c    +           iys,iye,iy_s,iy_e,iz_ss,iz_ee,myid,ncpu_s,numprocs)
c     call yd_mpi(vy(1,iys,iz_ss),trigx(1,2),yk(1),
c    +           nnx,nny,ixs,ixe,ix_s,ix_e,
c    +           iys,iye,iy_s,iy_e,iz_ss,iz_ee,myid,ncpu_s,numprocs)
c     call yd_mpi(wy(1,iys,iz_ss),trigx(1,2),yk(1),
c    +           nnx,nny,ixs,ixe,ix_s,ix_e,
c    +           iys,iye,iy_s,iy_e,iz_ss,iz_ee,myid,ncpu_s,numprocs)
c
      call yd_mpi(uy(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
      call yd_mpi(vy(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
      call yd_mpi(wy(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
c
      return
      end
      subroutine get_means(istage)
c
c ------------ get means for all variables
c              for use in iso, surfvis, comp1, compmn.
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
c
c -------- set e to minimum value 
c
      do iz=izs-1,ize+1
         do iy=iys,iye
         do ix=1,nnx
            e(ix,iy,iz )=amax1(e(ix,iy,iz ),sml_eg)
         enddo
         enddo
      enddo

      do iz=0,nnz+1
         u_mn(iz)   = 0.0
         v_mn(iz)   = 0.0
         w_mn(iz)   = 0.0
         e_mn(iz)   = 0.0
         engz(iz)   = 0.0
         engsbz(iz) = 0.0
         divz(iz)   = 0.0
         pxym(iz)   = 0.0
      enddo
      do iscl=1,nscl
         do iz=0,nnz+1
            t_mn(iz,iscl) = 0.0
            alphaC(iz,iscl) = 0.0
         enddo
      enddo
      iz_ee = ize
      iz_ss = izs
      if(ize .eq. nnz) iz_ee = nnzp1
      if(izs .eq. 1)   iz_ss = 0 
      do iz=iz_ss,iz_ee
         do iy=iys,iye
         do ix=1,nnx
            u_mn(iz) = u_mn(iz) + u(ix,iy,iz)
            v_mn(iz) = v_mn(iz) + v(ix,iy,iz)
            w_mn(iz) = w_mn(iz) + w(ix,iy,iz)
            e_mn(iz) = e_mn(iz) + e(ix,iy,iz)

         enddo
         enddo
         u_mn(iz) = u_mn(iz)*fnxy
         v_mn(iz) = v_mn(iz)*fnxy
         w_mn(iz) = w_mn(iz)*fnxy
         e_mn(iz) = e_mn(iz)*fnxy
         do iscl=1,nscl
            t_mn(iz,iscl) = 0.0
            alphaC(iz,iscl) = 0.0
            do iy=iys,iye
            do ix=1,nnx
               t_mn(iz,iscl) = t_mn(iz,iscl) + t(ix,iy,iscl,iz)
               alphaC(iz,iscl) = alphaC(iz,iscl)+vis_s(ix,iy,iscl,iz)
            enddo
            enddo
            t_mn(iz,iscl) = t_mn(iz,iscl)*fnxy
            alphaC(iz,iscl) = alphaC(iz,iscl)*fnxy
         enddo
      enddo
      call mpi_sum_z(u_mn(0),i_root,myid,nnzp1+1,1)
      call mpi_sum_z(v_mn(0),i_root,myid,nnzp1+1,1)
      call mpi_sum_z(w_mn(0),i_root,myid,nnzp1+1,1)
      call mpi_sum_z(e_mn(0),i_root,myid,nnzp1+1,1)
      do iscl=1,nscl
         call mpi_sum_z(t_mn(0,iscl),i_root,myid,nnzp1+1,1)
         call mpi_sum_z(alphaC(0,iscl),i_root,myid,nnzp1+1,1)
      enddo
c
c -------------- get terms which contribute to mean pressure
c                careful with the sum, get the mean p_star pressure
c
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            e_temp     =  0.5*(e(ix,iy,iz) + e(ix,iy,izm1))
            q_temp     =  0.5*((u(ix,iy,iz) + stokes(iz))**2 +
     +                          v(ix,iy,iz)*v(ix,iy,iz) +
     +                       0.5*(w(ix,iy,iz)*w(ix,iy,iz) +
     +                            w(ix,iy,izm1)*w(ix,iy,izm1)))
            engz(iz)   = engz(iz) + q_temp
            engsbz(iz) = engsbz(iz) + e_temp
            pxym(iz)   = pxym(iz) + p(ix,iy,iz) - (c23*e_temp + q_temp)
         enddo
         enddo
         engz(iz)   = engz(iz)*fnxy
         engsbz(iz) = engsbz(iz)*fnxy
         pxym(iz)   = pxym(iz)*fnxy
      enddo
      call mpi_sum_z(engz(1),i_root,myid,nnzp1,1)
      call mpi_sum_z(engsbz(1),i_root,myid,nnzp1,1)
      call mpi_sum_z(pxym(1),i_root,myid,nnz,1)
c
c ------------ save means and divergence for printout and compmn
c              all cpus have means over all z
c
      if(istage .eq. 1) then
        do iz=izs,ize
           izm1 = iz - 1
           do iy=iys,iye
           do ix=1,nnx
              divz(iz) = divz(iz) + 
     +                  (ux(ix,iy,iz)+vy(ix,iy,iz)+
     +                  (w(ix,iy,iz)-w(ix,iy,izm1))*dzw_i(iz))**2
           enddo
           enddo
           divz(iz) = divz(iz)*fnxy
        enddo
        call mpi_sum_z(divz(1),i_root,myid,nnz,1)
c
        do iz=0,nnz+1
           uxym(iz) = u_mn(iz)
           vxym(iz) = v_mn(iz)
           wxym(iz) = w_mn(iz)
        enddo
        do iscl=1,nscl
           do iz=0,nnz+1
              txym(iz,iscl) = t_mn(iz,iscl)
           enddo
        enddo
      endif
c
      return
      end
      subroutine xderivp(ax,trigx,xk,nnx,iys,iye)
c
c -------- get multiple x derivatives using fftpack routines
c          use fftpack storage a0, (a1,b1), (a2,b2),...,an
c          assumes that first wavenumber xk(1) = 0.0
c
c          assumes that wavenumbers are normalized by number of points
c
      real xk(nnx), trigx(2*nnx+15), ax(nnx,iys:iye)
c
c     fn = 1.0/float(nnx)
      do iy=iys,iye
         call rfftf(nnx,ax(1,iy),trigx)
         ii = 1
         ax(1,iy) = 0.0
         ax(nnx,iy) = 0.0
         do ix=2,nnx-1,2
            ii          = ii + 1
            temp        = ax(ix,iy)
            ax(ix,iy)   = -xk(ii)*ax(ix+1,iy)
            ax(ix+1,iy) = xk(ii)*temp
         enddo
         call rfftb(nnx,ax(1,iy),trigx)
      enddo
c
      return
      end
      subroutine fft2d_mpi(ax,at,trigx,trigc,nx,ny,
     +           jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           iz1,iz2,myid,ncpu,np,isgn)
c
c -------- get 2d fft using fftpack routines and parallel mpi
c          use fftpack storage a0, (a1,b1), (a2,b2),...,
c
c         isgn = -1 do forward transform, get coefficients
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
c
c         isgn = -2 do forward transform, get coefficients
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is at(ny,jxs:jxe,iz1:iz2)
c
c         isgn =  1 do inverse transform, move to physical space
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
c
c         isgn =  2 do inverse transform, move to physical space
c                   incoming array is at(ny,jxs:jxe,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
c
      real ax(nx+2,iys:iye,iz1:iz2), at(ny,jxs:jxe,iz1:iz2),
     +     trigx(2*nx+15), trigc(4*ny+15),
     +     a2d(2,ny), a_wrk(nx)
      integer jx_s(0:np-1), jx_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      nxp2 = nx + 2
      if(isgn .lt. 0) then
         fn   = 1.0/(float(nx)*float(ny))
c
c ------ 1d fft in x over [iys,iye] for all z
c
         do iz=iz1,iz2
            do iy=iys,iye
               do ix=1,nx
                  a_wrk(ix) = ax(ix,iy,iz)*fn
               enddo
               call rfftf(nx,a_wrk(1),trigx(1))
               ax(1,iy,iz) = a_wrk(1)
               ax(2,iy,iz) = 0.0
               do ix=2,nx
                  ax(ix+1,iy,iz) = a_wrk(ix)
               enddo
               ax(nx+2,iy,iz) = 0.0
            enddo
         enddo
         call xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +        iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c ------ 1d fft in y over [jxs,jxe] for all z
c
         do iz=iz1,iz2
            do ix=jxs,jxe,2
               do iy=1,ny
                  a2d(1,iy) = at(iy,ix,iz)
                  a2d(2,iy) = at(iy,ix+1,iz)
               enddo
               call cfftf(ny,a2d(1,1),trigc(1))
               do iy=1,ny
                  at(iy,ix,iz)   = a2d(1,iy)
                  at(iy,ix+1,iz) = a2d(2,iy)
               enddo
            enddo
         enddo
c
c ---- decide whether to transpose back or leave as is
c
         if(isgn .eq. -1) then
            call ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
         endif
c
      else
c
c ---- decide whether to first transpose or leave as is
c
         if(isgn .eq. 1) then
            call xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
         endif
c
c ------ 1d fft in y over [jxs,jxe] for all z
c
         do iz=iz1,iz2
            do ix=jxs,jxe,2
               do iy=1,ny
                  a2d(1,iy) = at(iy,ix,iz)
                  a2d(2,iy) = at(iy,ix+1,iz)
               enddo
               call cfftb(ny,a2d(1,1),trigc(1))
               do iy=1,ny
                  at(iy,ix,iz)   = a2d(1,iy)
                  at(iy,ix+1,iz) = a2d(2,iy)
               enddo
            enddo
         enddo
         call ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +        iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c ------  1d fft in x over [iys,iye] for all z
c
         do iz=iz1,iz2
            do iy=iys,iye
               a_wrk(1) = ax(1,iy,iz)
               do ix=2,nx
                  a_wrk(ix) = ax(ix+1,iy,iz)
               enddo
               call rfftb(nx,a_wrk(1),trigx(1))
               do ix=1,nx
                  ax(ix,iy,iz) = a_wrk(ix)
               enddo
            enddo
         enddo
      endif
c
      return
      end
      subroutine yderiv(ay,trigy,yk,nnx,nny)
c
c -------- get multiple y derivatives using fftpack routines
c          use fftpack storage a_0, (a1,b1), (a2,b2), ...,
c          assumes that first wavenumber yk(1) = 0.0
c
c          assumes that wavenumbers are normalized by number of points
c
      real yk(nny), trigy(2*nny+15), ay(nnx,nny)
      real a_trans(nny)
c
c     fn = 1.0/float(nny)
      do ix=1,nnx
         do iy=1,nny
            a_trans(iy) = ay(ix,iy)
         enddo
         call rfftf(nny,a_trans(1),trigy)
         ii = 1
         a_trans(1)   = 0.0
         a_trans(nny) = 0.0
         do iy=2,nny-1,2
            ii            = ii + 1
            temp          = a_trans(iy)
            a_trans(iy)   = -yk(ii)*a_trans(iy+1)
            a_trans(iy+1) = yk(ii)*temp
         enddo
         call rfftb(nny,a_trans(1),trigy)
         do iy=1,nny
            ay(ix,iy) = a_trans(iy)
         enddo
      enddo
c
      return
      end
      subroutine yd_mpi(ay,trigy,yk,
     +           nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c -------- get multiple y derivatives using fftpack routines and mpi
c          use fftpack storage a_0, (a1,b1), (a2,b2), ...,
c          assumes that first wavenumber yk(1) = 0.0
c          wavenumbers are normalized by number of points, ny
c
      real yk(ny), trigy(2*ny+15), ay(nx,iys:iye,iz1:iz2)
      real ayt(ny,ixs:ixe,iz1:iz2)
c
      integer ix_s(0:np-1), ix_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      call xtoy_trans(ay,ayt,nx,ny,ixs,ixe,ix_s,ix_e,
     +         iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c     fn = 1.0/float(nny)
      do iz=iz1,iz2
         do ix=ixs,ixe
            call rfftf(ny,ayt(1,ix,iz),trigy)
            ii = 1
            ayt(1,ix,iz)  = 0.0
            ayt(ny,ix,iz) = 0.0
            do iy=2,ny-1,2
               ii              = ii + 1
               temp            = ayt(iy,ix,iz)
               ayt(iy,ix,iz)   = -yk(ii)*ayt(iy+1,ix,iz)
               ayt(iy+1,ix,iz) = yk(ii)*temp
            enddo
            call rfftb(ny,ayt(1,ix,iz),trigy)
         enddo
      enddo
      call ytox_trans(ayt,ay,nx,ny,ixs,ixe,ix_s,ix_e,
     +         iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
      return
      end
      function rlim(d1,d2,d3)
c
c ------------- Cees's kappa=1/3 scheme
c
      r = (d1-d2+1.e-100)/(d2-d3-1.e-100)
      rlim = (d2-d3)*amax1(0.,amin1(r,amin1(1./6.+1./3.*r,1.)))
c
c ------------- Cees's kappa=-1 scheme
c
c     r = (d1-d2+1.e-100)/(d2-d3-1.e-100)
c     rlim = (d2-d3)*amin1(abs(r),0.5)
c
c ------------- first order upwind
c
c     rlim = 0.0
c
c ------------- QUICK scheme
c
c     rlim = -0.25*d2 - 0.125*d3 + 0.375*d1
c
      return
      end
      function ran1(idum)
c
c ----------- stolen from numerical recipes,p. 271
c
      integer idum, ia, im, iq, ir, ntab, ndiv
      real ran1, am, eps, rnmx
      parameter (ia=16807,im=2147483647,am=1.0/im,iq=127773,ir=2836.0,
     +           ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-07,rnmx=1.0-eps)
      integer j, k, iv(ntab), iy
      save iv, iy
      data iv /ntab*0/, iy /0/
      if(idum .le. 0 .or. iy .eq. 0) then
         idum = max(-idum,1)
         do j=ntab+8,1,-1
            k = idum/iq
            idum = ia*(idum - k*iq) - ir*k
            if(idum .lt. 0) idum = idum + im
            if(j .le. ntab) iv(j) = idum
         enddo
         iy = iv(1)
      endif
      k     = idum/iq
      idum  = ia*(idum - k*iq) - ir*k
      if(idum .lt. 0) idum = idum + im
      j     = 1 + iy/ndiv
      iy    = iv(j)
      iv(j) = idum
      ran1  = min(am*iy, rnmx)
c
      return
      end
      function ranf()
      data inc /1/
      save inc, ix, ia, m, fm
      if(inc.eq.1) then
        inc = 2
        m = 2**20
        fm = float(m)
        ix = 566387
        ia = 2**10 + 3
      endif
      ix = mod(ia*ix,m)
      fx = float(ix)
      ranf = fx/fm
      return
      end
      subroutine stokesv
c
c ----------- get stokes drift velocity for assumed wavelength stokesw
c             and wave amplitude stokesa. Changed sign for z.
c
c
      use pars
      use con_data
      use con_stats
      include 'mpif.h'
c
      if(iocean .eq. 1) then
c
c ----------- compute stokes velocity for ocean pbls
c
c        stokesw = pi2/20.0
         stokesw = pi2/76.5
c        ak      = 0.04
         ak      = 0.00
c        stokesa = 1.0
         stokesa = ak/stokesw
         sigma = sqrt(abs(grav)*stokesw)
         stokess = sigma*stokesw*stokesa**2
         do iz=1,nnzp1
            stokes(iz) = stokess*exp(2.0*stokesw*zz(iz))
         enddo
         if(l_root) then
            write(6,6000) (iz,zz(iz),stokes(iz),iz=1,nnz)
 6000       format(' iz ',10x,' zz',10x,' stokes',/,(1x,i3,2e12.4))
         endif
c
      else
c
c ----------------- set stokes velocity = 0 for atmos. pbls 
c
         do iz=1,nnzp1
            stokes(iz) = 0.0
         enddo
         stokess = 0.0
         udrift = 0.0
         vdrift = 0.0
      endif
c
      return
      end
      subroutine busngr(zeta,phim,phis,psim,psis)
c
c ---- Businger's version of similarity theory
c
      data pih /1.57079633/
      save pih
c
      if(zeta .lt. 0.) then
         x=(1.0 - 15.0*zeta)**0.25
         phim = 1.0/x
         psim = 2.0*alog((1.0+x)/2.0) + alog((1.0+x*x)/2.0) - 
     +          2.0*atan(x)+pih
         if(psim.gt.2.0)psim=2.0
         y = sqrt(1.0-9.0*zeta)
         phis = 0.74/y
         psis = alog((1.0+y)/2.0)*2.0
      else if(zeta .gt. 0) then
         phim = 1.0 + 4.7*zeta
         phis = 0.74 + 4.7*zeta
         psim = -4.7*zeta
         psis = -4.7*zeta
      else
         phim = 1.0
         phis = 0.74
         psim = 0.0
         psis = 0.0
      endif
      return
      end
      subroutine fzol(zeta,phim,phis,psim,psis)
c        estimate the stability functions for momentum, m
c                                         and scalars,  c
c        from input of the stability parameter zeta = z/L

      data c1/5./
      data a3,b3,a4,b4/1.258,8.382,-28.862,98.9545/
      data zetam,zetas/-0.2,-1.0/
      save c1, a3, b3, a4, b4, zetam, zetas
c
      psimu(Y)  = 1.571 + 2.0*(alog(0.5*(1.0 + Y)) - atan(Y)) + 
     +            alog(0.5 + 0.5*Y**2)
      psisu(Y)  = 2.0*alog(0.5 + 0.5*Y)
      psicu(Y,G)= (1.0 - G)*alog(abs(Y - 1.0))
     +          + 0.5*(G + 2.0)*alog(abs(Y**2 + Y + 1.0))
     +          - (2.0*G + 1.0) / sqrt(3.0) * 
     +            atan((Y + 0.5)*2.0/sqrt(3.0))
      Xm(zol)   = (1.0 - 16.0*zol)**0.25
      Xs(zol)   = sqrt(1.0 - 16.0*zol)
      Xc(zol,f) =  abs(1.0 - f*zol)**(4.0/3.0)/(1.0 - f*zol)
c
      if(zeta.ge.0.0)       then
c                                          STABLE
      if(zeta.le.1.0) then
        phim = 1.0 + c1 * zeta
        psim = - c1 * zeta
        phis = phim
        psis = psim
                      else
c                                   use limiting form
        phim = c1 + zeta
        psim = (1.0 - c1)*(1.0 + alog(zeta) ) - zeta
        phis = phim
        psis = psim
                      endif

                            else
c                                         UNSTABLE
c                                                  momentum         
       if(zeta.ge.zetam) then
         phim = 1.0 / Xm(zeta)
         psim = psimu(Xm(zeta))
                         else
c                            use convective limit for momentum
         X = (1.0 - b3/a3 * zeta)**(1.0/3.0)

         fm = a3**(-1.0/3.0)
         phim = fm / Xc(zeta,b3/a3) 
         psim = psimu(Xm(zetam))
     *        + psicu(Xc(zeta,b3/a3),fm)
     *        - psicu(Xc(zetam,b3/a3),fm)
                         endif
      
c                                         UNSTABLE scalars
       if(zeta.ge.zetas) then
         phis = 1.0/Xs(zeta)
         psis = psisu(Xs(zeta))
                         else
c                              use convective limit for scalars
         fs =   abs(a4)**(-1.0/3.0)*abs(a4)/a4
         phis = (a4 - b4*zeta)**(-1.0/3.0)
         psis = psisu(Xs(zetas))
     *        + psicu(Xc(zeta,b4/a4),fs)
     *        - psicu(Xc(zetas,b4/a4),fs)
                         endif
               
                            endif
       return
       end
      subroutine suft
c
c ---------- iterate for zeta = z/L using bisection method
c            either businger or large functions can be specified
c
c            isfc = 0, specified surface heat flux
c                 = 1, specified surface temperature
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
      implicit none

      real, parameter :: zeta_min=-6.0, zeta_max=3.0
      integer, parameter :: iter_mo = 30
 
      integer :: iz,izp1,izm1,ix,iy,iscl,iter,ierr
      real ::  buf(3+nscl)
      real :: ufree,tol,vsfc,utau2
      real :: zeta,zeta_mn,zeta_mx,zeta_a
      real :: f_con,f_new,d_theta,u_fac
      real :: phim,phis,psim,psis,t_fac
      real :: dnom,tep,thta
      real :: mod_magnus


c
c ---------- limiting value for wind
c
      ufree = 0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)
c
c ---- save old utau
c
      utausv = utau
      utau2  = utau*utau
c
      iz   = 1
      izp1 = iz + 1
      izm1 = iz - 1
c
      buf(1)  = 0.0
      buf(2)  = 0.0
      buf(3)  = 0.0
      tol     = 0.01
      do iy=iys,iye
      do ix=1,nnx
         buf(1) = buf(1) + u(ix,iy,iz)
         buf(2) = buf(2) + v(ix,iy,iz)
         wind(ix,iy) = sqrt((u(ix,iy,iz)+ugal)**2
     +                    +v(ix,iy,iz)*v(ix,iy,iz))
         buf(3) = buf(3) + wind(ix,iy)
      enddo
      enddo
      do iscl=1,nscl
         buf(3+iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            buf(3+iscl) = buf(3+iscl) + t(ix,iy,iscl,iz)
         enddo
         enddo
      enddo
c
c -------- get x-y slab sums
c
      call mpi_sum_xy(buf,myid,iss,ise,(3+nscl))
      u1xy  = buf(1)*fnxy + ugal
      v1xy  = buf(2)*fnxy
      windm = buf(3)*fnxy
      do iscl=1,nscl
         t1xy(iscl) = buf(3+iscl)*fnxy
      enddo

      vsfc  = sqrt(u1xy*u1xy+v1xy*v1xy)
      windm = amax1(windm,ufree)
      vsfc  = amax1(vsfc,ufree)

c
c ---------- limits for zeta
c
      zeta_mn = zeta_min
      zeta_mx = zeta_max
      if(isfc(1) .eq. 0) then
         f_con = z1*batag*vk*qstar(1)/((windm*vk)**3)
      else
         d_theta = vk74in*(tsfcc(1) - t1xy(1))
         f_con   = z1*batag*vk*d_theta/((windm*vk)**2)
      endif
c
c --------- iteration for zeta
c
      do iter=1,iter_mo
         zeta_a = 0.5*(zeta_mn + zeta_mx)
         if(ismlt .eq. 1) then
             call busngr(zeta_a,phim,phis,psim,psis)
         else
             call fzol(zeta_a,phim,phis,psim,psis)
         endif
         u_fac = (zody - psim)
         if(isfc(1) .eq. 0) then
            f_new =  zeta_a + f_con*u_fac**3
         else
            t_fac = 1.0/(zody - psis)
            f_new =  zeta_a + f_con*u_fac*u_fac*t_fac
         endif
         if(f_new .lt. 0.0) then
            zeta_mn = zeta_a
         else
            zeta_mx = zeta_a
         endif
c
c ----------- iteration details
c
         utau      = windm*vk/(zody-psim)
!         write(nprt,1000) iter, zeta_a, utau, phim, psim
 1000    format(' 1000 iter = ',i5,' zeta = ',e15.6,' u_* = ',e15.6,
     +          ' phim = ',e15.6,' psim = ',e15.6)
      enddo
c
c --------- check if neutral surface layer
c

      !if (ibuoy.eq.0 .or. abs(qstar(1)) .lt. 1.0e-8) then
      if (ibuoy.eq.0 ) then
          amonin    = 1000.
          zeta      = 0.0
          utau      = windm*vk/zody
          thstar(1) = 0.0
          t10xy(1)  = 0.0

          !tsfcc(1)  = t1xy(1)
          dnom = zosdy*vk74in
          thstar(1) = (t1xy(1) - tsfcc(1))/dnom
          t10xy(1)  = thstar(1)*dnom
          qstar(1)  = -utau*thstar(1)
          u10       = utau/vk*(alog(10.0/zo))
      else
         utau = windm*vk/(zody-psim)
         dnom = (zosdy-psis)*vk74in
         if(isfc(1) .eq. 0) then
            thstar(1) = -qstar(1)/utau
            tsfcc(1)  = t1xy(1)-thstar(1)*dnom
            t10xy(1)  = thstar(1)*dnom
         else
            thstar(1) = (t1xy(1) - tsfcc(1))/dnom
            t10xy(1)  = thstar(1)*dnom
            qstar(1)  = -utau*thstar(1)
         endif
         amonin = -utau**3/(batagk*qstar(1))
         zeta   = z1/amonin
         if (ismlt .eq. 1) then
            call busngr(zeta,phim,phis,psim,psis)
         else
            call fzol(zeta,phim,phis,psim,psis)
         endif
         u10 = utau/vk*(alog(10.0/zo) - psim)
      endif
c
c ------- surface details, for debug
c
!      if (myid.eq.0) then
!      write(*,2000) windm, utau, qstar(1), tsfcc(1), amonin, zeta,
!     +              z1, batag, vk, batagk, zo, zos
! 2000 format(' 2000 suft ',/,
!     +       '    windm = ',e15.6,' utau = ',e15.6,' qstar = ',e15.6,/,
!     +       '    tsfcc = ',e15.6,' MO L = ',e15.6,' z1/L = ',e15.6,/,
!     +       '    z1 = ',e15.6,' batag = ',e15.6,' vk = ',e15.6,/,
!     +       '    batagk = ',e15.6,' zo = ',e15.6,' zos = ',e15.6)
!      end if
c
      if (utau.gt.10.0) then
         write(6,9000)
         write(6,9200) utau,windm
         go to 9999
      endif
      if (t10xy(1).gt.0. .and. qstar(1) .gt. 0.) then
         write(6,9000)
         write(6,9300) u1xy,v1xy,t1xy(1),
     +                 tsfcc(1),amonin,utau,it
         go to 9999
      endif

!Now that temperature and momentum have been done, apply BCs on remaining scalars:

      !Treat humidity differently
      iscl=2

      !if (ibuoy.eq.0 .or. abs(qstar(1)) .lt. 1.0e-8) then
      if (ibuoy.eq.0 ) then
          
          thstar(iscl) = 0.0
          t10xy(iscl)  = 0.0
          tsfcc(iscl)  = t1xy(iscl)

          dnom = zosdy*vk74in
          tsfcc(2)=surf_RH/100.0*
     +    Mw/Ru/tsfcc(1)*mod_magnus(tsfcc(1))/surf_rho
            thstar(iscl) = (t1xy(iscl) - tsfcc(iscl))/dnom
            t10xy(iscl)  = thstar(iscl)*dnom
            qstar(iscl)  = -utau*thstar(iscl)
      else
         dnom = (zosdy-psis)*vk74in
         if(isfc(iscl) .eq. 0) then
            thstar(iscl) = -qstar(iscl)/utau
            tsfcc(iscl)  = t1xy(iscl)-thstar(iscl)*dnom
            t10xy(iscl)  = thstar(iscl)*dnom
         else

            !Humidity is a Dirichlet condition at surf_RH
            !Based on temperature in tsfcc(1)
            tsfcc(2)=surf_RH/100.0*
     +      Mw/Ru/tsfcc(1)*mod_magnus(tsfcc(1))/surf_rho

            thstar(iscl) = (t1xy(iscl) - tsfcc(iscl))/dnom
            t10xy(iscl)  = thstar(iscl)*dnom
            qstar(iscl)  = -utau*thstar(iscl)
         endif
         amonin = -utau**3/(batagk*qstar(iscl))
         zeta   = z1/amonin
      endif


      !Other scalars:
      do iscl=3,nscl

      if (ibuoy.eq.0 .or. qstar(1) .eq. 0.) then
          thstar(iscl) = 0.0
          t10xy(iscl)  = 0.0
          tsfcc(iscl)  = t1xy(iscl)
      else
         dnom = (zosdy-psis)*vk74in
         if(isfc(iscl) .eq. 0) then
            thstar(iscl) = -qstar(iscl)/utau
            tsfcc(iscl)  = t1xy(iscl)-thstar(iscl)*dnom
            t10xy(iscl)  = thstar(iscl)*dnom
         else
            thstar(iscl) = (t1xy(iscl) - tsfcc(iscl))/dnom
            t10xy(iscl)  = thstar(iscl)*dnom
            qstar(iscl)  = -utau*thstar(iscl)
         endif
         amonin = -utau**3/(batagk*qstar(iscl))
         zeta   = z1/amonin
      endif

      enddo

c ---------- examples of two other scalars
c
c     c
c     c **** get flux of b scalar, specified surface value
c     c
c           dnom      = (zody-psis)*vk74in
c           thstar(2) = (t1xy(2)-tsfcc(2))/dnom
c           qstar(2)  = -thstar(2)*utau
c           t10xy(2)  = thstar(2)*dnom
c           aut3m(2)  =  qstar(2)
c
c **** get surface value of c scalar, specified surface flux
c
c     dnom      = (zody-psis)*vk74in
c     thstar(2) = -qstar(2)/utau
c     tsfcc(2)  = t1xy(2) - dnom*thstar(2)
c     t10xy(2)  = thstar(2)*dnom
c     aut3m(2)  = qstar(2)
c
      zol = zeta
      hol = zol*zi/z1
c
c ---- note roundoff problem in angles if close to multiples of pi
c
      tep = u1xy/vsfc
      if(tep.gt.1.)  tep = 1.0
      if(tep.lt.-1.) tep = -1.0
      thta      = acos(tep)
      utau2     = utau*utau
      au13m     = -utau2*cos(thta)
      au23m     = -utau2*sin(thta)*sign(1.,v1xy)
      do iscl = 1,nscl
      aut3m(iscl)  =  qstar(iscl)
      end do

c
      return
c
c -------- iteration did not converge
c
 9999 continue
 9000 format(' Trouble in SR. suft')
 9200 format(' Stop because utau = ',e15.6,' windm = ',e15.6)
 9300 format(' ** CHECK SFC U = ',e15.6,' V=',e15.6,' T,TS = ',2e15.6,
     +       ' L =',e15.6,' U_* = ',e15.6,' AT IT = ',i5)
      call mpi_finalize(ierr)
      stop
      end
      subroutine sufto
c
      use pars
      use fields
      use con_data
      use con_stats
      real buf(3+nscl)
c
c ------- version of similarity theory adpated for ocean flows
c      option to use businger or large version of similarity theory
c
      iz    = 1
      izm1  = iz - 1
      izp1  = iz + 1
      z1_a  = abs(z1)
      buf(1)  = 0.0
      buf(2)  = 0.0
      buf(3)  = 0.0
      tol     = 0.01
      do iy=iys,iye
      do ix=1,nnx
         buf(1) = buf(1) + u(ix,iy,iz)
         buf(2) = buf(2) + v(ix,iy,iz)
         wind(ix,iy) = sqrt((u(ix,iy,iz)+ugal)**2
     +                    +v(ix,iy,iz)*v(ix,iy,iz))
         buf(3) = buf(3) + wind(ix,iy)
      enddo
      enddo
      do iscl=1,nscl
         buf(3+iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            buf(3+iscl) = buf(3+iscl) + t(ix,iy,iscl,iz)
         enddo
         enddo
      enddo
c
c -------- get x-y slab sums
c
      call mpi_sum_xy(buf,myid,iss,ise,(3+nscl))
      u1xy  = buf(1)*fnxy + ugal
      v1xy  = buf(2)*fnxy
      windm = buf(3)*fnxy
      do iscl=1,nscl
         t1xy(iscl) = buf(3+iscl)*fnxy
      enddo
      vsfc  = sqrt(u1xy*u1xy+v1xy*v1xy)
      windm = amax1(windm,ufree)
      vsfc  = amax1(vsfc,ufree)
c
      t10xy(1)=-qstar(1)/utau*zody*vk74in
c
c ---- check for temperature boundary condition
c
      if(isfc(iscl) .eq. 0 ) then
         tsfcc(1)=t1xy(1)-t10xy(1)
      endif
c     vsfc=sqrt(u1xy*u1xy+v1xy*v1xy)
c     if(windm.le.0.01)windm=0.01
c     if(vsfc .le.0.01)vsfc =0.01
c
c ----------- input surface wind stress (tau = 0.0184n/m*m)
c             density rho = 1000kg/m^3
c          
c     utau = 4.29e-03
c     utau = 6.10e-03
      utau = 7.00e-03
c
c **** save old utau
      utausv = utau
      utau2  = utau*utau
      if (ibuoy.eq.0 .or. qstar(1) .eq. 0.) then
          amonin    = 1000.
          zeta      = 0.
          thstar(1) = 0.0
          t10xy(1)  = 0.0
      else
          amonin = -utau2*utau/(batagk*qstar(1))
          zeta   = z1_a/amonin
      endif
      if (t10xy(1).lt.0. .and. qstar(1) .lt. 0.) then
         write(6,1234)u1xy,v1xy,t1xy(1),tsfcc(1),amonin,utau,it
 1234    format(' ** check sfc u=',e12.3,' v=',e12.3,' t,ts=',2f10.3,
     +     ' l=',e12.3,' u*=',e12.3,' at it=',i5)
         go to 9999
      endif
c
c -------- for stable,neutral and unstable pbl get drift velocity
c
      if(ismlt .eq. 1) then
          call busngr(zeta,phim,phis,psim,psis)
      else
          call fzol(zeta,phim,phis,psim,psis)
      endif
      udrift = windm + stokes(1) - stokess + utau*(zody-psim)*vkin
      vdrift = 0.0
      dnom      = (zody-psis)*vk74in
      if (isfc(iscl).eq.1) then
         thstar(1) = (t1xy(1) - tsfcc(1))/dnom
         t10xy(1)  = thstar(1)*dnom
         qstar(1)  = - utau*thstar(1)
      else
         thstar(1)  = -qstar(1)/utau
         tsfcc(1)   = t1xy(1)-thstar(1)*dnom
         t10xy(1)   = thstar(1)*dnom
      endif
      zol = zeta
      hol = zol*zi/z1
c
c ---------- examples of two other scalars
c
c     c
c     c **** get flux of b scalar, specified surface value
c     c
c           dnom      = (zody-psis)*vk74in
c           thstar(2) = (t1xy(2)-tsfcc(2))/dnom
c           qstar(2)  = -thstar(2)*utau
c           t10xy(2)  = thstar(2)*dnom
c           aut3m(2)  =  qstar(2)
c     c
c     c **** get surface value of c scalar, specified surface flux
c     c
c           dnom      = (zody-psis)*vk74in
c           thstar(3) = -qstar(3)/utau
c           tsfcc(3)  = t1xy(3) - dnom*thstar(3)
c           t10xy(3)  = thstar(3)*dnom
c           aut3m(3)  = qstar(3)
c
c **** note roundoff problem in angles are close to multiples of pi
c     tep=u1xy/vsfc
c     if(tep.gt.1.)tep=1.
c     if(tep.lt.-1.)tep=-1.
c     thta=acos(tep)
      utau2 = utau*utau
c     au13m = -utau2*cos(thta)
c     au23m = -utau2*sin(thta)*sign(1.,v1xy)
      au13m = utau2
      au23m = 0.0
      aut3m(1)= qstar(1)
c
      return
c
c --------- trouble in sl routine
c
 9999 continue
c
      write(nprt,9000)
 9000 format(' Trouble in SR. sufto')
      call mpi_finalize(ierr)
      stop
      end
      subroutine suft2(u_level1)
c
      use pars
      use fields
      use con_data
      use con_stats
c
      real u_level1(nnx,iys:iye,2+nscl)
c
c     u_level1(.,.,1) = u
c     u_level1(.,.,2) = v
c     u_level1(.,.,3) = theta
c     u_level1(.,.,4) = more scalars
c
      tol = 0.01
      ufree=0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)
      zeta_mn = -6.0
      zeta_mn_i = 1.0/zeta_mn
      iz   = 1
c     izm1 = iz - 1
c     izp1 = iz + 1
c
c      write(nprt,3131) myid, utau, zody, vk74in, batagk,
c    +               u_level1(jxs,1,3), u_level1(jxe,1,3),
c    +               u_level1(jxs,1,1), u_level1(jxe,1,1),
c    +               u_level1(jxs,1,2), u_level1(jxe,1,2)
c3131  format(' in suft2 myid = ',i4,/,
c    +        ' utau = ',e15.6,' zody = ',e15.6,/,
c    +        ' vk74in = ',e15.6,' batagk = ',e15.6,/,
c    +        ' t(jxs) = ',e15.6,' t(jxe) = ',e15.6,/,
c    +        ' u(jxs) = ',e15.6,' u(jxe) = ',e15.6,/,
c    +        ' v(jxs) = ',e15.6,' v(jxe) = ',e15.6)
c
      do iy=iys,iye
      do ix=mxs,mxe
c
c ----------------- first guess for utau
c
      utau = .001
c
      t10xy(1) = -qstar(1)/utau*zosdy*vk74in
      tsfcc(1) = u_level1(ix,iy,3) - t10xy(1)
      vsfc2    = u_level1(ix,iy,1)**2 + u_level1(ix,iy,2)**2
      vsfc     = sqrt(vsfc2)
      windm    = ufree+vsfc
      utausv   = utau
      utau2    = utau*utau
      amonin   = -utau2*utau/(batagk*qstar(1))
      if(amonin.eq.0.) then
            write(6,5050) ix,iy,it,utau,amonin
 5050       format(' 5050, sr. suft2, trouble at ',/,
     +             ' ix = ',i6,'iy = ',i6,' it = ',i6,' utau = ',e15.6,
     +             ' amonin = ',e15.6)
            stop
      endif
c
c ---- for unstable, free convection pbl
c
      iter = 0
 100  continue
c
c ----------------- limit the min (-l/z) change to accmmodate stable flow
c
      zeta_i = amin1(amonin/z1,zeta_mn_i)
      zeta_a = 1.0/zeta_i
c
      if(ismlt .eq. 1) then
          call busngr(zeta_a,phim,phis,psim,psis)
      else
          call fzol(zeta_a,phim,phis,psim,psis)
      endif
      utau     = windm*vk/(zody-psim)
      thstar(1)=-qstar(1)/utau
      amonold  = amonin
      amonin   = utau*utau/(batagk*thstar(1))
      diff     = abs(amonin - amonold)
c      write(nprt,5656)iter,psim,utau,zeta,amonin,dmonin,diff
c 5656 format(' iter=',i4,' phm=',e10.3,' utau=',e10.3,
c     1      ' zeta=',e10.3,' l=',e10.3,' diff = ',e12.4)
      iter = iter+1
      if(iter.gt.10)go to 1000
      if(diff.gt.abs(tol*amonin)) go to 100
 1000 continue
c
 2000 continue
c
      if (utau.gt.10.) then
         write(6,232)utau,windm
  232    format(' stop because utau=',e15.6,' windm=',e15.6)
         stop 9999
      endif
      t10xy(1) = -qstar(1)/utau*vk74in*(zosdy-psis)
      t_grnd(ix,iy,1) = u_level1(ix,iy,3) - t10xy(1)
c
      zol = zeta_a
      hol = zol*zi/z1
      tep = u_level1(ix,iy,1)/windm
      if(tep.gt.1.)  tep = 1.0
      if(tep.lt.-1.) tep = -1.0
      thta  = acos(tep)
      utau2 = utau*utau
c     au13m=-utau2*cos(thta)
c     au23m=-utau2*sin(thta)*sign(1.,u_level1(ix,iy,2))
c     aut3m(1)= qstar(1)
c
      tau13m(ix,iy)   = -utau2*cos(thta)
      tau23m(ix,iy)   = -utau2*sin(thta)*sign(1.,u_level1(ix,iy,2))
      taut3m(ix,iy,1) = qstar(1)
c
c **** get surface value of c scalar, specified surface flux
c
c     dnom      = (zosdy-psis)*vk74in
c     thstar(2) = -qstar(2)/utau
c     tsfcc(2)  = u_level1(ix,iy,4) - dnom*thstar(2)
c     t_grnd(ix,iy,2)  = u_level1(ix,iy,4) - dnom*thstar(2)
c     t10xy(2)  = thstar(2)*dnom
c     taut3m(ix,iy,2)  = qstar(2)
c
c
c ------- end of x-y loops
c
      enddo
      enddo
c
      return
      end
      subroutine init
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
c

      pi   = 4.0*atan(1.0)
      pi2  = 2.0*pi
      bfac = 1.0
      if(ibuoy.eq.0) bfac = 0.
c
c -------------------- case specific data
c
      if(iocean .eq. 1) then
         t00     = 283.
         t00b    = 5000.0
         cp      = 4.20e03
         gcp     = grav/cp
         batag   = bfac*grav/t00b
c        fcor    = 0.0
         fcor    = 1.39e-04
         fcor_h  = 0.0
c        wtsfc(1)=0.00
c        wtsfc(1)=4.96e-07
         wtsfc(1)=1.190476e-06
         qstar(1)=wtsfc(1)
c        dtdzf(1)=0.000
         dtdzf(1)=0.2548
         dtjump  = 0.
         divgls  = 0.
         zo      = 0.0001
         zi      = -5.
         izi     = 55
         xl      = 50.
         yl      = 50.
         zl      = -20.
c
c ---------- if stretched grid specify location of first point
c
         zw1 = -0.5
      else
         gcp     = grav/Cpa
         batag   = bfac*grav/t00
         fcor_h  = 0.0

         wtsfc = qstar
c
c
         !Gradients of temperature and scalars above inversion (upper BC)
         if (icase.eq.5) then  !Stratocumulus

         dtdzf(1)=(311.85-308.2)/(3000.0-2000.0)
         dtdzf(2)=(3.0e-3-4.2e-3)/(3000.0-2000.0)

         elseif (icase.eq.3) then !C-FOG

         dtdzf(1)=(285.5-284.0)/(80.0-30.0)
         dtdzf(2)=(0.006-0.00805)/(80.0-30.0)         

         elseif (icase.eq.6) then !FATIMA

         dtdzf(1)=(289.5-287.51)/(128.0-70.0)
         dtdzf(2)=(0.0094-0.0098)/(128.0-70.0)

         end if

         dtjump  = 0.0
         divgls  = 0.0

      endif
c
      time  = 0.0
c 
c ---------- outermost coarse grid  indicies are bounds of grid
c
      izlow = 1
      izup  = nnz
      dz    = zl/nnz
      dzg   = abs(dz)
      if(l_root) write(6,4040) zl,nnz,dzg
c
c --------------- generate z grids for particular mesh from
c                 iz = 0,1,...,nnz+1; this allows indexing
c                 to array elements z(0), etc.
c
      zwstrt = 0.0

c
c ------------ build z grid for w points
c
      if(iz_space .eq. 0) then
         do iz=0,nnz+1
            z(iz) = dz*float(iz) + zwstrt
         enddo
      elseif (iz_space .eq. 1) then
        call vgrid_channel(zw1,zi,zl,nnz,z(0),l_root,l_debug)
      elseif (iz_space .eq. 2) then
        call vgrid_channel_fstrm(zw1,zi,zl,nnz,z(0),l_root,l_debug)
      elseif (iz_space .eq. 3) then
        call vgrid(zw1,zi,zl,nnz,z(0),l_root,l_debug)
      endif
c
      call get_dz
c
      if(l_root) then
         write(6,8002) zwstrt
         write(6,8003) (iz,z(iz),zz(iz),iz=0,nnz+1)
      endif
c
      nnzm1 = nnz-1
      dx    = xl/nnx
      dy    = yl/nny
      fnxy  = 1./float(nxy)
      dzdz  = dzw(1)*dzw(1)
      z1    = zz(1)
c
      c23  = 2.0/3.0
      dsl  = (dx*1.5*dy*1.5*abs(dzw(1)))**(1./3.)
      dslg = dsl
      cs   = 0.2
c
      vk     = 0.4
      batagk = batag*vk
      vkin   = 1./vk
      ttmean = 0.
      zody   = alog(abs(z1/zo))
      zosdy   = alog(abs(z1/zos))
      write(nprt, 9901) z1,zo,zody,zosdy
 9901 format(' 9901 z1 = ',e15.6,' zo = ',e15.6,/,
     +       ' zody = ',e15.6,' zosdy = ',e15.6)
      zodyin = 1./zody
      wstar  = abs(batag*zi*wtsfc(1))**(1./3.)
      if(ismlt .eq. 1) then
c
c ---- set constants for businger similarity functions
c
         vk74   = vk*0.74
         vk74in = 0.74/vk
         zody74 = zody*0.74
      else 
c
c ---- set constants for large similarity functions
c
        vk74    = vk
        vk74in  = 1.0/vk
        zody74  = zody
      endif
      ugal   = 0.0
c      ugal   = ugcont*0.5
c     ugcont = ugcont - ugal
      cdbtm  = vk*vk/zody/zody
      if(iocean .eq. 1) then
c ----------- set surface friction velocity here and in sr. sufto
c        utau = 4.29e-03
         utau = 7.00e-03
      else
         ufree = 0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)
c
c ---- note : new estimate for utau !!!
c
         utau  = vk*(ufree+ugcont)/zody
c        utau  = vk*(ufree)/zody
      endif
      utau2    = utau*utau
      if(ibuoy .eq. 0 .or. qstar(1) .eq. 0.) then
        amonin = 1000.0
      else
        amonin = -utau2*utau/(batagk*qstar(1))
      endif
      hol   = abs(zi)/amonin
      zol   = abs(z1)/amonin
      uwsfc = -utau*utau
      vwsfc = -utau*utau
c
      if(l_root) then
         write(6,80)
         write(6,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +         ,cdbtm,ugcont
      endif
c
      if(l_debug) then
         write(nprt,80)
         write(nprt,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +         ,cdbtm,ugcont
      endif
c
      return
c ------------------------
   2  format(10x,' WT =',e12.4,',  U* =',e12.4,',  L =',e12.4,/,
     +       10x,' DTDZ FREE =',e12.4,',  ZODY=',e12.4,/,10x,
     +       ' ZO(BTM) =',e12.4,',  CDBTM=',e12.4,
     +       ',  UG = ',e12.4)
  80  format(///,' ***** SCRATCH RUN ***** ',//)
 4040 format(' zl = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 4043 format(' znest = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 8002 format(' zwstrt = ',e12.4)
 8003 format(' iz ',5x,' zw',5x,' zu ',5x,/,(i3,2e12.4))
      end
      subroutine vgrid(z1,zi,zl,nnz,z,l_root,ldebug)
c
      real z(0:nnz+1)
      logical l_root, l_debug
c
c ----------------- build grid up to zi first
c
      z_frst = z1
      z_cntr = zi*0.5
      !n_pbl  = nnz/2
      n_pbl  = (3*nnz)/4
      z_fac1 = z_cntr/z_frst
      z_fac2 = 1.0/float(n_pbl/2)
      z_fac  = 1.1
      knt = 0
      tol = 0.00001
   10 continue
        knt = knt + 1
        z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
        test   = abs(1.0 - z_facn/z_fac)
        if(knt .gt. 50) then
            if(l_root) write(6,9000) z_fac, z_facn, knt
 9000       format(' Cannot find stretching factor',/,
     +             ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 10
      if(l_root) write(6,9100) z_fac, z_cntr, z1, knt
 9100 format(' Stretching factor = ',e15.6,/,
     +       ' Match point       = ',e15.6,/,
     +       ' First z           = ',e15.6,/,
     +       ' Number of iters   = ',i4)
      z(1) = z_frst
      do iz=2,n_pbl/2-1
         z(iz) = z_frst*(z_fac**(float(iz)) - 1.0)/(z_fac - 1.0)
      enddo
      z(n_pbl/2) = z_cntr
      do iz=1,n_pbl/2 - 1
         z(n_pbl-iz) = zi - z(iz)
      enddo
      z(n_pbl) = zi
      z(0)   = 0.0
c
      if(l_root) write(6,5300) n_pbl
 5300 format(' n_pbl = ',i4)
c
c -------------- build grid from zi on up
c
      z_frst = z1
      z_cntr = zl - zi
      n_top  = nnz - n_pbl
      z_fac1 = z_cntr/z_frst
      z_fac2 = 1.0/float(n_top)
      z_fac  = 1.1
      knt = 0
      tol = 0.00001
   20 continue
        knt = knt + 1
        z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
        test   = abs(1.0 - z_facn/z_fac)
        if(knt .gt. 50) then
            if(l_root) write(6,8000) z_fac, z_facn, knt
 8000       format(' Cannot find stretching factor',/,
     +             ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 20
      if(l_root) write(6,8100) z_fac, z_cntr, z1, knt
 8100 format(' Stretching factor = ',e15.6,/,
     +       ' Match point       = ',e15.6,/,
     +       ' First z           = ',e15.6,/,
     +       ' Number of iters   = ',i4)
c
      z(n_pbl+1) = zi + z_frst
      do iz=n_pbl+2,nnz-1
         z(iz) = zi + z_frst*
     +           (z_fac**(float(iz-n_pbl)) - 1.0)/(z_fac - 1.0)
      enddo
      z(nnz) = zl
      z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))
c     if(l_root) write(6,5600) (iz,z(iz),iz=0,nnz+1)
 5600 format(' 5600 in vgrid ',/,
     +       ' iz ',5x,' zw ',/,(i3,e15.6))
c
c     write(1,2000)
c2000 format('#k ',/,
c    +       '#lw 0.5 ',/,
c    +       '#m 1',/,
c    +       '#x 0 100 50',/,
c    +       '#y -50 2100 500')
c     x1 = 30.0
c     x2 = 80.0
c     do iz=0,nnz+1
c        write(1,1000) x1,z(iz)
c1000    format('#k ',/,
c    +          (2e15.6))
c        write(1,1100) x2,z(iz)
c1100    format(2e15.6)
c     enddo
c
      return
      end
      subroutine vgrid_channel(z1,zi,zl,nnz,z,l_root,ldebug)
c
      real z(0:nnz+1)
      integer :: zidx
      logical l_root, l_debug
c
c ----------------- build grid up to zi first
c
      z_frst = z1
      z_cntr = zi*0.5
      n_pbl  = nnz
c     n_pbl  = (5*nnz)/8
      z_fac1 = z_cntr/z_frst
      z_fac2 = 1.0/float(n_pbl/2)
      z_fac  = 1.1
      knt = 0
      tol = 0.00001
   10 continue
        knt = knt + 1
        z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
        test   = abs(1.0 - z_facn/z_fac)
        if(knt .gt. 50) then
            if(l_root) write(6,9000) z_fac, z_facn, knt
 9000       format(' Cannot find stretching factor',/,
     +             ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 10
      if(l_root) write(6,9100) z_fac, z_cntr, z1, knt
 9100 format(' Stretching factor = ',e15.6,/,
     +       ' Match point       = ',e15.6,/,
     +       ' First z           = ',e15.6,/,
     +       ' Number of iters   = ',i4)
      z(1) = z_frst
      do iz=2,n_pbl/2-1
         z(iz) = z_frst*(z_fac**(float(iz)) - 1.0)/(z_fac - 1.0)
      enddo
      z(n_pbl/2) = z_cntr
      do iz=1,n_pbl/2 - 1
         z(n_pbl-iz) = zi - z(iz)
      enddo
      z(n_pbl) = zi
      z(0)   = 0.0
c
      if(l_root) write(6,5300) n_pbl
 5300 format(' n_pbl = ',i4)
c
c -------------- build grid from zi on up
!     For the channel, zi represents the channel centerline
!     Want the mesh to be a mirror image across this:
c
!      zidx = 1
!      do iz=n_pbl+1,nnz
!         z(iz) = zi + (zi - z(n_pbl-zidx))
!         zidx = zidx + 1
!      enddo
      z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))
c
      return
      end subroutine vgrid_channel
      subroutine vgrid_channel_fstrm(z1,zi,zl,nnz,z,l_root,ldebug)
c
      real z(0:nnz+1)
      real s(0:2*nnz+1)
      integer :: zidx
      logical l_root, l_debug
c      
                nnz = 2.0*nnz
                zi = 2.0*zi
                zl = 2.0*zl
c
c ----------------- build grid up to zi first
c
      z_frst = z1
      z_cntr = zi*0.5
      n_pbl  = nnz
c     n_pbl  = (5*nnz)/8
      z_fac1 = z_cntr/z_frst
      z_fac2 = 1.0/float(n_pbl/2)
      z_fac  = 1.1
      knt = 0
      tol = 0.00001
   10 continue
        knt = knt + 1
        z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
        test   = abs(1.0 - z_facn/z_fac)
        if(knt .gt. 50) then
            if(l_root) write(6,9000) z_fac, z_facn, knt
 9000       format(' Cannot find stretching factor',/,
     +             ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 10
      if(l_root) write(6,9100) z_fac, z_cntr, z1, knt
 9100 format(' Stretching factor = ',e15.6,/,
     +       ' Match point       = ',e15.6,/,
     +       ' First z           = ',e15.6,/,
     +       ' Number of iters   = ',i4)
      s(1) = z_frst
      do iz=2,n_pbl/2-1
         s(iz) = z_frst*(z_fac**(float(iz)) - 1.0)/(z_fac - 1.0)
      enddo
      s(n_pbl/2) = z_cntr
      do iz=1,n_pbl/2 - 1
         s(n_pbl-iz) = zi - s(iz)
      enddo
      s(n_pbl) = zi
      s(0)   = 0.0
c
      if(l_root) write(6,5300) n_pbl
 5300 format(' n_pbl = ',i4)
c
      s(nnz+1) = s(nnz) + (s(nnz) - s(nnz-1))

                nnz=nnz/2.0
                zi = zi/2.0
                zl = zl/2.0

                do iz=0,nnz
                        z(iz)=s(iz)
                enddo
                 z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))

      return
      end subroutine vgrid_channel_fstrm

      subroutine vgrid_uniform(z1,zi,zl,nnz,z,l_root,ldebug)
c
      real z(0:nnz+1),zdiff
      real s(0:2*nnz+1)
      integer :: zidx
      logical l_root, l_debug

      z_frst = z1
      n_pbl  = nnz
      z_fac  = 1
      zdiff = (zl-z_frst)/(nnz-1)
      s(1) = z_frst
      do iz=2,nnz
         s(iz) = z_frst+(zdiff*((float(iz)-1)))
      enddo
      s(0)   = 0.0
c
      s(nnz+1) = s(nnz) + (s(nnz) - s(nnz-1))

      do iz=0,nnz
         z(iz)=s(iz)
      enddo
         z(nnz+1) = z(nnz) + zdiff
      return
      end subroutine vgrid_uniform

      subroutine get_dz
c
c --------------- compute spacing for given vertical
c                 point distribution
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
      include 'mpif.h'
c
      do iz=1,nnz+1
         dzw(iz) = z(iz) - z(iz-1)
      enddo
      dzw(0)     = dzw(1)
      dzw(nnz+2) = dzw(nnz+1)
      do iz=0,nnz+2
         dzw_i(iz) = 1.0/dzw(iz)
      enddo
c
c ------------ build z grid for u points
c
      dzovr2 = dz*0.5
      do iz=1,nnz+1
         zz(iz) = 0.5*(z(iz) + z(iz-1))
      enddo
      zz(0) = - zz(1)
      do iz=1,nnz+1
         dzu(iz) = zz(iz) - zz(iz-1)
      enddo
      dzu(0)     = dzu(1)
      dzu(nnz+2) = dzu(nnz+1)
      do iz=0,nnz+2
         dzu_i(iz) = 1.0/dzu(iz)
      enddo
c
      return
      end
      subroutine random
c
c ----------- geostrophic winds designed for comparison case
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
      implicit none
 
      real ::  psi(nnx,iys:iye), psix(nnx,iys:iye),
     +     psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +     vyy(nnx,iys:iye,izs:izs),Ttmp

      real :: ampv,ampt,sum_psi,vmaxx,vmag,facv,ran1

      real :: slope_q,slope_t,slope_u,slope_v
      real :: slope_t1,slope_t2,slope_u1,slope_v1
      real :: z_switch,t_switch,q_switch,u_switch,v_switch
      real :: z_switch_q,z_switch_t,z_switch_v,z_switch_u
      real :: t_surf,q_surf,u_surf,v_surf
      real :: t_temp,q_temp,u_temp,v_temp
      real :: q_top,t_top,u_top,v_top,slope_q1,slope_q2

      integer :: ix,iy,iz
      integer :: idum


c ------------ note set nmatch in sr. iso so that
c              it is compatible with conditions here
c
      do iz=1,nnz
         ug(iz)   = ugcont
         vg(iz)   = vgcont
         divz(iz) = 0.0
      enddo

      if (icase.eq.5) then

      do iz=izs,ize

          ! MAGPIE
c         profile for potential temperature
            if (zz(iz) .le. 600.0) then

                  t_temp = 302.75   !301.75   !300.0   !301.75(special)   !302.5(ori)

            elseif (zz(iz) .gt. 600.0 .and. zz(iz) .le. 1900.0) then

                  t_surf   = 302.75   !301.75   !300.0   !301.75(special)   !302.5(ori)
                  t_switch = 308.0   !316.0   !315.5(special) !8   !313.0 !5   !308.0(ori)
                  slope_t  = (t_switch - t_surf) / (1900.0 - 600.0)
                  t_temp   = t_surf + (slope_t * (zz(iz) - 600.0))

            elseif (zz(iz) .gt. 1900.0 .and. zz(iz) .le. 2400.0) then

                  t_surf   = 308.0   !316.0   !315.5(special) !8   !313.0 !5   !308.0(ori)
                  t_switch = 315.0   !324.0   !323.0(special) !9   !321.0 !6   !315.0(ori)
                  slope_t  = (t_switch - t_surf) / (2400.0 - 1900.0)
                  t_temp   = t_surf + (slope_t * (zz(iz) - 1900.0))

            elseif (zz(iz) .gt. 2400.0 .and. zz(iz) .le. 3000.0) then

                  t_surf   = 315.0   !324.0   !323.0(special) !9   !321.0 !6   !315.0(ori)
                  t_switch = 317.25   !327.25   !326.25(special) !10   !324.25 !7  !317.25(ori)
                  slope_t  = (t_switch - t_surf) / (3000.0 - 2400.0)
                  t_temp   = t_surf + (slope_t * (zz(iz) - 2400.0))

            elseif (zz(iz) .ge. 3000.0) then

                  t_temp = 317.25   !327.25   !326.25(special) !10   !324.25 !7   !317.25(ori)

            endif

c         profile for qv
            if (zz(iz) .le. 600.0) then

                  q_surf   = 19.25e-3   !19.5e-3
                  q_switch = 18.0e-3   !17.05e-3
                  slope_q  = (q_switch - q_surf) / 600.0
                  q_temp   = q_surf + (slope_q * zz(iz))

            elseif (zz(iz) .gt. 600.0 .and. zz(iz) .le. 1800.0) then

                  q_surf   = 18.0e-3   !17.05e-3
                  q_switch = 11.5e-3   !8.5e-3
                  slope_q  = (q_switch - q_surf) / (1800.0 - 600.0)
                  q_temp   = q_surf + (slope_q * (zz(iz) - 600.0))

!            elseif (zz(iz) .gt. 1000.0 .and. zz(iz) .le. 2000.0) then

!                  q_surf   = 14.0e-3
!                  q_switch = 9.0e-3
!                  slope_q  = (q_switch - q_surf) / (2000.0 - 1000.0)
!                  q_temp   = q_surf + (slope_q * (zz(iz) - 1000.0))

            elseif (zz(iz) .gt. 1800.0 .and. zz(iz) .le. 2400.0) then

                  q_surf   = 11.5e-3   !8.5e-3
                  q_switch = 5.15e-3
                  slope_q  = (q_switch - q_surf) / (2400.0 - 1800.0)
                  q_temp   = q_surf + (slope_q * (zz(iz) - 1800.0))

            elseif (zz(iz) .gt. 2400.0 .and. zz(iz) .le. 3000.0) then

                  q_surf   = 5.15e-3
                  q_switch = 4.25e-3   !3.0e-3
                  slope_q  = (q_switch - q_surf) / (3000.0 - 2400.0)
                  q_temp   = q_surf + (slope_q * (zz(iz) - 2400.0))

            elseif (zz(iz) .ge. 3000.0) then

                  q_temp = 4.25e-3   !3.0e-3

            endif

c         profile for u
            if (zz(iz) .le. 3000.0) then

                  u_surf   = -5.0   !-4.0
                  u_switch = -15.0   !-7.25
                  slope_u  = (u_switch - u_surf) / 3000.0
                  u_temp   = u_surf + (slope_u * zz(iz))

!            elseif (zz(iz) .gt. 100.0 .and. zz(iz) .le. 700.0) then

!                  u_surf = -7.25
!                  u_switch = -8.0
!                  slope_u = (u_switch - u_surf) / (700.0 - 100.0)
!                  u_temp = u_surf + (slope_u * (zz(iz) - 100.0))

!            elseif (zz(iz) .gt. 700.0 .and. zz(iz) .le. 3000.0) then

!                  u_surf = -8.0
!                  u_switch = -15.5
!                  slope_u = (u_switch - u_surf) / (3000.0 - 700.0)
!                  u_temp = u_surf + (slope_u * (zz(iz) - 700.0))

            elseif (zz(iz) .ge. 3000.0) then

                  u_temp = -15.0   !-15.5

            endif

c         profile for v
            if (zz(iz) .le. 3000.0) then

                  v_surf   = -3.5   !-3.25
                  v_switch = 0.0   !-6.5
                  slope_v  = (v_switch - v_surf) / 3000.0
                  v_temp   = v_surf + (slope_v * zz(iz))

!            elseif (zz(iz) .gt. 250.0 .and. zz(iz) .le. 3000.0) then

!                  u_surf = -6.5
!                  u_switch = 0.0
!                  slope_u = (u_switch - u_surf) / (3000.0 - 250.0)
!                  u_temp = u_surf + (slope_u * (zz(iz) - 250.0))

            elseif (zz(iz) .ge. 3000.0) then

                  u_temp = 0.0

            endif



      !!! BOMEX initial profile
!      do iz=1,nnz
!            ug(iz)   = -10.0+1.8e-3*zz(iz)
!            vg(iz)   = 0.0
!            divz(iz) = 0.0
!      enddo       
!
!      do iz=izs,ize
!c         profile for pot temp
!            if (zz(iz).le.520.0) then
!                  t_temp = 298.7
!            elseif (zz(iz).gt.520.0 .and. zz(iz).le.1480.0) then
!                  t_surf = 298.7
!                  t_switch = 302.4
!                  slope_t = (t_switch-t_surf)/(1480.0-520.0)
!                  t_temp = t_surf + slope_t*(zz(iz)-520.0)
!            elseif (zz(iz).gt.1480.0 .and. zz(iz).le.2000.0) then
!                  t_surf = 302.4
!                  t_switch = 308.2
!                  slope_t = (t_switch-t_surf)/(2000.0-1480.0)
!                  t_temp = t_surf + slope_t*(zz(iz)-1480.0)
!            elseif (zz(iz).gt.2000.0 .and. zz(iz).le.3000.0) then
!                  t_surf = 308.2
!                  t_switch = 311.85
!                  slope_t = (t_switch-t_surf)/(3000.0-2000.0)
!                  t_temp = t_surf + slope_t*(zz(iz)-2000.0)
!            elseif (zz(iz).ge.3000.0) then
!                  t_temp = 311.85
!            endif
!
!c         profile for qv
!             if (zz(iz).le.520.0) then
!                  q_surf = 17.0e-3
!                  q_switch = 16.3e-3
!                  slope_q = (q_switch-q_surf)/520.0
!                  q_temp = q_surf + slope_q*zz(iz)
!            elseif (zz(iz).gt.520.0 .and. zz(iz).le.1480.0) then
!                  q_surf = 16.3e-3
!                  q_switch = 10.7e-3
!                  slope_q = (q_switch-q_surf)/(1480.0-520.0)
!                  q_temp = q_surf + slope_q*(zz(iz)-520.0)
!            elseif (zz(iz).gt.1480.0 .and. zz(iz).le.2000.0) then
!                  q_surf = 10.7e-3
!                  q_switch = 4.2e-3
!                  slope_q = (q_switch-q_surf)/(2000.0-1480.0)
!                  q_temp = q_surf + slope_q*(zz(iz)-1480.0)
!            elseif (zz(iz).gt.2000.0 .and. zz(iz).le.3000.0) then
!                  q_surf = 4.2e-3
!                  q_switch = 3.0e-3
!                  slope_q = (q_switch-q_surf)/(3000.0-2000.0)
!                  q_temp = q_surf + slope_q*(zz(iz)-2000.0)
!            elseif (zz(iz).ge.3000.0) then
!                  q_temp = 3.0e-3
!            endif
!
!c         profile for u
!             if (zz(iz).le.700.0) then
!                  u_temp = -8.75
!            elseif (zz(iz).gt.700.0 .and. zz(iz).le.3000.0) then
!                  u_surf = -8.75
!                  u_switch = -4.61
!                  slope_u = (u_switch-u_surf)/(3000.0-700.0)
!                  u_temp = u_surf + slope_u*(zz(iz)-700.0)
!            elseif (zz(iz).ge.3000.0) then
!                  u_temp = -4.61
!            endif

            do iy=iys,iye
            do ix=1,nnx
                  u(ix,iy,iz) = u_temp
                  v(ix,iy,iz) = 0.0
                  w(ix,iy,iz) = 0.0
                  e(ix,iy,iz) = 1.0-zz(iz)/3000.0

c                 pot temp
                  t(ix,iy,1,iz) = t_temp

c                 qv
                  t(ix,iy,2,iz) = q_temp

                  w(ix,iy,iz)   = 0.
                  r1(ix,iy,iz)  = 0.
                  r2(ix,iy,iz)  = 0.
                  r3(ix,iy,iz)  = 0.
                  r4(ix,iy,1,iz)= 0.
                  r4(ix,iy,2,iz)= 0.
                  r5(ix,iy,iz)  = 0.
             enddo
             enddo
       enddo

      elseif (icase.eq.6) then 
      !FATIMA fog setup

            !qv setup
            z_switch_q = 70.0
            q_surf = 0.0099
            q_switch = 0.0098
            q_top = 0.0094
            slope_q1 = (q_switch-q_surf)/z_switch_q
            slope_q2 = (q_top-q_switch)/(128.0-z_switch_q)

            !theta setup
            z_switch_t = 70.0
            t_surf = 287.5
            t_switch = 287.51
            t_top = 289.5
            slope_t1 = (t_switch-t_surf)/z_switch_t
            slope_t2 = (t_top-t_switch)/(128.0-z_switch_t)


            !u setup
            z_switch_u = 128.0
            u_surf = 3.8
            u_switch = 4.1
            u_top = 4.1
            slope_u1 = (u_top-u_surf)/128.0

            !v setup
            z_switch_v = 128.0
            v_surf = 1.4
            v_switch = 0.4
            v_top = 0.4
            slope_v1 = (v_top-v_surf)/128.0

            do iz=izs,ize
            do iy=iys,iye
            do ix=1,nnx

                ! u
                !u(ix,iy,iz) = u_surf + slope_u1*zz(iz)
                u(ix,iy,iz) = 4.0
                 ! v
                !v(ix,iy,iz) = v_surf + slope_v1*zz(iz)
                v(ix,iy,iz) = 0.0

                ! q
                if (zz(iz) .le. z_switch_q) then
                    t(ix,iy,2,iz) = q_surf + slope_q1*zz(iz)
                elseif (zz(iz) .gt. z_switch_q) then
                    t(ix,iy,2,iz) = q_switch +
     +                          slope_q2*(zz(iz)-z_switch_q)
                end if

                ! theta
                if (zz(iz) .le. z_switch_t) then
                    t(ix,iy,1,iz) = t_surf + slope_t1*zz(iz)
                elseif (zz(iz) .ge. z_switch_t) then
                    t(ix,iy,1,iz) = t_switch +
     +                          slope_t2*(zz(iz)-z_switch_t)
                end if

                w(ix,iy,iz) = 0.0
                e(ix,iy,iz) = 0.0
            enddo
            enddo
            enddo
     
      else  !switch the case


      if (icase.eq.3) then
      !!!! Use for the C-FOG fog LCM case:
      !Parameters from Charlotte's C-FOG profile
      z_switch = 30.0
      t_switch = 284.05
      q_switch = 0.00811
      t_surf = 284.0
      q_surf = 0.00823
      end if

   
      if (icase.eq.4) then
      z_switch = 600.0
      t_switch = 300.81
      t_surf = 300.81
      q_switch = 0.0215
      q_surf = 0.0215
      end if

      if (icase.eq.0) then
      t_surf = 285.0
      t_switch = 285.0
      q_surf = 0.0105
      q_switch = 0.0105
      slope_t = 0.0
      slope_q = 0.0
      z_switch = 1000.0
      end if

      !Slopes below the inversion
      slope_t = (t_switch-t_surf)/z_switch
      slope_q = (q_switch-q_surf)/z_switch

c
      do iz=izs,ize


         !Try and alleviate the abrupt IC on the velocity
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = ugcont
            v(ix,iy,iz) = vgcont
            w(ix,iy,iz) = 0.0
            e(ix,iy,iz) = 0.0
         enddo
         enddo

         if(zz(iz) .le. z_switch) then
           do iy=iys,iye
           do ix=1,nnx
              t(ix,iy,1,iz) = t_surf + slope_t*zz(iz)
              t(ix,iy,2,iz) = q_surf + slope_q*zz(iz)
           enddo
           enddo
         elseif(zz(iz) .ge. z_switch) then
           do iy=iys,iye
           do ix=1,nnx
              t(ix,iy,1,iz) = t_switch + (zz(iz) - z_switch)*dtdzf(1)
              t(ix,iy,2,iz) = q_switch + (zz(iz) - z_switch)*dtdzf(2)
           enddo
           enddo
         endif

         do iy=iys,iye
         do ix=1,nnx
            w(ix,iy,iz)   = 0.
            r1(ix,iy,iz)  = 0.
            r2(ix,iy,iz)  = 0.
            r3(ix,iy,iz)  = 0.
            r4(ix,iy,1,iz)= 0.
            r4(ix,iy,2,iz)= 0.
            r5(ix,iy,iz)  = 0.
         enddo 
         enddo 
      enddo

      end if
c
c ------------- set initial random field to be
c               divergence free
c
      idum = -1 - myid
      do iz=izs,ize
c
c ----------- ampv and ampt are max amplitudes of random 
c             velocity and temperature fields
c             make sure ampv is set if free convection so
c             that we have motions at first time step
c
         ampv = 0.0
         ampv = 0.001
         ampt = 0.10
c  
c ------- simple random field scaled between -0.5 and 0.5
c
         sum_psi = 0.0
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = ran1(idum)
            sum_psi = sum_psi + psi(ix,iy)
         enddo
         enddo
         sum_psi = sum_psi*fnxy
         call mpi_sum_xy(sum_psi,myid,iss,ise,1)
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = psi(ix,iy) - sum_psi
            psix(ix,iy)     = psi(ix,iy)
            psiy(ix,iy,izs) = psi(ix,iy)
         enddo
         enddo
         call xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         vmaxx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
            if(vmag .gt. vmaxx) vmaxx = vmag
         enddo
         enddo
         facv = ampv/vmaxx
c
         if (z(iz) .le. 50.0) then
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz)   = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
            v(ix,iy,iz)   = v(ix,iy,iz) + psix(ix,iy)*facv
            t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
         enddo
         enddo
         endif
c
         if(z(iz) .le. 250.0) then
         do iy=iys,iye
         do ix=1,nnx
            e(ix,iy,iz) = 0.4*(1.0 - z(iz)/250.0)**3
         enddo
         enddo
         endif
c
c ---------- check divergence of initial field
c
         do iy=iys,iye
         do ix=1,nnx
            uxx(ix,iy) = u(ix,iy,iz)
            vyy(ix,iy,izs) = v(ix,iy,iz)
         enddo
         enddo
         call xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         do iy=iys,iye
         do ix=1,nnx
            divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
         enddo
         enddo
         divz(iz) = divz(iz)*fnxy
c
c -------- end z loop
c
      enddo
c
      call mpi_sum_z(divz(1),i_root,myid,nnz,1)
c
      write(nprt,6000)
 6000 format(' check of divergence for initial state',/,
     +         ' iz ',5x,' divergence')
      write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100 format(i5,e15.6)
c
c ------------ fix for baroclinic and subsidence effects !!
c
c     do iz=izs,ize
c        ug(iz)=ugcont
c        vg(iz)=vgcont
c        if (.not.(ibrcl.eq.1)) go to 19988
c        if (.not.(iz.le.izi)) go to 19987
c        ug(iz)=0.
c        vg(iz)=0.
c 19987    continue
c 19988    continue
c        zz2=zz(iz)
c        wls(iz)=-divgls*zz2
c        if (.not.(iz.eq.1)) go to 19986
c        do ix=1,nnx
c        uls(ix)=divgls*(dx*float(ix-1)-xl*.5)
c        enddo
c     enddo
c     write(nprt,9)(uls(ix),ix=1,nnx)
c  9  format(1x,8e12.3)
c 19986 continue
c
      return
      end
      subroutine random_f
c
c ---------- example of using given (sparse) initial 
c            sounding profiles (FIX for ncpu_s).
c            
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +     psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +     vyy(nnx,iys:iye,izs:izs)
c
      parameter (nt=12, nz=11)
      real zg(nz), u_i(nz,nt), v_i(nz,nt), theta_i(nz,nt)
      real ui_temp(nz), vi_temp(nz), ti_temp(nz)
      real time_g(nt)
c
      data time_g /
     +  0.0000E+00,  0.3600E+04,  0.7200E+04,  0.1080E+05,  0.1440E+05,
     +  0.1800E+05,  0.2160E+05,  0.2520E+05,  0.2880E+05,  0.3240E+05,
     +  0.3600E+05,  0.3960E+05
     +/
      data zg /
     +  0.1000E+02,  0.3000E+02,  0.5500E+02,  0.9000E+02,  0.1400E+03,
     +  0.2150E+03,  0.3300E+03,  0.5000E+03,  0.7500E+03,  0.1100E+04,
     +  0.1600E+04
     +/
      data u_i /
     + -0.1510E+01, -0.1560E+01, -0.1580E+01, -0.1580E+01, -0.1560E+01,
     + -0.1530E+01, -0.1510E+01, -0.9000E+00, -0.1390E+01, -0.1220E+01,
     + -0.5100E+00,
     + -0.1090E+01, -0.1110E+01, -0.1120E+01, -0.1120E+01, -0.1030E+01,
     + -0.9900E+00, -0.9500E+00, -0.6200E+00, -0.1230E+01, -0.9400E+00,
     +  0.2800E+00,
     + -0.9100E+00, -0.9200E+00, -0.9100E+00, -0.9000E+00, -0.8800E+00,
     + -0.8400E+00, -0.8000E+00, -0.6500E+00, -0.1510E+01, -0.1070E+01,
     +  0.2400E+00,
     + -0.8900E+00, -0.8900E+00, -0.8900E+00, -0.8800E+00, -0.8700E+00,
     + -0.8500E+00, -0.8100E+00, -0.7000E+00, -0.1830E+01, -0.8400E+00,
     +  0.3500E+00,
     + -0.1250E+01, -0.1260E+01, -0.1260E+01, -0.1250E+01, -0.1240E+01,
     + -0.1220E+01, -0.1160E+01, -0.8800E+00, -0.1980E+01, -0.1900E+00,
     +  0.7500E+00,
     + -0.1800E+01, -0.1810E+01, -0.1820E+01, -0.1820E+01, -0.1800E+01,
     + -0.1780E+01, -0.1710E+01, -0.1150E+01, -0.1960E+01,  0.3900E+00,
     +  0.9200E+00,
     + -0.2110E+01, -0.2130E+01, -0.2140E+01, -0.2140E+01, -0.2130E+01,
     + -0.2110E+01, -0.2050E+01, -0.9300E+00, -0.1400E+01,  0.8800E+00,
     +  0.9600E+00,
     + -0.2250E+01, -0.2280E+01, -0.2290E+01, -0.2300E+01, -0.2290E+01,
     + -0.2260E+01, -0.2070E+01, -0.4000E-01, -0.1600E+00,  0.1440E+01,
     +  0.1190E+01,
     + -0.2160E+01, -0.2200E+01, -0.2220E+01, -0.2220E+01, -0.2220E+01,
     + -0.2190E+01, -0.1610E+01,  0.1470E+01,  0.1420E+01,  0.2050E+01,
     +  0.1610E+01,
     + -0.2230E+01, -0.2270E+01, -0.2290E+01, -0.2300E+01, -0.2300E+01,
     + -0.2260E+01, -0.1350E+01,  0.2480E+01,  0.2380E+01,  0.2320E+01,
     +  0.1740E+01,
     + -0.1890E+01, -0.1930E+01, -0.1950E+01, -0.1950E+01, -0.1940E+01,
     + -0.1890E+01, -0.1120E+01,  0.3010E+01,  0.3030E+01,  0.2800E+01,
     +  0.2000E+01,
     + -0.1210E+01, -0.1230E+01, -0.1240E+01, -0.1230E+01, -0.1210E+01,
     + -0.1140E+01, -0.4600E+00,  0.3320E+01,  0.3510E+01,  0.3420E+01,
     +  0.2340E+01
     +/
      data v_i /
     +  0.4800E+00,  0.5100E+00,  0.5300E+00,  0.5700E+00,  0.6900E+00,
     +  0.7300E+00,  0.7600E+00,  0.1410E+01, -0.4200E+00, -0.3060E+01,
     + -0.3500E+01,
     +  0.7800E+00,  0.8100E+00,  0.8400E+00,  0.8900E+00,  0.1060E+01,
     +  0.1110E+01,  0.1130E+01,  0.1190E+01, -0.1040E+01, -0.2900E+01,
     + -0.3440E+01,
     +  0.3000E+00,  0.3200E+00,  0.3400E+00,  0.3800E+00,  0.4800E+00,
     +  0.5300E+00,  0.5800E+00,  0.5300E+00, -0.1330E+01, -0.2040E+01,
     + -0.2830E+01,
     + -0.2700E+00, -0.2600E+00, -0.2400E+00, -0.2200E+00, -0.1800E+00,
     + -0.1300E+00, -0.5000E-01,  0.1000E+00, -0.1170E+01, -0.1100E+01,
     + -0.2370E+01,
     + -0.5500E+00, -0.5400E+00, -0.5300E+00, -0.5100E+00, -0.4800E+00,
     + -0.4100E+00, -0.2600E+00,  0.1700E+00, -0.4200E+00, -0.2200E+00,
     + -0.2080E+01,
     + -0.2700E+00, -0.2600E+00, -0.2500E+00, -0.2400E+00, -0.2100E+00,
     + -0.1600E+00, -0.1000E-01,  0.8500E+00,  0.9700E+00,  0.3500E+00,
     + -0.2250E+01,
     +  0.5300E+00,  0.5400E+00,  0.5600E+00,  0.5700E+00,  0.6000E+00,
     +  0.6500E+00,  0.7600E+00,  0.1960E+01,  0.2280E+01,  0.3600E+00,
     + -0.2590E+01,
     +  0.1590E+01,  0.1630E+01,  0.1650E+01,  0.1680E+01,  0.1720E+01,
     +  0.1780E+01,  0.2010E+01,  0.3260E+01,  0.3110E+01,  0.1600E+00,
     + -0.2580E+01,
     +  0.2560E+01,  0.2620E+01,  0.2660E+01,  0.2690E+01,  0.2740E+01,
     +  0.2830E+01,  0.3400E+01,  0.4030E+01,  0.3030E+01, -0.7000E-01,
     + -0.2320E+01,
     +  0.3500E+01,  0.3600E+01,  0.3650E+01,  0.3700E+01,  0.3750E+01,
     +  0.3860E+01,  0.4580E+01,  0.4100E+01,  0.2450E+01,  0.6000E-01,
     + -0.1770E+01,
     +  0.4500E+01,  0.4640E+01,  0.4700E+01,  0.4760E+01,  0.4830E+01,
     +  0.4930E+01,  0.5420E+01,  0.3960E+01,  0.2000E+01,  0.5000E+00,
     + -0.1150E+01,
     +  0.5290E+01,  0.5470E+01,  0.5550E+01,  0.5620E+01,  0.5690E+01,
     +  0.5790E+01,  0.6070E+01,  0.4000E+01,  0.1910E+01,  0.9700E+00,
     + -0.5600E+00
     +/
      data theta_i /
     +  0.2936E+03,  0.2936E+03,  0.2937E+03,  0.2937E+03,  0.2938E+03,
     +  0.2942E+03,  0.2948E+03,  0.2980E+03,  0.3027E+03,  0.3092E+03,
     +  0.3186E+03,
     +  0.2937E+03,  0.2937E+03,  0.2937E+03,  0.2938E+03,  0.2939E+03,
     +  0.2942E+03,  0.2946E+03,  0.2978E+03,  0.3024E+03,  0.3090E+03,
     +  0.3184E+03,
     +  0.2939E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,
     +  0.2941E+03,  0.2944E+03,  0.2976E+03,  0.3023E+03,  0.3089E+03,
     +  0.3182E+03,
     +  0.2940E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,
     +  0.2941E+03,  0.2943E+03,  0.2975E+03,  0.3022E+03,  0.3087E+03,
     +  0.3181E+03,
     +  0.2940E+03,  0.2940E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,
     +  0.2940E+03,  0.2942E+03,  0.2974E+03,  0.3021E+03,  0.3086E+03,
     +  0.3180E+03,
     +  0.2941E+03,  0.2940E+03,  0.2940E+03,  0.2940E+03,  0.2941E+03,
     +  0.2941E+03,  0.2941E+03,  0.2973E+03,  0.3019E+03,  0.3085E+03,
     +  0.3179E+03,
     +  0.2942E+03,  0.2942E+03,  0.2942E+03,  0.2942E+03,  0.2941E+03,
     +  0.2941E+03,  0.2941E+03,  0.2973E+03,  0.3020E+03,  0.3086E+03,
     +  0.3179E+03,
     +  0.2943E+03,  0.2943E+03,  0.2943E+03,  0.2943E+03,  0.2943E+03,
     +  0.2943E+03,  0.2943E+03,  0.2975E+03,  0.3022E+03,  0.3087E+03,
     +  0.3181E+03,
     +  0.2945E+03,  0.2945E+03,  0.2945E+03,  0.2945E+03,  0.2945E+03,
     +  0.2944E+03,  0.2946E+03,  0.2978E+03,  0.3025E+03,  0.3090E+03,
     +  0.3184E+03,
     +  0.2947E+03,  0.2947E+03,  0.2947E+03,  0.2947E+03,  0.2946E+03,
     +  0.2946E+03,  0.2949E+03,  0.2980E+03,  0.3027E+03,  0.3093E+03,
     +  0.3187E+03,
     +  0.2949E+03,  0.2949E+03,  0.2949E+03,  0.2948E+03,  0.2948E+03,
     +  0.2948E+03,  0.2950E+03,  0.2982E+03,  0.3028E+03,  0.3094E+03,
     +  0.3188E+03,
     +  0.2950E+03,  0.2950E+03,  0.2950E+03,  0.2950E+03,  0.2950E+03,
     +  0.2950E+03,  0.2950E+03,  0.2982E+03,  0.3029E+03,  0.3095E+03,
     +  0.3188E+03
     +/
c
      save time_g, zg, u_i, v_i, theta_i
c
c --------- find time location of initial profiles 
c
      call lterp(nt,time_g,t_factor,jt,jtp1,t_weit)
c
      do iz=1,nz
         ui_temp(iz) = u_i(iz,jt)*(1.0 - t_weit) +
     +                 u_i(iz,jtp1)*t_weit
         vi_temp(iz) = v_i(iz,jt)*(1.0 - t_weit) +
     +                 v_i(iz,jtp1)*t_weit
         ti_temp(iz) = theta_i(iz,jt)*(1.0 - t_weit) +
     +                 theta_i(iz,jtp1)*t_weit
      enddo
c
c ----------- interpolate vertically
c
      do iz=izs,ize
         call lterp(nz,zg,zz(iz),kk,kkp1,weit)
         u_temp = ui_temp(kk)*(1.0 - weit) +
     +            ui_temp(kkp1)*weit
         v_temp = vi_temp(kk)*(1.0 - weit) +
     +            vi_temp(kkp1)*weit
         theta_temp = ti_temp(kk)*(1.0 - weit) +
     +            ti_temp(kkp1)*weit
c
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz)   = u_temp
            v(ix,iy,iz)   = v_temp
            t(ix,iy,1,iz) = theta_temp
            w(ix,iy,iz)   = 0.
            r1(ix,iy,iz)  = 0.
            r2(ix,iy,iz)  = 0.
            r3(ix,iy,iz)  = 0.
            r4(ix,iy,1,iz)= 0.
            r5(ix,iy,iz)  = 0.
         enddo 
         enddo 
      enddo
c
c ------------- set initial random field to be
c               divergence free
c
      idum = -1
      do iz=izs,ize
         if (iz.le.8) then
c
c ----------- ampv and ampt are max amplitudes of random 
c             velocity and temperature fields
c
         ampv = 0.5
         ampt = 0.1
c  
c ------- simple random field scaled between 0 and 1
c
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = ran1(idum)
         enddo
         enddo
c
         do iy=iys,iye
         do ix=1,nnx
            psix(ix,iy) = psi(ix,iy)
            psiy(ix,iy,izs) = psi(ix,iy)
         enddo
         enddo
         call xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
c
         vmaxx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
            if(vmag .gt. vmaxx) vmaxx = vmag
         enddo
         enddo
         facv = ampv/vmaxx
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz)   = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
            v(ix,iy,iz)   = v(ix,iy,iz) + psix(ix,iy)*facv
            t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
            e(ix,iy,iz)   = 1.0
         enddo
         enddo
         endif
c
c ---------- check divergence of initial field
c
         do iy=iys,iye
         do ix=1,nnx
            uxx(ix,iy)     = u(ix,iy,iz)
            vyy(ix,iy,izs) = v(ix,iy,iz)
         enddo
         enddo
         call xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         do iy=iys,iye
         do ix=1,nnx
            divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
         enddo
         enddo
         divz(iz) = divz(iz)*fnxy
c
c -------- end z loop
c
      enddo
c
      call mpi_sum_z(divz(1),i_root,myid,nnz,1)
c
      write(nprt,6000)
 6000 format(' check of divergence for initial state',/,
     +         ' iz ',5x,' divergence')
      write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100 format(i5,e15.6)
c
c ------------ fix for baroclinic and subsidence effects !!
c
c     do iz=izs,ize
c        ug(iz)=ugcont
c        vg(iz)=vgcont
c        if (.not.(ibrcl.eq.1)) go to 19988
c        if (.not.(iz.le.izi)) go to 19987
c        ug(iz)=0.
c        vg(iz)=0.
c 19987    continue
c 19988    continue
c        zz2=zz(iz)
c        wls(iz)=-divgls*zz2
c        if (.not.(iz.eq.1)) go to 19986
c        do ix=1,nnx
c        uls(ix)=divgls*(dx*float(ix-1)-xl*.5)
c        enddo
c     enddo
c     write(nprt,9)(uls(ix),ix=1,nnx)
c  9  format(1x,8e12.3)
c 19986 continue
c
      return
      end
      subroutine randoc
c
c -------- random initial conditions for an
c          ocean simulation
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +     psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +     vyy(nnx,iys:iye,izs:izs)
c
      izi=(5*nnz)/20
      zi=z(izi)
      tmixed = 283.0
      do iz=izs,ize
         if (iz.le.izi) then
            do iy=iys,iye
            do ix=1,nnx
               u(ix,iy,iz)   = ugcont-ugal
               v(ix,iy,iz)   = vgcont
               w(ix,iy,iz)   = 0.0
               t(ix,iy,1,iz) = tmixed
               e(ix,iy,iz)   = 0.0
            enddo
            enddo
         endif
         if (iz.gt.izi) then
            do iy=iys,iye
            do ix=1,nnx
               u(ix,iy,iz)   = ugcont-ugal
               v(ix,iy,iz)   = vgcont
               w(ix,iy,iz)   = 0.0
               t(ix,iy,1,iz) = tmixed + dtdzf(1)*(zz(iz)-zi)
               e(ix,iy,iz)   = 0.0
            enddo
            enddo
         endif
         do iy=iys,iye
         do ix=1,nnx
            w(ix,iy,iz)    = 0.0
            r1(ix,iy,iz)   = 0.0
            r2(ix,iy,iz)   = 0.0
            r3(ix,iy,iz)   = 0.0
            r4(ix,iy,1,iz) = 0.0
            r5(ix,iy,iz)   = 0.0
         enddo
         enddo
      enddo
c
c ------------- set initial random field to be
c               divergence free
c
      idum = -1
      do iz=izs,ize
      if (iz.le.4) then
c
c ----------- ampv and ampt are max amplitudes of random 
c             velocity and temperature fields
c
         ampv = 0.01
c        ampt = 0.00
         ampt = 0.0001
c  
c ------- simple random field scaled between 0 and 1
c
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = ran1(idum)
         enddo
         enddo
c
         do iy=iys,iye
         do ix=1,nnx
            psix(ix,iy) = psi(ix,iy)
            psiy(ix,iy,izs) = psi(ix,iy)
         enddo
         enddo
         call xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
c
         vmaxx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
            if(vmag .gt. vmaxx) vmaxx = vmag
         enddo
         enddo
         facv = ampv/vmaxx
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
            v(ix,iy,iz) = v(ix,iy,iz) + psix(ix,iy)*facv
            t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
            e(ix,iy,iz) = 0.0001
         enddo
         enddo
      endif
c
c ---------- check divergence of initial field
c
         do iy=iys,iye
         do ix=1,nnx
            uxx(ix,iy) = u(ix,iy,iz)
            vyy(ix,iy,izs) = v(ix,iy,iz)
         enddo
         enddo
         call xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         do iy=iys,iye
         do ix=1,nnx
            divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
         enddo
         enddo
         divz(iz) = divz(iz)*fnxy
c
c -------- end z loop
c
      enddo
c
      call mpi_sum_z(divz(1),i_root,myid,nnz,1)
c
      write(nprt,6000)
 6000 format(' check of divergence for initial state',/,
     +         ' iz ',5x,' divergence')
      write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100 format(i5,e15.6)
c
      do iz=izs,ize
         ug(iz)=ugcont
         vg(iz)=vgcont
      enddo
c
      return
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
      subroutine forcing
c
c ----------- update surface temperature based on a 
c             constant cooling rate
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      tsfcc(1) = t_surf_i - c_rate*time
c
      return
      end
      subroutine set_sav(istart)
c
      use pars
      use fields
      use con_data
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
      subroutine pbltop(itop)
c
c ---------- get estimate of pbl top
c
c            method = 0, min of wt flux
c                        (good for buoyancy cases)
c            method = 1, uw flux less than critical value
c                        (good for ekman cases)
c            method = 2, running t average exceeds criterion
c                        (good for neutral cases with capping
c                         inversions)
c            method = 3, maximum gradient in temperature field
c                        (good for finding local zi see jas paper)
c                        with minimum search height (sr. setup)
c
c ------------ if method uses average statistics then only root
c              process need find zi
c
      use pars
      use fields
      use con_data
      use con_stats
      real trun(maxnz)
      include 'mpif.h'
      real gradloc(2,nnx,nny), gradmax(2,nnx,nny)
      external get_zi
c
      if(method .le. 2 .and. l_root) then
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = 1.0
      if (method .le. 0 .or. method .gt. 2) then
         itop=1
         wttot=wtle(1,1)+wtsb(1,1)
         wtmin=wttot*sgn
         do iz=2,nnz
            wttot=(wtle(iz,1)+wtsb(iz,1))*sgn
            if (wttot.le.wtmin) then
               itop=iz
               wtmin=wttot
            endif
         enddo
         zi=z(itop)
      else if (method .eq. 1) then
         itop = 1
         crit = 0.05
         uwsf = utau*utau
         do iz=1,nnzm1
               uwtot = (uwle(iz) + uwsb(iz))**2 +
     $                 (vwle(iz) + vwsb(iz))**2
               uwtot = sqrt(uwtot)
               if(uwtot/uwsf .gt. crit) then
                  itop=iz
               endif
         enddo
         zi=z(itop)
      else if (method .eq. 2) then
         trun(1) = txym(1,1)
         do iz=2,nnz
             weight = z(iz-1)/z(iz)
             trun(iz) = trun(iz-1)*weight + (1.0-weight)*txym(iz,1)
         enddo
         itop = 1
         tcrit = 0.25
         if(iocean .eq. 1) tcrit = 0.1
         do iz=2,nnz
                if(txym(iz,1) .gt. (trun(iz) + tcrit)) then
                  itop = iz
                  go to 320
                endif
         enddo
  320    continue
         zi=z(itop)
      endif
      do iy=1,nny
      do ix=1,nnx
         gradmax(2,ix,iy) = zi
      enddo
      enddo
c
c ----------- use gradient method, every process computes
c
      elseif(method .eq. 3) then
c
c ---------------- get local zi from gradient in temperaure field
c
c     dz_i = dzu_i(izs+1)
c     do iy=1,nny
c     do ix=1,nnx
c        gradloc(1,ix,iy) = (t(ix,iy,1,izs+1) - t(ix,iy,1,izs))*dz_i
c        gradloc(2,ix,iy) = z(izs)
c     enddo
c     enddo
c
c ------- similar to zeroing the stat array in sr. mean_stat
c
      do iy=1,nny
      do ix=1,nnx
         gradloc(1,ix,iy) = 0.0
         gradloc(2,ix,iy) = z(iz_min)
      enddo
      enddo
c
c ------------- now all z in this process
c
      if(iz_min .le. ize) then
      do iz=max(izs,iz_min),ize
         izp1 = iz + 1
         do iy=iys,iye
         do ix=1,nnx
            grad = (t(ix,iy,1,izp1) - t(ix,iy,1,iz))*dzu_i(izp1)
            if(grad .gt. gradloc(1,ix,iy)) then
               gradloc(1,ix,iy) = grad
               gradloc(2,ix,iy) = z(iz)
            endif
         enddo
         enddo
      enddo
      endif
c
c     call mpi_reduce(gradloc,gradmax,2*nnx*nny,mpi_real8,ziloc,
c    +                i_root,mpi_comm_world,ierror)
c
c ----------- alternate version using already defined function in mpi
c             passes 2 real8 variables
c
      call mpi_reduce(gradloc,gradmax,nnx*nny,mpi_2double_precision,
     +                mpi_maxloc,i_root,mpi_comm_world,ierror)
c
c ------------ get average on root process
c
      if(l_root) then
         zi_avg = 0.0
         do iy=1,nny
         do ix=1,nnx
            zi_avg = zi_avg + gradmax(2,ix,iy)
         enddo
         enddo
         zi = zi_avg*fnxy
c        itop = nint(zi/dz)
      endif
c
      endif
c
c -------- send average zi everywhere
c
      call mpi_bcast(zi,1,mpi_real8,
     +              i_root,mpi_comm_world,ierr)
c
      if(iocean .ne. 1) then
         do iz=1,nnz
            if(zi .ge. z(iz) .and.
     +         zi .lt. z(iz+1)) itop = iz
         enddo
      else
         do iz=1,nnz
            if(zi .le. z(iz) .and.
     +         zi .gt. z(iz+1)) itop = iz
         enddo
      endif
c
c     if(l_root) write(6,7001) myid,zi,itop
 7001 format(' 7001 in pbltop myid = ',i4,' zi = ',e15.6,
     +       ' itop = ',i3)
c
      return
      end
      subroutine get_zi(gradmax,gradout,len,itype)
c
      use pars
      real gradmax(*), gradout(*)
c
c     write(nprt,2001) myid, len
c2001 format(' 2001 in get_zi myid = ',i4,' len = ',i8)
c     write(nprt,2002) (i,gradmax(i),gradmax(i+1),i=1,len,2)
c2002 format(' i ',5x,' grad ',5x,' location ',/,
c    +      (i5,2e15.6))
c
      do i=1,len,2
         if(gradmax(i) .gt. gradout(i)) then
              gradout(i)   = gradmax(i)
              gradout(i+1) = gradmax(i+1)
         endif
      enddo
c
      return
      end
      subroutine print(lu,iz_strt,iz_end)
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      write(lu,4000)
 4000 format(30X,' --- SOLUTION ---')
      write(lu,4100) it,time,dt,zi,tsfcc(1),uwsfc,vwsfc,wtsfc(1),
     +              zol,hol,ucfl, vcfl, wcfl,
     +              t_ref
 4100 format(' IT=',I7,5x,'TIME (s) = ',e15.8,',  DT(s) = ',e15.6,/,
     +       10x,'ZTOP = ',e15.6,
     +       ',  TSFC = ',e15.6,
     +       ',  UW = ',e15.6,',  VW = ',e15.6,/,10x,
     +       'WT = ',e15.6,',  ZL =',e15.6,
     +       ',  HL = ',e15.6,/,10x,'U_cfl = ',e15.6,
     +       ',  V_cfl = ',e15.6,',  W_cfl = ',e15.6,/,10x,
     +       'Theta Ref = ',e15.6)
      write(lu,4200)
 4200 format(//,20x,'--------- HORIZONTAL MEAN VALUES ---------- ',
     +       //,2x,'IZ',4x,'T_MEAN',7x,'T2_MEAN',6x,
     + 'DIVG',8X,'LE_KE',6X,'SGS_KE',7X,'LE_WT',6X,
     + 'SGS_WT',7X,'SHRZ',8X,'BUOY')
      do 19999 iz=iz_end,iz_strt,-1
         !write(lu,4300)iz,txym(iz,1)-t_ref,divz(iz),
	 write(lu,4300)iz,txym(iz,1)-t_ref, txym(iz,2),divz(iz),
     +              englez(iz),eavg(iz),wtle(iz,1),
     +              wtsb(iz,1),shrz(iz),buyz(iz)
 4300    format(1X,I3,e12.4,8e12.4)
19999 continue
      write(lu,4400)tsfcc(1),wtsfc(1)
 4400 format('  SURFACE VALUE: TXYM=',F8.2,'               WTSB=',E9.2)
      if(iocean .eq. 1) then
         write(lu,4500) stokess,udrift,vdrift
 4500    format(/,' STOKESS = ',e12.4,' UDRIFT = ',e12.4,
     +          ' VDRIFT = ',e12.4)
      endif
      write(lu,4600) (iz,uxym(iz)+ugal,vxym(iz),uwle(iz),
     +       uwsb(iz),vwle(iz),vwsb(iz),iz=iz_strt,iz_end)
 4600 format(//,' IZ',5x,' UXYM + UGAL',8x,' VXYM',10x,' UWLE',10x,
     +          ' UWSB',10x,' VWLE',10x,' VWSB'
     +       ,/,(1x,i4,6(3x,e15.6)))
      if(ivis .eq. 1) then
         write(lu,4800) xksurf, nmatch, viscon, vise
 4800    format(//,' XKSURF = ',e15.6,' NMATCH = ',i4,/,
     +             ' VISCON = ',e15.6,' VISE = ',e15.6)
!         write(lu,4700) (iz,dfac(iz),iz=iz_strt,iz_end)
! 4700    format(//,'   IZ',5x,'  DFAC',/,(1x,i4,3x,e15.6))
      endif
c
c --------------- output additional scalars
c
      if(nscl .eq. 2) then
      write(lu,5005)tsfcc(2),wtsfc(2)
 5005 format(/,'  SURFACE VALUE: TXYM(2) =',e15.6,' WTSFC(2) = ',e15.6)
      write(lu,5100) (iz,txym(iz,2),wtle(iz,2),
     +              wtsb(iz,2),iz=iz_strt,iz_end)
 5100 format(//,' IZ',5x,' SCALAR-1 MEAN',8x,' WS1LE',10x,
     +          ' WS1SB',10x
     +       ,/,(1x,i4,3(3x,e15.6)))
      else if (nscl .eq. 3) then
!      write(lu,5205)tsfcc(2),wtsfc(2),tsfcc(3),wtsfc(3)
 5205 format(/,'  SURFACE VALUE: TXYM(2) =',e15.6,' WTSFC(2) = ',e15.6,
     +       /,'  SURFACE VALUE: TXYM(3) =',e15.6,' WTSFC(3) = ',e15.6)
      write(lu,5200) (iz,txym(iz,2),txym(iz,3),wtle(iz,2),
     +    wtsb(iz,2),wtle(iz,3),wtsb(iz,3),iz=iz_strt,iz_end)
 5200 format(//,' IZ',5x,' SCALAR-1 MEAN',8x,' SCALAR-2 MEAN',10x,
     +          ' WS1LE',10x,' WS1SB',10x,' WS2LE',10x,' WS1SB'
     +       ,/,(1x,i4,6(3x,e15.6)))
      endif
 
      return
      end
      subroutine xy_stats
c
c ------------ get statistics 
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
      implicit none
c
c ------- indices for indexing array stat(.,.)
c         js = number of non-scalar stats
c         ns = number of scalar stats
c
      integer, parameter :: js = 19, ns = 5, nstat = js + ns*nscl
      real :: stat(1:nnz,nstat)

      integer :: ix,iy,iz,izm1,izp1,izp2,m1,m2,m3,m4,m5,i,l

      real :: sgn
      real :: rlim
      real :: RHtmp,Ttmp,mod_magnus,exner,func_p_base,func_rho_base
      real :: rhoa,msqrRH
c
c -------- stat(.,1) = u*u = ups
c          stat(.,2) = v*v = vps
c          stat(.,3) = w*w = wps
c          stat(.,4) = w**3 = wcube
c          stat(.,5) = w**4 = wfour
c          stat(.,6) = resolved tke at zw = englez
c          stat(.,7) = sgs e at zu = engsbz
c          stat(.,8) = sgs e at zw = eavg
c          stat(.,9) = resolved uw at zw = uwle
c          stat(.,10) = resolved vw at zw = vwle
c          stat(.,11) = partsrc(1) 
c          stat(.,12) = partsrc(2) 
c          stat(.,13) = partsrc(3) 
c          stat(.,14) = Tpsrc
c          stat(.,15) = Hpsrc
c          stat(.,16) = TEpsrc
c          stat(.,17) = relative humidity
c          stat(.,18) = temperature (NOT potential temp)
c          stat(.,19) = RH^2
c          stat(.,m1) = resolved scalar flux wt at zw = wtle
c          stat(.,m2) = resolved scalar flux ut at zw = utle
c          stat(.,m3) = resolved scalar flux vt at zw = vtle
c          stat(.,m4) = scalar t*t at zu = tps
c          stat(.,m5) = scalar t*t*t at zu = tcube
c
c --------- use a trick with mpi reduce over all z to get averages
c           by setting stat array = 0 for all z on each process
c
      do i=1,nstat
      do iz=1,nnz
         stat(iz,i) = 0.0
      enddo
      enddo
c
c -------- indices for scalars
c
      m1 = js
      m2 = js + nscl
      m3 = js + 2*nscl
      m4 = js + 3*nscl
      m5 = js + 4*nscl
c
      sgn = 1.0
      if(iocean .eq. 1 .and. iupwnd .eq. 1) sgn = -1.0
c
      do iz=izs,ize
c
      izp2 = iz + 2
      izp1 = iz + 1
      izm1 = iz - 1

      !Set the density based on whether to take base state into account
      if (iexner .eq. 1) then
         rhoa = func_rho_base(surf_p,tsfcc(1),zz(iz))
      else
         rhoa = surf_rho
      end if


      do iy=iys,iye
      do ix=1,nnx
         stat(iz,1) = stat(iz,1) + (u(ix,iy,iz) - uxym(iz))**2
         stat(iz,2) = stat(iz,2) + (v(ix,iy,iz) - vxym(iz))**2
         stat(iz,3) = stat(iz,3) + (w(ix,iy,iz) - wxym(iz))**2
         stat(iz,4) = stat(iz,4) + (w(ix,iy,iz) - wxym(iz))**3
         stat(iz,5) = stat(iz,5) + (w(ix,iy,iz) - wxym(iz))**4
         stat(iz,6) = stat(iz,6) + 
     +                ((w(ix,iy,iz)-wxym(iz))**2 +
     +                (0.5*(u(ix,iy,iz)-uxym(iz) + 
     +                      u(ix,iy,izp1)-uxym(izp1)))**2 +
     +                (0.5*(v(ix,iy,iz)-vxym(iz) + 
     +                      v(ix,iy,izp1)-vxym(izp1)))**2)*0.5
c
         stat(iz,7) = stat(iz,7) + 0.5*(e(ix,iy,iz)+e(ix,iy,izm1))
         stat(iz,8) = stat(iz,8) + e(ix,iy,iz)
         stat(iz,9) = stat(iz,9) + (w(ix,iy,iz)-wxym(iz))*
     +              0.5*((u(ix,iy,iz)-uxym(iz))+
     +                   (u(ix,iy,izp1)-uxym(izp1)))
         stat(iz,10) = stat(iz,10) + (w(ix,iy,iz)-wxym(iz))*
     +              0.5*((v(ix,iy,iz)-vxym(iz))+
     +                   (v(ix,iy,izp1)-vxym(izp1)))
         stat(iz,11) = stat(iz,11) + partsrc(ix,iy,iz,1)
         stat(iz,12) = stat(iz,12) + partsrc(ix,iy,iz,2)
         stat(iz,13) = stat(iz,13) + partsrc(ix,iy,iz,3)
         stat(iz,14) = stat(iz,14) + partTsrc(ix,iy,iz)
         stat(iz,15) = stat(iz,15) + partHsrc(ix,iy,iz)
         stat(iz,16) = stat(iz,16) + partTEsrc(ix,iy,iz)

         if (iexner .eq. 1) then
         Ttmp = t(ix,iy,1,iz)*
     +   exner(surf_p,func_p_base(surf_p,tsfcc(1),zz(iz)))
         else
         Ttmp = t(ix,iy,1,iz)
         end if

         RHtmp = Ttmp*
     +   t(ix,iy,2,iz)*Ru/Mw/mod_magnus(Ttmp)*rhoa*100.0

         stat(iz,17) = stat(iz,17) + RHtmp
         stat(iz,18) = stat(iz,18) + Ttmp
         stat(iz,19) = stat(iz,19) + RHtmp**2

      enddo
      enddo

c
c ------------ get scalar resolved fluxes and variances
c
      do l=1,nscl
         if(iupwnd .ne. 1 .or. iz .eq. nnz) then
            do iy=iys,iye
            do ix=1,nnx
               stat(iz,m1+l)=stat(iz,m1+l) +
     +               (w(ix,iy,iz)-wxym(iz))*
     +               0.5*(t(ix,iy,l,iz)-txym(iz,l) +
     +                    t(ix,iy,l,izp1)-txym(izp1,l))
            enddo
            enddo
         else
c
c ------------------- monotone fluxes
c
           do iy=iys,iye
           do ix=1,nnx
              stat(iz,m1+l) = stat(iz,m1+l) +
     +    amax1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,l,iz) +
     + rlim(t(ix,iy,l,izp1),t(ix,iy,l,iz),t(ix,iy,l,izm1))) +
     +    amin1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,l,izp1) +
     + rlim(t(ix,iy,l,iz),t(ix,iy,l,izp1),t(ix,iy,l,izp2)))
           enddo
           enddo
         endif
         stat(iz,m1+l)= sgn*stat(iz,m1+l)
c
c ------------ get horizontal scalar resolved fluxes 
c
         do iy=iys,iye
         do ix=1,nnx
            stat(iz,m2+l) = stat(iz,m2+l)+
     +               (u(ix,iy,iz)-uxym(iz))*
     +               (t(ix,iy,l,iz)-txym(iz,l)) 
            stat(iz,m3+l) = stat(iz,m3+l)+
     +               (v(ix,iy,iz)-vxym(iz))*
     +               (t(ix,iy,l,iz)-txym(iz,l)) 
         enddo
         enddo
c
c ------------------- scalar variances & higher moments
c
         do iy=iys,iye
         do ix=1,nnx
            stat(iz,m4+l) = stat(iz,m4+l) + 
     +                (t(ix,iy,l,iz) - txym(iz,l))**2
            stat(iz,m5+l) = stat(iz,m5+l) + 
     +                (t(ix,iy,l,iz) - txym(iz,l))**3
         enddo
         enddo
c
c ------ end scalar loop
c
      enddo
c
c ------ end z loop
c
      enddo
c
c -------- add partial sums and send it to all
c
      call mpi_sum_z(stat(1,1),i_root,myid,nnz*nstat,1)
c
c ------ fill arrays for printout and constant file
c
      do iz=1,nnz
c 
      ups(iz)    = stat(iz,1)*fnxy
      vps(iz)    = stat(iz,2)*fnxy
      wps(iz)    = stat(iz,3)*fnxy
      wcube(iz)  = stat(iz,4)*fnxy
      wfour(iz)  = stat(iz,5)*fnxy
      englez(iz) = stat(iz,6)*fnxy
      engsbz(iz) = stat(iz,7)*fnxy
      eavg(iz)   = stat(iz,8)*fnxy
      uwle(iz)   = stat(iz,9)*fnxy
      vwle(iz)   = stat(iz,10)*fnxy

      m1src(iz) = stat(iz,11)*fnxy
      m2src(iz) = stat(iz,12)*fnxy
      m3src(iz) = stat(iz,13)*fnxy
      uw_tot(iz) = uwle(iz) + uwsb(iz)
      vw_tot(iz) = vwle(iz) + vwsb(iz)
      Tpsrc(iz) = stat(iz,14)*fnxy
      Hpsrc(iz) = stat(iz,15)*fnxy
      TEpsrc(iz) = stat(iz,16)*fnxy
      RHxym(iz) = stat(iz,17)*fnxy
      tempxym(iz) = stat(iz,18)*fnxy
      RHmsqr(iz) = stat(iz,19)*fnxy



c
c ------------ get scalar resolved fluxes and variances
c
      do l=1,nscl
         wtle(iz,l)   = stat(iz,m1+l)*fnxy
         utle(iz,l)   = stat(iz,m2+l)*fnxy
         vtle(iz,l)   = stat(iz,m3+l)*fnxy
         tps(iz,l)    = stat(iz,m4+l)*fnxy
         tcube(iz,l)  = stat(iz,m5+l)*fnxy
         wt_tot(iz,l) = wtle(iz,l) + wtsb(iz,l)
      enddo
      enddo

      !Append this with summing pflux and pfluxdiff
      call mpi_sum_z(pflux(0),i_root,myid,nnz+1,1)
      call mpi_sum_z(pfluxdiff(0),i_root,myid,nnz+1,1)
      call mpi_sum_z(pmassflux(0),i_root,myid,nnz+1,1)
      call mpi_sum_z(penegflux(0),i_root,myid,nnz+1,1)

      !Make these fluxes with units per area per time:
      pflux = pflux/dt/yl/xl
      pfluxdiff = pfluxdiff/dt/yl/xl
      pmassflux = pmassflux/dt/yl/xl
      penegflux = penegflux/dt/yl/xl

      !Do the volume-average RH calculation
      meanRH = 0.0
      msqrRH = 0.0
      do iz=1,nnz-1

         meanRH = meanRH + 
     +   0.5*(zz(iz+1)-zz(iz))*(RHxym(iz)+RHxym(iz+1))

         msqrRH = msqrRH + 
     +   0.5*(zz(iz+1)-zz(iz))*(RHmsqr(iz)+RHmsqr(iz+1))

      end do 

      meanRH = meanRH/(zz(nnz)-zz(1))  !Note: divide by the distance between the points, not the entire zl, to get the true average
      msqrRH = msqrRH/(zz(nnz)-zz(1))
      varRH = msqrRH - meanRH**2
      

c
      return
      end
      subroutine tke_budget

!NOTE::: THIS HAS NOT BEEN UPDATED TO REFLECT HUMIDITY IN THE BUOYANCY
!EQUATION   DHR 7/5/16

c
c -------- get terms in resolved scale tke budget
c          as in gabls writeup at w-points
c          at istage = 1. 
c
      use pars
      use particles
      use fields
      use con_data
      use con_stats
c
      real :: stat(1:nnz,5)
      real :: s11s,s22s,s33s,s12s,s13s,s23s,wz,wzp
      real :: s13,s23,s33
      real :: ufluc,ufluc_t,ufluc_b,vfluc,vfluc_t,vfluc_b,wfluc
      real :: ffluc1,ffluc2,ffluc3
      real :: ffluc1p,ffluc2p,ffluc3p
      real :: weit,weit1
      integer :: iz,i,j,izp1,izm1
c
c -------- stat(.,1) = tke transport  = wq
c          stat(.,2) = pressure transport  = wp
c          stat(.,3) = tke dissipation
c          stat(.,4) = tke dissipation
c          stat(.,5) = particle force correlation
c
      do iz=1,nnz
         stat(iz,1) = 0.0
         stat(iz,2) = 0.0
         stat(iz,3) = 0.0
         stat(iz,4) = 0.0
         stat(iz,5) = 0.0
      enddo

!Compute DNS dissipation, since there is no subgrid dissipation now:
      do iz=izs,ize
         izp1 = iz + 1
         izm1 = iz - 1
         weit = dzw(iz)/(dzw(iz) + dzw(izp1))
         weit1 = 1.0 - weit
!
! ---- get fluctuating strain rates:
!      here, sij = 1/2*(duidxj + dujdxi)
!      then t_diss = 2*nu*<sij sij>
! ---- NOTE: these are computed at the w-locations!  (not u,v locations)
!
         t_diss(iz) = 0.0
         do j=iys,iye
         do i=1,nnx

            !Things for dissipation - these are computed at w-locations since
            !there is no z-derivative after
            s11s = weit1*ux(i,j,iz)**2 + weit*ux(i,j,izp1)**2
            s22s = weit1*vy(i,j,iz)**2 + weit*vy(i,j,izp1)**2
            wz  = (w(i,j,iz)-w(i,j,izm1))*dzw_i(iz)
            wzp = (w(i,j,izp1)-w(i,j,iz))*dzw_i(izp1)
            s33s = weit*wzp**2 + weit1*wz**2
            s12s = weit1*(0.5*(uy(i,j,iz) + vx(i,j,iz)))**2 +
     +            weit*(0.5*(uy(i,j,izp1) + vx(i,j,izp1)))**2
            s13s = (0.5*((u(i,j,izp1) - u(i,j,iz) +
     +            uxym(iz) - uxym(izp1))*dzu_i(izp1) +
     +            wx(i,j,iz)))**2
            s23s = (0.5*((v(i,j,izp1) - v(i,j,iz) +
     +          vxym(iz) - vxym(izp1))*dzu_i(izp1) +
     +          wy(i,j,iz)))**2

         stat(iz,3) = stat(iz,3) + 2.0*vis_m(i,j,iz)*(
     +               s11s+s22s+s33s+2.0*(s12s+s13s+s23s))

            !Things for viscous transport - these are computed at u-locations since
            !a z-derivative is done to the average, which lands t_tau on w-locations
            
            ufluc_t   = u(i,j,izp1) - uxym(izp1)
            ufluc   = u(i,j,iz) - uxym(iz)
            ufluc_b   = u(i,j,izm1) - uxym(izm1)
            vfluc_t   = v(i,j,izp1) - vxym(izp1)
            vfluc   = v(i,j,iz) - vxym(iz)
            vfluc_b   = v(i,j,izm1) - vxym(izm1)
            wfluc = 0.5*( (w(i,j,iz)-wxym(iz)) 
     +                  + (w(i,j,izm1)-wxym(izm1)) )

            uz_t = (ufluc_t-ufluc)*dzu_i(izp1)
            uz_b = (ufluc-ufluc_b)*dzu_i(iz)
            vz_t = (vfluc_t-vfluc)*dzu_i(izp1)
            vz_b = (vfluc-vfluc_b)*dzu_i(iz)
            
            uz = 0.5*(uz_t+uz_b)
            vz = 0.5*(vz_t+vz_b)

            s13 = 0.5*(uz + 0.5*(wx(i,j,iz)+wx(i,j,izm1)))
            s23 = 0.5*(vz + 0.5*(wy(i,j,iz)+wy(i,j,izm1)))
            s33 = wz

         !Note: just uses vis_m(1,1,iz) since it's equal everywhere:
         stat(iz,4) = stat(iz,4) + 2.0*vis_m(i,j,iz)*(
     +               ufluc*s13 + vfluc*s23 + wfluc*s33)
                    
         !Finally get the particle force correlation term:
         ffluc1 = partsrc(i,j,iz,1)-m1src(iz)
         ffluc2 = partsrc(i,j,iz,2)-m2src(iz)
         ffluc3 = partsrc(i,j,iz,3)-m3src(iz)
         if (iz==nnz) then
         ffluc1p = 0.0
         ffluc2p = 0.0
         ffluc3p = 0.0
         else
         ffluc1p = partsrc(i,j,izp1,1)-m1src(izp1)
         ffluc2p = partsrc(i,j,izp1,2)-m2src(izp1)
         ffluc3p = partsrc(i,j,izp1,3)-m3src(izp1)
         end if
         stat(iz,5) = stat(iz,5) +
     +       weit*(ufluc_t*ffluc1p)+weit1*(ufluc*ffluc1)+
     +       weit*(vfluc_t*ffluc2p) + weit1*(vfluc*ffluc2) + 
     +       (w(i,j,iz)-wxym(iz))*(weit*ffluc3p+weit1*ffluc3)
         enddo
         enddo
         stat(iz,3) = stat(iz,3)*fnxy
         stat(iz,4) = stat(iz,4)*fnxy
         stat(iz,5) = stat(iz,5)*fnxy
      enddo
c
c --------------- get transport terms as vertical arrays
c
      do iz=izs,ize
c
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
c
c --------- get estimate of turbulent transport term
c
            ufluc   = u(ix,iy,iz) - uxym(iz)
            vfluc   = v(ix,iy,iz) - vxym(iz)
            wfluc   = w(ix,iy,iz) - wxym(iz)
            wfluc_l = w(ix,iy,izm1) - wxym(izm1)
            stat(iz,1)  = stat(iz,1) + 0.25*(wfluc + wfluc_l)*
     +             (ufluc**2 + vfluc**2 + 0.5*(wfluc_l**2 + wfluc**2))
c
c --------- get estimate of pressure transport term
c
            pfluc = p(ix,iy,iz) - pxym(iz)
     +           -0.5*((u(ix,iy,iz)+stokes(iz))**2 + 
     +                 v(ix,iy,iz)*v(ix,iy,iz) + 
     +      0.5*(w(ix,iy,iz)*w(ix,iy,iz)+w(ix,iy,izm1)*w(ix,iy,izm1)))
            stat(iz,2) = stat(iz,2) + pfluc*0.5*(wfluc_l + wfluc)
         enddo
         enddo
         stat(iz,1) = stat(iz,1)*fnxy
         stat(iz,2) = stat(iz,2)*fnxy
      enddo
      call mpi_sum_z(stat(1,1),i_root,myid,nnz*5,1)
c
c ------ we have all terms on all processors for all z, add them up
c        treat tr_tau at bottom special, tr_tau = 0
c
!      tr_tau(0) = 0.0
      do iz=1,nnz
c
         izp1 = iz + 1
         izm1 = iz - 1
         if(iz .eq. nnz) then
            t_tau(iz) = 0.0
            t_wp(iz)  = 0.0
            t_wq(iz)  = 0.0
         else
            t_tau(iz) = (stat(izp1,4) - stat(iz,4))*dzu_i(izp1) 
            t_wq(iz)  = -(stat(izp1,1) - stat(iz,1))*dzu_i(izp1)
            t_wp(iz)  = -(stat(izp1,2) - stat(iz,2))*dzu_i(izp1)
         endif
         dudz = (uxym(izp1) - uxym(iz))*dzu_i(izp1)
         dvdz = (vxym(izp1) - vxym(iz))*dzu_i(izp1)
c
c ------------- gather all the budget terms
c
         t_tran(iz)  = t_wq(iz) + t_wp(iz) + t_tau(iz)
         t_rprod(iz) = -(dudz*uwle(iz) + dvdz*vwle(iz))
         !Old t_sprod had subgrid stuff
         !t_sprod(iz) =  (dudz*uwsb(iz) + dvdz*vwsb(iz))
         !Now make t_sprod the spray tke term to reduce new variables:
         t_sprod(iz) = -stat(iz,5)
         t_buoy(iz)  =  batag*wtle(iz,1)
         t_diss(iz) = stat(iz,3)
c
      enddo
c
      return
      end
      subroutine extra_flux_terms
      use pars
      use particles
      use fields
      use con_data
      use con_stats
      use fftwk
      implicit none

      real :: stat(1:nnz,3)
      real :: weit,weit1
      real :: Tflucp1,Tflucm1,Tfluc,Tmean
      real :: qflucp1,qflucm1,qfluc,qmean,dqpdz
      real :: Sq,Sqp1
      real :: gradTp(3),dTpdz1,dTpdz,dTdz
      integer :: iz,i,j,izp1,izm1,iscl

!     Compute the "extra" enthalpy budget terms, located at w locations:
c -------- stat(.,1) = < q' T' w' > 
c          stat(.,2) = < T' Sq > where Sq is q source
c          stat(.,3) = Dv*< T' dq'/dz > 

       

      stat = 0.0
      do iz=izs,ize
         izp1 = iz + 1
         izm1 = iz - 1
         weit = dzw(iz)/(dzw(iz) + dzw(izp1))
         weit1 = 1.0 - weit

       do j=iys,iye
       do i=1,nnx

         if (iz==1)  then
           Tmean = 2.0*Tbot(1) - txym(iz,1)
           Tflucm1 = t(i,j,1,izm1)-Tmean
           Tflucp1 = t(i,j,1,izp1)-txym(izp1,1)
         elseif (iz==nnz) then
           Tmean = 2.0*Ttop(1) - txym(iz,1)
           Tflucp1 = t(i,j,1,izp1)-Tmean
           Tflucm1 = t(i,j,1,izm1)-txym(izm1,1)
         else
           Tflucp1 = t(i,j,1,izp1)-txym(izp1,1)
           Tflucm1 = t(i,j,1,izm1)-txym(izm1,1)
         end if
         Tfluc = t(i,j,1,iz)-txym(iz,1)

         if (iz==1)  then
           qmean = 2.0*Tbot(2) - txym(iz,2)
           qflucm1 = t(i,j,2,izm1)-qmean
           qflucp1 = t(i,j,2,izp1)-txym(izp1,2)
           Sqp1 = partHsrc(i,j,izp1)
         elseif (iz==nnz) then
           qmean = 2.0*Ttop(2) - txym(iz,2)
           qflucp1 = t(i,j,2,izp1)-qmean
           qflucm1 = t(i,j,2,izm1)-txym(izm1,2)
           Sqp1 = 0.0
         else
           qflucp1 = t(i,j,2,izp1)-txym(izp1,2)
           qflucm1 = t(i,j,2,izm1)-txym(izm1,2)
           Sqp1 = partHsrc(i,j,izp1)
         end if
         qfluc = t(i,j,2,iz)-txym(iz,2)
         Sq = partHsrc(i,j,iz)

         !Get dq'/dz at w-locations:
         dqpdz = (qflucp1 - qfluc)*dzu_i(izp1)

         !Then can get source #3:
         stat(iz,3) = stat(iz,3) + vis_s(i,j,2,iz)*
     +   (weit1*Tfluc + weit*Tflucp1)*dqpdz


         !Now source #1, which must be evaluated at w-location:

         stat(iz,1) = stat(iz,1) + w(i,j,iz)*
     +   (weit1*Tfluc + weit*Tflucp1)*
     +   (weit1*qfluc + weit*qflucp1)


         !Now source #2:
         stat(iz,2) = stat(iz,2) + 
     +   (weit1*Tfluc + weit*Tflucp1)*
     +   (weit1*Sq + weit*Sqp1)

     
       end do
       end do

         stat(iz,1) = stat(iz,1)*fnxy
         stat(iz,2) = stat(iz,2)*fnxy
         stat(iz,3) = stat(iz,3)*fnxy
      end do
         

      call mpi_sum_z(stat(1,1),i_root,myid,nnz*3,1)


c
c ------ we have all terms on all processors for all z, add them up
c
      do iz=1,nnz
c
c ------------- gather all the budget terms
c
         trip(iz) = stat(iz,1)
         TpSq(iz) = stat(iz,2)
         Tpdqp(iz) = stat(iz,3)
      enddo
c
      return
      end
      subroutine Tvar_budget(iscl)
      use pars
      use particles
      use fields
      use con_data
      use con_stats
      use fftwk
      implicit none

      real :: stat(1:nnz,4)
      real :: weit,weit1
      real :: tx_tmp(nnx,iys:iye), ty(nnx,iys:iye,izs-1:ize+1)
      real :: tx(nnx,iys:iye,izs-1:ize+1)
      real :: trans(izs:ize)
      real :: gradTp(3),Tfluc,dTpdz1,dTpdz,dTdz
      real :: Tflucp1,Tflucm1,qfluc,qflucp1,Tmean
      integer :: iz,i,j,izp1,izm1,iscl

c -------- stat(.,1) = Transport: -del.[U<T'2> + <u T'2> - alpha*grad(T'2)]
c          stat(.,2) = Dissipation: -2*alpha <grad(T').grad(T')> 
c          stat(.,3) = Particle: <T' Q'>

      !Need the y gradient of temp:
      do iz=izs-1,ize+1
      do j=iys,iye
      do i=1,nnx
         ty(i,j,iz)  = t(i,j,iscl,iz)
         tx_tmp(i,j)  = t(i,j,iscl,iz)
      enddo
      enddo
      call xderivp(tx_tmp(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      tx(1:nnx,iys:iye,iz) = tx_tmp(1:nnx,iys:iye)
      enddo

      call yd_mpi(ty(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
       

      stat = 0.0
      Tv_part2 = 0.0
      do iz=izs,ize
         izp1 = iz + 1
         izm1 = iz - 1
         weit = dzw(iz)/(dzw(iz) + dzw(izp1))
         weit1 = 1.0 - weit

       do j=iys,iye
       do i=1,nnx

         if (iz==1)  then
         Tmean = 2.0*Tbot(iscl) - txym(iz,iscl)
         Tflucm1 = t(i,j,iscl,izm1)-Tmean
         Tflucp1 = t(i,j,iscl,izp1)-txym(izp1,iscl)
         elseif (iz==nnz) then
         Tmean = 2.0*Ttop(iscl) - txym(iz,iscl)
         Tflucp1 = t(i,j,iscl,izp1)-Tmean
         Tflucm1 = t(i,j,iscl,izm1)-txym(izm1,iscl)
         else
         Tflucp1 = t(i,j,iscl,izp1)-txym(izp1,iscl)
         Tflucm1 = t(i,j,iscl,izm1)-txym(izm1,iscl)
         end if
         Tfluc = t(i,j,iscl,iz)-txym(iz,iscl)

         !First dissipation: 
         !Note that gradients of total T and T' in x,y directions are equal since d<T>/dx = d<T>/dy = 0
         gradTp(1) = weit1*tx(i,j,iz) + weit*tx(i,j,izp1)
         gradTp(2) = weit1*ty(i,j,iz) + weit*ty(i,j,izp1)

         !Now get dT'/dz:
         gradTp(3) = (Tflucp1 - Tfluc)*dzu_i(izp1)
         
         stat(iz,2) = stat(iz,2)  - 2.0*vis_s(i,j,iscl,iz)*
     +                (gradTp(1)**2+gradTp(2)**2+gradTp(3)**2)

         !Next transport

         !Store the transport sum at u-locations since z-derivative at the end

         stat(iz,1) = stat(iz,1) + w(i,j,iz)*Tfluc**2

         dTpdz1 = (Tflucp1**2-Tfluc**2)*dzu_i(izp1)
         dTpdz = (Tfluc**2-Tflucm1**2)*dzu_i(iz)
         stat(iz,1) = stat(iz,1) - vis_s(i,j,iscl,iz)*0.5*(dTpdz1+dTpdz)

          !Particle source:
          if (iscl == 1) then
          if (iTcouple ==1) then
            qfluc = partTsrc(i,j,iz)-Tpsrc(iz)
            if (iz==nnz) then
            qflucp1 = 0.0
            else
            qflucp1 = partTsrc(i,j,izp1)-Tpsrc(izp1)
            end if
            stat(iz,3) = stat(iz,3) + 
     +         weit*(qflucp1*Tflucp1) + weit1*(qfluc*Tfluc)
          endif
          endif


          if (iscl == 2) then
          if (iHcouple ==1) then
            qfluc = partHsrc(i,j,iz)-Hpsrc(iz)
            if (iz==nnz) then
            qflucp1 = 0.0
            else
            qflucp1 = partHsrc(i,j,izp1)-Hpsrc(izp1)
            end if
            stat(iz,3) = stat(iz,3) + 
     +         weit*(qflucp1*Tflucp1) + weit1*(qfluc*Tfluc)
          endif
          endif

          if (iscl .eq. 2 .and. iHcouple .eq. 1 ) then
              qfluc = partTEsrc(i,j,iz)-TEpsrc(iz)
              if (iz==nnz) then
              qflucp1 = 0.0
              else
              qflucp1 = partTEsrc(i,j,izp1)-TEpsrc(izp1)
              end if
              Tv_part2(iz) = Tv_part2(iz) + 
     +          weit*(qflucp1*Tflucp1) + weit1*(qfluc*Tfluc)
  
         endif
     
       end do
       end do

         stat(iz,1) = stat(iz,1)*fnxy
         stat(iz,2) = stat(iz,2)*fnxy
         stat(iz,3) = stat(iz,3)*fnxy
         if (iscl .eq. 2 .and. iHcouple .eq. 1 ) then
         Tv_part2(iz) = Tv_part2(iz)*fnxy
         endif
      end do
         

      call mpi_sum_z(stat(1,1),i_root,myid,nnz*3,1)
      if (iscl .eq. 2 .and. iHcouple .eq. 1 ) then
      call mpi_sum_z(Tv_part2(1),i_root,myid,nnz,1)
      endif


c
c ------ we have all terms on all processors for all z, add them up
c
      do iz=1,nnz
         izp1 = iz + 1
         izm1 = iz - 1
         if(iz .eq. nnz) then
            Tv_tran(iz,iscl) = 0.0
         else
            Tv_tran(iz,iscl) = -(stat(izp1,1)-stat(iz,1))*dzu_i(izp1)
         endif
         dTdz = (txym(izp1,iscl)-txym(iz,iscl))*dzu_i(izp1)
c
c ------------- gather all the budget terms
c
         Tv_prod(iz,iscl) = -2.0*wtle(iz,iscl)*dTdz
         Tv_diss(iz,iscl) = stat(iz,2)
         Tv_part1(iz,iscl) = -stat(iz,3)
         !Tv_part2 is already gathered
      enddo
c
      return
      end
      subroutine write_his(iloc)
c
c ----- write history file with global parameters
c       write tsfcc specially to preserve digits!
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
c
      divgmax = 0.0
      do iz=1,nnz
         divgmax = amax1(divgmax, divz(iz))
      enddo
c
      ziavg = zi
      holtop = hol
      wt_min = wtsb(iloc,1)
      wt_le  = wtle(iloc,1)
      krec = krec + 1
      mid = nnz/4
      write(nhis1,6000) time,dt,utau,ziavg,amonin,holtop,
     +         (tsfcc(1)-t_ref),uwsfc,vwsfc,divgmax, wt_min, wt_le,
     +         ucfl, vcfl, wcfl, wtsfc(1),
     +         ups(mid),vps(mid),wps(mid),tps(mid,1),
     +         uwle(mid),uwsb(mid),uw_tot(mid),
     +         vwle(mid),vwsb(mid),vw_tot(mid),
     +         wtle(mid,1),wtsb(mid,1),wt_tot(mid,1),
     +         englez(mid),eavg(mid), wabs,
     +         lwc,
     +         Rep_avg,radavg,radmin,radmax

 6000 format(37e17.8)
c
c -------------- write profile information
c
      call write_prof(nhisp,krec,isize,c_s%wwsb)
c
      return
      end
      subroutine write_prof(nhisp,krec,num,f)
      real f(num)
      real*4 f32(num)
c
c -------------- build special 32 bit arrays for profiles
c
      do i=1,num
         f32(i) = f(i)
      enddo
c
      write(nhisp,err=999,rec=krec) (f32(i),i=1,num)
c
      return
c --------------- errors
  999 continue
      write(6,9000) num,krec
 9000 format(' 9000, trouble in ',
     +       'SR. save_prof cannot write profile data ',/,
     +       ' num = ',i8, 'krec = ',i6)
      stop
      end
      subroutine close_histograms
c
c ---- close history files
c
      use pars
      logical there
c
c ---- root closes and checks the files
c
      close(nrad)
      close(nres)
      close(nactres)
      inquire(file=path_rad_hist,exist=there)
      if(.not.there) then
         write(6,8100) path_rad_hist
         call mpi_finalize(ierr)
         stop
      endif
      inquire(file=path_res_hist,exist=there)
      if(.not.there) then
         write(6,8200) path_res_hist
         call mpi_finalize(ierr)
         stop
      endif
      inquire(file=path_actres_hist,exist=there)
      if(.not.there) then
         write(6,8300) path_actres_hist
         call mpi_finalize(ierr)
         stop
      endif
      write(6,7100) path_rad_hist
      write(6,7200) path_res_hist
      write(6,7300) path_actres_hist
c
      return
c -------------------- process write errors
 7100 format(' RADHIST DATA IS WRITTEN IN FILE  ',a80)
 7200 format(' RESHIST DATA IS WRITTEN IN FILE  ',a80)
 7300 format(' ACTRESHIST DATA IS WRITTEN IN FILE  ',a80)
 8100 format(' SR. SAVE_RADPDF: Trouble rad histogram file =',a80)
 8200 format(' SR. SAVE_RESPDF: Trouble res histogram file =',a80)
 8300 format(' SR. SAVE_ACTRESPDF: Trouble res histogram file =',a80)
      end
      subroutine close_his
c
c ---- close history files
c
      use pars
      logical there
c
c ---- root closes and checks the files
c
      close(nhis1)
      close(nhisp)
      inquire(file=path_sav_h,exist=there)
      if(.not.there) then
         write(6,8000) path_sav_h
         call mpi_finalize(ierr)
         stop
      endif
      inquire(file=path_sav_hp,exist=there)
      if(.not.there) then
         write(6,8100) path_sav_hp
         call mpi_finalize(ierr)
         stop
      endif
      write(6,7000) path_sav_h
      write(6,7100) path_sav_hp
c
      return
c -------------------- process write errors
 7000 format(' HISTORY DATA IS WRITTEN IN FILE  ',a80)
 7100 format(' PROFILE HISTORY DATA IS WRITTEN IN FILE  ',a80)
 7200 format(' RADHIST DATA IS WRITTEN IN FILE  ',a80)
 8000 format(' SR. SAVE_HIS: Trouble history file not in path =',a80)
 8100 format(' SR. SAVE_HIS: Trouble profile history file',
     +       ' not in path =',a80)
 8200 format(' SR. SAVE_RADPDF: Trouble rad histogram file =',a80)
      end
      subroutine dealias
c
c --------- wave cutoff filter using 2d fft
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real wve(nny,jxs:jxe,izs:ize)
      real wves(nnxp2,iys:iye,izs:ize)
c
c --------- sharp spectral cutoff, specific to current 2dfft
c
      ix_cut   = 2*int(float(nnx)/3.) + 3
      iy_cut_l = int(float(nny)/3.) + 2
      iy_cut_u = nnyp2 - iy_cut_l
c
c ---------- u-equation
c
      call fft2d_mpi(u(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(u(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ---------- v-equation
c
      call fft2d_mpi(v(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(v(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ---------- w-equation
c
      call fft2d_mpi(w(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(w(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ---------- e-equation
c
      call fft2d_mpi(e(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(e(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ------------- scalars, not stored in correct order
c
      do iscl=1,nscl
         do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            wves(ix,iy,iz) = t(ix,iy,iscl,iz)
         enddo
         enddo
         enddo
         call fft2d_mpi(wves(1,iys,izs),wve(1,jxs,izs),trigx(1,1),
     +           trigc,nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
         call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
         call fft2d_mpi(wves(1,iys,izs),wve(1,jxs,izs),trigx(1,1),
     +           trigc,nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
         do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,iscl,iz) = wves(ix,iy,iz)
         enddo
         enddo
         enddo
      enddo
c
      return
      end
      subroutine sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
c
c --------- sharp cutoff filter for field wve stored
c           in 2d-fft order
c
      use pars
      real wve(nny,jxs:jxe,izs:ize)
c
      do iz=izs,ize
         do ix=jxs,jxe
         do iy=iy_cut_l,iy_cut_u
            wve(iy,ix,iz) = 0.0
         enddo
         enddo
      enddo
c
      if(jxe .lt. ix_cut) go to 999
c
         do iz=izs,ize
            do ix=max(jxs,ix_cut),jxe
            do iy=1,nny
               wve(iy,ix,iz) = 0.0
            enddo
            enddo
         enddo
c
  999 continue
c
      return
      end
      subroutine gridd
c
c ----------- allocate space and pass arrays using modules
c
      use pars
      use fields
      use particles 
      use fftwk
      use con_data
      use con_stats

      if (myid==0) write(6,5001) isize
 5001 format(' size of stats array = ',i8)

c
c
c ---------------- debug for arrays
c
      big = -99.0e+300
c
c ---------------- setup grid
c
      nnx = maxnx
      nny = maxny
      nnz = maxnz
c     izs = 1
c     ize = nnz
c
c
c ----------- make sure problem and cpu's match
c
      maxp   = numprocs-1
      ncpu_z = numprocs/ncpu_s
      if(mod(numprocs,ncpu_s) .ne. 0 .or.
     +   ncpu_z .gt. nnz) then
         go to 999
      endif
      if(l_root) write(6, 1100) ncpu_s, ncpu_z, numprocs,
     +                          maxp
      write(nprt,1100) ncpu_s, ncpu_z, numprocs, maxp
 1100 format(' Number of x-y slab cpus = ',i5,/,
     +       ' Number of z-level cpus  = ',i5,/,
     +       ' Total number of cpus    = ',i5,/,
     +       ' Max-p for index arrays  = ',i5)
c
c ---------------- allocate arrays for (i,j,k)-indexing on
c                  each processor (see set_range)
c
      allocate(ix_s(0:maxp), ix_e(0:maxp),
     +         jx_s(0:maxp), jx_e(0:maxp),
     +         kx_s(0:maxp), kx_e(0:maxp),
     +         mx_s(0:maxp), mx_e(0:maxp),
     +         iy_s(0:maxp), iy_e(0:maxp),
     +         jy_s(0:maxp), jy_e(0:maxp),
     +         is_s(0:maxp), is_e(0:maxp),
     +         iz_s(0:maxp), iz_e(0:maxp))
c
c ---------------- setup array sizes and variable dimensions
c
      nxy   = nnx*nny
      ncx   = nnx/2 + 1
      ncy   = nny/2 + 1
      nnxp1 = nnx + 1
      nnyp1 = nny + 1
      nnxp2 = nnx + 2
      nnyp2 = nny + 2
      nnzp1 = nnz + 1
      nnzm1 = nnz - 1
      ivis = ivis0
      fnxy  = 1.0/float(nnx*nny)
c
      write(nprt,7001) nnx,nny,nnz
 7001 format(' 7001 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)
c
      call set_range
c
      write(nprt,7002) nnx,nny,nnz
 7002 format(' 7002 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)
c
      num_y = iye + 1 - iys
c
c ------------- allocate solution arrays
c               account for nnxp2 for fields but not in rhs
c               and possible monotone for scalars
c
      allocate(u(nnxp2,iys:iye,izs-1:ize+1), 
     +         v(nnxp2,iys:iye,izs-1:ize+1), 
     +         w(nnxp2,iys:iye,izs-1:ize+1), 
     +         t(nnxp2,iys:iye,nscl,izs-2:ize+2), 
     +         e(nnxp2,iys:iye,izs-1:ize+1), 
     +         r1(nnx,iys:iye,izs-1:ize+1),
     +         r2(nnx,iys:iye,izs-1:ize+1),
     +         r3(nnx,iys:iye,izs-1:ize+1),
     +         r4(nnx,iys:iye,nscl,izs-1:ize+1),
     +         r5(nnx,iys:iye,izs-1:ize+1))
      u = 0.0
      v = 0.0
      w = 0.0
      t = 0.0
      e = 0.0
      r1 = 0.0
      r2 = 0.0
      r3 = 0.0
      r4 = 0.0
      r5 = 0.0

      allocate(p_base(0:nnz+1),
     +         rho_base(0:nnz+1),
     +         T_base(0:nnz+1),
     +         theta_base(0:nnz+1))
      p_base = 0.0
      rho_base = 0.0
      T_base = 0.0
      theta_base = 0.0

c   ------Allocate arrays for SFS stress and its derivatives
      allocate(sigm_s(nnxp2,iys:iye,izs-1:ize+1),
     +         sigm_sdx(nnxp2,iys:iye,izs-1:ize+1),
     +         sigm_sdy(nnxp2,iys:iye,izs-1:ize+1),
     +         sigm_sdz(nnxp2,iys:iye,izs-1:ize+1),
     +         vis_ss(nnxp2,iys:iye,izs-1:ize+1))

      allocate (sigm_sext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         sigm_sdxext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         sigm_sdyext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         sigm_sdzext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         vis_sext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3))

c
c ------------- allocate extended arrays for interpolation of
c               particle/spray location
c
      allocate(uext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3), 
     +         vext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3), 
     +         wext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         Text(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         T2ext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3)) 
      uext = 0.0
      vext = 0.0
      wext = 0.0
      Text = 0.0
      T2ext = 0.0

      !Transposed velocities to do the uf interpolation:
      allocate(u_t(0:nnz+1,iys:iye,mxs:mxe),
     +         v_t(0:nnz+1,iys:iye,mxs:mxe),
     +         w_t(0:nnz+1,iys:iye,mxs:mxe), 
     +         T_t(0:nnz+1,iys:iye,mxs:mxe),
     +         T2_t(0:nnz+1,iys:iye,mxs:mxe))
      u_t = 0.0
      v_t = 0.0
      w_t = 0.0
      T_t = 0.0
      T2_t = 0.0

      !Column-oriented array to calculate max. droplet relaxation time
      allocate(partsrc(nnx,iys:iye,izs-1:ize+1,1:3))
      partsrc = 0.0
      allocate(partsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1,1:3))
      partsrc_t = 0.0
      allocate(partTsrc(nnx,iys:iye,izs-1:ize+1))
      partTsrc = 0.0
      allocate(partTsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1))
      partTsrc_t = 0.0
      allocate(partHsrc(nnx,iys:iye,izs-1:ize+1))
      partHsrc = 0.0
      allocate(partHsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1))
      partHsrc_t = 0.0
      allocate(partTEsrc(nnx,iys:iye,izs-1:ize+1))
      partTEsrc = 0.0
      allocate(partTEsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1))
      partTEsrc_t = 0.0

      allocate(LWP(mxs:mxe,iys:iye))
      LWP = 0.0

      allocate(surf_precip(mxs:mxe,iys:iye))
      surf_precip = 0.0

      radsrc = 0.0
      viz_t_elapsed = 1.0e-10 !Initialize very small to prevent blowing up on 1st time step

c ------------- allocate space for boundary condition arrays
c               on top and bottom of domain
c
      allocate(ubc(nnx,iys:iye,2),
     +         vbc(nnx,iys:iye,2),
     +         wbc(nnx,iys:iye,2),
     +         tbc(nnx,iys:iye,nscl,2),
     +         ebc(nnx,iys:iye,2),
     +         pbc(nnx,iys:iye,2),
     +         pbc2(nnx,iys:iye,2))
c
c ------------ allocate space for wind and surface arrays
c
      allocate(wind(nnx,iys:iye), 
     +         tau13m(nnx,iys:iye), 
     +         tau23m(nnx,iys:iye), 
     +         taut3m(nnx,iys:iye,nscl), 
     +         t_grnd(nnx,iys:iye,nscl))
c
c ------------------- allocate space for derivative arrays
c
      allocate(ux(nnx,iys:iye,izs-1:ize+1),
     +         uy(nnx,iys:iye,izs-1:ize+1),
     +         vx(nnx,iys:iye,izs-1:ize+1),
     +         vy(nnx,iys:iye,izs-1:ize+1),
     +         wx(nnx,iys:iye,izs-1:ize+1),
     +         wy(nnx,iys:iye,izs-1:ize+1))
c
c ------------- allocate space for pressure, pressure bcs
c
      allocate(p(nnxp2,iys:iye,izs-1:ize+1),
     +         ptop(nnxp2,iys:iye,2))
      p = 0.0
      ptop = 0.0
c
c ------------- allocate space for viscosity and diffusivity
c
      allocate(vis_m(nnx,iys:iye,izs-1:ize+1),
     +         vis_s(nnx,iys:iye,nscl,izs-1:ize+1))
c
c ------------- allocate space for fft trig factors
c
      nq_trig = max(nnx,nny)
      allocate(trigx(2*nq_trig+15,2),
     +         trigc(4*nq_trig+15))
c
      return
  999 continue
c
      if(l_root) write(6,1000) numprocs, ncpu_s, mmz
      write(nprt,1000) numprocs, ncpu_s, nnz
 1000 format(' Gridd Trouble number of processors and grid',
     +          ' partitioning do not match!',/,
     +          ' Total num of cpus   = ',i5,
     +          ' Num cpu on x-y slab = ',i5,/,
     +          ' Num of z-levels     = ',i5)
      call mpi_finalize(ierr)
      end
      subroutine restart
c
c ----------- get restart file from local directory
c
      use pars
      use fields
      use con_data
      use con_stats
      character*80 path_res_c
      logical there
c
c --------------------- check if file is there
c
      inquire(file=path_res,exist=there)
      if(there) then
         if(l_root) write(6,6001) path_res
      else
         if(l_root) write(6,6005) path_res
         stop
      endif
c
c ------------------ get constant file
c
      iloc = index(path_res,' ')
      path_res_c = path_res(1:iloc-1)//'.con'
      inquire(file=path_res_c,exist=there)
      if(there) then
         if(l_root) write(6,6002) path_res_c
      else
         if(l_root) write(6,6006) path_res_c
         stop
      endif
      open(nvelc,err=200,file=path_res_c,form='unformatted',
     +        status='old')
c
      call read_res
c
      return
c ---------------------------- process errors
  100 continue
      write(6,9000) path_res, nvel
      call mpi_finalize(ierr)
      stop
c -----------------------
  200 continue
      write(6,9001) path_res_c, nvelc
      call mpi_finalize(ierr)
      stop
c -----------------------
 6001 format(' SR. RESTART: FILE READ = ',A80)
 6002 format(' SR. RESTART: CONSTANT FILE READ = ',A80)
 6005 format(' 6005, SR. RESTART: cannot find restart file = ',a80)
 6006 format(' 6005, SR. RESTART: cannot find constant file = ',a80)
 9000 format(' 9000, SR. RESTART: cannot open file =',a80,/,
     +       ' to unit number = ',i2)
 9001 format(' 9001, SR. RESTART: cannot open file =',a80,/,
     +       ' to unit number = ',i2)
      end
      subroutine read_res
c
c -------------- read restart file including constant file
c                changed for iys:iye
c
      use pars
      use fields
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
      allocate(temp(nvar,nnx,iys:iye))
c
c ---- open file
c
      call mpi_file_open(mpi_comm_world, path_res,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, nvel, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ------------ read 3d fields
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
c ------------ every mpi process reads constant file
c
      rewind(nvelc)
      read(nvelc,err=9993) c_c, c_s, c_i, case
      close(nvelc)
c
      if(l_root) write(6,4001) case
 4001 format(' 4001, SR. RESTART: case from restart = ',a3)
c
c ----- special restart conditions -------------------------------------
c
c -------- set case name to case input
c
      case   = case_inp
      if(l_root) write(6,4002) case_inp, utau, utausv
 4002 format(' 4002, SR. RESTART:',/,
     +       ' files will be saved with case name = ',a3,/,
     +       ' utau = ',e15.6,' utausv = ',e15.6)
c
c ------------------- if new vis model set match point for
c                     outer grid
      nmatch = 48
      utau = utausv
c
c -------- hand coded changes to restart if needed
c
c       qstars = 0.000
c       wtsfcs = 0.000
c
c
c ------   reset qstar and wtsfc for no heat flux
c              qstar(1) = qstars
c              wtsfc(1) = wtsfcs
c              qstar(2) = qstars
c              wtsfc(2) = wtsfcs
c ------   redefine case id to input value
c              case = cases
c
      if(l_root) write(6,4012) time
      if(l_root) write(6,4013) qstar(1) , nmatch, case
c
      call get_dz
c
      return
c ------------------------  process errors from read
c9991 continue
c     write(6,6000) nvel,iz
c6000 format(' SR. READ_RES: hit end of file on unit number = ',i2,/,
c    +       '               at iz = ',i4)
c     call mpi_finalize(ierr)
c     stop
c ---------------------
 9992 continue
      write(6,6100) nvel,iz
 6100 format(' SR. READ_RES: error reading file on unit number = ',i2,/,
     +       '               at iz = ',i4)
      call mpi_finalize(ierr)
      stop
c ---------------------
 9993 continue
      write(6,6200) nvelc
 6200 format(' SR. READ_RES:',/,
     +       '    error reading constant file on unit number = ',i2)
      call mpi_finalize(ierr)
      stop
c ---------------------
 4012 format(' SR. RESTART: restart completed at T=',e15.6)
 4013 format('    after restart qstar = ',e15.6,' nmatch = ',i5,
     +       ' case = ',a3)
      end
      subroutine save_viz
c
c --------------- save multiple (x-y), (x-z), (y-z), planes of data .
c                 modify recl in all open statements for more or less
c                 variables. 
c                 Constant - x, implies yz planes
c                 Constant - y, implies xz planes
c                 Constant - z, implies xy planes
c
c ------------- routine uses send/recv to get information in y-z planes
c
      use pars
      use fields
      use con_data
      use con_stats
      use fftwk
      use particles
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
c
      parameter(nvar_o = 6)
c
      integer ix_pick(maxnx),  iy_pick(maxny),  iz_pick(maxnz), 
     +        ix_order(maxnx), iy_order(maxny), iz_order(maxnz)
c
      integer istatus(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
c
      real(kind=4), dimension(nvar_o,nny,izs:ize) :: temp_x
      real(kind=4), dimension(nvar_o,nnx,izs:ize) :: temp_y
      real(kind=4), dimension(nvar_o,nnx,iys:iye) :: temp_z
      real, dimension(nvar_o,iys:iye,izs:ize)     :: buf_send
c
c ------------- don't touch
c
      data iviz_x,  iviz_y,  iviz_z  /0, 0, 0/
      data ionce_x, ionce_y, ionce_z, istuff /0, 0, 0, 0/
      data ix_pick, iy_pick, iz_pick /maxnx*0, maxny*0, maxnz*0/
      data ix_order, iy_order, iz_order /maxnx*0, maxny*0, maxnz*0/
      save iviz_x,  iviz_y,  iviz_z, 
     +     ix_pick, iy_pick, iz_pick, 
     +     ix_order, iy_order, iz_order,
     +     ionce_x, ionce_y, ionce_z, istuff,
     +     npln_x, npln_y, npln_z
c
c
c ----------- turn on z levels to save. Customize for your own use.
c             Set iz_pick(iz) = iz, ix_pick(ix) = ix, iy_pick(iy) = iy
c             Data is round-robin alternated in the data file for more than
c             1 plane for any particular view.
c

      iz_pick(12) = 12
      iz_pick(28) = 28
      iz_pick(64) = 64
      !iz_pick(20) = 20
      !iz_pick(45) = 45
      !iz_pick(60) = 60
c
c -------------- pick an x-z plane of data (can add more)
c
      iy_pick(nny/2) = nny/2
c     iy_pick(nny)   = nny
c
c -------------- pick a y-z plane of data (can add more)
c
      ix_pick(nnx/2) = nnx/2
c     ix_pick(nnx)   = nnx
c
c ------ find total number of z's turned on and open file once
c
      if (ionce_z .eq. 0) then
         npln_z = 0
         do k=1,nnz
            if(iz_pick(k) .eq. k) then
               npln_z = npln_z + 1
               iz_order(k) = npln_z
            endif
         enddo
         ionce_z = 1
         iviz_z =  -npln_z
         if (npln_z .ne. 0) then
      call mpi_barrier(mpi_comm_world,ierr)
            call mpi_file_open(mpi_comm_world, path_viz_xy,
     +                         mpi_mode_create+mpi_mode_rdwr,
     +                         mpi_info_null, nviz_z, ierr)
            disp = 0
            call mpi_file_set_view(nviz_z,disp,mpi_real4,mpi_real4,
     +                            'native',mpi_info_null,ierr)
         endif
      endif
c
c ------ find total number of y's turned on and open file once
c
      if (ionce_y .eq. 0) then
         npln_y = 0
         do j = 1,nny
            if(iy_pick(j) .eq. j) then
               npln_y = npln_y + 1
               iy_order(j) = npln_y
            endif
         enddo
         ionce_y = 1
         iviz_y  = -npln_y
         if (npln_y .ne. 0) then
            call mpi_file_open(mpi_comm_world, path_viz_xz,
     +                         mpi_mode_create+mpi_mode_rdwr,
     +                         mpi_info_null, nviz_y, ierr)
            disp = 0
            call mpi_file_set_view(nviz_y,disp,mpi_real4,mpi_real4,
     +                            'native',mpi_info_null,ierr)
         endif
      endif
c
c ------ find total number of x's turned on and open file once
c
      if (ionce_x .eq. 0) then
         npln_x = 0
         do i=1,nnx
            if(ix_pick(i) .eq. i) then
               npln_x = npln_x + 1
               ix_order(i) = npln_x
            endif
         enddo
         ionce_x = 1
         iviz_x  = -npln_x
         if (npln_x .ne. 0) then
            call mpi_file_open(mpi_comm_world, path_viz_yz,
     +                         mpi_mode_create+mpi_mode_rdwr,
     +                         mpi_info_null, nviz_x, ierr)
            disp = 0
            call mpi_file_set_view(nviz_x,disp,mpi_real4,mpi_real4,
     +                            'native',mpi_info_null,ierr)
         endif
      endif
c
      if(istuff .eq. 0 .and. l_root) then
         open(nviz_s,file=path_stuf)
         istuff = 1
      endif
c
c --------- write data, subtract t_ref to increase
c           resolution on 32 bit machines
c
c ---------- xy planes of data
c
      iviz_z  = iviz_z + npln_z
      nsize   = int(nvar_o,k8)*nnx*nny
      nsize2  = int(nvar_o,k8)*nnx*(iys-1)
      n_write = nvar_o*nnx*(iye+1-iys)
      do k=izs,ize
         if(iz_pick(k) .eq. k) then
            km1 = k - 1
            do j=iys,iye
            do i=1,nnx
               temp_z(1,i,j) = u(i,j,k)
               temp_z(2,i,j) = v(i,j,k)
               temp_z(3,i,j) = w(i,j,k)
               temp_z(4,i,j) = (t(i,j,1,k) - t_ref)
c
c ---------- get just the fluctuating pressure field
c
               temp_z(5,i,j) = p(i,j,k) - pxym(k)
     +                        -(e(i,j,k) + e(i,j,km1))/3.0
     +                        -0.5*((u(i,j,k) + stokes(k))**2 +
     +                               v(i,j,k)*v(i,j,k) +
     +                         0.5*(w(i,j,k)*w(i,j,k) + 
     +                              w(i,j,km1)*w(i,j,km1)))

               temp_z(6,i,j) = u(i,j,k)-uxym(k)
               !temp_z(6,i,j) = partsrc(i,j,k,1)
               !temp_z(7,i,j) = partsrc(i,j,k,2)
               !temp_z(8,i,j) = partsrc(i,j,k,3)
            enddo
            enddo
#if defined(SWAP)
            call byteswap(temp_z)
#endif
            offset = int((iviz_z + iz_order(k) - 1),k8)*nsize + nsize2
            call mpi_file_write_at(nviz_z,offset,temp_z,n_write,
     +                             mpi_real4,istatus,ierr)
            if (ierr .ne. 0) go to 9991
         endif
      enddo
c
c ---------- xz planes of data
c
      iviz_y = iviz_y + npln_y
      nsize  = int(nvar_o,k8)*nnx*nnz
      nsize2 = int(nvar_o,k8)*nnx*(izs-1)
      nwrite = nvar_o*nnx*(ize+1-izs)
      do j=iys,iye
         if(iy_pick(j) .eq. j) then
            do k=izs,ize
            km1 = k - 1
            do i=1,nnx
               temp_y(1,i,k) = u(i,j,k)
               temp_y(2,i,k) = v(i,j,k)
               temp_y(3,i,k) = w(i,j,k)
               temp_y(4,i,k) = (t(i,j,1,k) - t_ref)
c
c ---------- get just the fluctuating pressure field
c
               temp_y(5,i,k) =  p(i,j,k) - pxym(k)
     +                          -(e(i,j,k)+e(i,j,km1))/3.0
     +                          -0.5*((u(i,j,k)+stokes(k))**2 +
     +                               v(i,j,k)*v(i,j,k) +
     +                           0.5*(w(i,j,k)*w(i,j,k) + 
     +                                w(i,j,km1)*w(i,j,km1)))

               temp_y(6,i,k) = u(i,j,k) - uxym(k)
               !temp_y(6,i,k) = partsrc(i,j,k,1)
               !temp_y(7,i,k) = partsrc(i,j,k,2)
               !temp_y(8,i,k) = partsrc(i,j,k,3)
            enddo
            enddo
#if defined(SWAP)
            call byteswap(temp_y)
#endif
            offset = int((iviz_y + iy_order(j) - 1),k8)*nsize + nsize2
            call mpi_file_write_at(nviz_y,offset,temp_y,nwrite,
     +                                mpi_real4,istatus,ierr)
            if (ierr .ne. 0) goto 9992
         endif
      enddo
c
c ---------- yz planes that cut across all processors
c            just have root node on that slab write data
c
      iviz_x  = iviz_x + npln_x
      n_write = nvar_o*nny*(ize+1-izs)
      nsize   = int(nvar_o,k8)*nny*nnz
      nsize2  = int(nvar_o,k8)*nny*(izs-1)
      n_send  = nvar_o*(ize+1-izs)*(iye+1-iys)
      do i=1,nnx
         if(ix_pick(i) .eq. i) then
c
c ----------- build send buffer
c
            do k=izs,ize
            km1 = k - 1
            do j=iys,iye
               buf_send(1,j,k) = u(i,j,k)
               buf_send(2,j,k) = v(i,j,k)
               buf_send(3,j,k) = w(i,j,k)
               buf_send(4,j,k) = (t(i,j,1,k) - t_ref)
c
c ---------- get just the fluctuating pressure field
c
               buf_send(5,j,k) = p(i,j,k) - pxym(k)
     +                          -(e(i,j,k) + e(i,j,km1))/3.0
     +                          -0.5*((u(i,j,k) + stokes(k))**2 +
     +                               v(i,j,k)*v(i,j,k) +
     +                           0.5*(w(i,j,k)*w(i,j,k) + 
     +                                w(i,j,km1)*w(i,j,km1)))
              buf_send(6,j,k) = u(i,j,k)-uxym(k)
              !buf_send(6,j,k) = partsrc(i,j,k,1)
              !buf_send(7,j,k) = partsrc(i,j,k,2)
              !buf_send(8,j,k) = partsrc(i,j,k,3)
            enddo
            enddo
            if(myid .ne. iss) then
              call mpi_send(buf_send(1,iys,izs),n_send,
     +                      mpi_real8,iss,1,
     +                      mpi_comm_world,ierr)
            else
              do k=izs,ize
              do j=iys,iye
              do ii=1,nvar_o
                 temp_x(ii,j,k) = buf_send(ii,j,k)
              enddo
              enddo
              enddo
              do l=iss+1,ise
                 call recv_yz_var(temp_x,nvar_o,nny,
     +                            iy_s(l),iy_e(l),izs,ize,l)
              enddo
#if defined(SWAP)
              call byteswap(temp_x)
#endif
              offset = int((iviz_x + ix_order(i) - 1),k8)*nsize + nsize2
              call mpi_file_write_at(nviz_x,offset,temp_x,n_write,
     +                          mpi_real4,istatus,ierr)
              if (ierr .ne. 0) goto 9993
            endif
         endif
      enddo
c
c ------------- ascii file with facts in it that goes
c               with visualization
c
      if(l_root) then
         write(nviz_s,5000) time, amonin, zi, utau
 5000    format(4e20.8)
      endif
c
c ---- last time step close the files
c
!      if (it .eq. itmax) then
!         call mpi_file_close(nviz_z, ierr)
!         call mpi_file_close(nviz_y, ierr)
!         call mpi_file_close(nviz_x, ierr)
!         if (l_root) then
!            close(nviz_s)
!         endif
!      endif
       if (it .eq. itmax .or. mtape) then
        if(npln_z .ne. 0) then
            call mpi_file_close(nviz_z, ierr)
            ionce_z = 0
         endif
         if(npln_y .ne. 0) then
            call mpi_file_close(nviz_y, ierr)
            ionce_y = 0
         endif
         if(npln_x .ne. 0) then
            call mpi_file_close(nviz_x, ierr)
            ionce_x = 0
         endif
         if (l_root) then
            close(nviz_s)
            istuff = 0
         endif
      endif
c
      return
c --------------------------  errors in writing viz file
 9991 continue
      write(6,6000) nviz_z, iz
 6000 format(' SR. SAVE_VIS:',/,
     +       '    trouble cannot write xy viz file on unit = ',i2,/,
     +       '             at iz = ',i4)
      call mpi_finalize(ierr)
      stop
c --------------------------  errors in writing viz file
 9992 continue
      write(6,6100) nviz_y, iz, iviz_y
 6100 format(' SR. SAVE_VIS:',/,
     +       '    trouble cannot write xz viz file on unit = ',i2,/,
     +       '             at iz = ',i4,/,
     +       '            iviz_y = ',i8)
c --------------------------  errors in writing viz file
 9993 continue
      write(6,6200) nviz_x, iz, iviz_x
 6200 format(' SR. SAVE_VIS:',/,
     +       '    trouble cannot write yz viz file on unit = ',i2,/,
     +       '             at iz = ',i4,/,
     +       '            iviz_x = ',i8)
      call mpi_finalize(ierr)
      stop
      end
      subroutine recv_yz_var(temp_x,nvar,nny,iys,iye,izs,ize,ir)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real buf(nvar,iys:iye,izs:ize)
      real(kind=4), dimension(nvar,nny,izs:ize) :: temp_x
c
      num = nvar*(ize+1-izs)*(iye+1-iys)
      call mpi_recv(buf(1,iys,izs),num,mpi_real8,ir,1,
     +             mpi_comm_world,istatus,ierr)
      do k=izs,ize
      do j=iys,iye
      do ii=1,nvar
         temp_x(ii,j,k) = buf(ii,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine save_v
c
c --------------- save 3d fields
c
      use pars
      use fields
      use con_data
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
      logical there
c
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)                 nsize, nsize2
c
      real, allocatable, dimension(:,:,:) :: temp
      allocate(temp(nvar,nnx,iys:iye))
c
c ---- open file
c
      call mpi_file_open(mpi_comm_world, path_sav_v,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, nvel, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ---- write data
c
      nsize   = int(nvar,k8)*nnx*nny
      nsize2  = int(nvar,k8)*nnx*(iys-1)
      n_write = nvar*nnx*(iye+1-iys)
c
      do k=izs,ize
         do j = iys,iye
         do i = 1,nnx
            temp(1,i,j)    = u(i,j,k)
            temp(2,i,j)    = v(i,j,k)
            temp(3,i,j)    = w(i,j,k)
            temp(nvar,i,j) = e(i,j,k)
         enddo
         enddo
         do is = 1,nscl
            do j = iys,iye
            do i = 1,nnx
               temp(3+is,i,j) = t(i,j,is,k)
            enddo
            enddo
         enddo


#if defined(SWAP)
      call byteswap(temp)
#endif
c

         offset = int((k-1),k8)*nsize + nsize2
c         call mpi_file_write_at(nvel,offset,temp,n_write,
c     +                              mpi_real8,status,ierr)
         call mpi_file_write_at_all(nvel,offset,temp,n_write,
     +                              mpi_real8,status,ierr)
         if (ierr /= 0) goto 9991
c
      enddo

c
c ---- close file
c
      call mpi_file_close(nvel, ierr)

c
c ---- check file
c
      if (l_root) then
         inquire(file=path_sav_v,exist=there)
         if(.not.there) then
            write(6,8000) nvel,myid
            call mpi_finalize(ierr)
            stop
         endif
         write(6,7000) it,path_sav_v
      endif
c
      deallocate(temp)
c
      return
c --------------------------  errors in writing restart file
 9991 continue
      write(6,6000) nvel, iz
 6000 format(' SR. SAVE_V:',/,
     +       '    trouble cannot write restart file on unit = ',i2,/,
     +       '             at iz = ',i4)
      call mpi_finalize(ierr)
      stop
c --------------------
 7000 format(' **** DATA SET AT IT = ',I6,/,
     +       '      VELOCITY DATA IS WRITTEN IN FILE  ',a80)
 8000 format(' in SAVE_V: trouble writing file ',i5,'  myid = ',i5,
     +       ' at iz = ',i5)
      end
      subroutine save_c
c
c --------------- root process writes constant file
c                 sequential fortan binary
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      logical there
      character options*8, passwd*1
c
c ---- open file
c
      open(nvelc,err=9992,file=path_sav_c,form='unformatted',
     +                status='unknown')
      write(nvelc,err=9992) c_c, c_s, c_i, case
      close(nvelc)
c
        inquire(file=path_sav_c,exist=there)
        if(.not.there) then
           write(6,8001) path_sav_c
           call mpi_finalize(ierr)
           stop
        endif
c -----------------------------  output ok message
      write(6,7001) path_sav_c
c
      return
c --------------------------  errors in writing constant file
 9992 continue
      write(6,6100) nvelc
 6100 format(' SR. SAVE_V:',/,
     +  '    trouble cannot open/write constant file on unit = ',i2)
      call mpi_finalize(ierr)
      stop
c ---------------------
 7001 format('      CONSTANT DATA IS WRITTEN IN FILE  ',a80)
 8001 format(' SR. SAVE_C: Trouble constant file not in path =',a80)
      end
      subroutine save_p
c
c -------------- save pressure file
c
      use pars
      use fields
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
      logical there
c
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
c
      real temp(nnx,iys:iye)
c
c ---- open file
c
      call mpi_file_open(mpi_comm_world, path_sav_p,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, npre, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(npre,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ---- write data
c
      nsize   = int(nnx,k8)*nny
      nsize2  = int(nnx,k8)*(iys -1)
      n_write = nnx*(iye+1-iys)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            temp(i,j) = p(i,j,k)
         enddo
         enddo
#if defined(SWAP)
      call byteswap(temp)
#endif
         offset = int((k-1),k8)*nsize + nsize2
         call mpi_file_write_at_all(npre,offset,temp,n_write,
     +                              mpi_real8,status,ierr)
c         call mpi_file_write_at(npre,offset,temp,n_write,
c     +                              mpi_real8,status,ierr)
      enddo
c
c ---- close file
c
      call mpi_file_close(npre, ierr)
c
c ---- check file
c
      if (l_root) then
         inquire(file=path_sav_p,exist=there)
         if(.not.there) then
            write(6,8000) path_sav_p
            call mpi_finalize(ierr)
            stop
         endif
         write(6,7000) path_sav_p
      endif
c
      return
c -------------------- process write errors
 9991 continue
      write(6,6000) npre, iz
 6000 format(' SR. SAVE_P:',/,
     +       '    trouble cannot write pressure file on unit = ',i2,/,
     +       '             at iz = ',i4)
      call mpi_finalize(ierr)
      stop
c -----------------------
 7000 format('      PRESSURE DATA IS WRITTEN IN FILE  ',a80)
 8000 format(' SR. SAVE_P: Trouble pressure file not in path =',a80)
      end
      subroutine get_units
c
      use pars
c
c -------------- unit numbers for files
c
      nvel  = 20 
      npre  = 30
      nhis1 = 40
      nrad = 45
      nres = 46
      nactres = 47
      nvelc = 50
      nhisp = 60
      nviz_z = 80
      nviz_y = 82
      nviz_x = 84
      nviz_s = 90 
c
c ------------- unit number for standard print out
c               for each mpi task
c
      nprt = 1 
c
c ------------- open unit for standard printout
c
      path_prt = case_inp(1:3)//'.mp.xxxxx.out'
      write(path_prt(8:12),'(i5.5)') myid
      open(nprt,file=path_prt,form='formatted')
c
      return
      end
      subroutine get_output_filenames
c
c ----------- build file names for velocity, pressure, and constants
c
      use pars
      use con_data
      include 'mpif.h'
      character cgrid*10, num*3
c
c --------------- build character strings for file name
c
      cgrid = '.le.'
      write(num,'(i3.3)') itn
      iblnk = index(path_seed,' ')
      call blnk(path_sav_v)
      call blnk(path_sav_p)
      call blnk(path_sav_c)
      call blnk(path_sav_part)
      path_sav_v = path_seed(1:iblnk-1)//'/u'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)
      path_sav_p = path_seed(1:iblnk-1)//'/p'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)
      path_sav_c = path_seed(1:iblnk-1)//'/u'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)//'.con'
      path_sav_part = path_seed(1:iblnk-1)//'/part'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)
c
      return
      end
      subroutine open_histograms(istep)
c
c ------------------- open history files by root
c                     isize determined in sr. fill_cs
c
      use pars
      include 'mpif.h'
      character cgrid*4, iblks*16, num*3
      logical there
     
c
c --------------- build character strings for ascii history file name
c
      !Open PDF files
      cgrid = '.le.'

      call blnk(iblks)
      write(iblks(1:7),'(i7.7)') istep
      iblks(8:8) = '_'
      write(iblks(9:15),'(i7.7)') (istep+itape)
      iblnk = index(path_seed,' ')
      call blnk(path_rad_hist)
      call blnk(path_res_hist)


      path_rad_hist = path_seed(1:iblnk-1)//'/radhist'//
     +       cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'      

      path_res_hist = path_seed(1:iblnk-1)//'/reshist'//
     +       cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'      

      path_actres_hist = path_seed(1:iblnk-1)//'/actreshist'//
     +       cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'      


c
c ----------------- save data in directory
c
      if(l_root) then


      close(nrad)
      open(nrad,err=4500,file=path_rad_hist,form='formatted')

      close(nres)
      open(nres,err=4600,file=path_res_hist,form='formatted')

      close(nactres)
      open(nactres,err=4700,file=path_actres_hist,form='formatted')

      endif
c
      return
c ------------------- process errors
c-------------------
 4500 continue
      write(6,6303) nrad, path_rad_hist
 6303 format(' 6303, SR. OPEN_RADPDF:',/,
     +       '    cannot open radhist file on unit = ',i2,/,
     +       '    path = ',a80)
      stop
 4600 continue
      write(6,6303) nres, path_res_hist
 6304 format(' 6303, SR. OPEN_RESPDF:',/,
     +       '    cannot open reshist file on unit = ',i2,/,
     +       '    path = ',a80)
      stop
 4700 continue
      write(6,6305) nactres, path_actres_hist
 6305 format(' 6305, SR. OPEN_ACTRESPDF:',/,
     +       '    cannot open actreshist file on unit = ',i2,/,
     +       '    path = ',a80)
      stop
      end
      subroutine open_his(istep)
c
c ------------------- open history files by root
c                     isize determined in sr. fill_cs
c
      use pars
      include 'mpif.h'
      character cgrid*4, iblks*16, num*3
      logical there
     
c
c --------------- build character strings for ascii history file name
c
      cgrid = '.le.'
      call blnk(iblks)
      write(iblks(1:7),'(i7.7)') istep
      iblks(8:8) = '_'
      write(iblks(9:15),'(i7.7)') (istep + itape)
      iblnk = index(path_seed,' ')
      call blnk(path_sav_h)
      path_sav_h = path_seed(1:iblnk-1)//'/his'//
     +         cgrid(1:4)//case(1:3)//'.'//iblks(1:15)
c
c --------------- build character strings for ieee profile history file
c                 set record counter for direct access file = 0
c
      krec = 0
      cgrid = '.le.'
      call blnk(iblks)
      write(iblks(1:7),'(i7.7)') istep
      iblks(8:8) = '_'
      write(iblks(9:15),'(i7.7)') (istep + itape)
      iblnk = index(path_seed,' ')
      call blnk(path_sav_hp)
      path_sav_hp = path_seed(1:iblnk-1)//'/his'//
     +         cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.ieee'

c
c ----------------- save data in directory
c
      if(l_root) then

      close(nhis1)
      open(nhis1,err=3000,file=path_sav_h,form='formatted')
c
      close(nhisp)
      open(nhisp,err=4000,file=path_sav_hp,
     +        form='unformatted',access='direct',recl=isize*j_recl,
     +        status='unknown')

      end if

c
      return
c ------------------- process errors
 3000 continue
      write(6,6301) nhis1, path_sav_h
 6301 format(' 6301, SR. OPEN_HIS:',/,
     +       '    cannot open history1 file on unit = ',i2,/,
     +       '    path = ',a80)
      stop
c-------------------
 4000 continue
      write(6,6302) nhisp, path_sav_hp
 6302 format(' 6302, SR. OPEN_HIS:',/,
     +       '    cannot open history profile file on unit = ',i2,/,
     +       '    path = ',a80)
      stop

      end
      subroutine viz_output_filename(istep)
c
c ------------------- set visualization files,
c                     leaves files in scratch directory 
c
      use pars
      include 'mpif.h'
      character iblks*16

c
c --------------- build character strings for file names
c                 with time step
c
      call blnk(iblks)
      iblks(1:1) = '.'
      write(iblks(2:8),'(i7.7)') istep
      iblks(9:9) = '_'
      write(iblks(10:16),'(i7.7)') (istep + itape)
c
      iloc = index(path_seed,' ')
      path_viz_xy = path_seed(1:iloc-1)
     +         //'viz.'//case(1:3)//iblks(1:16)//'.xy.data'
c
      path_viz_xz = path_seed(1:iloc-1)
     +         //'viz.'//case(1:3)//iblks(1:16)//'.xz.data'
c
      path_viz_yz = path_seed(1:iloc-1)
     +         //'viz.'//case(1:3)//iblks(1:16)//'.yz.data'
c
      path_stuf = path_seed(1:iloc-1)
     +         //'stuff.'//case(1:3)//iblks(1:16)//'.data'
c
c     if(l_root) then
c        write(6,8001) path_viz_xy
c8001    format(' 8001: viz file = ',a80)
c        write(6,8001) path_viz_xz
c        write(6,8001) path_viz_yz
c        write(6,8001) path_stuf
c        write(6,8001) path_seed
c     endif
c
      return
      end
      subroutine open_viz
c
c ------------------- set visualization files,
c                     leaves files in scratch directory 
c
      use pars
      include 'mpif.h'
      character iblks*16
c
c --------------- build character strings for file names
c                 with time step
c
      call blnk(iblks)
      iblks(1:1) = '.'
      write(iblks(2:8),'(i7.7)') iti
      iblks(9:9) = '_'
      write(iblks(10:16),'(i7.7)') itmax
c
      iloc = index(path_viz_xy,' ')
      path_viz_xy = path_viz_xy(1:iloc-1)
     +      //'/viz.'//case(1:3)//iblks(1:16)//'.xy.data'
      iloc = index(path_viz_xz,' ')
      path_viz_xz = path_viz_xz(1:iloc-1)
     +      //'/viz.'//case(1:3)//iblks(1:16)//'.xz.data'
      iloc = index(path_viz_yz,' ')
      path_viz_yz = path_viz_yz(1:iloc-1)
     +      //'/viz.'//case(1:3)//iblks(1:16)//'.yz.data'
      path_stuf = path_stuf(1:iloc-1)
     +      //'/stuff.'//case(1:3)//iblks(1:16)//'.data'
      close(nviz_z)
      close(nviz_y)
      close(nviz_x)
      close(nviz_s)
c
c ----------- do not actually open the files here since
c             not all processors may have been picked and
c             its unknown how many variables are selected.
c             customized in sr. save_viz 
c
      return
      end
      subroutine range(n1,n2,nprocs,irank,ista,iend)
c
c ---------- the ibm range finder to balance load
c
      iwork1 = (n2 - n1 + 1)/nprocs
      iwork2 = mod(n2 - n1 +1, nprocs)
      ista = irank*iwork1 + n1 + min(irank,iwork2)
      iend = ista + iwork1 - 1
      if(iwork2 .gt. irank) iend = iend + 1
c
      return
      end
      subroutine set_range
c
c ---- build special x,y,z-ranges. dimensioned for 0:numprocs-1
c      indexed with myid
c
c      [ix_s:ix_e] x-range for computing y-derivatives nx-pts/ncpu_s
c                  in xtoy and ytox tranposes
c
c      [jx_s:jx_e] x-range for computing 2d fft (nx+2)-pts/ncpu_s
c                  must be even in each x-interval for complex fft in y
c
c      [kx_s:kx_e] x-range for pressure solver transpose (nx+2)-pts/ncpu_z
c                  nx+2 fourier coefficients for xtoz and ztox transposes
c
c      [mx_s:mx_e] x-range split across z cpus as nx-pts/ncpu_z
c                  for use in surface layer routines 
c
c      [is_s:is_e] starting and ending processor id's for a
c                  particular z-level
c
c      [iy_s:iy_e] y-range for computing y-derivatives ny-pts/ncpu_s
c                  in xtoy and ytox tranposes
c
c      [jy_s:jy_e] y-range for use in xtoz and ztox transposes 
c                  in pressure solution
c
c      [iz_s:iz_e] z-range for a particular vertical slab
c
c
      use pars
c
      write(nprt,7002) nnx,nny,nnz
 7002 format(' 7002 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)
c
      ii = -1
      do nn=0,ncpu_z-1
         call range(1,nnx+2,ncpu_z,nn,lx_s,lx_e)
         call range(1,nnx,ncpu_z,nn,nx_s,nx_e)
         call range(1,nny,ncpu_z,nn,ly_s,ly_e)
         call range(1,nnz,ncpu_z,nn,mz_s,mz_e)
         do mm=0,ncpu_s-1
            call range(1,nny,ncpu_s,mm,ny_s,ny_e)
            call range(1,nnx,ncpu_s,mm,nxy_s,nxy_e)
            call range(1,ncx,ncpu_s,mm,l2x_s,l2x_e)
            ii       = ii + 1
c
            ix_s(ii) = nxy_s
            ix_e(ii) = nxy_e
            jx_s(ii) = (l2x_s - 1)*2 + 1
            jx_e(ii) = l2x_e*2
            kx_s(ii) = lx_s
            kx_e(ii) = lx_e
            mx_s(ii) = nx_s
            mx_e(ii) = nx_e
c
            iy_s(ii) = ny_s
            iy_e(ii) = ny_e
            jy_s(ii) = ly_s
            jy_e(ii) = ly_e
c
            iz_s(ii) = mz_s
            iz_e(ii) = mz_e
c
            is_s(ii) = (ii/ncpu_s)*ncpu_s
            is_e(ii) = is_s(ii) + ncpu_s - 1
         enddo
      enddo
c
      iys = iy_s(myid)
      iye = iy_e(myid)
      jys = jy_s(myid)
      jye = jy_e(myid)
      ixs = ix_s(myid)
      ixe = ix_e(myid)
      jxs = jx_s(myid)
      jxe = jx_e(myid)
      kxs = kx_s(myid)
      kxe = kx_e(myid)
      mxs = mx_s(myid)
      mxe = mx_e(myid)
      izs = iz_s(myid)
      ize = iz_e(myid)
c
c ----------- get starting and  ending processor id's on each
c             vertical slab
c
      iss = is_s(myid)
      ise = is_e(myid)
c
c ------------ debug ranges 
c
      if(l_debug) then
         write(nprt,1200) myid, (nn, ix_s(nn), ix_e(nn), jx_s(nn),
     +                     jx_e(nn), kx_s(nn), kx_e(nn),
     +                     nn = 0,numprocs-1)
 1200    format(' myid =  ',i4,/,
     +       ' nn',5x,' ixs ',5x,' ixe ',5x,' jxs ',5x,' jxe '
     +       ,5x,' kxs ',5x,' kxe',/,(7i6))
c
         write(nprt,1213) myid, (nn, iy_s(nn), iy_e(nn),
     +                   jy_s(nn), jy_e(nn),
     +                   iz_s(nn), iz_e(nn), is_s(nn), is_e(nn),
     +                   nn=0,numprocs-1)
 1213    format(' myid = ',i3,/,
     +       ' nn ',3x,' iys ',5x,' iye ',5x,
     +       ' jys ',5x,' jye ',5x,
     +       ' izs ',5x,' ize',5x,' iss ',5x,' ise ',/,
     +       (9i6))
      endif
c
      return
      end
      subroutine mpi_sum_xy(f,myid,iss,ise,nsend)
c
c --------- get horizontal x-y sum over a set of proccessors [iss:ise]
c           for vector f(i). f(i) is overwritten. skip if single processor
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real work(nsend,iss:ise), f(nsend)
c
      if(iss .eq. ise) go to 999
c
      do j=1,nsend
         work(j,myid) = f(j)
         f(j)         = 0.0
      enddo
      do i=iss,ise
         if(i .ne. myid) then
            call mpi_sendrecv(work(1,myid),nsend,mpi_real8,i,1,
     +               work(1,i),nsend,mpi_real8,i,1,
     +           mpi_comm_world,istatus,ierr)
         endif
      enddo
      do i=iss,ise
      do j=1,nsend
         f(j) = f(j) + work(j,i)
      enddo
      enddo
c
  999 continue
c
      return
      end
      subroutine mpi_sum_z(f,i_root,myid,nsend,iall)
c
c --------- get sums on root or all processors
c           for all z for vector f(i)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real recv_b(nsend), f(nsend)
c
c -------- just root gets the result
c
      if(iall .ne. 1) then
         call mpi_reduce(f(1),recv_b(1),nsend,mpi_real8,mpi_sum,i_root,
     +                  mpi_comm_world,ierr)
         if(myid .eq. i_root) then
            do i=1,nsend
               f(i) = recv_b(i)
            enddo
         endif
      else
c
c -------- everyone gets the result
c
         call mpi_allreduce(f(1),recv_b(1),nsend,mpi_real8,mpi_sum,
     +                  mpi_comm_world,ierr)
         do i=1,nsend
            f(i) = recv_b(i)
         enddo
      endif
c
      return
      end
      subroutine mpi_sum_z_s(f,i_root,myid,nsend,nscl,iall)
c
c --------- get sums on root or all processors
c           for all z for vector f(i,nscl)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real recv_b(nsend,nscl), f(nsend,nscl)
c
      if(iall .ne. 1) then
         call mpi_reduce(f(1,1),recv_b(1,1),nsend*nscl,mpi_real8,
     +        mpi_sum,i_root,mpi_comm_world,ierr)
         if(myid .eq. i_root) then
            do iscl=1,nscl
            do i=1,nsend
               f(i,iscl) = recv_b(i,iscl)
            enddo
            enddo
         endif
      else
         call mpi_allreduce(f(1,1),recv_b(1,1),nsend*nscl,mpi_real8,
     +        mpi_sum, mpi_comm_world,ierr)
         do iscl=1,nscl
         do i=1,nsend
            f(i,iscl) = recv_b(i,iscl)
         enddo
         enddo
      endif
c
      return
      end
      subroutine mpi_gath_root(fs,fr,iz_s,iz_e,izs,ize,nz,myid,np,ns)
c
c ---------- gather results on root processors
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      integer iz_s(0:np-1), iz_e(0:np-1)
      real fs(izs:ize), fr(1:nz)
c
      if(np .eq. 1) go to 999
c
      irow_r = mod(myid,ns)
      if(myid .gt. ns) then
        call mpi_send(fs(izs),ize+1-izs,mpi_real8,irow_r,1,
     +       mpi_comm_world,ierr)
      else
        do l=irow_r+ns,np-1,ns
           ind = iz_s(l) + 1
           num = iz_e(l) + 1 - iz_s(l)
           call mpi_recv(fr(ind),num,mpi_real8,l,1,
     +       mpi_comm_world,istatus,ierr)
        enddo
      endif
c
  999 continue
c
      return
      end
      subroutine mpi_send_root(fs,num,myid,np,ns)
c
c ---------- send root results to other processors above it
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real fs(num)
c
      if(np .eq. 1) go to 999
c
      irow_r = mod(myid,ns)
      if(myid .ge. ns) then
        call mpi_recv(fs(1),num,mpi_real8,irow_r,1,
     +       mpi_comm_world,istatus,ierr)
      else
        do l=irow_r+ns,np-1,ns
           call mpi_send(fs(1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
  999 continue
c
      return
      end
      subroutine xtoy_trans(f,g,nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,
     +           myid,ncpu_s,np)
c 
c ------- transpose array  f(nx,iys:iye,iz1:iz2) ---> g(ny,ixs:ixe,iz1:iz2)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,iz1:iz2), 
     +     g(ny,ixs:ixe,iz1:iz2)
      real ft(nx*(iye+1-iys)*(iz2 - iz1 + 1)),
     +     gt(ny*(ixe+1-ixs)*(iz2 - iz1 + 1))
      integer ix_s(0:np-1), ix_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      jk = (iye - iys + 1)*(iz2 - iz1 + 1)
      ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)
c
c ----------- get cpus on slab for myid
c
      islab = myid/ncpu_s
      iss   = islab*ncpu_s
      ise   = iss + ncpu_s - 1
c
      do i=iss,ise
         nsend = (ix_e(i) - ix_s(i) + 1)*jk
         nrecv = (iy_e(i) - iy_s(i) + 1)*ik
         if(i .eq. myid) then
            call send_xtoy(f,gt(1),nx,ix_s(i),ix_e(i),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)
         else
            call send_xtoy(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)
            call mpi_sendrecv(
     +          ft(1),nsend,mpi_real8,i,1,
     +          gt(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_xtoy(g,gt(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(i),iy_e(i),iz1,iz2)
      enddo
c
      return
      end
      subroutine send_xtoy(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
c
c ------------- grab correct chunk of array to be sent
c
      real f(nx,iys:iye,izs:ize), ft(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         ft(i,j,k) = f(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_xtoy(g,gt,ny,ixs,ixe,iys,iye,izs,ize)
      real g(ny,ixs:ixe,izs:ize), gt(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         g(j,i,k) = gt(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine ytox_trans(g,f,nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,
     +           myid,ncpu_s,np)
c 
c ------- transpose array g(ny,ixs:ixe,iz1:iz2) ---> f(nx,iys:iye,iz1:iz2)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,iz1:iz2), 
     +     g(ny,ixs:ixe,iz1:iz2)
      real ft(nx*(iye+1-iys)*(iz2 - iz1 + 1)),
     +     gt(ny*(ixe+1-ixs)*(iz2 - iz1 + 1))
c
      integer ix_s(0:np-1), ix_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      jk = (iye - iys + 1)*(iz2 - iz1 + 1)
      ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)
c
c ----------- get cpus on slab for myid
c
      islab = myid/ncpu_s
      iss   = islab*ncpu_s
      ise   = iss + ncpu_s - 1
      do i=iss,ise
         nsend = (iy_e(i) - iy_s(i) + 1)*ik
         nrecv = (ix_e(i) - ix_s(i) + 1)*jk
         if(i .eq. myid) then
            call send_ytox(g,ft(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(i),iy_e(i),iz1,iz2)
         else
            call send_ytox(g,gt(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(i),iy_e(i),iz1,iz2)
            call mpi_sendrecv(
     +          gt(1),nsend,mpi_real8,i,1,
     +          ft(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_ytox(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)
      enddo
c
      return
      end
      subroutine send_ytox(g,gt,ny,ixs,ixe,iys,iye,izs,ize)
c
c ------------- grab correct chunk of array to be sent
c
      real g(ny,ixs:ixe,izs:ize), gt(iys:iye,ixs:ixe,izs:ize)
c
      do k=izs,ize
      do i=ixs,ixe
      do j=iys,iye
         gt(j,i,k) = g(j,i,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_ytox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
      real f(nx,iys:iye,izs:ize), ft(iys:iye,ixs:ixe,izs:ize)
c
      do k=izs,ize
      do i=ixs,ixe
      do j=iys,iye
         f(i,j,k) = ft(j,i,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine xtoz_trans(f,g,nx,nz,ixs,ixe,ix_s,ix_e,
     +           iys,iye,izs,ize,iz_s,iz_e,
     +           myid,ncpu_s,numprocs)
c
c ------- transpose array  f(nx,iys:iye,izs-1:ize+1) 
c                     ---> g(0:nz+1,iys:iye,ixs:ixe)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,izs-1:ize+1), g(0:nz+1,iys:iye,ixs:ixe)
      real ft(nx*(iye+1-iys)*(ize-izs+1)),
     +     gt(nz*(ixe+1-ixs)*(iye-iys+1))
      integer ix_s(0:numprocs-1), ix_e(0:numprocs-1),
     +        iz_s(0:numprocs-1), iz_e(0:numprocs-1)
c
      jk = (ize - izs + 1)*(iye - iys + 1)
      ij = (ixe - ixs + 1)*(iye - iys + 1)
c
c ----------- get starting location
c
      iss = myid - (myid/ncpu_s)*ncpu_s
c
      do i=iss,numprocs-1,ncpu_s
         nsend = (ix_e(i) - ix_s(i) + 1)*jk
         nrecv = (iz_e(i) - iz_s(i) + 1)*ij
         if(i .eq. myid) then
            call send_xtoz(f,gt(1),nx,ix_s(i),ix_e(i),
     +                  iys,iye,iz_s(myid),iz_e(myid))
         else
            call send_xtoz(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iys,iye,iz_s(myid),iz_e(myid))
            call mpi_sendrecv(
     +          ft(1),nsend,mpi_real8,i,1,
     +          gt(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_xtoz(g,gt(1),nz,ix_s(myid),ix_e(myid),
     +                  iys,iye,iz_s(i),iz_e(i))
      enddo
c
      return
      end
      subroutine send_xtoz(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
c
c ------- grab correct chunk of array to be sent and skip ghost points
c
      real f(nx,iys:iye,izs-1:ize+1), ft(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         ft(i,j,k) = f(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_xtoz(g,gt,nz,ixs,ixe,iys,iye,izs,ize)
      real g(0:nz+1,iys:iye,ixs:ixe), gt(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         g(k,j,i) = gt(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine ztox_trans(g,f,nx,nz,ixs,ixe,ix_s,ix_e,
     +           iys,iye,izs,ize,iz_s,iz_e,
     +           myid,ncpu_s,numprocs)
c
c ------- transpose array g(0:nz+1,iys:iye,ixs:ixe) 
c                    ---> f(nx,iys:iye,izs-1:ize+1)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,izs-1:ize+1), g(0:nz+1,iys:iye,ixs:ixe)
      real ft(nx*(iye+1-iys)*(ize-izs+3)),
     +     gt((nz+3)*(iye+1-iys)*(ixe-ixs+1))
c
      integer ix_s(0:numprocs-1), ix_e(0:numprocs-1),
     +        iz_s(0:numprocs-1), iz_e(0:numprocs-1)
c
      jk = (ize - izs + 3)*(iye - iys + 1)
      ij = (ixe - ixs + 1)*(iye - iys + 1)
c
c ------------- get starting location
c
      iss = myid - (myid/ncpu_s)*ncpu_s
c
      do i=iss,numprocs-1,ncpu_s
         nsend = (iz_e(i) - iz_s(i) + 3)*ij
         nrecv = (ix_e(i) - ix_s(i) + 1)*jk
         if(i .eq. myid) then
            call send_ztox(g,ft(1),nz,ix_s(myid),ix_e(myid),
     +                  iys,iye,iz_s(i),iz_e(i))
         else
            call send_ztox(g,gt(1),nz,ix_s(myid),ix_e(myid),
     +                  iys,iye,iz_s(i),iz_e(i))
            call mpi_sendrecv(
     +          gt(1),nsend,mpi_real8,i,1,
     +          ft(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_ztox(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iys,iye,iz_s(myid),iz_e(myid))
      enddo
c
      return
      end
      subroutine send_ztox(g,gt,nz,ixs,ixe,iys,iye,izs,ize)
c
c ------------- grab correct chunk of array to be sent,
c               account for ghost points
c
      real g(0:nz+1,iys:iye,ixs:ixe), gt(izs-1:ize+1,iys:iye,ixs:ixe)
c
      do j=iys,iye
      do i=ixs,ixe
      do k=izs-1,ize+1
         gt(k,j,i) = g(k,j,i)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_ztox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
      real f(nx,iys:iye,izs-1:ize+1), ft(izs-1:ize+1,iys:iye,ixs:ixe)
c
      do i=ixs,ixe
      do j=iys,iye
      do k=izs-1,ize+1
         f(i,j,k) = ft(k,j,i)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine exchange
c
c ------------- exchange ghost points with mpi,
c               nb and nt are the destination and
c               source nodes. Allows for 1z per cpu
c
      use pars
      use fields
c     use fftwk
      include 'mpif.h'
c
      real fs(nnx,iys:iye,(4+nscl)),fr(nnx,iys:iye,(4+nscl))
      integer istatus(mpi_status_size)
c
      nb = myid - ncpu_s
      nt = myid + ncpu_s
c
c ------------ account for endpoints
c
      if(iss .eq. 0) then
         nb = mpi_proc_null
      endif
      if(ise .eq. numprocs-1) then
         nt = mpi_proc_null
      endif
      nsend = nnx*(iye + 1 - iys)*(4+nscl)
      nrecv = nsend
c
c --------- send top of myid, receive bottom from myid - ncpu_s
c
      do iy=iys,iye
      do ix=1,nnx
         fs(ix,iy,1) = u(ix,iy,ize)
         fs(ix,iy,2) = v(ix,iy,ize)
         fs(ix,iy,3) = w(ix,iy,ize)
         fs(ix,iy,4) = e(ix,iy,ize)
      enddo
      enddo
      do iscl=1,nscl
         jloc = 4 + iscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,jloc) = t(ix,iy,iscl,ize)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nt,0,
     +     fr(1,iys,1),nrecv,mpi_real8,nb,0,
     +     mpi_comm_world,istatus,ierr)
c
      if(iss .ne. 0) then
         izm1 = izs-1
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,izm1) = fr(ix,iy,1)
            v(ix,iy,izm1) = fr(ix,iy,2)
            w(ix,iy,izm1) = fr(ix,iy,3)
            e(ix,iy,izm1) = fr(ix,iy,4)
         enddo
         enddo
         do iscl=1,nscl
            jloc = 4 + iscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izm1) = fr(ix,iy,jloc)
            enddo
            enddo
         enddo
      endif
c
c -------- send bottom of myid, receive bottom from myid + ncpu_s
c
      do iy=iys,iye
      do ix=1,nnx
         fs(ix,iy,1) = u(ix,iy,izs)
         fs(ix,iy,2) = v(ix,iy,izs)
         fs(ix,iy,3) = w(ix,iy,izs)
         fs(ix,iy,4) = e(ix,iy,izs)
      enddo
      enddo
      do iscl=1,nscl
         jloc = 4 + iscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,jloc) = t(ix,iy,iscl,izs)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nb,1,
     +     fr(1,iys,1),nrecv,mpi_real8,nt,1,
     +     mpi_comm_world,istatus,ierr)
c
      if(ise .ne. numprocs-1) then
         izp1 = ize+1
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,izp1) = fr(ix,iy,1)
            v(ix,iy,izp1) = fr(ix,iy,2)
            w(ix,iy,izp1) = fr(ix,iy,3)
            e(ix,iy,izp1) = fr(ix,iy,4)
         enddo
         enddo
         do iscl=1,nscl
            jloc = 4 + iscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izp1) = fr(ix,iy,jloc)
            enddo
            enddo
         enddo
      endif
c
c --------------- send extra scalar points 
c
      nsend = nnx*(iye + 1 - iys)*nscl
      nrecv = nsend
c
c -------------- send top of myid, receive bottom from myid - ncpu_s
c
      izm1 = ize-1
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,iscl) = t(ix,iy,iscl,izm1)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nt,0,
     +     fr(1,iys,1),nrecv,mpi_real8,nb,0,
     +     mpi_comm_world,istatus,ierr)
c
      if(iss .ne. 0) then
         izm2 = izs-2
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izm2) = fr(ix,iy,iscl)
            enddo
            enddo
         enddo
      endif
c
c -------------- send bottom of myid, receive bottom from myid + ncpu_s
c
      izp1 = izs+1
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,iscl) = t(ix,iy,iscl,izp1)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nb,1,
     +     fr(1,iys,1),nrecv,mpi_real8,nt,1,
     +     mpi_comm_world,istatus,ierr)
c
      if(ise .ne. numprocs-1) then
         izp2 = ize+2
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izp2) = fr(ix,iy,iscl)
            enddo
            enddo
         enddo
      endif
c
      return
      end
      subroutine bcast_pbc
c
c ---- send upper boundary conditions to other processors 
c      for fft solution of pressure
c
      use pars
      use fields
      include 'mpif.h'
      integer istatus(mpi_status_size),ierr
c
      if(numprocs .eq. 1) go to 999
c
      irow_r = mod(myid,ncpu_s)
      irow_t = is_s(numprocs-1) + irow_r
      num = nnx*(iye+1-iys)
c
c
      if(iss .ne. is_s(numprocs-1)) then
c
c ------ not in the top row, receive from top
c
        call mpi_recv(pbc(1,iys,1),num,mpi_real8,irow_t,1,
     +       mpi_comm_world,istatus,ierr)
      else
c
c ------ myid is in the top row, send to everyone below
c
        do l=irow_r,irow_t-ncpu_s,ncpu_s
           call mpi_send(pbc(1,iys,1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
c --------- same thing for another variable
c
      if(iss .ne. is_s(numprocs-1)) then
c
c ------ not in the top row, receive from top
c
        call mpi_recv(pbc2(1,iys,1),num,mpi_real8,irow_t,1,
     +       mpi_comm_world,istatus,ierr)
      else
c
c ------ in the top row, send to everyone below
c
        do l=irow_r,irow_t-ncpu_s,ncpu_s
           call mpi_send(pbc2(1,iys,1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
  999 continue
c
      return
      end
      subroutine bcast_lbc
c
c ---- send lower boundary condition information to everyone
c
      use pars
      use fields
      use con_data
      implicit none
      include 'mpif.h'
      integer :: ierr
      real :: upars(4) 


      upars(1) = utau
      upars(2) = u10
      upars(3) = qstar(1)
      upars(4) = qstar(2)

      call mpi_sum_xy(upars,myid,iss,ise,4)

      u10 = upars(2)/ncpu_s
      !qstar(1) = qstar(1)/ncpu_s
      !qstar(2) = qstar(2)/ncpu_s

      call mpi_bcast(u10,1,mpi_real8,
     +              0,mpi_comm_world,ierr)
      call mpi_bcast(qstar,2,mpi_real8,
     +              0,mpi_comm_world,ierr)


      return
      end

      subroutine calc_radiation

      !Calculate a simple radiative temperature source based on Larson
      !et al. (2007) MWR, used in Mellado et al. (2018) JAMES
      use pars
      use particles
      use con_stats
      use con_data
      implicit none
      include 'mpif.h'

      real :: LWP_tmp(nnz-1),LWP_rad(nnz-1)  !Have to calulate the integral from top down, then reverse ordering
      real :: rhoa,func_rho_base
      integer :: iz,lc

      !Do trapezoidal integration from the TOP down:
      !LWC only calculated on INTERIOR zw points

      lc = 1
      do iz=nnz-1,1,-1

         if (lc .eq. 1) then

         if (iexner.eq.1) then
            rhoa = func_rho_base(surf_p,tsfcc(1),zz(iz))
         else
            rhoa = surf_rho
         end if
         LWP_tmp(lc) = 0.5*rhoa*(ql(iz+1) + ql(iz))*(zz(iz+1)-zz(iz))

         else


         LWP_tmp(lc) = LWP_tmp(lc-1) +
     +   0.5*rhoa*(ql(iz+1) + ql(iz))*(zz(iz+1)-zz(iz))

         end if

         lc = lc + 1
      end do

      !Reverse the order so that LWP_rad(1) is the bottom
      !LWP_rad is stored at the interior zw points, since this is center of integral trapezoid
      do iz=1,nnz-1
         LWP_rad(iz) = LWP_tmp(nnz-iz)
      end do

      !Now calculate the radiative flux based on LWP_rad -- at zw points,
      !including top/bottom

      do iz=1,nnz-1

         radflux(iz) = rad_Fct*exp(-rad_kappa*LWP_rad(iz)) +
     +         rad_Fcb*exp(-rad_kappa*(LWP_rad(1) - LWP_rad(iz)))


      end do


      !Give no radiative flux divergence at top/bottom
      radflux(0) = radflux(1)
      radflux(nnz) = radflux(nnz-1)

      !Now finally calculate the "radsrc" tendency term for the temperature equation
      !Based on the divergence of the radflux

      do iz=1,nnz

         if (iexner.eq.1) then
            rhoa = func_rho_base(surf_p,tsfcc(1),zz(iz))
         else
            rhoa = surf_rho
         end if
         radsrc(iz) = -1/rhoa/Cpa*(radflux(iz)-radflux(iz-1))/dzw(iz)

      end do


      end subroutine calc_radiation

      subroutine humidity_control

      use pars
      use fields
      use fftwk
      use con_data
      use con_stats

      implicit none
      include 'mpif.h'

      real :: meanRH_curr,meanRH_target,msqrRH
      integer :: iz,ierr


      !Trapezoidal rule to get current volume mean RH:
      meanRH = 0.0
      do iz=1,nnz-1

         meanRH = meanRH + 
     +   0.5*(zz(iz+1)-zz(iz))*(RHxym(iz)+RHxym(iz+1))

      end do 
    
      meanRH = meanRH/(zz(nnz)-zz(1))  !Note: divide by the distance between the points, not the entire zl, to get the true average

      meanRH_curr = meanRH

      !!!!As an aside, also compute the variance of the RH:
      msqrRH = 0.0
      do iz=1,nnz-1
         msqrRH = msqrRH + 
     +   0.5*(zz(iz+1)-zz(iz))*(RHmsqr(iz)+RHmsqr(iz+1))

      end do 
     
      msqrRH = msqrRH/(zz(nnz)-zz(1))

      varRH = msqrRH - meanRH**2
      !!!! 


      if (ihumiditycontrol) then
          meanRH_target = 103.3
          Swall = Swall + (meanRH_target - meanRH_curr)*1.0e-9
      end if


      end subroutine humidity_control

      subroutine change_RH_bcs_to_q
      !This converts the RH provided in the params file to a value of q that code needs
      use pars
      use con_data
      implicit none
      real :: RHB,RHT,mod_magnus

      !Assuming tsfcc(2),Ttop(2),Tbot(2) in params.in are giving RH:
      if (iDNS .eq. 1) then

         RHT = Ttop(2)
         RHB = Tbot(2)

         !Convert RH given in input file into specific humidity
         Ttop(2) = RHT/100.0*Mw/Ru/Ttop(1)*mod_magnus(Ttop(1))/surf_rho
         Tbot(2) = RHB/100.0*Mw/Ru/Tbot(1)*mod_magnus(Tbot(1))/surf_rho

         if (myid==0) then
           write(*,*) 'SURFACE HUMIDITY: %RH, q = ',RHB,Tbot(2)
         end if

      else !(doing LES)

         RHB = tsfcc(2)

         tsfcc(2)=RHB/100.0*Mw/Ru/tsfcc(1)*mod_magnus(tsfcc(1))/surf_rho

         if (myid==0) then
           write(*,*) 'SURFACE HUMIDITY: %RH, q = ',RHB,tsfcc(2)
         end if

      end if



      end subroutine change_RH_bcs_to_q

      subroutine fill_base(T_surf)
      !Fill the base state profiles p_base,rho_base,T_base,theta_base
      !given the surface conditions T_surf and p_surf.
      !theta_base is assumed equal to T_surf (not some other constant
      !reference)
      use fields
      use pars
      use con_stats
      implicit none

      integer :: iz
      real, intent(in) :: T_surf
      real :: func_T_base,func_p_base,func_rho_base

      if (iexner.eq.1) then
         do iz=0,nnz+1
            theta_base(iz) = T_surf
            T_base(iz) = func_T_base(T_surf,zz(iz))
            p_base(iz) = func_p_base(surf_p,T_surf,zz(iz))
            rho_base(iz)=func_rho_base(surf_p,T_surf,zz(iz))
         end do
      else
         !Set entire arrays to constants
         theta_base = t00
         T_base = t00
         p_base = surf_p
         rho_base = surf_rho
      end if
         

      end subroutine fill_base

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
      function mod_magnus(T)
      implicit none

      !Take in T in Kelvin and return saturation vapor pressure using function of Alduchov and Eskridge, 1996
      real,intent(in) :: T
      real :: mod_magnus

      mod_magnus = 610.94 *exp((17.6257*(T-273.15))/(243.04+(T-273.15)))


      end function mod_magnus
      function exner(p0,p)
      use pars
      implicit none
      real, intent(in) :: p0,p
      real :: exner

      !Take in the reference pressure p0 and p(z), and compute exner = T/theta
      exner = (p0/p)**(-Rd/Cpa)

      end function exner

      !!!!!The dry adiabatic functions:
      function func_p_base(p_surf,T_surf,z)
      use pars, only: grav,Cpa,Rd
      implicit none
      real, intent(in) :: p_surf,T_surf,z
      real :: func_p_base

      func_p_base = p_surf*(1 - z*grav/Cpa/T_surf)**(Cpa/Rd)

      end function func_p_base
      function func_T_base(T_surf,z)
      use pars, only: grav,Cpa
      implicit none
      real, intent(in) :: T_surf,z
      real :: func_T_base
      
      func_T_base = T_surf - grav/Cpa*z

      end function func_T_base
      function func_rho_base(p_surf,T_surf,z)
      use pars, only: grav,Cpa,Rd
      implicit none
      real, intent(in) :: p_surf,T_surf,z
      real :: func_rho_base
      
      func_rho_base = p_surf*T_surf**(-Cpa/Rd)/Rd*
     +         (T_surf - grav*z/Cpa)**(Cpa/Rd-1.0)

      end function func_rho_base
      !!!!!

 
      function ran2(idum)
      integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real :: ran2,AM,EPS,RNMX
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
       real :: fac,gset,rsq,v1,v2,ran2
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
      real :: d,CDF,M,S,cdf_func_single

      cdf_func_single = 1.0/2.0*
     + (1.0 + erf((log(d) - M)/S/sqrt(2.0)))  - CDF

      end function cdf_func_single
      function cdf_func(d,CDF,factor,Ma,Sa,Mc,Sc)
      implicit none
      real :: d,CDF,factor,Ma,Sa,Mc,Sc,cdf_func

      cdf_func = 1.0/(1.0+factor)/2.0*
     + ( 1.0 + erf((log(d) - Ma)/Sa/sqrt(2.0)) +
     + factor*(1 + erf((log(d) - Mc)/Sc/sqrt(2.0))) ) - CDF

      end function cdf_func
