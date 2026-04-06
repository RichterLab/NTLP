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
      subroutine setcon
c
      use pars
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
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
      call rffti(nnx,trigx(:,1))
      call rffti(nny,trigx(:,2))
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
      use mod_solver
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
      use mod_mpi
      use mod_fft
      use mod_solver
      implicit none
 
      real :: psi(nnx, iys:iye), psix(nnx, iys:iye),
     +     psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +     vyy(nnx,iys:iye,izs:izs),Ttmp

      real :: ampv, ampt, sum_psi, vmaxx, vmag, facv

      real :: slope_q, slope_t, slope_u, slope_v
      real :: slope_t1, slope_t2, slope_u1, slope_v1
      real :: z_switch, t_switch, q_switch, u_switch, v_switch
      real :: z_switch_q, z_switch_t, z_switch_v, z_switch_u
      real :: t_surf, q_surf, u_surf, v_surf
      real :: t_temp, q_temp, u_temp, v_temp
      real :: q_top, t_top, u_top, v_top, slope_q1, slope_q2

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
         call xderivp(psix(1,iys),trigx(:,1),xk,nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(:,2),yk,
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
         call xderivp(uxx(1,iys),trigx(:,1),xk,nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(:,2),yk,
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
      call mpi_sum_z(divz,i_root,myid,nnz,1)
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
      use mod_mpi
      use mod_fft
      use mod_solver
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
         call xderivp(psix(1,iys),trigx(:,1),xk,nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(:,2),yk,
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
         call xderivp(uxx(1,iys),trigx(:,1),xk,nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(:,2),yk,
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
      call mpi_sum_z(divz,i_root,myid,nnz,1)
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
      use mod_mpi
      use mod_fft
      use mod_solver
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
         call xderivp(psix(1,iys),trigx(:,1),xk,nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(:,2),yk,
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
         call xderivp(uxx(1,iys),trigx(:,1),xk,nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(:,2),yk,
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
      call mpi_sum_z(divz,i_root,myid,nnz,1)
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
      use mod_mpi
      use mod_thermo
      use mod_solver
      implicit none
c
c ------- indices for indexing array stat(.,.)
c         js = number of non-scalar stats
c         ns = number of scalar stats
c
      integer, parameter :: js = 19, ns = 5, nstat = js + ns*nscl
      real :: stat(1:nnz, nstat)

      integer :: ix,iy,iz,izm1,izp1,izp2,m1,m2,m3,m4,m5,i,l

      real :: sgn
      real :: RHtmp, Ttmp
      real :: rhoa, msqrRH
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
     +      exner(surf_p,func_p_base(surf_p,tsfcc(1),zz(iz)))
            !rhoa = func_rho_base(surf_p,tsfcc(1),zz(iz))
	    rhoa = func_p_base(surf_p,tsfcc(1),zz(iz))/Rd/Ttmp
         else
            Ttmp = t(ix,iy,1,iz)
            rhoa = surf_rho
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
      call mpi_sum_z(pflux,i_root,myid,nnz+1,1)
      call mpi_sum_z(pfluxdiff,i_root,myid,nnz+1,1)
      call mpi_sum_z(pmassflux,i_root,myid,nnz+1,1)
      call mpi_sum_z(penegflux,i_root,myid,nnz+1,1)

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
      use mod_mpi
c
      real :: stat(1:nnz, 5)
      real :: s11s, s22s, s33s, s12s, s13s, s23s, wz, wzp
      real :: s13, s23, s33
      real :: ufluc, ufluc_t, ufluc_b, vfluc, vfluc_t, vfluc_b, wfluc
      real :: ffluc1, ffluc2, ffluc3
      real :: ffluc1p, ffluc2p, ffluc3p
      real :: weit, weit1
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
      subroutine Tvar_budget(iscl)
      use pars
      use particles
      use fields
      use con_data
      use con_stats
      use fftwk
      use mod_mpi
      use mod_fft
      implicit none

      real :: stat(1:nnz, 4)
      real :: weit, weit1
      real :: tx_tmp(nnx, iys:iye), ty(nnx, iys:iye, izs-1:ize+1)
      real :: tx(nnx, iys:iye, izs-1:ize+1)
      real :: trans(izs:ize)
      real :: gradTp(3), Tfluc, dTpdz1, dTpdz, dTdz
      real :: Tflucp1, Tflucm1, qfluc, qflucp1, Tmean
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
      call xderivp(tx_tmp(1,iys),trigx(:,1),xk,nnx,iys,iye)
      tx(1:nnx,iys:iye,iz) = tx_tmp(1:nnx,iys:iye)
      enddo

      call yd_mpi(ty(1,iys,izs-1),trigx(:,2),yk,
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
      call mpi_sum_z(Tv_part2,i_root,myid,nnz,1)
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
      use mod_mpi

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
      implicit none
      include 'mpif.h'
      character cgrid*10, num*3
      integer iblnk
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
