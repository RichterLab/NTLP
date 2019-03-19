      subroutine particle_update_rk3(it,istage)
      use pars
      use particles
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'

      integer :: istage,ierr,it
      real :: g(3)
      real :: denom,dtl,sigma,ttest
      integer :: ix,iy,iz,ind
      integer :: apr_count
      real :: Rep,diff(3),diffnorm,corrfac,myRep_avg
      real :: mydpmean, mydpmin, mydpmax
      real :: xtmp(3),vtmp(3),Tptmp,radiustmp
      real :: Nup,Shp,rhop,taup_i,estar,einf
      real :: myphip_sum,myphiw_sum,myphiv_sum,Volp      
      real :: TfC,TpC,Lv,Eff_C,Eff_S
      real :: t_s,t_f,t_s1,t_f1

      g(1:3) = (/0.0, 0.0, part_grav/)

      !First fill extended velocity field for interpolation
      !t_s = mpi_wtime()
      call fill_ext 
      !t_f = mpi_wtime()
      !call mpi_barrier(mpi_comm_world,ierr)
      !if (myid==5) write(*,*) 'time fill_ext:',t_f-t_s

      partcount_t     = 0.0
      partflux_t      = 0.0
      partcountU_t    = 0.0
      pcountrange_t   = 0.0
      vpsum_t         = 0.0
      upwp_t          = 0.0
      vpsqrsum_t      = 0.0
      Tpsum_t         = 0.0
      Tfsum_t         = 0.0
      qfsum_t         = 0.0
      radsum_t        = 0.0  
      rad2sum_t       = 0.0  
      mpsum_t         = 0.0
      mwsum_t         = 0.0
      Tpsqrsum_t      = 0.0
      wpTpsum_t       = 0.0
      myRep_avg       = 0.0
      myphip_sum      = 0.0
      myphiw_sum      = 0.0
      myphiv_sum      = 0.0
      qstarsum_t      = 0.0 
      mydpmean        = 0.0
      mydpmin         = 2.0*radius_max
      mydpmax         = 2.0*radius_min
      apr_count_t     = 0

      !t_s = mpi_wtime()

      
      !Loop over the linked list of particles:
      part => first_particle
      do while (associated(part))

      if(part%type_part .eq. 0)then !Aerosol or aircraft? 
        !First, interpolate to get the fluid velocity part%uf(1:3):
        call uf_interp
        if(it .LE. 1 )then 
          !part%xrhs(1:3) = part%vp(1:3)
          !part%xp(1:3) = xtmp(1:3) + dt*gama(istage)*part%xrhs(1:3)
          part%vp(1:3) = part%uf
          part%Tp      = part%Tf
        endif

        !Now advance the particle and position via RK3 (same as velocity)
         
        ! Calculate equilibrium radius
        call particle_eq_radius        ! partical equ. radius
 
        !Intermediate Values
        diff(1:3)  = part%vp - part%uf
        diffnorm   = sqrt(diff(1)**2 + diff(2)**2 + diff(3)**2)
        Rep        = 2.0*part%radius*diffnorm/nuf  
        Volp       = pi2*2.0/3.0*part%radius**3
        rhop       = (m_s+Volp*rhow)/Volp
        taup_i     = 18.0*rhoa*nuf/rhop/(2.0*part%radius)**2


c       Correction here
c       if(myid==0)then
c       write(*,*) 'Tau_p = ', taup_i
c       end if  

        mydpmean = mydpmean + 2.0*part%radius
        if(2.0*part%radius .lt. mydpmin)then
          mydpmin = 2.0*part%radius
        endif

        if(2.0*part%radius .gt. mydpmax)then
          mydpmax = 2.0*part%radius
        endif 


        myRep_avg  = myRep_avg + Rep
        corrfac    = (1.0 + 0.15*Rep**(0.687))
        myphip_sum = myphip_sum + Volp*rhop
        myphiw_sum = myphiw_sum + Volp*rhow
        myphiv_sum = myphiv_sum + Volp

        !Compute Nusselt number for particle:
        !Ranz-Marshall relation
        Nup = 2.0 + 0.6*Rep**(1.0/2.0)*Pra**(1.0/3.0)
        Shp = 2.0 + 0.6*Rep**(1.0/2.0)*Sc**(1.0/3.0)

        !Aircraft measurements based on particle location
        !exact position in the code is variable..
        if(istage .eq. 3)then
          if(sqrt((part%xp(1)-apr_loc(1))**2+(part%xp(2)-apr_loc(2))**2
     +             +(part%xp(3)-apr_loc(3))**2).lt.apr_rad)then
            apr_count = apr_count+1
          end if
        end if



        !Mass Transfer calculations
        TfC   = part%Tf-273.15
        einf  = 610.94*EXP(17.6257*TfC/(TfC+243.04))
        TpC   = part%Tp-273.15
        Lv    = (25.0 - 0.02274*26.0)*10.0**5

        Eff_C = 2.0*Mw*Gam/(Ru*rhow*part%radius*part%Tp)
        Eff_S = Ion*Os*m_s*Mw/Ms/(Volp*rhop-m_s)
        estar = einf*exp(Mw*Lv/Ru*(1.0/part%Tf-1.0/part%Tp)
     +                                        +Eff_C-Eff_S)
        part%qstar = Mw/Ru*estar/part%Tp/rhoa

  
        xtmp(1:3)  = part%xp(1:3) + dt*zetas(istage)*part%xrhs(1:3)
        vtmp(1:3)  = part%vp(1:3) + dt*zetas(istage)*part%vrhs(1:3) 
        !Tptmp = part%Tp + dt*zetas(istage)*part%Tprhs_s
        !Tptmp = Tptmp + dt*zetas(istage)*part%Tprhs_L
        radiustmp  = part%radius + dt*zetas(istage)*part%radrhs

        part%xrhs(1:3) = part%vp(1:3)
        part%vrhs(1:3) = corrfac*taup_i*(part%uf(1:3)-part%vp(1:3))
     +                               - g(1:3)

        if(ievap .eq. 1)then      
          part%radrhs = Shp/9.0/Sc*rhop/rhow*part%radius*
     +                  taup_i*(part%qinf-part%qstar) !assumes qinf=rhov/rhoa rather than rhov/rhom
        else
          part%radrhs = 0.0
        end if 
         ! Restrict the particle diameter when tw<dt    
        if(ievap .eq. 1 .and. part%radius .le. 20.0e-6)then
          part%radrhs = 0.0
          call particle_eq_radius
          part%radius = part%r_eq
          !else
          !part%radius = radiustmp + dt*gama(istage)*part%radrhs
        end if
 
        if(istage==1)then
          part%xp_i(1:3)  = part%xp(1:3)  !RK3 flux reference
          part%xp_ti(1:3) = part%xp(1:3) !Total flux reference
        end if
        
         !part%Tprhs_s = -Nup/3.0/Pra*CpaCpp*
!     +                     rhop/rhow*taup_i*(part%Tp-part%Tf)
        ! part%Tprhs_L = 3.0*Lv/Cpp/part%radius*part%radrhs
  
        part%xp(1:3) = xtmp(1:3) + dt*gama(istage)*part%xrhs(1:3)
        do ind = 1,2
          part%tot_dist = part%tot_dist + xtmp(ind)+
     +                    dt*gama(istage)*part%xrhs(ind)
        end do
        part%vp(1:3)    = vtmp(1:3) + dt*gama(istage)*part%vrhs(1:3)
        !part%Tp         = Tptmp + dt*gama(istage)*part%Tprhs_s
        !part%Tp         = part%Tp + dt*gama(istage)*part%Tprhs_L
        part%radius     = radiustmp + dt*gama(istage)*part%radrhs

        !--------Velocity of the particle not resolved
        part%vp(1:3)    = part%uf(1:3) - g(1:3)/taup_i
c       part%vp(1:3)     = part%uf(1:3)    ! make particles tracers


        ! part%Tp = part%Tf
        !-------
c        Calculate Equilibrium temperature for droplet radius evolution ! Indrajith 02/07/2017
c       ---------------
c        call particle_eq_Temp 
c        part%Tp = part%T_eq
c       --------------
c        check whehter the particle radius go below r_eq
c       -------------REMOVE FOR CENCENTRATION PROFILE CALCULATIONS
c        if(part%radius .lt.part%r_eq)then
c            part%radius = part%r_eq
c        end if

      else !Aircraft particles!
        if(istage==3)then
          call uf_interp
          part%xp(1:3) = part%xp(1:3) + dt*part%vp(1:3)
!          write(*,*)'current proc:',myid
!          write(*,*)'xmin xmax: ',xmin,xmax
          apr_loc(1:3) = part%xp(1:3)
        end if



      end if !For aircraft probes
        part => part%next
      end do
      !t_f1 = mpi_wtime()
      !write(*,*) 'proc,loop time: ',myid,t_f1-t_s
      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time loop:', t_f-t_s


      call mpi_bcast(apr_loc,3,mpi_real8,
     +               apr_id,mpi_comm_world,ierr)


c     JP 5-15-2018 Initialize xp_i after RK3 for subgrid fluxes
      part => first_particle
      do while (associated(part))
        if(part%type_part==0)then
          part%xp_si(1:3) = part%xp(1:3)
        end if
          part => part%next
      end do

      !Enforce nonperiodic bcs (either elastic or destroying particles)
      call particle_bcs_nonperiodic

      !Check to see if particles left processor
      !If they did, remove from one list and add to another
      !t_s = mpi_wtime()
      call particle_exchange
      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time exchg:', t_f - t_s

      !Now enforce periodic bcs 
      !just updates x,y locations if over xl,yl or under 0
      call particle_bcs_periodic

      !Now that particles are in their updated position, 
      !compute their contribution to the momentum coupling:
      call particle_coupling_update
      call particle_rk3_Uflux
 
      !Finally, now that coupling and statistics arrays are filled, 
      !Transpose them back to align with the velocities:
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,1),
     +                partsrc(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,2),
     +                partsrc(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,3),
     +                partsrc(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(partTsrc_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partTsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(partHsrc_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partHsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(partTEsrc_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partTEsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)


      !Try only calling these when the history data is being written:
      if(mtrans) then
      call ztox_trans(upwp_t(0:nnz+1,iys:iye,mxs:mxe),
     +               upwp(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,1),
     +               vpsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,2),
     +               vpsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,3),
     +               vpsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,1),
     +               vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,2),
     +               vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,3),
     +               vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(Tpsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +               Tpsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(Tpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +               Tpsqrsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(Tfsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +               Tfsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(qfsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +               qfsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(wpTpsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +               wpTpsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(partcount_t(0:nnz+1,iys:iye,mxs:mxe),
     +               partcount(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(radsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +               radsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(rad2sum_t(0:nnz+1,iys:iye,mxs:mxe),
     +               rad2sum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs) 

      call ztox_trans(mpsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +               mpsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(mwsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +               mwsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(qstarsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +               qstarsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(partflux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               partflux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(partcountU_t(0:nnz+1,iys:iye,mxs:mxe),
     +               partcountU(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pcountrange_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pcountrange(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +               mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)



      end if


      !t_s = mpi_wtime
      !Get particle count:
      numpart = 0
      part => first_particle
      do while (associated(part))
        numpart = numpart + 1
        part => part%next
      end do
      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time numpart: ', t_f - t_s
 
      !t_s = mpi_wtime()
      !Compute total number of particles
      call mpi_allreduce(numpart,tnumpart,1,mpi_integer,mpi_sum,
     +                   mpi_comm_world,ierr)
      !Compute average particle Reynolds number
      call mpi_allreduce(myRep_avg,Rep_avg,1,mpi_real8,mpi_sum,
     +                   mpi_comm_world,ierr)

      Rep_avg = Rep_avg/tnumpart

      call mpi_allreduce(mydpmean,dpmean,1,mpi_real8,mpi_sum,
     +                   mpi_comm_world,ierr)

      dpmean = dpmean/tnumpart

      call mpi_allreduce(mydpmax,dpmax,1,mpi_real8,mpi_max,
     +                   mpi_comm_world,ierr)
      call mpi_allreduce(mydpmin,dpmin,1,mpi_real8,mpi_min,
     +                   mpi_comm_world,ierr)

      call mpi_allreduce(myphip_sum,phip,1,mpi_real8,mpi_sum,
     +                   mpi_comm_world,ierr)

      call mpi_allreduce(myphiw_sum,phiw,1,mpi_real8,mpi_sum,
     +                   mpi_comm_world,ierr)

      call mpi_allreduce(myphiv_sum,phiv,1,mpi_real8,mpi_sum,
     +                   mpi_comm_world,ierr)

      call mpi_allreduce(apr_count,apr_count_t,1,mpi_integer,
     +                   mpi_sum,mpi_comm_world,ierr)
 
      phip = phip/xl/yl/zl/rhoa
      phiw = phiw/xl/yl/zl/rhoa
      phiv = phiv/xl/yl/zl

      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time mpi_allreduce: ', t_f - t_s

      end subroutine particle_update_rk3

      subroutine create_particle(xp,vp,Tp,qinfp,idx,rp_init,Pid,time,
     +                           type_p ) !Currently ignoring Tp and qinfp
      use particles
      use pars
      implicit none

      real :: xp(3),vp(3),Tp,qinfp  !Brian 7/29/14
      real :: rp_init,time
      integer :: idx,Pid,type_p

      if (.NOT. associated(first_particle)) then
        allocate(first_particle)
        part => first_particle
        nullify(part%next,part%prev)
      else
        !Add to beginning of list since it's more convenient
        part => first_particle
        allocate(part%prev)
        first_particle => part%prev
        part%prev%next => part
        part => first_particle
        nullify(part%prev)
      end if
  
      part%xp(1:3)    = xp(1:3)
      part%xp_i(1:3)  = xp(1:3)
      part%xp_ti(1:3) = xp(1:3)
      part%vp(1:3)    = vp(1:3)
      part%u_sub(1:3) = 0.0     !Indrajith 07/17 SFS vel
      part%sigm_s     = 0.0     !Indrajith 07/17 isotropic stress
      part%Tp         = Tp
      part%radius     = rp_init 
      part%uf(1:3)    = 0.0
      part%qinf       = qinfp  !Brian 7/14/14
      part%xrhs(1:3)  = 0.0
      part%vrhs(1:3)  = 0.0 
      part%Tprhs_s    = 0.0
      part%Tprhs_L    = 0.0
      part%radrhs     = 0.0
      part%pidx       = idx 
      part%procidx    = Pid
      part%birth_time = time 
      part%max_speed  = 0.0
      part%tot_dist   = 0.0
      part%z_highest  = 0.0
      part%type_part  = type_p !whether aerosol, aircraft, or w/e
      
      end subroutine create_particle

      function dFdr_uniform(r)
      implicit none
      real :: r,dFdr_uniform 

      dFdr_uniform = 0.5e-6

      end function dFdr_uniform

      subroutine particle_bcs_nonperiodic
      use particles
      use pars
      use con_data
      implicit none
      include 'mpif.h'  
      real :: top,bot
      integer :: values(8)
      integer :: bin_location
c      integer :: Npdx

      !Assumes domain goes from [0,xl),[0,yl),[0,zl]
      !Also maintain the number of particles on each proc
c      Npdx = 0 
      part => first_particle
      do while (associated(part))

      !perfectly elastic collisions on top, bottom walls
      !i.e. location is reflected, w-velocity is negated

c        top = zl  - delta_vis
c        bot = 0.0 + delta_vis

        !For particle radius bounce-back
        top = zl  - part%radius
        bot = 0.0 + part%radius

        if(part%xp(3) .GT. top)then
c          if (part%xp(3) .GT. zi) then
c            part%xp(3) = zi - (part%xp(3)-zi)     
          part%xp(3) = top - (part%xp(3)-top)
          part%vp(3) = -part%vp(3)
          part => part%next
        elseif(part%xp(3) .LT. bot)then
          if(inonperiodic .eq. 1)then
            if(ipolydisp .eq. 1)then
c             figure out which particle size this range is under.
              bin_location = minloc(class_maxr,DIM=1,
     +                              MASK=(part%radius<class_maxr))
              del_radclass(bin_location) = 
     +        del_radclass(bin_location) + 1
            end if
            part%death_time=time
            call get_units
            write(nparts+myid) part%birth_time, part%death_time,
     +          part%z_highest,part%tot_dist,part%max_speed,
     +          part%pidx 

            call destroy_particle 
            died = died + 1

!       This reinject is in order to take account into particle loss
!        reinject = reinject+1
          else
         !periodic vertical boundary condition
            part%xp(3) = bot+zwpi/2 !Half of plane flux
            part%vp(3) = -part%vp(3)
            if(isfs==1)part%u_sub(3) = -part%u_sub(3)
            part => part%next
          end if
        else
          part => part%next
        end if
      end do
      end subroutine particle_bcs_nonperiodic

      subroutine particle_bcs_periodic
      use particles
      use pars
      implicit none 

      !Assumes domain goes from [0,xl),[0,yl),[0,zl] 
      !Also maintain the number of particles on each proc
      
      part => first_particle
      do while (associated(part))

      !x,y periodic
   
      if (part%xp(1) .GT. xl) then
         part%xp(1) = part%xp(1)-xl
      elseif (part%xp(1) .LT. 0) then
         part%xp(1) = xl + part%xp(1)
      end if

      if (part%xp(2) .GT. yl) then
         part%xp(2) = part%xp(2)-yl
      elseif (part%xp(2) .LT. 0) then
         part%xp(2) = yl + part%xp(2)
      end if

      part => part%next

      end do
      end subroutine particle_bcs_periodic
      
      subroutine destroy_particle
      use particles
      implicit none

      type(particle), pointer :: tmp

      !Is it the first and last in the list?
      if (associated(part,first_particle) .AND. 
     +    (.NOT. associated(part%next)) ) then
          nullify(first_particle)
          deallocate(part)
      else
        if (associated(part,first_particle)) then !Is it the first particle?
           first_particle => part%next
           part => first_particle
           deallocate(part%prev)
        elseif (.NOT. associated(part%next)) then !Is it the last particle?
           nullify(part%prev%next)
           deallocate(part)
        else
           tmp => part
           part => part%next
           tmp%prev%next => tmp%next
           tmp%next%prev => tmp%prev
           deallocate(tmp)
        end if
      end if
   
      end subroutine destroy_particle
      subroutine particle_exchange
      use pars
      use particles
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'

      type(particle), pointer :: tmp
      integer :: idx,psum,csum
      integer :: ir,itr,itop,itl,il,ibl,ib,ibr
      integer :: istatus(mpi_status_size),ierr
      integer :: status_array(mpi_status_size,16),req(16)
      type(particle), allocatable :: rbuf_s(:),trbuf_s(:)
      type(particle), allocatable :: tbuf_s(:),tlbuf_s(:)
      type(particle), allocatable :: lbuf_s(:),blbuf_s(:)
      type(particle), allocatable :: bbuf_s(:),brbuf_s(:)
      type(particle), allocatable :: rbuf_r(:),trbuf_r(:)
      type(particle), allocatable :: tbuf_r(:),tlbuf_r(:)
      type(particle), allocatable :: lbuf_r(:),blbuf_r(:)
      type(particle), allocatable :: bbuf_r(:),brbuf_r(:)
      type(particle), allocatable :: totalbuf(:)
      integer :: aircraft_moved1,apr_id_old
        

c     Currently only for one aircraft, if multiple aircraft this will
c     not work!!!
      aircraft_moved1 = 0
      aircraft_moved2 = 0
      apr_id_old = apr_id
      !Zero out the counters for how many particles to send each dir.
      pr_s=0;ptr_s=0;pt_s=0;ptl_s=0;pl_s=0;pbl_s=0;pb_s=0;pbr_s=0
      
      !As soon as the location is updated, must check to see if it left the proc:
      !May be a better way of doing this, but it seems most reasonable:
      part => first_particle
      do while (associated(part))     




         !First get numbers being sent to all sides:
         if (part%xp(2) .GT. ymax) then 
           if (part%xp(1) .GT. xmax) then !top right (1)
             ptr_s = ptr_s + 1
               if(part%type_part == 1)then
                 aircraft_moved1=trproc 
               end if
           elseif (part%xp(1) .LT. xmin) then !bottom right (3)
             pbr_s = pbr_s + 1
               if(part%type_part == 1)then
                 aircraft_moved1=brproc
               end if
           else  !right (2)
             pr_s = pr_s + 1
               if(part%type_part == 1)then
                 aircraft_moved1=rproc
               end if
         end if

c_________________________
c            write(*,*)'ptr_s pbr_s pr_s'
c_________________________

         elseif (part%xp(2) .LT. ymin) then
           if (part%xp(1) .GT. xmax) then !top left (7)
             ptl_s = ptl_s + 1
               if(part%type_part == 1)then
                 aircraft_moved1=tlproc
               end if
           else if (part%xp(1) .LT. xmin) then !bottom left (5)
             pbl_s = pbl_s + 1
               if(part%type_part == 1)then
                 aircraft_moved1=blproc
               end if
           else  !left (6)
             pl_s = pl_s + 1
               if(part%type_part == 1)then
                 aircraft_moved1=lproc
               end if
           end if

c_________________________
c            write(*,*)'ptl_s pbl_s pl_s'
c_________________________

         elseif ( (part%xp(1) .GT. xmax) .AND.
     +           (part%xp(2) .LT. ymax) .AND.
     +           (part%xp(2) .GT. ymin) ) then !top (8)
            pt_s = pt_s + 1
               if(part%type_part == 1)then
                 aircraft_moved1=tproc
               end if
c_________________________
c            write(*,*)'pt_s'
c_________________________

         elseif ( (part%xp(1) .LT. xmin) .AND.
     +           (part%xp(2) .LT. ymax) .AND.
     +           (part%xp(2) .GT. ymin) ) then !bottom (4)
            pb_s = pb_s + 1
               if(part%type_part == 1)then
                 aircraft_moved1=bproc
               end if

c_________________________
c            write(*,*)'pb_s'
c_________________________

         end if
         
         part => part%next
      end do
      
      !Now allocate the send buffers based on these counts:
      allocate(rbuf_s(pr_s),trbuf_s(ptr_s),tbuf_s(pt_s),tlbuf_s(ptl_s))
      allocate(lbuf_s(pl_s),blbuf_s(pbl_s),bbuf_s(pb_s),brbuf_s(pbr_s))

      !Now loop back through the particles and fill the buffers:
      !NOTE: If it finds one, add it to buffer and REMOVE from list
      ir=1;itr=1;itop=1;itl=1;il=1;ibl=1;ib=1;ibr=1

      part => first_particle
      do while (associated(part))
         
         if (part%xp(2) .GT. ymax) then 
            if (part%xp(1) .GT. xmax) then !top right
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',trproc,'TOP RIGHT'
               trbuf_s(itr) = part
               call destroy_particle
               itr = itr + 1 
            elseif (part%xp(1) .LT. xmin) then !bottom right
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',brproc,'BOTTOM RIGHT'
               brbuf_s(ibr) = part
               call destroy_particle
               ibr = ibr + 1
            else   !right
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',rproc,'RIGHT'
               rbuf_s(ir) = part
               call destroy_particle
               ir = ir + 1
            end if
         elseif (part%xp(2) .LT. ymin) then
            if (part%xp(1) .GT. xmax) then !top left
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',tlproc,'TOP LEFT'
               tlbuf_s(itl) = part
               call destroy_particle
               itl = itl + 1
            else if (part%xp(1) .LT. xmin) then !bottom left
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',blproc,'BOTTOM LEFT'
               blbuf_s(ibl) = part
               call destroy_particle
               ibl = ibl + 1
            else  !left
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',lproc,'LEFT'
               lbuf_s(il) = part
               call destroy_particle
               il = il + 1
            end if
         elseif ( (part%xp(1) .GT. xmax) .AND.
     +           (part%xp(2) .LT. ymax) .AND.
     +           (part%xp(2) .GT. ymin) ) then !top
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',tproc,'TOP'
            tbuf_s(itop) = part
            call destroy_particle
            itop = itop + 1
         elseif ( (part%xp(1) .LT. xmin) .AND.
     +           (part%xp(2) .LT. ymax) .AND.
     +           (part%xp(2) .GT. ymin) ) then !bottom
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',bproc,'BOTTOM'
            bbuf_s(ib) = part
            call destroy_particle
            ib = ib + 1 
         else
         part => part%next
         end if 
         
      end do

      !Now everyone exchanges the counts with all neighbors:
      !Left/right:
      call MPI_Sendrecv(pr_s,1,mpi_integer,rproc,3,
     +        pl_r,1,mpi_integer,lproc,3,mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pl_s,1,mpi_integer,lproc,4,
     +        pr_r,1,mpi_integer,rproc,4,mpi_comm_world,istatus,ierr)

      !Top/bottom:
      call MPI_Sendrecv(pt_s,1,mpi_integer,tproc,5,
     +        pb_r,1,mpi_integer,bproc,5,mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pb_s,1,mpi_integer,bproc,6,
     +        pt_r,1,mpi_integer,tproc,6,mpi_comm_world,istatus,ierr)

      !Top right/bottom left:
      call MPI_Sendrecv(ptr_s,1,mpi_integer,trproc,7,
     +        pbl_r,1,mpi_integer,blproc,7,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pbl_s,1,mpi_integer,blproc,8,
     +        ptr_r,1,mpi_integer,trproc,8,
     +        mpi_comm_world,istatus,ierr)

       !Top left/bottom right:
      call MPI_Sendrecv(ptl_s,1,mpi_integer,tlproc,9,
     +        pbr_r,1,mpi_integer,brproc,9,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pbr_s,1,mpi_integer,brproc,10,
     +         ptl_r,1,mpi_integer,tlproc,10,
     +         mpi_comm_world,istatus,ierr)

      !Now everyone has the number of particles arriving from every neighbor
      !If the count is greater than zero, exchange:

      !Allocate room to receive from each side
      allocate(rbuf_r(pr_r),trbuf_r(ptr_r),tbuf_r(pt_r),tlbuf_r(ptl_r))
      allocate(lbuf_r(pl_r),blbuf_r(pbl_r),bbuf_r(pb_r),brbuf_r(pbr_r))
     
      !Send to right:
      if (pr_s .GT. 0) then
      call mpi_isend(rbuf_s,pr_s,particletype,rproc,11,
     +               mpi_comm_world,req(1),ierr)
      else
      req(1) = mpi_request_null
      end if

      !Receive from left:
      if (pl_r .GT. 0) then
      call mpi_irecv(lbuf_r,pl_r,particletype,lproc,11,
     +               mpi_comm_world,req(2),ierr)
      else
      req(2) = mpi_request_null
      end if

      !Send to left:
      if (pl_s .GT. 0) then
      call mpi_isend(lbuf_s,pl_s,particletype,lproc,12,
     +               mpi_comm_world,req(3),ierr)
      else
      req(3) = mpi_request_null
      end if

      !Receive from right:
      if (pr_r .GT. 0) then
      call mpi_irecv(rbuf_r,pr_r,particletype,rproc,12,
     +               mpi_comm_world,req(4),ierr)
      else
      req(4) = mpi_request_null
      end if

      !Send to top:
      if (pt_s .GT. 0) then
!      write(*,*)'particle information sent from proc: ',myid
      call mpi_isend(tbuf_s,pt_s,particletype,tproc,13,
     +                mpi_comm_world,req(5),ierr)
      else
      req(5) = mpi_request_null
      end if
      
      !Receive from bottom:
      if (pb_r .GT. 0) then
!      write(*,*)'particle information received to proc: ',myid
      call mpi_irecv(bbuf_r,pb_r,particletype,bproc,13,
     +                mpi_comm_world,req(6),ierr)
      else
      req(6) = mpi_request_null
      end if

      !Send to bottom:
      if (pb_s .GT. 0) then
      call mpi_isend(bbuf_s,pb_s,particletype,bproc,14,
     +                mpi_comm_world,req(7),ierr)
      else
      req(7) = mpi_request_null
      end if
      
      !Recieve from top:
      if (pt_r .GT. 0) then
      call mpi_irecv(tbuf_r,pt_r,particletype,tproc,14,
     +                mpi_comm_world,req(8),ierr)
      else
      req(8) = mpi_request_null
      end if

      !Send to top right:
      if (ptr_s .GT. 0) then
      call mpi_isend(trbuf_s,ptr_s,particletype,trproc,15,
     +                mpi_comm_world,req(9),ierr)
      else
      req(9) = mpi_request_null
      end if
     
      !Receive from bottom left:
      if (pbl_r .GT. 0) then
      call mpi_irecv(blbuf_r,pbl_r,particletype,blproc,15,
     +                mpi_comm_world,req(10),ierr)
      else 
      req(10) = mpi_request_null
      end if
    
      !Send to bottom left:
      if (pbl_s .GT. 0) then
      call mpi_isend(blbuf_s,pbl_s,particletype,blproc,16,
     +                mpi_comm_world,req(11),ierr)
      else
      req(11) = mpi_request_null
      end if
     
      !Receive from top right:
      if (ptr_r .GT. 0) then
      call mpi_irecv(trbuf_r,ptr_r,particletype,trproc,16,
     +                mpi_comm_world,req(12),ierr)
      else 
      req(12) = mpi_request_null
      end if

      !Send to top left:
      if (ptl_s .GT. 0) then
      call mpi_isend(tlbuf_s,ptl_s,particletype,tlproc,17,
     +                mpi_comm_world,req(13),ierr)
      else 
      req(13) = mpi_request_null
      end if
    
      !Receive from bottom right:
      if (pbr_r .GT. 0) then
      call mpi_irecv(brbuf_r,pbr_r,particletype,brproc,17,
     +                mpi_comm_world,req(14),ierr)
      else 
      req(14) = mpi_request_null
      end if
  
      !Send to bottom right:
      if (pbr_s .GT. 0) then
      call mpi_isend(brbuf_s,pbr_s,particletype,brproc,18,
     +                mpi_comm_world,req(15),ierr)
      else
      req(15) = mpi_request_null
      end if
  
      !Receive from top left:
      if (ptl_r .GT. 0) then
      call mpi_irecv(tlbuf_r,ptl_r,particletype,tlproc,18,
     +                mpi_comm_world,req(16),ierr)
      else
      req(16) = mpi_request_null
      end if

      call mpi_waitall(16,req,status_array,ierr)

      !Now add incoming particles to linked list:
      !NOTE: add them to beginning since it's easiest to access (first_particle)

      !Form one large buffer to loop through and add:
      psum = pr_r+ptr_r+pt_r+ptl_r+pl_r+pbl_r+pb_r+pbr_r
      csum = 0
      allocate(totalbuf(psum))
      if (pr_r .GT. 0) then 
         totalbuf(1:pr_r) = rbuf_r(1:pr_r)
         csum = csum + pr_r 
      end if
      if (ptr_r .GT. 0) then 
         totalbuf(csum+1:csum+ptr_r) = trbuf_r(1:ptr_r)
         csum = csum + ptr_r
      end if
      if (pt_r .GT. 0) then 
         totalbuf(csum+1:csum+pt_r) = tbuf_r(1:pt_r)
         csum = csum + pt_r
!         write(*,*)'pt_r, myid: ',pt_r,myid
      end if
      if (ptl_r .GT. 0) then 
         totalbuf(csum+1:csum+ptl_r) = tlbuf_r(1:ptl_r)
         csum = csum + ptl_r
      end if
      if (pl_r .GT. 0) then 
         totalbuf(csum+1:csum+pl_r) = lbuf_r(1:pl_r)
         csum = csum + pl_r
      end if
      if (pbl_r .GT. 0) then 
         totalbuf(csum+1:csum+pbl_r) = blbuf_r(1:pbl_r)
         csum = csum + pbl_r
      end if
      if (pb_r .GT. 0) then 
         totalbuf(csum+1:csum+pb_r) = bbuf_r(1:pb_r)
         csum = csum + pb_r
      end if
      if (pbr_r .GT. 0) then 
         totalbuf(csum+1:csum+pbr_r) = brbuf_r(1:pbr_r)
         csum = csum + pbr_r
      end if
      if(csum .GT. 0)then
!        write(*,*)'csum, myid: ',csum, myid
      end if
      do idx = 1,psum
        if (.NOT. associated(first_particle)) then
           allocate(first_particle)
           first_particle = totalbuf(idx)
           nullify(first_particle%next,first_particle%prev)
        else
           allocate(first_particle%prev)
           tmp => first_particle%prev
           tmp = totalbuf(idx)
           tmp%next => first_particle
           nullify(tmp%prev)
           first_particle => tmp
           nullify(tmp)
        end if
      end do  
      
      deallocate(rbuf_s,trbuf_s,tbuf_s,tlbuf_s)
      deallocate(lbuf_s,blbuf_s,bbuf_s,brbuf_s)
      deallocate(rbuf_r,trbuf_r,tbuf_r,tlbuf_r)
      deallocate(lbuf_r,blbuf_r,bbuf_r,brbuf_r)
      deallocate(totalbuf)

      call mpi_allreduce(aircraft_moved1,aircraft_moved2,1,mpi_integer,
     +                   mpi_sum,mpi_comm_world,ierr)
      if(myid==0)then
        apr_id = apr_id_old
        if(apr_loc(1) .GE. xmax .or. apr_loc(2) .GE. ymax)then !Outside
c       the bounds of root
          if(aircraft_moved2 .NE. 0)then
            apr_id = aircraft_moved2
          end if
        end if
      end if
     
      call mpi_bcast(apr_id,1,mpi_integer,0,mpi_comm_world,ierr)

      end subroutine particle_exchange

      subroutine fill_ext 
      use pars
      use particles
      use fields
      use con_stats
      use con_data
      implicit none
      include 'mpif.h'

      integer :: istatus(mpi_status_size),ierr
      integer :: ix,iy,iz
      !preceding letter: r=right,l=left,t=top,b=bot.
      !_s: buf of things to send TO r,l,t,b
      !_r: buf of things to recv FROM r,l,t,b 
c      real :: tbuf_s(nnz+2,iye-iys+1,2,4),tbuf_r(nnz+2,iye-iys+1,3,4) !Brian 7/14/14
c      real :: bbuf_s(nnz+2,iye-iys+1,3,4),bbuf_r(nnz+2,iye-iys+1,2,4) !Brian 7/14/14
c      real :: rbuf_s(nnz+2,2,mxe-mxs+1,4),rbuf_r(nnz+2,3,mxe-mxs+1,4) !Brian 7/14/14
c      real :: lbuf_s(nnz+2,3,mxe-mxs+1,4),lbuf_r(nnz+2,2,mxe-mxs+1,4) !Brian 7/14/14
      real :: tbuf_s(nnz+2,iye-iys+1,2,5),tbuf_r(nnz+2,iye-iys+1,3,5) !Brian 7/14/14
      real :: bbuf_s(nnz+2,iye-iys+1,3,5),bbuf_r(nnz+2,iye-iys+1,2,5) !Brian 7/14/14
      real :: rbuf_s(nnz+2,2,mxe-mxs+1,5),rbuf_r(nnz+2,3,mxe-mxs+1,5) !Brian 7/14/14
      real :: lbuf_s(nnz+2,3,mxe-mxs+1,5),lbuf_r(nnz+2,2,mxe-mxs+1,5) !Brian 7/14/14

      !Corners:
c      real :: trbuf_s(nnz+2,2,2,4),trbuf_r(nnz+2,3,3,4) !Brian 7/14/14
c      real :: brbuf_s(nnz+2,2,3,4),brbuf_r(nnz+2,3,2,4) !Brian 7/14/14
c      real :: blbuf_s(nnz+2,3,3,4),blbuf_r(nnz+2,2,2,4) !Brian 7/14/14
c      real :: tlbuf_s(nnz+2,3,2,4),tlbuf_r(nnx+2,2,3,4) !Brian 7/14/14
      real :: trbuf_s(nnz+2,2,2,5),trbuf_r(nnz+2,3,3,5) !Brian 7/14/14
      real :: brbuf_s(nnz+2,2,3,5),brbuf_r(nnz+2,3,2,5) !Brian 7/14/14
      real :: blbuf_s(nnz+2,3,3,5),blbuf_r(nnz+2,2,2,5) !Brian 7/14/14
      real :: tlbuf_s(nnz+2,3,2,5),tlbuf_r(nnz+2,2,3,5) !Brian 7/14/14

      !MPI send counts:
      integer :: rc_s,rc_r,trc_s,trc_r,tc_s,tc_r,tlc_s,tlc_r
      integer :: lc_s,lc_r,blc_s,blc_r,bc_s,bc_r,brc_s,brc_r

      !Debugging:
      real :: xv,yv,zv
 
      !To update the particle ODE in time, need the interpolated
      !velocity field
      !This requires filling uext,vext,wext from nearby procs
      uext = 0.0
      vext = 0.0
      wext = 0.0
      Text = 0.0
      T2ext = 0.0   ! Brian 7/14/14

      !FOR CHECKING INTERPOLATION: create an artificial velocity field
      !to interpolate from:
!      do iz=izs-1,ize+1
!        do iy=iys,iye
!          do ix=1,nnx
!            xv = dx*(ix-1)
!            yv = dy*(iy-1)
!            zv = z(iz)
!            u(ix,iy,iz) = 10.0 
!            v(ix,iy,iz) = 0.0 
!            w(ix,iy,iz) = 0.0 
!            u(ix,iy,iz) = 0.3*xv**6 + 0.5*xv*yv**3*zv + 3.0*xv**3*zv**2+
!     +                      0.7*yv**2*zv**3 + 0.5*zv**3 + 0.1*yv**2
!          end do
!        end do
!      end do

      !First fill the center, since this is just u,v,w on that proc:

      !In the column setup, need to tranpose u,v,w first into u_t,v_t,w_t:
      call xtoz_trans(u(1:nnx,iys:iye,izs-1:ize+1),u_t,nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(v(1:nnx,iys:iye,izs-1:ize+1),v_t,nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(w(1:nnx,iys:iye,izs-1:ize+1),w_t,nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(t(1:nnx,iys:iye,1,izs-1:ize+1),T_t,nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(t(1:nnx,iys:iye,2,izs-1:ize+1),T2_t,nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)   !Brian 7/14/14

      uext(0:nnz+1,iys:iye,mxs:mxe) = u_t(0:nnz+1,iys:iye,mxs:mxe)
      vext(0:nnz+1,iys:iye,mxs:mxe) = v_t(0:nnz+1,iys:iye,mxs:mxe)
      wext(0:nnz+1,iys:iye,mxs:mxe) = w_t(0:nnz+1,iys:iye,mxs:mxe)
      Text(0:nnz+1,iys:iye,mxs:mxe) = T_t(0:nnz+1,iys:iye,mxs:mxe)
      T2ext(0:nnz+1,iys:iye,mxs:mxe) = T2_t(0:nnz+1,iys:iye,mxs:mxe) !Brian 7/14/14

      !Recall that SR assign_nbrs assigned rproc,lproc, etc.

      !Going to call 6 sendrecv calls - one for each proc. nbr.:
      
      !Fill the send buffers:
      
c      !I know these are redundant, but so I can keep them straight...
c      tc_s = 4*(nnz+2)*2*(iye-iys+1) !Brian 7/14/14
c      tc_r = 4*(nnz+2)*3*(iye-iys+1) !Brian 7/14/14
c      trc_s = 4*(nnz+2)*2*2 !Brian 7/14/14
c      trc_r = 4*(nnz+2)*3*3 !Brian 7/14/14
c      rc_s = 4*(nnz+2)*(mxe-mxs+1)*2 !Brian 7/14/14
c      rc_r = 4*(nnx+2)*(mxe-mxs+1)*3 !Brian 7/14/14
c      tlc_s = 4*(nnz+2)*3*2 !Brian 7/14/14
c      tlc_r = 4*(nnz+2)*2*3 !Brian 7/14/14
c      bc_s = 4*(nnz+2)*3*(iye-iys+1) !Brian 7/14/14
c      bc_r = 4*(nnz+2)*2*(iye-iys+1) !Brian 7/14/14
c      blc_s = 4*(nnz+2)*3*3 !Brian 7/14/14
c      blc_r = 4*(nnz+2)*2*2 !Brian 7/14/14
c      lc_s = 4*(nnz+2)*(mxe-mxs+1)*3 !Brian 7/14/14
c      lc_r = 4*(nnz+2)*(mxe-mxs+1)*2 !Brian 7/14/14
c      brc_s = 4*(nnz+2)*2*3 !Brian 7/14/14
c      brc_r = 4*(nnz+2)*3*2 !Brian 7/14/14
      !I know these are redundant, but so I can keep them straight...
      tc_s = 5*(nnz+2)*2*(iye-iys+1) !Brian 7/14/14
      tc_r = 5*(nnz+2)*3*(iye-iys+1) !Brian 7/14/14
      trc_s = 5*(nnz+2)*2*2 !Brian 7/14/14
      trc_r = 5*(nnz+2)*3*3 !Brian 7/14/14
      rc_s = 5*(nnz+2)*(mxe-mxs+1)*2 !Brian 7/14/14
      rc_r = 5*(nnz+2)*(mxe-mxs+1)*3 !Brian 7/14/14
      tlc_s = 5*(nnz+2)*3*2 !Brian 7/14/14
      tlc_r = 5*(nnz+2)*2*3 !Brian 7/14/14
      bc_s = 5*(nnz+2)*3*(iye-iys+1) !Brian 7/14/14
      bc_r = 5*(nnz+2)*2*(iye-iys+1) !Brian 7/14/14
      blc_s = 5*(nnz+2)*3*3 !Brian 7/14/14
      blc_r = 5*(nnz+2)*2*2 !Brian 7/14/14
      lc_s = 5*(nnz+2)*(mxe-mxs+1)*3 !Brian 7/14/14
      lc_r = 5*(nnz+2)*(mxe-mxs+1)*2 !Brian 7/14/14
      brc_s = 5*(nnz+2)*2*3 !Brian 7/14/14
      brc_r = 5*(nnz+2)*3*2 !Brian 7/14/14
     
      !First u:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,1) = u_t(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,1) = u_t(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,1) = u_t(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,1) = u_t(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,1) = u_t(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,1) = u_t(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,1) = u_t(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,1) = u_t(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !v:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,2) = v_t(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,2) = v_t(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,2) = v_t(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,2) = v_t(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,2) = v_t(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,2) = v_t(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,2) = v_t(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,2) = v_t(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !w:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,3) = w_t(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,3) = w_t(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,3) = w_t(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,3) = w_t(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,3) = w_t(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,3) = w_t(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,3) = w_t(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,3) = w_t(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !T:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,4) = T_t(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,4) = T_t(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,4) = T_t(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,4) = T_t(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,4) = T_t(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,4) = T_t(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,4) = T_t(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,4) = T_t(0:nnz+1,iys:iys+2,mxe-1:mxe)
     
      !T2:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,5) =
     +      T2_t(0:nnz+1,iys:iye,mxe-1:mxe) !Brian 7/14/14
      trbuf_s(1:nnz+2,1:2,1:2,5) = 
     +      T2_t(0:nnz+1,iye-1:iye,mxe-1:mxe) !Brian 7/14/14
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,5) = 
     +      T2_t(0:nnz+1,iye-1:iye,mxs:mxe) !Brian 7/14/14
      brbuf_s(1:nnz+2,1:2,1:3,5) = 
     +      T2_t(0:nnz+1,iye-1:iye,mxs:mxs+2) !Brian 7/14/14
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,5) = 
     +      T2_t(0:nnz+1,iys:iye,mxs:mxs+2) !Brian 7/14/14
      blbuf_s(1:nnz+2,1:3,1:3,5) = 
     +      T2_t(0:nnz+1,iys:iys+2,mxs:mxs+2) !Brian 7/14/14
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,5) = 
     +      T2_t(0:nnz+1,iys:iys+2,mxs:mxe) !Brian 7/14/14
      tlbuf_s(1:nnz+2,1:3,1:2,5) = 
     +      T2_t(0:nnz+1,iys:iys+2,mxe-1:mxe) !Brian 7/14/14
     
      !Zero out recieve buffers 
      rbuf_r=0.0;trbuf_r=0.0;tbuf_r=0.0;tlbuf_r=0.0;lbuf_r=0.0
      blbuf_r=0.0;bbuf_r=0.0;brbuf_r=0.0

      !Left/right:
      call MPI_Sendrecv(rbuf_s,rc_s,mpi_real8,rproc,3,
     +        lbuf_r,lc_r,mpi_real8,lproc,3,mpi_comm_world,istatus,ierr)

      call mpi_barrier(mpi_comm_world,ierr)
      call MPI_Sendrecv(lbuf_s,lc_s,mpi_real8,lproc,4,
     +        rbuf_r,rc_r,mpi_real8,rproc,4,mpi_comm_world,istatus,ierr)

      !Top/bottom:
      call MPI_Sendrecv(tbuf_s,tc_s,mpi_real8,tproc,5,
     +        bbuf_r,bc_r,mpi_real8,bproc,5,mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(bbuf_s,bc_s,mpi_real8,bproc,6,
     +        tbuf_r,tc_r,mpi_real8,tproc,6,mpi_comm_world,istatus,ierr)

      !Top right/bottom left:
      call MPI_Sendrecv(trbuf_s,trc_s,mpi_real8,trproc,7,
     +        blbuf_r,blc_r,mpi_real8,blproc,7,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(blbuf_s,blc_s,mpi_real8,blproc,8,
     +        trbuf_r,trc_r,mpi_real8,trproc,8,
     +        mpi_comm_world,istatus,ierr)

       !Top left/bottom right:
      call MPI_Sendrecv(tlbuf_s,tlc_s,mpi_real8,tlproc,9,
     +        brbuf_r,brc_r,mpi_real8,brproc,9,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(brbuf_s,brc_s,mpi_real8,brproc,10,
     +         tlbuf_r,tlc_r,mpi_real8,tlproc,10,
     +         mpi_comm_world,istatus,ierr)

      !Now fill the ext arrays with the recieved buffers:
      uext(0:nnz+1,iys:iye,mxe+1:mxe+3) = 
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,1)
      uext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,1)
      uext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,1)
      uext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,1)
      uext(0:nnz+1,iys:iye,mxs-2:mxs-1) = 
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,1)
      uext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,1)
      uext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,1)
      uext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,1)
 
      vext(0:nnz+1,iys:iye,mxe+1:mxe+3) = 
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,2)
      vext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,2)
      vext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,2)
      vext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,2)
      vext(0:nnz+1,iys:iye,mxs-2:mxs-1) = 
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,2)
      vext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,2)
      vext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,2)
      vext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,2)

      wext(0:nnz+1,iys:iye,mxe+1:mxe+3) = 
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,3)
      wext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,3)
      wext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,3)
      wext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,3)
      wext(0:nnz+1,iys:iye,mxs-2:mxs-1) = 
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,3)
      wext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,3)
      wext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,3)
      wext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,3)

      Text(0:nnz+1,iys:iye,mxe+1:mxe+3) = 
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,4)
      Text(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,4)
      Text(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,4)
      Text(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,4)
      Text(0:nnz+1,iys:iye,mxs-2:mxs-1) = 
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,4)
      Text(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,4)
      Text(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,4)
      Text(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,4)


      T2ext(0:nnz+1,iys:iye,mxe+1:mxe+3) = 
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,5) !Brian 7/14/14
      T2ext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = 
     +     trbuf_r(1:nnz+2,1:3,1:3,5)!Brian 7/14/14
      T2ext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,5)!Brian 7/14/14
      T2ext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = 
     +     brbuf_r(1:nnz+2,1:3,1:2,5)!Brian 7/14/14
      T2ext(0:nnz+1,iys:iye,mxs-2:mxs-1) = 
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,5)!Brian 7/14/14
      T2ext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = 
     +     blbuf_r(1:nnz+2,1:2,1:2,5)!Brian 7/14/14
      T2ext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,5)!Brian 7/14/14
      T2ext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = 
     +     tlbuf_r(1:nnz+2,1:2,1:3,5)!Brian 7/14/14

      end subroutine fill_ext

      subroutine fill_extSFS ! Indrajith 06/201
c       This subroutine calculte the extented fields for SFS velocity
c       calculations

      use pars
      use particles
      use fields
      use con_stats
      use con_data
      implicit none
      include 'mpif.h'

      integer :: istatus(mpi_status_size),ierr
      integer :: ix,iy,iz
      real :: sigm_st(0:nnz+1,iys:iye,mxs:mxe)
      real :: sigm_sdxt(0:nnz+1,iys:iye,mxs:mxe)
      real :: sigm_sdyt(0:nnz+1,iys:iye,mxs:mxe)
      real :: sigm_sdzt(0:nnz+1,iys:iye,mxs:mxe)
      real :: vis_st(0:nnz+1,iys:iye,mxs:mxe)
      !preceding letter: r=right,l=left,t=top,b=bot.
      !_s: buf of things to send TO r,l,t,b
      !_r: buf of things to recv FROM r,l,t,b
      real :: tbuf_s(nnz+2,iye-iys+1,2,5),tbuf_r(nnz+2,iye-iys+1,3,5)
      real :: bbuf_s(nnz+2,iye-iys+1,3,5),bbuf_r(nnz+2,iye-iys+1,2,5)
      real :: rbuf_s(nnz+2,2,mxe-mxs+1,5),rbuf_r(nnz+2,3,mxe-mxs+1,5)
      real :: lbuf_s(nnz+2,3,mxe-mxs+1,5),lbuf_r(nnz+2,2,mxe-mxs+1,5)

      !Corners:
      real :: trbuf_s(nnz+2,2,2,5),trbuf_r(nnz+2,3,3,5) !Brian 7/14/14
      real :: brbuf_s(nnz+2,2,3,5),brbuf_r(nnz+2,3,2,5) !Brian 7/14/14
      real :: blbuf_s(nnz+2,3,3,5),blbuf_r(nnz+2,2,2,5) !Brian 7/14/14
      real :: tlbuf_s(nnz+2,3,2,5),tlbuf_r(nnz+2,2,3,5) !Brian 7/14/14


      !MPI send counts:
      integer :: rc_s,rc_r,trc_s,trc_r,tc_s,tc_r,tlc_s,tlc_r
      integer :: lc_s,lc_r,blc_s,blc_r,bc_s,bc_r,brc_s,brc_r
      sigm_sext = 0.0
      sigm_sdxext = 0.0
      sigm_sdyext = 0.0
      sigm_sdzext = 0.0
      vis_sext = 0.0
      vis_st = 0.0
      sigm_st = 0.0
      sigm_sdxt = 0.0
      sigm_sdyt = 0.0
      sigm_sdzt = 0.0
      
              !In the column setup, need to tranpose u,v,w first into
              !u_t,v_t,w_t:
      call xtoz_trans(sigm_s(1:nnx,iys:iye,izs-1:ize+1),sigm_st,nnx,
     +                nnz,mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(sigm_sdx(1:nnx,iys:iye,izs-1:ize+1),sigm_sdxt,nnx,
     +                nnz,mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(sigm_sdy(1:nnx,iys:iye,izs-1:ize+1),sigm_sdyt,nnx,
     +                nnz,mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(sigm_sdz(1:nnx,iys:iye,izs-1:ize+1),sigm_sdzt,nnx,
     +                nnz,mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)

      call xtoz_trans(vis_ss(1:nnx,iys:iye,izs-1:ize+1),vis_st,nnx,
     +                nnz,mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)


      sigm_sext(0:nnz+1,iys:iye,mxs:mxe) 
     +                              = sigm_st(0:nnz+1,iys:iye,mxs:mxe)
      sigm_sdxext(0:nnz+1,iys:iye,mxs:mxe) 
     +                              = sigm_sdxt(0:nnz+1,iys:iye,mxs:mxe)
      sigm_sdyext(0:nnz+1,iys:iye,mxs:mxe) 
     +                              = sigm_sdyt(0:nnz+1,iys:iye,mxs:mxe)
      sigm_sdzext(0:nnz+1,iys:iye,mxs:mxe) 
     +                              = sigm_sdzt(0:nnz+1,iys:iye,mxs:mxe)
      vis_sext(0:nnz+1,iys:iye,mxs:mxe)
     +                              = vis_st(0:nnz+1,iys:iye,mxs:mxe)


        
      !Fill the send buffers:
      tc_s = 5*(nnz+2)*2*(iye-iys+1) !Brian 7/14/14
      tc_r = 5*(nnz+2)*3*(iye-iys+1) !Brian 7/14/14
      trc_s = 5*(nnz+2)*2*2 !Brian 7/14/14
      trc_r = 5*(nnz+2)*3*3 !Brian 7/14/14
      rc_s = 5*(nnz+2)*(mxe-mxs+1)*2 !Brian 7/14/14
      rc_r = 5*(nnz+2)*(mxe-mxs+1)*3 !Brian 7/14/14
      tlc_s = 5*(nnz+2)*3*2 !Brian 7/14/14
      tlc_r = 5*(nnz+2)*2*3 !Brian 7/14/14
      bc_s = 5*(nnz+2)*3*(iye-iys+1) !Brian 7/14/14
      bc_r = 5*(nnz+2)*2*(iye-iys+1) !Brian 7/14/14
      blc_s = 5*(nnz+2)*3*3 !Brian 7/14/14
      blc_r = 5*(nnz+2)*2*2 !Brian 7/14/14
      lc_s = 5*(nnz+2)*(mxe-mxs+1)*3 !Brian 7/14/14
      lc_r = 5*(nnz+2)*(mxe-mxs+1)*2 !Brian 7/14/14
      brc_s = 5*(nnz+2)*2*3 !Brian 7/14/14
      brc_r = 5*(nnz+2)*3*2 !Brian 7/14/14
        
      !First sigm_s:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,1) 
     +                          = sigm_st(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,1) 
     +                          = sigm_st(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,1) 
     +                          = sigm_st(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,1) 
     +                          = sigm_st(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,1) 
     +                          = sigm_st(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,1) 
     +                          = sigm_st(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,1) 
     +                          = sigm_st(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,1)
     +                          = sigm_st(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !sigm_sdx:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,2) 
     +                          = sigm_sdxt(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,2) 
     +                          = sigm_sdxt(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,2) 
     +                          = sigm_sdxt(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,2) 
     +                          = sigm_sdxt(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,2) 
     +                          = sigm_sdxt(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,2) 
     +                          = sigm_sdxt(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,2) 
     +                          = sigm_sdxt(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,2) 
     +                          = sigm_sdxt(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !sigm_sdy:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,3) 
     +                         = sigm_sdyt(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,3)
     +                         = sigm_sdyt(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,3) 
     +                         = sigm_sdyt(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,3) 
     +                         = sigm_sdyt(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,3) 
     +                         = sigm_sdyt(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,3) 
     +                         = sigm_sdyt(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,3) 
     +                         = sigm_sdyt(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,3)
     +                         = sigm_sdyt(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !sigm_sdz
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,4) 
     +                          = sigm_sdzt(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,4) 
     +                          = sigm_sdzt(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,4) 
     +                          = sigm_sdzt(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,4) 
     +                          = sigm_sdzt(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,4) 
     +                          = sigm_sdzt(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,4) 
     +                          = sigm_sdzt(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,4) 
     +                          = sigm_sdzt(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,4) 
     +                          = sigm_sdzt(0:nnz+1,iys:iys+2,mxe-1:mxe)
 
      !vis_s
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,5)
     +                          = vis_st(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,5)
     +                          = vis_st(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,5)
     +                          = vis_st(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,5)
     +                          = vis_st(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,5)
     +                          = vis_st(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,5)
     +                          = vis_st(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,5)
     +                          = vis_st(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,5)
     +                          = vis_st(0:nnz+1,iys:iys+2,mxe-1:mxe)


      !Zero out recieve buffers
      rbuf_r=0.0;trbuf_r=0.0;tbuf_r=0.0;tlbuf_r=0.0;lbuf_r=0.0
      blbuf_r=0.0;bbuf_r=0.0;brbuf_r=0.0

      !Left/right:
      call MPI_Sendrecv(rbuf_s,rc_s,mpi_real8,rproc,3,
     +        lbuf_r,lc_r,mpi_real8,lproc,3,mpi_comm_world,istatus,ierr)

      call mpi_barrier(mpi_comm_world,ierr)
      call MPI_Sendrecv(lbuf_s,lc_s,mpi_real8,lproc,4,
     +        rbuf_r,rc_r,mpi_real8,rproc,4,mpi_comm_world,istatus,ierr)

      !Top/bottom:
      call MPI_Sendrecv(tbuf_s,tc_s,mpi_real8,tproc,5,
     +        bbuf_r,bc_r,mpi_real8,bproc,5,mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(bbuf_s,bc_s,mpi_real8,bproc,6,
     +        tbuf_r,tc_r,mpi_real8,tproc,6,mpi_comm_world,istatus,ierr)

      !Top right/bottom left:
      call MPI_Sendrecv(trbuf_s,trc_s,mpi_real8,trproc,7,
     +        blbuf_r,blc_r,mpi_real8,blproc,7,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(blbuf_s,blc_s,mpi_real8,blproc,8,
     +        trbuf_r,trc_r,mpi_real8,trproc,8,
     +        mpi_comm_world,istatus,ierr)

       !Top left/bottom right:
      call MPI_Sendrecv(tlbuf_s,tlc_s,mpi_real8,tlproc,9,
     +        brbuf_r,brc_r,mpi_real8,brproc,9,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(brbuf_s,brc_s,mpi_real8,brproc,10,
     +         tlbuf_r,tlc_r,mpi_real8,tlproc,10,
     +         mpi_comm_world,istatus,ierr)

c     Now fill the ext arrays with the recieved buffers:
      sigm_sext(0:nnz+1,iys:iye,mxe+1:mxe+3) =
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,1)
      sigm_sext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = 
     +     trbuf_r(1:nnz+2,1:3,1:3,1)
      sigm_sext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,1)
      sigm_sext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = 
     +     brbuf_r(1:nnz+2,1:3,1:2,1)
      sigm_sext(0:nnz+1,iys:iye,mxs-2:mxs-1) =
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,1)
      sigm_sext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = 
     +     blbuf_r(1:nnz+2,1:2,1:2,1)
      sigm_sext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,1)
      sigm_sext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = 
     +     tlbuf_r(1:nnz+2,1:2,1:3,1)

      sigm_sdxext(0:nnz+1,iys:iye,mxe+1:mxe+3) =
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,2)
      sigm_sdxext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = 
     +     trbuf_r(1:nnz+2,1:3,1:3,2)
      sigm_sdxext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,2)
      sigm_sdxext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = 
     +     brbuf_r(1:nnz+2,1:3,1:2,2)
      sigm_sdxext(0:nnz+1,iys:iye,mxs-2:mxs-1) =
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,2)
      sigm_sdxext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = 
     +     blbuf_r(1:nnz+2,1:2,1:2,2)
      sigm_sdxext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,2)
      sigm_sdxext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = 
     +     tlbuf_r(1:nnz+2,1:2,1:3,2)
        
      sigm_sdyext(0:nnz+1,iys:iye,mxe+1:mxe+3) =
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,3)
      sigm_sdyext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = 
     +     trbuf_r(1:nnz+2,1:3,1:3,3)
      sigm_sdyext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,3)
      sigm_sdyext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = 
     +     brbuf_r(1:nnz+2,1:3,1:2,3)
      sigm_sdyext(0:nnz+1,iys:iye,mxs-2:mxs-1) =
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,3)
      sigm_sdyext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = 
     +     blbuf_r(1:nnz+2,1:2,1:2,3)
      sigm_sdyext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,3)
      sigm_sdyext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = 
     +     tlbuf_r(1:nnz+2,1:2,1:3,3)

      sigm_sdzext(0:nnz+1,iys:iye,mxe+1:mxe+3) =
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,4)
      sigm_sdzext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = 
     +     trbuf_r(1:nnz+2,1:3,1:3,4)
      sigm_sdzext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,4)
      sigm_sdzext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = 
     +     brbuf_r(1:nnz+2,1:3,1:2,4)
      sigm_sdzext(0:nnz+1,iys:iye,mxs-2:mxs-1) =
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,4)
      sigm_sdzext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = 
     +     blbuf_r(1:nnz+2,1:2,1:2,4)
      sigm_sdzext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,4)
      sigm_sdzext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = 
     +     tlbuf_r(1:nnz+2,1:2,1:3,4)

     
      vis_sext(0:nnz+1,iys:iye,mxe+1:mxe+3) =
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,5)
      vis_sext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) =
     +     trbuf_r(1:nnz+2,1:3,1:3,5)
      vis_sext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,5)
      vis_sext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) =
     +     brbuf_r(1:nnz+2,1:3,1:2,5)
      vis_sext(0:nnz+1,iys:iye,mxs-2:mxs-1) =
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,5)
      vis_sext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) =
     +     blbuf_r(1:nnz+2,1:2,1:2,5)
      vis_sext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,5)
      vis_sext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) =
     +     tlbuf_r(1:nnz+2,1:2,1:3,5)

      end subroutine fill_extSFS

      subroutine uf_interp
      use pars
      use particles
      use fields
      use con_stats
      use con_data
      implicit none
       
      integer :: ix,iy,izuv,izw,iz,i,k,j
      integer :: first,last
      real :: xkval,xjval,pj,dxvec(2)
      integer :: ijpts(2,6),kuvpts(6),kwpts(6)
      real :: wt(4,6)
      real :: ran2
      
      !Following Orlando's code, get the "leftmost" node
      !This is just the minimum (i,j,k) on the volume 

      ijpts(1,3) = floor(part%xp(1)/dx) + 1 
      ijpts(2,3) = floor(part%xp(2)/dy) + 1
     
      !Fill in the neighbors:
      ijpts(1,2) = ijpts(1,3)-1
      ijpts(1,1) = ijpts(1,2)-1
      ijpts(1,4) = ijpts(1,3)+1
      ijpts(1,5) = ijpts(1,4)+1
      ijpts(1,6) = ijpts(1,5)+1

      ijpts(2,2) = ijpts(2,3)-1
      ijpts(2,1) = ijpts(2,2)-1
      ijpts(2,4) = ijpts(2,3)+1
      ijpts(2,5) = ijpts(2,4)+1
      ijpts(2,6) = ijpts(2,5)+1
     
      !Finding the k-lhnode is different since grid may be stretched
      !AND since (u,v) and w stored differently
      !Will get a k-index for (u,v) and one for w
      
      !Do (u,v) loop first:
      do iz=0,nnz+1
         if (zz(iz) .GT. part%xp(3)) go to 1199
      end do
1199  continue

      kuvpts(3) = iz-1
      !Then fill in the rest:
      kuvpts(4) = kuvpts(3)+1
      kuvpts(5) = kuvpts(4)+1
      kuvpts(6) = kuvpts(5)+1
      kuvpts(2) = kuvpts(3)-1
      kuvpts(1) = kuvpts(2)-1

      !Do again for w:
      do iz = 0,nnz+1 
         if (z(iz) .GT. part%xp(3)) go to 1120
      end do
1120  continue
   
      kwpts(3) = iz-1
      !Then fill in the rest:
      kwpts(4) = kwpts(3)+1
      kwpts(5) = kwpts(4)+1
      kwpts(6) = kwpts(5)+1
      kwpts(2) = kwpts(3)-1
      kwpts(1) = kwpts(2)-1

!---------
!      !As an aside, use ijpts,kwpts to update the particle numbers at each cell:
!
!      !call inputs  are the nodes associated with the volume that each part. lies in
!      !The +1 on kwpts is since stats are stored at uv (center) points in z
!      call particle_stats(ijpts(1,3),ijpts(2,3),kwpts(3)+1)
!   
!---------

      !Fill in the weights:
      !First for x and y since they are periodic:
      wt(1:4,1:6) = 0.0
      dxvec(1) = dx
      dxvec(2) = dy
      do iz = 1,2
      do j = 1,6
         xjval = dxvec(iz)*(ijpts(iz,j)-1)
         pj = 1.0
         do k = 1,6
            xkval = dxvec(iz)*(ijpts(iz,k)-1)
            if (j .NE. k) then
                  pj = pj*(part%xp(iz)-xkval)/(xjval-xkval)
            end if
         end do
         wt(iz,j) = pj
       end do
       end do
      

!      !Enforce periodicity in x by adjusting the indices:
!      !NOTE: doing this after having computed the weights
!      do i=1,6
!         if (ijpts(1,i) .LT. 1) then
!            ijpts(1,i) = ijpts(1,i)+nnx
!         elseif (ijpts(1,i) .GT. nnx) then
!            ijpts(1,i) = ijpts(1,i)-nnx
!         end if
!      end do
     
       !Now compute weights in z-dir
       !There are 2 sections: weights at (u,v) nodes (kuvpts) 
       !And weights computed at w nodes (kwpts)

       !Compute weights at kuvpts
       !Must check to see how close we are to a top/bot boundary
       if (kuvpts(3) == 1) then
          first = 3
          last = 4
          !Set these equal to 1 so uext(-1) won't be accessed
          !Note: the value doesn't matter since weight will be 0
          kuvpts(1) = 1
          kuvpts(2) = 1
       elseif (kuvpts(3) == 0) then
          first = 4
          last = 5
          kuvpts(1) = 1
          kuvpts(2) = 1
          kuvpts(3) = 1
       elseif (kuvpts(3) .LT. 0) then 
          first = 0
          last = 0
       elseif (kuvpts(3) == 2) then 
          first = 2
          last = 5
          kuvpts(1) = 1
          kuvpts(6) = 1
       !Between top cell center and the domain boundary
       elseif (kuvpts(3) == nnz) then
          first = 2
          last = 3
          kuvpts(4) = nnz
          kuvpts(5) = nnz
          kuvpts(6) = nnz
       elseif (kuvpts(3) .GT. nnz) then
          first = 0
          last = 0
       !Between 2nd to last and last cell center at top
          first = 3
          last = 4
          kuvpts(5) = nnz
          kuvpts(6) = nnz
       elseif (kuvpts(3) == nnz-2) then
          first = 2
          last = 5
       else
          first = 1
          last = 6
       end if

       !Recall that wt has been set to zero, so
       !weights will be zero if (first,last) isn't (1,6)
       do j = first,last
           xjval = zz(kuvpts(j))
           pj = 1.0
           do k = first,last
              xkval = zz(kuvpts(k))
              if (j .NE. k) then
                 pj = pj*(part%xp(3)-xkval)/(xjval-xkval)
              end if
           end do
           wt(3,j) = pj
      end do

       !Now compute weights at kwpts
       !Again must check to see how close we are to a top/bot boundary
       if (kwpts(3) == 0) then
          first = 3
          last = 4
          kwpts(1) = 1
          kwpts(2) = 1
       elseif (kwpts(3) .LT. 0) then 
          first = 0
          last = 0
          kwpts(1) = 1
          kwpts(2) = 1
          kwpts(3) = 1
       elseif (kwpts(3) == 1) then 
          first = 2
          last = 5
          kwpts(1) = 1
       elseif (kwpts(3) == nnz-1) then
          first = 3
          last = 4
          kwpts(5) = nnz
          kwpts(6) = nnz                   
       elseif (kwpts(3) .GE. nnz) then
          first = 0
          last = 0
          kwpts(3) = nnz
          kwpts(4) = nnz
          kwpts(5) = nnz
          kwpts(6) = nnz
       elseif (kwpts(3) == nnz-2) then
          first = 2
          last = 5
          kwpts(6) = nnz
       else
          first = 1
          last = 6
       end if

       !Recall that wt has been set to zero, so
       !weights will be zero if (first,last) isn't (1,6)
       do j = first,last
           xjval = z(kwpts(j))
           pj = 1.0
           do k = first,last
              xkval = z(kwpts(k))
              if (j .NE. k) then
                 pj = pj*(part%xp(3)-xkval)/(xjval-xkval)
              end if
           end do
           wt(4,j) = pj
      end do

      !Now we have the weights - compute the velocity at xp:
        part%uf(1:3) = 0.0
        part%Tf = 0.0
        part%qinf = 0.0 ! Brian 7/14/14 
        do k = 1,6
        do j = 1,6
        do i = 1,6
            ix = ijpts(1,i)
            iy = ijpts(2,j)
            izuv = kuvpts(k)
            izw = kwpts(k)
            part%uf(1) = part%uf(1)+uext(izuv,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(3,k) 
            part%uf(2) = part%uf(2)+vext(izuv,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(3,k) 
            part%uf(3) = part%uf(3)+wext(izw,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(4,k) 
            part%Tf = part%Tf+Text(izuv,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(3,k)
            part%qinf = part%qinf+T2ext(izuv,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(3,k)  ! Brian 7/14/14 
         end do
         end do 
         end do


      end subroutine uf_interp 
     
      subroutine sigm_interp(sigm_sdxp,sigm_sdyp,sigm_sdzp,
     +                          vis_sp,iz_part)
      use pars
      use particles
      use fields
      use con_stats
      use con_data
      implicit none

      integer :: ix,iy,izuv,izw,iz,i,k,j
      integer :: first,last
      real :: xkval,xjval,pj,dxvec(2)
      integer :: ijpts(2,6),kuvpts(6),kwpts(6),iz_part
      real :: wt(4,6)
      real :: ran2
      real :: sigm_sdxp,sigm_sdyp,sigm_sdzp,vis_sp
      !Following Orlando's code, get the "leftmost" node
      !This is just the minimum (i,j,k) on the volume

      ijpts(1,3) = floor(part%xp(1)/dx) + 1
      ijpts(2,3) = floor(part%xp(2)/dy) + 1

      !Fill in the neighbors:
      ijpts(1,2) = ijpts(1,3)-1
      ijpts(1,1) = ijpts(1,2)-1
      ijpts(1,4) = ijpts(1,3)+1
      ijpts(1,5) = ijpts(1,4)+1
      ijpts(1,6) = ijpts(1,5)+1

      ijpts(2,2) = ijpts(2,3)-1
      ijpts(2,1) = ijpts(2,2)-1
      ijpts(2,4) = ijpts(2,3)+1
      ijpts(2,5) = ijpts(2,4)+1
      ijpts(2,6) = ijpts(2,5)+1

      !Finding the k-lhnode is different since grid may be stretched
      !AND since (u,v) and w stored differently
      !Will get a k-index for (u,v) and one for w

      !Do (u,v) loop first:
      do iz=0,nnz+1
         if (zz(iz) .GT. part%xp(3)) go to 1199
      end do
1199  continue
c      iz_part = iz-1
      kuvpts(3) = iz-1
      !Then fill in the rest:
      kuvpts(4) = kuvpts(3)+1
      kuvpts(5) = kuvpts(4)+1
      kuvpts(6) = kuvpts(5)+1
      kuvpts(2) = kuvpts(3)-1
      kuvpts(1) = kuvpts(2)-1

      !Do again for w:
      do iz = 0,nnz+1
         if (z(iz) .GT. part%xp(3)) go to 1120
      end do
1120  continue
      iz_part = iz-1

      kwpts(3) = iz-1
      !Then fill in the rest:
      kwpts(4) = kwpts(3)+1
      kwpts(5) = kwpts(4)+1
      kwpts(6) = kwpts(5)+1
      kwpts(2) = kwpts(3)-1
      kwpts(1) = kwpts(2)-1
 
      !Fill in the weights:
      !First for x and y since they are periodic:
      wt(1:4,1:6) = 0.0
      dxvec(1) = dx
      dxvec(2) = dy
      do iz = 1,2
      do j = 1,6
         xjval = dxvec(iz)*(ijpts(iz,j)-1)
         pj = 1.0
         do k = 1,6
            xkval = dxvec(iz)*(ijpts(iz,k)-1)
            if (j .NE. k) then
                  pj = pj*(part%xp(iz)-xkval)/(xjval-xkval)
            end if
         end do
         wt(iz,j) = pj
       end do
       end do
       !Compute weights at kuvpts
       !Must check to see how close we are to a top/bot boundary
       if (kuvpts(3) == 1) then
          first = 3
          last = 4
          !Set these equal to 1 so uext(-1) won't be accessed
          !Note: the value doesn't matter since weight will be 0
          kuvpts(1) = 1
          kuvpts(2) = 1
       elseif (kuvpts(3) == 0) then
          first = 4
          last = 5
          kuvpts(1) = 1
          kuvpts(2) = 1
          kuvpts(3) = 1
       elseif (kuvpts(3) .LT. 0) then
          first = 0
          last = 0
       elseif (kuvpts(3) == 2) then
          first = 2
          last = 5
          kuvpts(1) = 1
          kuvpts(6) = 1
       !Between top cell center and the domain boundary
       elseif (kuvpts(3) == nnz) then
          first = 2
          last = 3
          kuvpts(4) = nnz
          kuvpts(5) = nnz
          kuvpts(6) = nnz
       elseif (kuvpts(3) .GT. nnz) then
          first = 0
          last = 0
       !Between 2nd to last and last cell center at top
       elseif (kuvpts(3) == nnz-1) then
          first = 3
          last = 4
          kuvpts(5) = nnz
          kuvpts(6) = nnz
       elseif (kuvpts(3) == nnz-2) then
          first = 2
          last = 5
          kuvpts(1) = nnz
          kuvpts(6) = nnz
       else
          first = 1
          last = 6
       end if

       !Recall that wt has been set to zero, so
       !weights will be zero if (first,last) isn't (1,6)
       do j = first,last
           xjval = zz(kuvpts(j))
           pj = 1.0
           do k = first,last
              xkval = zz(kuvpts(k))
              if (j .NE. k) then
                 pj = pj*(part%xp(3)-xkval)/(xjval-xkval)
              end if
           end do
           wt(3,j) = pj
      end do
       !Now compute weights at kwpts
       !Again must check to see how close we are to a top/bot boundary
       if (kwpts(3) == 0) then
          first = 3
          last = 4
          kwpts(1) = 1
          kwpts(2) = 1
       elseif (kwpts(3) .LT. 0) then
          first = 0
          last = 0
          kwpts(1) = 1
          kwpts(2) = 1
          kwpts(3) = 1
       elseif (kwpts(3) == 1) then
          first = 2
          last = 5
          kwpts(1) = 1
       elseif (kwpts(3) == nnz-1) then
          first = 3
          last = 4
          kwpts(5) = nnz
          kwpts(6) = nnz
       elseif (kwpts(3) .GE. nnz) then
          first = 0
          last = 0
          kwpts(3) = nnz
          kwpts(4) = nnz
          kwpts(5) = nnz
          kwpts(6) = nnz
       elseif (kwpts(3) == nnz-2) then
          first = 2
          last = 5
          kwpts(6) = nnz
       else
          first = 1
          last = 6
       end if
       !Recall that wt has been set to zero, so
       !weights will be zero if (first,last) isn't (1,6)
       do j = first,last
           xjval = z(kwpts(j))
           pj = 1.0
           do k = first,last
              xkval = z(kwpts(k))
              if (j .NE. k) then
                 pj = pj*(part%xp(3)-xkval)/(xjval-xkval)
              end if
           end do
           wt(4,j) = pj
      end do

      !Now we have the weights - compute the sigma_s and derivatives at xp:
c        part%uf(1:3) = 0.0
c        part%Tf = 0.0
c        part%qinf = 0.0 ! Brian 7/14/14
        part%sigm_s = 0.0
        sigm_sdxp = 0.0
        sigm_sdyp = 0.0
        sigm_sdzp = 0.0
        vis_sp = 0.0
        do k = 1,6
        do j = 1,6
        do i = 1,6
            ix = ijpts(1,i)
            iy = ijpts(2,j)
            izuv = kuvpts(k)
            izw = kwpts(k)
            part%sigm_s = part%sigm_s+sigm_sext(izw,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(4,k)
            sigm_sdxp = sigm_sdxp+ sigm_sdxext(izw,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(4,k)
            sigm_sdyp = sigm_sdyp+sigm_sdyext(izw,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(4,k)
            sigm_sdzp = sigm_sdzp+sigm_sdzext(izuv,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(3,k)

            vis_sp =vis_sp +vis_sext(izw,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(4,k)

         if (isnan(vis_sp)) then
         write(*,*) 'P:',vis_sext(izw,iy,ix),wt(1,i),wt(2,j),wt(4,k),
     +   part%xp(3),i,j,k,kwpts(3)
         end if

  
         end do
         end do
         end do


      end subroutine sigm_interp

      subroutine particle_init
      use particles
      use pars
      use con_data
      implicit none
      include 'mpif.h' 
      integer :: values(8)
      integer :: idx,ind
      integer :: ierr
      integer :: bin_location
      real :: xv,yv,zv,ran2,deltaz
      real :: maxx,maxy,maxz
      real :: rpmax,rpmin,rp_init !max radius, min radius, particle radiuis
      real :: clssize_log, log_ri
      real :: pi,init_apr_loc1
      !Create the seed for the random number generator:
     
      call date_and_time(VALUES=values)
      iseed = -(myid+values(8)+values(7)+values(6))

      !Initialize ngidx, the particle global index
      ngidx = 1
      Npdx  = 0 !set all Npdx to 0 at first particle init
      init_apr_loc1=0
      !Set the radius range Indrajith 01/03/2017
      rpmax = radius_max  ! maximum initial particle radius
      rpmin = radius_min  ! minimum initial particle radius
      !For the channel case, set the total number of particles:
      deltaz = zmax-zmin
c      write(*,*)'num-part',tnumpart     !DEBUG ONLY
      if(iPprod==0)then
        numpart=tnumpart/numprocs
        if(myid==0)then
          numpart = numpart + MOD(tnumpart,numprocs)
        end if
      else
c       Making p_l (point location) at the center of the grid
        p_l=xl/2.0
        ipointsource=0
        if(xmin .le. p_l .and. xmax .gt. p_l)then
          if(ymin .le. p_l .and. ymax .gt. p_l)then
            write(*,*)'proc_id',myid
            ipointsource=1
          end if
        end if
      end if

      !Initialize the linked list of particles:
      nullify(part,first_particle)
      
      
      ! For particle histogram calculations
c____________________________________________
      if(ipolydisp == 1)then
        allocate (hist_radclass(Nrp+2))
        allocate (del_radclass(Nrp+2))
        allocate (hist_req_init(Nrp+2))
        allocate (bin_radclass(Nrp+2))
        allocate (class_maxr(Nrp+2))
        allocate (class_minr(Nrp+2))
        allocate (init_radclass(Nrp+2))

        hist_radclass = 0
        del_radclass = 0
        hist_req_init = 0
        bin_radclass = 0.0
        class_maxr = 0.0
        class_minr = 0.0
        init_radclass = 0

        bin_radclass(1) = rpmin/2.0 ! dummy class

        !Checking params are fine..
        if(rpmin .gt. rpmax)then
          write(*,*)'how is your rpmin greater than rpmax'
          stop
        end if
        clssize_log = (log10(rpmax)-log10(rpmin))/Nrp
        do ind =2,Nrp + 2
          log_ri = log10(rpmin) +(ind-2)*clssize_log
          class_maxr(ind) = 10**(log_ri+clssize_log)
          class_minr(ind) = 10**(log_ri)
          bin_radclass(ind) = 10**(log_ri+clssize_log/2)
          if (myid .eq. 0)then
            write(*,*)'drop_class',ind,':'
            write(*,*)'bin center: ',bin_radclass(ind)
            write(*,*)'bin minimum radius: ',class_minr(ind)
            write(*,*)'bin maximum radius: ',class_maxr(ind)
          end if
        end do
        class_maxr(1)= class_minr(2)
      end if
c____________________________________________

      !Now initialize all particles with a random (or specified)  location on that processor
      maxx=0.0
      maxy=0.0
      maxz=0.0
      rp_init = 0.0 

!     Aircraft particle creation
      if(iaircraft .eq. 1)then
!       Restrict within the bounds of the processor
        if(   xmin<=apr_init(1) .and. apr_init(1)<=xmax
     +  .and. ymin<=apr_init(2) .and. apr_init(2)<ymax)then
          write(*,*)'myid and zmin zmax',myid,zmin,zmax
          write(*,*)'myid and xmin xmax',myid,xmin,xmax
          write(*,*)'myid and ymin ymax',myid,ymin,ymax
          Npdx = Npdx + 1
          init_apr_id = myid
          call create_particle(apr_init,apr_v,
     +         0,0,Npdx,radius_max,myid,time,1)
          write(*,*)'myid for aircraft: ',myid
          write(*,*)'aircraft velocity: ',apr_v(1),apr_v(2),apr_v(3)
          write(*,*)'Npdx: ',Npdx
        end if

        call mpi_allreduce(init_apr_id,apr_id,1,mpi_integer8,
     +                     mpi_sum,mpi_comm_world,ierr)
        call mpi_bcast(apr_id,1,mpi_integer8,0,mpi_comm_world,ierr)
      end if


      !Whether a point source
      if(iPprod==0)then
        do idx=1,numpart
          xv = ran2(iseed)*(xmax-xmin) + xmin
          yv = ran2(iseed)*(ymax-ymin) + ymin
          zv = ran2(iseed)*(zi-zw1) + zw1 
        ! Calculate particle radius
          if(ipolydisp .eq. 1)then
c           random distribution within the size bins
            call create_random_partsize(rp_init)
            bin_location = minloc(class_maxr,DIM=1,
     +                            MASK=(rp_init<class_maxr))  
            init_radclass(bin_location) = 
     +                    init_radclass(bin_location) + 1
          else
            rp_init = radius_max
          end if
          radius_mass = rp_init
          pi          = 4.0*atan(1.0)
          m_s         = Sal/1000.0*3.0/4.0*pi*radius_mass**3.0*rhow
          Npdx        = Npdx + 1        
          call create_particle((/xv,yv,zv/),vp_init,
     +         Tp_init,qf_init,Npdx,rp_init,myid,time,0.0)
        end do
      else
        if(ipointsource==1)then
c       Locations to where to single point inject particles.
          xv=xl/2.0
          yv=yl/2.0
c       Sourced pretty close to the surface
          zv=zw1
          numpart=tnumpart
          do idx=1,numpart
            if(ipolydisp==1)then
c             random distribution within the size bins
              call create_random_partsize(rp_init)
              bin_location = minloc(class_maxr,DIM=1,
     +                              MASK=(rp_init<class_maxr)) 
              init_radclass(bin_location) = 
     +                      init_radclass(bin_location) + 1
            else
              rp_init   = radius_max
            end if
            radius_mass = rp_init
            pi          = 4.0*atan(1.0)
            m_s         = Sal/1000.0*3.0/4.0*pi*radius_mass**3.0*rhow
            Npdx        = Npdx + 1
            call create_particle((/xv,yv,zv/),vp_init,
     +           Tp_init,qf_init,Npdx,rp_init,myid,time,0.0)

          end do
        end if
      end if

      partTsrc = 0.0
      partTsrc_t = 0.0
      partHsrc = 0.0
      partHsrc_t = 0.0
      partTEsrc = 0.0
      partTEsrc_t = 0.0

      end subroutine particle_init


      subroutine particle_setup

      use particles
      use pars
      implicit none 
      include 'mpif.h'

      integer :: blcts(3),types(3)
      integer :: ierr
      real :: pi
      integer(kind=MPI_ADDRESS_KIND) :: extent,lb
      integer(kind=MPI_ADDRESS_KIND) :: extent2,lb2,displs(3)
      integer :: num_reals

      !First set up the neighbors for the interpolation stage:
      call assign_nbrs

      !Also assign the x,y,z max and mins to track particles leaving
      call set_part_bounds

c      pi   = 4.0*atan(1.0)  
    

c      m_s = Sal/1000.0*4.0/3.0*pi*radius_mass**3*rhow

c     settling velocity for eulerian information based on radius_max
      ws = grav*(rhow*(2*radius_max)**2/18.0/nuf/rhoa)

c       **************************************************************
c      ws = 1e-6  !***************************************REMOVE LATER
c       *************************************************************

      if (myid==0) write(*,*) 'ws (Eulerian) part_setup = ',ws
      
c       Indrajith 03/13/2017
      if (myid==0) write(*,*) 'Init Particle max_radius = ',radius_max  
      if (myid==0) write(*,*) 'Init Particle min_radius = ',radius_min 
      if (myid==0) write(*,*) 'Particle Salinity = ',Sal     
      !if (myid==0) write(*,*) 'Particle density = ',rhop !rhop is not calculated yet
      !if (myid==0) write(*,*) 'Particle taup = ',1.0/taup_i
      if (myid==0) write(*,*) 'Particle gravity = ',part_grav
      if (myid==0) write(*,*) 'nuf = ',nuf
      if (myid==0) write(*,*) 'CpaCpp = ',CpaCpp
      if (myid==0) write(*,*) 'Pra = ',Pra

      !Initialize the linked list of particles:
      nullify(part,first_particle)

      !Set up MPI datatypes for sending particle information
      !MUST UPDATE IF THINGS ARE ADDED/REMOVED FROM PARTICLE STRUCTURE

       num_reals = 9*3+18
      blcts(1:3) = (/4,num_reals,2/)
      displs(1) = 0
      types(1) = mpi_integer
      call mpi_type_get_extent(mpi_integer,lb,extent,ierr)
      
      !Displace 4*size of mpi_integer (4 integer:
      !pidx,procidx,type_part,dummy)
      displs(2) = extent*4
      types(2) = mpi_real8
      call mpi_type_get_extent(mpi_real8,lb,extent,ierr)
      !Displace num_reals*size of mpi_real8
      displs(3) = displs(2) + extent*num_reals
      types(3) = mpi_integer8

      !Now define the type:
      call mpi_type_create_struct(3,blcts,
     +            displs,types,particletype,ierr)


       call mpi_type_get_true_extent(particletype,lb2,extent2,ierr)
       call mpi_type_get_extent(particletype,lb2,extent,ierr)
       if (extent .NE. sizeof(part) ) then
          if (myid==0) then
          write(*,*) 'WARNING: extent of particletype not equal
     +                  to sizeof(part):'
          write(*,*) 'sizeof(part) = ', sizeof(part)
!          write(*,*) 'sizeof(part%pidx) = ', sizeof(part%pidx)
          write(*,*) 'mpi_type_get_true_extent(particletype) = ',extent2
          write(*,*) 'mpi_type_get_extent(particletype) = ',extent
          end if
       end if
      
      !Need to compute any padding which may exist in particle struct:
      pad_diff = extent-extent2 
      if (myid==0) then
      write(*,*) 'mpi_get_extent = ',extent
      write(*,*) 'mpi_get_true_extent = ',extent2
      write(*,*) 'sizeof(part) = ',sizeof(part)
      write(*,*) 'DIFF = ',pad_diff
      end if
      if (pad_diff .LT. 0) then
        write(*,*) 'WARNING: mpi_get_extent - mpi_get_true_extent LT 0!'
        call mpi_finalize(ierr)
        stop
      end if
      
      if (myid==0) then
      write(*,*) 'huge(tnumpart) = ',huge(tnumpart)
      write(*,*) 'huge(part%pidx) = ',huge(part%pidx)
      end if


      call mpi_type_commit(particletype,ierr)

      end subroutine particle_setup
      subroutine assign_nbrs
        use pars
        use particles
        include 'mpif.h'
      !Figure out which processors lie to all sides: 
      !NOTE: For this updated case, where particles lie in columns not 
      !aligning with the velocity, there will be no MPI_PROC_NULL since
      !x and y are BOTH periodic
     
      !On right boundary:
      if ( mod(myid+1,ncpu_s) == 0 ) then
         !On the top:
         if ( myid .GE. ncpu_s*(ncpu_z-1) ) then
            rproc = myid-ncpu_s+1
            trproc = 0 
            tproc = ncpu_s-1 
            tlproc = ncpu_s-2 
            lproc = myid-1
            blproc = myid-ncpu_s-1
            bproc = myid-ncpu_s
            brproc = myid-ncpu_s - ncpu_s+1
         !On the bottom:
         elseif ( myid .LT. ncpu_s ) then
            rproc = myid-ncpu_s+1
            trproc = myid+1
            tproc = myid+ncpu_s
            tlproc = myid+ncpu_s-1
            lproc = myid-1
            blproc = myid+ncpu_s*(ncpu_z-1)-1 
            bproc = myid+ncpu_s*(ncpu_z-1) 
            brproc = ncpu_s*(ncpu_z-1) 
         !In the middle of right side:
         else 
            rproc = myid-ncpu_s+1
            trproc = myid+1
            tproc = myid+ncpu_s
            tlproc = myid+ncpu_s-1
            lproc = myid-1
            blproc = myid-ncpu_s-1
            bproc = myid-ncpu_s
            brproc = myid-ncpu_s - ncpu_s+1
         end if 

      !On the left boundary:
      elseif ( mod(myid,ncpu_s) == 0) then
         !On the top:
         if ( myid .GE. ncpu_s*(ncpu_z-1) ) then
            rproc = myid+1
            trproc = 1 
            tproc = 0 
            tlproc = ncpu_s-1
            lproc = myid+ncpu_s-1
            blproc = myid-1
            bproc = myid-ncpu_s
            brproc = myid-ncpu_s+1
         !On the bottom:
         elseif ( myid .LT. ncpu_s ) then
            rproc = myid+1
            trproc = myid+ncpu_s+1
            tproc = myid+ncpu_s
            tlproc = myid+ncpu_s+ncpu_s-1
            lproc = myid+ncpu_s-1
            blproc = numprocs-1 
            bproc = ncpu_s*(ncpu_z-1) 
            brproc = ncpu_s*(ncpu_z-1)+1 
         !In the middle of left side:
         else
            rproc = myid+1
            trproc = myid+ncpu_s+1
            tproc = myid+ncpu_s
            tlproc = myid+ncpu_s + ncpu_s-1
            lproc = myid+ncpu_s-1
            blproc = myid-1
            bproc = myid-ncpu_s
            brproc = myid-ncpu_s+1
         end if
      !On the top boundary
      elseif ( myid .GE. ncpu_s*(ncpu_z-1) ) then
         !Only check if in the middle:
         if ( .NOT. ( mod(myid,ncpu_s) == 0) ) then
            if ( .NOT. (mod(myid+1,ncpu_s) == 0) ) then
               rproc = myid+1
               trproc = myid-(ncpu_s*(ncpu_z-1))+1 
               tproc = myid-(ncpu_s*(ncpu_z-1)) 
               tlproc = myid-(ncpu_s*(ncpu_z-1))-1 
               lproc = myid-1
               blproc = myid-ncpu_s-1
               bproc = myid-ncpu_s
               brproc = myid-ncpu_s+1
            end if
         end if 
      !On the bottom boundary
      elseif ( myid .LT. ncpu_s) then
         if ( .NOT. ( mod(myid,ncpu_s) == 0) ) then
            if ( .NOT. (mod(myid+1,ncpu_s) == 0) ) then
               rproc = myid+1
               trproc = myid+ncpu_s+1
               tproc = myid+ncpu_s
               tlproc = myid+ncpu_s-1
               lproc = myid-1
               blproc = myid+ncpu_s*(ncpu_z-1)-1
               bproc = myid+ncpu_s*(ncpu_z-1) 
               brproc = myid+ncpu_s*(ncpu_z-1)+1 
            end if
         end if
      !Everywhere else:
      else 
         rproc = myid+1
         trproc = myid+ncpu_s+1
         tproc = myid+ncpu_s
         tlproc = myid+ncpu_s-1
         lproc = myid-1
         blproc = myid-ncpu_s-1
         bproc = myid-ncpu_s
         brproc = myid-ncpu_s+1
      end if

      return
      end subroutine assign_nbrs
      subroutine set_part_bounds  
      use pars
      use particles
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'
      integer :: ierr

      !Each processor must figure out at what ymin,ymax,zmin,zmax a particle leaves
      ymin        = dy*(iys-1)
      ymax        = dy*(iye)
!     Apparently this is [0 zl] for vertical extent..
      zmin        = z(izs-1)
      zmax        = z(ize)  
      xmin        = dx*(mxs-1)
      xmax        = dx*(mxe)
      if(myid==0)then
        root_xmin = xmin
        root_xmax = xmax
        root_ymin = ymin
        root_ymax = ymax
      end if
      call mpi_bcast(root_xmin,1,mpi_real8,
     +               0,mpi_comm_world,ierr)
      call mpi_bcast(root_xmax,1,mpi_real8,
     +               0,mpi_comm_world,ierr)
      call mpi_bcast(root_ymin,1,mpi_real8,
     +               0,mpi_comm_world,ierr)
      call mpi_bcast(root_ymax,1,mpi_real8,
     +               0,mpi_comm_world,ierr)
      if(myid==0)then
        write(*,*)'base xmin, xmax, ymin, ymax: '
        write(*,*) root_xmin, ' ',root_xmax, ' ',root_ymin,
     +             ' ',root_ymax
      end if


      end subroutine set_part_bounds
      subroutine save_particles
      use particles
      use pars
      implicit none
      include 'mpif.h'

      integer :: istatus(mpi_status_size), ierr, fh
      integer(kind=mpi_offset_kind) :: zoffset,offset
      integer :: pnum_vec(numprocs)
      integer :: iproc,i
      type(particle) :: writebuf(numpart),tmp

      !Do this with mpi_write_at_all
      !Need to figure out the displacements - need numpart from each proc
      call mpi_allgather(numpart,1,mpi_integer,pnum_vec,1,mpi_integer,
     +                   mpi_comm_world,ierr)

      !Package all the particles into writebuf:
      i = 1
      part => first_particle
      do while (associated(part))
      writebuf(i) = part
      !write(*,'a5,3e15.6') 'xp:',part%xp(1:3)
      part => part%next
      i = i + 1
      end do

      !Now only write to the file if you actually have particles
      !EXCEPTION: proc 0, which needs to write tnumpart regardless
      call mpi_file_open(mpi_comm_world, path_sav_part,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null,fh,ierr)

      zoffset = 0
      !Write tnumpart first:
      if (myid==0) then
      call mpi_file_write_at(fh,zoffset,tnumpart,1,mpi_integer,
     +                       istatus,ierr)
      write(*,*) 'wrote tnumpart = ',tnumpart
      end if

      zoffset = zoffset + 4
     
      !Now compute the offset (in bytes!):
      offset = zoffset 
      do iproc = 0,myid-1
         offset = offset + pnum_vec(iproc+1)*(sizeof(tmp)-pad_diff) 
      end do

      !Now everyone else write, ONLY if numpart > 0
      if (numpart .GT. 0) then
      call mpi_file_write_at(fh,offset,writebuf,numpart,
     +                       particletype,istatus,ierr)
!      write(*,*) 'wrote writebuf(1)%pidx = ',writebuf(1)%pidx
!      write(*,'a30,3e15.6') 'wrote writebuf(1)%vp =',writebuf(1)%vp(1:3)
!      write(*,'a30,3e15.6') 'wrote writebuf(1)%xp =',writebuf(1)%xp(1:3)
!      write(*,*) 'wrote writebuf(2)%pidx = ',writebuf(2)%pidx
!      write(*,'a30,3e15.6') 'wrote writebuf(2)%vp =',writebuf(2)%vp(1:3)
!      write(*,'a30,3e15.6') 'wrote writebuf(2)%xp =',writebuf(2)%xp(1:3)
      end if

      call mpi_file_close(fh,ierr)

      write(*,*) 'proc',myid,'wrote numpart = ',numpart

      if (myid==0) write(*,7000) path_sav_part
 7000 format(' PARTICLE DATA IS WRITTEN IN FILE  ',a80)

      end subroutine save_particles
      subroutine read_part_res
      use pars
      use particles
      implicit none
      include 'mpif.h'

      integer :: istatus(mpi_status_size), ierr, fh
      integer(kind=mpi_offset_kind) :: zoffset,offset
      integer :: myp,totalp 
      integer :: iproc,i,pidxmax
      type(particle), allocatable :: readbuf(:)

      if (myid==0) write(*,7000) path_part 
 7000 format(' READING PARTICLE DATA FROM  ',a80)


      call mpi_file_open(mpi_comm_world,path_part,
     +                   mpi_mode_rdonly,
     +                   mpi_info_null,fh,ierr)


      !Read in the total number of particles:
      offset = 0
      call mpi_file_read_at_all(fh,offset,tnumpart,1,
     +                      mpi_integer,istatus,ierr)
      if (myid==0) write(*,*) 'read tnumpart = ',tnumpart
    
      offset = 4
      allocate(readbuf(tnumpart))
      call mpi_file_read_at_all(fh,offset,readbuf,tnumpart,
     +                       particletype,istatus,ierr)

      do i = 1,tnumpart
        !Now - does it lie within this proc's bounds?
        if (readbuf(i)%xp(2) .GT. ymin .AND.
     +       readbuf(i)%xp(2) .LT. ymax .AND.
     +       readbuf(i)%xp(1) .GT. xmin .AND.
     +       readbuf(i)%xp(1) .LT. xmax) then 
            if (.NOT. associated(first_particle)) then
               allocate(first_particle)
               first_particle = readbuf(i)
               nullify(first_particle%prev,first_particle%next)
               part => first_particle
            else
               allocate(part%next)
               part%next = readbuf(i)
               part%next%prev => part
               part => part%next
               nullify(part%next)
            end if

        end if
      end do

      deallocate(readbuf)

      call mpi_file_close(fh,ierr)
      
      !Now just check how many each processor obtained:
      !At the same time, figure out max(pidx) and set ngidx 
      !to one plus this value:
      pidxmax = 0
      part => first_particle
      myp = 0
      do while (associated(part))
         myp = myp+1
         if (part%pidx .gt. pidxmax) pidxmax = part%pidx
         part => part%next
      end do

      !Set ngidx (the index for creating new particles) to 1+pidmax:
      ngidx = pidxmax + 1

      numpart = myp
     
      call mpi_allreduce(myp,totalp,1,mpi_integer,mpi_sum,
     +                   mpi_comm_world,ierr)

      write(*,*) 'proc',myid,'read in numpart:',myp
      if (myid==0) write(*,*) 'total number of particles read:',totalp

      end subroutine read_part_res
      
c      subroutine particle_stats(ipt,jpt,kpt,iptc,jptc,kptc)
      subroutine particle_stats(ipt,jpt,kpt)
      use particles
      use pars
      implicit none
      integer :: i,ipt,jpt,kpt
c      integer :: iptc,jptc,kptc ! Indrajit-for SFS vel test
      real :: rhop,pi


      pi   = 4.0*atan(1.0)

      rhop = (m_s+4.0/3.0*pi*part%radius**3*rhow)/
     +                           (4.0/3.0*pi*part%radius**3)


      !Takes in ipt,jpt,kpt as the node to the "bottom left" of the particle
      !(i.e. the node in the negative direction for x,y,z)
      !and computes quantities needed to get particle statistics

      partcount_t(kpt,jpt,ipt) = partcount_t(kpt,jpt,ipt) + 1.0 !ORIGINAL

c     "Particle Sizers", determine the range
      if(part%radius .gt. 7.5e-6 .and. part%radius .lt. 15.5e-6)then
        pcountrange_t(kpt,jpt,ipt) = pcountrange_t(kpt,jpt,ipt) + 1.0
      end if

c      ----CHANGED FOR PARTICLE SFS VEL TEST RMOVE LATER INDRAJITH
c     partcount_t(kptc,jptc,iptc) = partcount_t(kptc,jptc,iptc) + 1.0
c     ----------------------------------------------------------
      
      !Get su mean, mean-squared of particle velocities at each level
c      upwp_t(kpt,jpt,ipt) = upwp_t(kpt,jpt,ipt) + part%vp(1)*part%vp(3)
c      do i = 1,3
c      vpsum_t(kpt,jpt,ipt,i) = vpsum_t(kpt,jpt,ipt,i) + part%vp(i)
c      vpsqrsum_t(kpt,jpt,ipt,i)=vpsqrsum_t(kpt,jpt,ipt,i)+part%vp(i)**2
c      end do

c      Tpsum_t(kpt,jpt,ipt) = Tpsum_t(kpt,jpt,ipt) + part%Tp
c      Tpsqrsum_t(kpt,jpt,ipt) = Tpsqrsum_t(kpt,jpt,ipt) + part%Tp**2

c      Tfsum_t(kpt,jpt,ipt) = Tfsum_t(kpt,jpt,ipt) + part%Tf

c      qfsum_t(kpt,jpt,ipt) = qfsum_t(kpt,jpt,ipt) + part%qinf

c      wpTpsum_t(kpt,jpt,ipt) = wpTpsum_t(kpt,jpt,ipt) + 
c     +                          part%Tp*part%vp(3)


c      radsum_t(kpt,jpt,ipt) = radsum_t(kpt,jpt,ipt) + part%radius  !Brian 7/30/14

c      rad2sum_t(kpt,jpt,ipt) = rad2sum_t(kpt,jpt,ipt) + part%radius**2  

c      mpsum_t(kpt,jpt,ipt) = mpsum_t(kpt,jpt,ipt) + 
c     +                rhop*4.0/3.0*pi*part%radius**3  !Brian 7/30/14

c      mwsum_t(kpt,jpt,ipt) = mwsum_t(kpt,jpt,ipt) + 
c     +                rhow*4.0/3.0*pi*part%radius**3  !Brian 7/30/14

c      qstarsum_t(kpt,jpt,ipt) = qstarsum_t(kpt,jpt,ipt) + part%qstar



      end subroutine particle_stats

      subroutine particle_stats_U(ipt,jpt,kpt)
      use particles
      use pars
      implicit none
      integer :: i,ipt,jpt,kpt
      real :: rhop, pi

      pi   = 4.0*atan(1.0)

      partcountU_t(kpt,jpt,ipt) = partcountU_t(kpt,jpt,ipt) + 1.0

      !Get su mean, mean-squared of particle velocities at each level
      upwp_t(kpt,jpt,ipt) = upwp_t(kpt,jpt,ipt) + part%vp(1)*part%vp(3)
      do i = 1,3
      vpsum_t(kpt,jpt,ipt,i) = vpsum_t(kpt,jpt,ipt,i) + part%vp(i)
      vpsqrsum_t(kpt,jpt,ipt,i)=vpsqrsum_t(kpt,jpt,ipt,i)+part%vp(i)**2
      end do

      Tpsum_t(kpt,jpt,ipt) = Tpsum_t(kpt,jpt,ipt) + part%Tp
      Tpsqrsum_t(kpt,jpt,ipt) = Tpsqrsum_t(kpt,jpt,ipt) + part%Tp**2

      Tfsum_t(kpt,jpt,ipt) = Tfsum_t(kpt,jpt,ipt) + part%Tf

      qfsum_t(kpt,jpt,ipt) = qfsum_t(kpt,jpt,ipt) + part%qinf

      wpTpsum_t(kpt,jpt,ipt) = wpTpsum_t(kpt,jpt,ipt) +
     +                          part%Tp*part%vp(3)

      radsum_t(kpt,jpt,ipt) = radsum_t(kpt,jpt,ipt) + part%radius  !Brian 7/30/14

      rad2sum_t(kpt,jpt,ipt) = rad2sum_t(kpt,jpt,ipt) + part%radius**2

      mpsum_t(kpt,jpt,ipt) = mpsum_t(kpt,jpt,ipt) +
     +                rhop*4.0/3.0*pi*part%radius**3  !Brian 7/30/14

      mwsum_t(kpt,jpt,ipt) = mwsum_t(kpt,jpt,ipt) +
     +                rhow*4.0/3.0*pi*part%radius**3  !Brian 7/30/14

      qstarsum_t(kpt,jpt,ipt) = qstarsum_t(kpt,jpt,ipt) + part%qstar


      end subroutine particle_stats_U

      subroutine particle_coupling_update
      use particles
      use pars
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'
      real :: wtx,wty,wtz,wtt,dV
      real :: rhop,taup_i,partmass
      real :: xv,yv,zv
      real :: ctbuf_s(nnz+2,1:iye-iys+2,6),cbbuf_r(nnz+2,1:iye-iys+2,6)
      real :: crbuf_s(nnz+2,1:mxe-mxs+1,6),clbuf_r(nnz+2,1:mxe-mxs+1,6)
      integer :: i,j,k,ncount,ipt,jpt,kpt,kwpt,kUpt
      integer :: istatus(mpi_status_size),ierr
      integer :: ix,iy,iz
c      integer :: iptc,jptc,kwptc
      real :: g(3)
      real :: lref,Cvv,Cva,Cpv
      real :: zu_1,zu_nnz

      integer :: iz_beg,iz_end

      g(1:3) = (/0.0, 0.0, part_grav/)

c     Reset all quantities for new readings
      partsrc_t = 0.0
      partTsrc_t = 0.0
      partHsrc_t = 0.0
      partTEsrc_t = 0.0

      part => first_particle
      do while (associated(part))
      !First, as done in uf_interp, must find the "leftmost" node
      !of volume where particle belongs:
      !(must repeat since now particle locations have been updated)
     
        ipt = floor(part%xp(1)/dx) + 1
        jpt = floor(part%xp(2)/dy) + 1

!       Some particles may have died.. account only for those still
!       living
        if (part%xp(3) .gt. 0.0) then
          kUpt = floor((part%xp(3)-0.0)/dz_U) + 1
          call particle_stats_U(ipt,jpt,kUpt)


          do iz=0,nnz+1
            if (zz(iz) .GT. part%xp(3)) go to 1299
          end do
1299  continue
              kpt = iz-1

              !Do again for w:
              do iz = 0,nnz+1
                if (z(iz) .GT. part%xp(3)) go to 1120
              end do
1120  continue

                kwpt = iz

      !As an aside, use ipt,jpt,kwpt to update the particle numbers at each cell:

      !call inputs are the nodes associated with the volume that each part. lies in
c      call particle_stats(ipt,jpt,kwpt,iptc,jptc,kwptc)
                call particle_stats(ipt,jpt,kwpt)  
!---------

      !Add contribution to each of the 8 surrounding nodes:
      do i=0,1
      do j=0,1
      do k=0,1

      xv = dx*(i+ipt-1)
      yv = dy*(j+jpt-1)
      zv = zz(k+kpt)

      dV = dx*dy*dzu(kpt+1)

      wtx = (1.0 - abs(part%xp(1)-xv)/dx)
      wty = (1.0 - abs(part%xp(2)-yv)/dy)
      wtz = (1.0 - abs(part%xp(3)-zv)/dzu(kpt+1))
      wtt = wtx*wty*wtz

      rhop = (m_s+pi2*2.0/3.0*part%radius**3*rhow)/
     +                           (pi2*2.0/3.0*part%radius**3)
      partmass = rhop*2.0/3.0*pi2*(part%radius)**3
      taup_i = 18.0*rhoa*nuf/rhop/(2.0*part%radius)**2 !Brian 8/4/14

      !write(*,'a6,3i') 'node:',ipt+i,jpt+j,kpt+k
      !write(*,'a6,5e15.6') 'Vwt: ',dV,wtx,wty,wtz,wtt

      ix = ipt+i
!      !Account for periodicity:
!      if (ix == nnx+1) ix = 1
      iy = jpt+j
      iz = kpt+k

      !This error is saying that a particle is outside the bounds
      !of a CPU
      if (ix .gt. mxe+1) write(*,*) 'proc',myid,'has ix = ',ix,mxe+1
      if (ix .lt. mxs) write(*,*) 'proc',myid,'has ix = ',ix,mxs
      if (iy .gt. iye+1) write(*,*) 'proc',myid,'has iy = ',iy,iye+1
      if (iy .lt. iys) write(*,*) 'proc',myid,'has iy = ',iy,iys
      if (iz .gt. nnz+1) write(*,*) 'proc',myid,'has iz = ',iz,nnz+1
      if (iz .lt. 0) write(*,*) 'proc',myid,'has iz = ',iz

      !Recall to subtract minus g (add g) since momentum is extracted form
      !fluid only through drag term - NOT the gravity term as well
      if (icouple == 1) then
      !drag momentum coupling
        partsrc_t(iz,iy,ix,1:3) =
     +     partsrc_t(iz,iy,ix,1:3) +
     +     partmass/rhoa*(part%vrhs(1:3)+g(1:3))*wtt/dV

      !vapor momentum coupling
        partsrc_t(iz,iy,ix,1:3) =
     +     partsrc_t(iz,iy,ix,1:3) +
     +     rhow/rhoa*pi2*2*part%radius**2*
     +     part%radrhs*part%vp(1:3)*wtt/dV
      endif

      if (iTcouple == 1) then
        partTsrc_t(iz,iy,ix) =
     +     partTsrc_t(iz,iy,ix) +
     +     (part%Tprhs_s*6.0*rhow/rhop/CpaCpp/taup_i*
     +     (pi2/2.0)*part%radius*nuf)*wtt/dV
      endif

      if (iHcouple == 1) then
        partHsrc_t(iz,iy,ix) =     !Brian 8/5/14
     +     partHsrc_t(iz,iy,ix) +  !minus comes later
     +     rhow/rhoa*pi2*2*part%radius**2*part%radrhs*wtt/dV



!      http://eesc.columbia.edu/courses/ees/climate/SI_Units.pdf
!      or http://catalog.conveyorspneumatic.com/Asset/FLS%20Specific%20Heat%20Capacities%20of%20Gases.pdf
         Cvv = 1463.0 !J/degkg
         Cpv = 1952.0 !from E. Toolbox -> 1864.0 !J/degkg
       !Cva = 717.0 !J/degkg
       !lref = (25.0 - 0.02274*20.0)*10.0**5
         partTEsrc_t(iz,iy,ix) = 
     +     partTEsrc_t(iz,iy,ix) +  !minus comes later
     +         rhow/rhoa*pi2*2*part%radius**2*
     +                  part%radrhs*Cpv/Cpa*part%Tp*wtt/dV -
     +         rhow/rhoa*pi2*2*part%radius**2*
     +                  part%radrhs*Cpv/Cpa*part%Tf*wtt/dV

      endif

      !write(*,'a5,3i,3e15.6') 'bsrc:',ix,iy,iz,partsrc(ix,iy,iz,1:3)

      end do
      end do
      end do

      end if !Particle within domain?

      part => part%next
      end do

      !Now, partsrc and partTsrc have halos on each processor - give these to the rightful owner:
      crbuf_s=0.0;ctbuf_s=0.0
      clbuf_r=0.0;cbbuf_r=0.0
    
      !First send top: 
      !get send buffer ready:
      ctbuf_s(1:nnz+2,1:iye-iys+2,1:3)=partsrc_t(0:nnz+1,
     +                                iys:iye+1,mxe+1,1:3)
      ctbuf_s(1:nnz+2,1:iye-iys+2,4)=partTsrc_t(0:nnz+1,
     +                                iys:iye+1,mxe+1)

      ctbuf_s(1:nnz+2,1:iye-iys+2,5)=partHsrc_t(0:nnz+1,
     +                                iys:iye+1,mxe+1)

       ctbuf_s(1:nnz+2,1:iye-iys+2,6)=partTEsrc_t(0:nnz+1,
     +                                iys:iye+1,mxe+1)

      ncount = 6*(nnz+2)*(iye-iys+2) 
      call mpi_sendrecv(ctbuf_s,ncount,mpi_real8,tproc,1,
     +     cbbuf_r,ncount,mpi_real8,bproc,1,mpi_comm_world,istatus,ierr)

      !Now just add the contents of the receive buffer into the entire 
      !iys column of this proc:

      partsrc_t(0:nnz+1,iys:iye+1,mxs,1:3) = 
     +   partsrc_t(0:nnz+1,iys:iye+1,mxs,1:3) + 
     +    cbbuf_r(1:nnz+2,1:iye-iys+2,1:3)
      partTsrc_t(0:nnz+1,iys:iye+1,mxs) = 
     +   partTsrc_t(0:nnz+1,iys:iye+1,mxs) + 
     +     cbbuf_r(1:nnz+2,1:iye-iys+2,4)

      partHsrc_t(0:nnz+1,iys:iye+1,mxs) = 
     +   partHsrc_t(0:nnz+1,iys:iye+1,mxs) + 
     +     cbbuf_r(1:nnz+2,1:iye-iys+2,5)
      
    
      partTEsrc_t(0:nnz+1,iys:iye+1,mxs) = 
     +   partTEsrc_t(0:nnz+1,iys:iye+1,mxs) + 
     +   cbbuf_r(1:nnz+2,1:iye-iys+2,6)
      
      !Now get the right send buffer ready:
      crbuf_s(1:nnz+2,1:mxe-mxs+1,1:3)=
     +        partsrc_t(0:nnz+1,iye+1,mxs:mxe,1:3)
      crbuf_s(1:nnz+2,1:mxe-mxs+1,4)=
     +        partTsrc_t(0:nnz+1,iye+1,mxs:mxe)
      crbuf_s(1:nnz+2,1:mxe-mxs+1,5)=
     +        partHsrc_t(0:nnz+1,iye+1,mxs:mxe)

      crbuf_s(1:nnz+2,1:mxe-mxs+1,6)=
     +        partTEsrc_t(0:nnz+1,iye+1,mxs:mxe)

      !Now send to right:
      ncount = 6*(nnz+2)*(mxe-mxs+1)
      call mpi_sendrecv(crbuf_s,ncount,mpi_real8,rproc,2,
     +     clbuf_r,ncount,mpi_real8,lproc,2,mpi_comm_world,istatus,ierr)

      !And again add the contents to the top/bottom rows of partsrc:
      partsrc_t(0:nnz+1,iys,mxs:mxe,1:3) = 
     +           partsrc_t(0:nnz+1,iys,mxs:mxe,1:3) + 
     +                clbuf_r(1:nnz+2,1:mxe-mxs+1,1:3)

      partTsrc_t(0:nnz+1,iys,mxs:mxe) = 
     +           partTsrc_t(0:nnz+1,iys,mxs:mxe) + 
     +                 clbuf_r(1:nnz+2,1:mxe-mxs+1,4)
      partHsrc_t(0:nnz+1,iys,mxs:mxe) = 
     +           partHsrc_t(0:nnz+1,iys,mxs:mxe) + 
     +                 clbuf_r(1:nnz+2,1:mxe-mxs+1,5)  
      partTEsrc_t(0:nnz+1,iys,mxs:mxe) = 
     +           partTEsrc_t(0:nnz+1,iys,mxs:mxe) + 
     +                 clbuf_r(1:nnz+2,1:mxe-mxs+1,6)  

      !Now everyone should have the proper summed contributions from particles
      !lying within (i-1,i+1),(j-1,j+1),(k-1,k+1)

!      do ix=1,nnx
!      do iy=iys,iye
!      do iz=izs,ize
!
!      if ( (abs(partsrc(ix,iy,iz,1)) .gt. 1.0e-8) .or.
!     +      (abs(partsrc(ix,iy,iz,2)) .gt. 1.0e-8) .or. 
!     +      (abs(partsrc(ix,iy,iz,3)) .gt. 1.0e-8)) then
!
!      write(*,'a5,4i,3e15.6')'asrc:',myid,ix,iy,iz,partsrc(ix,iy,iz,1:3)
!      end if
!
!      end do
!      end do
!      end do

      end subroutine particle_coupling_update

      subroutine read_input_file
      use pars
      use particles
      use con_data
      implicit none

      character(48) :: label
      character(180) :: params_dir
      namelist /step_params/ iti,itmax,imean,ihst,itape,iphist,
     +                       itstr,it_his,it_viz,i_viz,itn,dt_con,
     +                       it_probe

      namelist /grid_params/ ncpu_s, Uo, Ttop, Tbot,
     +         qstar, tsfcc,
     +         zi, zl, xl, yl, zw1,zwpi

      namelist /path_names/ path_seed,path_part,path_res,
     +         path_sav,path_his,path_viz_xy,path_viz_xz,
     +         path_viz_yz,path_stuf,path_ran,part_hist,
     +         path_probe

      namelist /flags/ ismlt,ifree,isfc,iradup,
     +         iupwnd,ibuoy,ifilt,itcut,isubs,ibrcl,iocean,
     +         method,idebug,iz_space,ivis0,ifix_dt,new_vis,iDNS,
     +         ispray,icouple,iTcouple,iHcouple,ievap,istbl,isfs,
     +         iUprod,ipolydisp,ishift,ifields,iPprod,
     +         iSSGF, inonperiodic,
     +         iaircraft !Aircraft probe


      namelist /constants/ rhoa, nuf, Cpa, Pra, Sc,
     +         tnumpart,tnumpart_max,Nrp,numpart_dt,
     +         part_init,rhow,part_grav,
     +         Cpp,Mw,Ms,Ru,Gam,Ion,Os,Sal,
     +         radius_max,radius_min,delta_vis,
     +         grav, t00,fcor,
     +         dtdzf,zo,ugcont,vgcont,
     +         cond_shift_init, ugcont_new,
     +         vgcont_new, qstar_new,
     +         vp_init,Tp_init,qf_init,
     +         apr_v,apr_rad,apr_init

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

      dz_U = (zl-0.0)/maxnz



      end subroutine read_input_file

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

      subroutine particle_eq_Temp       ! Indrajith 02/07/2017
      use pars
      use particles
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'
      real :: esp_abs, esp_step
      real :: Tp1, Tp2, Tp3, FT1,FT2,FT3, T_eqb
      integer :: steps
      real :: FT
      ! Calculate Equilibrium temperature for droplet radius evolution 
      esp_abs = 10**-5
      esp_step = 10**-5
      T_eqb = part%Tf
      Tp1 = T_eqb - 10
      Tp2 = T_eqb
      FT1 = FT(Tp1)
      FT2 = FT(Tp2)
      steps = 0
      !write(*,*) 'Teq-calculation',Tp1,Tp2,part%Tf
      ! do while((Tp2-Tp1).ge.esp_step.or.abs(FT1).ge.esp_abs.and.
!      +       abs(FT2).ge.esp_abs.and.steps.le.10**4)
      do while((Tp2-Tp1).ge.esp_step.and.steps.le.100)
        Tp3 = (Tp1+Tp2)/2
        FT3 = FT(Tp3)
        if (FT3.eq.0.0)then
          exit
        else if (FT1*FT3.lt.0.0) then
          Tp2 = Tp3
        else
          Tp1 =Tp3
        end if
        FT1 = FT(Tp1)
        FT2 = FT(Tp2)
        steps = steps + 1
      end do
      T_eqb = (Tp1+Tp2)/2
      part%T_eq = T_eqb
      end subroutine particle_eq_Temp

      function FT(Tp)         ! Indrajith 02/07/2017
      use pars
      use particles
      use con_data
      use con_stats
      implicit none
      real :: Rep, diff(3),diffnorm
      real :: Nup, Shp, rhop, taup_i, estar, einf
      real :: Tfc, TpC, Lv, Eff_C, Eff_S, Volp
      real :: Temp_eq, FTtemp, Tp, FT
      ! This calculate the value of the function used to define the particle temperature
      diff(1:3) = part%vp - part%uf
      diffnorm = sqrt(diff(1)**2 + diff(2)**2 + diff(3)**2)
      Rep = 2.0*part%radius*diffnorm/nuf  
      Volp = pi2*2.0/3.0*part%radius**3
      rhop = (m_s+Volp*rhow)/Volp
      taup_i = 18.0*rhoa*nuf/rhop/(2.0*part%radius)**2 
         !Compute Nusselt number for particle:
         !Ranz-Marshall relation
      Nup = 2.0 + 0.6*Rep**(1.0/2.0)*Pra**(1.0/3.0)
      Shp = 2.0 + 0.6*Rep**(1.0/2.0)*Sc**(1.0/3.0)

      TfC = part%Tf-273.15
      einf = 610.94*EXP(17.6257*TfC/(TfC+243.04))
      TpC = part%Tp-273.15
      Lv =(25.0 - 0.02274*26.0)*10.0**5

      Eff_C = 2.0*Mw*Gam/(Ru*rhow*part%radius*part%Tp)
      Eff_S = Ion*Os*m_s*Mw/Ms/(Volp*rhop-m_s)
      estar = einf*exp(Mw*Lv/Ru*(1.0/part%Tf-1.0/part%Tp)
     +                                      +Eff_C-Eff_S)
      part%qstar = Mw/Ru*estar/part%Tp/rhoa
      Temp_eq = Shp/Sc*Pra/Nup* Lv/Cpa*(part%qinf-part%qstar)        !assumes qinf=rhov/rhoa rather than rhov/rhom
      FTtemp = part%Tf + Temp_eq
      FT = part%Tp -FTtemp 
      return
      end function FT


      subroutine particle_eq_radius       ! Indrajith 02/07/2017
      use pars
      use particles
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'
      real :: esp_abs, esp_step
      real :: rp1, rp2, rp3, FR1,FR2,FR3, R_eqb
      integer :: steps
      real :: FR
      ! Calculate Equilibrium temperature for droplet radius evolution 
      esp_abs  = 1e-7
      esp_step = 1e-7
      R_eqb    = part%radius
      rp1      = R_eqb
      rp2      = R_eqb/2
      FR1      = FR(rp1)
      FR2      = FR(rp2)
      steps    = 0
      do while(FR1*FR2.ge.0.0)
        rp2 = rp2/2
        FR2 = FR(rp2)
      end do
        
        ! do while((Tp2-Tp1).ge.esp_step.or.abs(FT1).ge.esp_abs.and.
!        +       abs(FT2).ge.esp_abs.and.steps.le.10**4)
       
      do while(ABS(rp1-rp2).gt.esp_step.and.steps.le.1000)
        rp3 = (rp1+rp2)/2
        FR3 = FR(rp3)
        if(FR3.eq.0.0)then
          rp1 = rp3
          rp2 = rp3
          exit
        else if (FR1*FR3.lt.0.0) then
          rp2 = rp3
        else
          rp1 =rp3
        end if
          FR1 = FR(rp1)
          FR2 = FR(rp2)
          steps = steps + 1
          if(steps.gt.1000)then
            write(*,*)'R_eq did not converge'
          end if
      end do
      R_eqb     = (rp1+rp2)/2
      part%r_eq = R_eqb
      !write(*,*) steps,(rp1-rp2),part%radius,R_eqb
      end subroutine particle_eq_radius

      function FR(rp)         ! Indrajith 02/25/2017
      use pars
      use particles
      use con_data
      use con_stats
      implicit none
      real :: rhop, einf
      real :: TfC, RHa, Lv, Eff_C, Eff_S, Volp
      real :: rp, q_inf, FR
      ! This calculate the value of the function used to define the particle equilibrium  radius                 
      Volp = pi2*2.0/3.0*rp**3
      rhop = (m_s+Volp*rhow)/Volp
      TfC = part%Tf-273.15
      einf = 610.94*EXP(17.6257*TfC/(TfC+243.04))
      q_inf = Mw/Ru*einf/part%Tf/rhoa
      RHa = part%qinf/q_inf
      Lv =(25.0 - 0.02274*26.0)*10.0**5
      Eff_C = 2.0*Mw*Gam/(Ru*rhow*rp*part%Tf)
      !Eff_S = Ion*Os*m_s*Mw/Ms/(Volp*rhop-m_s)
      Eff_S = Ion*Os*(Mw/Ms)*(Sal/1000.0)*part%radius**3/(rp**3) 
      FR =  Eff_S - Eff_C + log(RHa)
      return
      end function FR

      subroutine particle_histogram(istep)
      use pars
      use fields
      use con_data
      use con_stats
      use particles
      implicit none
      include 'mpif.h'
      integer:: i
      integer :: ierr
      integer,intent(in) :: istep
      integer :: iblnk 
      integer :: sumbuf_radclass(Nrp+2)
      integer :: sumbuf_iniradclass(Nrp+2)
      integer :: sumbuf_delradclass(Nrp+2)
      integer :: sumbuf_req_init(Nrp+2)
      real :: binval_radclass(Nrp+2)
       
      character*80 path_histogram
      character cgrid*10, num*3, iblks*16

      sumbuf_radclass    = 0
      sumbuf_iniradclass = 0
      sumbuf_delradclass = 0
      sumbuf_req_init    = 0
      binval_radclass    = bin_radclass

c       sum over each procesor
      call mpi_reduce(hist_radclass(:),sumbuf_radclass,Nrp+2,
     +     mpi_integer8,mpi_sum,0,mpi_comm_world,ierr)

      call mpi_reduce(init_radclass(:),sumbuf_iniradclass,Nrp+2,
     +     mpi_integer8,mpi_sum,0,mpi_comm_world,ierr)

      call mpi_reduce(del_radclass(:),sumbuf_delradclass,Nrp+2,
     +     mpi_integer8,mpi_sum,0,mpi_comm_world,ierr)

      call mpi_reduce(hist_req_init(:),sumbuf_req_init,Nrp+2,
     +     mpi_integer8,mpi_sum,0,mpi_comm_world,ierr)



c     build the character string for output file name
      cgrid = '.le.'
      write(num,'(i3.3)') itn
      iblnk = index(part_hist,' ')
      call blnk(iblks)
      write(iblks(1:7),'(i7.7)') istep
      iblks(8:8) = '_'
      write(iblks(9:15),'(i7.7)') (istep+iphist)
      call blnk(path_sav_hist)  
      
c     path_histogram = part_hist(1:iblnk-1)//'/particle'//
c     +               cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'
      path_histogram = part_hist(1:iblnk-1)//'/particle'//
     +                   '.'//iblks(1:15)//'.txt'

      if (myid==0)then
        open(unit =900, file = path_histogram,form ='formatted',
     +       access = 'sequential')
        do i = 1,Nrp+2
          write(900,'(e16.7,I8,I8,I8,I8)') bin_radclass(i),
     +          sumbuf_radclass(i),sumbuf_iniradclass(i),
     +          sumbuf_delradclass(i),sumbuf_req_init(i)
c     Writing histogram statistics on the front screen
          write(*,*)'Histogram: ',bin_radclass(i), sumbuf_radclass(i)
     +                  ,sumbuf_iniradclass(i),sumbuf_delradclass(i)
        end do
        close(900)         
      end if
      hist_radclass = 0
      hist_req_init = 0
      end subroutine particle_histogram
        

      subroutine particle_sizeclass
      use pars
      use fields
      use con_data
      use con_stats
      use particles
      implicit none
      integer :: i_class,i
      ! calculate the no of droplets in each calss
      !if (myid .eq.0)then
      !   do i = 2,Nrp+2
      !        write(*,*)'sizeclass',class_minr(i),class_maxr(i)
      !   end do
      !end if
      part => first_particle
      do while (associated(part))
c     particle radius is less than zero
        if(part%radius .lt. 0.0)then
          write(*,*) 'Error_particle_radius', part%radius
          stop
        end if

c       particle radius is larger than zero
c       particle radius is less than the smallest class

        if(part%radius .lt. class_maxr(1))then
          hist_radclass(1) = hist_radclass(1) + 1
        end if

c       particle radius is larger than largest class
!       JP 9/5/2018 maybe record this in another bin as well?
        if(part%radius .ge. class_maxr(Nrp+2))then
          write(*,*)'Large_drop. radius:',part%radius
        end if

c       Particle radius is with in the range

        do i_class = 2,Nrp+2
          if (part%radius .ge.class_minr(i_class).and.
     +        part%radius .lt.class_maxr(i_class)) then
            hist_radclass(i_class) = hist_radclass(i_class) + 1
          end if
        end do

        part=>part%next
      end do
      end subroutine particle_sizeclass

      subroutine req_sizeclass
c     this subroutine is only for checking the equilibrium radius distribution for initial radius 
      use pars
      use fields
      use con_data
      use con_stats
      use particles
      implicit none
      integer :: i_class,i
      ! calculate the no of droplets in each class
      !if (myid .eq.0)then
      !   do i = 2,Nrp+2
      !        write(*,*)'sizeclass',class_minr(i),class_maxr(i)
      !   end do
      !end if
      part => first_particle
      do while (associated(part))   
c       particle eq_radius is less than zero
        if(part%req_init .lt. 0.0)then
          write(*,*) 'Error_req_init', part%req_init
          stop
        end if
c       particle eq_radius is larger than zero
c       particle eq_radius is less than the smallest class
        if(part%req_init.lt.class_maxr(1))then
          hist_req_init(1) = hist_req_init(1) + 1
        end if

c       particle eq_radius is larger than largest class
!       JP 9/5/2018 maybe record this in another bin as well?
        if (part%req_init .ge.class_maxr(Nrp+2))then
          write(*,*)'Large_drop for req_init',part%req_init
        end if
c       Particle radius is with in the range

        do i_class = 2,Nrp+2
          if (part%req_init .ge.class_minr(i_class).and.
     +        part%req_init .lt.class_maxr(i_class)) then
            hist_req_init(i_class) = hist_req_init(i_class) + 1
          end if
        end do
        part=>part%next
      end do
      end subroutine req_sizeclass
        

        subroutine matrix_inverse(A,B)
c       This calculate inverse matrix of a 3x3 matrix
c	This requires in SFS velocity calculation for
c	Non isotropic case
        real :: A(3,3)
        real :: B(3,3)
        real :: det_i,det
        
        det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2))-
     +          A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1))+
     +          A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))

        det_i = 1/det        
        B(1,1) = det_i*(A(2,2)*A(3,3) - A(2,3)*A(3,2))
        B(1,2) = -det_i*(A(2,1)*A(3,3) - A(2,3)*A(3,1))
        B(1,3) = det_i*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
        B(2,1) = -det_i*(A(1,2)*A(3,3) - A(3,2)*A(3,1))
        B(2,2) = det_i*(A(1,1)*A(3,3) - A(1,3)*A(3,1))
        B(2,3) = -det_i*(A(1,1)*A(3,2) - A(1,2)*A(3,1))
        B(3,1) = det_i*(A(1,2)*A(2,3) - A(1,3)*A(2,2))
        B(3,2) = -det_i*(A(1,1)*A(2,3) - A(1,3)*A(2,1))
        B(3,3) = det_i*(A(1,1)*A(2,2) - A(1,2)*A(2,1))
        end subroutine matrix_inverse

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
