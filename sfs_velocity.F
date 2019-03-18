

      subroutine SFS_velocity
c     This subroutine calculate the SFS velocity for spray drops
c     This uses the Weil et al. (2004) formulation
c     Calculations assume isotropic turbulence
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
      implicit none
      include 'mpif.h'
      real :: sigm_sdxp,sigm_sdyp,sigm_sdzp,vis_sp
      real :: sigm_su,sigm_sl,us_ran,gasdev,tengz,englez_bar
      real :: engsbz_bar,sigm_w, sigm_ws
      real :: L_flt,epsn,fs,C0,a1,a2,a3,sigm_sprev,fs1
      real :: weit,weit1,weit3,weit4, T_lagr
      real :: us(3)
      real :: u_a1,u_a2,u_a3
      real :: v_a1,v_a2,v_a3
      real :: w_a1,w_a2,w_a3
      integer :: ix,iy,iz,izp1,izm1,ind,iz_part
      integer :: ierr
c     ---initialize (and reset all to zero)-------
      fs         = 0.0
      C0         = 0.0
      T_lagr     = 0.0
      l_flt      = 0.0
      epsn       = 0.0
      tengz      = 0.0
      englez_bar = 0.0
      engsbz_bar = 0.0
      sigm_s     = 0.0 
      sigm_sdx   = 0.0
      sigm_sdy   = 0.0
      sigm_sdz   = 0.0
      sigm_su    = 0.0
      sigm_sl    = 0.0
      us_ran     = 0.0
      us         = 0.0
      sigm_ws    = 0.0
      sigm_w     = 0.0
c	------------------
c       compute sigma squre (sigm_s) based on subgrid energy field
c	-----------------       
      do iz = izs,ize
        izp1  = iz+1
        izm1  = iz-1
        weit  = dzw(iz)/(dzw(iz)+dzw(izp1))
        weit1 = 1-weit
        weit3 = dzw(izm1)/(dzw(iz)+dzw(izm1))
        weit4 = 1-weit3
        do ix =1,nnx
          do iy = iys,iye
            sigm_s(ix,iy,iz)   = 2.0*e(ix,iy,iz)/3.0
            vis_ss(ix,iy,iz)   = vis_s(ix,iy,3,iz)
            sigm_sdx(ix,iy,iz) =  sigm_s(ix,iy,iz)       !for xderiv
            sigm_sdy(ix,iy,iz) =  sigm_s(ix,iy,iz)       !for yderiv
c	--------------------
c       calculate z derivative of sigma_s
c       this will be at the u-point!
c	--------------------
            sigm_sdz(ix,iy,iz)=
     +      (sigm_s(ix,iy,iz)-sigm_s(ix,iy,izm1))*dzw_i(iz)
          end do
        end do

c	-------------------
c       calculate x derivatives of sigma_s
c	-------------------
        call xderivp(sigm_sdx(1,iys,iz),trigx(1,1),xk(1),nnx,iys,iye) 

      end do
        
c	------------------
c       calculate y derivative of sigma_s
c	------------------
      call yd_mpi(sigm_sdy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c       ----------------
c      calculate extented fileds of sigm_s and its derivatives
c       ---------------
      call fill_extSFS
        !Loop over the linked list of particles:
      part => first_particle
      do while (associated(part))
      ! interpolate sigm_s and its derivative at particle location
        sigm_sprev  = part%sigm_s
        call sigm_interp(sigm_sdxp,sigm_sdyp,sigm_sdzp,vis_sp,iz_part)
        part%sigm_s = abs(part%sigm_s)  !Interpolation near surface can give small negative numbers
c	-----------------
c       calculate the subgrid velocity Weil et al ,2004-isotropic turb.
c	----------------
        l_flt = (2.25*dx*dy*dzw(iz_part+1))**(1.0/3.0)! filtered with
        epsn  = (0.93/l_flt)*(3.0*part%sigm_s/2.0)**(1.5) ! tur. dis.rt
        ! TKE : resolved + Subgrid at grid center
        !---------------------------------
c       tot_eng = (engsbz(iz_part+1) + engz(iz_part+1))
        englez_bar = 0.5*(englez(iz_part)+englez(iz_part+1))
        engsbz_bar = 0.5*(engsbz(iz_part)+engsbz(iz_part+1))
        tengz      = englez_bar + engsbz_bar
        !---------------------------------
        ! Calculate fs basd on w-componet of velocity
        sigm_w  = 0.5*(wps(iz_part)+wps(iz_part+1))
        sigm_ws = 0.5*(engsbz(iz_part)+engsbz(iz_part+1))/3.0

c       ---------write for single droplet ---------
c         write(*,*)'sigm_w:', sigm_w,sigm_ws
c       ------------------------------------------

        if(tengz.gt.0.0)then
c          fs = engsbz(iz_part+1)/(engsbz(iz_part+1) + engz(iz_part+1))
c          fs = engsbz_bar/tengz
          fs = sigm_ws/(sigm_w + sigm_ws)
        else
          fs = 0.0
        end if
        C0 = 6.0  ! Changed to 6.0 from 3.0 Indrjith 11-20-17

c       ---------Check for single-part----------
        T_lagr = 2*part%sigm_s/(C0*epsn)   ! Lagrangian time scale
c         write(*,*) 'L_time:',part%xp(3),T_lagr
c       -----------------------------------
c      ------------------
c        Calculate subgrid velocity components
c       -----------------
        us(1:3) = part%u_sub(1:3)
c      -----x component ----------------
        u_a1 = 0.0
        u_a2 = 0.0
        u_a3 = 0.0
        u_a1 = (-0.5)*fs*C0*epsn*part%u_sub(1)/part%sigm_s
        u_a2 = 0.5*(1.0/part%sigm_s)*part%u_sub(1)*
     +          (part%sigm_s - sigm_sprev)/dt
        u_a3 = 0.5*sigm_sdxp
        us_ran        = sqrt(fs*C0*epsn*dt)*gasdev(iseed)
        part%u_sub(1) = (u_a1+u_a2+u_a3)*dt + us_ran

c      ------FOR CHECK-INDRAJIT-11/17/2017-----
c        if(fs.gt.0.5)then
c        part%u_sub(1) = gasdev(iseed)*sqrt(2*vis_sp/dt)
c        end if
c      -----------------------------------------
        ! Check to see it blows up 
        if(isnan(part%u_sub(1)))then
          write(*,*)'u-nan',u_a1,u_a2,u_a3,us_ran,part%sigm_s
          write(*,*)'part-con',fs,epsn,dzw(iz_part+1),iz_part
          write(*,*)'part-nan',part%xp(1),part%xp(2),part%xp(3)
          write(*,*)'energy',engsbz(iz_part+1),engz(iz_part+1)
          write(*,*)'derev-nan',sigm_sdxp,sigm_sdyp,sigm_sdzp
          stop
        end if
c       ----------------------------------
c       -----y component ----------------
        v_a1 = 0.0
        v_a2 = 0.0
        v_a3 = 0.0
        v_a1 = (-0.5)*fs*C0*epsn*part%u_sub(2)/part%sigm_s
        v_a2 = 0.5*(1.0/part%sigm_s)*part%u_sub(2)*
     +          (part%sigm_s - sigm_sprev)/dt
        v_a3 = 0.5*sigm_sdyp
        us_ran        = sqrt(fs*C0*epsn*dt)*gasdev(iseed)
        part%u_sub(2) = (v_a1+v_a2+v_a3)*dt + us_ran

c      ------FOR CHECK-INDRAJIT-11/17/2017-----
c        if(fs.gt.0.5)then
c        part%u_sub(2) = gasdev(iseed)*sqrt(2*vis_sp/dt)
c        end if
c      -----------------------------------------

        if(isnan(part%u_sub(2)))then
          write(*,*)'v-nan',v_a1,v_a2,v_a3,us_ran,part%sigm_s
          write(*,*)'part-nan',part%xp(1),part%xp(2),part%xp(3)
          write(*,*)'derev-nan',sigm_sdxp,sigm_sdyp,sigm_sdzp
          stop
        end if

c       ----------------------------------
c       -----z component ----------------
        a1 = 0.0
        a2 = 0.0
        a3 = 0.0
        a1 = (-0.5)*fs*C0*epsn*part%u_sub(3)/part%sigm_s
        a2 = 0.5*(1.0/part%sigm_s)*part%u_sub(3)*
     +          (part%sigm_s - sigm_sprev)/dt
        a3 = 0.5*sigm_sdzp
        us_ran        = sqrt(fs*C0*epsn*dt)*gasdev(iseed)
        part%u_sub(3) = (w_a1+w_a2+w_a3)*dt + us_ran

!        if (myid==0 .and. part%xp(3) .lt. 8.0) then
!        write(*,'(a3,5e15.6,i)') 'T:',a3,us_ran,us_ran/a3,part%u_sub(3),
!     +  part%xp(3),iz_part
!        end if
        
c       ------- increase the u_sub close to surface -------
c        if(fs.gt.0.5)then
c        us_ran = sqrt(fs*C0*epsn*dt)*gasdev(iseed)*2.0
c        part%u_sub(3) = gasdev(iseed)*sqrt(2*vis_sp/dt)
c        end if
c       -------------------------------------------------
c        part%u_sub(3)= (a2+a3)*dt + us_ran
c         write(*,*)'fs1:',part%xp(3),fs,us_ran
c        write(*,*)'as1:',a1,a2,a3,us_ran
c      ------FOR CHECK-INDRAJIT-11/17/2017-----
c        if(part%xp(3).lt. 0.1*zi.and.fs.gt.0.5)then
c        if(part%xp(3).lt. zw1)then
c        part%u_sub(3) = gasdev(iseed)*sqrt(2*vis_sp/dt)
c        part%u_sub(3)= part%u_sub(3) + gasdev(iseed)
c        part%u_sub(3)= part%u_sub(3)+9.0*us_ran
c        end if
c      -----------------------------------------

        if(isnan(part%u_sub(3)))then
          write(*,*)'w-nan',w_a1,w_a2,w_a3,us_ran,part%sigm_s
          write(*,*)'part-nan',part%xp(1),part%xp(2),part%xp(3)
          write(*,*)'derev-nan',sigm_sdxp,sigm_sdyp,sigm_sdzp
          stop
        end if

c       some more checks because how the hell
c       are particles being teleported so far..
        if(part%u_sub(1).gt. 100.0 .or. part%u_sub(1) .lt. -100.0)then
          write(*,*)'u_sub(1) extrema: ',part%u_sub(1)
          write(*,*)'part loc: ',part%xp(1),part%xp(2),part%xp(3)
          write(*,*)'a1, a2, a3: ',u_a1,u_a2,u_a3
          write(*,*)'sigm_s: ',part%sigm_s
          write(*,*)'sigm_sprev: ',sigm_sprev
        end if
        if(part%u_sub(2).gt. 100.0 .or. part%u_sub(2) .lt. -100.0)then
          write(*,*)'u_sub(2) extrema: ',part%u_sub(2)
          write(*,*)'part loc: ',part%xp(1),part%xp(2),part%xp(3)
          write(*,*)'a1, a2, a3: ',v_a1,v_a2,v_a3
          write(*,*)'sigm_s: ',part%sigm_s
          write(*,*)'sigm_sprev: ',sigm_sprev
        end if
        if(part%u_sub(3).gt. 100.0 .or. part%u_sub(3) .lt. -100.0)then
          write(*,*)'u_sub(3) extrema: ',part%u_sub(3)
          write(*,*)'part loc: ',part%xp(1),part%xp(2),part%xp(3)
          write(*,*)'a1, a2, a3: ',w_a1,w_a2,w_a3
          write(*,*)'sigm_s: ',part%sigm_s
          write(*,*)'sigm_sprev: ',sigm_sprev
        end if








c       -------------------
c       Update particle location and velocity
c       ------------------
        do ind = 1,3
          part%xp(ind)  = part%xp(ind) + part%u_sub(ind)*dt
        end do
        do ind = 1,2
          part%tot_dist = part%tot_dist + part%u_sub(ind)*dt
        end do
c       ---------------------------------        
        part => part%next
      end do

      call particle_bcs_nonperiodic
      call particle_exchange
      call particle_bcs_periodic


      end subroutine SFS_velocity

      subroutine SFS_position
c     Use a more simplistic SFS treatment: Stochastic particle
c     position rather than velocity
c     Designed to be consistent with LES subgrid eddy diffusivity
c     Does not use Weil et al. 2004 formulation at all
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
      implicit none

      real :: sigm_sdxp,sigm_sdyp,sigm_sdzp,vis_sp
      real :: phim,phis,psim,psis,zeta
      real :: dadz,gasdev
      integer :: ix,iy,iz,izp1,izm1,ind,iz_part

c     ---initialize -------
      sigm_s = 0.0 

c     ------------------
c     compute sigma squre (sigm_s) based on subgrid energy field
c     -----------------       
      do iz =izs,ize
        izp1 = iz+1
        izm1 = iz-1
        do ix =1,nnx
          do iy = iys,iye
            sigm_s(ix,iy,iz) = 2.0*e(ix,iy,iz)/3.0
            vis_ss(ix,iy,iz) = vis_s(ix,iy,3,iz)
          end do
        end do
      end do

      call fill_extSFS

      !Loop over the linked list of particles:
      part => first_particle
      do while (associated(part))
        if(part%type_part == 0)then
          ! interpolate sigm_s and its derivative at particle location
          call sigm_interp(sigm_sdxp,sigm_sdyp,sigm_sdzp,vis_sp,iz_part)

          !Need vertical derivative of average vis_s:
          !Crude approximation: 0th order interpolation -- simply take
          !d(alphaC)/dz of the w-points surrounding the particle location
          if (part%xp(3) .lt. zw1) then
            !Beneath 1st zw point, use MO to approximate dadz:
            zeta = part%xp(3)/amonin
            call fzol(zeta,phim,phis,psim,psis)
            dadz = utau*vk/phis
          else
            dadz = (alphaC(iz_part+1,3)-alphaC(iz_part,3))*dzw_i(iz)
          end if

          !Now simply solve Langevin equation:
          part%xp(1) = part%xp(1)+gasdev(iseed)*sqrt(2.0*abs(vis_sp)*dt)
          part%xp(2) = part%xp(2)+gasdev(iseed)*sqrt(2.0*abs(vis_sp)*dt)
          part%xp(3) = part%xp(3)+gasdev(iseed)*sqrt(2.0*abs(vis_sp)*dt)
     +             + dadz*dt


          part => part%next
        else
          part => part%next
        end if
      end do

      call particle_bcs_nonperiodic
      call particle_exchange
      call particle_bcs_periodic


      end subroutine SFS_position
   

      subroutine particle_generation(it)
      use particles
      use pars
      use fields
      use con_data
      implicit none
      include 'mpif.h'
      real :: xv,yv,zv,ran2
      real :: rp_init
      real :: pi,xdummy,ydummy
      real :: src_rad
      integer :: idx,j,ierr
      integer :: it
      integer :: values(8),xy_found
      integer :: bin_location
      real, allocatable :: xpi(:)
       
c       If we want to reinject lost particles
c      call mpi_allreduce(reinject,tot_reinject,1,mpi_integer,mpi_sum,
c     +     mpi_comm_world,ierr)

c      if (myid == 0 .and. it .gt. part_init)then
c      write(*,*) 'total died at this timestep',tot_reinject
c      end if

c       Simply tracking how many died
      call mpi_allreduce(died,tot_died,1,mpi_integer,mpi_sum,
     +     mpi_comm_world,ierr)
      if(myid == 0 .and. inonperiodic == 1)then
        write(*,*) 'total died at this timestep', tot_died
        all_died=all_died+tot_died
      end if  
      allocate(xpi(2*(numpart_dt+tot_reinject)))

      if (it .eq. part_init)then
        t(:,:,3,:) = 0.0
      end if

      !Particle count
      numpart = 0
      part => first_particle
      do while (associated(part))
        numpart = numpart + 1
        part => part%next
      end do

      call mpi_allreduce(numpart,tnumpart,1,mpi_integer,mpi_sum,
     +                   mpi_comm_world,ierr)

c        not sure if this is currently needed       
c        if(iUprod .eq 1)then
c        else
        !Define our circular boundaries for heterogeneous source
        !production

c        end if

      if (it .gt. part_init .and. tnumpart_max .ge. tnumpart)then
c      hard code this for now
        src_rad = 267.62 

c     Homogeneous Particle creation?
        if(iUprod .eq. 1)then
          if(myid .eq. 0)then
            do j = 1,numpart_dt+tot_reinject
              xpi(j)                         = ran2(iseed)*xl
              xpi(j+numpart_dt+tot_reinject) = ran2(iseed)*yl
            end do
          end if
          call mpi_bcast(xpi,2*(numpart_dt+tot_reinject),mpi_real8,0,
     +                  mpi_comm_world,ierr)
          do idx = 1,numpart_dt+tot_reinject
            if(xpi(idx).ge.xmin .and. xpi(idx).lt.xmax .and.
     +         xpi(idx+numpart_dt+tot_reinject).ge.ymin .and.
     +         xpi(idx+numpart_dt+tot_reinject).le.ymax) then
              xv     = xpi(idx)
              yv     = xpi(idx+numpart_dt+tot_reinject)

              !Decide which heights particles need to be injected 
              zv     = zwpi
c              zv     = ran2(iseed)*(zwpi-zw1)+ zw1
              if(ipolydisp .eq. 1)then
                call create_random_partsize(rp_init)
                bin_location =minloc(class_maxr,DIM=1,
     +                               MASK=(rp_init<class_maxr))
                init_radclass(bin_location) = 
     +                        init_radclass(bin_location)+1
              else
                rp_init = radius_max
                
              end if
              radius_mass           = rp_init
              m_s                   = Sal/1000.0*3.0/4.0*
     +                                pi*radius_mass**3.0*rhow
              Npdx                  = Npdx + 1
        !Indrajith's initial random SFS velocity
              vp_init(1) = ran2(iseed)
              vp_init(2) = ran2(iseed)
              vp_init(3) = ran2(iseed)

              call create_particle((/xv,yv,zv/),vp_init,
     +        Tp_init,qf_init,Npdx,rp_init,myid,time,0.0)
              part => first_particle
              call particle_eq_radius

            end if
          end do
        else

c       Heterogeneous particle creation
          if(myid .eq. 0) then
c            call bw_change !Changing breaking wave statistics
            do j = 1,numpart_dt+tot_reinject
c       This is where the circle-based restriction will happen!
c       Theres a ton of extensions you can do with this, but just
c       get a base model working and continually add more stuff to it
c       (such as different geometries, number of sources, time-dependent
c       and spatially stochastic..)
              xy_found =0;
              do while(xy_found .eq. 0)
                xdummy=ran2(iseed)*xl
                ydummy=ran2(iseed)*yl
c       Is it within radius of the prescribed circles?
                if((xdummy-0.25*xl)**2+(ydummy-0.25*yl)**2.lt.src_rad**2
     +         .or.(xdummy-0.25*xl)**2+(ydummy-0.75*yl)**2.lt.src_rad**2
     +         .or.(xdummy-0.75*xl)**2+(ydummy-0.25*yl)**2.lt.src_rad**2
     +         .or.(xdummy-0.75*xl)**2+(ydummy-0.75*yl)**2.lt.src_rad**2
     +         )then
c       Is it within a currently occurring white-cap?
c                if(whitecap.active_wc==1)then
c                  xy_found = 1
c                end if

c             currently just finding its position
                  xy_found = 1
                end if
              end do
              xpi(j)                         = xdummy
              xpi(j+numpart_dt+tot_reinject) = ydummy
            end do
          end if

          call mpi_bcast(xpi,2*(numpart_dt+tot_reinject),mpi_real8,0,
     +                  mpi_comm_world,ierr)
          do idx = 1,numpart_dt+tot_reinject
            if(xpi(idx).ge.xmin .and. xpi(idx).lt.xmax .and.
     +         xpi(idx+numpart_dt+tot_reinject).ge.ymin .and.
     +         xpi(idx+numpart_dt+tot_reinject).le.ymax) then
              xv     = xpi(idx)
              yv     = xpi(idx+numpart_dt+tot_reinject)
              !Decide at which heights particles should be injected
              zv     = zwpi

              if(ipolydisp .eq. 1)then
                call create_random_partsize(rp_init)
                bin_location = minloc(class_maxr,DIM=1,
     +                                MASK=(rp_init<class_maxr))
                init_radclass(bin_location) = 
     +                        init_radclass(bin_location)+1
              else
                rp_init = radius_max
              end if
              radius_mass           = rp_init
              m_s                   = Sal/1000.0*3.0/4.0
     +                                *pi*radius_mass**3.0*rhow
              Npdx                  = Npdx + 1

        !Indrajith's initial random SFS velocity
              vp_init(1) = ran2(iseed)
              vp_init(2) = ran2(iseed)
              vp_init(3) = ran2(iseed)

              call create_particle((/xv,yv,zv/),vp_init,
     +          Tp_init,qf_init,Npdx,rp_init,myid,time,0.0)

              part => first_particle

              call particle_eq_radius

            end if
          end do

        end if !Whether heterogeneous or homogeneous?
      end if !If creating any particles?


      deallocate(xpi)

c        call mpi_allreduce(numpart,tnumpart,1,mpi_integer,mpi_sum,
c     +                     mpi_comm_world,ierr)
c        call mpi_allreduce(myRep_avg,Rep_avg,1,mpi_real8,mpi_sum,
c     +                     mpi_comm_world,ierr)
c        Rep_avg = Rep_avg/tnumpart


      reinject     = 0
      tot_reinject = 0
      died         = 0 
      tot_died     = 0


      !Particle count bro
      numpart = 0
      part => first_particle
      do while (associated(part))
        numpart = numpart + 1
        if(part%type_part==1)then
          numpart = numpart - 1 !Aircraft particle
        end if
        part => part%next
      end do

        !Compute total number of particles
      call mpi_allreduce(numpart,tnumpart,1,mpi_integer,mpi_sum,
     +                   mpi_comm_world,ierr)

      end subroutine particle_generation


!     3/18/2019, Now finally break it into 
!     advective and gravitational components.
      subroutine particle_rk3_Uflux
      use particles
      use pars
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'
      integer :: i,j,k,ipt,jpt,kpt,kwpt,kUpt
      integer :: ix,iy,iz
      real :: zu_1,zu_nnz
      integer :: iz_beg,iz_end
      part => first_particle
      do while (associated(part))
      if(part%type_part==0)then 
        ipt = floor(part%xp(1)/dx) + 1
        jpt = floor(part%xp(2)/dy) + 1
 
        !For particles within the domain
        if(part%xp(3) .gt. 0.0)then
          kUpt   = floor((part%xp(3)-0.0)/dz_U) + 1   
          zu_1   = 0.5*zl/nnz
          zu_nnz = zl - zu_1
 
          if(part%xp_i(3) .gt. 0.0)then
 
            if(part%xp_i(3) .lt. part%xp(3))then !Moved up
            
              iz_beg = floor((part%xp_i(3)-zu_1)/dz_U) + 2
              iz_end = floor((part%xp(3) - zu_1)/dz_U) + 1
            
              do iz=iz_beg,iz_end
                if (iz .lt. 1)then
                  write(*,*) 'up iz<1:',iz,part%xp(3),part%xp_i(3)
                end if
                if (iz .gt. nnz)then
                  write(*,*) 'up iz>nnz:',iz,part%xp(3),part%xp_i(3)
                end if 
                partflux_t(iz,jpt,ipt) = partflux_t(iz,jpt,ipt) + 1.0
              end do
            elseif (part%xp_i(3) .gt. part%xp(3))then !Moved down
              iz_beg = floor((part%xp(3)  -zu_1)/dz_U) + 2
              iz_end = floor((part%xp_i(3)-zu_1)/dz_U) + 1
              do iz=iz_beg,iz_end
                if (iz .lt. 1)then
                  write(*,*) 'dz iz<1:',iz,part%xp(3),part%xp_i(3)
                end if
                if (iz .gt. nnz)then
                  write(*,*) 'dn iz>nnz:',iz,part%xp(3),part%xp_i(3)
                end if
                partflux_t(iz,jpt,ipt) = partflux_t(iz,jpt,ipt) - 1.0
              end do
            end if
          end if
        end if
        part => part%next
      else
        part => part%next
      end if
      end do
      end subroutine particle_rk3_Uflux 

      subroutine SFS_flux
!     Same concept as particle_Uflux, now just add another thingamabob
!     to represent the subgrid flux instead of the uniform flux
      use particles
      use pars
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'
      integer :: i,j,k,ipt,jpt,kpt,kwpt,kUpt
      integer :: ix,iy,iz
      real :: zu_1,zu_nnz
      integer :: iz_beg,iz_end

      Spartflux_t = 0.0
      part => first_particle
      do while (associated(part))
      if(part%type_part==0)then  
        ipt = floor(part%xp(1)/dx) + 1
        jpt = floor(part%xp(2)/dy) + 1
 
        !For particles within the domain
        if (part%xp(3) .gt. 0.0) then
        kUpt   = floor((part%xp(3)-0.0)/dz_U) + 1 
        zu_1   = 0.5*zl/nnz
        zu_nnz = zl - zu_1
 
          if(part%xp_si(3) .gt. 0.0)then
         
            if(part%xp_si(3) .lt. part%xp(3))then !Up
            
              iz_beg = floor((part%xp_si(3)- zu_1)/dz_U) + 2
              iz_end = floor((part%xp(3)  - zu_1)/dz_U) + 1
              do iz=iz_beg,iz_end
                if(iz .lt. 1) then
                  write(*,*) 'up iz<1: ',iz,part%xp(3),part%xp_si(3)
                end if
                if(iz .gt. nnz) then
                  write(*,*) 'up iz>nnz: ',iz,part%xp(3),part%xp_si(3)
                end if
                Spartflux_t(iz,jpt,ipt) = Spartflux_t(iz,jpt,ipt) + 1.0
c              write(*,*)' (+1 vertical!)'
              end do
            elseif(part%xp_si(3) .gt. part%xp(3))then !Down
              iz_beg = floor((part%xp(3)  -zu_1)/dz_U) + 2
              iz_end = floor((part%xp_si(3)-zu_1)/dz_U) + 1
              do iz=iz_beg,iz_end
                if(iz .lt. 1)then
                  write(*,*) 'dz iz<1: ',iz,part%xp(3),part%xp_si(3)
                end if
                if(iz .gt. nnz)then
                  write(*,*) 'dn iz>nnz: ',iz,part%xp(3),part%xp_si(3)
                end if
                Spartflux_t(iz,jpt,ipt) = Spartflux_t(iz,jpt,ipt) - 1.0
c               write(*,*)' (-1 vertical!)'
              end do
            end if
          end if
        end if
        part => part%next
      else
        part => part%next
      end if
      end do

      call ztox_trans(Spartflux_t(0:nnz+1,iys:iye,mxs:mxe),
     +                Spartflux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      end subroutine SFS_flux

c     Currently this subroutine calculates the total flux
c     as well as the upward flux and downward flux
c     (the two latter quantities HAVE to add up to
c     the former mentioned total flux)

c     Now contains a subsection "probe" that reads
c     up to some dx dy

c     Now contains 16 square "probes" that reads
c     some dx dy, all of these should add to be
c     the TOTAL flux - 12/1/2018 JP
      subroutine total_flux 
      use particles
      use pars
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'
      integer :: i,j,k,ipt,jpt,kpt,kwpt,kUpt
      integer :: ix,iy,iz
      integer :: plane_d, plane_u
      real :: zu_1,zu_nnz
      integer :: iz_beg,iz_end

      tpartflux_t   = 0.0
      tupartflux_t  = 0.0
      tdpartflux_t  = 0.0
      xytpartflux_t = 0.0
      xyutpartflux_t= 0.0
      xydtpartflux_t= 0.0

      pp1flux_t     = 0.0
      pp2flux_t     = 0.0
      pp3flux_t     = 0.0
      pp4flux_t     = 0.0
      pp5flux_t     = 0.0
      pp6flux_t     = 0.0
      pp7flux_t     = 0.0
      pp8flux_t     = 0.0
      pp9flux_t     = 0.0
      pp10flux_t    = 0.0
      pp11flux_t    = 0.0
      pp12flux_t    = 0.0
      pp13flux_t    = 0.0
      pp14flux_t    = 0.0
      pp15flux_t    = 0.0
      pp16flux_t    = 0.0


!     for a discretized dx and dy of a *probe*'
      part => first_particle
      do while (associated(part))
      if(part%type_part==0)then
        ipt = floor(part%xp(1)/dx) + 1
        jpt = floor(part%xp(2)/dy) + 1
        !For particles within the domain
        if(part%xp(3) .gt. 0.0)then
            kUpt   = floor((part%xp(3)-0.0)/dz_U) + 1
            zu_1   = 0.5*zl/nnz
            zu_nnz = zl - zu_1
          if(part%xp_ti(3) .gt. 0.0)then
            if(part%xp_ti(3) .lt. part%xp(3))then !Moved up
              iz_beg = floor((part%xp_ti(3)-zu_1)/dz_U) + 2
              iz_end = floor((part%xp(3) - zu_1)/dz_U) + 1

              do iz=iz_beg,iz_end
                if(iz .lt. 1)then
                  write(*,*) 'up iz<1:',iz,part%xp(3),part%xp_ti(3)
                end if
                if(iz .gt. nnz)then
                  write(*,*) 'up iz>nnz:',iz,part%xp(3),part%xp_ti(3)
                end if
                tpartflux_t (iz,jpt,ipt) = tpartflux_t (iz,jpt,ipt) + 1.0
                tupartflux_t(iz,jpt,ipt) = tupartflux_t(iz,jpt,ipt) + 1.0

c               multiple up surface checks
                do plane_u=1,16
                  call pp_flux(plane_u,1,iz,jpt,ipt)
                end do
  


c               Probe at a set surface dx dy
                if(part%xp(1) .ge. p_l .and. part%xp(1) .le. 
     +            p_l+4.0*dx)then
                  if(part%xp(2) .ge. p_l .and. part%xp(2) .le. 
     +               p_l+4.0*dy)then
                    xytpartflux_t(iz,jpt,ipt)=xytpartflux_t(iz,jpt,ipt) 
     +                                        + 1.0
                    xyutpartflux_t(iz,jpt,ipt)=xyutpartflux_t
     +                                         (iz,jpt,ipt) + 1.0
                  end if
                end if


              end do
            elseif(part%xp_ti(3) .gt. part%xp(3))then !Moved down
              iz_beg = floor((part%xp(3)  -zu_1)/dz_U) + 2
              iz_end = floor((part%xp_ti(3)-zu_1)/dz_U) + 1
              do iz=iz_beg,iz_end
                if(iz .lt. 1)then
                  write(*,*) 'dz iz<1:',iz,part%xp(3),part%xp_ti(3)
                end if
                if(iz .gt. nnz)then
                  write(*,*) 'dn iz>nnz:',iz,part%xp(3),part%xp_ti(3)
                end if
                tpartflux_t (iz,jpt,ipt) = tpartflux_t (iz,jpt,ipt) - 1.0
                tdpartflux_t(iz,jpt,ipt) = tdpartflux_t(iz,jpt,ipt) + 1.0

c               multiple down surface checks
                do plane_d=1,16
                  call pp_flux(plane_d,2,iz,jpt,ipt)
                end do

c               Probe at a set surface dx dy
                if(part%xp(1) .ge. p_l .and. part%xp(1) .le. 
     +            p_l+4.0*dx)then
                  if(part%xp(2) .ge. p_l .and. part%xp(2) .le. 
     +               p_l+4.0*dy)then
                    xytpartflux_t(iz,jpt,ipt)=xytpartflux_t(iz,jpt,ipt)
     +                                        - 1.0
                    xydtpartflux_t(iz,jpt,ipt)=xydtpartflux_t
     +                                         (iz,jpt,ipt) + 1.0
                  end if
                end if


              end do
            end if
          end if
        end if
        part => part%next
      else
        part => part%next
      end if
      end do

!       This is where I write down probe data to get a
!       Lagrangian total flux
c        if(myid .eq. probe_loc)then
c          call probe_locator(probe_z_locations(1))
c          write(nprobe1) tpartflux_t(iploc,1,1)
c          call probe_locator(probe_z_locations(2))
c          write(nprobe2) tpartflux_t(iploc,1,1)
c          call probe_locator(probe_z_locations(3))
c          write(nprobe3) tpartflux_t(iploc,1,1)
c        end if

      call ztox_trans(tpartflux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               tpartflux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(tupartflux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               tupartflux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(tdpartflux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               tdpartflux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(xytpartflux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               xytpartflux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(xyutpartflux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               xyutpartflux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(xydtpartflux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               xydtpartflux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp1flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp1flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp2flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp2flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp3flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp3flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp4flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp4flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp5flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp5flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp6flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp6flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp7flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp7flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp8flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp8flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp9flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp9flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp10flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp10flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp11flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp11flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp12flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp12flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp13flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp13flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp14flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp14flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp15flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp15flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

      call ztox_trans(pp16flux_t(0:nnz+1,iys:iye,mxs:mxe),
     +               pp16flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +               ncpu_s,numprocs)

!      call ztox_trans(pp17flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp17flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp18flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp18flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp19flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp19flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp20flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp20flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp21flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp21flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp22flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp22flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp23flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp23flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp24flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp24flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp25flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp25flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp26flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp26flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp27flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp27flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp28flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp28flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp29flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp29flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp30flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp30flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp31flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp31flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp32flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp32flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp33flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp33flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp34flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp34flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp35flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp35flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp36flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp36flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp37flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp37flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp38flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp38flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp39flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp39flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp40flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp40flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp41flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp41flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp42flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp42flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp43flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp43flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp44flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp44flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp45flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp45flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp46flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp46flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp47flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp47flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp48flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp48flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp49flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp49flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp50flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp50flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp51flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp51flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp52flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp52flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp53flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp53flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp54flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp54flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp55flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp55flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp56flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp56flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp57flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp57flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp58flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp58flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp59flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp59flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp60flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp60flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp61flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp61flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp62flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp62flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp63flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp63flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)

!      call ztox_trans(pp64flux_t(0:nnz+1,iys:iye,mxs:mxe),
!     +               pp64flux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,
!     +               mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
!     +               ncpu_s,numprocs)




      end subroutine total_flux

c      subroutine bw_change
c      implicit none
c      include 'mpif.h'
c      !Have all processors contacted within this subroutine
c      !(MPI stuff)
c      integer :: stuff
c      real :: stuff2
c        
c      !Change spatial characteristics (Future)
c      !Change temporal characteristics (Right meow)
c        !-Active or steady?
cc     if(t_life .gt. bw%time) then
cc       bw%active = 0
cc     end
c
c      end subroutine bw_change

      subroutine write_probe_bit(nprobe)
      use pars
      use fields
      use con_data
      use con_stats
      use particles
      integer :: d_x,pcount
      real :: dummy,dummy2

c     only write if it matches the same iploc and height
      d_x = 2
c     dummy is particle count for 3 grid locations
c     dummy2 is the mean w (fluctuation) for 3 grid locations
      dummy=0
      do i=0,d_x
        dummy = dummy+partcount(ixs+i,1,iploc-1)
      enddo
c      write(*,*)'dummy:',dummy

      do i=0,d_x
        dummy2 = dummy2+w(ixs+i,1,iploc-1)
      enddo
      dummy2=dummy2/d_x !normalize by how many points added

      write(nprobe) u(ixs,iys,iploc),v(ixs,iys,iploc),
     + w(ixs,iys,iploc-1), w(ixs,iys,iploc), t(ixs,iys,3,iploc),
     + partcount(ixs,iys,iploc),
     + dummy, dummy2

c      write(*,*)'pcount(1,1,iploc-1):',partcount(ixs,1,iploc-1)
c      write(*,*)'pcount(1,1,iploc):',partcount(ixs+1,1,iploc)
c      write(*,*)'pcount(1,1,iploc+1):',partcount(ixs+2,1,iploc+1)
c      write(*,*)'w(ixs,iys,iploc-1)',w(ixs,iys,iploc-1)
c      write(*,*)'mean w of fluc',dummy2

      

      return
      end
