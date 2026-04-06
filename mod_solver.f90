! ----------------------------------------------------------------------
! mod_solver — combined SGS, RHS, and physics/boundary-condition routines
!
! Combines what would have been mod_sgs, mod_rhs, and mod_physics into a
! single module to avoid circular dependencies (the three groups call
! each other freely within the time-stepping loop).
!
! Public interface covers all subroutines called from the main program
! or from particles.f90. Helper functions and internal routines are
! private.
! ----------------------------------------------------------------------
module mod_solver
  use mpi
  use pars
  use fields
  use fftwk
  use con_data
  use con_stats
  use particles
  use mod_mpi
  use mod_fft
  implicit none
  private

  ! Time-step control
  public :: get_max, get_dt
  ! SGS turbulence
  public :: dns_vis, tke_vis, iso, surfvis, stokesv, busngr, fzol
  ! RHS / time-stepping core
  public :: comp1, rhs_uvw, rhs_uvw_DNS, rhs_scl, rhs_scl_dns
  public :: comp_p, comp2, pressure, solve_trid, tridv
  public :: get_derv, get_means, forcing, extra_flux_terms
  ! Boundary conditions and surface-layer physics
  public :: lower, lower_dns, upper_dns, lower_free, f_suft2, upper
  public :: suft, sufto, suft2, lterp
  ! Diagnostics
  public :: pbltop, get_zi
  ! Moisture / radiation
  public :: calc_radiation, humidity_control, change_RH_bcs_to_q, fill_base

contains

  ! ──────────────────────────────────────────────────────────────────
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
      end subroutine get_max

  ! ──────────────────────────────────────────────────────────────────
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
      end subroutine get_dt

  ! ──────────────────────────────────────────────────────────────────
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
      end subroutine lterp

  ! ──────────────────────────────────────────────────────────────────
      subroutine iso
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use mod_mpi
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
      end subroutine iso

  ! ──────────────────────────────────────────────────────────────────
      subroutine surfvis
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use mod_mpi
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
      end subroutine surfvis

  ! ──────────────────────────────────────────────────────────────────
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

      end subroutine dns_vis

  ! ──────────────────────────────────────────────────────────────────
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
      use mod_fft
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
      call xderivp(ex(1,iys),trigx(:,1),xk,nnx,iys,iye)
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
      call xderivp(fnt1(1,iys),trigx(:,1),xk,nnx,iys,iye)
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
      call yd_mpi(ey(1,iys,izs),trigx(:,2),yk,
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
      call yd_mpi(fnt2(1,iys,izs),trigx(:,2),yk,
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
      end subroutine tke_vis

  ! ──────────────────────────────────────────────────────────────────
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
      end subroutine stokesv

  ! ──────────────────────────────────────────────────────────────────
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
      end subroutine busngr

  ! ──────────────────────────────────────────────────────────────────
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
      end subroutine fzol

  ! ──────────────────────────────────────────────────────────────────
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
      use mod_mpi
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
      end subroutine comp1

  ! ──────────────────────────────────────────────────────────────────
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
      use mod_mpi
      use mod_fft
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
      call xderivp(fnt1(1,iys),trigx(:,1),xk,nnx,iys,iye)
      call xderivp(fnt2(1,iys),trigx(:,1),xk,nnx,iys,iye)
      call xderivp(fnt3(1,iys),trigx(:,1),xk,nnx,iys,iye)
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
      call yd_mpi(fntd(1,iys,izs),trigx(:,2),yk,
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
      call yd_mpi(fntd(1,iys,izs),trigx(:,2),yk,
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
      call yd_mpi(fntd(1,iys,izs),trigx(:,2),yk,
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
      end subroutine rhs_uvw

  ! ──────────────────────────────────────────────────────────────────
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
      use mod_mpi
      use mod_fft
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
      call xderivp(fnt1(1,iys),trigx(:,1),xk,nnx,iys,iye)
      call xderivp(fnt2(1,iys),trigx(:,1),xk,nnx,iys,iye)
      call xderivp(fnt3(1,iys),trigx(:,1),xk,nnx,iys,iye)
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
      call yd_mpi(fntd(1,iys,izs),trigx(:,2),yk,
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
      call yd_mpi(fntd(1,iys,izs),trigx(:,2),yk,
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
      call yd_mpi(fntd(1,iys,izs),trigx(:,2),yk,
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

  ! ──────────────────────────────────────────────────────────────────
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
      use mod_fft
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
      call xderivp(tx(1,iys),trigx(:,1),xk,nnx,iys,iye)
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

      call xderivp(fnt1(1,iys,iz),trigx(:,1),xk,nnx,iys,iye)
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
      call yd_mpi(ty(1,iys,izs),trigx(:,2),yk,
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
      call yd_mpi(fnt1(1,iys,izs),trigx(:,2),yk,
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
      end subroutine rhs_scl

  ! ──────────────────────────────────────────────────────────────────
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
      use mod_mpi
      use mod_fft
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
      call xderivp(tx(1,iys),trigx(:,1),xk,nnx,iys,iye)
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
      call xderivp(fnt1(1,iys,iz),trigx(:,1),xk,nnx,iys,iye)
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
      call yd_mpi(ty(1,iys,izs),trigx(:,2),yk,
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
      call yd_mpi(fnt1(1,iys,izs),trigx(:,2),yk,
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

  ! ──────────────────────────────────────────────────────────────────
      subroutine comp_p
c
c --------- setup pressure solver
c
      use pars
      use fftwk
      use fields
      use con_data
      use con_stats
      use mod_fft
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
         call xderivp(fnt1(1,iys,iz),trigx(:,1),xk,nnx,iys,iye)
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
      call yd_mpi(fnt1(1,iys,izs),trigx(:,2),yk,
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
      end subroutine comp_p

  ! ──────────────────────────────────────────────────────────────────
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
      use mod_mpi
      use mod_fft
      include 'mpif.h'
      real fnt1(nnx,iys:iye,izs:ize), fnt2(nnx,iys:iye)
      real r3_sum(1:nnz)
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
      call yd_mpi(fnt1(1,iys,izs),trigx(:,2),yk,
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
         call xderivp(fnt2(1,iys),trigx(:,1),xk,nnx,iys,iye)
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
      end subroutine comp2

  ! ──────────────────────────────────────────────────────────────────
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
      use mod_mpi
      use mod_fft
      include 'mpif.h'
      real pfft(nny,jxs:jxe,izs-1:ize+1)
      real pt(0:nnz+1,jxs:jxe,jys:jye)
      real ptopfft(nny,jxs:jxe,1:2)
      real psum(1:nnz)
      integer istatus(mpi_status_size)
c
c ------------ Fourier analyze the right hand side
c              at all iz = izs,ize. results are in pfft
c
c


      call fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(:,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)

c
c ------------ Fourier analyze the radiation bc arrays
c
      if(ibcu .eq. 1) then
        call fft2d_mpi(ptop(1,iys,1),ptopfft(1,jxs,1),
     +           trigx(:,1),trigc,
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
      call fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(:,1),trigc,
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
      end subroutine pressure

  ! ──────────────────────────────────────────────────────────────────
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
      end subroutine solve_trid

  ! ──────────────────────────────────────────────────────────────────
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
      end subroutine tridv

  ! ──────────────────────────────────────────────────────────────────
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
      use mod_mpi
      use mod_fft
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
         call xderivp(ux(1,iys,iz),trigx(:,1),xk,
     +                 nnx,iys,iye)
         call xderivp(vx(1,iys,iz),trigx(:,1),xk,
     +                 nnx,iys,iye)
         call xderivp(wx(1,iys,iz),trigx(:,1),xk,
     +                 nnx,iys,iye)
      enddo
c
c ---------- get y derivatives for (u,v,w)
c
c     call yd_mpi(uy(1,iys,iz_ss),trigx(:,2),yk,
c    +           nnx,nny,ixs,ixe,ix_s,ix_e,
c    +           iys,iye,iy_s,iy_e,iz_ss,iz_ee,myid,ncpu_s,numprocs)
c     call yd_mpi(vy(1,iys,iz_ss),trigx(:,2),yk,
c    +           nnx,nny,ixs,ixe,ix_s,ix_e,
c    +           iys,iye,iy_s,iy_e,iz_ss,iz_ee,myid,ncpu_s,numprocs)
c     call yd_mpi(wy(1,iys,iz_ss),trigx(:,2),yk,
c    +           nnx,nny,ixs,ixe,ix_s,ix_e,
c    +           iys,iye,iy_s,iy_e,iz_ss,iz_ee,myid,ncpu_s,numprocs)
c
      call yd_mpi(uy(1,iys,izs-1),trigx(:,2),yk,
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
      call yd_mpi(vy(1,iys,izs-1),trigx(:,2),yk,
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
      call yd_mpi(wy(1,iys,izs-1),trigx(:,2),yk,
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
c
      return
      end subroutine get_derv

  ! ──────────────────────────────────────────────────────────────────
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
      use mod_mpi
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
      call mpi_sum_z(engz,i_root,myid,nnzp1,1)
      call mpi_sum_z(engsbz,i_root,myid,nnzp1,1)
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
        call mpi_sum_z(divz,i_root,myid,nnz,1)
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
      end subroutine get_means

  ! ──────────────────────────────────────────────────────────────────
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
      end function rlim

  ! ──────────────────────────────────────────────────────────────────
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
      end function ran1

  ! ──────────────────────────────────────────────────────────────────
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
      end function ranf

  ! ──────────────────────────────────────────────────────────────────
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
      end subroutine forcing

  ! ──────────────────────────────────────────────────────────────────
      subroutine extra_flux_terms
      use pars
      use particles
      use fields
      use con_data
      use con_stats
      use fftwk
      use mod_mpi
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
      end subroutine extra_flux_terms

  ! ──────────────────────────────────────────────────────────────────
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
      use mod_mpi
      real sfc_flx(2+nscl)
      real upars(2)
c

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
      end subroutine lower

  ! ──────────────────────────────────────────────────────────────────
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

  ! ──────────────────────────────────────────────────────────────────
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

  ! ──────────────────────────────────────────────────────────────────
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
      use mod_mpi
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
      end subroutine lower_free

  ! ──────────────────────────────────────────────────────────────────
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
      end subroutine f_suft2

  ! ──────────────────────────────────────────────────────────────────
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
      use mod_mpi
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
      end subroutine upper

  ! ──────────────────────────────────────────────────────────────────
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
      use mod_mpi
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
      end subroutine suft

  ! ──────────────────────────────────────────────────────────────────
      subroutine sufto
c
      use pars
      use fields
      use con_data
      use con_stats
      use mod_mpi
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
      end subroutine sufto

  ! ──────────────────────────────────────────────────────────────────
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
      end subroutine suft2

  ! ──────────────────────────────────────────────────────────────────
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
      end subroutine pbltop

  ! ──────────────────────────────────────────────────────────────────
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
      end subroutine get_zi

  ! ──────────────────────────────────────────────────────────────────
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
      real :: rhoa,func_rho_base,Ttmp,exner,func_p_base
      integer :: iz,lc

      !Do trapezoidal integration from the TOP down:
      !LWC only calculated on INTERIOR zw points

      lc = 1
      do iz=nnz-1,1,-1

         if (lc .eq. 1) then

         if (iexner.eq.1) then
            Ttmp = txym(iz,1)*
     +      exner(surf_p,func_p_base(surf_p,tsfcc(1),zz(iz)))
            !rhoa = func_rho_base(surf_p,tsfcc(1),zz(iz))
	    rhoa = func_p_base(surf_p,tsfcc(1),zz(iz))/Rd/txym(iz,1)
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
            Ttmp = txym(iz,1)*
     +      exner(surf_p,func_p_base(surf_p,tsfcc(1),zz(iz)))
            !rhoa = func_rho_base(surf_p,tsfcc(1),zz(iz))
	    rhoa = func_p_base(surf_p,tsfcc(1),zz(iz))/Rd/txym(iz,1)
         else
            rhoa = surf_rho
         end if
         radsrc(iz) = -1/rhoa/Cpa*(radflux(iz)-radflux(iz-1))/dzw(iz)

      end do


      end subroutine calc_radiation


  ! ──────────────────────────────────────────────────────────────────
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


  ! ──────────────────────────────────────────────────────────────────
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


  ! ──────────────────────────────────────────────────────────────────
      subroutine fill_base(T_surf)
      !Fill the base state profiles p_base,rho_base,T_base,theta_base
      !given the surface conditions T_surf and p_surf.
      !theta_base is assumed equal to T_surf (not some other constant
      !reference)
      use fields
      use pars
      use con_stats
      use con_data
      implicit none

      integer :: iz
      real, intent(in) :: T_surf

          
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


  ! ──────────────────────────────────────────────────────────────────
      function mod_magnus(T)
      implicit none

      !Take in T in Kelvin and return saturation vapor pressure using function of Alduchov and Eskridge, 1996
      real,intent(in) :: T

      mod_magnus = 610.94 *exp((17.6257*(T-273.15))/(243.04+(T-273.15)))


      end function mod_magnus

  ! ──────────────────────────────────────────────────────────────────
      function exner(p0,p)
      use pars
      implicit none
      real, intent(in) :: p0,p

      !Take in the reference pressure p0 and p(z), and compute exner = T/theta
      exner = (p0/p)**(-Rd/Cpa)

      end function exner

      !!!!!The dry adiabatic functions:

  ! ──────────────────────────────────────────────────────────────────
      function func_p_base(p_surf,T_surf,z)
      use pars, only: grav,Cpa,Rd
      implicit none
      real, intent(in) :: p_surf,T_surf,z

      func_p_base = p_surf*(1 - z*grav/Cpa/T_surf)**(Cpa/Rd)

      end function func_p_base

  ! ──────────────────────────────────────────────────────────────────
      function func_T_base(T_surf,z)
      use pars, only: grav,Cpa
      implicit none
      real, intent(in) :: T_surf,z
      
      func_T_base = T_surf - grav/Cpa*z

      end function func_T_base

  ! ──────────────────────────────────────────────────────────────────
      function func_rho_base(p_surf,T_surf,z)
      use pars, only: grav,Cpa,Rd
      implicit none
      real, intent(in) :: p_surf,T_surf,z
      
      func_rho_base = p_surf*T_surf**(-Cpa/Rd)/Rd*
     +         (T_surf - grav*z/Cpa)**(Cpa/Rd-1.0)

      end function func_rho_base
      !!!!!

 

end module mod_solver
