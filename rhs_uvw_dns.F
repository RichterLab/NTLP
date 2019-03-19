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
         w_buy = batag*(t(ix,iy,1,iz)*weit1 +
     +                  t(ix,iy,1,izp1)*weit)
c
c ------------ geostrophic wind
c
         !u_geo = -fcor*vg(iz)
         !v_geo =  fcor*(ug(iz)-ugal)
         !Instead of geostrophic wind (which is a pressure gradient)
         !make u_geo and v_geo equal to my pressure gradient:

         u_geo = 0.0
         v_geo = 0.0

c
c ------------ totals
c
         r1(ix,iy,iz) = u_adv + u_cor + u_geo
         r2(ix,iy,iz) = v_adv + v_cor + v_geo
         r3(ix,iy,iz) = w_adv + w_cor + w_buy

         !Add particle momentum coupling
         if (icouple == 1) then
         r1(ix,iy,iz) = r1(ix,iy,iz) - partsrc(ix,iy,iz,1)
         r2(ix,iy,iz) = r2(ix,iy,iz) - partsrc(ix,iy,iz,2)
         !Note: partsrc(3,ix,iy,iz) is located at u,v-locations
         !Interpolate to w-location:
         r3(ix,iy,iz) = r3(ix,iy,iz) - (weit*partsrc(ix,iy,izp1,3)+
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
         fc_flx = 0.0
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




c -------- JP Sept 2017 David's Eulerian Concentration: Vertical
c advection due to gravitational settling (iscl == 3 only assuming 1,2
c are temp, humidity)
c      if(myid==0)then
c         write(*,*) "ws eul conc ",ws
c      endif

      if (iscl==3) then
      if (iz .eq. 1) then
      
      do iy=iys,iye
      do ix=1,nnx

         theta_u = ws*(weit1*t(ix,iy,iscl,iz) +
     +                weit*t(ix,iy,iscl,izp1))
         theta_l = 0.0
         r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz)
     +     +(theta_u - theta_l)*dzw_i(iz)

      enddo
      enddo

      else

      do iy=iys,iye
      do ix=1,nnx
         theta_u = weit1*t(ix,iy,iscl,iz) +
     +                weit*t(ix,iy,iscl,izp1)
         theta_l = weit3*t(ix,iy,iscl,iz) +
     +                weit4*t(ix,iy,iscl,izm1)
         r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz)
     +     +ws*(theta_u - theta_l)*dzw_i(iz)
      enddo
      enddo
      
      end if
      end if

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
