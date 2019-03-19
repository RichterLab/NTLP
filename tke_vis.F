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
c          vis_s(i,j,1,iz) = (1.+2.*alk(i,j,iz)/dslk)*vis_m(i,j,iz)
c          vis_s(i,j,2:nscl,iz) = vis_m(i,j,iz) 
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
