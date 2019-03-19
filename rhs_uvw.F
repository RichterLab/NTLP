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
         w_buy = batag*(t(ix,iy,1,iz)*weit1 +
     +                  t(ix,iy,1,izp1)*weit)
c
c ------------ geostrophic wind
c
         u_geo = -fcor*vg(iz)
         v_geo =  fcor*(ug(iz)-ugal)
c
c ------------ totals
c
         r1(ix,iy,iz) = u_adv + u_cor + u_geo
         r2(ix,iy,iz) = v_adv + v_cor + v_geo
         r3(ix,iy,iz) = w_adv + w_cor + w_buy
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
       
c   ------Debug only REMOVE BELOW SECTION LATER---------
c         if(ix.eq.32.and.iy.eq.32.and.iz.eq.32)then
c         write(*,*)'rhsuvw-tau 11,12,13,23',it,fnt1(ix,iy),fnt2(ix,iy),
c     +          tau13_u(ix,iy),tau23_u(ix,iy)
c        end if
c   -----REMOVE ABOVE AFTER DEBUG-----------------
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

c   ------Debug only REMOVE BELOW SECTION LATER---------
c         if(ix.eq.32.and.iy.eq.32.and.iz.eq.32)then
c         write(*,*)'rhsuvw-tau33',it,fnt2(ix,iy),fnt4(ix,iy)
c        end if
c   -----REMOVE ABOVE AFTER DEBUG-----------------

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

c   ------Debug only REMOVE BELOW SECTION LATER---------
c         if(ix.eq.32.and.iy.eq.32.and.iz.eq.32)then
c         write(*,*)'rhsuvw-tau 21',it,fntd(ix,iy,iz)
c        end if
c   -----REMOVE ABOVE AFTER DEBUG-----------------

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

c   ------Debug only REMOVE BELOW SECTION LATER---------
c         if(ix.eq.32.and.iy.eq.32.and.iz.eq.32)then
c         write(*,*)'rhsuvw-tau 22',it,fntd(ix,iy,iz)
c        end if
c   -----REMOVE ABOVE AFTER DEBUG-----------------

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
