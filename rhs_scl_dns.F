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
          r4(ix,iy,1,iz) = r4(ix,iy,1,iz) - partTsrc(ix,iy,iz)
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
          r4(ix,iy,1,iz) = r4(ix,iy,1,iz) - partTEsrc(ix,iy,iz)
          r4(ix,iy,2,iz) = r4(ix,iy,2,iz) - partHsrc(ix,iy,iz)
         end do
         end do
      end do
      end if
      end if 



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
!      vis_s(1:nnx,iys:iye,2,izs-1:ize+1) = nuf/Sc    !Dv =nu/Sc
      r5 = 0.0
      e = 0.0

      end
