      subroutine lower(it)
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
      real sfc_flx(2+nscl)
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
         call sufto(it)
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
         call suft(it)
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
c            dnom3=t10xy(iscl)*windm
c            if(dnom3 .ne. 0.) then
c               dnom_i = 1.0/dnom3
c               do iy=iys,iye
c               do ix=1,nnx
c                  taut3m(ix,iy,iscl)=aut3m(iscl)*
c     +                 (windm*(t(ix,iy,iscl,iz)-t1xy(iscl))+
c     +                  wind(ix,iy)*(t1xy(iscl)-tsfcc(iscl)))*dnom_i
c               enddo
c               enddo
c            else
               do iy=iys,iye
               do ix=1,nnx
                  taut3m(ix,iy,iscl) = aut3m(iscl)
               enddo
               enddo
c            endif
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
      enddo
c     write(nprt,2345) uwsfc, vwsfc, wtsfc(nscl), tsfcc(nscl)
 2345 format(' in lower 2345 uwsfc = ',e15.6,' vwsfc = ',e15.6,
     +       ' wtsfc = ',e15.6,' tsfcc = ',e15.6)
c
      do iy=iys,iye
      do ix=1,nnx
         dudz     = 2.*(u(ix,iy,iz) + ugal)*dz_i
         dvdz     = 2.*v(ix,iy,iz)*dz_i
         ubc(ix,iy,2) = u(ix,iy,iz) - dudz*dzu(iz)
         vbc(ix,iy,2) = v(ix,iy,iz) - dvdz*dzu(iz)
      enddo
      enddo
      do iscl=1,nsl
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
