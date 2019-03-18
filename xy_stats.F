      subroutine xy_stats
c
c ------------ get statistics 
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
c
c ------- indices for indexing array stat(.,.)
c         js = number of non-scalar stats
c         ns = number of scalar stats
c
      parameter(js = 61, ns = 5, nstat = js + ns*nscl+1)
      real stat(1:nnz,nstat)
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
c          stat(.,11) = particle number in each cell 
c          stat(.,12) = vpsum(1) 
c          stat(.,13) = vpsum(2) 
c          stat(.,14) = vpsum(3) 
c          stat(.,15) = vpsqrsum(1) 
c          stat(.,16) = vpsqrsum(2) 
c          stat(.,17) = vpsqrsum(3) 
c          stat(.,18) = partsrc(1) 
c          stat(.,19) = partsrc(2) 
c          stat(.,20) = partsrc(3) 
c          stat(.,21) = upwp - up'*wp'
c          stat(.,22) = Tpsum
c          stat(.,23) = Tpsqrsum
c          stat(.,24) = Tfsum
c          stat(.,25) = qfsum
c          stat(.,26) = wpTpsum
c          stat(.,27) = Tpsrc
c          stat(.,28) = radsum
c          stat(.,29) = rad2sum
c          stat(.,30) = mass_p_sum
c          stat(.,31) = mass_w_sum
c          stat(.,32) = Hpsrc
c          stat(.,33) = TEpsrc
c          stat(.,34) = qstarsum
c          stat(.,35) = total uniform particle count
c          stat(.,36) = (resolved for LES) advection and gravitational
c                       flux
c          stat(.,37) = subgrid particle flux
c          stat(.,38) = total particle flux
c          stat(.,39) = scalar viscosity..?
c          stat(.,40) = total up particle flux
c          stat(.,41) = total down particle flux
c          stat(.,42) = probe total particle flux
c
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
c
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
         stat(iz,11) = stat(iz,11) + partcount(ix,iy,iz)
         stat(iz,12) = stat(iz,12) + vpsum(ix,iy,iz,1)
         stat(iz,13) = stat(iz,13) + vpsum(ix,iy,iz,2)
         stat(iz,14) = stat(iz,14) + vpsum(ix,iy,iz,3)
         stat(iz,15) = stat(iz,15) + vpsqrsum(ix,iy,iz,1)
         stat(iz,16) = stat(iz,16) + vpsqrsum(ix,iy,iz,2)
         stat(iz,17) = stat(iz,17) + vpsqrsum(ix,iy,iz,3)
         stat(iz,18) = stat(iz,18) + partsrc(ix,iy,iz,1)
         stat(iz,19) = stat(iz,19) + partsrc(ix,iy,iz,2)
         stat(iz,20) = stat(iz,20) + partsrc(ix,iy,iz,3)
         stat(iz,21) = stat(iz,21) + upwp(ix,iy,iz)
         stat(iz,22) = stat(iz,22) + Tpsum(ix,iy,iz)
         stat(iz,23) = stat(iz,23) + Tpsqrsum(ix,iy,iz)
         stat(iz,24) = stat(iz,24) + Tfsum(ix,iy,iz)
         stat(iz,25) = stat(iz,25) + qfsum(ix,iy,iz)
         stat(iz,26) = stat(iz,26) + wpTpsum(ix,iy,iz)
         stat(iz,27) = stat(iz,27) + partTsrc(ix,iy,iz)
         stat(iz,28) = stat(iz,28) + radsum(ix,iy,iz) 
         stat(iz,29) = stat(iz,29) + rad2sum(ix,iy,iz) 
         stat(iz,30) = stat(iz,30) + mpsum(ix,iy,iz)
         stat(iz,31) = stat(iz,31) + mwsum(ix,iy,iz)
         stat(iz,32) = stat(iz,32) + partHsrc(ix,iy,iz)
         stat(iz,33) = stat(iz,33) + partTEsrc(ix,iy,iz)
         stat(iz,34) = stat(iz,34) + qstarsum(ix,iy,iz)
c         stat(iz,35) =stat(iz,35) +  partcount(ix,iy,iz)
         stat(iz,35) = stat(iz,35) + partcountU(ix,iy,iz)
         stat(iz,36) = stat(iz,36) + partflux(ix,iy,iz)  !This is RK3
         stat(iz,37) = stat(iz,37) + Spartflux(ix,iy,iz) !SFS
         stat(iz,38) = stat(iz,38) + tpartflux(ix,iy,iz) !=Total
         stat(iz,39) = stat(iz,39) + vis_s(ix,iy,3,iz)
         stat(iz,40) = stat(iz,40) + tupartflux(ix,iy,iz)
         stat(iz,41) = stat(iz,41) + tdpartflux(ix,iy,iz)
         stat(iz,42) = stat(iz,42) + xytpartflux(ix,iy,iz)
         stat(iz,43) = stat(iz,43) + xyutpartflux(ix,iy,iz)
         stat(iz,44) = stat(iz,44) + xydtpartflux(ix,iy,iz)
         stat(iz,45) = stat(iz,45) + pcountrange(ix,iy,iz)
         stat(iz,46) = stat(iz,46) + pp1flux(ix,iy,iz)
         stat(iz,47) = stat(iz,47) + pp2flux(ix,iy,iz)
         stat(iz,48) = stat(iz,48) + pp3flux(ix,iy,iz)
         stat(iz,49) = stat(iz,49) + pp4flux(ix,iy,iz)
         stat(iz,50) = stat(iz,50) + pp5flux(ix,iy,iz)
         stat(iz,51) = stat(iz,51) + pp6flux(ix,iy,iz)
         stat(iz,52) = stat(iz,52) + pp7flux(ix,iy,iz)
         stat(iz,53) = stat(iz,53) + pp8flux(ix,iy,iz)
         stat(iz,54) = stat(iz,54) + pp9flux(ix,iy,iz)
         stat(iz,55) = stat(iz,55) + pp10flux(ix,iy,iz)
         stat(iz,56) = stat(iz,56) + pp11flux(ix,iy,iz)
         stat(iz,57) = stat(iz,57) + pp12flux(ix,iy,iz)
         stat(iz,58) = stat(iz,58) + pp13flux(ix,iy,iz)
         stat(iz,59) = stat(iz,59) + pp14flux(ix,iy,iz)
         stat(iz,60) = stat(iz,60) + pp15flux(ix,iy,iz)
         stat(iz,61) = stat(iz,61) + pp16flux(ix,iy,iz)





      enddo                                              !
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
      ups(iz)        = stat(iz,1)*fnxy
      vps(iz)        = stat(iz,2)*fnxy
      wps(iz)        = stat(iz,3)*fnxy
      wcube(iz)      = stat(iz,4)*fnxy
      wfour(iz)      = stat(iz,5)*fnxy
      englez(iz)     = stat(iz,6)*fnxy
      engsbz(iz)     = stat(iz,7)*fnxy
      eavg(iz)       = stat(iz,8)*fnxy
      uwle(iz)       = stat(iz,9)*fnxy
      vwle(iz)       = stat(iz,10)*fnxy
      zconc(iz)      = stat(iz,11)/xl/yl/dzw(iz)
      zconc_U(iz)    = stat(iz,35)/xl/yl/dz_U
      zconc_range(iz)= stat(iz,45)/xl/yl/dzw(iz)
      vp1mean(iz)    = stat(iz,12)/stat(iz,11)
      vp2mean(iz)    = stat(iz,13)/stat(iz,11)
      vp3mean(iz)    = stat(iz,14)/stat(iz,11)
      vp1msqr(iz)    = sqrt(stat(iz,15)/stat(iz,11)-vp1mean(iz)**2)
      vp2msqr(iz)    = sqrt(stat(iz,16)/stat(iz,11)-vp2mean(iz)**2)
      vp3msqr(iz)    = sqrt(stat(iz,17)/stat(iz,11)-vp3mean(iz)**2)
      m1src(iz)      = stat(iz,18)*fnxy
      m2src(iz)      = stat(iz,19)*fnxy
      m3src(iz)      = stat(iz,20)*fnxy
      uw_tot(iz)     = uwle(iz) + uwsb(iz)
      vw_tot(iz)     = vwle(iz) + vwsb(iz)
      upwpm(iz)      = stat(iz,21)/stat(iz,11)-(vp1mean(iz)*vp3mean(iz))
      Tpmean(iz)     = stat(iz,22)/stat(iz,11)
      Tpmsqr(iz)     = sqrt(stat(iz,23)/stat(iz,11)-Tpmean(iz)**2)
      Tfmean(iz)     = stat(iz,24)/stat(iz,11)
      qfmean(iz)     = stat(iz,25)/stat(iz,11)
      wpTpm(iz)      = stat(iz,26)/stat(iz,11)-(Tpmean(iz)*vp3mean(iz))
      Tpsrc(iz)      = stat(iz,27)*fnxy
      radmean(iz)    = stat(iz,28)/stat(iz,11) 
      rad2mean(iz)   = stat(iz,29)/stat(iz,11)
 
      mpmean(iz)     = stat(iz,30)/stat(iz,11) !make sure stat matchs with phimp
      mwmean(iz)     = stat(iz,31)/stat(iz,11) !make sure stat matchs with phimw
      phimp(iz)      = stat(iz,30)/xl/yl/dzw(iz)/rhoa
      phimw(iz)      = stat(iz,31)/xl/yl/dzw(iz)/rhoa
 
      Hpsrc(iz)      = stat(iz,32)*fnxy
      TEpsrc(iz)     = stat(iz,33)*fnxy
      qstarm(iz)     = stat(iz,34)*fnxy

      rk3pflux(iz)   = stat(iz,36)
      Spflux(iz)     = stat(iz,37)
      Tpflux(iz)     = stat(iz,38)

      vismean(iz)    = stat(iz,39)
      TUpflux(iz)    = stat(iz,40)
      TDpflux(iz)    = stat(iz,41)
      probepflux(iz) = stat(iz,42)



      pp1(iz)        = stat(iz,46)
      pp2(iz)        = stat(iz,47)
      pp3(iz)        = stat(iz,48)
      pp4(iz)        = stat(iz,49)
      pp5(iz)        = stat(iz,50)
      pp6(iz)        = stat(iz,51)
      pp7(iz)        = stat(iz,52)
      pp8(iz)        = stat(iz,53)
      pp9(iz)        = stat(iz,54)
      pp10(iz)       = stat(iz,55)
      pp11(iz)       = stat(iz,56)
      pp12(iz)       = stat(iz,57)
      pp13(iz)       = stat(iz,58)
      pp14(iz)       = stat(iz,59)
      pp15(iz)       = stat(iz,60)
      pp16(iz)       = stat(iz,61)



      
     
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
c
      return
      end
