module mod_io
  use mpi
  use pars
  use fields
  use fftwk
  use con_data
  use con_stats
  use mod_mpi
  use mod_fft
  use mod_solver
  use mod_thermo
  use mod_init
  use particles
  private
  public :: nblnk, blnk, print, xy_stats, tke_budget, Tvar_budget, write_his, write_prof, close_histograms, close_his, save_viz, recv_yz_var, save_v, save_c, save_p, get_units, get_output_filenames, open_histograms, open_his, viz_output_filename, open_viz
!ontains
      subroutine nblnk(word)
      implicit real(a-h,o-z), integer(i-n)
      parameter (nmax=304)
      character wordt*304, word*(*)
      nchar = len(word)
      if(nchar .gt. nmax) then
         write(6,6000) nchar,nmax
 6000    format(' TROUBLE, IN SR. NBLNK : NCHAR = ',i6, &
                ' EXCEEDS NMAX = ',i6)
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
!
      return
      end subroutine nblnk

      subroutine blnk(word)
      implicit real(a-h,o-z), integer(i-n)
      character word*(*)
      nchar = len(word)
      do j=1,nchar
         word(j:j) = ' '
      enddo
!
      return
      end subroutine blnk

      subroutine print(lu,iz_strt,iz_end)
!
!
      implicit real(a-h,o-z), integer(i-n)
      write(lu,4000)
 4000 format(30X,' --- SOLUTION ---')
      write(lu,4100) it,time,dt,zi,tsfcc(1),uwsfc,vwsfc,wtsfc(1), &
                    zol,hol,ucfl, vcfl, wcfl, &
                    t_ref
 4100 format(' IT=',I7,5x,'TIME (s) = ',e15.8,',  DT(s) = ',e15.6,/, &
             10x,'ZTOP = ',e15.6, &
             ',  TSFC = ',e15.6, &
             ',  UW = ',e15.6,',  VW = ',e15.6,/,10x, &
             'WT = ',e15.6,',  ZL =',e15.6, &
             ',  HL = ',e15.6,/,10x,'U_cfl = ',e15.6, &
             ',  V_cfl = ',e15.6,',  W_cfl = ',e15.6,/,10x, &
             'Theta Ref = ',e15.6)
      write(lu,4200)
 4200 format(//,20x,'--------- HORIZONTAL MEAN VALUES ---------- ', &
             //,2x,'IZ',4x,'T_MEAN',7x,'T2_MEAN',6x, &
       'DIVG',8X,'LE_KE',6X,'SGS_KE',7X,'LE_WT',6X, &
       'SGS_WT',7X,'SHRZ',8X,'BUOY')
      do 19999 iz=iz_end,iz_strt,-1
         !write(lu,4300)iz,txym(iz,1)-t_ref,divz(iz),
         write(lu,4300)iz,txym(iz,1)-t_ref, txym(iz,2),divz(iz), &
                    englez(iz),eavg(iz),wtle(iz,1), &
                    wtsb(iz,1),shrz(iz),buyz(iz)
 4300    format(1X,I3,e12.4,8e12.4)
19999 continue
      write(lu,4400)tsfcc(1),wtsfc(1)
 4400 format('  SURFACE VALUE: TXYM=',F8.2,'               WTSB=',E9.2)
      if(iocean .eq. 1) then
         write(lu,4500) stokess,udrift,vdrift
 4500    format(/,' STOKESS = ',e12.4,' UDRIFT = ',e12.4, &
                ' VDRIFT = ',e12.4)
      endif
      write(lu,4600) (iz,uxym(iz)+ugal,vxym(iz),uwle(iz), &
             uwsb(iz),vwle(iz),vwsb(iz),iz=iz_strt,iz_end)
 4600 format(//,' IZ',5x,' UXYM + UGAL',8x,' VXYM',10x,' UWLE',10x, &
                ' UWSB',10x,' VWLE',10x,' VWSB' &
             ,/,(1x,i4,6(3x,e15.6)))
      if(ivis .eq. 1) then
         write(lu,4800) xksurf, nmatch, viscon, vise
 4800    format(//,' XKSURF = ',e15.6,' NMATCH = ',i4,/, &
                   ' VISCON = ',e15.6,' VISE = ',e15.6)
!         write(lu,4700) (iz,dfac(iz),iz=iz_strt,iz_end)
! 4700    format(//,'   IZ',5x,'  DFAC',/,(1x,i4,3x,e15.6))
      endif
!
! --------------- output additional scalars
!
      if(nscl .eq. 2) then
      write(lu,5005)tsfcc(2),wtsfc(2)
 5005 format(/,'  SURFACE VALUE: TXYM(2) =',e15.6,' WTSFC(2) = ',e15.6)
      write(lu,5100) (iz,txym(iz,2),wtle(iz,2), &
                    wtsb(iz,2),iz=iz_strt,iz_end)
 5100 format(//,' IZ',5x,' SCALAR-1 MEAN',8x,' WS1LE',10x, &
                ' WS1SB',10x &
             ,/,(1x,i4,3(3x,e15.6)))
      else if (nscl .eq. 3) then
!      write(lu,5205)tsfcc(2),wtsfc(2),tsfcc(3),wtsfc(3)
 5205 format(/,'  SURFACE VALUE: TXYM(2) =',e15.6,' WTSFC(2) = ',e15.6, &
             /,'  SURFACE VALUE: TXYM(3) =',e15.6,' WTSFC(3) = ',e15.6)
      write(lu,5200) (iz,txym(iz,2),txym(iz,3),wtle(iz,2), &
          wtsb(iz,2),wtle(iz,3),wtsb(iz,3),iz=iz_strt,iz_end)
 5200 format(//,' IZ',5x,' SCALAR-1 MEAN',8x,' SCALAR-2 MEAN',10x, &
                ' WS1LE',10x,' WS1SB',10x,' WS2LE',10x,' WS1SB' &
             ,/,(1x,i4,6(3x,e15.6)))
      endif
 
      return
      end subroutine print

      subroutine xy_stats
!
! ------------ get statistics 
!
      implicit none
!
! ------- indices for indexing array stat(.,.)
!         js = number of non-scalar stats
!         ns = number of scalar stats
!
      integer, parameter :: js = 19, ns = 5, nstat = js + ns*nscl
      real :: stat(1:nnz, nstat)

      integer :: ix,iy,iz,izm1,izp1,izp2,m1,m2,m3,m4,m5,i,l

      real :: sgn
      real :: RHtmp, Ttmp
      real :: rhoa, msqrRH
!
! -------- stat(.,1) = u*u = ups
!          stat(.,2) = v*v = vps
!          stat(.,3) = w*w = wps
!          stat(.,4) = w**3 = wcube
!          stat(.,5) = w**4 = wfour
!          stat(.,6) = resolved tke at zw = englez
!          stat(.,7) = sgs e at zu = engsbz
!          stat(.,8) = sgs e at zw = eavg
!          stat(.,9) = resolved uw at zw = uwle
!          stat(.,10) = resolved vw at zw = vwle
!          stat(.,11) = partsrc(1) 
!          stat(.,12) = partsrc(2) 
!          stat(.,13) = partsrc(3) 
!          stat(.,14) = Tpsrc
!          stat(.,15) = Hpsrc
!          stat(.,16) = TEpsrc
!          stat(.,17) = relative humidity
!          stat(.,18) = temperature (NOT potential temp)
!          stat(.,19) = RH^2
!          stat(.,m1) = resolved scalar flux wt at zw = wtle
!          stat(.,m2) = resolved scalar flux ut at zw = utle
!          stat(.,m3) = resolved scalar flux vt at zw = vtle
!          stat(.,m4) = scalar t*t at zu = tps
!          stat(.,m5) = scalar t*t*t at zu = tcube
!
! --------- use a trick with mpi reduce over all z to get averages
!           by setting stat array = 0 for all z on each process
!
      do i=1,nstat
      do iz=1,nnz
         stat(iz,i) = 0.0
      enddo
      enddo
!
! -------- indices for scalars
!
      m1 = js
      m2 = js + nscl
      m3 = js + 2*nscl
      m4 = js + 3*nscl
      m5 = js + 4*nscl
!
      sgn = 1.0
      if(iocean .eq. 1 .and. iupwnd .eq. 1) sgn = -1.0
!
      do iz=izs,ize
!
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
         stat(iz,6) = stat(iz,6) + &
                      ((w(ix,iy,iz)-wxym(iz))**2 + &
                      (0.5*(u(ix,iy,iz)-uxym(iz) + &
                            u(ix,iy,izp1)-uxym(izp1)))**2 + &
                      (0.5*(v(ix,iy,iz)-vxym(iz) + &
                            v(ix,iy,izp1)-vxym(izp1)))**2)*0.5
!
         stat(iz,7) = stat(iz,7) + 0.5*(e(ix,iy,iz)+e(ix,iy,izm1))
         stat(iz,8) = stat(iz,8) + e(ix,iy,iz)
         stat(iz,9) = stat(iz,9) + (w(ix,iy,iz)-wxym(iz))* &
                    0.5*((u(ix,iy,iz)-uxym(iz))+ &
                         (u(ix,iy,izp1)-uxym(izp1)))
         stat(iz,10) = stat(iz,10) + (w(ix,iy,iz)-wxym(iz))* &
                    0.5*((v(ix,iy,iz)-vxym(iz))+ &
                         (v(ix,iy,izp1)-vxym(izp1)))
         stat(iz,11) = stat(iz,11) + partsrc(ix,iy,iz,1)
         stat(iz,12) = stat(iz,12) + partsrc(ix,iy,iz,2)
         stat(iz,13) = stat(iz,13) + partsrc(ix,iy,iz,3)
         stat(iz,14) = stat(iz,14) + partTsrc(ix,iy,iz)
         stat(iz,15) = stat(iz,15) + partHsrc(ix,iy,iz)
         stat(iz,16) = stat(iz,16) + partTEsrc(ix,iy,iz)

         if (iexner .eq. 1) then
            Ttmp = t(ix,iy,1,iz)* &
            exner(surf_p,func_p_base(surf_p,tsfcc(1),zz(iz)))
            !rhoa = func_rho_base(surf_p,tsfcc(1),zz(iz))
            rhoa = func_p_base(surf_p,tsfcc(1),zz(iz))/Rd/Ttmp
         else
            Ttmp = t(ix,iy,1,iz)
            rhoa = surf_rho
         end if


         RHtmp = Ttmp* &
         t(ix,iy,2,iz)*Ru/Mw/mod_magnus(Ttmp)*rhoa*100.0

         stat(iz,17) = stat(iz,17) + RHtmp
         stat(iz,18) = stat(iz,18) + Ttmp
         stat(iz,19) = stat(iz,19) + RHtmp**2

      enddo
      enddo

!
! ------------ get scalar resolved fluxes and variances
!
      do l=1,nscl
         if(iupwnd .ne. 1 .or. iz .eq. nnz) then
            do iy=iys,iye
            do ix=1,nnx
               stat(iz,m1+l)=stat(iz,m1+l) + &
                     (w(ix,iy,iz)-wxym(iz))* &
                     0.5*(t(ix,iy,l,iz)-txym(iz,l) + &
                          t(ix,iy,l,izp1)-txym(izp1,l))
            enddo
            enddo
         else
!
! ------------------- monotone fluxes
!
           do iy=iys,iye
           do ix=1,nnx
              stat(iz,m1+l) = stat(iz,m1+l) + &
          amax1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,l,iz) + &
       rlim(t(ix,iy,l,izp1),t(ix,iy,l,iz),t(ix,iy,l,izm1))) + &
          amin1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,l,izp1) + &
       rlim(t(ix,iy,l,iz),t(ix,iy,l,izp1),t(ix,iy,l,izp2)))
           enddo
           enddo
         endif
         stat(iz,m1+l)= sgn*stat(iz,m1+l)
!
! ------------ get horizontal scalar resolved fluxes 
!
         do iy=iys,iye
         do ix=1,nnx
            stat(iz,m2+l) = stat(iz,m2+l)+ &
                     (u(ix,iy,iz)-uxym(iz))* &
                     (t(ix,iy,l,iz)-txym(iz,l)) 
            stat(iz,m3+l) = stat(iz,m3+l)+ &
                     (v(ix,iy,iz)-vxym(iz))* &
                     (t(ix,iy,l,iz)-txym(iz,l)) 
         enddo
         enddo
!
! ------------------- scalar variances & higher moments
!
         do iy=iys,iye
         do ix=1,nnx
            stat(iz,m4+l) = stat(iz,m4+l) + &
                      (t(ix,iy,l,iz) - txym(iz,l))**2
            stat(iz,m5+l) = stat(iz,m5+l) + &
                      (t(ix,iy,l,iz) - txym(iz,l))**3
         enddo
         enddo
!
! ------ end scalar loop
!
      enddo
!
! ------ end z loop
!
      enddo
!
! -------- add partial sums and send it to all
!
      call mpi_sum_z(stat(1,1),i_root,myid,nnz*nstat,1)
!
! ------ fill arrays for printout and constant file
!
      do iz=1,nnz
! 
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




!
! ------------ get scalar resolved fluxes and variances
!
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

         meanRH = meanRH + &
         0.5*(zz(iz+1)-zz(iz))*(RHxym(iz)+RHxym(iz+1))

         msqrRH = msqrRH + &
         0.5*(zz(iz+1)-zz(iz))*(RHmsqr(iz)+RHmsqr(iz+1))

      end do 

      meanRH = meanRH/(zz(nnz)-zz(1))  !Note: divide by the distance between the points, not the entire zl, to get the true average
      msqrRH = msqrRH/(zz(nnz)-zz(1))
      varRH = msqrRH - meanRH**2
      

!
      return
      end subroutine xy_stats

      subroutine tke_budget

!NOTE::: THIS HAS NOT BEEN UPDATED TO REFLECT HUMIDITY IN THE BUOYANCY
!EQUATION   DHR 7/5/16

!
! -------- get terms in resolved scale tke budget
!          as in gabls writeup at w-points
!          at istage = 1. 
!
!
      implicit real(a-h,o-z), integer(i-n)
      real :: stat(1:nnz, 5)
      real :: s11s, s22s, s33s, s12s, s13s, s23s, wz, wzp
      real :: s13, s23, s33
      real :: ufluc, ufluc_t, ufluc_b, vfluc, vfluc_t, vfluc_b, wfluc
      real :: ffluc1, ffluc2, ffluc3
      real :: ffluc1p, ffluc2p, ffluc3p
      real :: weit, weit1
      integer :: iz,i,j,izp1,izm1
!
! -------- stat(.,1) = tke transport  = wq
!          stat(.,2) = pressure transport  = wp
!          stat(.,3) = tke dissipation
!          stat(.,4) = tke dissipation
!          stat(.,5) = particle force correlation
!
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
            s12s = weit1*(0.5*(uy(i,j,iz) + vx(i,j,iz)))**2 + &
                  weit*(0.5*(uy(i,j,izp1) + vx(i,j,izp1)))**2
            s13s = (0.5*((u(i,j,izp1) - u(i,j,iz) + &
                  uxym(iz) - uxym(izp1))*dzu_i(izp1) + &
                  wx(i,j,iz)))**2
            s23s = (0.5*((v(i,j,izp1) - v(i,j,iz) + &
                vxym(iz) - vxym(izp1))*dzu_i(izp1) + &
                wy(i,j,iz)))**2

         stat(iz,3) = stat(iz,3) + 2.0*vis_m(i,j,iz)*( &
                     s11s+s22s+s33s+2.0*(s12s+s13s+s23s))

            !Things for viscous transport - these are computed at u-locations since
            !a z-derivative is done to the average, which lands t_tau on w-locations
            
            ufluc_t   = u(i,j,izp1) - uxym(izp1)
            ufluc   = u(i,j,iz) - uxym(iz)
            ufluc_b   = u(i,j,izm1) - uxym(izm1)
            vfluc_t   = v(i,j,izp1) - vxym(izp1)
            vfluc   = v(i,j,iz) - vxym(iz)
            vfluc_b   = v(i,j,izm1) - vxym(izm1)
            wfluc = 0.5*( (w(i,j,iz)-wxym(iz)) &
                        + (w(i,j,izm1)-wxym(izm1)) )

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
         stat(iz,4) = stat(iz,4) + 2.0*vis_m(i,j,iz)*( &
                     ufluc*s13 + vfluc*s23 + wfluc*s33)
                    
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
         stat(iz,5) = stat(iz,5) + &
             weit*(ufluc_t*ffluc1p)+weit1*(ufluc*ffluc1)+ &
             weit*(vfluc_t*ffluc2p) + weit1*(vfluc*ffluc2) + &
             (w(i,j,iz)-wxym(iz))*(weit*ffluc3p+weit1*ffluc3)
         enddo
         enddo
         stat(iz,3) = stat(iz,3)*fnxy
         stat(iz,4) = stat(iz,4)*fnxy
         stat(iz,5) = stat(iz,5)*fnxy
      enddo
!
! --------------- get transport terms as vertical arrays
!
      do iz=izs,ize
!
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
!
! --------- get estimate of turbulent transport term
!
            ufluc   = u(ix,iy,iz) - uxym(iz)
            vfluc   = v(ix,iy,iz) - vxym(iz)
            wfluc   = w(ix,iy,iz) - wxym(iz)
            wfluc_l = w(ix,iy,izm1) - wxym(izm1)
            stat(iz,1)  = stat(iz,1) + 0.25*(wfluc + wfluc_l)* &
                   (ufluc**2 + vfluc**2 + 0.5*(wfluc_l**2 + wfluc**2))
!
! --------- get estimate of pressure transport term
!
            pfluc = p(ix,iy,iz) - pxym(iz) &
                 -0.5*((u(ix,iy,iz)+stokes(iz))**2 + &
                       v(ix,iy,iz)*v(ix,iy,iz) + &
            0.5*(w(ix,iy,iz)*w(ix,iy,iz)+w(ix,iy,izm1)*w(ix,iy,izm1)))
            stat(iz,2) = stat(iz,2) + pfluc*0.5*(wfluc_l + wfluc)
         enddo
         enddo
         stat(iz,1) = stat(iz,1)*fnxy
         stat(iz,2) = stat(iz,2)*fnxy
      enddo
      call mpi_sum_z(stat(1,1),i_root,myid,nnz*5,1)
!
! ------ we have all terms on all processors for all z, add them up
!        treat tr_tau at bottom special, tr_tau = 0
!
!      tr_tau(0) = 0.0
      do iz=1,nnz
!
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
!
! ------------- gather all the budget terms
!
         t_tran(iz)  = t_wq(iz) + t_wp(iz) + t_tau(iz)
         t_rprod(iz) = -(dudz*uwle(iz) + dvdz*vwle(iz))
         !Old t_sprod had subgrid stuff
         !t_sprod(iz) =  (dudz*uwsb(iz) + dvdz*vwsb(iz))
         !Now make t_sprod the spray tke term to reduce new variables:
         t_sprod(iz) = -stat(iz,5)
         t_buoy(iz)  =  batag*wtle(iz,1)
         t_diss(iz) = stat(iz,3)
!
      enddo
!
      return
      end subroutine tke_budget

      subroutine Tvar_budget(iscl)
      implicit none

      real :: stat(1:nnz, 4)
      real :: weit, weit1
      real :: tx_tmp(nnx, iys:iye), ty(nnx, iys:iye, izs-1:ize+1)
      real :: tx(nnx, iys:iye, izs-1:ize+1)
      real :: trans(izs:ize)
      real :: gradTp(3), Tfluc, dTpdz1, dTpdz, dTdz
      real :: Tflucp1, Tflucm1, qfluc, qflucp1, Tmean
      integer :: iz,i,j,izp1,izm1,iscl

! -------- stat(.,1) = Transport: -del.[U<T'2> + <u T'2> - alpha*grad(T'2)]
!          stat(.,2) = Dissipation: -2*alpha <grad(T').grad(T')> 
!          stat(.,3) = Particle: <T' Q'>

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

      call yd_mpi(ty(1,iys,izs-1),trigx(:,2),yk, &
                 nnx,nny,ixs,ixe,ix_s,ix_e, &
                 iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
       

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
         
         stat(iz,2) = stat(iz,2)  - 2.0*vis_s(i,j,iscl,iz)* &
                      (gradTp(1)**2+gradTp(2)**2+gradTp(3)**2)

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
            stat(iz,3) = stat(iz,3) + &
               weit*(qflucp1*Tflucp1) + weit1*(qfluc*Tfluc)
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
            stat(iz,3) = stat(iz,3) + &
               weit*(qflucp1*Tflucp1) + weit1*(qfluc*Tfluc)
          endif
          endif

          if (iscl .eq. 2 .and. iHcouple .eq. 1 ) then
              qfluc = partTEsrc(i,j,iz)-TEpsrc(iz)
              if (iz==nnz) then
              qflucp1 = 0.0
              else
              qflucp1 = partTEsrc(i,j,izp1)-TEpsrc(izp1)
              end if
              Tv_part2(iz) = Tv_part2(iz) + &
                weit*(qflucp1*Tflucp1) + weit1*(qfluc*Tfluc)
  
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


!
! ------ we have all terms on all processors for all z, add them up
!
      do iz=1,nnz
         izp1 = iz + 1
         izm1 = iz - 1
         if(iz .eq. nnz) then
            Tv_tran(iz,iscl) = 0.0
         else
            Tv_tran(iz,iscl) = -(stat(izp1,1)-stat(iz,1))*dzu_i(izp1)
         endif
         dTdz = (txym(izp1,iscl)-txym(iz,iscl))*dzu_i(izp1)
!
! ------------- gather all the budget terms
!
         Tv_prod(iz,iscl) = -2.0*wtle(iz,iscl)*dTdz
         Tv_diss(iz,iscl) = stat(iz,2)
         Tv_part1(iz,iscl) = -stat(iz,3)
         !Tv_part2 is already gathered
      enddo
!
      return
      end subroutine Tvar_budget

      subroutine write_his(iloc)
!
! ----- write history file with global parameters
!       write tsfcc specially to preserve digits!
!
!
      implicit real(a-h,o-z), integer(i-n)
      divgmax = 0.0
      do iz=1,nnz
         divgmax = amax1(divgmax, divz(iz))
      enddo
!
      ziavg = zi
      holtop = hol
      wt_min = wtsb(iloc,1)
      wt_le  = wtle(iloc,1)
      krec = krec + 1
      mid = nnz/4
      write(nhis1,6000) time,dt,utau,ziavg,amonin,holtop, &
               (tsfcc(1)-t_ref),uwsfc,vwsfc,divgmax, wt_min, wt_le, &
               ucfl, vcfl, wcfl, wtsfc(1), &
               ups(mid),vps(mid),wps(mid),tps(mid,1), &
               uwle(mid),uwsb(mid),uw_tot(mid), &
               vwle(mid),vwsb(mid),vw_tot(mid), &
               wtle(mid,1),wtsb(mid,1),wt_tot(mid,1), &
               englez(mid),eavg(mid), wabs, &
               lwc, &
               Rep_avg,radavg,radmin,radmax

 6000 format(37e17.8)
!
! -------------- write profile information
!
      call write_prof(nhisp,krec,isize,c_s%wwsb)
!
      return
      end subroutine write_his

      subroutine write_prof(nhisp,krec,num,f)
      implicit real(a-h,o-z), integer(i-n)
      real f(num)
      real*4 f32(num)
!
! -------------- build special 32 bit arrays for profiles
!
      do i=1,num
         f32(i) = f(i)
      enddo
!
      write(nhisp,err=999,rec=krec) (f32(i),i=1,num)
!
      return
! --------------- errors
  999 continue
      write(6,9000) num,krec
 9000 format(' 9000, trouble in ', &
             'SR. save_prof cannot write profile data ',/, &
             ' num = ',i8, 'krec = ',i6)
      stop
      end subroutine write_prof

      subroutine close_histograms
!
! ---- close history files
!
      implicit real(a-h,o-z), integer(i-n)
      logical there
!
! ---- root closes and checks the files
!
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
!
      return
! -------------------- process write errors
 7100 format(' RADHIST DATA IS WRITTEN IN FILE  ',a80)
 7200 format(' RESHIST DATA IS WRITTEN IN FILE  ',a80)
 7300 format(' ACTRESHIST DATA IS WRITTEN IN FILE  ',a80)
 8100 format(' SR. SAVE_RADPDF: Trouble rad histogram file =',a80)
 8200 format(' SR. SAVE_RESPDF: Trouble res histogram file =',a80)
 8300 format(' SR. SAVE_ACTRESPDF: Trouble res histogram file =',a80)
      end subroutine close_histograms

      subroutine close_his
!
! ---- close history files
!
      implicit real(a-h,o-z), integer(i-n)
      logical there
!
! ---- root closes and checks the files
!
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
!
      return
! -------------------- process write errors
 7000 format(' HISTORY DATA IS WRITTEN IN FILE  ',a80)
 7100 format(' PROFILE HISTORY DATA IS WRITTEN IN FILE  ',a80)
 7200 format(' RADHIST DATA IS WRITTEN IN FILE  ',a80)
 8000 format(' SR. SAVE_HIS: Trouble history file not in path =',a80)
 8100 format(' SR. SAVE_HIS: Trouble profile history file', &
             ' not in path =',a80)
 8200 format(' SR. SAVE_RADPDF: Trouble rad histogram file =',a80)
      end subroutine close_his

      subroutine save_viz
!
! --------------- save multiple (x-y), (x-z), (y-z), planes of data .
!                 modify recl in all open statements for more or less
!                 variables. 
!                 Constant - x, implies yz planes
!                 Constant - y, implies xz planes
!                 Constant - z, implies xy planes
!
! ------------- routine uses send/recv to get information in y-z planes
!
#if defined(SWAP)
      use module_byteswap
#endif
!
      parameter(nvar_o = 6)
!
      integer ix_pick(maxnx),  iy_pick(maxny),  iz_pick(maxnz), &
      implicit real(a-h,o-z), integer(i-n)
              ix_order(maxnx), iy_order(maxny), iz_order(maxnz)
!
      integer istatus(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
!
      real(kind=4), dimension(nvar_o,nny,izs:ize) :: temp_x
      real(kind=4), dimension(nvar_o,nnx,izs:ize) :: temp_y
      real(kind=4), dimension(nvar_o,nnx,iys:iye) :: temp_z
      real, dimension(nvar_o,iys:iye,izs:ize)     :: buf_send
!
! ------------- don't touch
!
      data iviz_x,  iviz_y,  iviz_z  /0, 0, 0/
      data ionce_x, ionce_y, ionce_z, istuff /0, 0, 0, 0/
      data ix_pick, iy_pick, iz_pick /maxnx*0, maxny*0, maxnz*0/
      data ix_order, iy_order, iz_order /maxnx*0, maxny*0, maxnz*0/
      save iviz_x,  iviz_y,  iviz_z, &
           ix_pick, iy_pick, iz_pick, &
           ix_order, iy_order, iz_order, &
           ionce_x, ionce_y, ionce_z, istuff, &
           npln_x, npln_y, npln_z
!
!
! ----------- turn on z levels to save. Customize for your own use.
!             Set iz_pick(iz) = iz, ix_pick(ix) = ix, iy_pick(iy) = iy
!             Data is round-robin alternated in the data file for more than
!             1 plane for any particular view.
!

      iz_pick(12) = 12
      iz_pick(28) = 28
      iz_pick(64) = 64
      !iz_pick(20) = 20
      !iz_pick(45) = 45
      !iz_pick(60) = 60
!
! -------------- pick an x-z plane of data (can add more)
!
      iy_pick(nny/2) = nny/2
!     iy_pick(nny)   = nny
!
! -------------- pick a y-z plane of data (can add more)
!
      ix_pick(nnx/2) = nnx/2
!     ix_pick(nnx)   = nnx
!
! ------ find total number of z's turned on and open file once
!
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
            call mpi_file_open(mpi_comm_world, path_viz_xy, &
                               mpi_mode_create+mpi_mode_rdwr, &
                               mpi_info_null, nviz_z, ierr)
            disp = 0
            call mpi_file_set_view(nviz_z,disp,mpi_real4,mpi_real4, &
                                  'native',mpi_info_null,ierr)
         endif
      endif
!
! ------ find total number of y's turned on and open file once
!
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
            call mpi_file_open(mpi_comm_world, path_viz_xz, &
                               mpi_mode_create+mpi_mode_rdwr, &
                               mpi_info_null, nviz_y, ierr)
            disp = 0
            call mpi_file_set_view(nviz_y,disp,mpi_real4,mpi_real4, &
                                  'native',mpi_info_null,ierr)
         endif
      endif
!
! ------ find total number of x's turned on and open file once
!
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
            call mpi_file_open(mpi_comm_world, path_viz_yz, &
                               mpi_mode_create+mpi_mode_rdwr, &
                               mpi_info_null, nviz_x, ierr)
            disp = 0
            call mpi_file_set_view(nviz_x,disp,mpi_real4,mpi_real4, &
                                  'native',mpi_info_null,ierr)
         endif
      endif
!
      if(istuff .eq. 0 .and. l_root) then
         open(nviz_s,file=path_stuf)
         istuff = 1
      endif
!
! --------- write data, subtract t_ref to increase
!           resolution on 32 bit machines
!
! ---------- xy planes of data
!
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
!
! ---------- get just the fluctuating pressure field
!
               temp_z(5,i,j) = p(i,j,k) - pxym(k) &
                              -(e(i,j,k) + e(i,j,km1))/3.0 &
                              -0.5*((u(i,j,k) + stokes(k))**2 + &
                                     v(i,j,k)*v(i,j,k) + &
                               0.5*(w(i,j,k)*w(i,j,k) + &
                                    w(i,j,km1)*w(i,j,km1)))

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
            call mpi_file_write_at(nviz_z,offset,temp_z,n_write, &
                                   mpi_real4,istatus,ierr)
            if (ierr .ne. 0) go to 9991
         endif
      enddo
!
! ---------- xz planes of data
!
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
!
! ---------- get just the fluctuating pressure field
!
               temp_y(5,i,k) =  p(i,j,k) - pxym(k) &
                                -(e(i,j,k)+e(i,j,km1))/3.0 &
                                -0.5*((u(i,j,k)+stokes(k))**2 + &
                                     v(i,j,k)*v(i,j,k) + &
                                 0.5*(w(i,j,k)*w(i,j,k) + &
                                      w(i,j,km1)*w(i,j,km1)))

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
            call mpi_file_write_at(nviz_y,offset,temp_y,nwrite, &
                                      mpi_real4,istatus,ierr)
            if (ierr .ne. 0) goto 9992
         endif
      enddo
!
! ---------- yz planes that cut across all processors
!            just have root node on that slab write data
!
      iviz_x  = iviz_x + npln_x
      n_write = nvar_o*nny*(ize+1-izs)
      nsize   = int(nvar_o,k8)*nny*nnz
      nsize2  = int(nvar_o,k8)*nny*(izs-1)
      n_send  = nvar_o*(ize+1-izs)*(iye+1-iys)
      do i=1,nnx
         if(ix_pick(i) .eq. i) then
!
! ----------- build send buffer
!
            do k=izs,ize
            km1 = k - 1
            do j=iys,iye
               buf_send(1,j,k) = u(i,j,k)
               buf_send(2,j,k) = v(i,j,k)
               buf_send(3,j,k) = w(i,j,k)
               buf_send(4,j,k) = (t(i,j,1,k) - t_ref)
!
! ---------- get just the fluctuating pressure field
!
               buf_send(5,j,k) = p(i,j,k) - pxym(k) &
                                -(e(i,j,k) + e(i,j,km1))/3.0 &
                                -0.5*((u(i,j,k) + stokes(k))**2 + &
                                     v(i,j,k)*v(i,j,k) + &
                                 0.5*(w(i,j,k)*w(i,j,k) + &
                                      w(i,j,km1)*w(i,j,km1)))
              buf_send(6,j,k) = u(i,j,k)-uxym(k)
              !buf_send(6,j,k) = partsrc(i,j,k,1)
              !buf_send(7,j,k) = partsrc(i,j,k,2)
              !buf_send(8,j,k) = partsrc(i,j,k,3)
            enddo
            enddo
            if(myid .ne. iss) then
              call mpi_send(buf_send(1,iys,izs),n_send, &
                            mpi_real8,iss,1, &
                            mpi_comm_world,ierr)
            else
              do k=izs,ize
              do j=iys,iye
              do ii=1,nvar_o
                 temp_x(ii,j,k) = buf_send(ii,j,k)
              enddo
              enddo
              enddo
              do l=iss+1,ise
                 call recv_yz_var(temp_x,nvar_o,nny, &
                                  iy_s(l),iy_e(l),izs,ize,l)
              enddo
#if defined(SWAP)
              call byteswap(temp_x)
#endif
              offset = int((iviz_x + ix_order(i) - 1),k8)*nsize + nsize2
              call mpi_file_write_at(nviz_x,offset,temp_x,n_write, &
                                mpi_real4,istatus,ierr)
              if (ierr .ne. 0) goto 9993
            endif
         endif
      enddo
!
! ------------- ascii file with facts in it that goes
!               with visualization
!
      if(l_root) then
         write(nviz_s,5000) time, amonin, zi, utau
 5000    format(4e20.8)
      endif
!
! ---- last time step close the files
!
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
!
      return
! --------------------------  errors in writing viz file
 9991 continue
      write(6,6000) nviz_z, iz
 6000 format(' SR. SAVE_VIS:',/, &
             '    trouble cannot write xy viz file on unit = ',i2,/, &
             '             at iz = ',i4)
      call mpi_finalize(ierr)
      stop
! --------------------------  errors in writing viz file
 9992 continue
      write(6,6100) nviz_y, iz, iviz_y
 6100 format(' SR. SAVE_VIS:',/, &
             '    trouble cannot write xz viz file on unit = ',i2,/, &
             '             at iz = ',i4,/, &
             '            iviz_y = ',i8)
! --------------------------  errors in writing viz file
 9993 continue
      write(6,6200) nviz_x, iz, iviz_x
 6200 format(' SR. SAVE_VIS:',/, &
             '    trouble cannot write yz viz file on unit = ',i2,/, &
             '             at iz = ',i4,/, &
             '            iviz_x = ',i8)
      call mpi_finalize(ierr)
      stop
      end subroutine save_viz

      subroutine recv_yz_var(temp_x,nvar,nny,iys,iye,izs,ize,ir)
!
      implicit real(a-h,o-z), integer(i-n)
      integer istatus(mpi_status_size)
!
      real buf(nvar,iys:iye,izs:ize)
      real(kind=4), dimension(nvar,nny,izs:ize) :: temp_x
!
      num = nvar*(ize+1-izs)*(iye+1-iys)
      call mpi_recv(buf(1,iys,izs),num,mpi_real8,ir,1, &
                   mpi_comm_world,istatus,ierr)
      do k=izs,ize
      do j=iys,iye
      do ii=1,nvar
         temp_x(ii,j,k) = buf(ii,j,k)
      enddo
      enddo
      enddo
!
      return
      end subroutine recv_yz_var

      subroutine save_v
!
! --------------- save 3d fields
!
#if defined(SWAP)
      use module_byteswap
#endif
      logical there
!
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)                 nsize, nsize2
!
      real, allocatable, dimension(:,:,:) :: temp
      implicit real(a-h,o-z), integer(i-n)
      allocate(temp(nvar,nnx,iys:iye))
!
! ---- open file
!
      call mpi_file_open(mpi_comm_world, path_sav_v, &
                         mpi_mode_create+mpi_mode_rdwr, &
                         mpi_info_null, nvel, ierr)
!
! ---- set file view
!
      disp = 0
      call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8, &
                            'native',mpi_info_null,ierr)
!
! ---- write data
!
      nsize   = int(nvar,k8)*nnx*nny
      nsize2  = int(nvar,k8)*nnx*(iys-1)
      n_write = nvar*nnx*(iye+1-iys)
!
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
!

         offset = int((k-1),k8)*nsize + nsize2
!         call mpi_file_write_at(nvel,offset,temp,n_write,
!     +                              mpi_real8,status,ierr)
         call mpi_file_write_at_all(nvel,offset,temp,n_write, &
                                    mpi_real8,status,ierr)
         if (ierr /= 0) goto 9991
!
      enddo

!
! ---- close file
!
      call mpi_file_close(nvel, ierr)

!
! ---- check file
!
      if (l_root) then
         inquire(file=path_sav_v,exist=there)
         if(.not.there) then
            write(6,8000) nvel,myid
            call mpi_finalize(ierr)
            stop
         endif
         write(6,7000) it,path_sav_v
      endif
!
      deallocate(temp)
!
      return
! --------------------------  errors in writing restart file
 9991 continue
      write(6,6000) nvel, iz
 6000 format(' SR. SAVE_V:',/, &
             '    trouble cannot write restart file on unit = ',i2,/, &
             '             at iz = ',i4)
      call mpi_finalize(ierr)
      stop
! --------------------
 7000 format(' **** DATA SET AT IT = ',I6,/, &
             '      VELOCITY DATA IS WRITTEN IN FILE  ',a80)
 8000 format(' in SAVE_V: trouble writing file ',i5,'  myid = ',i5, &
             ' at iz = ',i5)
      end subroutine save_v

      subroutine save_c
!
! --------------- root process writes constant file
!                 sequential fortan binary
!
!
      implicit real(a-h,o-z), integer(i-n)
      logical there
      character options*8, passwd*1
!
! ---- open file
!
      open(nvelc,err=9992,file=path_sav_c,form='unformatted', &
                      status='unknown')
      write(nvelc,err=9992) c_c, c_s, c_i, case
      close(nvelc)
!
        inquire(file=path_sav_c,exist=there)
        if(.not.there) then
           write(6,8001) path_sav_c
           call mpi_finalize(ierr)
           stop
        endif
! -----------------------------  output ok message
      write(6,7001) path_sav_c
!
      return
! --------------------------  errors in writing constant file
 9992 continue
      write(6,6100) nvelc
 6100 format(' SR. SAVE_V:',/, &
        '    trouble cannot open/write constant file on unit = ',i2)
      call mpi_finalize(ierr)
      stop
! ---------------------
 7001 format('      CONSTANT DATA IS WRITTEN IN FILE  ',a80)
 8001 format(' SR. SAVE_C: Trouble constant file not in path =',a80)
      end subroutine save_c

      subroutine save_p
!
! -------------- save pressure file
!
#if defined(SWAP)
      use module_byteswap
#endif
      logical there
!
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
!
      real temp(nnx,iys:iye)
!
! ---- open file
!
      implicit real(a-h,o-z), integer(i-n)
      call mpi_file_open(mpi_comm_world, path_sav_p, &
                         mpi_mode_create+mpi_mode_rdwr, &
                         mpi_info_null, npre, ierr)
!
! ---- set file view
!
      disp = 0
      call mpi_file_set_view(npre,disp,mpi_real8,mpi_real8, &
                            'native',mpi_info_null,ierr)
!
! ---- write data
!
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
         call mpi_file_write_at_all(npre,offset,temp,n_write, &
                                    mpi_real8,status,ierr)
!         call mpi_file_write_at(npre,offset,temp,n_write,
!     +                              mpi_real8,status,ierr)
      enddo
!
! ---- close file
!
      call mpi_file_close(npre, ierr)
!
! ---- check file
!
      if (l_root) then
         inquire(file=path_sav_p,exist=there)
         if(.not.there) then
            write(6,8000) path_sav_p
            call mpi_finalize(ierr)
            stop
         endif
         write(6,7000) path_sav_p
      endif
!
      return
! -------------------- process write errors
 9991 continue
      write(6,6000) npre, iz
 6000 format(' SR. SAVE_P:',/, &
             '    trouble cannot write pressure file on unit = ',i2,/, &
             '             at iz = ',i4)
      call mpi_finalize(ierr)
      stop
! -----------------------
 7000 format('      PRESSURE DATA IS WRITTEN IN FILE  ',a80)
 8000 format(' SR. SAVE_P: Trouble pressure file not in path =',a80)
      end subroutine save_p

      subroutine get_units
!
!
! -------------- unit numbers for files
!
      implicit real(a-h,o-z), integer(i-n)
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
!
! ------------- unit number for standard print out
!               for each mpi task
!
      nprt = 1 
!
! ------------- open unit for standard printout
!
      path_prt = case_inp(1:3)//'.mp.xxxxx.out'
      write(path_prt(8:12),'(i5.5)') myid
      open(nprt,file=path_prt,form='formatted')
!
      return
      end subroutine get_units

      subroutine get_output_filenames
!
! ----------- build file names for velocity, pressure, and constants
!
      implicit none
      character cgrid*10, num*3
      integer iblnk
!
! --------------- build character strings for file name
!
      cgrid = '.le.'
      write(num,'(i3.3)') itn
      iblnk = index(path_seed,' ')
      call blnk(path_sav_v)
      call blnk(path_sav_p)
      call blnk(path_sav_c)
      call blnk(path_sav_part)
      path_sav_v = path_seed(1:iblnk-1)//'/u'// &
                       cgrid(1:4)//case(1:3)//num(1:3)
      path_sav_p = path_seed(1:iblnk-1)//'/p'// &
                       cgrid(1:4)//case(1:3)//num(1:3)
      path_sav_c = path_seed(1:iblnk-1)//'/u'// &
                       cgrid(1:4)//case(1:3)//num(1:3)//'.con'
      path_sav_part = path_seed(1:iblnk-1)//'/part'// &
                       cgrid(1:4)//case(1:3)//num(1:3)
!
      return
      end subroutine get_output_filenames

      subroutine open_histograms(istep)
!
! ------------------- open history files by root
!                     isize determined in sr. fill_cs
!
      implicit real(a-h,o-z), integer(i-n)
      character cgrid*4, iblks*16, num*3
      logical there
     
!
! --------------- build character strings for ascii history file name
!
      !Open PDF files
      cgrid = '.le.'

      call blnk(iblks)
      write(iblks(1:7),'(i7.7)') istep
      iblks(8:8) = '_'
      write(iblks(9:15),'(i7.7)') (istep+itape)
      iblnk = index(path_seed,' ')
      call blnk(path_rad_hist)
      call blnk(path_res_hist)


      path_rad_hist = path_seed(1:iblnk-1)//'/radhist'// &
             cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'      

      path_res_hist = path_seed(1:iblnk-1)//'/reshist'// &
             cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'      

      path_actres_hist = path_seed(1:iblnk-1)//'/actreshist'// &
             cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'      


!
! ----------------- save data in directory
!
      if(l_root) then


      close(nrad)
      open(nrad,err=4500,file=path_rad_hist,form='formatted')

      close(nres)
      open(nres,err=4600,file=path_res_hist,form='formatted')

      close(nactres)
      open(nactres,err=4700,file=path_actres_hist,form='formatted')

      endif
!
      return
! ------------------- process errors
!-------------------
 4500 continue
      write(6,6303) nrad, path_rad_hist
 6303 format(' 6303, SR. OPEN_RADPDF:',/, &
             '    cannot open radhist file on unit = ',i2,/, &
             '    path = ',a80)
      stop
 4600 continue
      write(6,6303) nres, path_res_hist
 6304 format(' 6303, SR. OPEN_RESPDF:',/, &
             '    cannot open reshist file on unit = ',i2,/, &
             '    path = ',a80)
      stop
 4700 continue
      write(6,6305) nactres, path_actres_hist
 6305 format(' 6305, SR. OPEN_ACTRESPDF:',/, &
             '    cannot open actreshist file on unit = ',i2,/, &
             '    path = ',a80)
      stop
      end subroutine open_histograms

      subroutine open_his(istep)
!
! ------------------- open history files by root
!                     isize determined in sr. fill_cs
!
      implicit real(a-h,o-z), integer(i-n)
      character cgrid*4, iblks*16, num*3
      logical there
     
!
! --------------- build character strings for ascii history file name
!
      cgrid = '.le.'
      call blnk(iblks)
      write(iblks(1:7),'(i7.7)') istep
      iblks(8:8) = '_'
      write(iblks(9:15),'(i7.7)') (istep + itape)
      iblnk = index(path_seed,' ')
      call blnk(path_sav_h)
      path_sav_h = path_seed(1:iblnk-1)//'/his'// &
               cgrid(1:4)//case(1:3)//'.'//iblks(1:15)
!
! --------------- build character strings for ieee profile history file
!                 set record counter for direct access file = 0
!
      krec = 0
      cgrid = '.le.'
      call blnk(iblks)
      write(iblks(1:7),'(i7.7)') istep
      iblks(8:8) = '_'
      write(iblks(9:15),'(i7.7)') (istep + itape)
      iblnk = index(path_seed,' ')
      call blnk(path_sav_hp)
      path_sav_hp = path_seed(1:iblnk-1)//'/his'// &
               cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.ieee'

!
! ----------------- save data in directory
!
      if(l_root) then

      close(nhis1)
      open(nhis1,err=3000,file=path_sav_h,form='formatted')
!
      close(nhisp)
      open(nhisp,err=4000,file=path_sav_hp, &
              form='unformatted',access='direct',recl=isize*j_recl, &
              status='unknown')

      end if

!
      return
! ------------------- process errors
 3000 continue
      write(6,6301) nhis1, path_sav_h
 6301 format(' 6301, SR. OPEN_HIS:',/, &
             '    cannot open history1 file on unit = ',i2,/, &
             '    path = ',a80)
      stop
!-------------------
 4000 continue
      write(6,6302) nhisp, path_sav_hp
 6302 format(' 6302, SR. OPEN_HIS:',/, &
             '    cannot open history profile file on unit = ',i2,/, &
             '    path = ',a80)
      stop

      end subroutine open_his

      subroutine viz_output_filename(istep)
!
! ------------------- set visualization files,
!                     leaves files in scratch directory 
!
      implicit real(a-h,o-z), integer(i-n)
      character iblks*16

!
! --------------- build character strings for file names
!                 with time step
!
      call blnk(iblks)
      iblks(1:1) = '.'
      write(iblks(2:8),'(i7.7)') istep
      iblks(9:9) = '_'
      write(iblks(10:16),'(i7.7)') (istep + itape)
!
      iloc = index(path_seed,' ')
      path_viz_xy = path_seed(1:iloc-1) &
               //'viz.'//case(1:3)//iblks(1:16)//'.xy.data'
!
      path_viz_xz = path_seed(1:iloc-1) &
               //'viz.'//case(1:3)//iblks(1:16)//'.xz.data'
!
      path_viz_yz = path_seed(1:iloc-1) &
               //'viz.'//case(1:3)//iblks(1:16)//'.yz.data'
!
      path_stuf = path_seed(1:iloc-1) &
               //'stuff.'//case(1:3)//iblks(1:16)//'.data'
!
!     if(l_root) then
!        write(6,8001) path_viz_xy
!8001    format(' 8001: viz file = ',a80)
!        write(6,8001) path_viz_xz
!        write(6,8001) path_viz_yz
!        write(6,8001) path_stuf
!        write(6,8001) path_seed
!     endif
!
      return
      end subroutine viz_output_filename

      subroutine open_viz
!
! ------------------- set visualization files,
!                     leaves files in scratch directory 
!
      implicit real(a-h,o-z), integer(i-n)
      character iblks*16
!
! --------------- build character strings for file names
!                 with time step
!
      call blnk(iblks)
      iblks(1:1) = '.'
      write(iblks(2:8),'(i7.7)') iti
      iblks(9:9) = '_'
      write(iblks(10:16),'(i7.7)') itmax
!
      iloc = index(path_viz_xy,' ')
      path_viz_xy = path_viz_xy(1:iloc-1) &
            //'/viz.'//case(1:3)//iblks(1:16)//'.xy.data'
      iloc = index(path_viz_xz,' ')
      path_viz_xz = path_viz_xz(1:iloc-1) &
            //'/viz.'//case(1:3)//iblks(1:16)//'.xz.data'
      iloc = index(path_viz_yz,' ')
      path_viz_yz = path_viz_yz(1:iloc-1) &
            //'/viz.'//case(1:3)//iblks(1:16)//'.yz.data'
      path_stuf = path_stuf(1:iloc-1) &
            //'/stuff.'//case(1:3)//iblks(1:16)//'.data'
      close(nviz_z)
      close(nviz_y)
      close(nviz_x)
      close(nviz_s)
!
! ----------- do not actually open the files here since
!             not all processors may have been picked and
!             its unknown how many variables are selected.
!             customized in sr. save_viz 
!
      return
      end subroutine open_viz

end module mod_io
