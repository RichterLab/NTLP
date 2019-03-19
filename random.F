      subroutine random
c
c ----------- geostrophic winds designed for comparison case
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
      real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +     psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +     vyy(nnx,iys:iye,izs:izs) ,T_real

      !Initialize partcount to 0:
      partcount     = 0.0
      partcount_t   = 0.0
c     --Particle flux JP--
      partflux      = 0.0
      partflux_t    = 0.0
      partcountU    = 0.0
      partcountU_t  = 0.0
      pcountrange   = 0.0
      pcountrange_t = 0.0
      Spartflux_t   = 0.0
      Spartflux     = 0.0
      tpartflux_t   = 0.0
      tpartflux     = 0.0
      tupartflux_t  = 0.0
      tupartflux    = 0.0
      tdpartflux_t  = 0.0
      tdpartflux    = 0.0
      xytpartflux_t = 0.0
      xytpartflux   = 0.0
      xyutpartflux_t= 0.0
      xydtpartflux  = 0.0
c     --------------------
      vpsum         = 0.0
      vpsum_t       = 0.0
      vpsqrsum      = 0.0
      vpsqrsum_t    = 0.0
      upwp_t        = 0.0
      upwp          = 0.0
      Tpsum         = 0.0
      Tpsum_t       = 0.0
      Tpsqrsum      = 0.0
      Tpsqrsum_t    = 0.0
      Tfsum         = 0.0
      Tfsum_t       = 0.0
      qfsum         = 0.0
      qfsum_t       = 0.0
      wpTpsum       = 0.0
      wpTpsum_t     = 0.0
      partsrc       = 0.0
      partsrc_t     = 0.0
      partTsrc      = 0.0
      partTsrc_t    = 0.0
      partHsrc      = 0.0
      partHsrc_t    = 0.0
      partTEsrc     = 0.0
      partTEsrc_t   = 0.0
      radsum        = 0.0 
      radsum_t      = 0.0 
      rad2sum       = 0.0 
      rad2sum_t     = 0.0 
      mpsum         = 0.0 
      mpsum_t       = 0.0 
      mwsum         = 0.0 
      mwsum_t       = 0.0 
      qstarsum      = 0.0
      qstarsum_t    = 0.0 

c ------------ note set nmatch in sr. iso so that
c              it is compatible with conditions here
c
      do iz=1,nnz
c        ug(iz)   = ugcont*(zz(iz)/zl)
         ug(iz)   = ugcont
         vg(iz)   = vgcont
         divz(iz) = 0.0
      enddo
c
c      izi = (100*nnz)/100  !initiall 50,changed 100 by Indrajith 04-06/17
      izi = (50*nnz)/100 
      zi  = z(izi)
c
      z_lower = zi - 50.0
      t_lower = 283.15
      z_upper = zi + 50.0
      ! intially t_upper =308, changed to make stronger inversion,
      ! Indrajith 09/08/2017
      t_upper = 310.0
      slope   = (t_upper - t_lower)/(z_upper - z_lower)
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = ugcont-ugal
c           u(ix,iy,iz) = ug(iz) - ugal
            v(ix,iy,iz) = vgcont
            w(ix,iy,iz) = 0.0
            e(ix,iy,iz) = 0.0
         enddo
         enddo
c-----------------------------------------------------------------
c       Commented by Indrajith 09-22-17 for MABL-concentration
c       comparison
c         if(z(iz) .le. z_lower) then
c           do iy=iys,iye
c           do ix=1,nnx
c             t(ix,iy,1,iz) = t_lower
c              !t(ix,iy,1,iz) = tsfcc(1) - 5.0 !Brian 5/28/14

c           enddo
c           enddo
c         elseif(z(iz) .ge. z_upper) then
c           do iy=iys,iye
c           do ix=1,nnx
c              t(ix,iy,1,iz) = t_upper + (zz(iz+1) - z_upper)*dtdzf(1)
c              !t(ix,iy,1,iz) = tsfcc(1) - 5.0
c           enddo
c           enddo
c         else
c           do iy=iys,iye
c           do ix=1,nnx
c              t(ix,iy,1,iz) = t_lower + slope*(zz(iz+1) - z_lower)
              !t(ix,iy,1,iz) = tsfcc(1) - 5.0 
c           enddo
c           enddo
c         endif
c------------------------------------------------------------------
c       Added by Indrajith 09/22/17 for MABL-concentration comparison
c        
        if(z(iz) .le. zi) then
           do iy=iys,iye
           do ix=1,nnx
             t(ix,iy,1,iz) = t_lower
           enddo
           enddo
        else
           do iy=iys,iye
           do ix=1,nnx
              t(ix,iy,1,iz) = t_lower + (zz(iz+1) - zi)*dtdzf(1)
           enddo
           enddo
        end if
c-----------------------------------------------------------------
         do iy=iys,iye
         do ix=1,nnx
            !t(ix,iy,2,iz) = 0.0  !For now just set the initial humidity field to zero: Changed to 80% RH by Indrajith 01/02/17
            T_real = t(ix,iy,1,iz)*(1.0-9.81*z(iz)/(Cpa*tsfcc(1)))    ! real temperature
            if (z(iz).le.z_lower) then
            t(ix,iy,2,iz) =  0.8*Mw/Ru/T_real*610.94*
     +               EXP(17.6257*(T_real-273.15)/(T_real-
     +               273.15+243.04))/rhoa       ! RH set to 80%
            else
              
           t(ix,iy,2,iz) =  0.05*Mw/Ru/T_real*610.94*
     +               EXP(17.6257*(T_real-273.15)/(T_real-
     +               273.15+243.04))/rhoa       ! RH set to 5%
          ! t(ix,iy,2,iz) = 0.0 ! make the RH 0 for check , Indrajith
            end if
          ! t(ix,iy,2,iz) = 0.0
           ! if (iy .eq. 64 .and. ix .eq.64) then
           ! write(*,*) 'D5EBUG:',T_real,t(ix,iy,2,iz),Mw,Ru,rhoa,Cpa

           ! end if 
            w(ix,iy,iz)   = 0.
            r1(ix,iy,iz)  = 0.
            r2(ix,iy,iz)  = 0.
            r3(ix,iy,iz)  = 0.
            r4(ix,iy,1,iz)= 0.
            !r4(ix,iy,2,iz)= 0.
            r5(ix,iy,iz)  = 0.
         enddo 
         enddo 
      enddo
c
c ------------- set initial random field to be
c               divergence free
c
      idum = -1 - myid
      do iz=izs,ize
c
c ----------- ampv and ampt are max amplitudes of random 
c             velocity and temperature fields
c             make sure ampv is set if free convection so
c             that we have motions at first time step
c
         ampv = 0.0
         ampv = 0.001
         ampt = 0.10
c  
c ------- simple random field scaled between -0.5 and 0.5
c
         sum_psi = 0.0
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = ran1(idum)
            sum_psi = sum_psi + psi(ix,iy)
         enddo
         enddo
         sum_psi = sum_psi*fnxy
         call mpi_sum_xy(sum_psi,myid,iss,ise,1)
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = psi(ix,iy) - sum_psi
            psix(ix,iy)     = psi(ix,iy)
            psiy(ix,iy,izs) = psi(ix,iy)
         enddo
         enddo
         call xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         vmaxx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
            if(vmag .gt. vmaxx) vmaxx = vmag
         enddo
         enddo
         facv = ampv/vmaxx
c
         if (z(iz) .le. 50.0) then
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz)   = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
            v(ix,iy,iz)   = v(ix,iy,iz) + psix(ix,iy)*facv
            t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
         enddo
         enddo
         endif
c
         if(z(iz) .le. 250.0) then
         do iy=iys,iye
         do ix=1,nnx
            e(ix,iy,iz) = 0.4*(1.0 - z(iz)/250.0)**3
         enddo
         enddo
         endif
c
c ---------- check divergence of initial field
c
         do iy=iys,iye
         do ix=1,nnx
            uxx(ix,iy) = u(ix,iy,iz)
            vyy(ix,iy,izs) = v(ix,iy,iz)
         enddo
         enddo
         call xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         do iy=iys,iye
         do ix=1,nnx
            divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
         enddo
         enddo
         divz(iz) = divz(iz)*fnxy
c
c -------- end z loop
c
      enddo
c
      call mpi_sum_z(divz(1),i_root,myid,nnz,1)
c
      write(nprt,6000)
 6000 format(' check of divergence for initial state',/,
     +         ' iz ',5x,' divergence')
      write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100 format(i5,e15.6)
c
c ------------ fix for baroclinic and subsidence effects !!
c
c     do iz=izs,ize
c        ug(iz)=ugcont
c        vg(iz)=vgcont
c        if (.not.(ibrcl.eq.1)) go to 19988
c        if (.not.(iz.le.izi)) go to 19987
c        ug(iz)=0.
c        vg(iz)=0.
c 19987    continue
c 19988    continue
c        zz2=zz(iz)
c        wls(iz)=-divgls*zz2
c        if (.not.(iz.eq.1)) go to 19986
c        do ix=1,nnx
c        uls(ix)=divgls*(dx*float(ix-1)-xl*.5)
c        enddo
c     enddo
c     write(nprt,9)(uls(ix),ix=1,nnx)
c  9  format(1x,8e12.3)
c 19986 continue
c
      return
      end
