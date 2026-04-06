module mod_init
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
  use particles
  implicit none
  private
  public :: setcon, setup, init, vgrid, vgrid_channel, vgrid_channel_fstrm, vgrid_uniform, get_dz, random, random_f, randoc, gridd, restart, read_res
contains

! --- extracted from les.F: setcon ---
      subroutine setcon
      implicit real(a-h,o-z), integer(i-n)
!
!
! ----------------- get machine type, can erad datadir also
!
      open(99,file='mach.file')
      read(99,9000) imach
 9000 format(i1)
!
      close(99)
!
! ----------- initialize fft
!
      call rffti(nnx,trigx(:,1))
      call rffti(nny,trigx(:,2))
      call cffti(nny,trigc(1))
!
! ----------- start step for history files
!
      if (iti.eq.0) then
         it_his = 1
         it_viz = 1
      else
         it_his = iti+2
         it_viz = iti+2
      endif

      it_his_nxt = it_his
      it_viz_nxt = it_viz
!
! ---------------- set min value of e
!
      if(iocean .eq. 1) then
!
         smal_e = 0.0
         smal_e = 1.0e-12
!        smal_e = 6.0e-03
      else
         smal_e = 1.0e-09
!        smal_e = 0.0
      endif
!
! ---------------------- set constants in eddy viscosity model
!
      ck       = 0.1
      ceps     = 0.93
      csmag    = sqrt(ck*sqrt(ck/ceps))
      stab_c   = 0.76
!
! ----------------- set stability constant
!
      stabmin = 1.0e-12
!
! ---------------- minimum dsl length constant
! 
      almin_c = 0.0001
!
! -------------- initialize grid restart flag
!
      igrdr = 1
!
! -------------- create mpi operation to find max and location
!                using local gradient method
!
      call mpi_op_create(get_zi,.true.,ziloc,ierror)
!
! ------------------- define coefficients for 3-order runge-kutta
!                     time stepping scheme, borrowed from Spalart,
!                     Moser and Rogers, J. Comp. Physics 3/21/90
!                     Note this is a simplier version since all terms
!                     are lumped in the non-linear terms.
!                     cfl number is for an entire runge-kutta step
!                     in this case three stages. cfl = max(u)*dt/dx
!
      zetas(1) = 0.0
      zetas(2) = -17.0/60.0
      zetas(3) = -5.0/12.0
      gama(1)  = 8.0/15.0
      gama(2)  = 5.0/12.0
      gama(3)  = 3.0/4.0
      etas(1)  = -1.0
      etas(2)  = -1.0 + 8.0/15.0
      etas(3)  = -1.0 + 2.0/3.0
!
! ----------- a full step, at the new time
!
      etas(4)  =  0.0
!
!
      return
!
      end

! --- extracted from les.F: setup ---
      subroutine setup
      implicit real(a-h,o-z), integer(i-n)
!
!
!
! ------------ turn on new sgs model at a particular step
!
      if(it .ge. new_vis .and. ivis0 .eq. 1) then
          ivis = 1
      else
          ivis = 0
      endif
!
      if(igrdr .eq. 3) then
         if(l_root) then
            write(6,6)iti,utau,tsfcc(1) ,qstar(1)
            write(6,510)
            write(6,520) wtsfc(1),utau,amonin,dtdzf(1),zody,zo &
      ,cdbtm,ugcont
            call print(6,1,nnz)
         endif
          if(l_debug) then
            write(nprt,6)iti,utau,tsfcc(1) ,qstar(1)
            write(nprt,510)
            write(nprt,520) wtsfc(1),utau,amonin,dtdzf(1),zody,zo &
      ,cdbtm,ugcont
            call print(nprt,izs,ize)
          endif
      endif
!     if(ifilt.eq.1)call filter
      if(l_root) then
         write(6,1) nnx,nny,nnz,ismlt,ifilt,iti,itmax, &
      iupwnd,ibuoy,itcut, &
      dt,zo,tsfcc(1),isubs,ibrcl, &
      method, iocean, ivis
      endif
      if(l_debug) then
         write(nprt,1) nnx,nny,nnz,ismlt,ifilt,iti,itmax, &
      iupwnd,ibuoy,itcut, &
      dt,zo,tsfcc(1),isubs,ibrcl, &
      method, iocean, ivis
      endif
!
! -------------- boundary condition flags 
!
      ibcu = iradup
!     ibcu = 0
      ibcl = 0
!
! -------------------- wavenumbers, introduce a normalized
!                      set of wavenumbers to eliminate computation
!                      in derivatives , xderiv, yderiv
!
      do i=1,nnx
         xkn(i) = float(i-1)*pi2/xl
         if(i.gt.ncx)xkn(i) = -float(nnx-i+1)*pi2/xl
      enddo
      fn = 1.0/float(nnx)
      do i=1,nnx
         xk(i) = xkn(i)*fn
      enddo
      do i=1,nny
         ykn(i) = float(i-1)*pi2/yl
         if(i.gt.ncy)ykn(i) = -float(nny-i+1)*pi2/yl
      enddo
      fn = 1.0/float(nny)
      do i=1,nny
         yk(i) = ykn(i)*fn
      enddo
      ii = -1
      do i=1,ncx
         ii = ii + 2
         temp = xkn(i)**2
         do j=1,nny
            temp1       = temp + ykn(j)**2
            xks(ii,j)   = temp1
            xks(ii+1,j) = temp1
         enddo
      enddo
      xnn = abs(batag*dtdzf(1))
!
! ----------- choose correct sign so gravity waves
!             propagate out of the domain
!
      sgn = 1.0
      if(iocean .eq. 1) sgn = -1.0
      if(ibcu.eq.1) then
         do iy=1,nny
         do ix=1,nnxp2
            if(xks(ix,iy) .le. 0.) then
              wavexy(ix,iy) = 0.0
            else
              wavexy(ix,iy) = sgn*sqrt(xnn/xks(ix,iy))
            endif
         enddo
         enddo
      endif
!
! -------------------- set length scale for SGS model
!
      if(iz_space .eq. 0) then
!
! ------------- uniform vertical spacing
!
      dx32 = dx*3./2.
      dy32 = dy*3./2.
      dsl  = (abs(dx32*dy32*dzw(1)))**(1./3.)
      dslg = dsl
      if(l_root)  write(6,2000) dsl
      if(l_debug) write(nprt,2000) dsl
!
! --------------------- create dsl array for easy indexing in comp1
!
      do iz=0,nnzp1
         dsl_z(iz) = dslg
      enddo
!
! ------------- variable vertical spacing
!
      else
!
! ----------- just estimate dsl for average spacing
!
         dx32 = dx*3./2.
         dy32 = dy*3./2.
!
         dsl_max = (abs(dx32*dy32*dzw(0)))**(1./3.)
         do iz=0,nnzp1
            dsl_z(iz) = (abs(dx32*dy32*dzw(iz)))**(1./3.)
            if(dsl_z(iz) .gt. dsl_max) dsl_max = dsl_z(iz)
         enddo
!        do iz=0,nnzp1
!           dsl_z(iz) = dsl_max
!        enddo
         dsl  = dsl_max
         dslg = dsl
      endif
!
      gridr = 1.0
      sml_eg = smal_e*gridr
! -------------------- set viscosity model parameters 
      if(ivis .ne. 1) then
        viscon = 0.0
        xksurf = 0.0
        nmatch = -1
        myid_newvis = 0
        do iz=1,nnz
           dfac(iz) = 1.0
        enddo
      endif
! ------------------- set stokes velocity for atmos/oceanic flow
      call stokesv
!
! --------- can add a time factor so as to skip into any part of
!           the specified geostrophic arrays. time factor in seconds
!
      t_factor = 7200.0
!
! ---------- for print out to get more digits
!
      t_ref = 300.0
!
! -------------------- specify cooling rate and initial
!                      temperature even for restarts
!
      c_rate   = 0.25/3600.0
      t_surf_i = 265.0
!
! -------------------- do not look for zi below zi_min
!
      zi_min = 30.0
      if(iocean .eq. 1) zi_min = -5.0
      iz_min = 1
      do iz=1,nnz-1
         if(zz(iz) .lt. zi_min .and. &
      zz(iz+1) .ge. zi_min) iz_min = iz
      enddo
      if(l_root) then
         write(6,9000) zi_min, iz_min
      endif
!
 9998 continue
      return
! --------------------------- format statements
    6 format(///,' DATA FROM RESTART FILE AT STEP =',I5, &
      ' U_* = ',e15.6,' TS = ',e15.6,' Q_* = ',e15.6,///)
  510 format(' RESTART ***** CASE WITH : ******',/)
  520 format(' WT = ',e12.4,', U_* = ',e12.4,', L = ',e12.4, &
      ', DTDZ FREE = ',e12.4,', ZODY = ',e12.4,/,10x, &
      '  ZO(BTM) = ',e12.4,', CDBTM = ',e12.4, &
      ', UG = ',e12.4)
    1 format(10x,' NNX = ',i5,',  NNY = ',i5, &
      ',  NNZ = ',i5,/,10x,' SFC SMLT = ',i1, &
      ',  FILTER = ',i1, &
      ',  ITI = ',i6,',  ITMAX = ',i8,/,10x, &
      ' IUPWIND = ',i1,',  BUYNCY = ',i1, &
      ',  ITCUT = ',i1,/,10x, &
      ' DT = ',e15.6,',  ZO = ',e15.6,',  TS = ',e15.6, &
      ',  SUBSD = ',i1,/, &
      10x,' BRCLICITY = ',i1,',  METHOD = ',i1,',  IOCEAN = ',i1, &
      ',  IVIS = ',i1)
 2000 format(10x,' DSL = ',e15.6)
 9000 format(' Search for zi above the height = ',e15.6,/, &
      ' iz_min = ',i5)
      end

! --- extracted from les.F: init ---
      subroutine init
      implicit real(a-h,o-z), integer(i-n)
!
!

      pi   = 4.0*atan(1.0)
      pi2  = 2.0*pi
      bfac = 1.0
      if(ibuoy.eq.0) bfac = 0.
!
! -------------------- case specific data
!
      if(iocean .eq. 1) then
         t00     = 283.
         t00b    = 5000.0
         cp      = 4.20e03
         gcp     = grav/cp
         batag   = bfac*grav/t00b
!        fcor    = 0.0
         fcor    = 1.39e-04
         fcor_h  = 0.0
!        wtsfc(1)=0.00
!        wtsfc(1)=4.96e-07
         wtsfc(1)=1.190476e-06
         qstar(1)=wtsfc(1)
!        dtdzf(1)=0.000
         dtdzf(1)=0.2548
         dtjump  = 0.
         divgls  = 0.
         zo      = 0.0001
         zi      = -5.
         izi     = 55
         xl      = 50.
         yl      = 50.
         zl      = -20.
!
! ---------- if stretched grid specify location of first point
!
         zw1 = -0.5
      else
         gcp     = grav/Cpa
         batag   = bfac*grav/t00
         fcor_h  = 0.0

         wtsfc = qstar
!
!
         !Gradients of temperature and scalars above inversion (upper BC)
         if (icase.eq.5) then  !Stratocumulus

         dtdzf(1)=(311.85-308.2)/(3000.0-2000.0)
         dtdzf(2)=(3.0e-3-4.2e-3)/(3000.0-2000.0)

         elseif (icase.eq.3) then !C-FOG

         dtdzf(1)=(285.5-284.0)/(80.0-30.0)
         dtdzf(2)=(0.006-0.00805)/(80.0-30.0)         

         elseif (icase.eq.6) then !FATIMA

         dtdzf(1)=(289.5-287.51)/(128.0-70.0)
         dtdzf(2)=(0.0094-0.0098)/(128.0-70.0)

         end if

         dtjump  = 0.0
         divgls  = 0.0

      endif
!
      time  = 0.0
! 
! ---------- outermost coarse grid  indicies are bounds of grid
!
      izlow = 1
      izup  = nnz
      dz    = zl/nnz
      dzg   = abs(dz)
      if(l_root) write(6,4040) zl,nnz,dzg
!
! --------------- generate z grids for particular mesh from
!                 iz = 0,1,...,nnz+1; this allows indexing
!                 to array elements z(0), etc.
!
      zwstrt = 0.0

!
! ------------ build z grid for w points
!
      if(iz_space .eq. 0) then
         do iz=0,nnz+1
            z(iz) = dz*float(iz) + zwstrt
         enddo
      elseif (iz_space .eq. 1) then
        call vgrid_channel(zw1,zi,zl,nnz,z(0:),l_root,l_debug)
      elseif (iz_space .eq. 2) then
        call vgrid_channel_fstrm(zw1,zi,zl,nnz,z(0:),l_root,l_debug)
      elseif (iz_space .eq. 3) then
        call vgrid(zw1,zi,zl,nnz,z(0:),l_root,l_debug)
      endif
!
      call get_dz
!
      if(l_root) then
         write(6,8002) zwstrt
         write(6,8003) (iz,z(iz),zz(iz),iz=0,nnz+1)
      endif
!
      nnzm1 = nnz-1
      dx    = xl/nnx
      dy    = yl/nny
      fnxy  = 1./float(nxy)
      dzdz  = dzw(1)*dzw(1)
      z1    = zz(1)
!
      c23  = 2.0/3.0
      dsl  = (dx*1.5*dy*1.5*abs(dzw(1)))**(1./3.)
      dslg = dsl
      cs   = 0.2
!
      vk     = 0.4
      batagk = batag*vk
      vkin   = 1./vk
      ttmean = 0.
      zody   = alog(abs(z1/zo))
      zosdy   = alog(abs(z1/zos))
      write(nprt, 9901) z1,zo,zody,zosdy
 9901 format(' 9901 z1 = ',e15.6,' zo = ',e15.6,/, &
      ' zody = ',e15.6,' zosdy = ',e15.6)
      zodyin = 1./zody
      wstar  = abs(batag*zi*wtsfc(1))**(1./3.)
      if(ismlt .eq. 1) then
!
! ---- set constants for businger similarity functions
!
         vk74   = vk*0.74
         vk74in = 0.74/vk
         zody74 = zody*0.74
      else 
!
! ---- set constants for large similarity functions
!
        vk74    = vk
        vk74in  = 1.0/vk
        zody74  = zody
      endif
      ugal   = 0.0
!      ugal   = ugcont*0.5
!     ugcont = ugcont - ugal
      cdbtm  = vk*vk/zody/zody
      if(iocean .eq. 1) then
! ----------- set surface friction velocity here and in sr. sufto
!        utau = 4.29e-03
         utau = 7.00e-03
      else
         ufree = 0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)
!
! ---- note : new estimate for utau !!!
!
         utau  = vk*(ufree+ugcont)/zody
!        utau  = vk*(ufree)/zody
      endif
      utau2    = utau*utau
      if(ibuoy .eq. 0 .or. qstar(1) .eq. 0.) then
        amonin = 1000.0
      else
        amonin = -utau2*utau/(batagk*qstar(1))
      endif
      hol   = abs(zi)/amonin
      zol   = abs(z1)/amonin
      uwsfc = -utau*utau
      vwsfc = -utau*utau
!
      if(l_root) then
         write(6,80)
         write(6,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo &
      ,cdbtm,ugcont
      endif
!
      if(l_debug) then
         write(nprt,80)
         write(nprt,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo &
      ,cdbtm,ugcont
      endif
!
      return
! ------------------------
   2  format(10x,' WT =',e12.4,',  U* =',e12.4,',  L =',e12.4,/, &
      10x,' DTDZ FREE =',e12.4,',  ZODY=',e12.4,/,10x, &
      ' ZO(BTM) =',e12.4,',  CDBTM=',e12.4, &
      ',  UG = ',e12.4)
  80  format(///,' ***** SCRATCH RUN ***** ',//)
 4040 format(' zl = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 4043 format(' znest = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 8002 format(' zwstrt = ',e12.4)
 8003 format(' iz ',5x,' zw',5x,' zu ',5x,/,(i3,2e12.4))
      end

! --- extracted from les.F: vgrid ---
      subroutine vgrid(z1,zi,zl,nnz,z,l_root,ldebug)
      implicit real(a-h,o-z), integer(i-n)
!
      real z(0:nnz+1)
      logical l_root, ldebug
!
! ----------------- build grid up to zi first
!
      z_frst = z1
      z_cntr = zi*0.5
      !n_pbl  = nnz/2
      n_pbl  = (3*nnz)/4
      z_fac1 = z_cntr/z_frst
      z_fac2 = 1.0/float(n_pbl/2)
      z_fac  = 1.1
      knt = 0
      tol = 0.00001
   10 continue
        knt = knt + 1
        z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
        test   = abs(1.0 - z_facn/z_fac)
        if(knt .gt. 50) then
            if(l_root) write(6,9000) z_fac, z_facn, knt
 9000       format(' Cannot find stretching factor',/, &
      ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 10
      if(l_root) write(6,9100) z_fac, z_cntr, z1, knt
 9100 format(' Stretching factor = ',e15.6,/, &
      ' Match point       = ',e15.6,/, &
      ' First z           = ',e15.6,/, &
      ' Number of iters   = ',i4)
      z(1) = z_frst
      do iz=2,n_pbl/2-1
         z(iz) = z_frst*(z_fac**(float(iz)) - 1.0)/(z_fac - 1.0)
      enddo
      z(n_pbl/2) = z_cntr
      do iz=1,n_pbl/2 - 1
         z(n_pbl-iz) = zi - z(iz)
      enddo
      z(n_pbl) = zi
      z(0)   = 0.0
!
      if(l_root) write(6,5300) n_pbl
 5300 format(' n_pbl = ',i4)
!
! -------------- build grid from zi on up
!
      z_frst = z1
      z_cntr = zl - zi
      n_top  = nnz - n_pbl
      z_fac1 = z_cntr/z_frst
      z_fac2 = 1.0/float(n_top)
      z_fac  = 1.1
      knt = 0
      tol = 0.00001
   20 continue
        knt = knt + 1
        z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
        test   = abs(1.0 - z_facn/z_fac)
        if(knt .gt. 50) then
            if(l_root) write(6,8000) z_fac, z_facn, knt
 8000       format(' Cannot find stretching factor',/, &
      ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 20
      if(l_root) write(6,8100) z_fac, z_cntr, z1, knt
 8100 format(' Stretching factor = ',e15.6,/, &
      ' Match point       = ',e15.6,/, &
      ' First z           = ',e15.6,/, &
      ' Number of iters   = ',i4)
!
      z(n_pbl+1) = zi + z_frst
      do iz=n_pbl+2,nnz-1
         z(iz) = zi + z_frst* &
      (z_fac**(float(iz-n_pbl)) - 1.0)/(z_fac - 1.0)
      enddo
      z(nnz) = zl
      z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))
!     if(l_root) write(6,5600) (iz,z(iz),iz=0,nnz+1)
 5600 format(' 5600 in vgrid ',/, &
      ' iz ',5x,' zw ',/,(i3,e15.6))
!
!     write(1,2000)
!2000 format('#k ',/,
!    +       '#lw 0.5 ',/,
!    +       '#m 1',/,
!    +       '#x 0 100 50',/,
!    +       '#y -50 2100 500')
!     x1 = 30.0
!     x2 = 80.0
!     do iz=0,nnz+1
!        write(1,1000) x1,z(iz)
!1000    format('#k ',/,
!    +          (2e15.6))
!        write(1,1100) x2,z(iz)
!1100    format(2e15.6)
!     enddo
!
      return
      end

! --- extracted from les.F: vgrid_channel ---
      subroutine vgrid_channel(z1,zi,zl,nnz,z,l_root,ldebug)
      implicit real(a-h,o-z), integer(i-n)
!
      real z(0:nnz+1)
      integer :: zidx
      logical l_root, ldebug
!
! ----------------- build grid up to zi first
!
      z_frst = z1
      z_cntr = zi*0.5
      n_pbl  = nnz
!     n_pbl  = (5*nnz)/8
      z_fac1 = z_cntr/z_frst
      z_fac2 = 1.0/float(n_pbl/2)
      z_fac  = 1.1
      knt = 0
      tol = 0.00001
   10 continue
        knt = knt + 1
        z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
        test   = abs(1.0 - z_facn/z_fac)
        if(knt .gt. 50) then
            if(l_root) write(6,9000) z_fac, z_facn, knt
 9000       format(' Cannot find stretching factor',/, &
      ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 10
      if(l_root) write(6,9100) z_fac, z_cntr, z1, knt
 9100 format(' Stretching factor = ',e15.6,/, &
      ' Match point       = ',e15.6,/, &
      ' First z           = ',e15.6,/, &
      ' Number of iters   = ',i4)
      z(1) = z_frst
      do iz=2,n_pbl/2-1
         z(iz) = z_frst*(z_fac**(float(iz)) - 1.0)/(z_fac - 1.0)
      enddo
      z(n_pbl/2) = z_cntr
      do iz=1,n_pbl/2 - 1
         z(n_pbl-iz) = zi - z(iz)
      enddo
      z(n_pbl) = zi
      z(0)   = 0.0
!
      if(l_root) write(6,5300) n_pbl
 5300 format(' n_pbl = ',i4)
!
! -------------- build grid from zi on up
!     For the channel, zi represents the channel centerline
!     Want the mesh to be a mirror image across this:
!
!      zidx = 1
!      do iz=n_pbl+1,nnz
!         z(iz) = zi + (zi - z(n_pbl-zidx))
!         zidx = zidx + 1
!      enddo
      z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))
!
      return
      end subroutine vgrid_channel

! --- extracted from les.F: vgrid_channel_fstrm ---
      subroutine vgrid_channel_fstrm(z1,zi,zl,nnz,z,l_root,ldebug)
      implicit real(a-h,o-z), integer(i-n)
!
      real z(0:nnz+1)
      real s(0:2*nnz+1)
      integer :: zidx
      logical l_root, ldebug
!      
                nnz = 2.0*nnz
                zi = 2.0*zi
                zl = 2.0*zl
!
! ----------------- build grid up to zi first
!
      z_frst = z1
      z_cntr = zi*0.5
      n_pbl  = nnz
!     n_pbl  = (5*nnz)/8
      z_fac1 = z_cntr/z_frst
      z_fac2 = 1.0/float(n_pbl/2)
      z_fac  = 1.1
      knt = 0
      tol = 0.00001
   10 continue
        knt = knt + 1
        z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
        test   = abs(1.0 - z_facn/z_fac)
        if(knt .gt. 50) then
            if(l_root) write(6,9000) z_fac, z_facn, knt
 9000       format(' Cannot find stretching factor',/, &
      ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 10
      if(l_root) write(6,9100) z_fac, z_cntr, z1, knt
 9100 format(' Stretching factor = ',e15.6,/, &
      ' Match point       = ',e15.6,/, &
      ' First z           = ',e15.6,/, &
      ' Number of iters   = ',i4)
      s(1) = z_frst
      do iz=2,n_pbl/2-1
         s(iz) = z_frst*(z_fac**(float(iz)) - 1.0)/(z_fac - 1.0)
      enddo
      s(n_pbl/2) = z_cntr
      do iz=1,n_pbl/2 - 1
         s(n_pbl-iz) = zi - s(iz)
      enddo
      s(n_pbl) = zi
      s(0)   = 0.0
!
      if(l_root) write(6,5300) n_pbl
 5300 format(' n_pbl = ',i4)
!
      s(nnz+1) = s(nnz) + (s(nnz) - s(nnz-1))

                nnz=nnz/2.0
                zi = zi/2.0
                zl = zl/2.0

                do iz=0,nnz
                        z(iz)=s(iz)
                enddo
                 z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))

      return
      end subroutine vgrid_channel_fstrm

! --- extracted from les.F: vgrid_uniform ---
      subroutine vgrid_uniform(z1,zi,zl,nnz,z,l_root,ldebug)
      implicit real(a-h,o-z), integer(i-n)
!
      real z(0:nnz+1),zdiff
      real s(0:2*nnz+1)
      integer :: zidx
      logical l_root, ldebug

      z_frst = z1
      n_pbl  = nnz
      z_fac  = 1
      zdiff = (zl-z_frst)/(nnz-1)
      s(1) = z_frst
      do iz=2,nnz
         s(iz) = z_frst+(zdiff*((float(iz)-1)))
      enddo
      s(0)   = 0.0
!
      s(nnz+1) = s(nnz) + (s(nnz) - s(nnz-1))

      do iz=0,nnz
         z(iz)=s(iz)
      enddo
         z(nnz+1) = z(nnz) + zdiff
      return
      end subroutine vgrid_uniform

! --- extracted from les.F: get_dz ---
      subroutine get_dz
      implicit real(a-h,o-z), integer(i-n)
!
! --------------- compute spacing for given vertical
!                 point distribution
!
!
      do iz=1,nnz+1
         dzw(iz) = z(iz) - z(iz-1)
      enddo
      dzw(0)     = dzw(1)
      dzw(nnz+2) = dzw(nnz+1)
      do iz=0,nnz+2
         dzw_i(iz) = 1.0/dzw(iz)
      enddo
!
! ------------ build z grid for u points
!
      dzovr2 = dz*0.5
      do iz=1,nnz+1
         zz(iz) = 0.5*(z(iz) + z(iz-1))
      enddo
      zz(0) = - zz(1)
      do iz=1,nnz+1
         dzu(iz) = zz(iz) - zz(iz-1)
      enddo
      dzu(0)     = dzu(1)
      dzu(nnz+2) = dzu(nnz+1)
      do iz=0,nnz+2
         dzu_i(iz) = 1.0/dzu(iz)
      enddo
!
      return
      end

! --- extracted from les.F: random ---
      subroutine random
!
! ----------- geostrophic winds designed for comparison case
!
      implicit none

      real :: psi(nnx, iys:iye), psix(nnx, iys:iye), &
      psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye), &
      vyy(nnx,iys:iye,izs:izs),Ttmp

      real :: ampv, ampt, sum_psi, vmaxx, vmag, facv

      real :: slope_q, slope_t, slope_u, slope_v
      real :: slope_t1, slope_t2, slope_u1, slope_v1
      real :: z_switch, t_switch, q_switch, u_switch, v_switch
      real :: z_switch_q, z_switch_t, z_switch_v, z_switch_u
      real :: t_surf, q_surf, u_surf, v_surf
      real :: t_temp, q_temp, u_temp, v_temp
      real :: q_top, t_top, u_top, v_top, slope_q1, slope_q2

      integer :: ix,iy,iz
      integer :: idum


! ------------ note set nmatch in sr. iso so that
!              it is compatible with conditions here
!
      do iz=1,nnz
         ug(iz)   = ugcont
         vg(iz)   = vgcont
         divz(iz) = 0.0
      enddo

      if (icase.eq.5) then

      do iz=izs,ize

          ! MAGPIE
!         profile for potential temperature
            if (zz(iz) .le. 600.0) then

                  t_temp = 302.75   !301.75   !300.0   !301.75(special)   !302.5(ori)

            elseif (zz(iz) .gt. 600.0 .and. zz(iz) .le. 1900.0) then

                  t_surf   = 302.75   !301.75   !300.0   !301.75(special)   !302.5(ori)
                  t_switch = 308.0   !316.0   !315.5(special) !8   !313.0 !5   !308.0(ori)
                  slope_t  = (t_switch - t_surf) / (1900.0 - 600.0)
                  t_temp   = t_surf + (slope_t * (zz(iz) - 600.0))

            elseif (zz(iz) .gt. 1900.0 .and. zz(iz) .le. 2400.0) then

                  t_surf   = 308.0   !316.0   !315.5(special) !8   !313.0 !5   !308.0(ori)
                  t_switch = 315.0   !324.0   !323.0(special) !9   !321.0 !6   !315.0(ori)
                  slope_t  = (t_switch - t_surf) / (2400.0 - 1900.0)
                  t_temp   = t_surf + (slope_t * (zz(iz) - 1900.0))

            elseif (zz(iz) .gt. 2400.0 .and. zz(iz) .le. 3000.0) then

                  t_surf   = 315.0   !324.0   !323.0(special) !9   !321.0 !6   !315.0(ori)
                  t_switch = 317.25   !327.25   !326.25(special) !10   !324.25 !7  !317.25(ori)
                  slope_t  = (t_switch - t_surf) / (3000.0 - 2400.0)
                  t_temp   = t_surf + (slope_t * (zz(iz) - 2400.0))

            elseif (zz(iz) .ge. 3000.0) then

                  t_temp = 317.25   !327.25   !326.25(special) !10   !324.25 !7   !317.25(ori)

            endif

!         profile for qv
            if (zz(iz) .le. 600.0) then

                  q_surf   = 19.25e-3   !19.5e-3
                  q_switch = 18.0e-3   !17.05e-3
                  slope_q  = (q_switch - q_surf) / 600.0
                  q_temp   = q_surf + (slope_q * zz(iz))

            elseif (zz(iz) .gt. 600.0 .and. zz(iz) .le. 1800.0) then

                  q_surf   = 18.0e-3   !17.05e-3
                  q_switch = 11.5e-3   !8.5e-3
                  slope_q  = (q_switch - q_surf) / (1800.0 - 600.0)
                  q_temp   = q_surf + (slope_q * (zz(iz) - 600.0))

!            elseif (zz(iz) .gt. 1000.0 .and. zz(iz) .le. 2000.0) then

!                  q_surf   = 14.0e-3
!                  q_switch = 9.0e-3
!                  slope_q  = (q_switch - q_surf) / (2000.0 - 1000.0)
!                  q_temp   = q_surf + (slope_q * (zz(iz) - 1000.0))

            elseif (zz(iz) .gt. 1800.0 .and. zz(iz) .le. 2400.0) then

                  q_surf   = 11.5e-3   !8.5e-3
                  q_switch = 5.15e-3
                  slope_q  = (q_switch - q_surf) / (2400.0 - 1800.0)
                  q_temp   = q_surf + (slope_q * (zz(iz) - 1800.0))

            elseif (zz(iz) .gt. 2400.0 .and. zz(iz) .le. 3000.0) then

                  q_surf   = 5.15e-3
                  q_switch = 4.25e-3   !3.0e-3
                  slope_q  = (q_switch - q_surf) / (3000.0 - 2400.0)
                  q_temp   = q_surf + (slope_q * (zz(iz) - 2400.0))

            elseif (zz(iz) .ge. 3000.0) then

                  q_temp = 4.25e-3   !3.0e-3

            endif

!         profile for u
            if (zz(iz) .le. 3000.0) then

                  u_surf   = -5.0   !-4.0
                  u_switch = -15.0   !-7.25
                  slope_u  = (u_switch - u_surf) / 3000.0
                  u_temp   = u_surf + (slope_u * zz(iz))

!            elseif (zz(iz) .gt. 100.0 .and. zz(iz) .le. 700.0) then

!                  u_surf = -7.25
!                  u_switch = -8.0
!                  slope_u = (u_switch - u_surf) / (700.0 - 100.0)
!                  u_temp = u_surf + (slope_u * (zz(iz) - 100.0))

!            elseif (zz(iz) .gt. 700.0 .and. zz(iz) .le. 3000.0) then

!                  u_surf = -8.0
!                  u_switch = -15.5
!                  slope_u = (u_switch - u_surf) / (3000.0 - 700.0)
!                  u_temp = u_surf + (slope_u * (zz(iz) - 700.0))

            elseif (zz(iz) .ge. 3000.0) then

                  u_temp = -15.0   !-15.5

            endif

!         profile for v
            if (zz(iz) .le. 3000.0) then

                  v_surf   = -3.5   !-3.25
                  v_switch = 0.0   !-6.5
                  slope_v  = (v_switch - v_surf) / 3000.0
                  v_temp   = v_surf + (slope_v * zz(iz))

!            elseif (zz(iz) .gt. 250.0 .and. zz(iz) .le. 3000.0) then

!                  u_surf = -6.5
!                  u_switch = 0.0
!                  slope_u = (u_switch - u_surf) / (3000.0 - 250.0)
!                  u_temp = u_surf + (slope_u * (zz(iz) - 250.0))

            elseif (zz(iz) .ge. 3000.0) then

                  u_temp = 0.0

            endif



      !!! BOMEX initial profile
!      do iz=1,nnz
!            ug(iz)   = -10.0+1.8e-3*zz(iz)
!            vg(iz)   = 0.0
!            divz(iz) = 0.0
!      enddo       
!
!      do iz=izs,ize
!c         profile for pot temp
!            if (zz(iz).le.520.0) then
!                  t_temp = 298.7
!            elseif (zz(iz).gt.520.0 .and. zz(iz).le.1480.0) then
!                  t_surf = 298.7
!                  t_switch = 302.4
!                  slope_t = (t_switch-t_surf)/(1480.0-520.0)
!                  t_temp = t_surf + slope_t*(zz(iz)-520.0)
!            elseif (zz(iz).gt.1480.0 .and. zz(iz).le.2000.0) then
!                  t_surf = 302.4
!                  t_switch = 308.2
!                  slope_t = (t_switch-t_surf)/(2000.0-1480.0)
!                  t_temp = t_surf + slope_t*(zz(iz)-1480.0)
!            elseif (zz(iz).gt.2000.0 .and. zz(iz).le.3000.0) then
!                  t_surf = 308.2
!                  t_switch = 311.85
!                  slope_t = (t_switch-t_surf)/(3000.0-2000.0)
!                  t_temp = t_surf + slope_t*(zz(iz)-2000.0)
!            elseif (zz(iz).ge.3000.0) then
!                  t_temp = 311.85
!            endif
!
!c         profile for qv
!             if (zz(iz).le.520.0) then
!                  q_surf = 17.0e-3
!                  q_switch = 16.3e-3
!                  slope_q = (q_switch-q_surf)/520.0
!                  q_temp = q_surf + slope_q*zz(iz)
!            elseif (zz(iz).gt.520.0 .and. zz(iz).le.1480.0) then
!                  q_surf = 16.3e-3
!                  q_switch = 10.7e-3
!                  slope_q = (q_switch-q_surf)/(1480.0-520.0)
!                  q_temp = q_surf + slope_q*(zz(iz)-520.0)
!            elseif (zz(iz).gt.1480.0 .and. zz(iz).le.2000.0) then
!                  q_surf = 10.7e-3
!                  q_switch = 4.2e-3
!                  slope_q = (q_switch-q_surf)/(2000.0-1480.0)
!                  q_temp = q_surf + slope_q*(zz(iz)-1480.0)
!            elseif (zz(iz).gt.2000.0 .and. zz(iz).le.3000.0) then
!                  q_surf = 4.2e-3
!                  q_switch = 3.0e-3
!                  slope_q = (q_switch-q_surf)/(3000.0-2000.0)
!                  q_temp = q_surf + slope_q*(zz(iz)-2000.0)
!            elseif (zz(iz).ge.3000.0) then
!                  q_temp = 3.0e-3
!            endif
!
!c         profile for u
!             if (zz(iz).le.700.0) then
!                  u_temp = -8.75
!            elseif (zz(iz).gt.700.0 .and. zz(iz).le.3000.0) then
!                  u_surf = -8.75
!                  u_switch = -4.61
!                  slope_u = (u_switch-u_surf)/(3000.0-700.0)
!                  u_temp = u_surf + slope_u*(zz(iz)-700.0)
!            elseif (zz(iz).ge.3000.0) then
!                  u_temp = -4.61
!            endif

            do iy=iys,iye
            do ix=1,nnx
                  u(ix,iy,iz) = u_temp
                  v(ix,iy,iz) = 0.0
                  w(ix,iy,iz) = 0.0
                  e(ix,iy,iz) = 1.0-zz(iz)/3000.0

!                 pot temp
                  t(ix,iy,1,iz) = t_temp

!                 qv
                  t(ix,iy,2,iz) = q_temp

                  w(ix,iy,iz)   = 0.
                  r1(ix,iy,iz)  = 0.
                  r2(ix,iy,iz)  = 0.
                  r3(ix,iy,iz)  = 0.
                  r4(ix,iy,1,iz)= 0.
                  r4(ix,iy,2,iz)= 0.
                  r5(ix,iy,iz)  = 0.
             enddo
             enddo
       enddo

      elseif (icase.eq.6) then 
      !FATIMA fog setup

            !qv setup
            z_switch_q = 70.0
            q_surf = 0.0099
            q_switch = 0.0098
            q_top = 0.0094
            slope_q1 = (q_switch-q_surf)/z_switch_q
            slope_q2 = (q_top-q_switch)/(128.0-z_switch_q)

            !theta setup
            z_switch_t = 70.0
            t_surf = 287.5
            t_switch = 287.51
            t_top = 289.5
            slope_t1 = (t_switch-t_surf)/z_switch_t
            slope_t2 = (t_top-t_switch)/(128.0-z_switch_t)


            !u setup
            z_switch_u = 128.0
            u_surf = 3.8
            u_switch = 4.1
            u_top = 4.1
            slope_u1 = (u_top-u_surf)/128.0

            !v setup
            z_switch_v = 128.0
            v_surf = 1.4
            v_switch = 0.4
            v_top = 0.4
            slope_v1 = (v_top-v_surf)/128.0

            do iz=izs,ize
            do iy=iys,iye
            do ix=1,nnx

                ! u
                !u(ix,iy,iz) = u_surf + slope_u1*zz(iz)
                u(ix,iy,iz) = 4.0
                 ! v
                !v(ix,iy,iz) = v_surf + slope_v1*zz(iz)
                v(ix,iy,iz) = 0.0

                ! q
                if (zz(iz) .le. z_switch_q) then
                    t(ix,iy,2,iz) = q_surf + slope_q1*zz(iz)
                elseif (zz(iz) .gt. z_switch_q) then
                    t(ix,iy,2,iz) = q_switch + &
      slope_q2*(zz(iz)-z_switch_q)
                end if

                ! theta
                if (zz(iz) .le. z_switch_t) then
                    t(ix,iy,1,iz) = t_surf + slope_t1*zz(iz)
                elseif (zz(iz) .ge. z_switch_t) then
                    t(ix,iy,1,iz) = t_switch + &
      slope_t2*(zz(iz)-z_switch_t)
                end if

                w(ix,iy,iz) = 0.0
                e(ix,iy,iz) = 0.0
            enddo
            enddo
            enddo

      else  !switch the case


      if (icase.eq.3) then
      !!!! Use for the C-FOG fog LCM case:
      !Parameters from Charlotte's C-FOG profile
      z_switch = 30.0
      t_switch = 284.05
      q_switch = 0.00811
      t_surf = 284.0
      q_surf = 0.00823
      end if


      if (icase.eq.4) then
      z_switch = 600.0
      t_switch = 300.81
      t_surf = 300.81
      q_switch = 0.0215
      q_surf = 0.0215
      end if

      if (icase.eq.0) then
      t_surf = 285.0
      t_switch = 285.0
      q_surf = 0.0105
      q_switch = 0.0105
      slope_t = 0.0
      slope_q = 0.0
      z_switch = 1000.0
      end if

      !Slopes below the inversion
      slope_t = (t_switch-t_surf)/z_switch
      slope_q = (q_switch-q_surf)/z_switch

!
      do iz=izs,ize


         !Try and alleviate the abrupt IC on the velocity
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = ugcont
            v(ix,iy,iz) = vgcont
            w(ix,iy,iz) = 0.0
            e(ix,iy,iz) = 0.0
         enddo
         enddo

         if(zz(iz) .le. z_switch) then
           do iy=iys,iye
           do ix=1,nnx
              t(ix,iy,1,iz) = t_surf + slope_t*zz(iz)
              t(ix,iy,2,iz) = q_surf + slope_q*zz(iz)
           enddo
           enddo
         elseif(zz(iz) .ge. z_switch) then
           do iy=iys,iye
           do ix=1,nnx
              t(ix,iy,1,iz) = t_switch + (zz(iz) - z_switch)*dtdzf(1)
              t(ix,iy,2,iz) = q_switch + (zz(iz) - z_switch)*dtdzf(2)
           enddo
           enddo
         endif

         do iy=iys,iye
         do ix=1,nnx
            w(ix,iy,iz)   = 0.
            r1(ix,iy,iz)  = 0.
            r2(ix,iy,iz)  = 0.
            r3(ix,iy,iz)  = 0.
            r4(ix,iy,1,iz)= 0.
            r4(ix,iy,2,iz)= 0.
            r5(ix,iy,iz)  = 0.
         enddo 
         enddo 
      enddo

      end if
!
! ------------- set initial random field to be
!               divergence free
!
      idum = -1 - myid
      do iz=izs,ize
!
! ----------- ampv and ampt are max amplitudes of random 
!             velocity and temperature fields
!             make sure ampv is set if free convection so
!             that we have motions at first time step
!
         ampv = 0.0
         ampv = 0.001
         ampt = 0.10
!  
! ------- simple random field scaled between -0.5 and 0.5
!
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
         call xderivp(psix(1,iys),trigx(:,1),xk,nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(:,2),yk, &
      nnx,nny,ixs,ixe,ix_s,ix_e, &
      iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         vmaxx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
            if(vmag .gt. vmaxx) vmaxx = vmag
         enddo
         enddo
         facv = ampv/vmaxx
!
         if (z(iz) .le. 50.0) then
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz)   = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
            v(ix,iy,iz)   = v(ix,iy,iz) + psix(ix,iy)*facv
            t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
         enddo
         enddo
         endif
!
         if(z(iz) .le. 250.0) then
         do iy=iys,iye
         do ix=1,nnx
            e(ix,iy,iz) = 0.4*(1.0 - z(iz)/250.0)**3
         enddo
         enddo
         endif
!
! ---------- check divergence of initial field
!
         do iy=iys,iye
         do ix=1,nnx
            uxx(ix,iy) = u(ix,iy,iz)
            vyy(ix,iy,izs) = v(ix,iy,iz)
         enddo
         enddo
         call xderivp(uxx(1,iys),trigx(:,1),xk,nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(:,2),yk, &
      nnx,nny,ixs,ixe,ix_s,ix_e, &
      iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         do iy=iys,iye
         do ix=1,nnx
            divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
         enddo
         enddo
         divz(iz) = divz(iz)*fnxy
!
! -------- end z loop
!
      enddo
!
      call mpi_sum_z(divz,i_root,myid,nnz,1)
!
      write(nprt,6000)
 6000 format(' check of divergence for initial state',/, &
      ' iz ',5x,' divergence')
      write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100 format(i5,e15.6)
!
! ------------ fix for baroclinic and subsidence effects !!
!
!     do iz=izs,ize
!        ug(iz)=ugcont
!        vg(iz)=vgcont
!        if (.not.(ibrcl.eq.1)) go to 19988
!        if (.not.(iz.le.izi)) go to 19987
!        ug(iz)=0.
!        vg(iz)=0.
! 19987    continue
! 19988    continue
!        zz2=zz(iz)
!        wls(iz)=-divgls*zz2
!        if (.not.(iz.eq.1)) go to 19986
!        do ix=1,nnx
!        uls(ix)=divgls*(dx*float(ix-1)-xl*.5)
!        enddo
!     enddo
!     write(nprt,9)(uls(ix),ix=1,nnx)
!  9  format(1x,8e12.3)
! 19986 continue
!
      return
      end

! --- extracted from les.F: random_f ---
      subroutine random_f
      implicit real(a-h,o-z), integer(i-n)
!
! ---------- example of using given (sparse) initial 
!            sounding profiles (FIX for ncpu_s).
!            
      real psi(nnx,iys:iye), psix(nnx,iys:iye), &
      psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye), &
      vyy(nnx,iys:iye,izs:izs)
!
      parameter (nt=12, nz=11)
      real zg(nz), u_i(nz,nt), v_i(nz,nt), theta_i(nz,nt)
      real ui_temp(nz), vi_temp(nz), ti_temp(nz)
      real time_g(nt)
!
      data time_g / &
      0.0000E+00,  0.3600E+04,  0.7200E+04,  0.1080E+05,  0.1440E+05, &
      0.1800E+05,  0.2160E+05,  0.2520E+05,  0.2880E+05,  0.3240E+05, &
      0.3600E+05,  0.3960E+05 &
      /
      data zg / &
      0.1000E+02,  0.3000E+02,  0.5500E+02,  0.9000E+02,  0.1400E+03, &
      0.2150E+03,  0.3300E+03,  0.5000E+03,  0.7500E+03,  0.1100E+04, &
      0.1600E+04 &
      /
      data u_i / &
      -0.1510E+01, -0.1560E+01, -0.1580E+01, -0.1580E+01, -0.1560E+01, &
      -0.1530E+01, -0.1510E+01, -0.9000E+00, -0.1390E+01, -0.1220E+01, &
      -0.5100E+00, &
      -0.1090E+01, -0.1110E+01, -0.1120E+01, -0.1120E+01, -0.1030E+01, &
      -0.9900E+00, -0.9500E+00, -0.6200E+00, -0.1230E+01, -0.9400E+00, &
      0.2800E+00, &
      -0.9100E+00, -0.9200E+00, -0.9100E+00, -0.9000E+00, -0.8800E+00, &
      -0.8400E+00, -0.8000E+00, -0.6500E+00, -0.1510E+01, -0.1070E+01, &
      0.2400E+00, &
      -0.8900E+00, -0.8900E+00, -0.8900E+00, -0.8800E+00, -0.8700E+00, &
      -0.8500E+00, -0.8100E+00, -0.7000E+00, -0.1830E+01, -0.8400E+00, &
      0.3500E+00, &
      -0.1250E+01, -0.1260E+01, -0.1260E+01, -0.1250E+01, -0.1240E+01, &
      -0.1220E+01, -0.1160E+01, -0.8800E+00, -0.1980E+01, -0.1900E+00, &
      0.7500E+00, &
      -0.1800E+01, -0.1810E+01, -0.1820E+01, -0.1820E+01, -0.1800E+01, &
      -0.1780E+01, -0.1710E+01, -0.1150E+01, -0.1960E+01,  0.3900E+00, &
      0.9200E+00, &
      -0.2110E+01, -0.2130E+01, -0.2140E+01, -0.2140E+01, -0.2130E+01, &
      -0.2110E+01, -0.2050E+01, -0.9300E+00, -0.1400E+01,  0.8800E+00, &
      0.9600E+00, &
      -0.2250E+01, -0.2280E+01, -0.2290E+01, -0.2300E+01, -0.2290E+01, &
      -0.2260E+01, -0.2070E+01, -0.4000E-01, -0.1600E+00,  0.1440E+01, &
      0.1190E+01, &
      -0.2160E+01, -0.2200E+01, -0.2220E+01, -0.2220E+01, -0.2220E+01, &
      -0.2190E+01, -0.1610E+01,  0.1470E+01,  0.1420E+01,  0.2050E+01, &
      0.1610E+01, &
      -0.2230E+01, -0.2270E+01, -0.2290E+01, -0.2300E+01, -0.2300E+01, &
      -0.2260E+01, -0.1350E+01,  0.2480E+01,  0.2380E+01,  0.2320E+01, &
      0.1740E+01, &
      -0.1890E+01, -0.1930E+01, -0.1950E+01, -0.1950E+01, -0.1940E+01, &
      -0.1890E+01, -0.1120E+01,  0.3010E+01,  0.3030E+01,  0.2800E+01, &
      0.2000E+01, &
      -0.1210E+01, -0.1230E+01, -0.1240E+01, -0.1230E+01, -0.1210E+01, &
      -0.1140E+01, -0.4600E+00,  0.3320E+01,  0.3510E+01,  0.3420E+01, &
      0.2340E+01 &
      /
      data v_i / &
      0.4800E+00,  0.5100E+00,  0.5300E+00,  0.5700E+00,  0.6900E+00, &
      0.7300E+00,  0.7600E+00,  0.1410E+01, -0.4200E+00, -0.3060E+01, &
      -0.3500E+01, &
      0.7800E+00,  0.8100E+00,  0.8400E+00,  0.8900E+00,  0.1060E+01, &
      0.1110E+01,  0.1130E+01,  0.1190E+01, -0.1040E+01, -0.2900E+01, &
      -0.3440E+01, &
      0.3000E+00,  0.3200E+00,  0.3400E+00,  0.3800E+00,  0.4800E+00, &
      0.5300E+00,  0.5800E+00,  0.5300E+00, -0.1330E+01, -0.2040E+01, &
      -0.2830E+01, &
      -0.2700E+00, -0.2600E+00, -0.2400E+00, -0.2200E+00, -0.1800E+00, &
      -0.1300E+00, -0.5000E-01,  0.1000E+00, -0.1170E+01, -0.1100E+01, &
      -0.2370E+01, &
      -0.5500E+00, -0.5400E+00, -0.5300E+00, -0.5100E+00, -0.4800E+00, &
      -0.4100E+00, -0.2600E+00,  0.1700E+00, -0.4200E+00, -0.2200E+00, &
      -0.2080E+01, &
      -0.2700E+00, -0.2600E+00, -0.2500E+00, -0.2400E+00, -0.2100E+00, &
      -0.1600E+00, -0.1000E-01,  0.8500E+00,  0.9700E+00,  0.3500E+00, &
      -0.2250E+01, &
      0.5300E+00,  0.5400E+00,  0.5600E+00,  0.5700E+00,  0.6000E+00, &
      0.6500E+00,  0.7600E+00,  0.1960E+01,  0.2280E+01,  0.3600E+00, &
      -0.2590E+01, &
      0.1590E+01,  0.1630E+01,  0.1650E+01,  0.1680E+01,  0.1720E+01, &
      0.1780E+01,  0.2010E+01,  0.3260E+01,  0.3110E+01,  0.1600E+00, &
      -0.2580E+01, &
      0.2560E+01,  0.2620E+01,  0.2660E+01,  0.2690E+01,  0.2740E+01, &
      0.2830E+01,  0.3400E+01,  0.4030E+01,  0.3030E+01, -0.7000E-01, &
      -0.2320E+01, &
      0.3500E+01,  0.3600E+01,  0.3650E+01,  0.3700E+01,  0.3750E+01, &
      0.3860E+01,  0.4580E+01,  0.4100E+01,  0.2450E+01,  0.6000E-01, &
      -0.1770E+01, &
      0.4500E+01,  0.4640E+01,  0.4700E+01,  0.4760E+01,  0.4830E+01, &
      0.4930E+01,  0.5420E+01,  0.3960E+01,  0.2000E+01,  0.5000E+00, &
      -0.1150E+01, &
      0.5290E+01,  0.5470E+01,  0.5550E+01,  0.5620E+01,  0.5690E+01, &
      0.5790E+01,  0.6070E+01,  0.4000E+01,  0.1910E+01,  0.9700E+00, &
      -0.5600E+00 &
      /
      data theta_i / &
      0.2936E+03,  0.2936E+03,  0.2937E+03,  0.2937E+03,  0.2938E+03, &
      0.2942E+03,  0.2948E+03,  0.2980E+03,  0.3027E+03,  0.3092E+03, &
      0.3186E+03, &
      0.2937E+03,  0.2937E+03,  0.2937E+03,  0.2938E+03,  0.2939E+03, &
      0.2942E+03,  0.2946E+03,  0.2978E+03,  0.3024E+03,  0.3090E+03, &
      0.3184E+03, &
      0.2939E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03, &
      0.2941E+03,  0.2944E+03,  0.2976E+03,  0.3023E+03,  0.3089E+03, &
      0.3182E+03, &
      0.2940E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03, &
      0.2941E+03,  0.2943E+03,  0.2975E+03,  0.3022E+03,  0.3087E+03, &
      0.3181E+03, &
      0.2940E+03,  0.2940E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03, &
      0.2940E+03,  0.2942E+03,  0.2974E+03,  0.3021E+03,  0.3086E+03, &
      0.3180E+03, &
      0.2941E+03,  0.2940E+03,  0.2940E+03,  0.2940E+03,  0.2941E+03, &
      0.2941E+03,  0.2941E+03,  0.2973E+03,  0.3019E+03,  0.3085E+03, &
      0.3179E+03, &
      0.2942E+03,  0.2942E+03,  0.2942E+03,  0.2942E+03,  0.2941E+03, &
      0.2941E+03,  0.2941E+03,  0.2973E+03,  0.3020E+03,  0.3086E+03, &
      0.3179E+03, &
      0.2943E+03,  0.2943E+03,  0.2943E+03,  0.2943E+03,  0.2943E+03, &
      0.2943E+03,  0.2943E+03,  0.2975E+03,  0.3022E+03,  0.3087E+03, &
      0.3181E+03, &
      0.2945E+03,  0.2945E+03,  0.2945E+03,  0.2945E+03,  0.2945E+03, &
      0.2944E+03,  0.2946E+03,  0.2978E+03,  0.3025E+03,  0.3090E+03, &
      0.3184E+03, &
      0.2947E+03,  0.2947E+03,  0.2947E+03,  0.2947E+03,  0.2946E+03, &
      0.2946E+03,  0.2949E+03,  0.2980E+03,  0.3027E+03,  0.3093E+03, &
      0.3187E+03, &
      0.2949E+03,  0.2949E+03,  0.2949E+03,  0.2948E+03,  0.2948E+03, &
      0.2948E+03,  0.2950E+03,  0.2982E+03,  0.3028E+03,  0.3094E+03, &
      0.3188E+03, &
      0.2950E+03,  0.2950E+03,  0.2950E+03,  0.2950E+03,  0.2950E+03, &
      0.2950E+03,  0.2950E+03,  0.2982E+03,  0.3029E+03,  0.3095E+03, &
      0.3188E+03 &
      /
!
      save time_g, zg, u_i, v_i, theta_i
!
! --------- find time location of initial profiles 
!
      call lterp(nt,time_g,t_factor,jt,jtp1,t_weit)
!
      do iz=1,nz
         ui_temp(iz) = u_i(iz,jt)*(1.0 - t_weit) + &
      u_i(iz,jtp1)*t_weit
         vi_temp(iz) = v_i(iz,jt)*(1.0 - t_weit) + &
      v_i(iz,jtp1)*t_weit
         ti_temp(iz) = theta_i(iz,jt)*(1.0 - t_weit) + &
      theta_i(iz,jtp1)*t_weit
      enddo
!
! ----------- interpolate vertically
!
      do iz=izs,ize
         call lterp(nz,zg,zz(iz),kk,kkp1,weit)
         u_temp = ui_temp(kk)*(1.0 - weit) + &
      ui_temp(kkp1)*weit
         v_temp = vi_temp(kk)*(1.0 - weit) + &
      vi_temp(kkp1)*weit
         theta_temp = ti_temp(kk)*(1.0 - weit) + &
      ti_temp(kkp1)*weit
!
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz)   = u_temp
            v(ix,iy,iz)   = v_temp
            t(ix,iy,1,iz) = theta_temp
            w(ix,iy,iz)   = 0.
            r1(ix,iy,iz)  = 0.
            r2(ix,iy,iz)  = 0.
            r3(ix,iy,iz)  = 0.
            r4(ix,iy,1,iz)= 0.
            r5(ix,iy,iz)  = 0.
         enddo 
         enddo 
      enddo
!
! ------------- set initial random field to be
!               divergence free
!
      idum = -1
      do iz=izs,ize
         if (iz.le.8) then
!
! ----------- ampv and ampt are max amplitudes of random 
!             velocity and temperature fields
!
         ampv = 0.5
         ampt = 0.1
!  
! ------- simple random field scaled between 0 and 1
!
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = ran1(idum)
         enddo
         enddo
!
         do iy=iys,iye
         do ix=1,nnx
            psix(ix,iy) = psi(ix,iy)
            psiy(ix,iy,izs) = psi(ix,iy)
         enddo
         enddo
         call xderivp(psix(1,iys),trigx(:,1),xk,nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(:,2),yk, &
      nnx,nny,ixs,ixe,ix_s,ix_e, &
      iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
!
         vmaxx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
            if(vmag .gt. vmaxx) vmaxx = vmag
         enddo
         enddo
         facv = ampv/vmaxx
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz)   = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
            v(ix,iy,iz)   = v(ix,iy,iz) + psix(ix,iy)*facv
            t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
            e(ix,iy,iz)   = 1.0
         enddo
         enddo
         endif
!
! ---------- check divergence of initial field
!
         do iy=iys,iye
         do ix=1,nnx
            uxx(ix,iy)     = u(ix,iy,iz)
            vyy(ix,iy,izs) = v(ix,iy,iz)
         enddo
         enddo
         call xderivp(uxx(1,iys),trigx(:,1),xk,nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(:,2),yk, &
      nnx,nny,ixs,ixe,ix_s,ix_e, &
      iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         do iy=iys,iye
         do ix=1,nnx
            divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
         enddo
         enddo
         divz(iz) = divz(iz)*fnxy
!
! -------- end z loop
!
      enddo
!
      call mpi_sum_z(divz,i_root,myid,nnz,1)
!
      write(nprt,6000)
 6000 format(' check of divergence for initial state',/, &
      ' iz ',5x,' divergence')
      write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100 format(i5,e15.6)
!
! ------------ fix for baroclinic and subsidence effects !!
!
!     do iz=izs,ize
!        ug(iz)=ugcont
!        vg(iz)=vgcont
!        if (.not.(ibrcl.eq.1)) go to 19988
!        if (.not.(iz.le.izi)) go to 19987
!        ug(iz)=0.
!        vg(iz)=0.
! 19987    continue
! 19988    continue
!        zz2=zz(iz)
!        wls(iz)=-divgls*zz2
!        if (.not.(iz.eq.1)) go to 19986
!        do ix=1,nnx
!        uls(ix)=divgls*(dx*float(ix-1)-xl*.5)
!        enddo
!     enddo
!     write(nprt,9)(uls(ix),ix=1,nnx)
!  9  format(1x,8e12.3)
! 19986 continue
!
      return
      end

! --- extracted from les.F: randoc ---
      subroutine randoc
      implicit real(a-h,o-z), integer(i-n)
!
! -------- random initial conditions for an
!          ocean simulation
!
      real psi(nnx,iys:iye), psix(nnx,iys:iye), &
      psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye), &
      vyy(nnx,iys:iye,izs:izs)
!
      izi=(5*nnz)/20
      zi=z(izi)
      tmixed = 283.0
      do iz=izs,ize
         if (iz.le.izi) then
            do iy=iys,iye
            do ix=1,nnx
               u(ix,iy,iz)   = ugcont-ugal
               v(ix,iy,iz)   = vgcont
               w(ix,iy,iz)   = 0.0
               t(ix,iy,1,iz) = tmixed
               e(ix,iy,iz)   = 0.0
            enddo
            enddo
         endif
         if (iz.gt.izi) then
            do iy=iys,iye
            do ix=1,nnx
               u(ix,iy,iz)   = ugcont-ugal
               v(ix,iy,iz)   = vgcont
               w(ix,iy,iz)   = 0.0
               t(ix,iy,1,iz) = tmixed + dtdzf(1)*(zz(iz)-zi)
               e(ix,iy,iz)   = 0.0
            enddo
            enddo
         endif
         do iy=iys,iye
         do ix=1,nnx
            w(ix,iy,iz)    = 0.0
            r1(ix,iy,iz)   = 0.0
            r2(ix,iy,iz)   = 0.0
            r3(ix,iy,iz)   = 0.0
            r4(ix,iy,1,iz) = 0.0
            r5(ix,iy,iz)   = 0.0
         enddo
         enddo
      enddo
!
! ------------- set initial random field to be
!               divergence free
!
      idum = -1
      do iz=izs,ize
      if (iz.le.4) then
!
! ----------- ampv and ampt are max amplitudes of random 
!             velocity and temperature fields
!
         ampv = 0.01
!        ampt = 0.00
         ampt = 0.0001
!  
! ------- simple random field scaled between 0 and 1
!
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = ran1(idum)
         enddo
         enddo
!
         do iy=iys,iye
         do ix=1,nnx
            psix(ix,iy) = psi(ix,iy)
            psiy(ix,iy,izs) = psi(ix,iy)
         enddo
         enddo
         call xderivp(psix(1,iys),trigx(:,1),xk,nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(:,2),yk, &
      nnx,nny,ixs,ixe,ix_s,ix_e, &
      iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
!
         vmaxx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
            if(vmag .gt. vmaxx) vmaxx = vmag
         enddo
         enddo
         facv = ampv/vmaxx
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
            v(ix,iy,iz) = v(ix,iy,iz) + psix(ix,iy)*facv
            t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
            e(ix,iy,iz) = 0.0001
         enddo
         enddo
      endif
!
! ---------- check divergence of initial field
!
         do iy=iys,iye
         do ix=1,nnx
            uxx(ix,iy) = u(ix,iy,iz)
            vyy(ix,iy,izs) = v(ix,iy,iz)
         enddo
         enddo
         call xderivp(uxx(1,iys),trigx(:,1),xk,nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(:,2),yk, &
      nnx,nny,ixs,ixe,ix_s,ix_e, &
      iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         do iy=iys,iye
         do ix=1,nnx
            divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
         enddo
         enddo
         divz(iz) = divz(iz)*fnxy
!
! -------- end z loop
!
      enddo
!
      call mpi_sum_z(divz,i_root,myid,nnz,1)
!
      write(nprt,6000)
 6000 format(' check of divergence for initial state',/, &
      ' iz ',5x,' divergence')
      write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100 format(i5,e15.6)
!
      do iz=izs,ize
         ug(iz)=ugcont
         vg(iz)=vgcont
      enddo
!
      return
      end

! --- extracted from les.F: gridd ---
      subroutine gridd
      implicit real(a-h,o-z), integer(i-n)
!
! ----------- allocate space and pass arrays using modules
!

      if (myid==0) write(6,5001) isize
 5001 format(' size of stats array = ',i8)

!
!
! ---------------- debug for arrays
!
      big = -99.0e+300
!
! ---------------- setup grid
!
      nnx = maxnx
      nny = maxny
      nnz = maxnz
!     izs = 1
!     ize = nnz
!
!
! ----------- make sure problem and cpu's match
!
      maxp   = numprocs-1
      ncpu_z = numprocs/ncpu_s
      if(mod(numprocs,ncpu_s) .ne. 0 .or. &
      ncpu_z .gt. nnz) then
         go to 999
      endif
      if(l_root) write(6, 1100) ncpu_s, ncpu_z, numprocs, &
      maxp
      write(nprt,1100) ncpu_s, ncpu_z, numprocs, maxp
 1100 format(' Number of x-y slab cpus = ',i5,/, &
      ' Number of z-level cpus  = ',i5,/, &
      ' Total number of cpus    = ',i5,/, &
      ' Max-p for index arrays  = ',i5)
!
! ---------------- allocate arrays for (i,j,k)-indexing on
!                  each processor (see set_range)
!
      allocate(ix_s(0:maxp), ix_e(0:maxp), &
      jx_s(0:maxp), jx_e(0:maxp), &
      kx_s(0:maxp), kx_e(0:maxp), &
      mx_s(0:maxp), mx_e(0:maxp), &
      iy_s(0:maxp), iy_e(0:maxp), &
      jy_s(0:maxp), jy_e(0:maxp), &
      is_s(0:maxp), is_e(0:maxp), &
      iz_s(0:maxp), iz_e(0:maxp))
!
! ---------------- setup array sizes and variable dimensions
!
      nxy   = nnx*nny
      ncx   = nnx/2 + 1
      ncy   = nny/2 + 1
      nnxp1 = nnx + 1
      nnyp1 = nny + 1
      nnxp2 = nnx + 2
      nnyp2 = nny + 2
      nnzp1 = nnz + 1
      nnzm1 = nnz - 1
      ivis = ivis0
      fnxy  = 1.0/float(nnx*nny)
!
      write(nprt,7001) nnx,nny,nnz
 7001 format(' 7001 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)
!
      call set_range
!
      write(nprt,7002) nnx,nny,nnz
 7002 format(' 7002 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)
!
      num_y = iye + 1 - iys
!
! ------------- allocate solution arrays
!               account for nnxp2 for fields but not in rhs
!               and possible monotone for scalars
!
      allocate(u(nnxp2,iys:iye,izs-1:ize+1),  &
      v(nnxp2,iys:iye,izs-1:ize+1),  &
      w(nnxp2,iys:iye,izs-1:ize+1),  &
      t(nnxp2,iys:iye,nscl,izs-2:ize+2),  &
      e(nnxp2,iys:iye,izs-1:ize+1),  &
      r1(nnx,iys:iye,izs-1:ize+1), &
      r2(nnx,iys:iye,izs-1:ize+1), &
      r3(nnx,iys:iye,izs-1:ize+1), &
      r4(nnx,iys:iye,nscl,izs-1:ize+1), &
      r5(nnx,iys:iye,izs-1:ize+1))
      u = 0.0
      v = 0.0
      w = 0.0
      t = 0.0
      e = 0.0
      r1 = 0.0
      r2 = 0.0
      r3 = 0.0
      r4 = 0.0
      r5 = 0.0

      allocate(p_base(0:nnz+1), &
      rho_base(0:nnz+1), &
      T_base(0:nnz+1), &
      theta_base(0:nnz+1))
      p_base = 0.0
      rho_base = 0.0
      T_base = 0.0
      theta_base = 0.0

!   ------Allocate arrays for SFS stress and its derivatives
      allocate(sigm_s(nnxp2,iys:iye,izs-1:ize+1), &
      sigm_sdx(nnxp2,iys:iye,izs-1:ize+1), &
      sigm_sdy(nnxp2,iys:iye,izs-1:ize+1), &
      sigm_sdz(nnxp2,iys:iye,izs-1:ize+1), &
      vis_ss(nnxp2,iys:iye,izs-1:ize+1))

      allocate (sigm_sext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3), &
      sigm_sdxext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3), &
      sigm_sdyext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3), &
      sigm_sdzext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3), &
      vis_sext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3))

!
! ------------- allocate extended arrays for interpolation of
!               particle/spray location
!
      allocate(uext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),  &
      vext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),  &
      wext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3), &
      Text(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3), &
      T2ext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3)) 
      uext = 0.0
      vext = 0.0
      wext = 0.0
      Text = 0.0
      T2ext = 0.0

      !Transposed velocities to do the uf interpolation:
      allocate(u_t(0:nnz+1,iys:iye,mxs:mxe), &
      v_t(0:nnz+1,iys:iye,mxs:mxe), &
      w_t(0:nnz+1,iys:iye,mxs:mxe),  &
      T_t(0:nnz+1,iys:iye,mxs:mxe), &
      T2_t(0:nnz+1,iys:iye,mxs:mxe))
      u_t = 0.0
      v_t = 0.0
      w_t = 0.0
      T_t = 0.0
      T2_t = 0.0

      !Column-oriented array to calculate max. droplet relaxation time
      allocate(partsrc(nnx,iys:iye,izs-1:ize+1,1:3))
      partsrc = 0.0
      allocate(partsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1,1:3))
      partsrc_t = 0.0
      allocate(partTsrc(nnx,iys:iye,izs-1:ize+1))
      partTsrc = 0.0
      allocate(partTsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1))
      partTsrc_t = 0.0
      allocate(partHsrc(nnx,iys:iye,izs-1:ize+1))
      partHsrc = 0.0
      allocate(partHsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1))
      partHsrc_t = 0.0
      allocate(partTEsrc(nnx,iys:iye,izs-1:ize+1))
      partTEsrc = 0.0
      allocate(partTEsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1))
      partTEsrc_t = 0.0

      allocate(LWP(mxs:mxe,iys:iye))
      LWP = 0.0

      allocate(surf_precip(mxs:mxe,iys:iye))
      surf_precip = 0.0

      radsrc = 0.0
      viz_t_elapsed = 1.0e-10 !Initialize very small to prevent blowing up on 1st time step

! ------------- allocate space for boundary condition arrays
!               on top and bottom of domain
!
      allocate(ubc(nnx,iys:iye,2), &
      vbc(nnx,iys:iye,2), &
      wbc(nnx,iys:iye,2), &
      tbc(nnx,iys:iye,nscl,2), &
      ebc(nnx,iys:iye,2), &
      pbc(nnx,iys:iye,2), &
      pbc2(nnx,iys:iye,2))
!
! ------------ allocate space for wind and surface arrays
!
      allocate(wind(nnx,iys:iye),  &
      tau13m(nnx,iys:iye),  &
      tau23m(nnx,iys:iye),  &
      taut3m(nnx,iys:iye,nscl),  &
      t_grnd(nnx,iys:iye,nscl))
!
! ------------------- allocate space for derivative arrays
!
      allocate(ux(nnx,iys:iye,izs-1:ize+1), &
      uy(nnx,iys:iye,izs-1:ize+1), &
      vx(nnx,iys:iye,izs-1:ize+1), &
      vy(nnx,iys:iye,izs-1:ize+1), &
      wx(nnx,iys:iye,izs-1:ize+1), &
      wy(nnx,iys:iye,izs-1:ize+1))
!
! ------------- allocate space for pressure, pressure bcs
!
      allocate(p(nnxp2,iys:iye,izs-1:ize+1), &
      ptop(nnxp2,iys:iye,2))
      p = 0.0
      ptop = 0.0
!
! ------------- allocate space for viscosity and diffusivity
!
      allocate(vis_m(nnx,iys:iye,izs-1:ize+1), &
      vis_s(nnx,iys:iye,nscl,izs-1:ize+1))
!
! ------------- allocate space for fft trig factors
!
      nq_trig = max(nnx,nny)
      allocate(trigx(2*nq_trig+15,2), &
      trigc(4*nq_trig+15))
!
      return
  999 continue
!
      if(l_root) write(6,1000) numprocs, ncpu_s, mmz
      write(nprt,1000) numprocs, ncpu_s, nnz
 1000 format(' Gridd Trouble number of processors and grid', &
      ' partitioning do not match!',/, &
      ' Total num of cpus   = ',i5, &
      ' Num cpu on x-y slab = ',i5,/, &
      ' Num of z-levels     = ',i5)
      call mpi_finalize(ierr)
      end

! --- extracted from les.F: restart ---
      subroutine restart
      implicit real(a-h,o-z), integer(i-n)
!
! ----------- get restart file from local directory
!
      character*80 path_res_c
      logical there
!
! --------------------- check if file is there
!
      inquire(file=path_res,exist=there)
      if(there) then
         if(l_root) write(6,6001) path_res
      else
         if(l_root) write(6,6005) path_res
         stop
      endif
!
! ------------------ get constant file
!
      iloc = index(path_res,' ')
      path_res_c = path_res(1:iloc-1)//'.con'
      inquire(file=path_res_c,exist=there)
      if(there) then
         if(l_root) write(6,6002) path_res_c
      else
         if(l_root) write(6,6006) path_res_c
         stop
      endif
      open(nvelc,err=200,file=path_res_c,form='unformatted', &
      status='old')
!
      call read_res
!
      return
! ---------------------------- process errors
  100 continue
      write(6,9000) path_res, nvel
      call mpi_finalize(ierr)
      stop
! -----------------------
  200 continue
      write(6,9001) path_res_c, nvelc
      call mpi_finalize(ierr)
      stop
! -----------------------
 6001 format(' SR. RESTART: FILE READ = ',A80)
 6002 format(' SR. RESTART: CONSTANT FILE READ = ',A80)
 6005 format(' 6005, SR. RESTART: cannot find restart file = ',a80)
 6006 format(' 6005, SR. RESTART: cannot find constant file = ',a80)
 9000 format(' 9000, SR. RESTART: cannot open file =',a80,/, &
      ' to unit number = ',i2)
 9001 format(' 9001, SR. RESTART: cannot open file =',a80,/, &
      ' to unit number = ',i2)
      end

! --- extracted from les.F: read_res ---
      subroutine read_res
      implicit real(a-h,o-z), integer(i-n)
!
! -------------- read restart file including constant file
!                changed for iys:iye
!
#if defined(SWAP)
      use module_byteswap
#endif
!
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
      real, allocatable, dimension(:,:,:) :: temp
      allocate(temp(nvar,nnx,iys:iye))
!
! ---- open file
!
      call mpi_file_open(mpi_comm_world, path_res, &
      mpi_mode_create+mpi_mode_rdwr, &
      mpi_info_null, nvel, ierr)
!
! ---- set file view
!
      disp = 0
      call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8, &
      'native',mpi_info_null,ierr)
!
! ------------ read 3d fields
!
      nsize  = int(nvar,k8)*nnx*nny
      nsize2 = int(nvar,k8)*nnx*(iys-1)
      n_read = nvar*nnx*(iye+1-iys)
!
      do k=izs,ize
         offset = int((k-1),k8)*nsize + nsize2
         call mpi_file_read_at_all(nvel,offset,temp,n_read, &
      mpi_real8,status,ierr)
         if (ierr /= 0) goto 9992
#if defined(SWAP)
         call byteswap(temp)
#endif
         do j=iys,iye
         do i=1,nnx
            u(i,j,k) = temp(1,i,j) 
            v(i,j,k) = temp(2,i,j)
            w(i,j,k) = temp(3,i,j)
            e(i,j,k) = temp(nvar,i,j)
         enddo
         enddo
         do is = 1,nscl
            do j = iys,iye
            do i = 1,nnx
               t(i,j,is,k) = temp(3+is,i,j)
            enddo
            enddo
         enddo
!
      enddo
!
! ---- close file
!
      call mpi_file_close(nvel, ierr)
!
      deallocate(temp)
!
! ------------ every mpi process reads constant file
!
      rewind(nvelc)
      read(nvelc,err=9993) c_c, c_s, c_i, case
      close(nvelc)
!
      if(l_root) write(6,4001) case
 4001 format(' 4001, SR. RESTART: case from restart = ',a3)
!
! ----- special restart conditions -------------------------------------
!
! -------- set case name to case input
!
      case   = case_inp
      if(l_root) write(6,4002) case_inp, utau, utausv
 4002 format(' 4002, SR. RESTART:',/, &
      ' files will be saved with case name = ',a3,/, &
      ' utau = ',e15.6,' utausv = ',e15.6)
!
! ------------------- if new vis model set match point for
!                     outer grid
      nmatch = 48
      utau = utausv
!
! -------- hand coded changes to restart if needed
!
!       qstars = 0.000
!       wtsfcs = 0.000
!
!
! ------   reset qstar and wtsfc for no heat flux
!              qstar(1) = qstars
!              wtsfc(1) = wtsfcs
!              qstar(2) = qstars
!              wtsfc(2) = wtsfcs
! ------   redefine case id to input value
!              case = cases
!
      if(l_root) write(6,4012) time
      if(l_root) write(6,4013) qstar(1) , nmatch, case
!
      call get_dz
!
      return
! ------------------------  process errors from read
!9991 continue
!     write(6,6000) nvel,iz
!6000 format(' SR. READ_RES: hit end of file on unit number = ',i2,/,
!    +       '               at iz = ',i4)
!     call mpi_finalize(ierr)
!     stop
! ---------------------
 9992 continue
      write(6,6100) nvel,iz
 6100 format(' SR. READ_RES: error reading file on unit number = ',i2,/, &
      '               at iz = ',i4)
      call mpi_finalize(ierr)
      stop
! ---------------------
 9993 continue
      write(6,6200) nvelc
 6200 format(' SR. READ_RES:',/, &
      '    error reading constant file on unit number = ',i2)
      call mpi_finalize(ierr)
      stop
! ---------------------
 4012 format(' SR. RESTART: restart completed at T=',e15.6)
 4013 format('    after restart qstar = ',e15.6,' nmatch = ',i5, &
      ' case = ',a3)
      end

end module mod_init
