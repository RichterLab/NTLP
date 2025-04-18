c     ----------------------------------------------------------------------
      module pars
c     ----------------------------------------------------------------------
      integer :: iti,itmax,imean,itape
      integer :: itstr
      real :: i_viz,ihst

      integer :: ismlt,ifree,iradup,
     +     iupwnd,ibuoy,ifilt,itcut,isubs,ibrcl,iocean,method,
     +     idebug,iz_space,ivis0,ifix_dt,new_vis,iDNS,
     +     icouple,iTcouple,iHcouple,ievap,ifields,ilin,ineighbor,
     +     icoalesce,ipart_method,ipartdiff,isfs,iexner,ilongwave,
     +     inetcdf,iviznetcdf,ihumiditycontrol,ireintro,itrajout,
     +     inewpart,ihurr,icase,ikernel,irestart

      real :: max_time
c     
c     -------- j_recl = record length in "bytes" for history file
c     k8     = kind parameter for integers in mpi_io routines
c     
      integer, parameter :: j_recl=4, k8=8

      integer, parameter :: nscl = 2, nvar = (4+nscl)
      integer, parameter :: maxnx=128,maxny=128,maxnz=128
      integer, parameter :: maxnz1=maxnz+1,
     + maxnz2=maxnz+2,maxnx2=maxnx+2,maxny2=maxny+2

      !Boundary conditions and other parameters
      real :: Uo, Ttop(nscl), Tbot(nscl),qstar(nscl),tsfcc(nscl)
      real :: grav,t00,fcor,zi,zl,xl,yl,zw1,dpdx
      real :: ugcont,vgcont,zo,zos,surf_RH,surf_p,surf_rho
      real :: meanRH,varRH
      real :: nprime
      real :: dvdr,hurr_rad
      integer :: isfc(nscl)
      real :: viz_t_elapsed

      !Radiation parameters
      real :: rad_kappa,rad_Fct,rad_Fcb
      real :: radflux(0:maxnz),radsrc(maxnz)

      !Material parameters, both air and particles
      real :: Lv = (25.0 - 0.02274*26.0)*10.0**5
      real :: Sal,Ms,Gam,Ion
      real :: rhow,nuf,rhos
      real :: CpaCpp,Cpa,Cpp,Pra,Sc,Mw,Ru,Rd

c     
c     ------------ leave following definitions as is
c     
c     ----------------------------------------------------------------------
      integer ::    nnx, nny, nnz, nxy, ncx, nnxp1, nnyp1, ncy,
     +     nnxp2, nnyp2, nnzp1, ivis, nnzm1, isize, krec,
     +     izs, ize, ixs, ixe, jxs, jxe, kxs, kxe,
     +     mxs, mxe, iss, ise, iys, iye, jys, jye
c     ----------------------------------------------------------------------
      character case*3
c     ----------------------------------------------------------------------
      integer  ::   nvel, npre, nhis1, nprt,nrad,nres,nactres,
     +     nhisp, nvelc, ntraj,
     +     nviz_z, nviz_y, 
     +     nviz_x, nviz_s,
     +     kfile, jfile, ibcl, ibcu,
     +     igrdr, imach, it_his_nxt, it_viz_nxt,
     +     it_his,it_viz
      logical ::    mnout, micut, mtape, mhis, msave,
     +     l_root, l_debug, msave_v, mviz
c     ----------------------------------------------------------------------
      real    ::    windm,u1xy,v1xy,t1xy(nscl),
     +     t10xy(nscl),au13m,au23m,aut3m(nscl),tsfcm(nscl),
     +     thstar(nscl), eavg(maxnz), tr_tau(0:maxnz),
     +     pxym(0:maxnz1), zi_min
      integer ::    izi, iz_min
      real, allocatable ::
     +     wind(:,:), tau13m(:,:), tau23m(:,:), 
     +     taut3m(:,:,:), t_grnd(:,:,:)
c     ----------------------------------------------------------------------
      real ::       u_mn(0:maxnz1), v_mn(0:maxnz1),
     +     w_mn(0:maxnz1), t_mn(0:maxnz1,nscl), alphaC(0:maxnz1,nscl),
     +     e_mn(0:maxnz1)
c     ----------------------------------------------------------------------
      real ::       dzw(0:maxnz2), dzu(0:maxnz2),
     +     dzw_i(0:maxnz2), dzu_i(0:maxnz2)
c     ----------------------------------------------------------------------
      real ::       t_factor, t_ref, c_rate, t_surf_i
c     ----------------------------------------------------------------------
      real ::       dfac(maxnz), dsl_z(0:maxnz1),
     +     xksurf, viscon, vise, almin_c,stabmin,
     +     ck,ceps,csmag,stab_c,vis_mean(0:maxnz)
      integer ::    nmatch
c     ----------------------------------------------------------------------
      real ::       zetas(3), gama(3), etas(4), dt_new,
     +     umax,vmax,wmax, wabs, vismax,
     +     cfl, tzero,
     +     ucfl, vcfl, wcfl
c     ----------------------------------------------------------------------
      character*80  path_res, path_prt,
     +     path_hp, path_sav_hp, path_part, path_ran,
     +     path_v, path_c, path_p, path_h,
     +     path_sav_v, path_sav_c,
     +     path_sav_p, path_sav_h, path_sav_part,
     +     path_rad_hist, path_res_hist, path_actres_hist,
     +     path_traj,
     +     bad_news
      character case_inp*3
      character*80 path_viz_xy, path_viz_xz, path_viz_yz, path_stuf,
     +             path_seed
c     ----------------------------------------------------------------------
      integer ::    myid, numprocs, i_root,
     +     ziloc, myid_newvis, ncpu_s, ncpu_z, maxp
      integer, allocatable, dimension(:) :: 
     +     ix_s, ix_e, jx_s, jx_e,
     +     kx_s, kx_e, mx_s, mx_e,
     +     iy_s, iy_e, jy_s, jy_e,
     +     is_s, is_e, iz_s, iz_e
      end module pars
c     ----------------------------------------------------------------------
c     --------------------------------------------------------------------- 
      module fields
      real, allocatable :: 
     +     u(:,:,:), v(:,:,:), w(:,:,:), t(:,:,:,:), e(:,:,:), 
     +     r1(:,:,:), r2(:,:,:), r3(:,:,:), r4(:,:,:,:), r5(:,:,:)
      real, allocatable :: 
     +     ux(:,:,:), uy(:,:,:), vx(:,:,:), vy(:,:,:), 
     +     wx(:,:,:), wy(:,:,:),
     +     p(:,:,:), ptop(:,:,:), vis_m(:,:,:), vis_s(:,:,:,:) 
      real, allocatable :: 
     +     ubc(:,:,:), vbc(:,:,:), wbc(:,:,:), tbc(:,:,:,:), 
     +     ebc(:,:,:), pbc(:,:,:), pbc2(:,:,:)
      real, allocatable ::
     +     p_base(:),rho_base(:),T_base(:),theta_base(:)
      end module fields
c     ----------------------------------------------------------------------
      module fftwk
      real, allocatable :: trigx(:,:), trigc(:)
      end module fftwk
c     ----------------------------------------------------------------------
      module con_data
c     ----------------------------------------------------------------------
      use pars, only : nscl
      type con_d
      sequence
      real ::  vk, vkin, vk74, vk74in, 
     +     gcp, fcor_h, pi2,
     +     batagk, batag, bfac,
     +     cdbtm, dtdzf(nscl), dtjump, ugal, divgls,
     +     z1, utausv, dx, dy, dz, dt, 
     +     fnxy, dzdz, dsl, c23, dtgama, dtzeta, xkmax,
     +     time, zody, zody74, zosdy,
     +     utau, wtsfc(nscl),
     +     uwsfc, vwsfc, amonin,
     +     zol, hol, smal_e, sml_eg, u10, Swall
      end type con_d
      type(con_d), target :: c_c
      real, pointer ::
     +     vk, vkin, vk74, vk74in, 
     +     gcp, fcor_h, pi2,
     +     batagk, batag, bfac,
     +     cdbtm, dtdzf(:), dtjump, ugal, divgls,
     +     z1, utausv, dx, dy, dz, dt,
     +     fnxy, dzdz, dsl, c23, dtgama, dtzeta, xkmax, 
     +     time, zody, zody74, zosdy,
     +     utau, wtsfc(:), 
     +     uwsfc, vwsfc, amonin,
     +     zol, hol, smal_e, sml_eg, u10, Swall

      !Integers
      type con_i
      sequence
      integer :: his_counter,histog_counter,viz_counter,itn,it,
     + tnum_destroy_accum
      end type con_i

      type(con_i), target :: c_i
      integer, pointer :: 
     + his_counter,histog_counter,viz_counter,itn,it,
     + tnum_destroy_accum
       
   
      contains
      subroutine fill_cc
c     
c     --------------- pointer associations for constant variables
c     
      vk     => c_c%vk 
      vkin   => c_c%vkin
      vk74   => c_c%vk74
      vk74in => c_c%vk74in
      gcp    => c_c%gcp
      fcor_h => c_c%fcor_h
      pi2    => c_c%pi2
      batagk => c_c%batagk
      batag  => c_c%batag
      bfac   => c_c%bfac
      cdbtm  => c_c%cdbtm
      dtdzf  => c_c%dtdzf
      dtjump => c_c%dtjump
      ugal   => c_c%ugal
      divgls => c_c%divgls
      z1     => c_c%z1
      utausv => c_c%utausv
      dx     => c_c%dx
      dy     => c_c%dy
      dz     => c_c%dz
      dt     => c_c%dt
      fnxy   => c_c%fnxy
      dzdz   => c_c%dzdz
      dsl    => c_c%dsl
      c23    => c_c%c23
      dtgama => c_c%dtgama
      dtzeta => c_c%dtzeta
      xkmax  => c_c%xkmax
      time   => c_c%time
      zody   => c_c%zody
      zody74 => c_c%zody74
      zosdy   => c_c%zosdy
      utau   => c_c%utau
      wtsfc  => c_c%wtsfc
      uwsfc  => c_c%uwsfc
      vwsfc  => c_c%vwsfc
      amonin => c_c%amonin
      zol    => c_c%zol
      hol    => c_c%hol
      smal_e => c_c%smal_e
      sml_eg => c_c%sml_eg
      u10    => c_c%u10
      Swall  => c_c%Swall
      return
      end subroutine fill_cc
      subroutine fill_ci

      his_counter    => c_i%his_counter
      histog_counter => c_i%histog_counter
      viz_counter    => c_i%viz_counter
      itn            => c_i%itn
      it             => c_i%it
      tnum_destroy_accum   => c_i%tnum_destroy_accum

      end subroutine fill_ci
      end module con_data
c ----------------------------------------------------------------------
      module con_stats
        use pars
        type con_s
        sequence
        real ::  wwsb(maxnz),engz(0:maxnz1),
     +           engsbz(0:maxnz1),
     +           englez(maxnz),uxym(0:maxnz1),
     +           vxym(0:maxnz1),wxym(0:maxnz1),
     +           txym(0:maxnz1,nscl),divz(0:maxnz1),
     +           utle(maxnz,nscl), utsb(maxnz,nscl),
     +           vtle(maxnz,nscl), vtsb(maxnz,nscl),
     +           wtle(maxnz,nscl), wtsb(0:maxnz,nscl),
     +           wt_tot(maxnz,nscl),
     +           z(0:maxnz1),zz(0:maxnz1),
     +           shrz(maxnz),buyz(maxnz),
     +           triz(maxnz),
     +           uwsb(0:maxnz),vwsb(0:maxnz),
     +           uwle(maxnz),vwle(maxnz),
     +           uw_tot(maxnz),vw_tot(maxnz),
     +           wcube(maxnz), wfour(maxnz),
     +           tcube(maxnz,nscl),
     +           ups(maxnz), vps(maxnz),
     +           wps(maxnz), tps(maxnz,nscl),
     +           t_rprod(maxnz), t_wq(maxnz),
     +           t_wp(maxnz), t_tau(maxnz),
     +           t_tran(maxnz), t_buoy(maxnz),
     +           t_diss(maxnz), t_sprod(maxnz),
     +           zconc(maxnz),
     +           vp1mean(maxnz),vp2mean(maxnz),vp3mean(maxnz),
     +           vp1msqr(maxnz),vp2msqr(maxnz),vp3msqr(maxnz),
     +           uf1mean(maxnz),uf2mean(maxnz),uf3mean(maxnz),
     +           uf1msqr(maxnz),uf2msqr(maxnz),uf3msqr(maxnz),
     +           m1src(maxnz),m2src(maxnz),m3src(maxnz),
     +           upwpm(maxnz),
     +           Tpmean(maxnz),Tpmsqr(maxnz),
     +           Tfmean(maxnz),qfmean(maxnz),
     +           wpTpm(maxnz),Tpsrc(maxnz),
     +           Tv_tran(maxnz,nscl),Tv_prod(maxnz,nscl),
     +           Tv_diss(maxnz,nscl),
     +           Tv_part1(maxnz,nscl),Tv_part2(maxnz),
     +           radmean(maxnz),rad2mean(maxnz),
     +           multmean(maxnz),mwmean(maxnz), 
     +           Nc(maxnz),ql(maxnz),  
     +           Hpsrc(maxnz),TEpsrc(maxnz),
     +           qstarm(maxnz),
     +           trip(maxnz),TpSq(maxnz),Tpdqp(maxnz),
     +           radtend(maxnz),
     +           RHxym(maxnz),RHmsqr(maxnz),tempxym(maxnz),
     +           pflux(0:maxnz),pfluxdiff(0:maxnz),
     +           pmassflux(0:maxnz),penegflux(0:maxnz)
        real ::  xkn(maxnx),ykn(maxny),
     +           xk(maxnx),yk(maxny),
     +           xks(maxnx2,maxny),wavexy(maxnx2,maxny)
        real ::  ug(maxnz),vg(maxnz),
     +           wls(maxnz),uls(maxnx)
        real ::  udrift,vdrift,
     +           stokesw,stokesa,
     +           stokess,stokes(maxnz1)
        real ::  dtg, dslg, dzg
        real ::  lwc,phiw,phiv 

        end type con_s
        type(con_s), target :: c_s
        real, pointer ::
     +           wwsb(:), engz(:), engsbz(:),
     +           englez(:), uxym(:), vxym(:), wxym(:),
     +           txym(:,:), divz(:), utle(:,:), utsb(:,:),
     +           vtle(:,:), vtsb(:,:), wtle(:,:), wtsb(:,:),
     +           wt_tot(:,:), z(:), zz(:), shrz(:), buyz(:),
     +           triz(:), uwsb(:), vwsb(:), uwle(:), vwle(:),
     +           uw_tot(:), vw_tot(:), wcube(:), wfour(:),
     +           tcube(:,:), ups(:), vps(:),
     +           wps(:), tps(:,:), t_rprod(:), t_wq(:),
     +           t_wp(:), t_tau(:), t_tran(:), t_buoy(:),
     +           t_diss(:), t_sprod(:),
     +           zconc(:),
     +           vp1mean(:),vp2mean(:),vp3mean(:),
     +           vp1msqr(:),vp2msqr(:),vp3msqr(:),
     +           uf1mean(:),uf2mean(:),uf3mean(:),
     +           uf1msqr(:),uf2msqr(:),uf3msqr(:),
     +           m1src(:),m2src(:),m3src(:),
     +           upwpm(:),
     +           Tpmean(:),Tpmsqr(:),
     +           Tfmean(:),qfmean(:),
     +           wpTpm(:),Tpsrc(:),
     +           Tv_tran(:,:),Tv_prod(:,:),Tv_diss(:,:),
     +           Tv_part1(:,:),Tv_part2(:),
     +           radmean(:),rad2mean(:),
     +           multmean(:),mwmean(:),
     +           Nc(:),ql(:), 
     +           Hpsrc(:), TEpsrc(:),
     +           qstarm(:),
     +           trip(:),TpSq(:),Tpdqp(:),
     +           radtend(:),
     +           RHxym(:),RHmsqr(:),tempxym(:),
     +           pflux(:),pfluxdiff(:),
     +           pmassflux(:),penegflux(:)
        real, pointer ::  
     +           xkn(:), ykn(:), xk(:), yk(:), xks(:,:), wavexy(:,:)
        real, pointer ::  
     +           ug(:), vg(:), wls(:), uls(:)
        real, pointer ::  
     +           udrift, vdrift, stokesw, stokesa,
     +           stokess, stokes(:)
        real, pointer ::  
     +           dtg, dslg, dzg
        real, pointer ::  
     +           lwc,phiw,phiv 
      contains
         subroutine fill_cs
c
c -------------- pointer association for stat arrays and get size
c                of stat arrays isize for history files
c
             isize = 0 
             wwsb    => c_s%wwsb     ; isize = isize + size(wwsb)
             engz    => c_s%engz     ; isize = isize + size(engz)
             engsbz  => c_s%engsbz   ; isize = isize + size(engsbz)
             englez  => c_s%englez   ; isize = isize + size(englez)
             uxym    => c_s%uxym     ; isize = isize + size(uxym)
             vxym    => c_s%vxym     ; isize = isize + size(vxym)
             wxym    => c_s%wxym     ; isize = isize + size(wxym)
             txym    => c_s%txym     ; isize = isize + size(txym)
             divz    => c_s%divz     ; isize = isize + size(divz)
             utle    => c_s%utle     ; isize = isize + size(utle)
             utsb    => c_s%utsb     ; isize = isize + size(utsb)
             vtle    => c_s%vtle     ; isize = isize + size(vtle)
             vtsb    => c_s%vtsb     ; isize = isize + size(vtsb)
             wtle    => c_s%wtle     ; isize = isize + size(wtle)
             wtsb    => c_s%wtsb     ; isize = isize + size(wtsb)
             wt_tot  => c_s%wt_tot   ; isize = isize + size(wt_tot)
             z       => c_s%z        ; isize = isize + size(z)
             zz      => c_s%zz       ; isize = isize + size(zz)
             shrz    => c_s%shrz     ; isize = isize + size(shrz)
             buyz    => c_s%buyz     ; isize = isize + size(buyz)
             triz    => c_s%triz     ; isize = isize + size(triz)
             uwsb    => c_s%uwsb     ; isize = isize + size(uwsb)
             vwsb    => c_s%vwsb     ; isize = isize + size(vwsb)
             uwle    => c_s%uwle     ; isize = isize + size(uwle)
             vwle    => c_s%vwle     ; isize = isize + size(vwle)
             uw_tot  => c_s%uw_tot   ; isize = isize + size(uw_tot)
             vw_tot  => c_s%vw_tot   ; isize = isize + size(vw_tot)
             wcube   => c_s%wcube    ; isize = isize + size(wcube)
             wfour   => c_s%wfour    ; isize = isize + size(wfour)
             tcube   => c_s%tcube    ; isize = isize + size(tcube)
             ups     => c_s%ups      ; isize = isize + size(ups)
             vps     => c_s%vps      ; isize = isize + size(vps)
             wps     => c_s%wps      ; isize = isize + size(wps)
             tps     => c_s%tps      ; isize = isize + size(tps)
             t_rprod => c_s%t_rprod  ; isize = isize + size(t_rprod)
             t_wq    => c_s%t_wq     ; isize = isize + size(t_wq)
             t_wp    => c_s%t_wp     ; isize = isize + size(t_wp)
             t_tau   => c_s%t_tau    ; isize = isize + size(t_tau)
             t_tran  => c_s%t_tran   ; isize = isize + size(t_tran)
             t_buoy  => c_s%t_buoy   ; isize = isize + size(t_buoy)
             t_diss  => c_s%t_diss   ; isize = isize + size(t_diss)
             t_sprod => c_s%t_sprod  ; isize = isize + size(t_sprod)
             zconc   => c_s%zconc    ; isize = isize + size(zconc)
             vp1mean  => c_s%vp1mean ; isize = isize + size(vp1mean)
             vp2mean  => c_s%vp2mean ; isize = isize + size(vp2mean)
             vp3mean  => c_s%vp3mean ; isize = isize + size(vp3mean)
             vp1msqr  => c_s%vp1msqr ; isize = isize + size(vp1msqr)
             vp2msqr  => c_s%vp2msqr ; isize = isize + size(vp2msqr)
             vp3msqr  => c_s%vp3msqr ; isize = isize + size(vp3msqr)
             uf1mean  => c_s%uf1mean ; isize = isize + size(uf1mean)
             uf2mean  => c_s%uf2mean ; isize = isize + size(uf2mean)
             uf3mean  => c_s%uf3mean ; isize = isize + size(uf3mean)
             uf1msqr  => c_s%uf1msqr ; isize = isize + size(uf1msqr)
             uf2msqr  => c_s%uf2msqr ; isize = isize + size(uf2msqr)
             uf3msqr  => c_s%uf3msqr ; isize = isize + size(uf3msqr)
             m1src   => c_s%m1src    ; isize = isize + size(m1src)
             m2src   => c_s%m2src    ; isize = isize + size(m2src)
             m3src   => c_s%m3src    ; isize = isize + size(m3src)
             upwpm   => c_s%upwpm    ; isize = isize + size(upwpm)
             Tpmean  => c_s%Tpmean   ; isize = isize + size(Tpmean)
             Tpmsqr  => c_s%Tpmsqr   ; isize = isize + size(Tpmsqr)
             Tfmean  => c_s%Tfmean   ; isize = isize + size(Tfmean)
             qfmean  => c_s%qfmean   ; isize = isize + size(qfmean)
             wpTpm   => c_s%wpTpm    ; isize = isize + size(wpTpm)
             Tpsrc   => c_s%Tpsrc    ; isize = isize + size(Tpsrc)
             Tv_tran => c_s%Tv_tran  ; isize = isize + size(Tv_tran)
             Tv_prod => c_s%Tv_prod  ; isize = isize + size(Tv_prod)
             Tv_diss => c_s%Tv_diss  ; isize = isize + size(Tv_diss)
             Tv_part1 => c_s%Tv_part1  ; isize = isize + size(Tv_part1)
             Tv_part2 => c_s%Tv_part2  ; isize = isize + size(Tv_part2)
             radmean => c_s%radmean  ; isize = isize + size(radmean)
             rad2mean => c_s%rad2mean  ; isize = isize + size(rad2mean)
             multmean => c_s%multmean ; isize = isize + size(multmean)
             mwmean => c_s%mwmean    ; isize = isize + size(mwmean) 
             Nc => c_s%Nc    ; isize = isize + size(Nc) 
             ql => c_s%ql    ; isize = isize + size(ql)
             Hpsrc   => c_s%Hpsrc    ; isize = isize + size(Hpsrc)
             TEpsrc  => c_s%TEpsrc  ; isize = isize + size(TEpsrc)
             qstarm => c_s%qstarm  ; isize = isize + size(qstarm)
             trip => c_s%trip; isize = isize + size(trip)
             TpSq => c_s%TpSq; isize = isize + size(TpSq)
             Tpdqp => c_s%Tpdqp; isize = isize + size(Tpdqp)
             radtend => c_s%radtend; isize = isize + size(radtend)
             RHxym => c_s%RHxym; isize = isize + size(RHxym)
             RHmsqr => c_s%RHmsqr; isize = isize + size(RHmsqr)
             tempxym => c_s%tempxym; isize = isize + size(tempxym)
             pflux => c_s%pflux; isize = isize + size(pflux)
             pfluxdiff => c_s%pfluxdiff; isize = isize + size(pfluxdiff)
             pmassflux => c_s%pmassflux; isize = isize + size(pmassflux)
             penegflux => c_s%penegflux; isize = isize + size(penegflux)


             xkn     => c_s%xkn
             ykn     => c_s%ykn
             xk      => c_s%xk
             yk      => c_s%yk
             xks     => c_s%xks 
             wavexy  => c_s%wavexy
             ug      => c_s%ug
             vg      => c_s%vg
             wls     => c_s%wls
             uls     => c_s%uls
             udrift  => c_s%udrift
             vdrift  => c_s%vdrift
             stokesw => c_s%stokesw
             stokesa => c_s%stokesa
             stokess => c_s%stokess
             stokes  => c_s%stokes
             dtg     => c_s%dtg
             dslg    => c_s%dslg 
             dzg     => c_s%dzg
             lwc    => c_s%lwc
             phiw    => c_s%phiw 
             phiv    => c_s%phiv
         return
         end subroutine fill_cs
      end module con_stats
