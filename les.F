c     ----------------------------------------------------------------------
      module pars
c     ----------------------------------------------------------------------
      integer :: iti,itmax,imean,ihst,itape
      integer :: itstr,it_his,it_viz,i_viz

      integer :: ismlt,ifree,iradup,
     +     iupwnd,ibuoy,ifilt,itcut,isubs,ibrcl,iocean,method,
     +     idebug,iz_space,ivis0,ifix_dt,new_vis,iDNS,
     +     icouple,iTcouple,iHcouple,ievap,ifields,ilin,ineighbor,
     +     icoalesce,ipart_method,ipartdiff,isfs,iexner,ilongwave
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
      real :: ugcont,vgcont,zo,surf_RH,psurf,surf_temp
      integer :: isfc(nscl)

      !Radiation parameters
      real :: rad_kappa,rad_Fct,rad_Fcb
      real :: radflux(0:maxnz),radsrc(maxnz)

      !Material parameters, both air and particles
      real :: Lv = (25.0 - 0.02274*26.0)*10.0**5
      real :: Sal,Ms,Gam,Ion
      real :: rhow,nuf,rhos
      real :: CpaCpp,Cpa,Cpp,Pra,Sc,Mw,Ru,rhoa,Rd

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
      integer  ::   nvel, npre, nhis1, nprt,nrad,
     +     nhisp, nvelc, 
     +     nviz_z, nviz_y, 
     +     nviz_x, nviz_s,
     +     kfile, jfile, ibcl, ibcu,
     +     igrdr, imach, itn, it_his_nxt, it_viz_nxt
      logical ::    mnout, micut, mtape, mhis, msave,mtrans,
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
     +     w_mn(0:maxnz1), t_mn(0:maxnz1,nscl), alphaC(0:maxnz1,nscl)
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
      character*80  path_res, path_sav, path_his, path_prt,
     +     path_hp, path_sav_hp, path_part, path_ran,
     +     path_v, path_c, path_p, path_h,
     +     path_sav_v, path_sav_c,
     +     path_sav_p, path_sav_h, path_sav_part,
     +     path_histog, path_sav_hist,
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
      module particles
      integer :: rproc,trproc,tproc,tlproc,lproc,blproc,bproc,brproc
      integer :: pr_r,pl_r,pt_r,pb_r,ptr_r,ptl_r,pbl_r,pbr_r
      integer :: pr_s,pl_s,pt_s,pb_s,ptr_s,ptl_s,pbl_s,pbr_s
      real :: ymin,ymax,zmin,zmax,xmax,xmin
      real, allocatable :: uext(:,:,:), vext(:,:,:), wext(:,:,:)
      real, allocatable :: u_t(:,:,:), v_t(:,:,:), w_t(:,:,:)
      real, allocatable :: Text(:,:,:),T_t(:,:,:)
      real, allocatable :: T2ext(:,:,:),T2_t(:,:,:)
      real, allocatable :: partTsrc(:,:,:),partTsrc_t(:,:,:)
      real, allocatable :: partHsrc(:,:,:),partHsrc_t(:,:,:) 
      real, allocatable :: partTEsrc(:,:,:),partTEsrc_t(:,:,:) 
      real, allocatable :: partcount_t(:,:,:),partsrc_t(:,:,:,:)
      real, allocatable :: vpsum_t(:,:,:,:),vpsqrsum_t(:,:,:,:)
      real, allocatable :: upwp_t(:,:,:),upwp(:,:,:)
      real, allocatable :: partcount(:,:,:),partsrc(:,:,:,:)
      real, allocatable :: vpsum(:,:,:,:),vpsqrsum(:,:,:,:)
      real, allocatable :: Tpsum(:,:,:),Tpsum_t(:,:,:)
      real, allocatable :: Tpsqrsum(:,:,:),Tpsqrsum_t(:,:,:)
      real, allocatable :: Tfsum(:,:,:),Tfsum_t(:,:,:)
      real, allocatable :: qfsum(:,:,:),qfsum_t(:,:,:)
      real, allocatable :: wpTpsum(:,:,:),wpTpsum_t(:,:,:)
      real, allocatable :: radsum(:,:,:),radsum_t(:,:,:) 
      real, allocatable :: rad2sum(:,:,:),rad2sum_t(:,:,:) 
      real, allocatable :: multcount(:,:,:),multcount_t(:,:,:) 
      real, allocatable :: mwsum(:,:,:),mwsum_t(:,:,:)
      real, allocatable :: qstarsum(:,:,:),qstarsum_t(:,:,:) 

c     --- SFS velocity calculation ---------
      real, allocatable :: sigm_s(:,:,:),sigm_sdx(:,:,:),sigm_sdy(:,:,:)
      real, allocatable :: sigm_sdz(:,:,:),sigm_sext(:,:,:)
      real, allocatable :: sigm_sdxext(:,:,:),sigm_sdyext(:,:,:)
      real, allocatable :: sigm_sdzext(:,:,:)
      real, allocatable :: vis_ss(:,:,:),vis_sext(:,:,:) 

      integer :: particletype,pad_diff
      integer :: numpart,tnumpart,ngidx
      integer :: iseed
      integer :: num100=0, num1000=0, numimpos=0
      integer :: tnum100, tnum1000, tnumimpos
      integer :: denum, actnum, tdenum, tactnum
      integer :: num_destroy=0,tnum_destroy=0

      real :: Rep_avg,part_grav(3)
      real :: radavg,radmin,radmax,tempmin,tempmax,qmin,qmax
      real :: vp_init(3),Tp_init,radius_init
      real :: pdf_factor
      integer*8 :: mult_init,mult_factor,mult_a,mult_c

      real :: avgres=0,tavgres=0

      integer, parameter :: histbins = 512
      real :: hist_rad(histbins+2)
      real :: bins_rad(histbins+2)      

      
      !REMEMBER: IF ADDING ANYTHING, MUST UPDATE MPI DATATYPE!
      type :: particle
      integer :: pidx,procidx,nbr_pidx,nbr_procidx
      real :: vp(3),xp(3),uf(3),xrhs(3),vrhs(3),Tp,Tprhs_s
      real :: Tprhs_L,Tf,radius,radrhs,qinf,qstar,dist
      real :: res,m_s,Os
      real :: u_sub(3),sigm_s
      integer*8 :: mult
      type(particle), pointer :: prev,next
      end type particle

      type(particle), pointer :: part,first_particle
      end module particles
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
     +     batagk, batag,
     +     cdbtm, dtdzf(nscl), dtjump, ugal, divgls,
     +     z1, utausv, dx, dy, dz, dt, 
     +     fnxy, dzdz, dsl, c23, dtgama, dtzeta, xkmax,
     +     time, zody, zody74,
     +     utau, wtsfc(nscl),
     +     uwsfc, vwsfc, amonin,
     +     zol, hol, smal_e, sml_eg
      end type con_d
      type(con_d), target :: c_c
      real, pointer ::
     +     vk, vkin, vk74, vk74in, 
     +     gcp, fcor_h, pi2,
     +     batagk, batag,
     +     cdbtm, dtdzf(:), dtjump, ugal, divgls,
     +     z1, utausv, dx, dy, dz, dt,
     +     fnxy, dzdz, dsl, c23, dtgama, dtzeta, xkmax, 
     +     time, zody, zody74,
     +     utau, wtsfc(:), 
     +     uwsfc, vwsfc, amonin,
     +     zol, hol, smal_e, sml_eg
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
      utau   => c_c%utau
      wtsfc  => c_c%wtsfc
      uwsfc  => c_c%uwsfc
      vwsfc  => c_c%vwsfc
      amonin => c_c%amonin
      zol    => c_c%zol
      hol    => c_c%hol
      smal_e => c_c%smal_e
      sml_eg => c_c%sml_eg
      return
      end subroutine fill_cc
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
     +           radtend(maxnz)
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
     +           radtend(:)
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
c ----------------------------------------------------------------------
      program les_mpi_2d
c
      use pars
      use fields
      use particles
      use con_data
      use con_stats
      include 'mpif.h'
c
c ------------- definition of internal flags
c
c
c       iDNS    = 0; call the subgrid computation of vis_m and vis_s
c               = 1; call the molecular viscosity and diffusivity
c
c       igrdr   =  3; data comes from restart file
c               =  2; data comes from initialization (random)
c               =  1; data comes from coarser grid (or otherwise)
c
c       ibcu    = 1 ; upper boundary condition set by radiation bc
c               = 0 ; fixed value = 0.
c               = -1; value defined by coarser mesh for all variables
c
c       ibcl    = 0 ; lower boundary condition set by similarity theory (sr. setup)
c               = -1; value defined by coarser mesh for all variables
c
c       ifix_dt = 0 ; variable time step with fixed cfl number in setcon
c               = 1 ; fixed time step set in sr. get_dt
c
c       ifree   = 0 ; use spatially averaged surface conditions for MO (call lower)
c               = 1 ; use point-by-point conditions for MO free convection (call lower_free)
c
c       ihst    = nn; frequency at which global variables are output in history file
c               < 0 ; no history files
c
c       it_his  = time step where history files start, incremented by itape
c
c       it_viz  = time step where viz files start, incremented by itape
c
c       ismlt   = 1 ; use businger formulas in MO 
c                 0 ; use large and everyone elses formulas in MO 
c
c       iupwnd  = 0;  use skew symmetric formulas for all derivatives
c                     in scalar equations
c               = 1;  use hybrid upwind scheme for all derivatives
c                     in scalar equations
c
c       ivis0   = 0; old eddy viscosity model 
c               = 1; new eddy viscosity model 
c
c       new_vis = step; the iteration step for which the new model
c                       is turned on when ivis0=1
c               < 0; new model is on at all steps for ivis0=1
c
c       nscl  .ge. 1   number of scalars to be followed set in parameter statements
c                      change entries in sr. init, and sr. suft for surface bc's
c
c -------------------------------------------------------------------------------
c
c ---------- initialize MPI, get myid, numprocs, 
c            test if on root process
c
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
      call mpi_comm_size(mpi_comm_world,numprocs,ierr)
c
      i_root = 0
      l_root = .false.
      if(myid .eq. i_root) l_root = .true.
c
      l_debug = .false.
      if(idebug .eq. 1) l_debug = .true.
c
      ts_mpi = mpi_wtime()

!----- Read the input file for all necessary parameters
      call read_input_file

c
c -------- set number of x-y slab cpus
c
!      ncpu_s = 8
c
      case_inp = 'cou'
c
      call get_units
      call gridd
      call setcon
c
c -------------- scratch run
c
      if (iti.eq.0)  then
         igrdr = 2
         case = case_inp
         call init
         call setup(it)

         call particle_setup
         call particle_init 
         !call read_part_res !comment out usually 

c
c ---------- choose routine for getting initial guess
c
         if(iocean .eq. 1) then
            call randoc
         else
            if (ifields .eq. 1) then
            call get_fields
            else
            call random   
            end if
         endif
         !Call dns_vis even when doing LES since vis_m needs to be initialized with something
         call dns_vis
         call get_max
      else

         igrdr = 3
         call restart
         call setup(it)

         !Call dns_vis even when doing LES since vis_m needs to be initialized with something
         call dns_vis
         call get_max

 
         call particle_setup
         call read_part_res
      endif


c
c --------------- time loop ------------
c
      tzero = time
      call get_dt
 9000 continue
      call set_sav(it,iti)
      if (myid==0) then
      write(*,*) 'Starting time loop'
      write(*,*) 'it,time = ',it,time
      end if
      if (it == 1) numpart = 0
      part => first_particle
      do while (associated(part))
      if (part%pidx == 1 .and. part%procidx == 0) then
      write(*,'(a7,4e15.6)') 'xp1:  ',time,part%xp(1:3)
      write(*,'(a7,4e15.6)') 'vp1:  ',time,part%vp(1:3)
      write(*,'(a7,4e15.6)') 'uf1:  ',time,part%uf(1:3)
      write(*,'(a7,2e15.6)') 'Tp1:  ',time,part%Tp
      write(*,'(a7,2e15.6)') 'Tf1:  ',time,part%Tf
      write(*,'(a7,2e15.6)') 'rad1: ',time,part%radius 
      write(*,'(a7,2e15.6)') 'qinf1:',time,part%qinf 
      write(*,'(a7,2e15.6)') 'qstr1:',time,part%qstar 
      write(*,'(a7,e15.6,i)') 'mult1:',time,part%mult

      end if

      part => part%next

      if (it == 1) numpart = numpart + 1
      end do

      if (myid==0) then
       write(*,*) 'time,tnumpart:',time,tnumpart
       write(*,*) 'time,tnum_destroy:',time,tnum_destroy
       write(*,'(a15,3e15.6)') 'radmin,radmax:',time,radmin,radmax
       write(*,'(a15,3e15.6)') 'tempmin,tempmax:',time,tempmin,tempmax
       write(*,'(a15,3e15.6)') 'qmin,qmax:',time,qmin,qmax
       write(*,'(a15,3e15.6)') 'time,radavg:',time,radavg
       write(*,'(a15,e15.6,3i)'),'time,100,1000,impos:',time,tnum100,
     + tnum1000,tnumimpos
      !if (myid==0) write(*,*), 'time,actnum',time,tactnum
      !if (myid==0) write(*,*), 'time,denum',time,tdenum
      end if


      if(it .ge. new_vis .and. ivis0 .eq. 1) then
          ivis = 1
      else
          ivis = 0
      endif


c
c ---------------- 3 stage runge-kutta time stepping
c
      t_stage_s = mpi_wtime()
      do  8999 istage=1,3
c
      dtzeta = dt*zetas(istage)
      dtgama = dt*gama(istage)
c
c ---------- compute derivatives of (u,v,w)
c

      call exchange
c

      call get_derv

c
c --------- new eddy viscosity, and bcs
c
      if(iss .eq. 0 .and. ifree .eq. 0) then
         if (iDNS .eq. 1) then
            call lower_dns(it)
         else
            call lower(it)
         end if
      elseif(ifree .eq. 1) then
         if (iDNS .eq. 1) then
            call lower_dns(it)
         else
            call lower_free(it)
         end if
      endif
      if(ise .eq. numprocs-1) then
         if (iDNS .eq. 1) then
            call upper_dns
         else
            call upper
         end if
      endif
      call bcast_pbc
      call get_means(istage)
      if(ivis .eq. 1) then
         call iso(it)
         call surfvis(it)
      endif
      if(istage .eq. 1)then
        call xy_stats
        call tke_budget
        do iscl=1,nscl
        call Tvar_budget(iscl)
        enddo
        call extra_flux_terms
        call pbltop(itop)
      endif


c
c ------------ save velocity field
c
      if(msave .and. istage .eq. 1) then
         call save_v(it)
         call save_particles
      endif
      if(msave_v .and. istage .eq. 1) then
         call save_viz(it)
      endif
c
c ------------ save pressure field
c
      if(msave .and. istage .eq. 1) then
         call save_p
      endif


      if(mhis  .and. istage .eq. 1) then
          call write_histograms
      end if

c
c --------- get rhs for all equations
c
      !t_s = mpi_wtime()
      call comp1(istage,it)
      !t_f = mpi_wtime()
      !write(*,*) 'proc, comp1: ',myid,t_f-t_s
      !call mpi_barrier(mpi_comm_world,ierr)
      !if (myid==5) write(*,*) 'time comp1: ',t_f-t_s
      if(istage .eq. 1) then
         if(msave .and. l_root) call save_c(it)
      endif

     
c
c --------- solve for pressure
c
      !t_s = mpi_wtime()
      call comp_p
      !t_f = mpi_wtime()
      !write(*,*) 'proc, comp_p: ',myid,t_f-t_s
      !call mpi_barrier(mpi_comm_world,ierr)
      !if (myid==5) write(*,*) 'time comp_p: ',t_f-t_s
c
c --------- add pressure gradient and dealias
c
      !t_s = mpi_wtime()
      call comp2
      !t_f = mpi_wtime()
      !!write(*,*) 'proc, comp2: ',myid,t_f-t_s
      !call mpi_barrier(mpi_comm_world,ierr)
      !if (myid==5) write(*,*) 'time comp2: ',t_f-t_s
      if(micut) then
         call dealias
      endif




c      
c -------- update particles
c
         if (ipart_method .eq. 1) then
         t_s = mpi_wtime()
         !call particle_update_cn(it,istage) 
         call particle_update_rk3(it,istage) 

         t_f = mpi_wtime()
         call mpi_barrier(mpi_comm_world,ierr)
         if (myid==5) write(*,*) 'time part: ',t_f-t_s
         end if


      if(mnout .and. istage .eq. 1)  then
          if(l_debug) then
             call print(nprt,it,izs,ize)
          endif
          if(l_root) call print(6,it,1,nnz)
      endif
      if(l_root) then
         if(mhis  .and. istage .eq. 1)  call write_his(itop)
         if(mhis  .and. istage .eq. 1 .and. mtape) call close_his
      endif
c
 8999 continue
      t_stage_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time stage: ',t_stage_f - t_stage_s

      if (ipart_method .eq. 2) then
         t_s = mpi_wtime()

         call particle_update_BE(it)

         t_f = mpi_wtime()
         call mpi_barrier(mpi_comm_world,ierr)
         if (myid==5) write(*,*) 'time part: ',t_f-t_s
      end if

      !Call coalescence outside of RK loop since it's not appropriate as
      !part of RK scheme
      if (icoalesce) then
      t_s = mpi_wtime()

      call particle_coalesce

      call mpi_barrier(mpi_comm_world,ierr)
      t_f = mpi_wtime()
      if (myid==5) write(*,*) 'time coalescence:', t_f - t_s
      end if

      if (ipartdiff) then
      t_s = mpi_wtime()

        if (isfs == 2) then
           !Call Weil et al. (2004) Lagrangian stochastic model
           call SFS_velocity
        elseif(isfs == 1) then
           !Call stochastic model for the position based on vis_s
           call SFS_position
        end if

      call mpi_barrier(mpi_comm_world,ierr)
      t_f = mpi_wtime()
      if (myid==5) write(*,*) 'time partdiff:', t_f - t_s
      end if

      call radius_histogram

      call get_max
      call get_dt
      if (it.ge.itmax) go to 99000
      go to 9000
c
99000 continue
      te_mpi = mpi_wtime()
      write(6,9997) (te_mpi - ts_mpi)
 9997 format(' Job Execution Time = ',e15.6)
c
 9998 continue
      call mpi_finalize(ierr)
c
      stop
      end
      subroutine get_max
c
c --------- routine computes max velocities as sweep through
c           the velocity field 
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      real u_send(5), u_recv(5)
c

      dx_i = 1.0/dx
      dy_i = 1.0/dy
c
      u_temp = 0.0
      v_temp = 0.0
      w_temp = 0.0
      vis_temp = 0.0
      vis_temp = 0.0
      do iz=izs,ize
c
        u_xy = 0.0
        v_xy = 0.0
        w_xy = 0.0
        vis_xym = 0.0
        vis_xys = 0.0
        vis_xy = 0.0
        do iy=iys,iye
        do ix=1,nnx
           u_xy = amax1(u_xy,abs(u(ix,iy,iz)+stokes(iz)))
           v_xy = amax1(v_xy,abs(v(ix,iy,iz)))
           w_xy = amax1(w_xy,abs(w(ix,iy,iz)))
           vis_xym = amax1(vis_xym,vis_m(ix,iy,iz))
           vis_xys = amax1(vis_xys,vis_s(ix,iy,1,iz))
           vis_xy = amax1(vis_xys,vis_xym)
        enddo
        enddo
        u_xy   = u_xy*dx_i
        v_xy   = v_xy*dy_i
        wsav   = w_xy
        w_xy   = w_xy/abs(dzw(iz))
        vis_xy = vis_xy/amin1(dx,dy,dzw(iz))**2
c
        u_temp = amax1(u_xy,u_temp)
        v_temp = amax1(v_xy,v_temp)
        w_temp = amax1(w_xy,w_temp)
        vis_temp = amax1(vis_xy,vis_temp)
c
c       if(iz .le. 15) then
c         write(6,6000) iz, wmax
c6000     format(' in get_dt iz = ',i3,' wmax = ',e15.6)
c       endif
c
      enddo
      u_send(1) = u_temp
      u_send(2) = v_temp
      u_send(3) = w_temp
      u_send(4) = wsav
      u_send(5) = vis_temp 
  
c
      call mpi_allreduce(u_send,u_recv,5,mpi_real8,
     +     mpi_max,mpi_comm_world,ierror)
c
      umax = u_recv(1)
      vmax = u_recv(2)
      wmax = u_recv(3)
      wabs = u_recv(4)
      vismax = u_recv(5)


c
      return
      end
      subroutine get_dt
c
c ---------- routine computes max time step for given cfl number
c            from max's found previously
c
      use pars
      use con_data
      use con_stats
c

      ucfl = umax
      vcfl = vmax
      wcfl = wmax
      ucflm = ucfl
      vcflm = vcfl
      wcflm = wcfl
      vel_max = wcflm
      vel_max = amax1(ucflm,vel_max)
      vel_max = amax1(vcflm,vel_max)

      if(vel_max .le. 0.0) then
          write(6,6000) ucflm, vcflm,wcflm, vel_max
 6000     format('6000, sr. get_dt bad news, umax = ',e15.6,/,
     +           ' vmax = ',e15.6,' wmax = ',e15.6,/,
     +           ' vel_max = ',e15.5,/,
     +           ' infinite time step !!!')
          stop
      endif

c
c ---------------- choose fixed or variable time step
c
      if(ifix_dt .ne. 0) then
c
c ------------- if used, change to fit your problem
c
!        dt_new = 0.5  !Now in input file
      else
c
c ------------------- new estimate of best time step
c                     from cfl constraint
c
      dt_new  = cfl/vel_max
      dt_new = amin1(dt_new, 5.0)
c     dt_new = amin1(dt_new, 10.0)
      endif

c ---------------- compare against viscous stability limit
c
      if (ifix_dt .eq. 0) then
      if(vismax*dt_new .gt. 0.5) then
         dt_cfl = dt_new
         dt_new = 0.5/vismax
         if(l_root) then
            write(6,6200) dt_new, dt_cfl, vismax
 6200       format(' 6200 get_dt: cfl time step too large',/,
     +      '   viscous time step = ',e15.6,
     +      ' cfl time step = ',e15.6,' vismax = ',e15.6)
         endif
      endif
      end if
c
c -------- for safety if restart set timestep = saved timestep in
c          restart file
c
      if(it_step .eq. iti .and. iti .ne. 0) then
        dt_new = dt
      endif

c
      return
      end
      subroutine lterp(n,zary,zpt,i,ip1,ratio)
c
c ---- linear interpolation for zpt in zary, where zary is 
c      monotonic increasing or decreasing function
c
      dimension zary(*)
      nm1 = n-1
      if(n.le.1) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
      endif
      if(zary(1) .lt. zary(2)) go to 1
                               go to 101
    1 continue
c
c **** monotonic increasing array
c
        if(zpt .lt. zary(1)) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
        else if(zpt .gt. zary(n)) then
          i = n
          ip1 = n
          ratio = 1.0
          go to 999
        endif
        do j=1,nm1
              if(zpt .ge. zary(j) .and.
     $           zpt .le. zary(j+1)) then
                 i = j
                 ip1 = j+1
                 ratio = (zpt - zary(i))/(zary(ip1) - zary(i))
                 go to 999
              endif
        enddo
c
c **** decreasing array
c
  101 continue 
        if(zpt .gt. zary(1)) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
        else if(zpt .lt. zary(n)) then
          i = n
          ip1 = n
          ratio = 1.0
          go to 999
        endif
        do j=1,nm1
              if(zpt .le. zary(j) .and.
     $           zpt .ge. zary(j+1)) then
                 i = j
                 ip1 = j+1
                 ratio = (zpt - zary(i))/(zary(ip1) - zary(I))
                 go to 999
              endif
        enddo
  999 continue
      return
      end
      subroutine setcon
c
      use pars
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      external get_zi
c
c ----------------- get machine type, can erad datadir also
c
      open(99,file='mach.file')
      read(99,9000) imach
 9000 format(i1)
c
      close(99)
c
c ----------- initialize fft
c
      call rffti(nnx,trigx(1,1))
      call rffti(nny,trigx(1,2))
      call cffti(nny,trigc(1))
c
c ----------- start step for history files
c
      it_his_nxt = it_his
      it_viz_nxt = it_viz
c
c ---------------- set min value of e
c
      if(iocean .eq. 1) then
c
         smal_e = 0.0
         smal_e = 1.0e-12
c        smal_e = 6.0e-03
      else
         smal_e = 1.0e-09
c        smal_e = 0.0
      endif
c
c ---------------------- set constants in eddy viscosity model
c
      ck       = 0.1
      ceps     = 0.93
      csmag    = sqrt(ck*sqrt(ck/ceps))
      stab_c   = 0.76
c
c ----------------- set stability constant
c
      stabmin = 1.0e-12
c
c ---------------- minimum dsl length constant
c 
      almin_c = 0.0001
c
c -------------- initialize grid restart flag
c
      igrdr = 1
c
c -------------- create mpi operation to find max and location
c                using local gradient method
c
      call mpi_op_create(get_zi,.true.,ziloc,ierror)
c
c ------------------- define coefficients for 3-order runge-kutta
c                     time stepping scheme, borrowed from Spalart,
c                     Moser and Rogers, J. Comp. Physics 3/21/90
c                     Note this is a simplier version since all terms
c                     are lumped in the non-linear terms.
c                     cfl number is for an entire runge-kutta step
c                     in this case three stages. cfl = max(u)*dt/dx
c
      zetas(1) = 0.0
      zetas(2) = -17.0/60.0
      zetas(3) = -5.0/12.0
      gama(1)  = 8.0/15.0
      gama(2)  = 5.0/12.0
      gama(3)  = 3.0/4.0
      etas(1)  = -1.0
      etas(2)  = -1.0 + 8.0/15.0
      etas(3)  = -1.0 + 2.0/3.0
c
c ----------- a full step, at the new time
c
      etas(4)  =  0.0
c
c
      return
c
      end

      subroutine setup(it)
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      it = iti
c
c ------------ turn on new sgs model at a particular step
c
      if(it .ge. new_vis .and. ivis0 .eq. 1) then
          ivis = 1
      else
          ivis = 0
      endif
c
      if(igrdr . eq. 3) then
         if(l_root) then
            write(6,6)iti,utau,tsfcc(1) ,qstar(1)
            write(6,510)
            write(6,520) wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +            ,cdbtm,ugcont
            call print(6,it,1,nnz)
         endif
          if(l_debug) then
            write(nprt,6)iti,utau,tsfcc(1) ,qstar(1)
            write(nprt,510)
            write(nprt,520) wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +            ,cdbtm,ugcont
            call print(nprt,it,izs,ize)
          endif
      endif
c     if(ifilt.eq.1)call filter
      if(l_root) then
         write(6,1) nnx,nny,nnz,ismlt,ifilt,iti,itmax,
     +             iupwnd,ibuoy,itcut,
     +             dt,zo,tsfcc(1),isubs,ibrcl,
     +             method, iocean, ivis
      endif
      if(l_debug) then
         write(nprt,1) nnx,nny,nnz,ismlt,ifilt,iti,itmax,
     +             iupwnd,ibuoy,itcut,
     +             dt,zo,tsfcc(1),isubs,ibrcl,
     +             method, iocean, ivis
      endif
c
c -------------- boundary condition flags 
c
      ibcu = iradup
c     ibcu = 0
      ibcl = 0
c
c -------------------- wavenumbers, introduce a normalized
c                      set of wavenumbers to eliminate computation
c                      in derivatives , xderiv, yderiv
c
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
c
c ----------- choose correct sign so gravity waves
c             propagate out of the domain
c
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
c
c -------------------- set length scale for SGS model
c
      if(iz_space .eq. 0) then
c
c ------------- uniform vertical spacing
c
      dx32 = dx*3./2.
      dy32 = dy*3./2.
      dsl  = (abs(dx32*dy32*dzw(1)))**(1./3.)
      dslg = dsl
      if(l_root)  write(6,2000) dsl
      if(l_debug) write(nprt,2000) dsl
c
c --------------------- create dsl array for easy indexing in comp1
c
      do iz=0,nnzp1
         dsl_z(iz) = dslg
      enddo
c
c ------------- variable vertical spacing
c
      else
c
c ----------- just estimate dsl for average spacing
c
         dx32 = dx*3./2.
         dy32 = dy*3./2.
c
         dsl_max = (abs(dx32*dy32*dzw(0)))**(1./3.)
         do iz=0,nnzp1
            dsl_z(iz) = (abs(dx32*dy32*dzw(iz)))**(1./3.)
            if(dsl_z(iz) .gt. dsl_max) dsl_max = dsl_z(iz)
         enddo
c        do iz=0,nnzp1
c           dsl_z(iz) = dsl_max
c        enddo
         dsl  = dsl_max
         dslg = dsl
      endif
c
      gridr = 1.0
      sml_eg = smal_e*gridr
c -------------------- set viscosity model parameters 
      if(ivis .ne. 1) then
        viscon = 0.0
        xksurf = 0.0
        nmatch = -1
        myid_newvis = 0
        do iz=1,nnz
           dfac(iz) = 1.0
        enddo
      endif
c ------------------- set stokes velocity for atmos/oceanic flow
      call stokesv
c
c --------- can add a time factor so as to skip into any part of
c           the specified geostrophic arrays. time factor in seconds
c
      t_factor = 7200.0
c
c ---------- for print out to get more digits
c
      t_ref = 300.0
c
c -------------------- specify cooling rate and initial
c                      temperature even for restarts
c
      c_rate   = 0.25/3600.0
      t_surf_i = 265.0
c
c -------------------- do not look for zi below zi_min
c
      zi_min = 30.0
      if(iocean .eq. 1) zi_min = -5.0
      iz_min = 1
      do iz=1,nnz-1
         if(zz(iz) .lt. zi_min .and.
     +      zz(iz+1) .ge. zi_min) iz_min = iz
      enddo
      if(l_root) then
         write(6,9000) zi_min, iz_min
      endif
c
 9998 continue
      return
c --------------------------- format statements
    6 format(///,' DATA FROM RESTART FILE AT STEP =',I5,
     +       ' U_* = ',e15.6,' TS = ',e15.6,' Q_* = ',e15.6,///)
  510 format(' RESTART ***** CASE WITH : ******',/)
  520 format(' WT = ',e12.4,', U_* = ',e12.4,', L = ',e12.4,
     +       ', DTDZ FREE = ',e12.4,', ZODY = ',e12.4,/,10x,
     +       '  ZO(BTM) = ',e12.4,', CDBTM = ',e12.4,
     +       ', UG = ',e12.4)
    1 format(10x,' NNX = ',i5,',  NNY = ',i5,
     + ',  NNZ = ',i5,/,10x,' SFC SMLT = ',i1,
     + ',  FILTER = ',i1,
     + ',  ITI = ',i6,',  ITMAX = ',i6,/,10x,
     + ' IUPWIND = ',i1,',  BUYNCY = ',i1,
     + ',  NO ALISNG = ',i1,',  ITCUT = ',i1,/,10x,
     + ' DT = ',e15.6,',  ZO = ',e15.6,',  TS = ',e15.6,
     + ',  SUBSD = ',i1,/,
     + 10x,' BRCLICITY = ',i1,',  METHOD = ',i1,',  IOCEAN = ',i1,
     + ',  IVIS = ',i1)
 2000 format(10x,' DSL = ',e15.6)
 9000 format(' Search for zi above the height = ',e15.6,/,
     +       ' iz_min = ',i5)
      end
      subroutine nblnk(word)
      parameter (nmax=304)
      character wordt*304, word*(*)
      nchar = len(word)
      if(nchar .gt. nmax) then
         write(6,6000) nchar,nmax
 6000    format(' TROUBLE, IN SR. NBLNK : NCHAR = ',i6,
     +          ' EXCEEDS NMAX = ',i6)
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
c
      return
      end
      subroutine blnk(word)
      character word*(*)
      nchar = len(word)
      do j=1,nchar
         word(j:j) = ' '
      enddo
c
      return
      end
      subroutine iso(it)
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      real sfk(1:nnz)
c
c ---- get isotropy factor and scale it to match at the matching
c      height. uses boundary conditions from lower and upper. 
c
      do iz=1,nnz
c        dfac(iz) = 1.0
         dfac(iz) = 0.0
         sfk(iz)  = 0.0
      enddo
      do iz=izs,ize
         dfac(iz) = 1.0
      enddo
c
c ------ set nmatch equal to fraction of initial zi in sr. random
c
c     nmatch = izi/2
c     nmatch = 16
      nmatch = 48
      do i=0,numprocs-1,ncpu_s
         if(nmatch .ge. iz_s(i) .and.
     +      nmatch .le. iz_e(i)) myid_newvis = i
      enddo
c
      do iz=izs,min(ize,nmatch)
         izp1 = iz + 1
         izm1 = iz - 1
         weit = dzw(iz)/(dzw(iz) + dzw(izp1))
         weit1 = 1.0 - weit
c
c ---- get fluctuating strains
c
         do j=iys,iye
         do i=1,nnx
            s11 = weit1*ux(i,j,iz)**2 + weit*ux(i,j,izp1)**2
            s22 = weit1*vy(i,j,iz)**2 + weit*vy(i,j,izp1)**2
            wz  = (w(i,j,iz)-w(i,j,izm1))*dzw_i(iz)
            wzp = (w(i,j,izp1)-w(i,j,iz))*dzw_i(izp1)
            s33 = weit*wzp**2 + weit1*wz**2
            s12 = weit1*(uy(i,j,iz) + vx(i,j,iz))**2 +
     +            weit*(uy(i,j,izp1) + vx(i,j,izp1))**2
            s13 = (((u(i,j,izp1) - u(i,j,iz) +
     +            u_mn(iz) - u_mn(izp1))*dzu_i(izp1) +
     +            wx(i,j,iz)))**2
            s23 = (((v(i,j,izp1) - v(i,j,iz) +
     +          v_mn(iz) - v_mn(izp1))*dzu_i(izp1) +
     +          wy(i,j,iz)))**2
            sfk(iz) = sfk(iz) + 2.0*(s11 + s22 + s33) +
     +                       s12 + s13 + s23
         enddo
         enddo
         sfk(iz) = sfk(iz)*fnxy
      enddo
      call mpi_sum_z(sfk,i_root,myid,nnz,1)
c
      do iz=izs,min(ize,nmatch)
         izp1 = iz + 1
         izm1 = iz - 1
c
         sfk(iz) = sqrt(sfk(iz))
         smk = sqrt((u_mn(izp1)-u_mn(iz))**2 +
     +              (v_mn(izp1)-v_mn(iz))**2)*abs(dzu_i(izp1))
         if(sfk(iz) .le. 0. .and. smk .le. 0.) then
           dfac(iz) = 1.0
         else
           dfac(iz) = sfk(iz)/(sfk(iz) + smk)
         endif
c     if(l_root) write(6,6001) iz,sfk(iz),smk,dfac(iz)
 6001 format(' iz = ',i3,' sfk = ',e15.6,
     +       ' smk = ',e15.6,' dfac = ',e15.6)
      enddo
c
c
c ---- rescale ratio to give unity at match height
c      and if nested grid match value at upper boundary
c      of coarser grid
c
      if(myid .eq. myid_newvis) then
         dfacm = dfac(nmatch)
      endif
c
      call mpi_bcast(dfacm,1,mpi_real8,
     +              myid_newvis,mpi_comm_world,ierr)
c
      do iz=izs,min(ize,nmatch)
         dfac(iz) = dfac(iz)/dfacm
         dfac(iz) = amax1(dfac(iz), 0.1)
         dfac(iz) = amin1(dfac(iz), 1.0)
      enddo
c
c --------- gather dfac on all processes for printing and use in tke_vis
c           use reduce and divide by number of slab cpus
c
      call mpi_sum_z(dfac,i_root,myid,nnz,1)
      fncpu_s = 1.0/float(ncpu_s)
      do iz=1,nnz
         dfac(iz) = dfac(iz)*fncpu_s
      enddo
c
c     call mpi_gath_root(dfac(izs),dfac(1),iz_s,iz_e,izs,ize,nnz,
c    +                   myid,numprocs,ncpu_s)
c
c     if(l_root) write(6,6000) nmatch,ivis,(iz,dfac(iz),iz=1,nnz)
 6000 format(' in sr. iso, nmatch = ',i3,/,
     +       ' ivis = ',i3,'iz',5x,'dfac',/,(i3,1x,e15.6))
c     write(nprt,3001) (iz,dfac(iz),iz=1,nnz)
 3001 format(' iz ',5x,' dfac ',/,(i5,e15.6))
      return
      end
      subroutine surfvis(it)
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      real xkvis(nnx,iys:iye), alwk(nnx,iys:iye)
c
      real send(3), buf(3)
c
      xksurf = 0.0
      viscon = 0.0
      vise   = 0.0
c
c ----------- only root process(es) compute 
c
      if(iss .eq. 0) then

c     ck = 0.1
c     csmag = 0.18
c     xkmax  = dzdz/dt/5.
      iz   = 1
      izm1 = iz - 1
      izp1 = iz + 1
c     xkmax  = dzu(izp1)*dzu(izp1)/(5.0*dt)
      dz_i = dzu_i(izp1)
      if(iocean .eq. 1) then
         call sufto(it)
      else
         call suft(it)
      endif
      if(qstar(1) .eq. 0.) then
         zeta = 0.0
      else
         zeta = abs(z(1))/amonin
      endif
      if(ismlt .eq. 1) then
          call busngr(zeta,phim,phis,psim,psis)
      else
          call fzol(zeta,phim,phis,psim,psis)
      endif
      viscon = vk*abs(z(1))/(utau*phim)
      vise   = utau*vk*abs(z(1))/phim
c
c ---- get special value at z1 to match with surface layer
c
      uws = 0.0
      vws = 0.0
      do iy=iys,iye
      do ix=1,nnx
         uws = uws + 0.5*(u(ix,iy,iz)-u_mn(iz) + 
     +         u(ix,iy,izp1) - u_mn(izp1))*(w(ix,iy,iz)-w_mn(iz))
         vws = vws + 0.5*(v(ix,iy,iz)-v_mn(iz) + 
     +         v(ix,iy,izp1) - v_mn(izp1))*(w(ix,iy,iz)-w_mn(iz))
      enddo
      enddo
      uws = uws*fnxy
      vws = vws*fnxy !Brian corrected from fxy 1/15/15
c
c ---- get average fluctuating eddy viscsoity
c
      do iy=iys,iye
      do ix=1,nnx
         e(ix,iy,iz)=amax1(e(ix,iy,iz),sml_eg)
      enddo
      enddo
      dslk = amin1(dsl,vk*abs(z(iz))/csmag)
c     stabmin = 1.e-12
c     almin = 0.0001*dsl
      almin = almin_c*dsl_z(iz)
      do iy=iys,iye
      do ix=1,nnx
         alwk(ix,iy)=dslk
c
c --------no stability corrected length scales when
c         new eddy viscosity is on
c
c         stab=batag*(t(ix,iy,1,izp1)-t(ix,iy,1,iz))*dz_i
c         if(stab.gt.stabmin) then
c           als = stab_c*sqrt(e(ix,iy,iz)/stab)
c           alwk(ix,iy) = amin1(dslk,als)
c         endif
c         alwk(ix,iy)=amax1(almin,alwk(ix,iy))
         xkvis(ix,iy)=ck*alwk(ix,iy)*sqrt(e(ix,iy,iz))*dfac(1)
c        xkvis(ix,iy)=amin1(xkvis(ix,iy),xkmax)
      enddo
      enddo
c
c ---- get average viscosity
c
      xkavg = 0.0
      do iy=iys,iye
      do ix=1,nnx
         xkavg = xkavg + xkvis(ix,iy)
      enddo
      enddo
      xkavg = xkavg*fnxy
c
      buf(1) = uws
      buf(2) = vws
      buf(3) = xkavg
      call mpi_sum_xy(buf,myid,iss,ise,3)
c
      uws   = buf(1)
      vws   = buf(2)
      xkavg = buf(3)
c
      xkz1 = vise - sqrt(uws**2 + vws**2)*viscon
      xksurf =  xkz1 - xkavg
      xksurf = amax1(xksurf,0.0)
      xksurf = amin1(xksurf,vise)
c     if(l_root) write(6,6000) dfac(1), xkavg, xkz1, vise, xksurf
 6000 format(' dfac = ',e12.4,' xkavg = ',e12.4,' xkz1 = ',e12.4,/,
     +       ' vise = ',e12.4,' xksurf = ',e12.4)
c
      endif
c
c ---------- broadcast values to other processes
c
      send(1) = xksurf
      send(2) = viscon
      send(3) = vise
c
      call mpi_bcast(send,3,mpi_real8,
     +              i_root,mpi_comm_world,ierr)
c
      xksurf = send(1)
      viscon = send(2)
      vise   = send(3)
c
      return
      end
      subroutine comp1(istep,it)
c
c ----- 3-order runge-kutta time stepping and monotone scalar fluxes in x,y,z.
c       designed to use mpi in x & y directions.
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      parameter(js = 6, ns = 3, nstat = js + ns*nscl)
      real stat(0:nnz,nstat)
c
c ------ temp arrays to hold rhs from step n-1 and
c        field variables from step n 
c
      real urhs(nnx,iys:iye,izs:ize), 
     +     vrhs(nnx,iys:iye,izs:ize),
     +     wrhs(nnx,iys:iye,izs:ize),
     +     erhs(nnx,iys:iye,izs:ize),
     +     trhs(nnx,iys:iye,nscl,izs:ize)
c

      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            urhs(ix,iy,iz) = u(ix,iy,iz) + dtzeta*r1(ix,iy,iz)
            vrhs(ix,iy,iz) = v(ix,iy,iz) + dtzeta*r2(ix,iy,iz)
            wrhs(ix,iy,iz) = w(ix,iy,iz) + dtzeta*r3(ix,iy,iz)
            erhs(ix,iy,iz) = e(ix,iy,iz) + dtzeta*r5(ix,iy,iz)
         enddo
         enddo
      enddo
      do iz=izs,ize
         do l=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            trhs(ix,iy,l,iz) = t(ix,iy,l,iz) + dtzeta*r4(ix,iy,l,iz)
         enddo
         enddo
         enddo
      enddo
c
c --------- get viscosity and rhs of (e,u,v,w)-equations
c           at next step
c
      if (iDNS .eq. 1) then
         call dns_vis
         call rhs_uvw_DNS(istep)
      else
         call tke_vis(istep)   ! WARNING: this has not been checked for iscl=2
         call rhs_uvw(istep)   ! WARNING: this has not been checked for iscl=2
      end if

c
c -------- evaluate rhs of scalar equations
c
      do l=1,nscl
         if (iDNS .eq. 1) then
            call rhs_scl_dns(istep,l,it) 
         else
            call rhs_scl(istep,l,it) 
         end if
      enddo

c
c ---------- gather stat sums on root processor
c            using mpi_reduction over all processors
c
      if(istep .eq. 1) then
c
        do j=1,nstat
        do iz=0,nnz
           stat(iz,j) = 0.0
        enddo
        enddo
        do iz=izs,ize
           stat(iz,1) = uwsb(iz)
           stat(iz,2) = vwsb(iz)
           stat(iz,3) = wwsb(iz)
!           stat(iz,4) = tr_tau(iz)
           stat(iz,5) = triz(iz)
           stat(iz,6) = shrz(iz)
           !stat(iz,7) = t_diss(iz)
        enddo
        m1 = js
        m2 = js + nscl
        m3 = js + 2*nscl
        do l=1,nscl
           do iz=izs,ize
              stat(iz,m1+l) = utsb(iz,l)
              stat(iz,m2+l) = vtsb(iz,l)
              stat(iz,m3+l) = wtsb(iz,l)
           enddo
        enddo

        !Populate the iz=0 locations if processor lies at the bottom -- only for certain quantities
        if (iss .eq. 0) then
           stat(0,1) = uwsb(0)
           stat(0,2) = vwsb(0)
           do l=1,nscl
              stat(0,m3+l) = wtsb(0,l)
           enddo
        end if

        call mpi_sum_z(stat(0,1),i_root,myid,nstat*(nnz+1),1)
        do iz=0,nnz
           uwsb(iz)   = stat(iz,1)
           vwsb(iz)   = stat(iz,2)
        enddo
        do iz=1,nnz
           wwsb(iz)   = stat(iz,3)
!           tr_tau(iz) = stat(iz,4)
           triz(iz)   = stat(iz,5)
           shrz(iz)   = stat(iz,6)
           !t_diss(iz) = stat(iz,7)
        enddo
        do l=1,nscl
           do iz=1,nnz
              utsb(iz,l) = stat(iz,m1+l)
              vtsb(iz,l) = stat(iz,m2+l)
           enddo
           do iz=0,nnz
              wtsb(iz,l) = stat(iz,m3+l)
           enddo
        enddo
        do iz=1,nnz
           buyz(iz) = batag*wtsb(iz,1)
        enddo

c
c -------- end if block
c
      endif
c
c ------- save old rhs in field variables for RK-advancement
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = urhs(ix,iy,iz)
            v(ix,iy,iz) = vrhs(ix,iy,iz)
            w(ix,iy,iz) = wrhs(ix,iy,iz)
            e(ix,iy,iz) = erhs(ix,iy,iz)
         enddo
         enddo
      enddo
      do iz=izs,ize
         do l=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,l,iz) = trhs(ix,iy,l,iz)
         enddo
         enddo
         enddo
      enddo
c
      return
      end
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
      use particles
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

!MODIFY BUOYANCY TERM TO INCLUDE HUMIDITY
!         w_buy = batag*(t(ix,iy,1,iz)*weit1 +
!     +                  t(ix,iy,1,izp1)*weit)

         w_buy = batag*(t(ix,iy,1,iz)*weit1 +
     +                  t(ix,iy,1,izp1)*weit)*
     +   (1.0+0.61*(weit1*t(ix,iy,2,iz) + 
     +    weit*t(ix,iy,2,izp1)))


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

         u_geo = -dpdx
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
         sfc_flx = 0.0
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


!---------add on simple radiation scheme:
      if (ilongwave == 1) then
      if (iscl == 1) then

         call calc_radiation

         do iz=izs,ize
           do iy=iys,iye
           do ix=1,nnx
           !calc_radiation fills "radsrc" vector
           r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) + radsrc(iz)

           end do
           end do
         end do


      end if
      end if

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
      vis_s(1:nnx,iys:iye,2,izs-1:ize+1) = nuf/Sc    !Dv =nu/Sc
      r5 = 0.0
      e = 0.0

      end
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
      use particles
      real sfc_flx(2+nscl)
c

          if (time .gt. 1800.0) then
             tsfcc(1) = 282.0
          else
             tsfcc(1)  = 283.9
          end if


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
            dnom3=t10xy(iscl)*windm
            if(dnom3 .ne. 0.) then
               dnom_i = 1.0/dnom3
               do iy=iys,iye
               do ix=1,nnx
                  taut3m(ix,iy,iscl)=aut3m(iscl)*
     +                 (windm*(t(ix,iy,iscl,iz)-t1xy(iscl))+
     +                  wind(ix,iy)*(t1xy(iscl)-tsfcc(iscl)))*dnom_i
               enddo
               enddo
            else
               do iy=iys,iye
               do ix=1,nnx
                  taut3m(ix,iy,iscl) = aut3m(iscl)
               enddo
               enddo
            endif
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
      do iscl=1,nscl
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
      subroutine lower_dns

      !Make lower BC by setting w = 0, u,v equal to mirror of interior points
      !Also set the scalar boundary condition based on either Neumann or Dirichlet

      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      implicit none

      integer :: iz,izm1,iscl,ix,iy
      real :: dz_i,dtdz,flux

      iz   = 1
      izm1 = iz - 1
      dz_i = dzu_i(1)

c ------------ initialize u, v, w, t and derivatives at izm1
c
      do iy=iys,iye
      do ix=1,nnx
         u(ix,iy,izm1)  = -2.0*Uo-u(ix,iy,iz)
         v(ix,iy,izm1)  = -v(ix,iy,iz) 
         w(ix,iy,izm1)  = 0.0 
         r3(ix,iy,izm1) = 0.0
         e(ix,iy,izm1)  = 0.0 
         ux(ix,iy,izm1) = 0.0
         uy(ix,iy,izm1) = 0.0
         vx(ix,iy,izm1) = 0.0
         vy(ix,iy,izm1) = 0.0
         wx(ix,iy,izm1) = 0.0 
         wy(ix,iy,izm1) = 0.0 
           wbc(ix,iy,2) = 0.0
           ebc(ix,iy,2) = 0.0
           ubc(ix,iy,2) = u(ix,iy,izm1)
           vbc(ix,iy,2) = v(ix,iy,izm1)
           pbc(ix,iy,2) = 0.0
           pbc2(ix,iy,2)= 0.0
         
      enddo
      enddo


!NOTE: sign convention is that wtsfc is flux upwards INTO domain (not just in vertical direction)
! Set the scalar boundary condition based on isfc
c            isfc = 0, specified surface heat flux (through qstar)
c                 = 1, specified surface temperature (Tbot)

      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            if (isfc(iscl)==1) then
            t(ix,iy,iscl,izm1) = 2.0*Tbot(iscl) - t(ix,iy,iscl,iz) 
            end if

            if (isfc(iscl)==0) then 
            flux = dzu(0)*wtsfc(iscl)/vis_s(ix,iy,iscl,izm1)
            t(ix,iy,iscl,izm1) = t(ix,iy,iscl,iz)     
     +       + flux  
            end if
         enddo
         enddo
      enddo

      end subroutine lower_dns
      subroutine upper_dns

      !Make upper BC by setting w = 0, u,v equal to mirror of interior points
      !Also set the scalar boundary condition based on either Neumann or Dirichlet

      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      implicit none

      integer :: iz,izp1,iscl,ix,iy
      real :: dz_i,dtdz,flux

      iz   = nnz
      izp1 = iz + 1

c
c ------------ initialize u, v, w, t and derivatives at izp1
c
      do iy=iys,iye
      do ix=1,nnx
         u(ix,iy,izp1)  = 2.0*Uo-u(ix,iy,iz)
         v(ix,iy,izp1)  = -v(ix,iy,iz) 
         w(ix,iy,izp1)  = 0.0 
         r3(ix,iy,izp1) = 0.0
         e(ix,iy,izp1)  = 0.0 
         ux(ix,iy,izp1) = 0.0
         uy(ix,iy,izp1) = 0.0
         vx(ix,iy,izp1) = 0.0
         vy(ix,iy,izp1) = 0.0
         wx(ix,iy,izp1) = 0.0 
         wy(ix,iy,izp1) = 0.0 
           wbc(ix,iy,1) = 0.0
           ebc(ix,iy,1) = 0.0
           ubc(ix,iy,1) = u(ix,iy,izp1)
           vbc(ix,iy,1) = v(ix,iy,izp1)
           pbc(ix,iy,1) = 0.0
           pbc2(ix,iy,1)= 0.0
      enddo
      enddo

! Set the scalar boundary condition based on isfc
c            isfc = 0, specified surface heat flux (through qstar)
c                 = 1, specified surface temperature (Ttop)

!NOTE: sign convention is that wtsfc is flux upwards INTO domain (not just in vertical direction)

      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            if (isfc(iscl)==1) then
            t(ix,iy,iscl,izp1) = 2.0*Ttop(iscl)-t(ix,iy,iscl,iz) 
            end if
            if (isfc(iscl)==0) then 
            flux = dzu(izp1)*wtsfc(iscl)/vis_s(ix,iy,iscl,izp1)
            t(ix,iy,iscl,izp1) = t(ix,iy,iscl,iz)      
     +       - flux 
            end if
         enddo
         enddo
      enddo

      end subroutine upper_dns
      subroutine lower_free(it)
c
c --------------- setup lower boundary condition for free
c                 convection where each processor applies
c                 log-law at several (ix,iy) for iz = 1
c
c                 index f(.,.,2)  indicates lower
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      real u_level1(nnx,iys:iye,2+nscl), buf(2+2*nscl)
      real sbuf(2+2*nscl,mxs:mxe,iys:iye)
      real rbuf((2+2*nscl)*nnx*(iye+1-iys))
c
c -------------- broadcast level 1 data everywhere
c
      if(iss .eq. 0) then
         do iy=iys,iye
         do ix=1,nnx
            u_level1(ix,iy,1) = u(ix,iy,1)
            u_level1(ix,iy,2) = v(ix,iy,1)
         enddo
         enddo
         do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            u_level1(ix,iy,2+iscl) = t(ix,iy,iscl,1)
         enddo
         enddo
         enddo
      endif
      num = nnx*(iye + 1 - iys)*(2+nscl)
c
c ------ send all of root data to other processors
c
      call mpi_send_root(u_level1(1,iys,1),
     +             num,myid,numprocs,ncpu_s)
c
c --------- every task gets their own fluxes and surface scalars
c
      call suft2(u_level1,it)
c
c --------- send surface scalars and momentum fluxes
c           back to root(s)
c
      if(numprocs .eq. 1) go to 999
c
      do iy=iys,iye
      do ix=mxs,mxe
         sbuf(1,ix,iy)  = tau13m(ix,iy)
         sbuf(2,ix,iy)  = tau23m(ix,iy)
      enddo
      enddo
      do iscl=1,nscl
      do iy=iys,iye
      do ix=mxs,mxe
         sbuf(2+iscl,ix,iy)      = taut3m(ix,iy,iscl)
         sbuf(2+nscl+iscl,ix,iy) = t_grnd(ix,iy,iscl)
      enddo
      enddo
      enddo
c
      irow_r = mod(myid,ncpu_s)
      if(myid .ge. ncpu_s) then
        num = (2+2*nscl)*(mxe+1-mxs)*(iye+1-iys)
        call mpi_send(sbuf(1,mxs,iys),num,mpi_real8,irow_r,1,
     +       mpi_comm_world,ierr)
      else
        do l=irow_r+ncpu_s,numprocs-1,ncpu_s
           num = (2+2*nscl)*(mx_e(l)+1-mx_s(l))*(iye+1-iys)
           call mpi_recv(rbuf(1),num,mpi_real8,l,1,
     +          mpi_comm_world,istatus,ierr)
c          call f_suft2(rbuf,maxnx,maxny,mx_s(l),mx_e(l),iys,iye,nscl,
           call f_suft2(rbuf,nnx,mx_s(l),mx_e(l),iys,iye,nscl,
     +                  tau13m,tau23m,taut3m,t_grnd)
        enddo
      endif
c
  999 continue
c
c ------------ only for root row = 0
c              get sums of surface conditions
c              and set surface boundary conditions
c
      if(iss .eq. 0) then
c
         buf(1) = 0.0
         buf(2) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            buf(1) = buf(1) + tau13m(ix,iy)
            buf(2) = buf(2) + tau23m(ix,iy)
         enddo
         enddo
         do iscl=1,nscl
            buf(2+iscl)      = 0.
            buf(2+nscl+iscl) = 0.
            do iy=iys,iye
            do ix=1,nnx
               buf(2+iscl)      = buf(2+iscl) + taut3m(ix,iy,iscl)
               buf(2+nscl+iscl) = buf(2+nscl+iscl) + t_grnd(ix,iy,iscl)
            enddo
            enddo
         enddo
c
         call mpi_sum_xy(buf,myid,iss,ise,2+2*nscl)
         uwsfc = buf(1)*fnxy
         vwsfc = buf(2)*fnxy
         do iscl=1,nscl
            wtsfc(iscl) = buf(2+iscl)*fnxy
            tsfcc(iscl) = buf(2+nscl+iscl)*fnxy
         enddo
c
         iz   = 1
         izm1 = iz - 1
         dz_i = dzu_i(iz)
c
         do iy=iys,iye
         do ix=1,nnx
            ebc(ix,iy,2)=amax1(e(ix,iy,iz),sml_eg)
            wbc(ix,iy,2)= 0.0
            pbc(ix,iy,2) = 0.0
            pbc2(ix,iy,2) = 0.0
         enddo
         enddo
c
         do iy=iys,iye
         do ix=1,nnx
            dudz     = 2.*u(ix,iy,iz)*dz_i
            dvdz     = 2.*v(ix,iy,iz)*dz_i
            ubc(ix,iy,2) = u(ix,iy,iz) - dudz*dzu(iz)
            vbc(ix,iy,2) = v(ix,iy,iz) - dvdz*dzu(iz)
         enddo
         enddo
         do iscl=1,nscl
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
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izm1) = tbc(ix,iy,iscl,2)
            enddo
            enddo
         enddo
c
c ----- end of if block for root row
c
      endif
c
 7999 continue
c
      return
      end
      subroutine f_suft2(rbuf,nnx,mxs,mxe,iys,iye,nscl,
     +                  tau13m,tau23m,taut3m,t_grnd)
c
c ------ fill surface arrays on root processors
c
      real rbuf(2+2*nscl,mxs:mxe,iys:iye)
      real tau13m(nnx,iys:iye), tau23m(nnx,iys:iye),
     +     taut3m(nnx,iys:iye,nscl), t_grnd(nnx,iys:iye,nscl)
c
      do iy=iys,iye
      do ix=mxs,mxe
         tau13m(ix,iy) = rbuf(1,ix,iy)
         tau23m(ix,iy) = rbuf(2,ix,iy)
      enddo
      enddo
      do iscl=1,nscl
         do iy=iys,iye
         do ix=mxs,mxe
            taut3m(ix,iy,iscl) = rbuf(2+iscl,ix,iy)
            t_grnd(ix,iy,iscl) = rbuf(2+nscl+iscl,ix,iy)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine upper
c
c ---- set boundary condition on upper boundary iz=nnz
c      option for special radiation boundary condition
c                 index f(.,.,1)  indicates upper. 
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      iz   = nnz
      izm1 = iz - 1
      izm2 = iz - 2
      izp1 = iz + 1
      izp2 = iz + 2
c
      if(ibcu .eq. 0) then
c
c --------- boundary conditions are gradient conditions
c
c            dudzbc = 0.0
c            dvdzbc = 0.0
c            dtdzbc = dtdzf
c            wbc    = 0.0
c            ebc    = 0.0
c
        do iy=iys,iye
        do ix=1,nnx
           wbc(ix,iy,1) = 0.0
           ebc(ix,iy,1) = 0.0
           ubc(ix,iy,1) = u(ix,iy,iz)
           vbc(ix,iy,1) = v(ix,iy,iz)
           pbc(ix,iy,1) = 0.0
           pbc2(ix,iy,1)= 0.0
        enddo
        enddo
        do iscl=1,nscl
c
c ---------- first get average scalar gradient
c
           dtdzf(iscl) = 0.0
           do iy=iys,iye
           do ix=1,nnx
              dtdzf(iscl) = dtdzf(iscl) + (t(ix,iy,iscl,nnz) -
     +                      t(ix,iy,iscl,nnz-1))*dzu_i(nnz)
           enddo
           enddo
           dtdzf(iscl) = dtdzf(iscl)*fnxy
        enddo
c
        call mpi_sum_xy(dtdzf,myid,iss,ise,nscl)
c
        do iscl=1,nscl
           do iy=iys,iye
           do ix=1,nnx
              tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) + 
     +                            dtdzf(iscl)*dzu(nnzp1)
           enddo
           enddo
        enddo
      else if(ibcu .eq. 1) then
c
c ------------- special if iradup boundary condition
c               get estimate of w from continuity and 
c               linearized relation for pressure
c
      xmeanp = 0.0
      grad_ug = ug(nnz) - ug((nnz-1))
      do iy=iys,iye
      do ix=1,nnx
         wbc(ix,iy,1) = w(ix,iy,izm1)-
     +                  (ux(ix,iy,iz)+vy(ix,iy,iz))*dzw(iz)
         pbc(ix,iy,1) = .5*(w(ix,iy,izm1)+wbc(ix,iy,1))
         ebc(ix,iy,1) = 0.0
         ubc(ix,iy,1) = u(ix,iy,iz) + grad_ug
         vbc(ix,iy,1) = v(ix,iy,iz)
         pbc2(ix,iy,1)=0.5*(u(ix,iy,iz)**2 + v(ix,iy,iz)**2) +
     +              0.25*(w(ix,iy,izm1)**2 + wbc(ix,iy,1)**2)
         xmeanp = xmeanp + pbc2(ix,iy,1)
      enddo
      enddo
      call mpi_sum_xy(xmeanp,myid,iss,ise,1)
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) + 
     +                          dtdzf(iscl)*dzu(nnzp1)
         enddo
         enddo
      enddo
      xmeanp = xmeanp*fnxy
      do iy=iys,iye
      do ix=1,nnx
         pbc2(ix,iy,1) = pbc2(ix,iy,1) - xmeanp
      enddo
      enddo
c
c ---------- end if block
c
      endif
c
      do iy=iys,iye
      do ix=1,nnx
         w(ix,iy,iz)   = wbc(ix,iy,1)
         e(ix,iy,iz)   = ebc(ix,iy,1)
         r3(ix,iy,iz)  = 0.0
         r5(ix,iy,iz)  = 0.0
         u(ix,iy,izp1) = ubc(ix,iy,1)
         v(ix,iy,izp1) = vbc(ix,iy,1)
c ------------- note w and e nnz+1 values are not needed
         w(ix,iy,izp1) = wbc(ix,iy,1)
         e(ix,iy,izp1) = ebc(ix,iy,1)
         r3(ix,iy,izp1)= 0.0
         r5(ix,iy,izp1)= 0.0
c
c ---------- set derivatives at top of box (wx,wy not needed)
c            ux,uy,vx,vy are used in e production, but neglect
c            at top of box becuase of bc
c
         wx(ix,iy,izp1) = 0.0
         wy(ix,iy,izp1) = 0.0
         ux(ix,iy,izp1) = 0.0
         uy(ix,iy,izp1) = 0.0
         vx(ix,iy,izp1) = 0.0
         vy(ix,iy,izp1) = 0.0
c        ux(ix,iy,izp1) = ubc(ix,iy,1)
c        uy(ix,iy,izp1) = ubc(ix,iy,1)
c        vx(ix,iy,izp1) = vbc(ix,iy,1)
c        vy(ix,iy,izp1) = vbc(ix,iy,1)
      enddo
      enddo
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,iscl,izp1) = tbc(ix,iy,iscl,1)
            t(ix,iy,iscl,izp2) = tbc(ix,iy,iscl,1)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine comp_p
c
c --------- setup pressure solver
c
      use pars
      use fftwk
      use fields
      use con_data
      use con_stats
      include 'mpif.h'
      real fnt1(nnx,iys:iye,izs:ize)
      real fs(nnx,iys:iye,2), fr(nnx,iys:iye,2)
      integer istatus(mpi_status_size)
c
      gami = 1.0/dtgama
c
      nb = myid - ncpu_s
      nt = myid + ncpu_s
c
c ------------ Send both r3 and updated w (from comp1)
c              to processor above the current myid.
c
      if(iss .eq. 0) then
         nb = mpi_proc_null
      endif
      if(ise .eq. numprocs-1) then
         nt = mpi_proc_null
      endif
      nsend = 2*nnx*(iye + 1 - iys)
      nrecv = nsend
      do iy=iys,iye
      do ix=1,nnx
         fs(ix,iy,1) = r3(ix,iy,ize)
         fs(ix,iy,2) = w(ix,iy,ize)
      enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nt,2,
     +     fr(1,iys,1),nrecv,mpi_real8,nb,2,
     +     mpi_comm_world,istatus,ierr)
      if(iss .ne. 0) then
         do iy=iys,iye
         do ix=1,nnx
            r3(ix,iy,izs-1) = fr(ix,iy,1)
            w(ix,iy,izs-1)  = fr(ix,iy,2)
         enddo
         enddo
      endif
c
c ----------- setup general pressure calculation
c             relies on rhs from step n-1 being included 
c             in velocity-arrays already
c
      do iz=izs,ize
         izm1 = iz -1
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = u(ix,iy,iz)*gami + r1(ix,iy,iz)
         enddo
         enddo
         call xderivp(fnt1(1,iys,iz),trigx(1,1),xk(1),nnx,iys,iye)
c
         if(iz .eq. 1) then
            do iy=iys,iye
            do ix=1,nnx
                p(ix,iy,iz) = fnt1(ix,iy,iz) +  
     +                     ((w(ix,iy,iz) -wbc(ix,iy,2))*gami +
     +                       r3(ix,iy,iz))*dzw_i(iz)
            enddo
            enddo
         else if(iz .eq. nnz) then
            do iy=iys,iye
            do ix=1,nnx
                p(ix,iy,iz) = fnt1(ix,iy,iz) + 
     +                     ((wbc(ix,iy,1) - w(ix,iy,izm1))*gami -
     +                      r3(ix,iy,izm1))*dzw_i(iz)
            enddo
            enddo
         else 
            do iy=iys,iye
            do ix=1,nnx
                p(ix,iy,iz) = fnt1(ix,iy,iz) + 
     +                    ((w(ix,iy,iz)  - w(ix,iy,izm1))*gami +
     +                      r3(ix,iy,iz) - r3(ix,iy,izm1))*dzw_i(iz)
            enddo
            enddo
         endif
c
c --------- end z loop
c
      enddo
c
c ----------- check for radiation boundary condition, all processors
c
      if(ibcu .eq. 1) then
        do iy=iys,iye
        do ix=1,nnx
           ptop(ix,iy,1) = pbc(ix,iy,1)
           ptop(ix,iy,2) = pbc2(ix,iy,1)
        enddo
        enddo
      endif
c
c --------- now y contribution
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = v(ix,iy,iz)*gami + r2(ix,iy,iz)
         enddo
         enddo
      enddo
c
      call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
             p(ix,iy,iz) = p(ix,iy,iz) + fnt1(ix,iy,iz) 
         enddo
         enddo
      enddo
c
      call pressure
c
      return
      end
      subroutine comp2
c
c ------- add p gradients to rhs. Use already defined p
c         at ize+1 to get w (see sr. pressure).
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real fnt1(nnx,iys:iye,izs:ize), fnt2(nnx,iys:iye)
      real r3_sum(1:nnz)
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      do iz=1,nnz
         r3_sum(iz) = 0.0
      enddo
c
c --------- dp/dy at all z
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = p(ix,iy,iz)
         enddo
         enddo
      enddo
c
      call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
      do iz=izs,ize
c
         izm1  = iz - 1
         izp1  = iz + 1
c
         do iy=iys,iye
         do ix=1,nnx
            fnt2(ix,iy) = p(ix,iy,iz)
         enddo
         enddo
         call xderivp(fnt2(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         do iy=iys,iye
         do ix=1,nnx
            r1(ix,iy,iz) = r1(ix,iy,iz) - fnt2(ix,iy)
            r2(ix,iy,iz) = r2(ix,iy,iz) - fnt1(ix,iy,iz)
         enddo
         enddo
         if (iz.ne.nnz) then
            do iy=iys,iye
            do ix=1,nnx
               r3(ix,iy,iz) = r3(ix,iy,iz) -
     +            (p(ix,iy,izp1)-p(ix,iy,iz))*dzu_i(izp1)
               r3_sum(iz) = r3_sum(iz) + r3(ix,iy,iz)
            enddo
            enddo
            r3_sum(iz) = r3_sum(iz)*fnxy
         endif
c
c ------------------------ time stepping with 3-order rk method
c                          first w variables
c
      if(iz .ne. nnz) then
         do iy=iys,iye
         do ix=1,nnx
c           w(ix,iy,iz)  = w(ix,iy,iz)+dtgama*r3(ix,iy,iz)
            e(ix,iy,iz)  = e(ix,iy,iz)+dtgama*r5(ix,iy,iz)
         enddo
         enddo
      else
c
c --------- update wout and eout by setting = to bc values
c
         do iy=iys,iye
         do ix=1,nnx
            w(ix,iy,iz)  = wbc(ix,iy,1)
            e(ix,iy,iz)  = ebc(ix,iy,1)
            r3(ix,iy,iz) = 0.0
            r5(ix,iy,iz) = 0.0
         enddo
         enddo
      endif
c
c -------- now all u-variables
c
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = u(ix,iy,iz)+dtgama*r1(ix,iy,iz)
            v(ix,iy,iz) = v(ix,iy,iz)+dtgama*r2(ix,iy,iz)
         enddo
         enddo
         do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,iscl,iz)  = t(ix,iy,iscl,iz)+
     +                          dtgama*r4(ix,iy,iscl,iz)
         enddo
         enddo
         enddo
c
c -------- end z loop
c
      enddo
c
c ---------- gather partial sums for w computation
c
      call mpi_sum_z(r3_sum,i_root,myid,nnz,1)
c
      do iz=izs,min(ize,nnz-1)
         do iy=iys,iye
         do ix=1,nnx
            r3(ix,iy,iz) = r3(ix,iy,iz) - r3_sum(iz)
            w(ix,iy,iz)  = w(ix,iy,iz) + dtgama*r3(ix,iy,iz)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine pressure
c
c -------- solve for pressure using a matrix transpose
c          across mpi tasks and tridiagonal solver. 
c          The transposed array
c          is dimensioned (0:nnz+1). Values 
c          (0 & nnz+1) are not needed but are useful in the 
c          matrix transpose when we return (see send_ztox).
c          On exit p is defined at all [izs-1:ize+1].
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real pfft(nny,jxs:jxe,izs-1:ize+1)
      real pt(0:nnz+1,jxs:jxe,jys:jye)
      real ptopfft(nny,jxs:jxe,1:2)
      real psum(1:nnz)
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
c ------------ Fourier analyze the right hand side
c              at all iz = izs,ize. results are in pfft
c
c


      call fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)

c
c ------------ Fourier analyze the radiation bc arrays
c
      if(ibcu .eq. 1) then
        call fft2d_mpi(ptop(1,iys,1),ptopfft(1,jxs,1),
     +           trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           1,2,myid,ncpu_s,numprocs,-2)
      endif
c
c ---------- transpose first and last index of array
c            the order of pfft is (y,x,z)
c
      call xtoz_trans(pfft,pt,nny,nnz,jys,jye,jy_s,jy_e,
     +                jxs,jxe,izs,ize,iz_s,iz_e,myid,ncpu_s,
     +                numprocs)


      call solve_trid(pt, ptopfft)

c
c ------------- transpose back
c
      call ztox_trans(pt,pfft,nny,nnz,jys,jye,jy_s,jy_e,
     +                jxs,jxe,izs,ize,iz_s,iz_e,myid,ncpu_s,
     +                numprocs)
c
      iz_ee = ize+1
      if(ise .eq. numprocs-1) then
         iz_ee = ize
      endif
c
c --------- inverse fft at all iz=izs,iz_ee to get p
c           see z indices
c
      call fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,iz_ee,myid,ncpu_s,numprocs,2)
c
c -------- partial sums for pressure
c
      do iz=1,nnz
         psum(iz) = 0.0
      enddo
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            psum(iz) = psum(iz) + p(ix,iy,iz)
         enddo
         enddo
         psum(iz) = psum(iz)*fnxy
      enddo
      call mpi_sum_z(psum,i_root,myid,nnz,1)
c
      do iz=izs,iz_ee
c        psum(iz) = -psum(iz) + engz(iz) + c23*engsbz(iz)
         do iy=iys,iye
         do ix=1,nnx
            p(ix,iy,iz) = p(ix,iy,iz) - psum(iz)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine solve_trid(pt, ptop)
c 
c --------- tridiagonal solver. odd order for ptop, ptop2
c           because of 2d-fft
c
      use pars
      use con_data
      use con_stats
c
      real ptop(nny,jxs:jxe,1:2)
      real pt(0:nnz+1,jxs:jxe,jys:jye)
      real aa(nnz,jxs:jxe),bb(nnz,jxs:jxe),
     +     dd(nnz,jxs:jxe),rh(nnz,jxs:jxe)
      real fac_u(nnz), fac_l(nnz), fac_a(nnz)
c
      do iz=1,nnz
         fac_u(iz) = 1.0/(dzw(iz)*dzu(iz+1))
         fac_l(iz) = 1.0/(dzw(iz)*dzu(iz))
         fac_a(iz) = fac_l(iz) + fac_u(iz)
      enddo
c
      do kp=jys,jye    
         do lp=jxs,jxe
         do iz=2,nnz-1
            bb(iz,lp)  = fac_l(iz)
            aa(iz,lp)  = fac_u(iz)
            dd(iz,lp)  = -xks(lp,kp) - fac_a(iz)
            rh(iz,lp)  = pt(iz,lp,kp)
         enddo
         enddo
c
c --------------- lower boundary, fill exterior pressure (not used)
c
         do lp=jxs,jxe
            bb(1,lp)  = 1.0
            aa(1,lp)  = fac_u(1)
            dd(1,lp)  = -xks(lp,kp) - fac_u(1)
            rh(1,lp)  = pt(1,lp,kp)
            pt(0,lp,kp) = 0.0
         enddo
c
c --------------- upper boundary, fill exterior pressure (not used)
c
         if(ibcu .eq. 1) then
            do lp=jxs,jxe
              bb(nnz,lp) = 0.0
              aa(nnz,lp) = 0.0
              dd(nnz,lp) = 1.0
              rh(nnz,lp) = ptop(kp,lp,1)*wavexy(lp,kp) + ptop(kp,lp,2)
              pt(nnz+1,lp,kp) = 0.0
            enddo
         else
            do lp=jxs,jxe
               bb(nnz,lp) = fac_l(nnz)
               aa(nnz,lp) = 1.0
               dd(nnz,lp) = -xks(lp,kp) - fac_l(nnz)
               rh(nnz,lp) = pt(nnz,lp,kp)
               pt(nnz+1,lp,kp) = 0.0
            enddo
         endif
c
c ---------------- special situation for zeroth mode
c                  makes mean pressure = 0
c
         if(kp .eq. 1 .and. jxs .eq. 1) then
           do iz=1,nnz
              dd(iz,1) = 1.0
              rh(iz,1) = 0.0
              aa(iz,1) = 0.0
              bb(iz,1) = 0.0
              dd(iz,2) = 1.0
              rh(iz,2) = 0.0
              aa(iz,2) = 0.0
              bb(iz,2) = 0.0
           enddo
         endif
c
c --------------- solve system
c


         call tridv(bb,dd,aa,rh,nnz,jxs,jxe)
         do lp=jxs,jxe
         do iz=1,nnz
            pt(iz,lp,kp) = rh(iz,lp)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine tridv(b,d,a,r,n,j1,j2)
c
c --- tridiagonal matrix solver with multiple vectors
c     (note j and i loops are reversed from cray version)
c
c --- input:   n   size of a,b,d and r
c              b   below diagonal elements (b(1) not used)
c              d   diagonal elements
c              a   above diagonal elements (a(n) not used)
c              r   right hand side
c              j1:j2  range of input vectors
c
c --- output:  r   solution vector
c
      real b(n,j1:j2), d(n,j1:j2), a(n,j1:j2), r(n,j1:j2)
c

      if(n .le. 1 ) then
         do j=j1,j2
            r(1,j) = r(1,j)/d(1,j)
         enddo
         go to 999
      endif
      do j=j1,j2
         d(1,j) = 1.0/d(1,j)
      enddo
      do j=j1,j2
      do i=2,n
         fac = b(i,j)*d(i-1,j)
         d(i,j) = 1.0/(d(i,j) - fac*a(i-1,j))
         r(i,j) = r(i,j) - fac*r(i-1,j)
      enddo
      enddo
      do j=j1,j2
         r(n,j) = r(n,j)*d(n,j)
      enddo
      do j=j1,j2
      do i=n-1,1,-1
         r(i,j) = d(i,j)*(r(i,j) - a(i,j)*r(i+1,j))
      enddo
      enddo
  999 continue
c
      return
      end
      subroutine get_derv
c
c ------- get ux,uy,vx,vy at all z for this node
c         using parallel fft. can be improved (?)
c         by using exchange to send derivatives
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      iz_ss = izs-1
      iz_ee = ize+1
      if(iss .eq. 0) then
         iz_ss = izs 
      endif
      if(ise .eq. numprocs-1) then
         iz_ee = ize
      endif
c
c ------- make sure <w> = 0
c
      do iz=izs-1,ize+1
         w_sum = 0.0
         do iy=iys,iye
         do ix=1,nnx
            w_sum = w_sum + w(ix,iy,iz)
         enddo
         enddo
         w_sum = w_sum*fnxy
         call mpi_sum_xy(w_sum,myid,iss,ise,1)
         do iy=iys,iye
         do ix=1,nnx
            w(ix,iy,iz) = w(ix,iy,iz) - w_sum
         enddo
         enddo
      enddo
c
c     do iz=iz_ss,iz_ee
      do iz=izs-1,ize+1
c        if(iz .eq. izs-1 .or. iz .eq. ize+1) then
c        do iy=iys,iye
c        do ix=1,nnx
c           ux(ix,iy,iz) = 0.0
c           vx(ix,iy,iz) = 0.0
c           wx(ix,iy,iz) = 0.0
c           uy(ix,iy,iz) = 0.0
c           vy(ix,iy,iz) = 0.0
c           wy(ix,iy,iz) = 0.0
c        enddo
c        enddo
c        else
         do iy=iys,iye
         do ix=1,nnx
            ux(ix,iy,iz) = u(ix,iy,iz)
            vx(ix,iy,iz) = v(ix,iy,iz)
            wx(ix,iy,iz) = w(ix,iy,iz)
            uy(ix,iy,iz) = u(ix,iy,iz)
            vy(ix,iy,iz) = v(ix,iy,iz)
            wy(ix,iy,iz) = w(ix,iy,iz)
         enddo
         enddo
c        endif
         call xderivp(ux(1,iys,iz),trigx(1,1),xk(1),
     +                 nnx,iys,iye)
         call xderivp(vx(1,iys,iz),trigx(1,1),xk(1),
     +                 nnx,iys,iye)
         call xderivp(wx(1,iys,iz),trigx(1,1),xk(1),
     +                 nnx,iys,iye)
      enddo
c
c ---------- get y derivatives for (u,v,w)
c
c     call yd_mpi(uy(1,iys,iz_ss),trigx(1,2),yk(1),
c    +           nnx,nny,ixs,ixe,ix_s,ix_e,
c    +           iys,iye,iy_s,iy_e,iz_ss,iz_ee,myid,ncpu_s,numprocs)
c     call yd_mpi(vy(1,iys,iz_ss),trigx(1,2),yk(1),
c    +           nnx,nny,ixs,ixe,ix_s,ix_e,
c    +           iys,iye,iy_s,iy_e,iz_ss,iz_ee,myid,ncpu_s,numprocs)
c     call yd_mpi(wy(1,iys,iz_ss),trigx(1,2),yk(1),
c    +           nnx,nny,ixs,ixe,ix_s,ix_e,
c    +           iys,iye,iy_s,iy_e,iz_ss,iz_ee,myid,ncpu_s,numprocs)
c
      call yd_mpi(uy(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
      call yd_mpi(vy(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
      call yd_mpi(wy(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
c
      return
      end
      subroutine get_means(istage)
c
c ------------ get means for all variables
c              for use in iso, surfvis, comp1, compmn.
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      do iz=0,nnz+1
         u_mn(iz)   = 0.0
         v_mn(iz)   = 0.0
         w_mn(iz)   = 0.0
         engz(iz)   = 0.0
         engsbz(iz) = 0.0
         divz(iz)   = 0.0
         pxym(iz)   = 0.0
      enddo
      do iscl=1,nscl
         do iz=0,nnz+1
            t_mn(iz,iscl) = 0.0
            alphaC(iz,iscl) = 0.0
         enddo
      enddo
      iz_ee = ize
      iz_ss = izs
      if(ize .eq. nnz) iz_ee = nnzp1
      if(izs .eq. 1)   iz_ss = 0 
      do iz=iz_ss,iz_ee
         do iy=iys,iye
         do ix=1,nnx
            u_mn(iz) = u_mn(iz) + u(ix,iy,iz)
            v_mn(iz) = v_mn(iz) + v(ix,iy,iz)
            w_mn(iz) = w_mn(iz) + w(ix,iy,iz)

         enddo
         enddo
         u_mn(iz) = u_mn(iz)*fnxy
         v_mn(iz) = v_mn(iz)*fnxy
         w_mn(iz) = w_mn(iz)*fnxy
         do iscl=1,nscl
            t_mn(iz,iscl) = 0.0
            alphaC(iz,iscl) = 0.0
            do iy=iys,iye
            do ix=1,nnx
               t_mn(iz,iscl) = t_mn(iz,iscl) + t(ix,iy,iscl,iz)
               alphaC(iz,iscl) = alphaC(iz,iscl)+vis_s(ix,iy,iscl,iz)
            enddo
            enddo
            t_mn(iz,iscl) = t_mn(iz,iscl)*fnxy
            alphaC(iz,iscl) = alphaC(iz,iscl)*fnxy
         enddo
      enddo
      call mpi_sum_z(u_mn(0),i_root,myid,nnzp1+1,1)
      call mpi_sum_z(v_mn(0),i_root,myid,nnzp1+1,1)
      call mpi_sum_z(w_mn(0),i_root,myid,nnzp1+1,1)
      do iscl=1,nscl
         call mpi_sum_z(t_mn(0,iscl),i_root,myid,nnzp1+1,1)
         call mpi_sum_z(alphaC(0,iscl),i_root,myid,nnzp1+1,1)
      enddo
c
c -------- set e to minimum value 
c
      do iz=izs-1,ize+1
         do iy=iys,iye
         do ix=1,nnx
            e(ix,iy,iz )=amax1(e(ix,iy,iz ),sml_eg)
         enddo
         enddo
      enddo
c
c -------------- get terms which contribute to mean pressure
c                careful with the sum, get the mean p_star pressure
c
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            e_temp     =  0.5*(e(ix,iy,iz) + e(ix,iy,izm1))
            q_temp     =  0.5*((u(ix,iy,iz) + stokes(iz))**2 +
     +                          v(ix,iy,iz)*v(ix,iy,iz) +
     +                       0.5*(w(ix,iy,iz)*w(ix,iy,iz) +
     +                            w(ix,iy,izm1)*w(ix,iy,izm1)))
            engz(iz)   = engz(iz) + q_temp
            engsbz(iz) = engsbz(iz) + e_temp
            pxym(iz)   = pxym(iz) + p(ix,iy,iz) - (c23*e_temp + q_temp)
         enddo
         enddo
         engz(iz)   = engz(iz)*fnxy
         engsbz(iz) = engsbz(iz)*fnxy
         pxym(iz)   = pxym(iz)*fnxy
      enddo
      call mpi_sum_z(engz(1),i_root,myid,nnzp1,1)
      call mpi_sum_z(engsbz(1),i_root,myid,nnzp1,1)
      call mpi_sum_z(pxym(1),i_root,myid,nnz,1)
c
c ------------ save means and divergence for printout and compmn
c              all cpus have means over all z
c
      if(istage .eq. 1) then
        do iz=izs,ize
           izm1 = iz - 1
           do iy=iys,iye
           do ix=1,nnx
              divz(iz) = divz(iz) + 
     +                  (ux(ix,iy,iz)+vy(ix,iy,iz)+
     +                  (w(ix,iy,iz)-w(ix,iy,izm1))*dzw_i(iz))**2
           enddo
           enddo
           divz(iz) = divz(iz)*fnxy
        enddo
        call mpi_sum_z(divz(1),i_root,myid,nnz,1)
c
        do iz=0,nnz+1
           uxym(iz) = u_mn(iz)
           vxym(iz) = v_mn(iz)
           wxym(iz) = w_mn(iz)
        enddo
        do iscl=1,nscl
           do iz=0,nnz+1
              txym(iz,iscl) = t_mn(iz,iscl)
           enddo
        enddo
      endif
c
      return
      end
      subroutine xderivp(ax,trigx,xk,nnx,iys,iye)
c
c -------- get multiple x derivatives using fftpack routines
c          use fftpack storage a0, (a1,b1), (a2,b2),...,an
c          assumes that first wavenumber xk(1) = 0.0
c
c          assumes that wavenumbers are normalized by number of points
c
      real xk(nnx), trigx(2*nnx+15), ax(nnx,iys:iye)
c
c     fn = 1.0/float(nnx)
      do iy=iys,iye
         call rfftf(nnx,ax(1,iy),trigx)
         ii = 1
         ax(1,iy) = 0.0
         ax(nnx,iy) = 0.0
         do ix=2,nnx-1,2
            ii          = ii + 1
            temp        = ax(ix,iy)
            ax(ix,iy)   = -xk(ii)*ax(ix+1,iy)
            ax(ix+1,iy) = xk(ii)*temp
         enddo
         call rfftb(nnx,ax(1,iy),trigx)
      enddo
c
      return
      end
      subroutine fft2d_mpi(ax,at,trigx,trigc,nx,ny,
     +           jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           iz1,iz2,myid,ncpu,np,isgn)
c
c -------- get 2d fft using fftpack routines and parallel mpi
c          use fftpack storage a0, (a1,b1), (a2,b2),...,
c
c         isgn = -1 do forward transform, get coefficients
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
c
c         isgn = -2 do forward transform, get coefficients
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is at(ny,jxs:jxe,iz1:iz2)
c
c         isgn =  1 do inverse transform, move to physical space
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
c
c         isgn =  2 do inverse transform, move to physical space
c                   incoming array is at(ny,jxs:jxe,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
c
      real ax(nx+2,iys:iye,iz1:iz2), at(ny,jxs:jxe,iz1:iz2),
     +     trigx(2*nx+15), trigc(4*ny+15),
     +     a2d(2,ny), a_wrk(nx)
      integer jx_s(0:np-1), jx_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      nxp2 = nx + 2
      if(isgn .lt. 0) then
         fn   = 1.0/(float(nx)*float(ny))
c
c ------ 1d fft in x over [iys,iye] for all z
c
         do iz=iz1,iz2
            do iy=iys,iye
               do ix=1,nx
                  a_wrk(ix) = ax(ix,iy,iz)*fn
               enddo
               call rfftf(nx,a_wrk(1),trigx(1))
               ax(1,iy,iz) = a_wrk(1)
               ax(2,iy,iz) = 0.0
               do ix=2,nx
                  ax(ix+1,iy,iz) = a_wrk(ix)
               enddo
               ax(nx+2,iy,iz) = 0.0
            enddo
         enddo
         call xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +        iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c ------ 1d fft in y over [jxs,jxe] for all z
c
         do iz=iz1,iz2
            do ix=jxs,jxe,2
               do iy=1,ny
                  a2d(1,iy) = at(iy,ix,iz)
                  a2d(2,iy) = at(iy,ix+1,iz)
               enddo
               call cfftf(ny,a2d(1,1),trigc(1))
               do iy=1,ny
                  at(iy,ix,iz)   = a2d(1,iy)
                  at(iy,ix+1,iz) = a2d(2,iy)
               enddo
            enddo
         enddo
c
c ---- decide whether to transpose back or leave as is
c
         if(isgn .eq. -1) then
            call ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
         endif
c
      else
c
c ---- decide whether to first transpose or leave as is
c
         if(isgn .eq. 1) then
            call xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
         endif
c
c ------ 1d fft in y over [jxs,jxe] for all z
c
         do iz=iz1,iz2
            do ix=jxs,jxe,2
               do iy=1,ny
                  a2d(1,iy) = at(iy,ix,iz)
                  a2d(2,iy) = at(iy,ix+1,iz)
               enddo
               call cfftb(ny,a2d(1,1),trigc(1))
               do iy=1,ny
                  at(iy,ix,iz)   = a2d(1,iy)
                  at(iy,ix+1,iz) = a2d(2,iy)
               enddo
            enddo
         enddo
         call ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +        iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c ------  1d fft in x over [iys,iye] for all z
c
         do iz=iz1,iz2
            do iy=iys,iye
               a_wrk(1) = ax(1,iy,iz)
               do ix=2,nx
                  a_wrk(ix) = ax(ix+1,iy,iz)
               enddo
               call rfftb(nx,a_wrk(1),trigx(1))
               do ix=1,nx
                  ax(ix,iy,iz) = a_wrk(ix)
               enddo
            enddo
         enddo
      endif
c
      return
      end
      subroutine yderiv(ay,trigy,yk,nnx,nny)
c
c -------- get multiple y derivatives using fftpack routines
c          use fftpack storage a_0, (a1,b1), (a2,b2), ...,
c          assumes that first wavenumber yk(1) = 0.0
c
c          assumes that wavenumbers are normalized by number of points
c
      real yk(nny), trigy(2*nny+15), ay(nnx,nny)
      real a_trans(nny)
c
c     fn = 1.0/float(nny)
      do ix=1,nnx
         do iy=1,nny
            a_trans(iy) = ay(ix,iy)
         enddo
         call rfftf(nny,a_trans(1),trigy)
         ii = 1
         a_trans(1)   = 0.0
         a_trans(nny) = 0.0
         do iy=2,nny-1,2
            ii            = ii + 1
            temp          = a_trans(iy)
            a_trans(iy)   = -yk(ii)*a_trans(iy+1)
            a_trans(iy+1) = yk(ii)*temp
         enddo
         call rfftb(nny,a_trans(1),trigy)
         do iy=1,nny
            ay(ix,iy) = a_trans(iy)
         enddo
      enddo
c
      return
      end
      subroutine yd_mpi(ay,trigy,yk,
     +           nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c -------- get multiple y derivatives using fftpack routines and mpi
c          use fftpack storage a_0, (a1,b1), (a2,b2), ...,
c          assumes that first wavenumber yk(1) = 0.0
c          wavenumbers are normalized by number of points, ny
c
      real yk(ny), trigy(2*ny+15), ay(nx,iys:iye,iz1:iz2)
      real ayt(ny,ixs:ixe,iz1:iz2)
c
      integer ix_s(0:np-1), ix_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      call xtoy_trans(ay,ayt,nx,ny,ixs,ixe,ix_s,ix_e,
     +         iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c     fn = 1.0/float(nny)
      do iz=iz1,iz2
         do ix=ixs,ixe
            call rfftf(ny,ayt(1,ix,iz),trigy)
            ii = 1
            ayt(1,ix,iz)  = 0.0
            ayt(ny,ix,iz) = 0.0
            do iy=2,ny-1,2
               ii              = ii + 1
               temp            = ayt(iy,ix,iz)
               ayt(iy,ix,iz)   = -yk(ii)*ayt(iy+1,ix,iz)
               ayt(iy+1,ix,iz) = yk(ii)*temp
            enddo
            call rfftb(ny,ayt(1,ix,iz),trigy)
         enddo
      enddo
      call ytox_trans(ayt,ay,nx,ny,ixs,ixe,ix_s,ix_e,
     +         iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
      return
      end
      function rlim(d1,d2,d3)
c
c ------------- Cees's kappa=1/3 scheme
c
      r = (d1-d2+1.e-100)/(d2-d3-1.e-100)
      rlim = (d2-d3)*amax1(0.,amin1(r,amin1(1./6.+1./3.*r,1.)))
c
c ------------- Cees's kappa=-1 scheme
c
c     r = (d1-d2+1.e-100)/(d2-d3-1.e-100)
c     rlim = (d2-d3)*amin1(abs(r),0.5)
c
c ------------- first order upwind
c
c     rlim = 0.0
c
c ------------- QUICK scheme
c
c     rlim = -0.25*d2 - 0.125*d3 + 0.375*d1
c
      return
      end
      function ran1(idum)
c
c ----------- stolen from numerical recipes,p. 271
c
      integer idum, ia, im, iq, ir, ntab, ndiv
      real ran1, am, eps, rnmx
      parameter (ia=16807,im=2147483647,am=1.0/im,iq=127773,ir=2836.0,
     +           ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-07,rnmx=1.0-eps)
      integer j, k, iv(ntab), iy
      save iv, iy
      data iv /ntab*0/, iy /0/
      if(idum .le. 0 .or. iy .eq. 0) then
         idum = max(-idum,1)
         do j=ntab+8,1,-1
            k = idum/iq
            idum = ia*(idum - k*iq) - ir*k
            if(idum .lt. 0) idum = idum + im
            if(j .le. ntab) iv(j) = idum
         enddo
         iy = iv(1)
      endif
      k     = idum/iq
      idum  = ia*(idum - k*iq) - ir*k
      if(idum .lt. 0) idum = idum + im
      j     = 1 + iy/ndiv
      iy    = iv(j)
      iv(j) = idum
      ran1  = min(am*iy, rnmx)
c
      return
      end
      function ranf()
      data inc /1/
      save inc, ix, ia, m, fm
      if(inc.eq.1) then
        inc = 2
        m = 2**20
        fm = float(m)
        ix = 566387
        ia = 2**10 + 3
      endif
      ix = mod(ia*ix,m)
      fx = float(ix)
      ranf = fx/fm
      return
      end
      subroutine stokesv
c
c ----------- get stokes drift velocity for assumed wavelength stokesw
c             and wave amplitude stokesa. Changed sign for z.
c
c
      use pars
      use con_data
      use con_stats
      include 'mpif.h'
c
      if(iocean .eq. 1) then
c
c ----------- compute stokes velocity for ocean pbls
c
c        stokesw = pi2/20.0
         stokesw = pi2/76.5
c        ak      = 0.04
         ak      = 0.00
c        stokesa = 1.0
         stokesa = ak/stokesw
         sigma = sqrt(abs(grav)*stokesw)
         stokess = sigma*stokesw*stokesa**2
         do iz=1,nnzp1
            stokes(iz) = stokess*exp(2.0*stokesw*zz(iz))
         enddo
         if(l_root) then
            write(6,6000) (iz,zz(iz),stokes(iz),iz=1,nnz)
 6000       format(' iz ',10x,' zz',10x,' stokes',/,(1x,i3,2e12.4))
         endif
c
      else
c
c ----------------- set stokes velocity = 0 for atmos. pbls 
c
         do iz=1,nnzp1
            stokes(iz) = 0.0
         enddo
         stokess = 0.0
         udrift = 0.0
         vdrift = 0.0
      endif
c
      return
      end
      subroutine busngr(zeta,phim,phis,psim,psis)
c
c ---- Businger's version of similarity theory
c
      data pih /1.57079633/
      save pih
c
      if(zeta .lt. 0.) then
         x=(1.0 - 15.0*zeta)**0.25
         phim = 1.0/x
         psim = 2.0*alog((1.0+x)/2.0) + alog((1.0+x*x)/2.0) - 
     +          2.0*atan(x)+pih
         if(psim.gt.2.0)psim=2.0
         y = sqrt(1.0-9.0*zeta)
         phis = 0.74/y
         psis = alog((1.0+y)/2.0)*2.0
      else if(zeta .gt. 0) then
         phim = 1.0 + 4.7*zeta
         phis = 0.74 + 4.7*zeta
         psim = -4.7*zeta
         psis = -4.7*zeta
      else
         phim = 1.0
         phis = 0.74
         psim = 0.0
         psis = 0.0
      endif
      return
      end
      subroutine fzol(zeta,phim,phis,psim,psis)
c        estimate the stability functions for momentum, m
c                                         and scalars,  c
c        from input of the stability parameter zeta = z/L

      data c1/5./
      data a3,b3,a4,b4/1.258,8.382,-28.862,98.9545/
      data zetam,zetas/-0.2,-1.0/
      save c1, a3, b3, a4, b4, zetam, zetas
c
      psimu(Y)  = 1.571 + 2.0*(alog(0.5*(1.0 + Y)) - atan(Y)) + 
     +            alog(0.5 + 0.5*Y**2)
      psisu(Y)  = 2.0*alog(0.5 + 0.5*Y)
      psicu(Y,G)= (1.0 - G)*alog(abs(Y - 1.0))
     +          + 0.5*(G + 2.0)*alog(abs(Y**2 + Y + 1.0))
     +          - (2.0*G + 1.0) / sqrt(3.0) * 
     +            atan((Y + 0.5)*2.0/sqrt(3.0))
      Xm(zol)   = (1.0 - 16.0*zol)**0.25
      Xs(zol)   = sqrt(1.0 - 16.0*zol)
      Xc(zol,f) =  abs(1.0 - f*zol)**(4.0/3.0)/(1.0 - f*zol)
c
      if(zeta.ge.0.0)       then
c                                          STABLE
      if(zeta.le.1.0) then
        phim = 1.0 + c1 * zeta
        psim = - c1 * zeta
        phis = phim
        psis = psim
                      else
c                                   use limiting form
        phim = c1 + zeta
        psim = (1.0 - c1)*(1.0 + alog(zeta) ) - zeta
        phis = phim
        psis = psim
                      endif

                            else
c                                         UNSTABLE
c                                                  momentum         
       if(zeta.ge.zetam) then
         phim = 1.0 / Xm(zeta)
         psim = psimu(Xm(zeta))
                         else
c                            use convective limit for momentum
         X = (1.0 - b3/a3 * zeta)**(1.0/3.0)

         fm = a3**(-1.0/3.0)
         phim = fm / Xc(zeta,b3/a3) 
         psim = psimu(Xm(zetam))
     *        + psicu(Xc(zeta,b3/a3),fm)
     *        - psicu(Xc(zetam,b3/a3),fm)
                         endif
      
c                                         UNSTABLE scalars
       if(zeta.ge.zetas) then
         phis = 1.0/Xs(zeta)
         psis = psisu(Xs(zeta))
                         else
c                              use convective limit for scalars
         fs =   abs(a4)**(-1.0/3.0)*abs(a4)/a4
         phis = (a4 - b4*zeta)**(-1.0/3.0)
         psis = psisu(Xs(zetas))
     *        + psicu(Xc(zeta,b4/a4),fs)
     *        - psicu(Xc(zetas,b4/a4),fs)
                         endif
               
                            endif
       return
       end
      subroutine suft(it)
c
c ---------- iterate for zeta = z/L using bisection method
c            either businger or large functions can be specified
c
c            isfc = 0, specified surface heat flux
c                 = 1, specified surface temperature
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
      real buf(3+nscl)
      real :: mod_Magnus
c
      parameter (iter_mo = 30, zeta_min = -6.0, zeta_max = 3.0)


c
c ---------- limiting value for wind
c
      ufree = 0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)
c
c ---- save old utau
c
      utausv = utau
      utau2  = utau*utau
c
      iz   = 1
      izp1 = iz + 1
      izm1 = iz - 1
c
      buf(1)  = 0.0
      buf(2)  = 0.0
      buf(3)  = 0.0
      tol     = 0.01
      do iy=iys,iye
      do ix=1,nnx
         buf(1) = buf(1) + u(ix,iy,iz)
         buf(2) = buf(2) + v(ix,iy,iz)
         wind(ix,iy) = sqrt((u(ix,iy,iz)+ugal)**2
     +                    +v(ix,iy,iz)*v(ix,iy,iz))
         buf(3) = buf(3) + wind(ix,iy)
      enddo
      enddo
      do iscl=1,nscl
         buf(3+iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            buf(3+iscl) = buf(3+iscl) + t(ix,iy,iscl,iz)
         enddo
         enddo
      enddo
c
c -------- get x-y slab sums
c
      call mpi_sum_xy(buf,myid,iss,ise,(3+nscl))
      u1xy  = buf(1)*fnxy + ugal
      v1xy  = buf(2)*fnxy
      windm = buf(3)*fnxy
      do iscl=1,nscl
         t1xy(iscl) = buf(3+iscl)*fnxy
      enddo
      vsfc  = sqrt(u1xy*u1xy+v1xy*v1xy)
      windm = amax1(windm,ufree)
      vsfc  = amax1(vsfc,ufree)
c
c ---------- limits for zeta
c
      zeta_mn = zeta_min
      zeta_mx = zeta_max
      if(isfc(1) .eq. 0) then
         f_con = z1*batag*vk*qstar(1)/((windm*vk)**3)
      else
         d_theta = vk74in*(tsfcc(1) - t1xy(1))
         f_con   = z1*batag*vk*d_theta/((windm*vk)**2)
      endif
c
c --------- iteration for zeta
c
      do iter=1,iter_mo
         zeta_a = 0.5*(zeta_mn + zeta_mx)
         if(ismlt .eq. 1) then
             call busngr(zeta_a,phim,phis,psim,psis)
         else
             call fzol(zeta_a,phim,phis,psim,psis)
         endif
         u_fac = (zody - psim)
         if(isfc(1) .eq. 0) then
            f_new =  zeta_a + f_con*u_fac**3
         else
            t_fac = 1.0/(zody - psis)
            f_new =  zeta_a + f_con*u_fac*u_fac*t_fac
         endif
         if(f_new .lt. 0.0) then
            zeta_mn = zeta_a
         else
            zeta_mx = zeta_a
         endif
c
c ----------- iteration details
c
         utau      = windm*vk/(zody-psim)
!         write(nprt,1000) iter, zeta_a, utau, phim, psim
 1000    format(' 1000 iter = ',i5,' zeta = ',e15.6,' u_* = ',e15.6,
     +          ' phim = ',e15.6,' psim = ',e15.6)
      enddo
c
c --------- check if neutral surface layer
c

      if (ibuoy.eq.0 .or. abs(qstar(1)) .lt. 1.0e-6) then
          amonin    = 1000.
          zeta      = 0.0
          utau      = windm*vk/zody
          thstar(1) = 0.0
          t10xy(1)  = 0.0

          !tsfcc(1)  = t1xy(1)
            dnom = zody*vk74in
            thstar(1) = (t1xy(1) - tsfcc(1))/dnom
            t10xy(1)  = thstar(1)*dnom
            qstar(1)  = -utau*thstar(1)
            write(*,*) 'DHR INSIDE',thstar(1),t10xy(1),qstar(1)
      else
         utau = windm*vk/(zody-psim)
         dnom = (zody-psis)*vk74in
         if(isfc(1) .eq. 0) then
            thstar(1) = -qstar(1)/utau
            tsfcc(1)  = t1xy(1)-thstar(1)*dnom
            t10xy(1)  = thstar(1)*dnom
         else
            thstar(1) = (t1xy(1) - tsfcc(1))/dnom
            t10xy(1)  = thstar(1)*dnom
            qstar(1)  = -utau*thstar(1)
         endif
         amonin = -utau**3/(batagk*qstar(1))
         zeta   = z1/amonin
      endif
c
c ------- surface details, for debug
c
!      write(nprt,2000) windm, utau, qstar(1), tsfcc(1), amonin, zeta,
!     +              z1, batag, vk, batagk, zo
 2000 format(' 2000 suft ',/,
     +       '    windm = ',e15.6,' utau = ',e15.6,' qstar = ',e15.6,/,
     +       '    tsfcc = ',e15.6,' MO L = ',e15.6,' z1/L = ',e15.6,/,
     +       '    z1 = ',e15.6,' batag = ',e15.6,' vk = ',e15.6,/,
     +       '    batagk = ',e15.6,' zo = ',e15.6)
c
      if (utau.gt.10.0) then
         write(6,9000)
         write(6,9200) utau,windm
         go to 9999
      endif
      if (t10xy(1).gt.0. .and. qstar(1) .gt. 0.) then
         write(6,9000)
         write(6,9300) u1xy,v1xy,t1xy(1),
     +                 tsfcc(1),amonin,utau,it
         go to 9999
      endif

!Now that temperature and momentum have been done, apply BCs on remaining scalars:

      !Treat humidity differently
      iscl=2

      if (ibuoy.eq.0 .or. abs(qstar(1)) .lt. 1.0e-6) then
          
          !thstar(iscl) = 0.0
          !t10xy(iscl)  = 0.0
          !tsfcc(iscl)  = t1xy(iscl)

          dnom = zody*vk74in
          tsfcc(2)=surf_RH/100.0*
     +    Mw/Ru/tsfcc(1)*mod_Magnus(tsfcc(1))/rhoa
            thstar(iscl) = (t1xy(iscl) - tsfcc(iscl))/dnom
            t10xy(iscl)  = thstar(iscl)*dnom
            qstar(iscl)  = -utau*thstar(iscl)
      else
         dnom = (zody-psis)*vk74in
         if(isfc(iscl) .eq. 0) then
            thstar(iscl) = -qstar(iscl)/utau
            tsfcc(iscl)  = t1xy(iscl)-thstar(iscl)*dnom
            t10xy(iscl)  = thstar(iscl)*dnom
         else
            !Humidity is always a Dirichlet condition at surf_RH
            !Based on temperature in tsfcc(1)
            tsfcc(2)=surf_RH/100.0*
     +      Mw/Ru/tsfcc(1)*mod_Magnus(tsfcc(1))/rhoa

            thstar(iscl) = (t1xy(iscl) - tsfcc(iscl))/dnom
            t10xy(iscl)  = thstar(iscl)*dnom
            qstar(iscl)  = -utau*thstar(iscl)
         endif
         amonin = -utau**3/(batagk*qstar(iscl))
         zeta   = z1/amonin
      endif


      !Other scalars:
      do iscl=3,nscl

      if (ibuoy.eq.0 .or. qstar(1) .eq. 0.) then
          thstar(iscl) = 0.0
          t10xy(iscl)  = 0.0
          tsfcc(iscl)  = t1xy(iscl)
      else
         dnom = (zody-psis)*vk74in
         if(isfc(iscl) .eq. 0) then
            thstar(iscl) = -qstar(iscl)/utau
            tsfcc(iscl)  = t1xy(iscl)-thstar(iscl)*dnom
            t10xy(iscl)  = thstar(iscl)*dnom
         else
            thstar(iscl) = (t1xy(iscl) - tsfcc(iscl))/dnom
            t10xy(iscl)  = thstar(iscl)*dnom
            qstar(iscl)  = -utau*thstar(iscl)
         endif
         amonin = -utau**3/(batagk*qstar(iscl))
         zeta   = z1/amonin
      endif

      enddo

c ---------- examples of two other scalars
c
c     c
c     c **** get flux of b scalar, specified surface value
c     c
c           dnom      = (zody-psis)*vk74in
c           thstar(2) = (t1xy(2)-tsfcc(2))/dnom
c           qstar(2)  = -thstar(2)*utau
c           t10xy(2)  = thstar(2)*dnom
c           aut3m(2)  =  qstar(2)
c
c **** get surface value of c scalar, specified surface flux
c
c     dnom      = (zody-psis)*vk74in
c     thstar(2) = -qstar(2)/utau
c     tsfcc(2)  = t1xy(2) - dnom*thstar(2)
c     t10xy(2)  = thstar(2)*dnom
c     aut3m(2)  = qstar(2)
c
      zol = zeta
      hol = zol*zi/z1
c
c ---- note roundoff problem in angles if close to multiples of pi
c
      tep = u1xy/vsfc
      if(tep.gt.1.)  tep = 1.0
      if(tep.lt.-1.) tep = -1.0
      thta      = acos(tep)
      utau2     = utau*utau
      au13m     = -utau2*cos(thta)
      au23m     = -utau2*sin(thta)*sign(1.,v1xy)
      do iscl = 1,nscl
      aut3m(iscl)  =  qstar(iscl)
      end do
c
      return
c
c -------- iteration did not converge
c
 9999 continue
 9000 format(' Trouble in SR. suft')
 9200 format(' Stop because utau = ',e15.6,' windm = ',e15.6)
 9300 format(' ** CHECK SFC U = ',e15.6,' V=',e15.6,' T,TS = ',2e15.6,
     +       ' L =',e15.6,' U_* = ',e15.6,' AT IT = ',i5)
      call mpi_finalize(ierr)
      stop
      end
      subroutine sufto(it)
c
      use pars
      use fields
      use con_data
      use con_stats
      real buf(3+nscl)
c
c ------- version of similarity theory adpated for ocean flows
c      option to use businger or large version of similarity theory
c
      iz    = 1
      izm1  = iz - 1
      izp1  = iz + 1
      z1_a  = abs(z1)
      buf(1)  = 0.0
      buf(2)  = 0.0
      buf(3)  = 0.0
      tol     = 0.01
      do iy=iys,iye
      do ix=1,nnx
         buf(1) = buf(1) + u(ix,iy,iz)
         buf(2) = buf(2) + v(ix,iy,iz)
         wind(ix,iy) = sqrt((u(ix,iy,iz)+ugal)**2
     +                    +v(ix,iy,iz)*v(ix,iy,iz))
         buf(3) = buf(3) + wind(ix,iy)
      enddo
      enddo
      do iscl=1,nscl
         buf(3+iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            buf(3+iscl) = buf(3+iscl) + t(ix,iy,iscl,iz)
         enddo
         enddo
      enddo
c
c -------- get x-y slab sums
c
      call mpi_sum_xy(buf,myid,iss,ise,(3+nscl))
      u1xy  = buf(1)*fnxy + ugal
      v1xy  = buf(2)*fnxy
      windm = buf(3)*fnxy
      do iscl=1,nscl
         t1xy(iscl) = buf(3+iscl)*fnxy
      enddo
      vsfc  = sqrt(u1xy*u1xy+v1xy*v1xy)
      windm = amax1(windm,ufree)
      vsfc  = amax1(vsfc,ufree)
c
      t10xy(1)=-qstar(1)/utau*zody*vk74in
c
c ---- check for temperature boundary condition
c
      if(isfc(iscl) .eq. 0 ) then
         tsfcc(1)=t1xy(1)-t10xy(1)
      endif
c     vsfc=sqrt(u1xy*u1xy+v1xy*v1xy)
c     if(windm.le.0.01)windm=0.01
c     if(vsfc .le.0.01)vsfc =0.01
c
c ----------- input surface wind stress (tau = 0.0184n/m*m)
c             density rho = 1000kg/m^3
c          
c     utau = 4.29e-03
c     utau = 6.10e-03
      utau = 7.00e-03
c
c **** save old utau
      utausv = utau
      utau2  = utau*utau
      if (ibuoy.eq.0 .or. qstar(1) .eq. 0.) then
          amonin    = 1000.
          zeta      = 0.
          thstar(1) = 0.0
          t10xy(1)  = 0.0
      else
          amonin = -utau2*utau/(batagk*qstar(1))
          zeta   = z1_a/amonin
      endif
      if (t10xy(1).lt.0. .and. qstar(1) .lt. 0.) then
         write(6,1234)u1xy,v1xy,t1xy(1),tsfcc(1),amonin,utau,it
 1234    format(' ** check sfc u=',e12.3,' v=',e12.3,' t,ts=',2f10.3,
     +     ' l=',e12.3,' u*=',e12.3,' at it=',i5)
         go to 9999
      endif
c
c -------- for stable,neutral and unstable pbl get drift velocity
c
      if(ismlt .eq. 1) then
          call busngr(zeta,phim,phis,psim,psis)
      else
          call fzol(zeta,phim,phis,psim,psis)
      endif
      udrift = windm + stokes(1) - stokess + utau*(zody-psim)*vkin
      vdrift = 0.0
      dnom      = (zody-psis)*vk74in
      if (isfc(iscl).eq.1) then
         thstar(1) = (t1xy(1) - tsfcc(1))/dnom
         t10xy(1)  = thstar(1)*dnom
         qstar(1)  = - utau*thstar(1)
      else
         thstar(1)  = -qstar(1)/utau
         tsfcc(1)   = t1xy(1)-thstar(1)*dnom
         t10xy(1)   = thstar(1)*dnom
      endif
      zol = zeta
      hol = zol*zi/z1
c
c ---------- examples of two other scalars
c
c     c
c     c **** get flux of b scalar, specified surface value
c     c
c           dnom      = (zody-psis)*vk74in
c           thstar(2) = (t1xy(2)-tsfcc(2))/dnom
c           qstar(2)  = -thstar(2)*utau
c           t10xy(2)  = thstar(2)*dnom
c           aut3m(2)  =  qstar(2)
c     c
c     c **** get surface value of c scalar, specified surface flux
c     c
c           dnom      = (zody-psis)*vk74in
c           thstar(3) = -qstar(3)/utau
c           tsfcc(3)  = t1xy(3) - dnom*thstar(3)
c           t10xy(3)  = thstar(3)*dnom
c           aut3m(3)  = qstar(3)
c
c **** note roundoff problem in angles are close to multiples of pi
c     tep=u1xy/vsfc
c     if(tep.gt.1.)tep=1.
c     if(tep.lt.-1.)tep=-1.
c     thta=acos(tep)
      utau2 = utau*utau
c     au13m = -utau2*cos(thta)
c     au23m = -utau2*sin(thta)*sign(1.,v1xy)
      au13m = utau2
      au23m = 0.0
      aut3m(1)= qstar(1)
c
      return
c
c --------- trouble in sl routine
c
 9999 continue
c
      write(nprt,9000)
 9000 format(' Trouble in SR. sufto')
      call mpi_finalize(ierr)
      stop
      end
      subroutine suft2(u_level1,it)
c
      use pars
      use fields
      use con_data
      use con_stats
c
      real u_level1(nnx,iys:iye,2+nscl)
c
c     u_level1(.,.,1) = u
c     u_level1(.,.,2) = v
c     u_level1(.,.,3) = theta
c     u_level1(.,.,4) = more scalars
c
      tol = 0.01
      ufree=0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)
      zeta_mn = -6.0
      zeta_mn_i = 1.0/zeta_mn
      iz   = 1
c     izm1 = iz - 1
c     izp1 = iz + 1
c
c      write(nprt,3131) myid, utau, zody, vk74in, batagk,
c    +               u_level1(jxs,1,3), u_level1(jxe,1,3),
c    +               u_level1(jxs,1,1), u_level1(jxe,1,1),
c    +               u_level1(jxs,1,2), u_level1(jxe,1,2)
c3131  format(' in suft2 myid = ',i4,/,
c    +        ' utau = ',e15.6,' zody = ',e15.6,/,
c    +        ' vk74in = ',e15.6,' batagk = ',e15.6,/,
c    +        ' t(jxs) = ',e15.6,' t(jxe) = ',e15.6,/,
c    +        ' u(jxs) = ',e15.6,' u(jxe) = ',e15.6,/,
c    +        ' v(jxs) = ',e15.6,' v(jxe) = ',e15.6)
c
      do iy=iys,iye
      do ix=mxs,mxe
c
c ----------------- first guess for utau
c
      utau = .001
c
      t10xy(1) = -qstar(1)/utau*zody*vk74in
      tsfcc(1) = u_level1(ix,iy,3) - t10xy(1)
      vsfc2    = u_level1(ix,iy,1)**2 + u_level1(ix,iy,2)**2
      vsfc     = sqrt(vsfc2)
      windm    = ufree+vsfc
      utausv   = utau
      utau2    = utau*utau
      amonin   = -utau2*utau/(batagk*qstar(1))
      if(amonin.eq.0.) then
            write(6,5050) ix,iy,it,utau,amonin
 5050       format(' 5050, sr. suft2, trouble at ',/,
     +             ' ix = ',i6,'iy = ',i6,' it = ',i6,' utau = ',e15.6,
     +             ' amonin = ',e15.6)
            stop
      endif
c
c ---- for unstable, free convection pbl
c
      iter = 0
 100  continue
c
c ----------------- limit the min (-l/z) change to accmmodate stable flow
c
      zeta_i = amin1(amonin/z1,zeta_mn_i)
      zeta_a = 1.0/zeta_i
c
      if(ismlt .eq. 1) then
          call busngr(zeta_a,phim,phis,psim,psis)
      else
          call fzol(zeta_a,phim,phis,psim,psis)
      endif
      utau     = windm*vk/(zody-psim)
      thstar(1)=-qstar(1)/utau
      amonold  = amonin
      amonin   = utau*utau/(batagk*thstar(1))
      diff     = abs(amonin - amonold)
c      write(nprt,5656)iter,psim,utau,zeta,amonin,dmonin,diff
c 5656 format(' iter=',i4,' phm=',e10.3,' utau=',e10.3,
c     1      ' zeta=',e10.3,' l=',e10.3,' diff = ',e12.4)
      iter = iter+1
      if(iter.gt.10)go to 1000
      if(diff.gt.abs(tol*amonin)) go to 100
 1000 continue
c
 2000 continue
c
      if (utau.gt.10.) then
         write(6,232)utau,windm
  232    format(' stop because utau=',e15.6,' windm=',e15.6)
         stop 9999
      endif
      t10xy(1) = -qstar(1)/utau*vk74in*(zody-psis)
      t_grnd(ix,iy,1) = u_level1(ix,iy,3) - t10xy(1)
c
      zol = zeta_a
      hol = zol*zi/z1
      tep = u_level1(ix,iy,1)/windm
      if(tep.gt.1.)  tep = 1.0
      if(tep.lt.-1.) tep = -1.0
      thta  = acos(tep)
      utau2 = utau*utau
c     au13m=-utau2*cos(thta)
c     au23m=-utau2*sin(thta)*sign(1.,u_level1(ix,iy,2))
c     aut3m(1)= qstar(1)
c
      tau13m(ix,iy)   = -utau2*cos(thta)
      tau23m(ix,iy)   = -utau2*sin(thta)*sign(1.,u_level1(ix,iy,2))
      taut3m(ix,iy,1) = qstar(1)
c
c **** get surface value of c scalar, specified surface flux
c
c     dnom      = (zody-psis)*vk74in
c     thstar(2) = -qstar(2)/utau
c     tsfcc(2)  = u_level1(ix,iy,4) - dnom*thstar(2)
c     t_grnd(ix,iy,2)  = u_level1(ix,iy,4) - dnom*thstar(2)
c     t10xy(2)  = thstar(2)*dnom
c     taut3m(ix,iy,2)  = qstar(2)
c
c
c ------- end of x-y loops
c
      enddo
      enddo
c
      return
      end
      subroutine init
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
c
      real :: RHT,RHB,mod_Magnus

      pi   = 4.0*atan(1.0)
      pi2  = 2.0*pi
      bfac = 1.0
      if(ibuoy.eq.0) bfac = 0.
c
c -------------------- case specific data
c
      if(iocean .eq. 1) then
         t00     = 283.
         t00b    = 5000.0
         cp      = 4.20e03
         gcp     = grav/cp
         batag   = bfac*grav/t00b
c        fcor    = 0.0
         fcor    = 1.39e-04
         fcor_h  = 0.0
c        wtsfc(1)=0.00
c        wtsfc(1)=4.96e-07
         wtsfc(1)=1.190476e-06
         qstar(1)=wtsfc(1)
c        dtdzf(1)=0.000
         dtdzf(1)=0.2548
         dtjump  = 0.
         divgls  = 0.
         zo      = 0.0001
         zi      = -5.
c        izi     = (55*nnz)/100
c        izi     = nnz
         izi     = 55
         xl      = 50.
         yl      = 50.
         zl      = -20.
c
c ---------- if stretched grid specify location of first point
c
         zw1 = -0.5
      else
         gcp     = grav/Cpa
         batag   = bfac*grav/t00
         fcor_h  = 0.0

         !!!If doing DNS, BCs are set somewhat differently than when doing LES
         if (iDNS .eq. 1) then

         !Get BCs of temperature in Celsius for use in specific expressions
         topC=(Ttop(1)-273.15)
         botC=(Tbot(1)-273.15)  
         
         !Normalize the RH given in input file
         RHT = Ttop(2)/100.0 
         RHB = Tbot(2)/100.0

         !Convert RH given in input file into specific humidity
         Ttop(2) =  RHT*Mw/Ru/Ttop(1)*610.94*
     +               EXP(17.6257*topC/(topC+243.04))/rhoa
         Tbot(2) =  RHB*Mw/Ru/Tbot(1)*610.94*
     +               EXP(17.6257*botC/(botC+243.04))/rhoa

         wtsfc(1) = qstar(1)
         wtsfc(2) = qstar(2)

         else !(doing LES)
c
         !qstar specified in input file for isfc = 0
        
         wtsfc = qstar

         !tsfcc specified in input file for isfc = 1

         !Assuming tsfcc(2) is giving RH:

         RHB = surf_RH/100.0

         tsfcc(2) = RHB*Mw/Ru/tsfcc(1)*mod_Magnus(tsfcc(1))/rhoa

         if (myid==0) then
         write(*,*) 'SURFACE HUMIDITY: %RH, q = ',surf_RH,tsfcc(2)
         end if



         end if
c
c
         !Gradients of temperature and scalars above inversion (upper BC)
         dtdzf(1)=(285.4-284.6)/(127.0-93.5)
         dtdzf(2)=(0.00663-0.007992)/(117.7-93.5)
c
         dtjump  = 0.0
         divgls  = 0.0

  
c
c ---------- if stretched grid specify location of first point
c
         !Particle channel:
         izi     = (60*nnz)/100
      endif
c
      time  = 0.0
c 
c ---------- outermost coarse grid  indicies are bounds of grid
c
      izlow = 1
      izup  = nnz
      dz    = zl/nnz
      dzg   = abs(dz)
      if(l_root) write(6,4040) zl,nnz,dzg
c
c --------------- generate z grids for particular mesh from
c                 iz = 0,1,...,nnz+1; this allows indexing
c                 to array elements z(0), etc.
c
      zwstrt = 0.0
c
c ------------ if uniform vertical spacing then
c
      if(iz_space .eq. 0) then
c
c ------------ build z grid for w points
c
         do iz=0,nnz+1
            z(iz) = dz*float(iz) + zwstrt
         enddo
      else
        !call vgrid(zw1,zi,zl,nnz,z(0),l_root,l_debug)
        !call vgrid_channel(zw1,zi,zl,nnz,z(0),l_root,l_debug)
        call vgrid_channel_fstrm(zw1,zi,zl,nnz,z(0),l_root,l_debug)
      endif
c
      call get_dz
c
      if(l_root) then
         write(6,8002) zwstrt
         write(6,8003) (iz,z(iz),zz(iz),iz=0,nnz+1)
      endif
c
      nnzm1 = nnz-1
      dx    = xl/nnx
      dy    = yl/nny
      fnxy  = 1./float(nxy)
      dzdz  = dzw(1)*dzw(1)
      z1    = zz(1)
c
      c23  = 2.0/3.0
      dsl  = (dx*1.5*dy*1.5*abs(dzw(1)))**(1./3.)
      dslg = dsl
      cs   = 0.2
c
      vk     = 0.4
      batagk = batag*vk
      vkin   = 1./vk
      ttmean = 0.
      zody   = alog(abs(z1/zo))
      write(nprt, 9901) z1,zo,zody
 9901 format(' 9901 z1 = ',e15.6,' zo = ',e15.6,/,
     +       ' zody = ',e15.6)
      zodyin = 1./zody
      wstar  = abs(batag*zi*wtsfc(1))**(1./3.)
      if(ismlt .eq. 1) then
c
c ---- set constants for businger similarity functions
c
         vk74   = vk*0.74
         vk74in = 0.74/vk
         zody74 = zody*0.74
      else 
c
c ---- set constants for large similarity functions
c
        vk74    = vk
        vk74in  = 1.0/vk
        zody74  = zody
      endif
      ugal   = 0.0
c      ugal   = ugcont*0.5
c     ugcont = ugcont - ugal
      cdbtm  = vk*vk/zody/zody
      if(iocean .eq. 1) then
c ----------- set surface friction velocity here and in sr. sufto
c        utau = 4.29e-03
         utau = 7.00e-03
      else
         ufree = 0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)
c
c ---- note : new estimate for utau !!!
c
         utau  = vk*(ufree+ugcont)/zody
c        utau  = vk*(ufree)/zody
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
c
      if(l_root) then
         write(6,80)
         write(6,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +         ,cdbtm,ugcont
      endif
c
      if(l_debug) then
         write(nprt,80)
         write(nprt,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +         ,cdbtm,ugcont
      endif
c
      return
c ------------------------
   2  format(10x,' WT =',e12.4,',  U* =',e12.4,',  L =',e12.4,/,
     +       10x,' DTDZ FREE =',e12.4,',  ZODY=',e12.4,/,10x,
     +       ' ZO(BTM) =',e12.4,',  CDBTM=',e12.4,
     +       ',  UG = ',e12.4)
  80  format(///,' ***** SCRATCH RUN ***** ',//)
 4040 format(' zl = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 4043 format(' znest = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 8002 format(' zwstrt = ',e12.4)
 8003 format(' iz ',5x,' zw',5x,' zu ',5x,/,(i3,2e12.4))
      end
      subroutine vgrid(z1,zi,zl,nnz,z,l_root,ldebug)
c
      real z(0:nnz+1)
      logical l_root, l_debug
c
c ----------------- build grid up to zi first
c
      z_frst = z1
      z_cntr = zi*0.5
      n_pbl  = nnz/2
c     n_pbl  = (5*nnz)/8
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
 9000       format(' Cannot find stretching factor',/,
     +             ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 10
      if(l_root) write(6,9100) z_fac, z_cntr, z1, knt
 9100 format(' Stretching factor = ',e15.6,/,
     +       ' Match point       = ',e15.6,/,
     +       ' First z           = ',e15.6,/,
     +       ' Number of iters   = ',i4)
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
c
      if(l_root) write(6,5300) n_pbl
 5300 format(' n_pbl = ',i4)
c
c -------------- build grid from zi on up
c
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
 8000       format(' Cannot find stretching factor',/,
     +             ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 20
      if(l_root) write(6,8100) z_fac, z_cntr, z1, knt
 8100 format(' Stretching factor = ',e15.6,/,
     +       ' Match point       = ',e15.6,/,
     +       ' First z           = ',e15.6,/,
     +       ' Number of iters   = ',i4)
c
      z(n_pbl+1) = zi + z_frst
      do iz=n_pbl+2,nnz-1
         z(iz) = zi + z_frst*
     +           (z_fac**(float(iz-n_pbl)) - 1.0)/(z_fac - 1.0)
      enddo
      z(nnz) = zl
      z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))
c     if(l_root) write(6,5600) (iz,z(iz),iz=0,nnz+1)
 5600 format(' 5600 in vgrid ',/,
     +       ' iz ',5x,' zw ',/,(i3,e15.6))
c
c     write(1,2000)
c2000 format('#k ',/,
c    +       '#lw 0.5 ',/,
c    +       '#m 1',/,
c    +       '#x 0 100 50',/,
c    +       '#y -50 2100 500')
c     x1 = 30.0
c     x2 = 80.0
c     do iz=0,nnz+1
c        write(1,1000) x1,z(iz)
c1000    format('#k ',/,
c    +          (2e15.6))
c        write(1,1100) x2,z(iz)
c1100    format(2e15.6)
c     enddo
c
      return
      end
      subroutine vgrid_channel(z1,zi,zl,nnz,z,l_root,ldebug)
c
      real z(0:nnz+1)
      integer :: zidx
      logical l_root, l_debug
c
c ----------------- build grid up to zi first
c
      z_frst = z1
      z_cntr = zi*0.5
      n_pbl  = nnz
c     n_pbl  = (5*nnz)/8
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
 9000       format(' Cannot find stretching factor',/,
     +             ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 10
      if(l_root) write(6,9100) z_fac, z_cntr, z1, knt
 9100 format(' Stretching factor = ',e15.6,/,
     +       ' Match point       = ',e15.6,/,
     +       ' First z           = ',e15.6,/,
     +       ' Number of iters   = ',i4)
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
c
      if(l_root) write(6,5300) n_pbl
 5300 format(' n_pbl = ',i4)
c
c -------------- build grid from zi on up
!     For the channel, zi represents the channel centerline
!     Want the mesh to be a mirror image across this:
c
!      zidx = 1
!      do iz=n_pbl+1,nnz
!         z(iz) = zi + (zi - z(n_pbl-zidx))
!         zidx = zidx + 1
!      enddo
      z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))
c
      return
      end subroutine vgrid_channel
      subroutine vgrid_channel_fstrm(z1,zi,zl,nnz,z,l_root,ldebug)
c
      real z(0:nnz+1)
      real s(0:2*nnz+1)
      integer :: zidx
      logical l_root, l_debug
c      
                nnz = 2.0*nnz
                zi = 2.0*zi
                zl = 2.0*zl
c
c ----------------- build grid up to zi first
c
      z_frst = z1
      z_cntr = zi*0.5
      n_pbl  = nnz
c     n_pbl  = (5*nnz)/8
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
 9000       format(' Cannot find stretching factor',/,
     +             ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 10
      if(l_root) write(6,9100) z_fac, z_cntr, z1, knt
 9100 format(' Stretching factor = ',e15.6,/,
     +       ' Match point       = ',e15.6,/,
     +       ' First z           = ',e15.6,/,
     +       ' Number of iters   = ',i4)
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
c
      if(l_root) write(6,5300) n_pbl
 5300 format(' n_pbl = ',i4)
c
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

      subroutine get_dz
c
c --------------- compute spacing for given vertical
c                 point distribution
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
      include 'mpif.h'
c
      do iz=1,nnz+1
         dzw(iz) = z(iz) - z(iz-1)
      enddo
      dzw(0)     = dzw(1)
      dzw(nnz+2) = dzw(nnz+1)
      do iz=0,nnz+2
         dzw_i(iz) = 1.0/dzw(iz)
      enddo
c
c ------------ build z grid for u points
c
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
c
      return
      end
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
     +     vyy(nnx,iys:iye,izs:izs),mod_Magnus,Ttmp

      real :: slope_q,slope_t
      real :: z_switch,t_switch,q_switch
      real :: t_surf,q_surf

      !Initialize partcount to 0:
      partcount = 0.0
      partcount_t = 0.0
      vpsum = 0.0
      vpsum_t = 0.0
      vpsqrsum = 0.0
      vpsqrsum_t = 0.0
      upwp_t = 0.0
      upwp = 0.0
      Tpsum = 0.0
      Tpsum_t = 0.0
      Tpsqrsum = 0.0
      Tpsqrsum_t = 0.0
      Tfsum = 0.0
      Tfsum_t = 0.0
      qfsum = 0.0
      qfsum_t = 0.0
      wpTpsum = 0.0
      wpTpsum_t = 0.0
      partsrc = 0.0
      partsrc_t = 0.0
      partTsrc = 0.0
      partTsrc_t = 0.0
      partHsrc = 0.0
      partHsrc_t = 0.0
      partTEsrc = 0.0
      partTEsrc_t = 0.0
      radsum = 0.0 
      radsum_t = 0.0 
      rad2sum = 0.0 
      rad2sum_t = 0.0 
      multcount = 0.0
      multcount_t = 0.0
      mwsum = 0.0
      mwsum_t = 0.0
      qstarsum = 0.0
      qstarsum_t = 0.0 
      u = 0.0
      v = 0.0
      w = 0.0
      e = 0.0
      t = 0.0
      r1 = 0.0
      r2 = 0.0
      r3 = 0.0
      r4 = 0.0
      r5 = 0.0

c ------------ note set nmatch in sr. iso so that
c              it is compatible with conditions here
c
      do iz=1,nnz
         ug(iz)   = ugcont
         vg(iz)   = vgcont
         divz(iz) = 0.0
      enddo

      !Parameters from Charlotte's C-FOG profile
      z_switch = 93.5
      t_switch = 284.5
      q_switch = 0.0079916
      t_surf = 284.0
      q_surf = 0.0081199

      !Slopes below the inversion
      slope_t = (t_switch-t_surf)/z_switch
      slope_q = (q_switch-q_surf)/z_switch
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = ugcont-ugal
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

         !Make a constant q profile at RH based on surface temp:
!            t(ix,iy,2,iz) = surf_RH/100.0*
!     +      Mw/Ru/t_lower*mod_Magnus(t_lower)/rhoa

         !Make a desired profile of RH:
!            !Convert theta to temperature:
!            Ttmp = 
!     +      t(ix,iy,1,iz)*(psurf/(psurf - zz(iz)*rhoa*grav))**(-Rd/Cpa)
!
!            !Get temperature based on desired RH
!            t(ix,iy,2,iz) = surf_RH/100.0*
!     +      Mw/Ru/Ttmp*mod_Magnus(Ttmp)/rhoa


!     Try to copy Charlotte's IC from the fog case (overwriting the
!     temperature profile from above)

!From the first version of the paper:
!            t(ix,iy,1,iz) = 286.0 + (1.0/128.0)*zz(iz)
!            t(ix,iy,2,iz) = 0.0093 - (0.0005/128.0)*zz(iz)
!
!            t(ix,iy,1,iz) = 286.0  !Remember this is the potential temp
!            !Actual temp
!            Ttmp = 
!     +      t(ix,iy,1,iz)*(psurf/(psurf - zz(iz)*rhoa*grav))**(-Rd/Cpa)
!
!            t(ix,iy,2,iz) = 99.0/100.0*
!     +      Mw/Ru/Ttmp*mod_Magnus(Ttmp)/rhoa


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
      subroutine random_f
c
c ---------- example of using given (sparse) initial 
c            sounding profiles (FIX for ncpu_s).
c            
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +     psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +     vyy(nnx,iys:iye,izs:izs)
c
      parameter (nt=12, nz=11)
      real zg(nz), u_i(nz,nt), v_i(nz,nt), theta_i(nz,nt)
      real ui_temp(nz), vi_temp(nz), ti_temp(nz)
      real time_g(nt)
c
      data time_g /
     +  0.0000E+00,  0.3600E+04,  0.7200E+04,  0.1080E+05,  0.1440E+05,
     +  0.1800E+05,  0.2160E+05,  0.2520E+05,  0.2880E+05,  0.3240E+05,
     +  0.3600E+05,  0.3960E+05
     +/
      data zg /
     +  0.1000E+02,  0.3000E+02,  0.5500E+02,  0.9000E+02,  0.1400E+03,
     +  0.2150E+03,  0.3300E+03,  0.5000E+03,  0.7500E+03,  0.1100E+04,
     +  0.1600E+04
     +/
      data u_i /
     + -0.1510E+01, -0.1560E+01, -0.1580E+01, -0.1580E+01, -0.1560E+01,
     + -0.1530E+01, -0.1510E+01, -0.9000E+00, -0.1390E+01, -0.1220E+01,
     + -0.5100E+00,
     + -0.1090E+01, -0.1110E+01, -0.1120E+01, -0.1120E+01, -0.1030E+01,
     + -0.9900E+00, -0.9500E+00, -0.6200E+00, -0.1230E+01, -0.9400E+00,
     +  0.2800E+00,
     + -0.9100E+00, -0.9200E+00, -0.9100E+00, -0.9000E+00, -0.8800E+00,
     + -0.8400E+00, -0.8000E+00, -0.6500E+00, -0.1510E+01, -0.1070E+01,
     +  0.2400E+00,
     + -0.8900E+00, -0.8900E+00, -0.8900E+00, -0.8800E+00, -0.8700E+00,
     + -0.8500E+00, -0.8100E+00, -0.7000E+00, -0.1830E+01, -0.8400E+00,
     +  0.3500E+00,
     + -0.1250E+01, -0.1260E+01, -0.1260E+01, -0.1250E+01, -0.1240E+01,
     + -0.1220E+01, -0.1160E+01, -0.8800E+00, -0.1980E+01, -0.1900E+00,
     +  0.7500E+00,
     + -0.1800E+01, -0.1810E+01, -0.1820E+01, -0.1820E+01, -0.1800E+01,
     + -0.1780E+01, -0.1710E+01, -0.1150E+01, -0.1960E+01,  0.3900E+00,
     +  0.9200E+00,
     + -0.2110E+01, -0.2130E+01, -0.2140E+01, -0.2140E+01, -0.2130E+01,
     + -0.2110E+01, -0.2050E+01, -0.9300E+00, -0.1400E+01,  0.8800E+00,
     +  0.9600E+00,
     + -0.2250E+01, -0.2280E+01, -0.2290E+01, -0.2300E+01, -0.2290E+01,
     + -0.2260E+01, -0.2070E+01, -0.4000E-01, -0.1600E+00,  0.1440E+01,
     +  0.1190E+01,
     + -0.2160E+01, -0.2200E+01, -0.2220E+01, -0.2220E+01, -0.2220E+01,
     + -0.2190E+01, -0.1610E+01,  0.1470E+01,  0.1420E+01,  0.2050E+01,
     +  0.1610E+01,
     + -0.2230E+01, -0.2270E+01, -0.2290E+01, -0.2300E+01, -0.2300E+01,
     + -0.2260E+01, -0.1350E+01,  0.2480E+01,  0.2380E+01,  0.2320E+01,
     +  0.1740E+01,
     + -0.1890E+01, -0.1930E+01, -0.1950E+01, -0.1950E+01, -0.1940E+01,
     + -0.1890E+01, -0.1120E+01,  0.3010E+01,  0.3030E+01,  0.2800E+01,
     +  0.2000E+01,
     + -0.1210E+01, -0.1230E+01, -0.1240E+01, -0.1230E+01, -0.1210E+01,
     + -0.1140E+01, -0.4600E+00,  0.3320E+01,  0.3510E+01,  0.3420E+01,
     +  0.2340E+01
     +/
      data v_i /
     +  0.4800E+00,  0.5100E+00,  0.5300E+00,  0.5700E+00,  0.6900E+00,
     +  0.7300E+00,  0.7600E+00,  0.1410E+01, -0.4200E+00, -0.3060E+01,
     + -0.3500E+01,
     +  0.7800E+00,  0.8100E+00,  0.8400E+00,  0.8900E+00,  0.1060E+01,
     +  0.1110E+01,  0.1130E+01,  0.1190E+01, -0.1040E+01, -0.2900E+01,
     + -0.3440E+01,
     +  0.3000E+00,  0.3200E+00,  0.3400E+00,  0.3800E+00,  0.4800E+00,
     +  0.5300E+00,  0.5800E+00,  0.5300E+00, -0.1330E+01, -0.2040E+01,
     + -0.2830E+01,
     + -0.2700E+00, -0.2600E+00, -0.2400E+00, -0.2200E+00, -0.1800E+00,
     + -0.1300E+00, -0.5000E-01,  0.1000E+00, -0.1170E+01, -0.1100E+01,
     + -0.2370E+01,
     + -0.5500E+00, -0.5400E+00, -0.5300E+00, -0.5100E+00, -0.4800E+00,
     + -0.4100E+00, -0.2600E+00,  0.1700E+00, -0.4200E+00, -0.2200E+00,
     + -0.2080E+01,
     + -0.2700E+00, -0.2600E+00, -0.2500E+00, -0.2400E+00, -0.2100E+00,
     + -0.1600E+00, -0.1000E-01,  0.8500E+00,  0.9700E+00,  0.3500E+00,
     + -0.2250E+01,
     +  0.5300E+00,  0.5400E+00,  0.5600E+00,  0.5700E+00,  0.6000E+00,
     +  0.6500E+00,  0.7600E+00,  0.1960E+01,  0.2280E+01,  0.3600E+00,
     + -0.2590E+01,
     +  0.1590E+01,  0.1630E+01,  0.1650E+01,  0.1680E+01,  0.1720E+01,
     +  0.1780E+01,  0.2010E+01,  0.3260E+01,  0.3110E+01,  0.1600E+00,
     + -0.2580E+01,
     +  0.2560E+01,  0.2620E+01,  0.2660E+01,  0.2690E+01,  0.2740E+01,
     +  0.2830E+01,  0.3400E+01,  0.4030E+01,  0.3030E+01, -0.7000E-01,
     + -0.2320E+01,
     +  0.3500E+01,  0.3600E+01,  0.3650E+01,  0.3700E+01,  0.3750E+01,
     +  0.3860E+01,  0.4580E+01,  0.4100E+01,  0.2450E+01,  0.6000E-01,
     + -0.1770E+01,
     +  0.4500E+01,  0.4640E+01,  0.4700E+01,  0.4760E+01,  0.4830E+01,
     +  0.4930E+01,  0.5420E+01,  0.3960E+01,  0.2000E+01,  0.5000E+00,
     + -0.1150E+01,
     +  0.5290E+01,  0.5470E+01,  0.5550E+01,  0.5620E+01,  0.5690E+01,
     +  0.5790E+01,  0.6070E+01,  0.4000E+01,  0.1910E+01,  0.9700E+00,
     + -0.5600E+00
     +/
      data theta_i /
     +  0.2936E+03,  0.2936E+03,  0.2937E+03,  0.2937E+03,  0.2938E+03,
     +  0.2942E+03,  0.2948E+03,  0.2980E+03,  0.3027E+03,  0.3092E+03,
     +  0.3186E+03,
     +  0.2937E+03,  0.2937E+03,  0.2937E+03,  0.2938E+03,  0.2939E+03,
     +  0.2942E+03,  0.2946E+03,  0.2978E+03,  0.3024E+03,  0.3090E+03,
     +  0.3184E+03,
     +  0.2939E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,
     +  0.2941E+03,  0.2944E+03,  0.2976E+03,  0.3023E+03,  0.3089E+03,
     +  0.3182E+03,
     +  0.2940E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,
     +  0.2941E+03,  0.2943E+03,  0.2975E+03,  0.3022E+03,  0.3087E+03,
     +  0.3181E+03,
     +  0.2940E+03,  0.2940E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,
     +  0.2940E+03,  0.2942E+03,  0.2974E+03,  0.3021E+03,  0.3086E+03,
     +  0.3180E+03,
     +  0.2941E+03,  0.2940E+03,  0.2940E+03,  0.2940E+03,  0.2941E+03,
     +  0.2941E+03,  0.2941E+03,  0.2973E+03,  0.3019E+03,  0.3085E+03,
     +  0.3179E+03,
     +  0.2942E+03,  0.2942E+03,  0.2942E+03,  0.2942E+03,  0.2941E+03,
     +  0.2941E+03,  0.2941E+03,  0.2973E+03,  0.3020E+03,  0.3086E+03,
     +  0.3179E+03,
     +  0.2943E+03,  0.2943E+03,  0.2943E+03,  0.2943E+03,  0.2943E+03,
     +  0.2943E+03,  0.2943E+03,  0.2975E+03,  0.3022E+03,  0.3087E+03,
     +  0.3181E+03,
     +  0.2945E+03,  0.2945E+03,  0.2945E+03,  0.2945E+03,  0.2945E+03,
     +  0.2944E+03,  0.2946E+03,  0.2978E+03,  0.3025E+03,  0.3090E+03,
     +  0.3184E+03,
     +  0.2947E+03,  0.2947E+03,  0.2947E+03,  0.2947E+03,  0.2946E+03,
     +  0.2946E+03,  0.2949E+03,  0.2980E+03,  0.3027E+03,  0.3093E+03,
     +  0.3187E+03,
     +  0.2949E+03,  0.2949E+03,  0.2949E+03,  0.2948E+03,  0.2948E+03,
     +  0.2948E+03,  0.2950E+03,  0.2982E+03,  0.3028E+03,  0.3094E+03,
     +  0.3188E+03,
     +  0.2950E+03,  0.2950E+03,  0.2950E+03,  0.2950E+03,  0.2950E+03,
     +  0.2950E+03,  0.2950E+03,  0.2982E+03,  0.3029E+03,  0.3095E+03,
     +  0.3188E+03
     +/
c
      save time_g, zg, u_i, v_i, theta_i
c
c --------- find time location of initial profiles 
c
      call lterp(nt,time_g,t_factor,jt,jtp1,t_weit)
c
      do iz=1,nz
         ui_temp(iz) = u_i(iz,jt)*(1.0 - t_weit) +
     +                 u_i(iz,jtp1)*t_weit
         vi_temp(iz) = v_i(iz,jt)*(1.0 - t_weit) +
     +                 v_i(iz,jtp1)*t_weit
         ti_temp(iz) = theta_i(iz,jt)*(1.0 - t_weit) +
     +                 theta_i(iz,jtp1)*t_weit
      enddo
c
c ----------- interpolate vertically
c
      do iz=izs,ize
         call lterp(nz,zg,zz(iz),kk,kkp1,weit)
         u_temp = ui_temp(kk)*(1.0 - weit) +
     +            ui_temp(kkp1)*weit
         v_temp = vi_temp(kk)*(1.0 - weit) +
     +            vi_temp(kkp1)*weit
         theta_temp = ti_temp(kk)*(1.0 - weit) +
     +            ti_temp(kkp1)*weit
c
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
c
c ------------- set initial random field to be
c               divergence free
c
      idum = -1
      do iz=izs,ize
         if (iz.le.8) then
c
c ----------- ampv and ampt are max amplitudes of random 
c             velocity and temperature fields
c
         ampv = 0.5
         ampt = 0.1
c  
c ------- simple random field scaled between 0 and 1
c
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = ran1(idum)
         enddo
         enddo
c
         do iy=iys,iye
         do ix=1,nnx
            psix(ix,iy) = psi(ix,iy)
            psiy(ix,iy,izs) = psi(ix,iy)
         enddo
         enddo
         call xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
c
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
c
c ---------- check divergence of initial field
c
         do iy=iys,iye
         do ix=1,nnx
            uxx(ix,iy)     = u(ix,iy,iz)
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
      subroutine randoc
c
c -------- random initial conditions for an
c          ocean simulation
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +     psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +     vyy(nnx,iys:iye,izs:izs)
c
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
c
c ------------- set initial random field to be
c               divergence free
c
      idum = -1
      do iz=izs,ize
      if (iz.le.4) then
c
c ----------- ampv and ampt are max amplitudes of random 
c             velocity and temperature fields
c
         ampv = 0.01
c        ampt = 0.00
         ampt = 0.0001
c  
c ------- simple random field scaled between 0 and 1
c
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = ran1(idum)
         enddo
         enddo
c
         do iy=iys,iye
         do ix=1,nnx
            psix(ix,iy) = psi(ix,iy)
            psiy(ix,iy,izs) = psi(ix,iy)
         enddo
         enddo
         call xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
c
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
      do iz=izs,ize
         ug(iz)=ugcont
         vg(iz)=vgcont
      enddo
c
      return
      end
      subroutine get_fields
c
c ----------- special routine to read just 3d fields
c             as an initial guess, easy to customize
c             if missing data, etc..
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
c
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
      real, allocatable, dimension(:,:,:) :: temp
      logical there
c
      allocate(temp(nvar,nnx,iys:iye))
c
c ---------- input file to read from
c
!      path_ran = 'XXXXXXXXX/u.le.cou000'
c
c --------------------- get restart file from local directory
c                       reuse unit number
c
      close(nvel)
c
      inquire(file=path_ran,exist=there)
      if(there) then
         if(l_root) write(6,6001) path_ran
      else
         if(l_root) write(6,6005) path_ran
         stop
      endif
c
c ---- open file
c
      call mpi_file_open(mpi_comm_world, path_ran,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, nvel, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ------------ read 3d fields, make rhs*8
c
      nsize  = int(nvar,k8)*nnx*nny
      nsize2 = int(nvar,k8)*nnx*(iys-1)
      n_read = nvar*nnx*(iye+1-iys)
c
      do k=izs,ize
         offset = int((k-1),k8)*nsize + nsize2
         call mpi_file_read_at_all(nvel,offset,temp,n_read,
     +                              mpi_real8,status,ierr)
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
c
      enddo
c
c ---- close file
c
      call mpi_file_close(nvel, ierr)
c
      deallocate(temp)
c
      do iz=izs,ize
c
         ug(iz) = 0.0
         vg(iz) = 0.0
c
c ---------------- initial guess for pressure
c
         do iy=iys,iye
         do ix=1,nnx
            p(ix,iy,iz) = 0.0
         enddo
         enddo
      enddo
c
      return
c ---------------------------- process errors
  100 continue
      write(6,9000) path_ran, nvel
      call mpi_finalize(ierr)
      stop
c
 9992 continue
      write(6,6100) nvel,iz
      call mpi_finalize(ierr)
      stop
c ---------------------
 6001 format(' SR. GET_FIELDS: FILE READ FOR INITIALIZATION = ',a80)
 6005 format(' 6005, SR. GET_FIELDS: cannot find restart file = ',a80)
 6100 format(' SR. GET_FIELDS: file read error on unit number = ',i2,/,
     +       '               at iz = ',i4)
 9000 format(' 9000, SR. GET_FIELDS: cannot open file =',a80,/,
     +       ' to unit number = ',i2)
      end
      subroutine forcing
c
c ----------- update surface temperature based on a 
c             constant cooling rate
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      tsfcc(1) = t_surf_i - c_rate*time
c
      return
      end
      subroutine set_sav(it,istart)
c
      use pars
      use fields
      use con_data
      use con_stats
c
      data ionce /0/
      save ionce
c
      if(it .ne. istart) then
c
c ------------------- increment time if not first time through
c
         time=time+dt
      endif
c
      it=it+1
c
      dt    = dt_new
      dtg   = dt
      mnout = (mod(it,imean).eq.0).or. (it.eq.1)
      mtape = (mod(it,itape).eq.0)
      micut = (mod(it,itcut).eq.0)
      mviz  = (mod(it,i_viz).eq.0)
      if(ihst .lt. 0) then
         mhis = .false.
      else
         mhis = (mod(it,ihst).eq.0 .and. it .ge. it_his)
      endif
       mtrans = (mod(it+1,ihst).eq.0)

      if (i_viz .lt. 0) then
         msave_v = .false.
      else
         msave_v = (mod(it,i_viz).eq.0 .and. it .ge. it_viz)
      endif
c
c ---------- decide whether velocity fields are saved
c
      msave = .false.
      if(it .ge. itstr .and. mtape) then
         itn=itn+1
         msave = .true.
         call get_output_filenames
      endif
c
c ---------- decide whether viz fields are saved
c
      !msave_v = .false.
      !if(it .ge. itstr .and. mviz .and. i_viz .gt. 0) then
      !   msave_v = .true.
      !   if(ionce .eq. 0) then
      !      ionce = 1
      !      call open_viz
      !   endif
      !endif
      if((i_viz .gt. 0) .and. (it .ge. it_viz_nxt)) then
        call viz_output_filename(it)
        it_viz_nxt = it_viz_nxt + itape
      endif
c
c --------- decide whether history files are to be saved
c
      if((ihst .gt. 0) .and. (it .ge. it_his_nxt)) then
         call open_his(it)
         it_his_nxt = it_his_nxt + itape
      endif
c
      return
      end
      subroutine pbltop(itop)
c
c ---------- get estimate of pbl top
c
c            method = 0, min of wt flux
c                        (good for buoyancy cases)
c            method = 1, uw flux less than critical value
c                        (good for ekman cases)
c            method = 2, running t average exceeds criterion
c                        (good for neutral cases with capping
c                         inversions)
c            method = 3, maximum gradient in temperature field
c                        (good for finding local zi see jas paper)
c                        with minimum search height (sr. setup)
c
c ------------ if method uses average statistics then only root
c              process need find zi
c
      use pars
      use fields
      use con_data
      use con_stats
      real trun(maxnz)
      include 'mpif.h'
      real gradloc(2,nnx,nny), gradmax(2,nnx,nny)
      external get_zi
c
      if(method .le. 2 .and. l_root) then
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = 1.0
      if (method .le. 0 .or. method .gt. 2) then
         itop=1
         wttot=wtle(1,1)+wtsb(1,1)
         wtmin=wttot*sgn
         do iz=2,nnz
            wttot=(wtle(iz,1)+wtsb(iz,1))*sgn
            if (wttot.le.wtmin) then
               itop=iz
               wtmin=wttot
            endif
         enddo
         zi=z(itop)
      else if (method .eq. 1) then
         itop = 1
         crit = 0.05
         uwsf = utau*utau
         do iz=1,nnzm1
               uwtot = (uwle(iz) + uwsb(iz))**2 +
     $                 (vwle(iz) + vwsb(iz))**2
               uwtot = sqrt(uwtot)
               if(uwtot/uwsf .gt. crit) then
                  itop=iz
               endif
         enddo
         zi=z(itop)
      else if (method .eq. 2) then
         trun(1) = txym(1,1)
         do iz=2,nnz
             weight = z(iz-1)/z(iz)
             trun(iz) = trun(iz-1)*weight + (1.0-weight)*txym(iz,1)
         enddo
         itop = 1
         tcrit = 0.25
         if(iocean .eq. 1) tcrit = 0.1
         do iz=2,nnz
                if(txym(iz,1) .gt. (trun(iz) + tcrit)) then
                  itop = iz
                  go to 320
                endif
         enddo
  320    continue
         zi=z(itop)
      endif
      do iy=1,nny
      do ix=1,nnx
         gradmax(2,ix,iy) = zi
      enddo
      enddo
c
c ----------- use gradient method, every process computes
c
      elseif(method .eq. 3) then
c
c ---------------- get local zi from gradient in temperaure field
c
c     dz_i = dzu_i(izs+1)
c     do iy=1,nny
c     do ix=1,nnx
c        gradloc(1,ix,iy) = (t(ix,iy,1,izs+1) - t(ix,iy,1,izs))*dz_i
c        gradloc(2,ix,iy) = z(izs)
c     enddo
c     enddo
c
c ------- similar to zeroing the stat array in sr. mean_stat
c
      do iy=1,nny
      do ix=1,nnx
         gradloc(1,ix,iy) = 0.0
         gradloc(2,ix,iy) = z(iz_min)
      enddo
      enddo
c
c ------------- now all z in this process
c
      if(iz_min .le. ize) then
      do iz=max(izs,iz_min),ize
         izp1 = iz + 1
         do iy=iys,iye
         do ix=1,nnx
            grad = (t(ix,iy,1,izp1) - t(ix,iy,1,iz))*dzu_i(izp1)
            if(grad .gt. gradloc(1,ix,iy)) then
               gradloc(1,ix,iy) = grad
               gradloc(2,ix,iy) = z(iz)
            endif
         enddo
         enddo
      enddo
      endif
c
c     call mpi_reduce(gradloc,gradmax,2*nnx*nny,mpi_real8,ziloc,
c    +                i_root,mpi_comm_world,ierror)
c
c ----------- alternate version using already defined function in mpi
c             passes 2 real8 variables
c
      call mpi_reduce(gradloc,gradmax,nnx*nny,mpi_2double_precision,
     +                mpi_maxloc,i_root,mpi_comm_world,ierror)
c
c ------------ get average on root process
c
      if(l_root) then
         zi_avg = 0.0
         do iy=1,nny
         do ix=1,nnx
            zi_avg = zi_avg + gradmax(2,ix,iy)
         enddo
         enddo
         zi = zi_avg*fnxy
c        itop = nint(zi/dz)
      endif
c
      endif
c
c -------- send average zi everywhere
c
      call mpi_bcast(zi,1,mpi_real8,
     +              i_root,mpi_comm_world,ierr)
c
      if(iocean .ne. 1) then
         do iz=1,nnz
            if(zi .ge. z(iz) .and.
     +         zi .lt. z(iz+1)) itop = iz
         enddo
      else
         do iz=1,nnz
            if(zi .le. z(iz) .and.
     +         zi .gt. z(iz+1)) itop = iz
         enddo
      endif
c
c     if(l_root) write(6,7001) myid,zi,itop
 7001 format(' 7001 in pbltop myid = ',i4,' zi = ',e15.6,
     +       ' itop = ',i3)
c
      return
      end
      subroutine get_zi(gradmax,gradout,len,itype)
c
      use pars
      real gradmax(*), gradout(*)
c
c     write(nprt,2001) myid, len
c2001 format(' 2001 in get_zi myid = ',i4,' len = ',i8)
c     write(nprt,2002) (i,gradmax(i),gradmax(i+1),i=1,len,2)
c2002 format(' i ',5x,' grad ',5x,' location ',/,
c    +      (i5,2e15.6))
c
      do i=1,len,2
         if(gradmax(i) .gt. gradout(i)) then
              gradout(i)   = gradmax(i)
              gradout(i+1) = gradmax(i+1)
         endif
      enddo
c
      return
      end
      subroutine print(lu,it,iz_strt,iz_end)
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      write(lu,4000)
 4000 format(30X,' --- SOLUTION ---')
      write(lu,4100) it,time,dt,zi,tsfcc(1),uwsfc,vwsfc,wtsfc(1),
     +              zol,hol,ucfl, vcfl, wcfl,
     +              t_ref
 4100 format(' IT=',I7,5x,'TIME (s) = ',e15.8,',  DT(s) = ',e15.6,/,
     +       10x,'ZTOP = ',e15.6,
     +       ',  TSFC = ',e15.6,
     +       ',  UW = ',e15.6,',  VW = ',e15.6,/,10x,
     +       'WT = ',e15.6,',  ZL =',e15.6,
     +       ',  HL = ',e15.6,/,10x,'U_cfl = ',e15.6,
     +       ',  V_cfl = ',e15.6,',  W_cfl = ',e15.6,/,10x,
     +       'Theta Ref = ',e15.6)
      write(lu,4200)
 4200 format(//,20x,'--------- HORIZONTAL MEAN VALUES ---------- ',
     +       //,2x,'IZ',4x,'T_MEAN',7x,'T2_MEAN',6x,
     + 'DIVG',8X,'LE_KE',6X,'SGS_KE',7X,'LE_WT',6X,
     + 'SGS_WT',7X,'SHRZ',8X,'BUOY')
      do 19999 iz=iz_end,iz_strt,-1
         !write(lu,4300)iz,txym(iz,1)-t_ref,divz(iz),
	 write(lu,4300)iz,txym(iz,1)-tref, txym(iz,2),divz(iz),
     +              englez(iz),eavg(iz),wtle(iz,1),
     +              wtsb(iz,1),shrz(iz),buyz(iz)
 4300    format(1X,I3,e12.4,8e12.4)
19999 continue
      write(lu,4400)tsfcc(1),wtsfc(1)
 4400 format('  SURFACE VALUE: TXYM=',F8.2,'               WTSB=',E9.2)
      if(iocean .eq. 1) then
         write(lu,4500) stokess,udrift,vdrift
 4500    format(/,' STOKESS = ',e12.4,' UDRIFT = ',e12.4,
     +          ' VDRIFT = ',e12.4)
      endif
      write(lu,4600) (iz,uxym(iz)+ugal,vxym(iz),uwle(iz),
     +       uwsb(iz),vwle(iz),vwsb(iz),iz=iz_strt,iz_end)
 4600 format(//,' IZ',5x,' UXYM + UGAL',8x,' VXYM',10x,' UWLE',10x,
     +          ' UWSB',10x,' VWLE',10x,' VWSB'
     +       ,/,(1x,i4,6(3x,e15.6)))
      if(ivis .eq. 1) then
         write(lu,4800) xksurf, nmatch, viscon, vise
 4800    format(//,' XKSURF = ',e15.6,' NMATCH = ',i4,/,
     +             ' VISCON = ',e15.6,' VISE = ',e15.6)
!         write(lu,4700) (iz,dfac(iz),iz=iz_strt,iz_end)
! 4700    format(//,'   IZ',5x,'  DFAC',/,(1x,i4,3x,e15.6))
      endif
c
c --------------- output additional scalars
c
      if(nscl .eq. 2) then
      write(lu,5005)tsfcc(2),wtsfc(2)
 5005 format(/,'  SURFACE VALUE: TXYM(2) =',e15.6,' WTSFC(2) = ',e15.6)
      write(lu,5100) (iz,txym(iz,2),wtle(iz,2),
     +              wtsb(iz,2),iz=iz_strt,iz_end)
 5100 format(//,' IZ',5x,' SCALAR-1 MEAN',8x,' WS1LE',10x,
     +          ' WS1SB',10x
     +       ,/,(1x,i4,3(3x,e15.6)))
      else if (nscl .eq. 3) then
!      write(lu,5205)tsfcc(2),wtsfc(2),tsfcc(3),wtsfc(3)
 5205 format(/,'  SURFACE VALUE: TXYM(2) =',e15.6,' WTSFC(2) = ',e15.6,
     +       /,'  SURFACE VALUE: TXYM(3) =',e15.6,' WTSFC(3) = ',e15.6)
      write(lu,5200) (iz,txym(iz,2),txym(iz,3),wtle(iz,2),
     +    wtsb(iz,2),wtle(iz,3),wtsb(iz,3),iz=iz_strt,iz_end)
 5200 format(//,' IZ',5x,' SCALAR-1 MEAN',8x,' SCALAR-2 MEAN',10x,
     +          ' WS1LE',10x,' WS1SB',10x,' WS2LE',10x,' WS1SB'
     +       ,/,(1x,i4,6(3x,e15.6)))
      endif
 
      return
      end
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
      parameter(js = 34, ns = 5, nstat = js + ns*nscl)   !Brian 7/30/14  js=25->26
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
         stat(iz,30) = stat(iz,30) + multcount(ix,iy,iz)
         stat(iz,31) = stat(iz,31) + mwsum(ix,iy,iz)
         stat(iz,32) = stat(iz,32) + partHsrc(ix,iy,iz)
         stat(iz,33) = stat(iz,33) + partTEsrc(ix,iy,iz)
         stat(iz,34) = stat(iz,34) + qstarsum(ix,iy,iz)
      enddo
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
      zconc(iz)  = stat(iz,11)/xl/yl/dzw(iz)
      vp1mean(iz) = stat(iz,12)/stat(iz,11)
      vp2mean(iz) = stat(iz,13)/stat(iz,11)
      vp3mean(iz) = stat(iz,14)/stat(iz,11)
      vp1msqr(iz) = sqrt(stat(iz,15)/stat(iz,11)-vp1mean(iz)**2)
      vp2msqr(iz) = sqrt(stat(iz,16)/stat(iz,11)-vp2mean(iz)**2)
      vp3msqr(iz) = sqrt(stat(iz,17)/stat(iz,11)-vp3mean(iz)**2)
      m1src(iz) = stat(iz,18)*fnxy
      m2src(iz) = stat(iz,19)*fnxy
      m3src(iz) = stat(iz,20)*fnxy
      uw_tot(iz) = uwle(iz) + uwsb(iz)
      vw_tot(iz) = vwle(iz) + vwsb(iz)
      upwpm(iz) = stat(iz,21)/stat(iz,11)-(vp1mean(iz)*vp3mean(iz))
      Tpmean(iz) = stat(iz,22)/stat(iz,11)
      Tpmsqr(iz) = sqrt(stat(iz,23)/stat(iz,11)-Tpmean(iz)**2)
      Tfmean(iz) = stat(iz,24)/stat(iz,11)
      qfmean(iz) = stat(iz,25)/stat(iz,11)
      wpTpm(iz) = stat(iz,26)/stat(iz,11)-(Tpmean(iz)*vp3mean(iz))
      Tpsrc(iz) = stat(iz,27)*fnxy
      radmean(iz) = stat(iz,28)/stat(iz,11) 
      rad2mean(iz)= stat(iz,29)/stat(iz,11)
 
      multmean(iz) = stat(iz,30)/stat(iz,11)
      mwmean(iz) = stat(iz,31)/stat(iz,11)
      Nc(iz) = stat(iz,30)/xl/yl/dzw(iz)
      ql(iz) = stat(iz,31)/xl/yl/dzw(iz)/rhoa
 
      Hpsrc(iz) = stat(iz,32)*fnxy
      TEpsrc(iz) = stat(iz,33)*fnxy
      qstarm(iz) = stat(iz,34)/stat(iz,11)

      !Doesn't need averaging if basing radiation source on average ql:
      radtend(iz) = radsrc(iz)

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
      subroutine tke_budget

!NOTE::: THIS HAS NOT BEEN UPDATED TO REFLECT HUMIDITY IN THE BUOYANCY
!EQUATION   DHR 7/5/16

c
c -------- get terms in resolved scale tke budget
c          as in gabls writeup at w-points
c          at istage = 1. 
c
      use pars
      use particles
      use fields
      use con_data
      use con_stats
c
      real :: stat(1:nnz,5)
      real :: s11s,s22s,s33s,s12s,s13s,s23s,wz,wzp
      real :: s13,s23,s33
      real :: ufluc,ufluc_t,ufluc_b,vfluc,vfluc_t,vfluc_b,wfluc
      real :: ffluc1,ffluc2,ffluc3
      real :: ffluc1p,ffluc2p,ffluc3p
      real :: weit,weit1
      integer :: iz,i,j,izp1,izm1
c
c -------- stat(.,1) = tke transport  = wq
c          stat(.,2) = pressure transport  = wp
c          stat(.,3) = tke dissipation
c          stat(.,4) = tke dissipation
c          stat(.,5) = particle force correlation
c
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
            s12s = weit1*(0.5*(uy(i,j,iz) + vx(i,j,iz)))**2 +
     +            weit*(0.5*(uy(i,j,izp1) + vx(i,j,izp1)))**2
            s13s = (0.5*((u(i,j,izp1) - u(i,j,iz) +
     +            uxym(iz) - uxym(izp1))*dzu_i(izp1) +
     +            wx(i,j,iz)))**2
            s23s = (0.5*((v(i,j,izp1) - v(i,j,iz) +
     +          vxym(iz) - vxym(izp1))*dzu_i(izp1) +
     +          wy(i,j,iz)))**2

         stat(iz,3) = stat(iz,3) + 2.0*vis_m(i,j,iz)*(
     +               s11s+s22s+s33s+2.0*(s12s+s13s+s23s))

            !Things for viscous transport - these are computed at u-locations since
            !a z-derivative is done to the average, which lands t_tau on w-locations
            
            ufluc_t   = u(i,j,izp1) - uxym(izp1)
            ufluc   = u(i,j,iz) - uxym(iz)
            ufluc_b   = u(i,j,izm1) - uxym(izm1)
            vfluc_t   = v(i,j,izp1) - vxym(izp1)
            vfluc   = v(i,j,iz) - vxym(iz)
            vfluc_b   = v(i,j,izm1) - vxym(izm1)
            wfluc = 0.5*( (w(i,j,iz)-wxym(iz)) 
     +                  + (w(i,j,izm1)-wxym(izm1)) )

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
         stat(iz,4) = stat(iz,4) + 2.0*vis_m(i,j,iz)*(
     +               ufluc*s13 + vfluc*s23 + wfluc*s33)
                    
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
         stat(iz,5) = stat(iz,5) +
     +       weit*(ufluc_t*ffluc1p)+weit1*(ufluc*ffluc1)+
     +       weit*(vfluc_t*ffluc2p) + weit1*(vfluc*ffluc2) + 
     +       (w(i,j,iz)-wxym(iz))*(weit*ffluc3p+weit1*ffluc3)
         enddo
         enddo
         stat(iz,3) = stat(iz,3)*fnxy
         stat(iz,4) = stat(iz,4)*fnxy
         stat(iz,5) = stat(iz,5)*fnxy
      enddo
c
c --------------- get transport terms as vertical arrays
c
      do iz=izs,ize
c
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
c
c --------- get estimate of turbulent transport term
c
            ufluc   = u(ix,iy,iz) - uxym(iz)
            vfluc   = v(ix,iy,iz) - vxym(iz)
            wfluc   = w(ix,iy,iz) - wxym(iz)
            wfluc_l = w(ix,iy,izm1) - wxym(izm1)
            stat(iz,1)  = stat(iz,1) + 0.25*(wfluc + wfluc_l)*
     +             (ufluc**2 + vfluc**2 + 0.5*(wfluc_l**2 + wfluc**2))
c
c --------- get estimate of pressure transport term
c
            pfluc = p(ix,iy,iz) - pxym(iz)
     +           -0.5*((u(ix,iy,iz)+stokes(iz))**2 + 
     +                 v(ix,iy,iz)*v(ix,iy,iz) + 
     +      0.5*(w(ix,iy,iz)*w(ix,iy,iz)+w(ix,iy,izm1)*w(ix,iy,izm1)))
            stat(iz,2) = stat(iz,2) + pfluc*0.5*(wfluc_l + wfluc)
         enddo
         enddo
         stat(iz,1) = stat(iz,1)*fnxy
         stat(iz,2) = stat(iz,2)*fnxy
      enddo
      call mpi_sum_z(stat(1,1),i_root,myid,nnz*5,1)
c
c ------ we have all terms on all processors for all z, add them up
c        treat tr_tau at bottom special, tr_tau = 0
c
!      tr_tau(0) = 0.0
      do iz=1,nnz
c
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
c
c ------------- gather all the budget terms
c
         t_tran(iz)  = t_wq(iz) + t_wp(iz) + t_tau(iz)
         t_rprod(iz) = -(dudz*uwle(iz) + dvdz*vwle(iz))
         !Old t_sprod had subgrid stuff
         !t_sprod(iz) =  (dudz*uwsb(iz) + dvdz*vwsb(iz))
         !Now make t_sprod the spray tke term to reduce new variables:
         t_sprod(iz) = -stat(iz,5)
         t_buoy(iz)  =  batag*wtle(iz,1)
         t_diss(iz) = stat(iz,3)
c
      enddo
c
      return
      end
      subroutine extra_flux_terms
      use pars
      use particles
      use fields
      use con_data
      use con_stats
      use fftwk
      implicit none

      real :: stat(1:nnz,3)
      real :: weit,weit1
      real :: Tflucp1,Tflucm1,Tfluc,Tmean
      real :: qflucp1,qflucm1,qfluc,qmean,dqpdz
      real :: Sq,Sqp1
      real :: gradTp(3),dTpdz1,dTpdz,dTdz
      integer :: iz,i,j,izp1,izm1,iscl

!     Compute the "extra" enthalpy budget terms, located at w locations:
c -------- stat(.,1) = < q' T' w' > 
c          stat(.,2) = < T' Sq > where Sq is q source
c          stat(.,3) = Dv*< T' dq'/dz > 

       

      stat = 0.0
      do iz=izs,ize
         izp1 = iz + 1
         izm1 = iz - 1
         weit = dzw(iz)/(dzw(iz) + dzw(izp1))
         weit1 = 1.0 - weit

       do j=iys,iye
       do i=1,nnx

         if (iz==1)  then
           Tmean = 2.0*Tbot(1) - txym(iz,1)
           Tflucm1 = t(i,j,1,izm1)-Tmean
           Tflucp1 = t(i,j,1,izp1)-txym(izp1,1)
         elseif (iz==nnz) then
           Tmean = 2.0*Ttop(1) - txym(iz,1)
           Tflucp1 = t(i,j,1,izp1)-Tmean
           Tflucm1 = t(i,j,1,izm1)-txym(izm1,1)
         else
           Tflucp1 = t(i,j,1,izp1)-txym(izp1,1)
           Tflucm1 = t(i,j,1,izm1)-txym(izm1,1)
         end if
         Tfluc = t(i,j,1,iz)-txym(iz,1)

         if (iz==1)  then
           qmean = 2.0*Tbot(2) - txym(iz,2)
           qflucm1 = t(i,j,2,izm1)-qmean
           qflucp1 = t(i,j,2,izp1)-txym(izp1,2)
           Sqp1 = -partHsrc(i,j,izp1)
         elseif (iz==nnz) then
           qmean = 2.0*Ttop(2) - txym(iz,2)
           qflucp1 = t(i,j,2,izp1)-qmean
           qflucm1 = t(i,j,2,izm1)-txym(izm1,2)
           Sqp1 = 0.0
         else
           qflucp1 = t(i,j,2,izp1)-txym(izp1,2)
           qflucm1 = t(i,j,2,izm1)-txym(izm1,2)
         end if
         qfluc = t(i,j,2,iz)-txym(iz,2)
         Sq = -partHsrc(i,j,iz)

         !Get dq'/dz at w-locations:
         dqpdz = (qflucp1 - qfluc)*dzu_i(izp1)

         !Then can get source #3:
         stat(iz,3) = stat(iz,3) + vis_s(i,j,2,iz)*
     +   (weit1*Tfluc + weit*Tflucp1)*dqpdz


         !Now source #1, which must be evaluated at w-location:

         stat(iz,1) = stat(iz,1) + w(i,j,iz)*
     +   (weit1*Tfluc + weit*Tflucp1)*
     +   (weit1*qfluc + weit*qflucp1)


         !Now source #2:
         stat(iz,2) = stat(iz,2) + 
     +   (weit1*Tfluc + weit*Tflucp1)*
     +   (weit1*Sq + weit*Sqp1)

     
       end do
       end do

         stat(iz,1) = stat(iz,1)*fnxy
         stat(iz,2) = stat(iz,2)*fnxy
         stat(iz,3) = stat(iz,3)*fnxy
      end do
         

      call mpi_sum_z(stat(1,1),i_root,myid,nnz*3,1)


c
c ------ we have all terms on all processors for all z, add them up
c
      do iz=1,nnz
c
c ------------- gather all the budget terms
c
         trip(iz) = stat(iz,1)
         TpSq(iz) = stat(iz,2)
         Tpdqp(iz) = stat(iz,3)
      enddo
c
      return
      end
      subroutine Tvar_budget(iscl)
      use pars
      use particles
      use fields
      use con_data
      use con_stats
      use fftwk
      implicit none

      real :: stat(1:nnz,4)
      real :: weit,weit1
      real :: tx_tmp(nnx,iys:iye), ty(nnx,iys:iye,izs-1:ize+1)
      real :: tx(nnx,iys:iye,izs-1:ize+1)
      real :: trans(izs:ize)
      real :: gradTp(3),Tfluc,dTpdz1,dTpdz,dTdz
      real :: Tflucp1,Tflucm1,qfluc,qflucp1,Tmean
      integer :: iz,i,j,izp1,izm1,iscl

c -------- stat(.,1) = Transport: -del.[U<T'2> + <u T'2> - alpha*grad(T'2)]
c          stat(.,2) = Dissipation: -2*alpha <grad(T').grad(T')> 
c          stat(.,3) = Particle: <T' Q'>

      !Need the y gradient of temp:
      do iz=izs-1,ize+1
      do j=iys,iye
      do i=1,nnx
         ty(i,j,iz)  = t(i,j,iscl,iz)
         tx_tmp(i,j)  = t(i,j,iscl,iz)
      enddo
      enddo
      call xderivp(tx_tmp(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      tx(1:nnx,iys:iye,iz) = tx_tmp(1:nnx,iys:iye)
      enddo

      call yd_mpi(ty(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
       

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
         
         stat(iz,2) = stat(iz,2)  - 2.0*vis_s(i,j,iscl,iz)*
     +                (gradTp(1)**2+gradTp(2)**2+gradTp(3)**2)

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
            stat(iz,3) = stat(iz,3) + 
     +         weit*(qflucp1*Tflucp1) + weit1*(qfluc*Tfluc)
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
            stat(iz,3) = stat(iz,3) + 
     +         weit*(qflucp1*Tflucp1) + weit1*(qfluc*Tfluc)
          endif
          endif

          if (iscl .eq. 2 .and. iHcouple .eq. 1 ) then
              qfluc = partTEsrc(i,j,iz)-TEpsrc(iz)
              if (iz==nnz) then
              qflucp1 = 0.0
              else
              qflucp1 = partTEsrc(i,j,izp1)-TEpsrc(izp1)
              end if
              Tv_part2(iz) = Tv_part2(iz) + 
     +          weit*(qflucp1*Tflucp1) + weit1*(qfluc*Tfluc)
  
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
      call mpi_sum_z(Tv_part2(1),i_root,myid,nnz,1)
      endif


c
c ------ we have all terms on all processors for all z, add them up
c
      do iz=1,nnz
         izp1 = iz + 1
         izm1 = iz - 1
         if(iz .eq. nnz) then
            Tv_tran(iz,iscl) = 0.0
         else
            Tv_tran(iz,iscl) = -(stat(izp1,1)-stat(iz,1))*dzu_i(izp1)
         endif
         dTdz = (txym(izp1,iscl)-txym(iz,iscl))*dzu_i(izp1)
c
c ------------- gather all the budget terms
c
         Tv_prod(iz,iscl) = -2.0*wtle(iz,iscl)*dTdz
         Tv_diss(iz,iscl) = stat(iz,2)
         Tv_part1(iz,iscl) = -stat(iz,3)
         !Tv_part2 is already gathered
      enddo
c
      return
      end
      subroutine write_his(iloc)
c
c ----- write history file with global parameters
c       write tsfcc specially to preserve digits!
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
c
      divgmax = 0.0
      do iz=1,nnz
         divgmax = amax1(divgmax, divz(iz))
      enddo
c
      ziavg = zi
      holtop = hol
      wt_min = wtsb(iloc,1)
      wt_le  = wtle(iloc,1)
      krec = krec + 1
      mid = nnz/4
      write(nhis1,6000) time,dt,utau,ziavg,amonin,holtop,
     +         (tsfcc(1)-t_ref),uwsfc,vwsfc,divgmax, wt_min, wt_le,
     +         ucfl, vcfl, wcfl, wtsfc(1),
     +         ups(mid),vps(mid),wps(mid),tps(mid,1),
     +         uwle(mid),uwsb(mid),uw_tot(mid),
     +         vwle(mid),vwsb(mid),vw_tot(mid),
     +         wtle(mid,1),wtsb(mid,1),wt_tot(mid,1),
     +         englez(mid),eavg(mid), wabs,
     +         lwc,
     +         Rep_avg,radavg,radmin,radmax

 6000 format(37e17.8)
c
c -------------- write profile information
c
      call write_prof(nhisp,krec,isize,c_s%wwsb)
c
      return
      end
      subroutine write_prof(nhisp,krec,num,f)
      real f(num)
      real*4 f32(num)
c
c -------------- build special 32 bit arrays for profiles
c
      do i=1,num
         f32(i) = f(i)
      enddo
c
      write(nhisp,err=999,rec=krec) (f32(i),i=1,num)
c
      return
c --------------- errors
  999 continue
      write(6,9000) num,krec
 9000 format(' 9000, trouble in ',
     +       'SR. save_prof cannot write profile data ',/,
     +       ' num = ',i8, 'krec = ',i6)
      stop
      end
      subroutine close_his
c
c ---- close history files
c
      use pars
      logical there
c
c ---- root closes and checks the files
c
      close(nhis1)
      close(nhisp)
      close(nrad)
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
      inquire(file=path_sav_hist,exist=there)
      if(.not.there) then
         write(6,8200) path_sav_hist
         call mpi_finalize(ierr)
         stop
      endif
      write(6,7000) path_sav_h
      write(6,7100) path_sav_hp
      write(6,7200) path_sav_hist
c
      return
c -------------------- process write errors
 7000 format(' HISTORY DATA IS WRITTEN IN FILE  ',a80)
 7100 format(' PROFILE HISTORY DATA IS WRITTEN IN FILE  ',a80)
 7200 format(' RADHIST DATA IS WRITTEN IN FILE  ',a80)
 8000 format(' SR. SAVE_HIS: Trouble history file not in path =',a80)
 8100 format(' SR. SAVE_HIS: Trouble profile history file',
     +       ' not in path =',a80)
 8200 format(' SR. SAVE_RADPDF: Trouble rad histogram file =',a80)
      end
      subroutine dealias
c
c --------- wave cutoff filter using 2d fft
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real wve(nny,jxs:jxe,izs:ize)
      real wves(nnxp2,iys:iye,izs:ize)
c
c --------- sharp spectral cutoff, specific to current 2dfft
c
      ix_cut   = 2*int(float(nnx)/3.) + 3
      iy_cut_l = int(float(nny)/3.) + 2
      iy_cut_u = nnyp2 - iy_cut_l
c
c ---------- u-equation
c
      call fft2d_mpi(u(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(u(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ---------- v-equation
c
      call fft2d_mpi(v(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(v(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ---------- w-equation
c
      call fft2d_mpi(w(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(w(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ---------- e-equation
c
      call fft2d_mpi(e(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(e(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ------------- scalars, not stored in correct order
c
      do iscl=1,nscl
         do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            wves(ix,iy,iz) = t(ix,iy,iscl,iz)
         enddo
         enddo
         enddo
         call fft2d_mpi(wves(1,iys,izs),wve(1,jxs,izs),trigx(1,1),
     +           trigc,nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
         call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
         call fft2d_mpi(wves(1,iys,izs),wve(1,jxs,izs),trigx(1,1),
     +           trigc,nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
         do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,iscl,iz) = wves(ix,iy,iz)
         enddo
         enddo
         enddo
      enddo
c
      return
      end
      subroutine sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
c
c --------- sharp cutoff filter for field wve stored
c           in 2d-fft order
c
      use pars
      real wve(nny,jxs:jxe,izs:ize)
c
      do iz=izs,ize
         do ix=jxs,jxe
         do iy=iy_cut_l,iy_cut_u
            wve(iy,ix,iz) = 0.0
         enddo
         enddo
      enddo
c
      if(jxe .lt. ix_cut) go to 999
c
         do iz=izs,ize
            do ix=max(jxs,ix_cut),jxe
            do iy=1,nny
               wve(iy,ix,iz) = 0.0
            enddo
            enddo
         enddo
c
  999 continue
c
      return
      end
      subroutine gridd
c
c ----------- allocate space and pass arrays using modules
c
      use pars
      use fields
      use particles 
      use fftwk
      use con_data
      use con_stats
c
c ------------- establish association between pointers
c               and data structures
c
      call fill_cc
      call fill_cs

      if (myid==0) write(6,5001) isize
 5001 format(' size of stats array = ',i8)

c
c
c ---------------- debug for arrays
c
      big = -99.0e+300
c
c ---------------- setup grid
c
      nnx = maxnx
      nny = maxny
      nnz = maxnz
c     izs = 1
c     ize = nnz
c
c
c ----------- make sure problem and cpu's match
c
      maxp   = numprocs-1
      ncpu_z = numprocs/ncpu_s
      if(mod(numprocs,ncpu_s) .ne. 0 .or.
     +   ncpu_z .gt. nnz) then
         go to 999
      endif
      if(l_root) write(6, 1100) ncpu_s, ncpu_z, numprocs,
     +                          maxp
      write(nprt,1100) ncpu_s, ncpu_z, numprocs, maxp
 1100 format(' Number of x-y slab cpus = ',i5,/,
     +       ' Number of z-level cpus  = ',i5,/,
     +       ' Total number of cpus    = ',i5,/,
     +       ' Max-p for index arrays  = ',i5)
c
c ---------------- allocate arrays for (i,j,k)-indexing on
c                  each processor (see set_range)
c
      allocate(ix_s(0:maxp), ix_e(0:maxp),
     +         jx_s(0:maxp), jx_e(0:maxp),
     +         kx_s(0:maxp), kx_e(0:maxp),
     +         mx_s(0:maxp), mx_e(0:maxp),
     +         iy_s(0:maxp), iy_e(0:maxp),
     +         jy_s(0:maxp), jy_e(0:maxp),
     +         is_s(0:maxp), is_e(0:maxp),
     +         iz_s(0:maxp), iz_e(0:maxp))
c
c ---------------- setup array sizes and variable dimensions
c
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
c
      write(nprt,7001) nnx,nny,nnz
 7001 format(' 7001 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)
c
      call set_range
c
      write(nprt,7002) nnx,nny,nnz
 7002 format(' 7002 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)
c
      num_y = iye + 1 - iys
c
c ------------- allocate solution arrays
c               account for nnxp2 for fields but not in rhs
c               and possible monotone for scalars
c
      allocate(u(nnxp2,iys:iye,izs-1:ize+1), 
     +         v(nnxp2,iys:iye,izs-1:ize+1), 
     +         w(nnxp2,iys:iye,izs-1:ize+1), 
     +         t(nnxp2,iys:iye,nscl,izs-2:ize+2), 
     +         e(nnxp2,iys:iye,izs-1:ize+1), 
     +         r1(nnx,iys:iye,izs-1:ize+1),
     +         r2(nnx,iys:iye,izs-1:ize+1),
     +         r3(nnx,iys:iye,izs-1:ize+1),
     +         r4(nnx,iys:iye,nscl,izs-1:ize+1),
     +         r5(nnx,iys:iye,izs-1:ize+1))

c   ------Allocate arrays for SFS stress and its derivatives
      allocate(sigm_s(nnxp2,iys:iye,izs-1:ize+1),
     +         sigm_sdx(nnxp2,iys:iye,izs-1:ize+1),
     +         sigm_sdy(nnxp2,iys:iye,izs-1:ize+1),
     +         sigm_sdz(nnxp2,iys:iye,izs-1:ize+1),
     +         vis_ss(nnxp2,iys:iye,izs-1:ize+1))

      allocate (sigm_sext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         sigm_sdxext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         sigm_sdyext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         sigm_sdzext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         vis_sext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3))

c
c ------------- allocate extended arrays for interpolation of
c               particle/spray location
c
      allocate(uext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3), 
     +         vext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3), 
     +         wext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         Text(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         T2ext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3)) 
      !Transposed velocities to do the uf interpolation:
      allocate(u_t(0:nnz+1,iys:iye,mxs:mxe),
     +         v_t(0:nnz+1,iys:iye,mxs:mxe),
     +         w_t(0:nnz+1,iys:iye,mxs:mxe), 
     +         T_t(0:nnz+1,iys:iye,mxs:mxe),
     +         T2_t(0:nnz+1,iys:iye,mxs:mxe))

      allocate(partcount(nnx,iys:iye,izs-1:ize+1))
      allocate(partcount_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(upwp_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(upwp(nnx,iys:iye,izs-1:ize+1))
      allocate(vpsum(nnx,iys:iye,izs-1:ize+1,1:3))
      allocate(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,1:3))
      allocate(vpsqrsum(nnx,iys:iye,izs-1:ize+1,1:3))
      allocate(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,1:3))
      allocate(Tpsum(nnx,iys:iye,izs-1:ize+1))
      allocate(Tpsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(Tpsqrsum(nnx,iys:iye,izs-1:ize+1))
      allocate(Tpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(Tfsum(nnx,iys:iye,izs-1:ize+1))
      allocate(Tfsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(qfsum(nnx,iys:iye,izs-1:ize+1))
      allocate(qfsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(wpTpsum(nnx,iys:iye,izs-1:ize+1))
      allocate(wpTpsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(partsrc(nnx,iys:iye,izs-1:ize+1,1:3))
      allocate(partsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1,1:3))
      allocate(partTsrc(nnx,iys:iye,izs-1:ize+1))
      allocate(partTsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1))
      allocate(radsum(nnx,iys:iye,izs-1:ize+1))
      allocate(radsum_t(0:nnz+1,iys:iye,mxs:mxe)) 
      allocate(rad2sum(nnx,iys:iye,izs-1:ize+1)) 
      allocate(rad2sum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(multcount(nnx,iys:iye,izs-1:ize+1)) 
      allocate(multcount_t(0:nnz+1,iys:iye,mxs:mxe)) 
      allocate(mwsum(nnx,iys:iye,izs-1:ize+1)) 
      allocate(mwsum_t(0:nnz+1,iys:iye,mxs:mxe)) 
      allocate(partHsrc(nnx,iys:iye,izs-1:ize+1))
      allocate(partHsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1))
      allocate(partTEsrc(nnx,iys:iye,izs-1:ize+1))
      allocate(partTEsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1))
      allocate(qstarsum(nnx,iys:iye,izs-1:ize+1)) 
      allocate(qstarsum_t(0:nnz+1,iys:iye,mxs:mxe)) 


c ------------- allocate space for boundary condition arrays
c               on top and bottom of domain
c
      allocate(ubc(nnx,iys:iye,2),
     +         vbc(nnx,iys:iye,2),
     +         wbc(nnx,iys:iye,2),
     +         tbc(nnx,iys:iye,nscl,2),
     +         ebc(nnx,iys:iye,2),
     +         pbc(nnx,iys:iye,2),
     +         pbc2(nnx,iys:iye,2))
c
c ------------ allocate space for wind and surface arrays
c
      allocate(wind(nnx,iys:iye), 
     +         tau13m(nnx,iys:iye), 
     +         tau23m(nnx,iys:iye), 
     +         taut3m(nnx,iys:iye,nscl), 
     +         t_grnd(nnx,iys:iye,nscl))
c
c ------------------- allocate space for derivative arrays
c
      allocate(ux(nnx,iys:iye,izs-1:ize+1),
     +         uy(nnx,iys:iye,izs-1:ize+1),
     +         vx(nnx,iys:iye,izs-1:ize+1),
     +         vy(nnx,iys:iye,izs-1:ize+1),
     +         wx(nnx,iys:iye,izs-1:ize+1),
     +         wy(nnx,iys:iye,izs-1:ize+1))
c
c ------------- allocate space for pressure, pressure bcs
c
      allocate(p(nnxp2,iys:iye,izs-1:ize+1),
     +         ptop(nnxp2,iys:iye,2))
c
c ------------- allocate space for viscosity and diffusivity
c
      allocate(vis_m(nnx,iys:iye,izs-1:ize+1),
     +         vis_s(nnx,iys:iye,nscl,izs-1:ize+1))
c
c ------------- allocate space for fft trig factors
c
      nq_trig = max(nnx,nny)
      allocate(trigx(2*nq_trig+15,2),
     +         trigc(4*nq_trig+15))
c
      return
  999 continue
c
      if(l_root) write(6,1000) numprocs, ncpu_s, mmz
      write(nprt,1000) numprocs, ncpu_s, nnz
 1000 format(' Gridd Trouble number of processors and grid',
     +          ' partitioning do not match!',/,
     +          ' Total num of cpus   = ',i5,
     +          ' Num cpu on x-y slab = ',i5,/,
     +          ' Num of z-levels     = ',i5)
      call mpi_finalize(ierr)
      end
      subroutine restart
c
c ----------- get restart file from local directory
c
      use pars
      use fields
      use con_data
      use con_stats
      character*80 path_res_c
      logical there
c
c --------------------- check if file is there
c
      inquire(file=path_res,exist=there)
      if(there) then
         if(l_root) write(6,6001) path_res
      else
         if(l_root) write(6,6005) path_res
         stop
      endif
c
c ------------------ get constant file
c
      iloc = index(path_res,' ')
      path_res_c = path_res(1:iloc-1)//'.con'
      inquire(file=path_res_c,exist=there)
      if(there) then
         if(l_root) write(6,6002) path_res_c
      else
         if(l_root) write(6,6006) path_res_c
         stop
      endif
      open(nvelc,err=200,file=path_res_c,form='unformatted',
     +        status='old')
c
      call read_res
c
      return
c ---------------------------- process errors
  100 continue
      write(6,9000) path_res, nvel
      call mpi_finalize(ierr)
      stop
c -----------------------
  200 continue
      write(6,9001) path_res_c, nvelc
      call mpi_finalize(ierr)
      stop
c -----------------------
 6001 format(' SR. RESTART: FILE READ = ',A80)
 6002 format(' SR. RESTART: CONSTANT FILE READ = ',A80)
 6005 format(' 6005, SR. RESTART: cannot find restart file = ',a80)
 6006 format(' 6005, SR. RESTART: cannot find constant file = ',a80)
 9000 format(' 9000, SR. RESTART: cannot open file =',a80,/,
     +       ' to unit number = ',i2)
 9001 format(' 9001, SR. RESTART: cannot open file =',a80,/,
     +       ' to unit number = ',i2)
      end
      subroutine read_res
c
c -------------- read restart file including constant file
c                changed for iys:iye
c
      use pars
      use fields
      use con_data
      use con_stats
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
c
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
      real, allocatable, dimension(:,:,:) :: temp
      allocate(temp(nvar,nnx,iys:iye))
c
c ---- open file
c
      call mpi_file_open(mpi_comm_world, path_res,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, nvel, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ------------ read 3d fields
c
      nsize  = int(nvar,k8)*nnx*nny
      nsize2 = int(nvar,k8)*nnx*(iys-1)
      n_read = nvar*nnx*(iye+1-iys)
c
      do k=izs,ize
         offset = int((k-1),k8)*nsize + nsize2
         call mpi_file_read_at_all(nvel,offset,temp,n_read,
     +                              mpi_real8,status,ierr)
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
c
      enddo
c
c ---- close file
c
      call mpi_file_close(nvel, ierr)
c
      deallocate(temp)
c
c ------------ every mpi process reads constant file
c
      rewind(nvelc)
      read(nvelc,err=9993) c_c, c_s, case
      close(nvelc)
c
      if(l_root) write(6,4001) case
 4001 format(' 4001, SR. RESTART: case from restart = ',a3)
c
c ----- special restart conditions -------------------------------------
c
c -------- set case name to case input
c
      case   = case_inp
      if(l_root) write(6,4002) case_inp, utau, utausv
 4002 format(' 4002, SR. RESTART:',/,
     +       ' files will be saved with case name = ',a3,/,
     +       ' utau = ',e15.6,' utausv = ',e15.6)
c
c ------------------- if new vis model set match point for
c                     outer grid
      nmatch = 48
      utau = utausv
c
c -------- hand coded changes to restart if needed
c
c       qstars = 0.000
c       wtsfcs = 0.000
c
c
c ------   reset qstar and wtsfc for no heat flux
c              qstar(1) = qstars
c              wtsfc(1) = wtsfcs
c              qstar(2) = qstars
c              wtsfc(2) = wtsfcs
c ------   redefine case id to input value
c              case = cases
c
      if(l_root) write(6,4012) time
      if(l_root) write(6,4013) qstar(1) , nmatch, case
c
      call get_dz
c
      return
c ------------------------  process errors from read
c9991 continue
c     write(6,6000) nvel,iz
c6000 format(' SR. READ_RES: hit end of file on unit number = ',i2,/,
c    +       '               at iz = ',i4)
c     call mpi_finalize(ierr)
c     stop
c ---------------------
 9992 continue
      write(6,6100) nvel,iz
 6100 format(' SR. READ_RES: error reading file on unit number = ',i2,/,
     +       '               at iz = ',i4)
      call mpi_finalize(ierr)
      stop
c ---------------------
 9993 continue
      write(6,6200) nvelc
 6200 format(' SR. READ_RES:',/,
     +       '    error reading constant file on unit number = ',i2)
      call mpi_finalize(ierr)
      stop
c ---------------------
 4012 format(' SR. RESTART: restart completed at T=',e15.6)
 4013 format('    after restart qstar = ',e15.6,' nmatch = ',i5,
     +       ' case = ',a3)
      end
      subroutine save_viz(it)
c
c --------------- save multiple (x-y), (x-z), (y-z), planes of data .
c                 modify recl in all open statements for more or less
c                 variables. 
c                 Constant - x, implies yz planes
c                 Constant - y, implies xz planes
c                 Constant - z, implies xy planes
c
c ------------- routine uses send/recv to get information in y-z planes
c
      use pars
      use fields
      use con_data
      use con_stats
      use fftwk
      use particles
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
c
      parameter(nvar_o = 6)
c
      integer ix_pick(maxnx),  iy_pick(maxny),  iz_pick(maxnz), 
     +        ix_order(maxnx), iy_order(maxny), iz_order(maxnz)
c
      integer istatus(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
c
      real(kind=4), dimension(nvar_o,nny,izs:ize) :: temp_x
      real(kind=4), dimension(nvar_o,nnx,izs:ize) :: temp_y
      real(kind=4), dimension(nvar_o,nnx,iys:iye) :: temp_z
      real, dimension(nvar_o,iys:iye,izs:ize)     :: buf_send
c
c ------------- don't touch
c
      data iviz_x,  iviz_y,  iviz_z  /0, 0, 0/
      data ionce_x, ionce_y, ionce_z, istuff /0, 0, 0, 0/
      data ix_pick, iy_pick, iz_pick /maxnx*0, maxny*0, maxnz*0/
      data ix_order, iy_order, iz_order /maxnx*0, maxny*0, maxnz*0/
      save iviz_x,  iviz_y,  iviz_z, 
     +     ix_pick, iy_pick, iz_pick, 
     +     ix_order, iy_order, iz_order,
     +     ionce_x, ionce_y, ionce_z, istuff,
     +     npln_x, npln_y, npln_z
c
c
c ----------- turn on z levels to save. Customize for your own use.
c             Set iz_pick(iz) = iz, ix_pick(ix) = ix, iy_pick(iy) = iy
c             Data is round-robin alternated in the data file for more than
c             1 plane for any particular view.
c
      iz_pick(12) = 12
      iz_pick(28) = 28
      iz_pick(64) = 64
      !iz_pick(20) = 20
      !iz_pick(45) = 45
      !iz_pick(60) = 60
c
c -------------- pick an x-z plane of data (can add more)
c
      iy_pick(nny/2) = nny/2
c     iy_pick(nny)   = nny
c
c -------------- pick a y-z plane of data (can add more)
c
      ix_pick(nnx/2) = nnx/2
c     ix_pick(nnx)   = nnx
c
c ------ find total number of z's turned on and open file once
c
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
            call mpi_file_open(mpi_comm_world, path_viz_xy,
     +                         mpi_mode_create+mpi_mode_rdwr,
     +                         mpi_info_null, nviz_z, ierr)
            disp = 0
            call mpi_file_set_view(nviz_z,disp,mpi_real4,mpi_real4,
     +                            'native',mpi_info_null,ierr)
         endif
      endif
c
c ------ find total number of y's turned on and open file once
c
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
            call mpi_file_open(mpi_comm_world, path_viz_xz,
     +                         mpi_mode_create+mpi_mode_rdwr,
     +                         mpi_info_null, nviz_y, ierr)
            disp = 0
            call mpi_file_set_view(nviz_y,disp,mpi_real4,mpi_real4,
     +                            'native',mpi_info_null,ierr)
         endif
      endif
c
c ------ find total number of x's turned on and open file once
c
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
            call mpi_file_open(mpi_comm_world, path_viz_yz,
     +                         mpi_mode_create+mpi_mode_rdwr,
     +                         mpi_info_null, nviz_x, ierr)
            disp = 0
            call mpi_file_set_view(nviz_x,disp,mpi_real4,mpi_real4,
     +                            'native',mpi_info_null,ierr)
         endif
      endif
c
      if(istuff .eq. 0 .and. l_root) then
         open(nviz_s,file=path_stuf)
         istuff = 1
      endif
c
c --------- write data, subtract t_ref to increase
c           resolution on 32 bit machines
c
c ---------- xy planes of data
c
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
c
c ---------- get just the fluctuating pressure field
c
               temp_z(5,i,j) = p(i,j,k) - pxym(k)
     +                        -(e(i,j,k) + e(i,j,km1))/3.0
     +                        -0.5*((u(i,j,k) + stokes(k))**2 +
     +                               v(i,j,k)*v(i,j,k) +
     +                         0.5*(w(i,j,k)*w(i,j,k) + 
     +                              w(i,j,km1)*w(i,j,km1)))

               temp_z(6,i,j) = u(i,j,k)-uxym(k)
               !temp_z(6,i,j) = partsrc(i,j,k,1)
               !temp_z(7,i,j) = partsrc(i,j,k,2)
               !temp_z(8,i,j) = partsrc(i,j,k,3)
               !temp_z(9,i,j) = partcount(i,j,k)/dx/dy/dzw(k)
            enddo
            enddo
#if defined(SWAP)
            call byteswap(temp_z)
#endif
            offset = int((iviz_z + iz_order(k) - 1),k8)*nsize + nsize2
            call mpi_file_write_at(nviz_z,offset,temp_z,n_write,
     +                             mpi_real4,istatus,ierr)
            if (ierr .ne. 0) go to 9991
         endif
      enddo
c
c ---------- xz planes of data
c
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
c
c ---------- get just the fluctuating pressure field
c
               temp_y(5,i,k) =  p(i,j,k) - pxym(k)
     +                          -(e(i,j,k)+e(i,j,km1))/3.0
     +                          -0.5*((u(i,j,k)+stokes(k))**2 +
     +                               v(i,j,k)*v(i,j,k) +
     +                           0.5*(w(i,j,k)*w(i,j,k) + 
     +                                w(i,j,km1)*w(i,j,km1)))

               temp_y(6,i,k) = u(i,j,k) - uxym(k)
               !temp_y(6,i,k) = partsrc(i,j,k,1)
               !temp_y(7,i,k) = partsrc(i,j,k,2)
               !temp_y(8,i,k) = partsrc(i,j,k,3)
               !temp_y(9,i,k) = partcount(i,j,k)/dx/dy/dzw(k)
            enddo
            enddo
#if defined(SWAP)
            call byteswap(temp_y)
#endif
            offset = int((iviz_y + iy_order(j) - 1),k8)*nsize + nsize2
            call mpi_file_write_at(nviz_y,offset,temp_y,nwrite,
     +                                mpi_real4,istatus,ierr)
            if (ierr .ne. 0) goto 9992
         endif
      enddo
c
c ---------- yz planes that cut across all processors
c            just have root node on that slab write data
c
      iviz_x  = iviz_x + npln_x
      n_write = nvar_o*nny*(ize+1-izs)
      nsize   = int(nvar_o,k8)*nny*nnz
      nsize2  = int(nvar_o,k8)*nny*(izs-1)
      n_send  = nvar_o*(ize+1-izs)*(iye+1-iys)
      do i=1,nnx
         if(ix_pick(i) .eq. i) then
c
c ----------- build send buffer
c
            do k=izs,ize
            km1 = k - 1
            do j=iys,iye
               buf_send(1,j,k) = u(i,j,k)
               buf_send(2,j,k) = v(i,j,k)
               buf_send(3,j,k) = w(i,j,k)
               buf_send(4,j,k) = (t(i,j,1,k) - t_ref)
c
c ---------- get just the fluctuating pressure field
c
               buf_send(5,j,k) = p(i,j,k) - pxym(k)
     +                          -(e(i,j,k) + e(i,j,km1))/3.0
     +                          -0.5*((u(i,j,k) + stokes(k))**2 +
     +                               v(i,j,k)*v(i,j,k) +
     +                           0.5*(w(i,j,k)*w(i,j,k) + 
     +                                w(i,j,km1)*w(i,j,km1)))
              buf_send(6,j,k) = u(i,j,k)-uxym(k)
              !buf_send(6,j,k) = partsrc(i,j,k,1)
              !buf_send(7,j,k) = partsrc(i,j,k,2)
              !buf_send(8,j,k) = partsrc(i,j,k,3)
              !buf_send(9,j,k) = partcount(i,j,k)/dx/dy/dzw(k)
            enddo
            enddo
            if(myid .ne. iss) then
              call mpi_send(buf_send(1,iys,izs),n_send,
     +                      mpi_real8,iss,1,
     +                      mpi_comm_world,ierr)
            else
              do k=izs,ize
              do j=iys,iye
              do ii=1,nvar_o
                 temp_x(ii,j,k) = buf_send(ii,j,k)
              enddo
              enddo
              enddo
              do l=iss+1,ise
                 call recv_yz_var(temp_x,nvar_o,nny,
     +                            iy_s(l),iy_e(l),izs,ize,l)
              enddo
#if defined(SWAP)
              call byteswap(temp_x)
#endif
              offset = int((iviz_x + ix_order(i) - 1),k8)*nsize + nsize2
              call mpi_file_write_at(nviz_x,offset,temp_x,n_write,
     +                          mpi_real4,istatus,ierr)
              if (ierr .ne. 0) goto 9993
            endif
         endif
      enddo
c
c ------------- ascii file with facts in it that goes
c               with visualization
c
      if(l_root) then
         write(nviz_s,5000) time, amonin, zi, utau
 5000    format(4e20.8)
      endif
c
c ---- last time step close the files
c
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
c
      return
c --------------------------  errors in writing viz file
 9991 continue
      write(6,6000) nviz_z, iz
 6000 format(' SR. SAVE_VIS:',/,
     +       '    trouble cannot write xy viz file on unit = ',i2,/,
     +       '             at iz = ',i4)
      call mpi_finalize(ierr)
      stop
c --------------------------  errors in writing viz file
 9992 continue
      write(6,6100) nviz_y, iz, iviz_y
 6100 format(' SR. SAVE_VIS:',/,
     +       '    trouble cannot write xz viz file on unit = ',i2,/,
     +       '             at iz = ',i4,/,
     +       '            iviz_y = ',i8)
c --------------------------  errors in writing viz file
 9993 continue
      write(6,6200) nviz_x, iz, iviz_x
 6200 format(' SR. SAVE_VIS:',/,
     +       '    trouble cannot write yz viz file on unit = ',i2,/,
     +       '             at iz = ',i4,/,
     +       '            iviz_x = ',i8)
      call mpi_finalize(ierr)
      stop
      end
      subroutine recv_yz_var(temp_x,nvar,nny,iys,iye,izs,ize,ir)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real buf(nvar,iys:iye,izs:ize)
      real(kind=4), dimension(nvar,nny,izs:ize) :: temp_x
c
      num = nvar*(ize+1-izs)*(iye+1-iys)
      call mpi_recv(buf(1,iys,izs),num,mpi_real8,ir,1,
     +             mpi_comm_world,istatus,ierr)
      do k=izs,ize
      do j=iys,iye
      do ii=1,nvar
         temp_x(ii,j,k) = buf(ii,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine save_v(it)
c
c --------------- save 3d fields
c
      use pars
      use fields
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
      logical there
c
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)                 nsize, nsize2
c
      real, allocatable, dimension(:,:,:) :: temp
      allocate(temp(nvar,nnx,iys:iye))
c
c ---- open file
c
      call mpi_file_open(mpi_comm_world, path_sav_v,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, nvel, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ---- write data
c
      nsize   = int(nvar,k8)*nnx*nny
      nsize2  = int(nvar,k8)*nnx*(iys-1)
      n_write = nvar*nnx*(iye+1-iys)
c
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
c

         offset = int((k-1),k8)*nsize + nsize2
c         call mpi_file_write_at(nvel,offset,temp,n_write,
c     +                              mpi_real8,status,ierr)
         call mpi_file_write_at_all(nvel,offset,temp,n_write,
     +                              mpi_real8,status,ierr)
         if (ierr /= 0) goto 9991
c
      enddo

c
c ---- close file
c
      call mpi_file_close(nvel, ierr)

c
c ---- check file
c
      if (l_root) then
         inquire(file=path_sav_v,exist=there)
         if(.not.there) then
            write(6,8000) nvel,myid
            call mpi_finalize(ierr)
            stop
         endif
         write(6,7000) it,path_sav_v
      endif
c
      deallocate(temp)
c
      return
c --------------------------  errors in writing restart file
 9991 continue
      write(6,6000) nvel, iz
 6000 format(' SR. SAVE_V:',/,
     +       '    trouble cannot write restart file on unit = ',i2,/,
     +       '             at iz = ',i4)
      call mpi_finalize(ierr)
      stop
c --------------------
 7000 format(' **** DATA SET AT IT = ',I6,/,
     +       '      VELOCITY DATA IS WRITTEN IN FILE  ',a80)
 8000 format(' in SAVE_V: trouble writing file ',i5,'  myid = ',i5,
     +       ' at iz = ',i5)
      end
      subroutine save_c(it)
c
c --------------- root process writes constant file
c                 sequential fortan binary
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      logical there
      character options*8, passwd*1
c
c ---- open file
c
      open(nvelc,err=9992,file=path_sav_c,form='unformatted',
     +                status='unknown')
      write(nvelc,err=9992) c_c, c_s, case
      close(nvelc)
c
        inquire(file=path_sav_c,exist=there)
        if(.not.there) then
           write(6,8001) path_sav_c
           call mpi_finalize(ierr)
           stop
        endif
c -----------------------------  output ok message
      write(6,7001) path_sav_c
c
      return
c --------------------------  errors in writing constant file
 9992 continue
      write(6,6100) nvelc
 6100 format(' SR. SAVE_V:',/,
     +  '    trouble cannot open/write constant file on unit = ',i2)
      call mpi_finalize(ierr)
      stop
c ---------------------
 7001 format('      CONSTANT DATA IS WRITTEN IN FILE  ',a80)
 8001 format(' SR. SAVE_C: Trouble constant file not in path =',a80)
      end
      subroutine save_p
c
c -------------- save pressure file
c
      use pars
      use fields
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
      logical there
c
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
c
      real temp(nnx,iys:iye)
c
c ---- open file
c
      call mpi_file_open(mpi_comm_world, path_sav_p,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, npre, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(npre,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ---- write data
c
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
         call mpi_file_write_at_all(npre,offset,temp,n_write,
     +                              mpi_real8,status,ierr)
c         call mpi_file_write_at(npre,offset,temp,n_write,
c     +                              mpi_real8,status,ierr)
      enddo
c
c ---- close file
c
      call mpi_file_close(npre, ierr)
c
c ---- check file
c
      if (l_root) then
         inquire(file=path_sav_p,exist=there)
         if(.not.there) then
            write(6,8000) path_sav_p
            call mpi_finalize(ierr)
            stop
         endif
         write(6,7000) path_sav_p
      endif
c
      return
c -------------------- process write errors
 9991 continue
      write(6,6000) npre, iz
 6000 format(' SR. SAVE_P:',/,
     +       '    trouble cannot write pressure file on unit = ',i2,/,
     +       '             at iz = ',i4)
      call mpi_finalize(ierr)
      stop
c -----------------------
 7000 format('      PRESSURE DATA IS WRITTEN IN FILE  ',a80)
 8000 format(' SR. SAVE_P: Trouble pressure file not in path =',a80)
      end
      subroutine get_units
c
      use pars
c
c -------------- unit numbers for files
c
      nvel  = 20 
      npre  = 30
      nhis1 = 40
      nrad = 45
      nvelc = 50
      nhisp = 60
      nviz_z = 80
      nviz_y = 82
      nviz_x = 84
      nviz_s = 90 
c
c ------------- unit number for standard print out
c               for each mpi task
c
      nprt = 1 
c
c ------------- open unit for standard printout
c
      path_prt = case_inp(1:3)//'.mp.xxxxx.out'
      write(path_prt(8:12),'(i5.5)') myid
      open(nprt,file=path_prt,form='formatted')
c
      return
      end
      subroutine get_output_filenames
c
c ----------- build file names for velocity, pressure, and constants
c
      use pars
      include 'mpif.h'
      character cgrid*10, num*3
c
c --------------- build character strings for file name
c
      cgrid = '.le.'
      write(num,'(i3.3)') itn
      iblnk = index(path_sav,' ')
      call blnk(path_sav_v)
      call blnk(path_sav_p)
      call blnk(path_sav_c)
      call blnk(path_sav_part)
      path_sav_v = path_sav(1:iblnk-1)//'/u'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)
      path_sav_p = path_sav(1:iblnk-1)//'/p'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)
      path_sav_c = path_sav(1:iblnk-1)//'/u'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)//'.con'
      path_sav_part = path_sav(1:iblnk-1)//'/part'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)
c
      return
      end
      subroutine open_his(istep)
c
c ------------------- open history files by root
c                     isize determined in sr. fill_cs
c
      use pars
      include 'mpif.h'
      character cgrid*4, iblks*16, num*3
      logical there
     
c
c --------------- build character strings for ascii history file name
c
      cgrid = '.le.'
      call blnk(iblks)
      write(iblks(1:7),'(i7.7)') istep
      iblks(8:8) = '_'
      write(iblks(9:15),'(i7.7)') (istep + itape)
      iblnk = index(path_his,' ')
      call blnk(path_sav_h)
      path_sav_h = path_his(1:iblnk-1)//'/his'//
     +         cgrid(1:4)//case(1:3)//'.'//iblks(1:15)
c
c --------------- build character strings for ieee profile history file
c                 set record counter for direct access file = 0
c
      krec = 0
      cgrid = '.le.'
      call blnk(iblks)
      write(iblks(1:7),'(i7.7)') istep
      iblks(8:8) = '_'
      write(iblks(9:15),'(i7.7)') (istep + itape)
      iblnk = index(path_his,' ')
      call blnk(path_sav_hp)
      path_sav_hp = path_his(1:iblnk-1)//'/his'//
     +         cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.ieee'

      !Do it again for radius pdf file:
      cgrid = '.le.'
      write(num,'(i3.3)') itn

      iblnk = index(path_histog,' ')
      call blnk(iblks)
      write(iblks(1:7),'(i7.7)') istep
      iblks(8:8) = '_'
      write(iblks(9:15),'(i7.7)') (istep+itape)
      call blnk(path_sav_hist)


      path_sav_hist = path_histog(1:iblnk-1)//'/radhist'//
     +          cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'      
c
c ----------------- save data in directory
c
      if(l_root) then

      close(nhis1)
      open(nhis1,err=3000,file=path_sav_h,form='formatted')
c
      close(nhisp)
      open(nhisp,err=4000,file=path_sav_hp,
     +        form='unformatted',access='direct',recl=isize*j_recl,
     +        status='unknown')

      close(nrad)
      open(nrad,err=4500,file=path_sav_hist,form='formatted')
      endif
c
      return
c ------------------- process errors
 3000 continue
      write(6,6301) nhis1, path_sav_h
 6301 format(' 6301, SR. OPEN_HIS:',/,
     +       '    cannot open history1 file on unit = ',i2,/,
     +       '    path = ',a80)
      stop
c-------------------
 4000 continue
      write(6,6302) nhisp, path_sav_hp
 6302 format(' 6302, SR. OPEN_HIS:',/,
     +       '    cannot open history profile file on unit = ',i2,/,
     +       '    path = ',a80)
      stop

 4500 continue
      write(6,6303) nrad, path_sav_hist
 6303 format(' 6303, SR. OPEN_RADPDF:',/,
     +       '    cannot open radhist file on unit = ',i2,/,
     +       '    path = ',a80)
      stop
      end
      subroutine viz_output_filename(istep)
c
c ------------------- set visualization files,
c                     leaves files in scratch directory 
c
      use pars
      include 'mpif.h'
      character iblks*16

c
c --------------- build character strings for file names
c                 with time step
c
      call blnk(iblks)
      iblks(1:1) = '.'
      write(iblks(2:8),'(i7.7)') istep
      iblks(9:9) = '_'
      write(iblks(10:16),'(i7.7)') (istep + itape)
c
      iloc = index(path_seed,' ')
      path_viz_xy = path_seed(1:iloc-1)
     +         //'/viz.'//case(1:3)//iblks(1:16)//'.xy.data'
c
      path_viz_xz = path_seed(1:iloc-1)
     +         //'/viz.'//case(1:3)//iblks(1:16)//'.xz.data'
c
      path_viz_yz = path_seed(1:iloc-1)
     +         //'/viz.'//case(1:3)//iblks(1:16)//'.yz.data'
c
      path_stuf = path_seed(1:iloc-1)
     +         //'/stuff.'//case(1:3)//iblks(1:16)//'.data'
c
c     if(l_root) then
c        write(6,8001) path_viz_xy
c8001    format(' 8001: viz file = ',a80)
c        write(6,8001) path_viz_xz
c        write(6,8001) path_viz_yz
c        write(6,8001) path_stuf
c        write(6,8001) path_seed
c     endif
c
      return
      end
      subroutine open_viz
c
c ------------------- set visualization files,
c                     leaves files in scratch directory 
c
      use pars
      include 'mpif.h'
      character iblks*16
c
c --------------- build character strings for file names
c                 with time step
c
      call blnk(iblks)
      iblks(1:1) = '.'
      write(iblks(2:8),'(i7.7)') iti
      iblks(9:9) = '_'
      write(iblks(10:16),'(i7.7)') itmax
c
      iloc = index(path_viz_xy,' ')
      path_viz_xy = path_viz_xy(1:iloc-1)
     +      //'/viz.'//case(1:3)//iblks(1:16)//'.xy.data'
      iloc = index(path_viz_xz,' ')
      path_viz_xz = path_viz_xz(1:iloc-1)
     +      //'/viz.'//case(1:3)//iblks(1:16)//'.xz.data'
      iloc = index(path_viz_yz,' ')
      path_viz_yz = path_viz_yz(1:iloc-1)
     +      //'/viz.'//case(1:3)//iblks(1:16)//'.yz.data'
      path_stuf = path_stuf(1:iloc-1)
     +      //'/stuff.'//case(1:3)//iblks(1:16)//'.data'
      close(nviz_z)
      close(nviz_y)
      close(nviz_x)
      close(nviz_s)
c
c ----------- do not actually open the files here since
c             not all processors may have been picked and
c             its unknown how many variables are selected.
c             customized in sr. save_viz 
c
      return
      end
      subroutine range(n1,n2,nprocs,irank,ista,iend)
c
c ---------- the ibm range finder to balance load
c
      iwork1 = (n2 - n1 + 1)/nprocs
      iwork2 = mod(n2 - n1 +1, nprocs)
      ista = irank*iwork1 + n1 + min(irank,iwork2)
      iend = ista + iwork1 - 1
      if(iwork2 .gt. irank) iend = iend + 1
c
      return
      end
      subroutine set_range
c
c ---- build special x,y,z-ranges. dimensioned for 0:numprocs-1
c      indexed with myid
c
c      [ix_s:ix_e] x-range for computing y-derivatives nx-pts/ncpu_s
c                  in xtoy and ytox tranposes
c
c      [jx_s:jx_e] x-range for computing 2d fft (nx+2)-pts/ncpu_s
c                  must be even in each x-interval for complex fft in y
c
c      [kx_s:kx_e] x-range for pressure solver transpose (nx+2)-pts/ncpu_z
c                  nx+2 fourier coefficients for xtoz and ztox transposes
c
c      [mx_s:mx_e] x-range split across z cpus as nx-pts/ncpu_z
c                  for use in surface layer routines 
c
c      [is_s:is_e] starting and ending processor id's for a
c                  particular z-level
c
c      [iy_s:iy_e] y-range for computing y-derivatives ny-pts/ncpu_s
c                  in xtoy and ytox tranposes
c
c      [jy_s:jy_e] y-range for use in xtoz and ztox transposes 
c                  in pressure solution
c
c      [iz_s:iz_e] z-range for a particular vertical slab
c
c
      use pars
c
      write(nprt,7002) nnx,nny,nnz
 7002 format(' 7002 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)
c
      ii = -1
      do nn=0,ncpu_z-1
         call range(1,nnx+2,ncpu_z,nn,lx_s,lx_e)
         call range(1,nnx,ncpu_z,nn,nx_s,nx_e)
         call range(1,nny,ncpu_z,nn,ly_s,ly_e)
         call range(1,nnz,ncpu_z,nn,mz_s,mz_e)
         do mm=0,ncpu_s-1
            call range(1,nny,ncpu_s,mm,ny_s,ny_e)
            call range(1,nnx,ncpu_s,mm,nxy_s,nxy_e)
            call range(1,ncx,ncpu_s,mm,l2x_s,l2x_e)
            ii       = ii + 1
c
            ix_s(ii) = nxy_s
            ix_e(ii) = nxy_e
            jx_s(ii) = (l2x_s - 1)*2 + 1
            jx_e(ii) = l2x_e*2
            kx_s(ii) = lx_s
            kx_e(ii) = lx_e
            mx_s(ii) = nx_s
            mx_e(ii) = nx_e
c
            iy_s(ii) = ny_s
            iy_e(ii) = ny_e
            jy_s(ii) = ly_s
            jy_e(ii) = ly_e
c
            iz_s(ii) = mz_s
            iz_e(ii) = mz_e
c
            is_s(ii) = (ii/ncpu_s)*ncpu_s
            is_e(ii) = is_s(ii) + ncpu_s - 1
         enddo
      enddo
c
      iys = iy_s(myid)
      iye = iy_e(myid)
      jys = jy_s(myid)
      jye = jy_e(myid)
      ixs = ix_s(myid)
      ixe = ix_e(myid)
      jxs = jx_s(myid)
      jxe = jx_e(myid)
      kxs = kx_s(myid)
      kxe = kx_e(myid)
      mxs = mx_s(myid)
      mxe = mx_e(myid)
      izs = iz_s(myid)
      ize = iz_e(myid)
c
c ----------- get starting and  ending processor id's on each
c             vertical slab
c
      iss = is_s(myid)
      ise = is_e(myid)
c
c ------------ debug ranges 
c
      if(l_debug) then
         write(nprt,1200) myid, (nn, ix_s(nn), ix_e(nn), jx_s(nn),
     +                     jx_e(nn), kx_s(nn), kx_e(nn),
     +                     nn = 0,numprocs-1)
 1200    format(' myid =  ',i4,/,
     +       ' nn',5x,' ixs ',5x,' ixe ',5x,' jxs ',5x,' jxe '
     +       ,5x,' kxs ',5x,' kxe',/,(7i6))
c
         write(nprt,1213) myid, (nn, iy_s(nn), iy_e(nn),
     +                   jy_s(nn), jy_e(nn),
     +                   iz_s(nn), iz_e(nn), is_s(nn), is_e(nn),
     +                   nn=0,numprocs-1)
 1213    format(' myid = ',i3,/,
     +       ' nn ',3x,' iys ',5x,' iye ',5x,
     +       ' jys ',5x,' jye ',5x,
     +       ' izs ',5x,' ize',5x,' iss ',5x,' ise ',/,
     +       (9i6))
      endif
c
      return
      end
      subroutine mpi_sum_xy(f,myid,iss,ise,nsend)
c
c --------- get horizontal x-y sum over a set of proccessors [iss:ise]
c           for vector f(i). f(i) is overwritten. skip if single processor
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real work(nsend,iss:ise), f(nsend)
c
      if(iss .eq. ise) go to 999
c
      do j=1,nsend
         work(j,myid) = f(j)
         f(j)         = 0.0
      enddo
      do i=iss,ise
         if(i .ne. myid) then
            call mpi_sendrecv(work(1,myid),nsend,mpi_real8,i,1,
     +               work(1,i),nsend,mpi_real8,i,1,
     +           mpi_comm_world,istatus,ierr)
         endif
      enddo
      do i=iss,ise
      do j=1,nsend
         f(j) = f(j) + work(j,i)
      enddo
      enddo
c
  999 continue
c
      return
      end
      subroutine mpi_sum_z(f,i_root,myid,nsend,iall)
c
c --------- get sums on root or all processors
c           for all z for vector f(i)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real recv_b(nsend), f(nsend)
c
c -------- just root gets the result
c
      if(iall .ne. 1) then
         call mpi_reduce(f(1),recv_b(1),nsend,mpi_real8,mpi_sum,i_root,
     +                  mpi_comm_world,ierr)
         if(myid .eq. i_root) then
            do i=1,nsend
               f(i) = recv_b(i)
            enddo
         endif
      else
c
c -------- everyone gets the result
c
         call mpi_allreduce(f(1),recv_b(1),nsend,mpi_real8,mpi_sum,
     +                  mpi_comm_world,ierr)
         do i=1,nsend
            f(i) = recv_b(i)
         enddo
      endif
c
      return
      end
      subroutine mpi_sum_z_s(f,i_root,myid,nsend,nscl,iall)
c
c --------- get sums on root or all processors
c           for all z for vector f(i,nscl)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real recv_b(nsend,nscl), f(nsend,nscl)
c
      if(iall .ne. 1) then
         call mpi_reduce(f(1,1),recv_b(1,1),nsend*nscl,mpi_real8,
     +        mpi_sum,i_root,mpi_comm_world,ierr)
         if(myid .eq. i_root) then
            do iscl=1,nscl
            do i=1,nsend
               f(i,iscl) = recv_b(i,iscl)
            enddo
            enddo
         endif
      else
         call mpi_allreduce(f(1,1),recv_b(1,1),nsend*nscl,mpi_real8,
     +        mpi_sum, mpi_comm_world,ierr)
         do iscl=1,nscl
         do i=1,nsend
            f(i,iscl) = recv_b(i,iscl)
         enddo
         enddo
      endif
c
      return
      end
      subroutine mpi_gath_root(fs,fr,iz_s,iz_e,izs,ize,nz,myid,np,ns)
c
c ---------- gather results on root processors
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      integer iz_s(0:np-1), iz_e(0:np-1)
      real fs(izs:ize), fr(1:nz)
c
      if(np .eq. 1) go to 999
c
      irow_r = mod(myid,ns)
      if(myid .gt. ns) then
        call mpi_send(fs(izs),ize+1-izs,mpi_real8,irow_r,1,
     +       mpi_comm_world,ierr)
      else
        do l=irow_r+ns,np-1,ns
           ind = iz_s(l) + 1
           num = iz_e(l) + 1 - iz_s(l)
           call mpi_recv(fr(ind),num,mpi_real8,l,1,
     +       mpi_comm_world,istatus,ierr)
        enddo
      endif
c
  999 continue
c
      return
      end
      subroutine mpi_send_root(fs,num,myid,np,ns)
c
c ---------- send root results to other processors above it
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real fs(num)
c
      if(np .eq. 1) go to 999
c
      irow_r = mod(myid,ns)
      if(myid .ge. ns) then
        call mpi_recv(fs(1),num,mpi_real8,irow_r,1,
     +       mpi_comm_world,istatus,ierr)
      else
        do l=irow_r+ns,np-1,ns
           call mpi_send(fs(1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
  999 continue
c
      return
      end
      subroutine xtoy_trans(f,g,nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,
     +           myid,ncpu_s,np)
c 
c ------- transpose array  f(nx,iys:iye,iz1:iz2) ---> g(ny,ixs:ixe,iz1:iz2)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,iz1:iz2), 
     +     g(ny,ixs:ixe,iz1:iz2)
      real ft(nx*(iye+1-iys)*(iz2 - iz1 + 1)),
     +     gt(ny*(ixe+1-ixs)*(iz2 - iz1 + 1))
      integer ix_s(0:np-1), ix_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      jk = (iye - iys + 1)*(iz2 - iz1 + 1)
      ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)
c
c ----------- get cpus on slab for myid
c
      islab = myid/ncpu_s
      iss   = islab*ncpu_s
      ise   = iss + ncpu_s - 1
c
      do i=iss,ise
         nsend = (ix_e(i) - ix_s(i) + 1)*jk
         nrecv = (iy_e(i) - iy_s(i) + 1)*ik
         if(i .eq. myid) then
            call send_xtoy(f,gt(1),nx,ix_s(i),ix_e(i),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)
         else
            call send_xtoy(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)
            call mpi_sendrecv(
     +          ft(1),nsend,mpi_real8,i,1,
     +          gt(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_xtoy(g,gt(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(i),iy_e(i),iz1,iz2)
      enddo
c
      return
      end
      subroutine send_xtoy(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
c
c ------------- grab correct chunk of array to be sent
c
      real f(nx,iys:iye,izs:ize), ft(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         ft(i,j,k) = f(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_xtoy(g,gt,ny,ixs,ixe,iys,iye,izs,ize)
      real g(ny,ixs:ixe,izs:ize), gt(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         g(j,i,k) = gt(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine ytox_trans(g,f,nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,
     +           myid,ncpu_s,np)
c 
c ------- transpose array g(ny,ixs:ixe,iz1:iz2) ---> f(nx,iys:iye,iz1:iz2)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,iz1:iz2), 
     +     g(ny,ixs:ixe,iz1:iz2)
      real ft(nx*(iye+1-iys)*(iz2 - iz1 + 1)),
     +     gt(ny*(ixe+1-ixs)*(iz2 - iz1 + 1))
c
      integer ix_s(0:np-1), ix_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      jk = (iye - iys + 1)*(iz2 - iz1 + 1)
      ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)
c
c ----------- get cpus on slab for myid
c
      islab = myid/ncpu_s
      iss   = islab*ncpu_s
      ise   = iss + ncpu_s - 1
      do i=iss,ise
         nsend = (iy_e(i) - iy_s(i) + 1)*ik
         nrecv = (ix_e(i) - ix_s(i) + 1)*jk
         if(i .eq. myid) then
            call send_ytox(g,ft(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(i),iy_e(i),iz1,iz2)
         else
            call send_ytox(g,gt(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(i),iy_e(i),iz1,iz2)
            call mpi_sendrecv(
     +          gt(1),nsend,mpi_real8,i,1,
     +          ft(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_ytox(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)
      enddo
c
      return
      end
      subroutine send_ytox(g,gt,ny,ixs,ixe,iys,iye,izs,ize)
c
c ------------- grab correct chunk of array to be sent
c
      real g(ny,ixs:ixe,izs:ize), gt(iys:iye,ixs:ixe,izs:ize)
c
      do k=izs,ize
      do i=ixs,ixe
      do j=iys,iye
         gt(j,i,k) = g(j,i,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_ytox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
      real f(nx,iys:iye,izs:ize), ft(iys:iye,ixs:ixe,izs:ize)
c
      do k=izs,ize
      do i=ixs,ixe
      do j=iys,iye
         f(i,j,k) = ft(j,i,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine xtoz_trans(f,g,nx,nz,ixs,ixe,ix_s,ix_e,
     +           iys,iye,izs,ize,iz_s,iz_e,
     +           myid,ncpu_s,numprocs)
c
c ------- transpose array  f(nx,iys:iye,izs-1:ize+1) 
c                     ---> g(0:nz+1,iys:iye,ixs:ixe)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,izs-1:ize+1), g(0:nz+1,iys:iye,ixs:ixe)
      real ft(nx*(iye+1-iys)*(ize-izs+1)),
     +     gt(nz*(ixe+1-ixs)*(iye-iys+1))
      integer ix_s(0:numprocs-1), ix_e(0:numprocs-1),
     +        iz_s(0:numprocs-1), iz_e(0:numprocs-1)
c
      jk = (ize - izs + 1)*(iye - iys + 1)
      ij = (ixe - ixs + 1)*(iye - iys + 1)
c
c ----------- get starting location
c
      iss = myid - (myid/ncpu_s)*ncpu_s
c
      do i=iss,numprocs-1,ncpu_s
         nsend = (ix_e(i) - ix_s(i) + 1)*jk
         nrecv = (iz_e(i) - iz_s(i) + 1)*ij
         if(i .eq. myid) then
            call send_xtoz(f,gt(1),nx,ix_s(i),ix_e(i),
     +                  iys,iye,iz_s(myid),iz_e(myid))
         else
            call send_xtoz(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iys,iye,iz_s(myid),iz_e(myid))
            call mpi_sendrecv(
     +          ft(1),nsend,mpi_real8,i,1,
     +          gt(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_xtoz(g,gt(1),nz,ix_s(myid),ix_e(myid),
     +                  iys,iye,iz_s(i),iz_e(i))
      enddo
c
      return
      end
      subroutine send_xtoz(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
c
c ------- grab correct chunk of array to be sent and skip ghost points
c
      real f(nx,iys:iye,izs-1:ize+1), ft(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         ft(i,j,k) = f(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_xtoz(g,gt,nz,ixs,ixe,iys,iye,izs,ize)
      real g(0:nz+1,iys:iye,ixs:ixe), gt(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         g(k,j,i) = gt(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine ztox_trans(g,f,nx,nz,ixs,ixe,ix_s,ix_e,
     +           iys,iye,izs,ize,iz_s,iz_e,
     +           myid,ncpu_s,numprocs)
c
c ------- transpose array g(0:nz+1,iys:iye,ixs:ixe) 
c                    ---> f(nx,iys:iye,izs-1:ize+1)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,izs-1:ize+1), g(0:nz+1,iys:iye,ixs:ixe)
      real ft(nx*(iye+1-iys)*(ize-izs+3)),
     +     gt((nz+3)*(iye+1-iys)*(ixe-ixs+1))
c
      integer ix_s(0:numprocs-1), ix_e(0:numprocs-1),
     +        iz_s(0:numprocs-1), iz_e(0:numprocs-1)
c
      jk = (ize - izs + 3)*(iye - iys + 1)
      ij = (ixe - ixs + 1)*(iye - iys + 1)
c
c ------------- get starting location
c
      iss = myid - (myid/ncpu_s)*ncpu_s
c
      do i=iss,numprocs-1,ncpu_s
         nsend = (iz_e(i) - iz_s(i) + 3)*ij
         nrecv = (ix_e(i) - ix_s(i) + 1)*jk
         if(i .eq. myid) then
            call send_ztox(g,ft(1),nz,ix_s(myid),ix_e(myid),
     +                  iys,iye,iz_s(i),iz_e(i))
         else
            call send_ztox(g,gt(1),nz,ix_s(myid),ix_e(myid),
     +                  iys,iye,iz_s(i),iz_e(i))
            call mpi_sendrecv(
     +          gt(1),nsend,mpi_real8,i,1,
     +          ft(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_ztox(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iys,iye,iz_s(myid),iz_e(myid))
      enddo
c
      return
      end
      subroutine send_ztox(g,gt,nz,ixs,ixe,iys,iye,izs,ize)
c
c ------------- grab correct chunk of array to be sent,
c               account for ghost points
c
      real g(0:nz+1,iys:iye,ixs:ixe), gt(izs-1:ize+1,iys:iye,ixs:ixe)
c
      do j=iys,iye
      do i=ixs,ixe
      do k=izs-1,ize+1
         gt(k,j,i) = g(k,j,i)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_ztox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
      real f(nx,iys:iye,izs-1:ize+1), ft(izs-1:ize+1,iys:iye,ixs:ixe)
c
      do i=ixs,ixe
      do j=iys,iye
      do k=izs-1,ize+1
         f(i,j,k) = ft(k,j,i)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine exchange
c
c ------------- exchange ghost points with mpi,
c               nb and nt are the destination and
c               source nodes. Allows for 1z per cpu
c
      use pars
      use fields
c     use fftwk
      include 'mpif.h'
c
      real fs(nnx,iys:iye,(4+nscl)),fr(nnx,iys:iye,(4+nscl))
      integer istatus(mpi_status_size)
c
      nb = myid - ncpu_s
      nt = myid + ncpu_s
c
c ------------ account for endpoints
c
      if(iss .eq. 0) then
         nb = mpi_proc_null
      endif
      if(ise .eq. numprocs-1) then
         nt = mpi_proc_null
      endif
      nsend = nnx*(iye + 1 - iys)*(4+nscl)
      nrecv = nsend
c
c --------- send top of myid, receive bottom from myid - ncpu_s
c
      do iy=iys,iye
      do ix=1,nnx
         fs(ix,iy,1) = u(ix,iy,ize)
         fs(ix,iy,2) = v(ix,iy,ize)
         fs(ix,iy,3) = w(ix,iy,ize)
         fs(ix,iy,4) = e(ix,iy,ize)
      enddo
      enddo
      do iscl=1,nscl
         jloc = 4 + iscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,jloc) = t(ix,iy,iscl,ize)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nt,0,
     +     fr(1,iys,1),nrecv,mpi_real8,nb,0,
     +     mpi_comm_world,istatus,ierr)
c
      if(iss .ne. 0) then
         izm1 = izs-1
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,izm1) = fr(ix,iy,1)
            v(ix,iy,izm1) = fr(ix,iy,2)
            w(ix,iy,izm1) = fr(ix,iy,3)
            e(ix,iy,izm1) = fr(ix,iy,4)
         enddo
         enddo
         do iscl=1,nscl
            jloc = 4 + iscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izm1) = fr(ix,iy,jloc)
            enddo
            enddo
         enddo
      endif
c
c -------- send bottom of myid, receive bottom from myid + ncpu_s
c
      do iy=iys,iye
      do ix=1,nnx
         fs(ix,iy,1) = u(ix,iy,izs)
         fs(ix,iy,2) = v(ix,iy,izs)
         fs(ix,iy,3) = w(ix,iy,izs)
         fs(ix,iy,4) = e(ix,iy,izs)
      enddo
      enddo
      do iscl=1,nscl
         jloc = 4 + iscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,jloc) = t(ix,iy,iscl,izs)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nb,1,
     +     fr(1,iys,1),nrecv,mpi_real8,nt,1,
     +     mpi_comm_world,istatus,ierr)
c
      if(ise .ne. numprocs-1) then
         izp1 = ize+1
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,izp1) = fr(ix,iy,1)
            v(ix,iy,izp1) = fr(ix,iy,2)
            w(ix,iy,izp1) = fr(ix,iy,3)
            e(ix,iy,izp1) = fr(ix,iy,4)
         enddo
         enddo
         do iscl=1,nscl
            jloc = 4 + iscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izp1) = fr(ix,iy,jloc)
            enddo
            enddo
         enddo
      endif
c
c --------------- send extra scalar points 
c
      nsend = nnx*(iye + 1 - iys)*nscl
      nrecv = nsend
c
c -------------- send top of myid, receive bottom from myid - ncpu_s
c
      izm1 = ize-1
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,iscl) = t(ix,iy,iscl,izm1)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nt,0,
     +     fr(1,iys,1),nrecv,mpi_real8,nb,0,
     +     mpi_comm_world,istatus,ierr)
c
      if(iss .ne. 0) then
         izm2 = izs-2
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izm2) = fr(ix,iy,iscl)
            enddo
            enddo
         enddo
      endif
c
c -------------- send bottom of myid, receive bottom from myid + ncpu_s
c
      izp1 = izs+1
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,iscl) = t(ix,iy,iscl,izp1)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nb,1,
     +     fr(1,iys,1),nrecv,mpi_real8,nt,1,
     +     mpi_comm_world,istatus,ierr)
c
      if(ise .ne. numprocs-1) then
         izp2 = ize+2
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izp2) = fr(ix,iy,iscl)
            enddo
            enddo
         enddo
      endif
c
      return
      end
      subroutine bcast_pbc
c
c ---- send upper boundary conditions to other processors 
c      for fft solution of pressure
c
      use pars
      use fields
      include 'mpif.h'
      integer istatus(mpi_status_size),ierr
c
      if(numprocs .eq. 1) go to 999
c
      irow_r = mod(myid,ncpu_s)
      irow_t = is_s(numprocs-1) + irow_r
      num = nnx*(iye+1-iys)
c
c
      if(iss .ne. is_s(numprocs-1)) then
c
c ------ not in the top row, receive from top
c
        call mpi_recv(pbc(1,iys,1),num,mpi_real8,irow_t,1,
     +       mpi_comm_world,istatus,ierr)
      else
c
c ------ myid is in the top row, send to everyone below
c
        do l=irow_r,irow_t-ncpu_s,ncpu_s
           call mpi_send(pbc(1,iys,1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
c --------- same thing for another variable
c
      if(iss .ne. is_s(numprocs-1)) then
c
c ------ not in the top row, receive from top
c
        call mpi_recv(pbc2(1,iys,1),num,mpi_real8,irow_t,1,
     +       mpi_comm_world,istatus,ierr)
      else
c
c ------ in the top row, send to everyone below
c
        do l=irow_r,irow_t-ncpu_s,ncpu_s
           call mpi_send(pbc2(1,iys,1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
  999 continue
c
      return
      end

      subroutine calc_radiation

      !Calculate a simple radiative temperature source based on Larson
      !et al. (2007) MWR, used in Mellado et al. (2018) JAMES
      use pars
      use particles
      use con_stats
      use con_data
      implicit none
      include 'mpif.h'

      real :: LWP_tmp(nnz-1),LWP(nnz-1)  !Have to calulate the integral from top down, then reverse ordering
      integer :: iz,lc

      !Do trapezoidal integration from the TOP down:
      !LWC only calculated on INTERIOR zw points

      lc = 1
      do iz=nnz-1,1,-1

         if (lc .eq. 1) then

         LWP_tmp(lc) = 0.5*rhoa*(ql(iz+1) + ql(iz))*(zz(iz+1)-zz(iz))

         else


         LWP_tmp(lc) = LWP_tmp(lc-1) +
     +   0.5*rhoa*(ql(iz+1) + ql(iz))*(zz(iz+1)-zz(iz))

         end if

         lc = lc + 1
      end do

      !Reverse the order so that LWP(1) is the bottom
      !LWP is stored at the interior zw points, since this is center of integral trapezoid
      do iz=1,nnz-1
         LWP(iz) = LWP_tmp(nnz-iz)
      end do

      !Now calculate the radiative flux based on LWP -- at zw points,
      !including top/bottom

      do iz=1,nnz-1

         radflux(iz) = rad_Fct*exp(-rad_kappa*LWP(iz)) +
     +         rad_Fcb*exp(-rad_kappa*(LWP(1) - LWP(iz)))


      end do


      !Give no radiative flux divergence at top/bottom
      radflux(0) = radflux(1)
      radflux(nnz) = radflux(nnz-1)

      !Now finally calculate the "radsrc" tendency term for the temperature equation
      !Based on the divergence of the radflux

      do iz=1,nnz

         radsrc(iz) = -1/rhoa/Cpa*(radflux(iz)-radflux(iz-1))/dzw(iz)


      end do


      end subroutine calc_radiation



      subroutine particle_update_rk3(it,istage)
      use pars
      use particles
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'

      integer :: istage,ierr,it
      real :: pi
      real :: denom,dtl,sigma
      integer :: ix,iy,iz
      real :: Rep,diff(3),diffnorm,corrfac,myRep_avg
      real :: xtmp(3),vtmp(3),Tptmp,radiustmp
      real :: Nup,Shp,rhop,taup_i,estar,einf
      real :: mylwc_sum,myphiw_sum,myphiv_sum,Volp      
      real :: Eff_C,Eff_S
      real :: t_s,t_f,t_s1,t_f1
      real :: mod_Magnus


      !First fill extended velocity field for interpolation
      !t_s = mpi_wtime()
      call fill_ext 
      !t_f = mpi_wtime()
      !call mpi_barrier(mpi_comm_world,ierr)
      !if (myid==5) write(*,*) 'time fill_ext:',t_f-t_s


      partcount_t = 0.0
      vpsum_t = 0.0
      upwp_t = 0.0
      vpsqrsum_t = 0.0
      Tpsum_t = 0.0
      Tfsum_t = 0.0
      qfsum_t = 0.0
      radsum_t = 0.0  
      rad2sum_t = 0.0  
      multcount_t = 0.0
      mwsum_t = 0.0
      Tpsqrsum_t = 0.0
      wpTpsum_t = 0.0
      myRep_avg = 0.0
      mylwc_sum = 0.0
      myphiw_sum = 0.0
      myphiv_sum = 0.0
      qstarsum_t = 0.0 

      t_s = mpi_wtime()


      !If you want, you can have the particles calculate nearest neighbor
      !Brute is there for checking, but WAY slower
      if (ineighbor) then
      !t_s = mpi_wtime()

      call particle_neighbor_search_kd
      !call particle_neighbor_search_brute

      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time neighbor:', t_f - t_s
      end if


      !Loop over the linked list of particles:
      part => first_particle
      do while (associated(part))     
         !First, interpolate to get the fluid velocity part%uf(1:3):
         if (ilin .eq. 1) then
            call uf_interp_lin   !Use trilinear interpolation
         else
            call uf_interp       !Use 6th order Lagrange interpolation
         end if

         if (iexner .eq. 1) then
             part%Tf = part%Tf*
     +  (psurf/(psurf-part%xp(3)*rhoa*grav))**(-Rd/Cpa)
         end if


         if (it .LE. 1 ) then 
            !part%xrhs(1:3) = part%vp(1:3)
            !part%xp(1:3) = xtmp(1:3) + dt*gama(istage)*part%xrhs(1:3)
            part%vp(1:3) = part%uf
            part%Tp = part%Tf
         endif

         !Now advance the particle and position via RK3 (same as velocity)
        
         !Intermediate Values
         pi = 4.0*atan(1.0)  
         diff(1:3) = part%vp - part%uf
         diffnorm = sqrt(diff(1)**2 + diff(2)**2 + diff(3)**2)
         Rep = 2.0*part%radius*diffnorm/nuf  
         Volp = pi2*2.0/3.0*part%radius**3
         rhop = (part%m_s+Volp*rhow)/Volp
         taup_i = 18.0*rhoa*nuf/rhop/(2.0*part%radius)**2 

         myRep_avg = myRep_avg + Rep
         corrfac = (1.0 + 0.15*Rep**(0.687))
         mylwc_sum = mylwc_sum + Volp*rhop*real(part%mult)
         myphiw_sum = myphiw_sum + Volp*rhow
         myphiv_sum = myphiv_sum + Volp


         !Compute Nusselt number for particle:
         !Ranz-Marshall relation
         Nup = 2.0 + 0.6*Rep**(1.0/2.0)*Pra**(1.0/3.0)
         Shp = 2.0 + 0.6*Rep**(1.0/2.0)*Sc**(1.0/3.0)


         !Mass Transfer calculations
         einf = mod_Magnus(part%Tf)
         Eff_C = 2.0*Mw*Gam/(Ru*rhow*part%radius*part%Tp)
         Eff_S = Ion*part%Os*part%m_s*Mw/Ms/(Volp*rhop-part%m_s)
         estar = einf*exp(Mw*Lv/Ru*(1.0/part%Tf-1.0/part%Tp)
     +                                           +Eff_C-Eff_S)
         part%qstar = Mw/Ru*estar/part%Tp/rhoa

  
         xtmp(1:3) = part%xp(1:3) + dt*zetas(istage)*part%xrhs(1:3)
         vtmp(1:3) = part%vp(1:3) + dt*zetas(istage)*part%vrhs(1:3) 
         Tptmp = part%Tp + dt*zetas(istage)*part%Tprhs_s
         Tptmp = Tptmp + dt*zetas(istage)*part%Tprhs_L
         radiustmp = part%radius + dt*zetas(istage)*part%radrhs

         part%xrhs(1:3) = part%vp(1:3)
         part%vrhs(1:3) = corrfac*taup_i*(part%uf(1:3)-part%vp(1:3))
     +                               - part_grav(1:3)

         if (ievap .EQ. 1) then      
            part%radrhs = Shp/9.0/Sc*rhop/rhow*part%radius*
     +                       taup_i*(part%qinf-part%qstar) !assumes qinf=rhov/rhoa rather than rhov/rhom
         else
            part%radrhs = 0.0
         end if


         part%Tprhs_s = -Nup/3.0/Pra*CpaCpp*
     +                     rhop/rhow*taup_i*(part%Tp-part%Tf)
         part%Tprhs_L = 3.0*Lv/Cpp/part%radius*part%radrhs



  
         part%xp(1:3) = xtmp(1:3) + dt*gama(istage)*part%xrhs(1:3)
         part%vp(1:3) = vtmp(1:3) + dt*gama(istage)*part%vrhs(1:3)
         part%Tp = Tptmp + dt*gama(istage)*part%Tprhs_s
         part%Tp = part%Tp + dt*gama(istage)*part%Tprhs_L
         part%radius = radiustmp + dt*gama(istage)*part%radrhs


         if (istage .eq. 3) part%res = part%res + dt


        part => part%next
      end do



      !t_f1 = mpi_wtime()
      !write(*,*) 'proc,loop time: ',myid,t_f1-t_s
      call mpi_barrier(mpi_comm_world,ierr)
      t_f = mpi_wtime()
      if (myid==5) write(*,*) 'time loop:', t_f-t_s

      !Enforce nonperiodic bcs (either elastic or destroying particles)
      !t_s = mpi_wtime()
      call particle_bcs_nonperiodic
      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time bc_non:', t_f - t_s

      !Check to see if particles left processor
      !If they did, remove from one list and add to another
      t_s = mpi_wtime()
      call particle_exchange
      call mpi_barrier(mpi_comm_world,ierr)
      t_f = mpi_wtime()
      if (myid==5) write(*,*) 'time exchg:', t_f - t_s

      !Now enforce periodic bcs 
      !just updates x,y locations if over xl,yl or under 0
      !t_s = mpi_wtime()
      call particle_bcs_periodic
      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time bc_per:', t_f - t_s

      part => first_particle
      do while (associated(part))     


      part => part%next
      end do

      !Now that particles are in their updated position, 
      !compute their contribution to the momentum coupling:
      !t_s = mpi_wtime()
      call particle_coupling_update
      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time cpl: ', t_f - t_s

      call particle_coupling_exchange

 
      !Finally, now that coupling and statistics arrays are filled, 
      !Transpose them back to align with the velocities:
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,1),
     +                partsrc(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,2),
     +                partsrc(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,3),
     +                partsrc(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(partTsrc_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partTsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(partHsrc_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partHsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(partTEsrc_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partTEsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(mwsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                mwsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      !Try only calling these when the history data is being written:
      if(mtrans  .and. istage .eq. 3) then
      call ztox_trans(upwp_t(0:nnz+1,iys:iye,mxs:mxe),
     +                upwp(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,1),
     +                vpsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,2),
     +                vpsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,3),
     +                vpsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,1),
     +                vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,2),
     +                vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,3),
     +                vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(Tpsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                Tpsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(Tpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                Tpsqrsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(Tfsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                Tfsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(qfsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                qfsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(wpTpsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                wpTpsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(partcount_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partcount(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(radsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                radsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(rad2sum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                rad2sum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs) 

      call ztox_trans(multcount_t(0:nnz+1,iys:iye,mxs:mxe),
     +                multcount(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)


      call ztox_trans(qstarsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                qstarsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      end if


      !t_s = mpi_wtime
      !Get particle count:
      numpart = 0
      part => first_particle
      do while (associated(part))
      numpart = numpart + 1
      part => part%next
      end do
      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time numpart: ', t_f - t_s
 
      !t_s = mpi_wtime()
      !Compute total number of particles
      call mpi_allreduce(numpart,tnumpart,1,mpi_integer,mpi_sum,
     +                   mpi_comm_world,ierr)
      !Compute average particle Reynolds number
      call mpi_allreduce(myRep_avg,Rep_avg,1,mpi_real8,mpi_sum,
     +                   mpi_comm_world,ierr)

      Rep_avg = Rep_avg/tnumpart

      call mpi_allreduce(mylwc_sum,lwc,1,mpi_real8,mpi_sum,
     +                   mpi_comm_world,ierr)

      call mpi_allreduce(myphiw_sum,phiw,1,mpi_real8,mpi_sum,
     +                   mpi_comm_world,ierr)

      call mpi_allreduce(myphiv_sum,phiv,1,mpi_real8,mpi_sum,
     +                   mpi_comm_world,ierr)

      
      phiw = phiw/xl/yl/zl/rhoa
      phiv = phiv/xl/yl/zl

      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time mpi_allreduce: ', t_f - t_s

      end subroutine particle_update_rk3
      subroutine particle_update_BE(it)
      use pars
      use particles
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'

      integer :: ierr,it
      real :: tmpbuf(5),tmpbuf_rec(5),tmp_reintro
      real :: myradavg,myradmax,myradmin,mytempmax,mytempmin
      real :: myqmin,myqmax
      real :: denom,dtl,sigma
      integer :: ix,iy,iz,im,flag,mflag
      real :: Rep,diff(3),diffnorm,corrfac,myRep_avg
      real :: xtmp(3),vtmp(3),Tptmp,radiustmp
      real :: Nup,Shp,rhop,taup_i,estar,einf
      real :: mylwc_sum,myphiw_sum,myphiv_sum,Volp
      real :: Eff_C,Eff_S
      real :: t_s,t_f,t_s1,t_f1
      real :: rt_start(2)
      real :: rt_zeroes(2)
      real :: taup0, dt_taup0, temp_r, temp_t, guess
      real :: mod_Magnus,tmp_coeff



      !First fill extended velocity field for interpolation
      call fill_ext

      partcount_t = 0.0
      vpsum_t = 0.0
      upwp_t = 0.0
      vpsqrsum_t = 0.0
      Tpsum_t = 0.0
      Tfsum_t = 0.0
      qfsum_t = 0.0
      radsum_t = 0.0
      rad2sum_t = 0.0
      multcount_t = 0.0
      mwsum_t = 0.0
      Tpsqrsum_t = 0.0
      wpTpsum_t = 0.0
      myRep_avg = 0.0
      mylwc_sum = 0.0
      myphiw_sum = 0.0
      myphiv_sum = 0.0
      qstarsum_t = 0.0

      partsrc_t = 0.0
      partTsrc_t = 0.0
      partHsrc_t = 0.0
      partTEsrc_t = 0.0

      denum = 0
      actnum = 0
      num100 = 0
      num1000 = 0
      numimpos = 0

      !loop over the linked list of particles
      part => first_particle
      do while (associated(part))


         !First, interpolate to get the fluid velocity part%uf(1:3):
         if (ilin .eq. 1) then
            call uf_interp_lin   !Use trilinear interpolation
         else
            call uf_interp       !Use 6th order Lagrange interpolation
         end if


        if (it .LE. 1) then
           part%vp(1:3) = part%uf
        end if

         if (iexner .eq. 1) then
             part%Tf = part%Tf*
     +  (psurf/(psurf-part%xp(3)*rhoa*grav))**(-Rd/Cpa)
         end if

        !Debug:
        part%sigm_s = part%radius

        diff(1:3) = part%vp - part%uf
        diffnorm = sqrt(diff(1)**2 + diff(2)**2 + diff(3)**2)
        Volp = pi2*2.0/3.0*part%radius**3
        rhop = (part%m_s+Volp*rhow)/Volp
        taup_i = 18.0*rhoa*nuf/rhop/(2.0*part%radius)**2
        Rep = 2.0*part%radius*diffnorm/nuf
        corrfac = (1.0 + 0.15*Rep**(0.687))

        !implicitly calculates next velocity and position
        part%xp(1:3) = part%xp(1:3) + dt*part%vp(1:3)
        part%vp(1:3) = (part%vp(1:3)+taup_i*dt*corrfac*part%uf(1:3)
     +   +dt*part_grav(1:3))/(1+dt*corrfac*taup_i)

        ! non-dimensionalizes particle radius and temperature before
        ! iteratively solving for next radius and temperature

        taup0 = (((part%m_s)/((2./3.)*pi2*radius_init**3) + rhow) *
     +   (radius_init*2)**2)/(18*rhoa*nuf)

        dt_taup0 = dt/taup0

        if (ievap .EQ. 1) then

               !Gives initial guess into nonlinear solver
               !mflag = 0, has equilibrium radius; mflag = 1, no
               !equilibrium (uses itself as initial guess)
               call rad_solver2(guess,mflag)

               if (mflag == 0) then
                rt_start(1) = guess/part%radius
                rt_start(2) = part%Tf/part%Tp
               else
                rt_start(1) = 1.0
                rt_start(2) = 1.0
               end if

               call gauss_newton_2d(part%vp,dt_taup0,
     +                           rt_start, rt_zeroes,flag)

               if (flag==1) then
               num100 = num100+1

               call LV_solver(part%vp,dt_taup0,
     +                              rt_start, rt_zeroes,flag)

               end if

               if (flag == 1) num1000 = num1000 + 1

               if (isnan(rt_zeroes(1)).OR.(rt_zeroes(1)<0)
     +    .OR.isnan(rt_zeroes(2)).OR.(rt_zeroes(2)<0)) then
                numimpos = numimpos + 1  !How many have failed?
                !If they failed (should be very small number), radius,
                !temp remain unchanged
                rt_zeroes(1) = 1.0
                rt_zeroes(2) = part%Tf/part%Tp
               end if


               !Crude way of counting activation. Better would be
               !computing critical radius
               if (part%radius > 1.0e-6 .AND.
     +             part%radius*rt_zeroes(1) < 1.0e-6) then
                   denum = denum + 1

               elseif (part%radius < 1.0e-6 .AND.
     +             part%radius*rt_zeroes(1) > 1.0e-6) then
                   actnum = actnum + 1

               endif

               !Redimensionalize
               part%radius = rt_zeroes(1)*part%radius
               part%Tp = rt_zeroes(2)*part%Tp
        end if

         if (part%radius .gt. 1.0e-2) then
         write(*,'(a15,12e15.6)') 'WARNING: BIG DROPLET', 
     +    part%radius,part%qinf,part%Tp,part%Tf,part%xp(3),
     +    part%Os,part%m_s,part%vp(1),part%vp(2),part%vp(3),
     +    part%res,part%sigm_s
         end if


         !Intermediate Values
         diff(1:3) = part%vp - part%uf
         diffnorm = sqrt(diff(1)**2 + diff(2)**2 + diff(3)**2)
         Rep = 2.0*part%radius*diffnorm/nuf

         myRep_avg = myRep_avg + Rep
         corrfac = (1.0 + 0.15*Rep**(0.687))
         mylwc_sum = mylwc_sum + Volp*rhop*real(part%mult)
         myphiw_sum = myphiw_sum + Volp*rhow
         myphiv_sum = myphiv_sum + Volp

         !Compute Nusselt number for particle:
         !Ranz-Marshall relation
         Nup = 2.0 + 0.6*Rep**(1.0/2.0)*Pra**(1.0/3.0)
         Shp = 2.0 + 0.6*Rep**(1.0/2.0)*Sc**(1.0/3.0)

         !Mass Transfer calculations
         einf = mod_Magnus(part%Tf)

         Eff_C = 2.0*Mw*Gam/(Ru*rhow*part%radius*part%Tp)
         Eff_S = Ion*part%Os*part%m_s*Mw/Ms/(Volp*rhop-part%m_s)
         estar = einf*exp(Mw*Lv/Ru*(1.0/part%Tf-1.0/part%Tp)
     +                                           +Eff_C-Eff_S)
         part%qstar = Mw/Ru*estar/part%Tp/rhoa

        if (ievap .EQ. 1) then
            part%radrhs = Shp/9.0/Sc*rhop/rhow*part%radius*
     +                       taup_i*(part%qinf-part%qstar) !assumes qinf=rhov/rhoa rather than rhov/rhom

 
        else

            part%radrhs = 0.0

            !Also update the temperature directly using BE:
            tmp_coeff = -Nup/3.0/Pra*CpaCpp*rhop/rhow*taup_i
            part%Tp = (part%Tp + tmp_coeff*dt*part%Tf)/(1+dt*tmp_coeff)
        end if

        part%Tprhs_s = -Nup/3.0/Pra*CpaCpp*
     +                     rhop/rhow*taup_i*(part%Tp-part%Tf)
        part%Tprhs_L = 3.0*Lv/Cpp/part%radius*part%radrhs

        part%xrhs(1:3) = part%vp(1:3)
        part%vrhs(1:3) = corrfac*taup_i*(part%uf(1:3)-part%vp(1:3))
     +                               + part_grav(1:3)

        part%res = part%res + dt


      part => part%next
      end do


      !Loop through and make sure to remove the rare case where droplet
      !explodes due to BE procedure
      !Still not 100% certain why this happens, but only occurs under
      !very specific and rare circumstances
      part => first_particle
      do while (associated(part))

       if (part%radius .gt. 1.0) then
         write(*,'(a16,12e15.6)') 'DHR BIG DROPLET: ',
     +    part%radius,part%qinf,part%Tp,part%Tf,part%xp(3),
     +    part%Os,part%m_s,part%vp(1),part%vp(2),part%vp(3),
     +    part%res
          call destroy_particle
       else
           part=> part%next
        end if
      end do

      !Enforce nonperiodic bcs (either elastic or destroying particles)
      call particle_bcs_nonperiodic

      !Check to see if particles left processor
      !If they did, remove from one list and add to another

      call particle_exchange

      !Now enforce periodic bcs 
      !just updates x,y locations if over xl,yl or under 0
      call particle_bcs_periodic

      call particle_coupling_update

      call particle_coupling_exchange
      !Finally, now that coupling and statistics arrays are filled, 
      !Transpose them back to align with the velocities:
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,1),
     +                partsrc(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,2),
     +                partsrc(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,3),
     +                partsrc(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(partTsrc_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partTsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(partHsrc_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partHsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(partTEsrc_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partTEsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(mwsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                mwsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)


      !Try only calling these when the history data is being written:
      if(mtrans) then
      call ztox_trans(upwp_t(0:nnz+1,iys:iye,mxs:mxe),
     +                upwp(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,1),
     +                vpsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,2),
     +                vpsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,3),
     +                vpsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,1),
     +                vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,2),
     +                vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,3),
     +                vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(Tpsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                Tpsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(Tpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                Tpsqrsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(Tfsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                Tfsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(qfsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                qfsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(wpTpsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                wpTpsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(partcount_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partcount(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(radsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                radsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(rad2sum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                rad2sum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(multcount_t(0:nnz+1,iys:iye,mxs:mxe),
     +                multcount(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(qstarsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                qstarsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      end if

      !Get particle count:
      numpart = 0
      myradavg = 0.0
      myradmin=1000.0
      myradmax = 0.0
      mytempmin = 1000.0
      mytempmax = 0.0
      myqmin = 1000.0
      myqmax = 0.0
      part => first_particle
      do while (associated(part))
      numpart = numpart + 1
      myradavg = myradavg + part%radius
      if (part%radius .gt. myradmax) myradmax = part%radius
      if (part%radius .lt. myradmin) myradmin = part%radius
      if (part%Tp .gt. mytempmax) mytempmax = part%Tp
      if (part%Tp .lt. mytempmin) mytempmin = part%Tp
      if (part%qstar .gt. myqmax) myqmax = part%qinf
      if (part%qstar .lt. myqmin) myqmin = part%qinf
      part => part%next
      end do


      !Compute total number of particles
      call mpi_allreduce(numpart,tnumpart,1,mpi_integer,mpi_sum,
     +                   mpi_comm_world,ierr)



      call mpi_allreduce(denum,tdenum,1,mpi_integer,mpi_sum,
     +                   mpi_comm_world,ierr)

      call mpi_allreduce(actnum,tactnum,1,mpi_integer,mpi_sum,
     +                   mpi_comm_world,ierr)

      tmpbuf(1) = myRep_avg
      tmpbuf(2) = mylwc_sum
      tmpbuf(3) = myphiw_sum
      tmpbuf(4) = myphiv_sum
      tmpbuf(5) = myradavg


      !calculate average particle residence time
      call mpi_allreduce(avgres,tavgres,1,mpi_real8,mpi_sum,
     +    mpi_comm_world,ierr)

      !Combine all reals that are being summed:
      call mpi_allreduce(tmpbuf,tmpbuf_rec,5,mpi_real8,mpi_sum,
     +    mpi_comm_world,ierr)

      call mpi_allreduce(num_destroy,tnum_destroy,1,mpi_integer,mpi_sum,
     +   mpi_comm_world,ierr)

      call mpi_allreduce(num100,tnum100,1,mpi_integer,mpi_sum,
     +    mpi_comm_world,ierr)
      call mpi_allreduce(num1000,tnum1000,1,mpi_integer,mpi_sum,
     +    mpi_comm_world,ierr)
      call mpi_allreduce(numimpos,tnumimpos,1,mpi_integer,mpi_sum,
     +    mpi_comm_world,ierr)

      Rep_avg = tmpbuf_rec(1)
      lwc = tmpbuf(2)
      phiw = tmpbuf_rec(3)
      phiv = tmpbuf_rec(4)
      radavg = tmpbuf_rec(5)

      phiw = phiw/xl/yl/zl/rhoa
      phiv = phiv/xl/yl/zl
      Rep_avg = Rep_avg/tnumpart
      radavg = radavg/tnumpart
      tavgres = tavgres/tnum_destroy

      !Min and max radius
      call mpi_allreduce(myradmin,radmin,1,mpi_real8,mpi_min,
     +                   mpi_comm_world,ierr)
      call mpi_allreduce(myradmax,radmax,1,mpi_real8,mpi_max,
     +                   mpi_comm_world,ierr)

      call mpi_allreduce(mytempmin,tempmin,1,mpi_real8,mpi_min,
     +                   mpi_comm_world,ierr)
      call mpi_allreduce(mytempmax,tempmax,1,mpi_real8,mpi_max,
     +                   mpi_comm_world,ierr)

      call mpi_allreduce(myqmin,qmin,1,mpi_real8,mpi_min,
     +                   mpi_comm_world,ierr)
      call mpi_allreduce(myqmax,qmax,1,mpi_real8,mpi_min,
     +                   mpi_comm_world,ierr)



      end subroutine particle_update_BE
      subroutine gauss_newton_2d(vnext,h,vec1,vec2,flag)
        use particles
        implicit none

        real, intent(in) :: vnext(3), h, vec1(2)
        real, intent(out) :: vec2(2)
        integer, intent(out) :: flag
        real :: error = 1E-8, fv1(2), fv2(2), v1(2), v_output(3), rel
        real :: diff, temp1(2), temp2(2), relax, coeff, correct(2)
        real, dimension(1:2, 1:2) :: J, fancy, inv, finalJ
        integer :: iterations, neg, counts

        iterations = 0
        flag = 0

        v1 = vec1
        fv2 = (/1., 1./)
        coeff = 0.1
        do while ((sqrt(dot_product(fv2, fv2)) > error).AND.
     +              (iterations<1000))


                !print*, 'conditions: ', part%Tf, part%qinf, 'error', 
                !+          sqrt(dot_product(fv2, fv2)), 'vel', part%uf

                iterations = iterations + 1
                !print*, 'vnext', vnext, 'v1', v1, 'v_output', v_output,
                !+          'h', h

                call ie_vrt_nd(vnext, v1(1), v1(2), v_output, fv1, h)
                !print*,'fv1', fv1
                call jacob_approx_2d(vnext, v1(1), v1(2), h, J)

                !print*,'error ', sqrt(dot_product(fv1, fv1))

                fancy = matmul(transpose(J), J)

                call inverse_finder_2d(fancy, inv)

                finalJ = matmul(inv, transpose(J))

                !print*,'corrective term', matmul(finalJ, fv1)
                correct = matmul(finalJ,fv1)
                vec2 = v1 - correct

                call ie_vrt_nd(vnext, v1(1), v1(2),v_output,temp1,h)
                call ie_vrt_nd(vnext, vec2(1), vec2(2),v_output,temp2,h)

                diff = sqrt(dot_product(temp1,temp1))-
     +                 sqrt(dot_product(temp2,temp2))

                if (sqrt(dot_product(correct,correct))<1E-8) then
                        !flag = 1
                        EXIT
                end if

                relax = 1.0
                counts = 0

               do while ((diff<0).OR.(vec2(1)<0).OR.(vec2(2)<0)
     +                   .OR.isnan(vec2(1)))
                        counts = counts + 1
                        coeff = 0.5
                        relax = relax * coeff
                        vec2 = v1-matmul(finalJ,fv1)*relax
                call ie_vrt_nd(vnext, vec2(1), vec2(2),v_output,temp2,h)
                        diff = sqrt(dot_product(temp1,temp1))-
     +                  sqrt(dot_product(temp2,temp2))

                        if (counts>10) EXIT
                end do

                v1 = vec2

                call ie_vrt_nd(vnext, vec2(1), vec2(2), v_output, fv2,h)
        end do
      if (iterations == 100) flag = 1
      if (isnan(vec2(1)) .OR. vec2(1)<0 .OR. isnan(vec2(2))
     + .OR. vec2(2)<0) flag = 1

      end subroutine gauss_newton_2d
      subroutine LV_solver(vnext,h,vec1,vec2,flag)
        use particles
        implicit none

        real, intent(in) :: vnext(3),h, vec1(2)
        real, intent(out) :: vec2(2)
        integer, intent(out) :: flag
        real :: error = 1E-8, fv1(2), fv2(2), v1(2), v_output(3), rel
        real :: diff, lambda,lup,ldown
        real :: C(2), newC(2), gradC(2), correct(2)
        real, dimension(1:2, 1:2) :: J,I,g,invg
        integer :: iterations, neg

        I = reshape((/1, 0, 0, 1/),shape(I))
        iterations = 0
        flag = 0
        v1 = vec1
        fv2 = (/1., 1./)

        lambda = 0.001
        lup = 2.0
        ldown = 2.0

        do while ((sqrt(dot_product(fv2, fv2)) > error).AND.
     +              (iterations<1000))

        iterations = iterations + 1
        call jacob_approx_2d(vnext, v1(1), v1(2), h,J)

        call ie_vrt_nd(vnext, v1(1), v1(2),v_output,fv1,h)
        g = matmul(transpose(J),J)+lambda*I
        gradC = matmul(transpose(J),fv1)
        C = 0.5*fv1*fv1

        call inverse_finder_2d(g, invg)
        correct = matmul(invg, gradC)
        if (sqrt(dot_product(correct,correct)) < 1E-12) then
                !flag = 1
                EXIT
        end if

        vec2 = v1 - correct
        call ie_vrt_nd(vnext, vec2(1), vec2(2),v_output,fv2,h)
        newC = 0.5*fv2*fv2

        if (sqrt(dot_product(newC,newC))<sqrt(dot_product(C,C))) then
                v1 = vec2
                lambda = lambda/ldown
        else
                lambda = lambda*lup
        end if

        end do

        if (iterations==1000) then
                flag = 1
                !print*,'not converged ', sqrt(dot_product(fv2, fv2))
        end if

        if (vec2(1) < 0 .OR. vec2(2) < 0) then
                flag = 1
                !print*, 'negative', sqrt(dot_product(fv2,fv2))
        end if


      end subroutine LV_solver
      subroutine jacob_approx_2d(vnext, rnext, tnext, h, J)
        use particles
        implicit none
        integer :: n

        real, intent(in) :: vnext(3), rnext, tnext, h
        real, intent(out), dimension(1:2, 1:2) :: J
        real :: diff = 0, v_output(3), rt_output(2),xper(2),
     +          fxper(2), ynext(2),xper2(2),fxper2(2)

        diff = 1E-12

        ynext(1) = rnext
        ynext(2) = tnext

        call ie_vrt_nd(vnext, rnext, tnext, v_output, rt_output, h)

        xper = ynext
        xper2 = ynext

        do n=1, 2
                xper(n) = xper(n) + diff
                xper2(n) = xper2(n) - diff
                call ie_vrt_nd(vnext, xper(1), xper(2),v_output,
     +           fxper,h)
                call ie_vrt_nd(vnext, xper2(1), xper2(2),v_output,
     +           fxper2,h)
                !              if (.NOT. 
                !+
                !sqrt(dot_product((fxper-rt_output),(fxper-rt_output)))
                !> 0)
                J(:, n) = (fxper-rt_output)/diff
                xper(n) = ynext(n)
                xper2(n) = ynext(n)
        end do

      end subroutine jacob_approx_2d
      subroutine inverse_finder_2d(C, invC)
        implicit none
        real :: det
        real, dimension(1:2, 1:2), intent(in) :: C
        real, dimension(1:2, 1:2), intent(out) :: invC

        det = C(1, 1) * C(2, 2) - C(1, 2) * C(2, 1)

        invC = reshape((/C(2, 2), -C(2,1), -C(1, 2), C(1, 1)/),
     +   shape(invC))
        invC = (1./det)*invC

      end subroutine inverse_finder_2d
      subroutine ie_vrt_nd(vnext, tempr, tempt, v_output,rt_output, h)
      use pars
      use particles
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'

      real, intent(in) :: vnext(3), tempr, tempt, h
      real, intent(out) :: v_output(3), rT_output(2)

      real :: esa, dnext,  m_w, rhop, Rep, taup,
     +         vprime(3), rprime, Tprime, qstr, Shp, Nup, dp, VolP,
     +         diff(3), diffnorm, Tnext, rnext, T
      real :: taup0, g(3), mod_Magnus


        taup0 = (((part%m_s)/((2./3.)*pi2*radius_init**3) + rhow) *
     +   (radius_init*2)**2)/(18*rhoa*nuf)
        g(1:3) = part_grav(1:3)

        ! quantities come in already non-dimensionalized, so must be
        ! converted back;
        ! velocity is not non-dimensionalized so no need to change
        rnext = tempr * part%radius
        Tnext = tempt * part%Tp
        dnext = rnext * 2.

        esa = mod_Magnus(part%Tf)
        VolP = (2./3.)*pi2*rnext**3
        rhop = (part%m_s + VolP*rhow) / VolP

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!! Velocity !!!
        diff(1:3) = part%uf - vnext
        diffnorm = sqrt(diff(1)**2 + diff(2)**2 + diff(3)**2)
        Rep = dnext * diffnorm/nuf
        taup = (rhop * dnext**2)/(18.0*rhoa*nuf)
        vprime(1:3) = (1. + 0.15 * (Rep**0.687)) * (1./taup) *
     +  diff(1:3) - g(1:3)
        vprime(1:3) = vprime(1:3) * taup0 ** 2
        !!!!!!!!!!!!!!!!

        !!! Humidity !!!
        qstr = (Mw/(Ru*Tnext*rhoa)) * esa *
     +  exp(((Lv*Mw/Ru)*((1./part%Tf) - (1./Tnext))) +
     +  ((2.*Mw*Gam)/(Ru*rhow*rnext*Tnext)) -
     +  ((Ion*part%Os*part%m_s*(Mw/Ms))/(Volp*rhop-part%m_s)))
        !!!!!!!!!!!!!!!!!!

        !!! Radius !!!
        Shp = 2. + 0.6 * Rep**(1./2.) * Sc**(1./3.)
        rprime = (1./9.) * (Shp/Sc) * (rhop/rhow) *
     +  (rnext/taup) * (part%qinf - qstr)
        rprime = rprime * (taup0/part%radius)
        !!!!!!!!!!!!!!!!!

        !!! Temperature !!!
        Nup = 2. + 0.6*Rep**(1./2.)*Pra**(1./3.);

        Tprime = -(1./3.)*(Nup/Pra)*CpaCpp*(rhop/rhow)*
     +           (1./taup)*(Tnext-part%Tf) +
     +           3.*Lv*(1./(rnext*Cpp))*rprime*(part%radius/taup0)
        Tprime = Tprime * (taup0/part%Tp)
        !!!!!!!!!!!!!!!!!

        ! velocity is not non-dimensionalized so it does not need to be
        ! changed back
        v_output(1:3) = vnext(1:3) - part%vp(1:3) - h * vprime(1:3)
        rT_output(1) = rnext/part%radius - 1.0  - h*rprime
        rT_output(2) = Tnext/part%Tp - 1.0  - h*Tprime

      end subroutine ie_vrt_nd
      subroutine rad_solver2(guess,mflag)
      use pars
      use particles
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'

      real, intent(OUT) :: guess
      integer, intent(OUT) :: mflag
      real :: a, c, esa, Q, R, M, val, theta, S, T
      real :: mod_Magnus

      mflag = 0
      esa = mod_Magnus(part%Tf)

      a = -(2*Mw*Gam)/(Ru*rhow*part%Tf)/
     +    LOG((Ru*part%Tf*rhoa*part%qinf)/(Mw*esa))
      c = (Ion*part%Os*part%m_s*(Mw/Ms))/((2.0/3.0)*pi2*rhow)/
     +    LOG((Ru*part%Tf*rhoa*part%qinf)/(Mw*esa))

      Q = (a**2.0)/9.0
      R = (2.0*a**3.0+27.0*c)/54.0
      M = R**2.0-Q**3.0
      val = (R**2.0)/(Q**3.0)

      if (M<0) then
        theta = acos(R/sqrt(Q**3.0))
        guess = -(2*sqrt(Q)*cos((theta-pi2)/3.0))-a/3.0

        if (guess < 0) then
        guess = -(2*sqrt(Q)*cos((theta+pi2)/3.0))-a/3.0
        end if

      else
        S = -(R/abs(R))*(abs(R)+sqrt(M))**(1.0/3.0)
        T = Q/S
        guess = S + T - a/3.0

        if (guess < 0) then
                guess = part%radius
                mflag = 1
        end if
      end if

      end subroutine rad_solver2

      subroutine create_particle(xp,vp,Tp,m_s,Os,mult,rad_init,idx)
      use particles
      use pars
      implicit none

      real :: xp(3),vp(3),Tp,qinfp,rad_init,pi,m_s,Os
      integer :: idx
      integer*8 :: mult

      if (.NOT. associated(first_particle)) then
         allocate(first_particle)
         part => first_particle
         nullify(part%next,part%prev)
      else
         !Add to beginning of list since it's more convenient
         part => first_particle
         allocate(part%prev)
         first_particle => part%prev
         part%prev%next => part
         part => first_particle
         nullify(part%prev)
      end if

      pi = 4.0*atan(1.0)
  
      part%xp(1:3) = xp(1:3)
      part%vp(1:3) = vp(1:3)
      part%Tp = Tp
      part%radius = rad_init
      part%uf(1:3) = vp(1:3)
      part%qinf = tsfcc(2)
      part%qstar = 0.008
      part%Tf = Tp
      part%xrhs(1:3) = 0.0
      part%vrhs(1:3) = 0.0 
      part%Tprhs_s = 0.0
      part%Tprhs_L = 0.0
      part%radrhs = 0.0
      part%pidx = idx 
      part%procidx = myid
      part%nbr_pidx = -1
      part%nbr_procidx = -1
      part%mult = mult
      part%res = 0.0
      part%m_s = m_s
      part%Os = Os
      part%dist = 0.0
      part%u_sub(1:3) = 0.0
      part%sigm_s = 0.0

      
      end subroutine create_particle
      function dFdr_uniform(r)
      implicit none
      real :: r,dFdr_uniform 

      dFdr_uniform = 0.5e-6

      end function dFdr_uniform
      subroutine particle_bcs_nonperiodic
      use particles
      use pars
      implicit none
      real :: top,bot,ran2
      real :: xp_init(3),m_s,Os,rad_init,prob
      real :: S,M
      integer :: idx
      integer*8 :: mult

      !Assumes domain goes from [0,xl),[0,yl),[0,zl]
      !Also maintain the number of particles on each proc

      part => first_particle
      do while (associated(part))

      !perfectly elastic collisions on top, bottom walls
      !i.e. location is reflected, w-velocity is negated

      top = zl
      !bot = 0.0 + part%radius
      bot = zw1

      if (part%xp(3) .GT. top) then
         part%xp(3) = top - (part%xp(3)-top)
         part%vp(3) = -part%vp(3)
         part => part%next
      elseif (part%xp(3) .LT. bot) then
         !part%xp(3) = bot + (bot-part%xp(3))
         !part%vp(3) = -part%vp(3)
         !part => part%next

         call destroy_particle

         !Immediately create a new one
         prob = pdf_factor/(1 + pdf_factor)

         !Use this to choose which to generate:

         if (ran2(iseed) .gt. prob) then   !It's accumulation mode
            S = 0.5
            M = -1.95
            Os = 0.5
            mult = mult_a

            !With these parameters, get m_s and rad_init from distribution
            call lognormal_dist(rad_init,m_s,Os,M,S)

          else  !It's coarse mode

            S = 0.45
            M = 0.0
            Os = 1.0
            mult = mult_c

            !With these parameters, get m_s and rad_init from distribution
            call lognormal_dist(rad_init,m_s,Os,M,S)

          end if

         xp_init(1) = ran2(iseed)*(xmax-xmin) + xmin
         xp_init(2) = ran2(iseed)*(ymax-ymin) + ymin
         xp_init(3) = ran2(iseed)*(zi-zw1) + zw1

         idx = 2

         call create_particle(xp_init,vp_init,
     +         Tp_init,m_s,Os,mult,rad_init,idx)

      else
         part => part%next
      end if

      end do


      end subroutine particle_bcs_nonperiodic
      subroutine particle_bcs_periodic
      use particles
      use pars
      implicit none 

      !Assumes domain goes from [0,xl),[0,yl),[0,zl] 
      !Also maintain the number of particles on each proc
      
      part => first_particle
      do while (associated(part))

      !x,y periodic
   
      if (part%xp(1) .GT. xl) then
         part%xp(1) = part%xp(1)-xl
      elseif (part%xp(1) .LT. 0) then
         part%xp(1) = xl + part%xp(1)
      end if

      if (part%xp(2) .GT. yl) then
         part%xp(2) = part%xp(2)-yl
      elseif (part%xp(2) .LT. 0) then
         part%xp(2) = yl + part%xp(2)
      end if

      part => part%next

      end do


      end subroutine particle_bcs_periodic
      subroutine particle_coalesce
      use particles
      use pars
      use kd_tree
      use pars
      use con_data
      implicit none 

      type(particle), pointer :: part_tmp
      type(tree_master_record), pointer :: tree

      real, allocatable :: xp_data(:,:),distances(:),rad_data(:)
      real, allocatable :: vel_data(:,:)
      integer, allocatable :: index_data(:,:),indexes(:)
      integer*8, allocatable :: mult_data(:)
      integer, allocatable :: destroy_data(:)
      integer, allocatable :: coal_data(:),ran_nq(:)

      integer :: i,nq,coal_idx,j,ran_idx,tmp_int,gm,gam_til
      integer :: ns,k_idx,j_idx
      integer*8 :: mult_tmp_j,mult_tmp_k,xi_j,xi_k
      real :: qv(3),dist_tmp,xdist,ydist,zdist,ran2
      real :: phi,K,Pjk,veldiff,E,dV,p_alpha,pvol_j,pvol_k,golovin_b
      real :: rad_j_tmp,rad_k_tmp

      !For whatever reason, the kd-search sometimes misses edge cases,
      !and you should do nq+1 if you actually want nq
      nq = 11

      allocate(xp_data(numpart,3),index_data(numpart,2))
      allocate(vel_data(numpart,3))
      allocate(distances(nq),indexes(nq),ran_nq(2:nq-1))
      allocate(rad_data(numpart),mult_data(numpart),coal_data(numpart))
      allocate(destroy_data(numpart))

      !Loop over particles and fill arrays
      i = 1
      part_tmp => first_particle
      do while (associated(part_tmp))

         xp_data(i,1:3) = part_tmp%xp(1:3) 
         index_data(i,1) = part_tmp%pidx
         index_data(i,2) = part_tmp%procidx

         rad_data(i) = part_tmp%radius
         mult_data(i) = part_tmp%mult
         coal_data(i) = 0

         vel_data(i,1:3) = part_tmp%vp(1:3)

      i = i+1   
      part_tmp => part_tmp%next
      end do

      !Build the kd-tree
      tree => create_tree(xp_data) 

      !Do the search for each of the particles
      part_tmp => first_particle
      do i=1,numpart

         qv(1:3) = xp_data(i,1:3)
         call n_nearest_to(tree,qv,nq,indexes,distances)
         !call n_nearest_to_brute_force(tree,qv,nq,indexes,distances)
 
         !Go back through and assign the shortest distance and index
         !NOTE: Must use 2nd one since it finds itself as nearest neighbor
         !Keep this turned on despite "ineighbor" flag -- consider it a bonus
         part_tmp%dist = sqrt(distances(2))
         part_tmp%nbr_pidx = index_data(indexes(2),1)
         part_tmp%nbr_procidx = index_data(indexes(2),2)

         !Okay, now the particle knows who the nearest nq particles are
         !--> pick one at random and apply coalescence rules
         !Loop over all nq until you find one that hasn't coalesced
         !Don't coalesce if all nq nearby have already done so

         !Set up an array 2,3,...,nq-1
         do j=2,nq-1
            ran_nq(j) = j
         end do

         !Get a random permutation of this array
         do j=2,nq-1
            ran_idx = floor(ran2(iseed)*(nq-3)) + 2 !Get a number between 2 and nq-1
            tmp_int = ran_nq(j) 
            ran_nq(j) = ran_nq(ran_idx)
            ran_nq(ran_idx) = tmp_int
         end do

         !Now loop through these coalescence candidates and take first one that hasn't already coalesced
         coal_idx = -1
         if (coal_data(i) .eq. 0) then
         do j=2,nq-1

             if (.not. coal_data(indexes(ran_nq(j)))) then        !Found one that has not already coalesced
                coal_idx = indexes(ran_nq(j))              !The index of the coalescence candidate in the arrays
                goto 101
             end if
         
         end do
         end if !coal_data = 0
   
101   continue


      !Now apply the coalescence rules to the pair (i,coal_idx) assuming a coal_idx was found (.ne. -1)
      if (coal_idx .ge. 1) then

         phi = ran2(iseed) 
         dV = 2*pi2/3.0*(sqrt(distances(nq-1)))**3  !The volume will be the sphere formed by outermost droplet considered
         E = 1.0  !Collision efficiency -- obviously this needs to be updated from 1.0
         
         veldiff = sqrt( (vel_data(i,1)-vel_data(coal_idx,1))**2 + 
     +                   (vel_data(i,2)-vel_data(coal_idx,2))**2 +  
     +                   (vel_data(i,3)-vel_data(coal_idx,3))**2 )

         if (mult_data(i) .ge. mult_data(coal_idx)) then
            xi_j = mult_data(i)
            xi_k = mult_data(coal_idx)
            j_idx = i
            k_idx = coal_idx
         else
            xi_j = mult_data(coal_idx)
            xi_k = mult_data(i)
            j_idx = coal_idx
            k_idx = i
         end if

         !Choose the kernel:
         K = pi2/2.0*E*veldiff*(rad_data(i) + rad_data(coal_idx))**2

         !Golovin (1963) kernel:
         !pvol_j = pi2*2.0/3.0*rad_data(j_idx)**3.0
         !pvol_k = pi2*2.0/3.0*rad_data(k_idx)**3.0
         !golovin_b = 1.5e3
         !K = golovin_b*(pvol_j + pvol_k)

         Pjk = K*dt/dV*xi_j

         !TESTING: cheat here and hard-code a different dt and dV than the flow
         !since there are issues
         !Pjk = K*1.0/1.0e6*xi_j
         !ns = tnumpart

         !ns = nq-1   !This would be the number of particles in the "cell" according to Shima et al. 2009
         p_alpha = Pjk*(real(ns)*(real(ns)-1.0)/2.0)/(real(ns)/2.0)

         !if (p_alpha .gt. 1) write(*,*) 'WARNING: p_alpha > 1'

         if (phi .lt. p_alpha-floor(p_alpha)) then
            gm = floor(p_alpha) + 1
         else
            gm = floor(p_alpha)
         end if


         if (gm .gt. 0) then  !Only update radii and multiplicities if the coin flip indicates
       
            gam_til = min(gm,floor(real(xi_j)/real(xi_k)))


            if (xi_j - gam_til*xi_k .gt. 0) then

            !Update particle j's multiplicity
            mult_data(j_idx) = mult_data(j_idx)-gam_til*mult_data(k_idx)
               
            !Update particle k's radius
            rad_data(k_idx) = (gam_til*rad_data(j_idx)**3 + 
     +                         rad_data(k_idx)**3)**(1.0/3.0)


            elseif (xi_j - gam_til*xi_k .eq. 0) then

            mult_tmp_j = floor(real(mult_data(k_idx))/2.0)
            mult_tmp_k = mult_data(k_idx) -
     +                   floor(real(mult_data(k_idx))/2.0)
       
            mult_data(j_idx) = mult_tmp_j
            mult_data(k_idx) = mult_tmp_k

            rad_j_tmp = (gam_til*rad_data(j_idx)**3 +
     +                         rad_data(k_idx)**3)**(1.0/3.0)
            rad_k_tmp = (gam_til*rad_data(j_idx)**3 +
     +                         rad_data(k_idx)**3)**(1.0/3.0)

            rad_data(j_idx) = rad_j_tmp
            rad_data(k_idx) = rad_k_tmp

            end if


          end if !gm .gt. 0

       !Now exclude both of these from checking again
       coal_data(coal_idx) = 1
       coal_data(i) = 1
      end if  !coal_idx .gt. 1
         
         

      part_tmp => part_tmp%next
      end do

      !Now finally update the particle linked list
      i = 1
      part_tmp => first_particle
      do while (associated(part_tmp))

         !Only things which should change are radius and multiplicity
         part_tmp%radius = rad_data(i)
         part_tmp%mult = mult_data(i)

      i = i+1   
      part_tmp => part_tmp%next
      end do

      !Finally remove dead particles from coalescence
      i = 1
      numpart = 0
      part => first_particle
      do while (associated(part))
         if (mult_data(i) .eq. 0) then
            call destroy_particle
         else
            numpart = numpart + 1
            part => part%next
         end if

      i = i+1
      end do



      call destroy_tree(tree)
      deallocate(xp_data,index_data)
      deallocate(vel_data,distances,indexes,ran_nq)
      deallocate(rad_data,mult_data,coal_data,destroy_data)


      end subroutine particle_coalesce
      subroutine SFS_velocity
c     This subroutine calculate the SFS velocity for spray drops
c     This use Weil et al. (2004) fromulation
c     Calculations assume isotropic turbulece
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
      implicit none
      include 'mpif.h'
      real :: sigm_sdxp,sigm_sdyp,sigm_sdzp,vis_sp
      real :: sigm_su,sigm_sl,us_ran,gasdev,tengz,englez_bar
      real :: engsbz_bar,sigm_w, sigm_ws
      real :: L_flt,epsn,fs,C0,a1,a2,a3,sigm_sprev,fs1
      real :: weit,weit1,weit3,weit4, T_lagr
      real :: us(3)
      integer :: ix,iy,iz,izp1,izm1,ind,iz_part,ierr
c       ---initialize -------
        fs = 0.0
        C0 = 0.0
        T_lagr = 0.0
        l_flt = 0.0
        epsn = 0.0
        tengz = 0.0
        englez_bar =0.0
        engsbz_bar = 0.0
        sigm_s = 0.0
        sigm_sdx = 0.0
        sigm_sdy = 0.0
        sigm_sdz = 0.0
        sigm_su = 0.0
        sigm_sl  = 0.0
        us_ran = 0.0
        us = 0.0
        sigm_ws = 0.0
        sigm_w = 0.0
c       ------------------
c       compute sigma squre (sigm_s) based on subgrid energy field
c       -----------------       
        do iz =izs,ize
        izp1 = iz+1
        izm1 = iz-1
        weit = dzw(iz)/(dzw(iz)+dzw(izp1))
        weit1 = 1-weit
        weit3 = dzw(izm1)/(dzw(iz)+dzw(izm1))
        weit4 = 1-weit3
        do ix =1,nnx
        do iy = iys,iye
           sigm_s(ix,iy,iz) = 2.0*e(ix,iy,iz)/3.0
           vis_ss(ix,iy,iz) = vis_s(ix,iy,3,iz)
           sigm_sdx(ix,iy,iz) =  sigm_s(ix,iy,iz)       !for xderiv
           sigm_sdy(ix,iy,iz) =  sigm_s(ix,iy,iz)       !for yderiv
c       --------------------
c       calculate z derivative of sigma_s
c       this will be at the u-point!
c       --------------------
           sigm_sdz(ix,iy,iz)=
     +     (sigm_s(ix,iy,iz)-sigm_s(ix,iy,izm1))*dzw_i(iz)
        end do
        end do

c       -------------------
c       calculate x derivatives of sigma_s
c       -------------------
         call xderivp(sigm_sdx(1,iys,iz),trigx(1,1),xk(1),nnx,iys,iye)

        end do

c       ------------------
c       calculate y derivative of sigma_s
c       ------------------
        call yd_mpi(sigm_sdy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c       ----------------
c      calculate extented fileds of sigm_s and its derivatives
c       ---------------
        call fill_extSFS
        !Loop over the linked list of particles:
        part => first_particle
        do while (associated(part))

        ! interpolate sigm_s and its derivative at particle location
        sigm_sprev = part%sigm_s
        call sigm_interp(sigm_sdxp,sigm_sdyp,sigm_sdzp,vis_sp,iz_part)

        part%sigm_s = abs(part%sigm_s)  !Interpolation near surface can give small negative numbers

c       -----------------
c       calculate the subgrid velocity Weil et al ,2004-isotropic turb.
c       ----------------
         l_flt = (2.25*dx*dy*dzw(iz_part+1))**(1.0/3.0)! filtered with
         epsn = (0.93/l_flt)*(3.0*part%sigm_s/2.0)**(1.5) ! tur. dis.rt
         ! TKE : resolved + Subgrid at grid center
        !---------------------------------
c        tot_eng = (engsbz(iz_part+1) + engz(iz_part+1))
         englez_bar = 0.5*(englez(iz_part)+englez(iz_part+1))
         engsbz_bar = 0.5*(engsbz(iz_part)+engsbz(iz_part+1))
         tengz = englez_bar + engsbz_bar
        !---------------------------------
        ! Calculate fs basd on w-componet of velocity
         sigm_w = 0.5*(wps(iz_part)+wps(iz_part+1))
         sigm_ws = 0.5*(engsbz(iz_part)+engsbz(iz_part+1))/3.0

c       ---------write for single droplet ---------
c         write(*,*)'sigm_w:', sigm_w,sigm_ws
c       ------------------------------------------
         if(tengz.gt.0.0)then
c        fs = engsbz(iz_part+1)/(engsbz(iz_part+1) + engz(iz_part+1))
c        fs = engsbz_bar/tengz
         fs = sigm_ws/(sigm_w + sigm_ws)
        else
         fs =0.0
        end if
         C0 = 6.0  ! Changed to 6.0 from 3.0 Indrjith 11-20-17

c       ---------Check for single-part----------
         T_lagr = 2*part%sigm_s/(C0*epsn)   ! Lagrangian time scale
c         write(*,*) 'L_time:',part%xp(3),T_lagr
c       -----------------------------------
c      ------------------
c        Calculate subgrid velocity components
c       -----------------
         us(1:3) = part%u_sub(1:3)
c      -----x component ----------------
         a1 = 0.0
         a2 =0.0
         a3 = 0.0
         a1 =(-0.5)*fs*C0*epsn*part%u_sub(1)/part%sigm_s
         a2 =0.5*(1.0/part%sigm_s)*part%u_sub(1)*
     +          (part%sigm_s - sigm_sprev)/dt
         a3 = 0.5*sigm_sdxp
         us_ran = sqrt(fs*C0*epsn*dt)*gasdev(iseed)
         part%u_sub(1) = (a1+a2+a3)*dt + us_ran

        ! Check to see it blows up 
        if(isnan(part%u_sub(1)))then
        write(*,*)'u-nan',a1,a2,a3,us_ran,part%sigm_s
        write(*,*)'part-con',fs,epsn,dzw(iz_part+1),iz_part
        write(*,*)'part-nan',part%xp(1),part%xp(2),part%xp(3)
        write(*,*)'energy',engsbz(iz_part+1),engz(iz_part+1)
        write(*,*)'derev-nan',sigm_sdxp,sigm_sdyp,sigm_sdzp
        stop
        end if
c       ----------------------------------
c       -----y component ----------------
        a1 = 0.0
        a2 =0.0
        a3 = 0.0
        a1 =(-0.5)*fs*C0*epsn*part%u_sub(2)/part%sigm_s
        a2 =0.5*(1.0/part%sigm_s)*part%u_sub(2)*
     +          (part%sigm_s - sigm_sprev)/dt
        a3 = 0.5*sigm_sdyp
        us_ran = sqrt(fs*C0*epsn*dt)*gasdev(iseed)
        part%u_sub(2) = (a1+a2+a3)*dt + us_ran


        if(isnan(part%u_sub(2)))then
        write(*,*)'v-nan',a1,a2,a3,us_ran,part%sigm_s
        write(*,*)'part-nan',part%xp(1),part%xp(2),part%xp(3)
        write(*,*)'derev-nan',sigm_sdxp,sigm_sdyp,sigm_sdzp
        stop
        end if

c       ----------------------------------
c       -----z component ----------------
        a1 = 0.0
        a2 =0.0
        a3 = 0.0
        a1 =(-0.5)*fs*C0*epsn*part%u_sub(3)/part%sigm_s
        a2 =0.5*(1.0/part%sigm_s)*part%u_sub(3)*
     +          (part%sigm_s - sigm_sprev)/dt
        a3 = 0.5*sigm_sdzp
        us_ran = sqrt(fs*C0*epsn*dt)*gasdev(iseed)
        part%u_sub(3) = (a1+a2+a3)*dt + us_ran

        if(isnan(part%u_sub(3)))then
        write(*,*)'w-nan',a1,a2,a3,us_ran,part%sigm_s
        write(*,*)'part-nan',part%xp(1),part%xp(2),part%xp(3)
        write(*,*)'derev-nan',sigm_sdxp,sigm_sdyp,sigm_sdzp
        stop
        end if

c       -------------------
c       Update particle location and velocity
c       ------------------
        do ind = 1,3
        part%xp(ind) = part%xp(ind) + part%u_sub(ind)*dt
        end do
c       ---------------------------------        
        part => part%next
        end do


        call particle_bcs_nonperiodic
        call particle_exchange
        call particle_bcs_periodic

      numpart = 0
      part => first_particle
      do while (associated(part))
         numpart = numpart + 1
      part => part%next
      end do

      !Compute total number of particles
      call mpi_allreduce(numpart,tnumpart,1,mpi_integer,mpi_sum,
     +                   mpi_comm_world,ierr)


      end subroutine SFS_velocity

      subroutine SFS_position
c     Use a more simplistic SFS treatment: Stochastic particle
c     position rather than velocity
c     Designed to be consistent with LES subgrid eddy diffusivity
c     Does not use Weil et al. 2004 formulation at all
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
      implicit none
      include 'mpif.h'

      real :: sigm_sdxp,sigm_sdyp,sigm_sdzp,vis_sp
      real :: phim,phis,psim,psis,zeta
      real :: dadz,gasdev
      integer :: ix,iy,iz,izp1,izm1,ind,iz_part,ierr

c       ---initialize -------
        sigm_s = 0.0

c       ------------------
c       compute sigma squre (sigm_s) based on subgrid energy field
c       -----------------       
        do iz =izs,ize
        izp1 = iz+1
        izm1 = iz-1
        do ix =1,nnx
        do iy = iys,iye
           sigm_s(ix,iy,iz) = 2.0*e(ix,iy,iz)/3.0
           vis_ss(ix,iy,iz) = vis_s(ix,iy,1,iz)
        end do
        end do

        end do

        call fill_extSFS

        !Loop over the linked list of particles:
        part => first_particle
        do while (associated(part))

        ! interpolate sigm_s and its derivative at particle location
        call sigm_interp(sigm_sdxp,sigm_sdyp,sigm_sdzp,vis_sp,iz_part)

        !Need vertical derivative of average vis_s:
        !Crude approximation: 0th order interpolation -- simply take
        !d(alphaC)/dz of the w-points surrounding the particle location
        if (part%xp(3) .lt. zw1) then
           !Beneath 1st zw point, use MO to approximate dadz:
           zeta = part%xp(3)/amonin
           call fzol(zeta,phim,phis,psim,psis)
           dadz = utau*vk/phis
        else
           !Assuming diffusivity of temperature
           dadz = (alphaC(iz_part+1,1)-alphaC(iz_part,1))*dzw_i(iz)
        end if


        !Now simply solve Langevin equation:
        part%xp(1) = part%xp(1) + gasdev(iseed)*sqrt(2.0*abs(vis_sp)*dt)
        part%xp(2) = part%xp(2) + gasdev(iseed)*sqrt(2.0*abs(vis_sp)*dt)
        part%xp(3) = part%xp(3) + gasdev(iseed)*sqrt(2.0*abs(vis_sp)*dt)
     +             + dadz*dt


        part => part%next
        end do

        call particle_bcs_nonperiodic
        call particle_exchange
        call particle_bcs_periodic

      numpart = 0
      part => first_particle
      do while (associated(part))
         numpart = numpart + 1
      part => part%next
      end do

      !Compute total number of particles
      call mpi_allreduce(numpart,tnumpart,1,mpi_integer,mpi_sum,
     +                   mpi_comm_world,ierr)


      end subroutine SFS_position

      subroutine particle_neighbor_search_kd
      use particles
      use pars
      use kd_tree
      implicit none 

      type(particle), pointer :: part_tmp
      type(tree_master_record), pointer :: tree

      real, allocatable :: xp_data(:,:),distances(:)
      integer, allocatable :: index_data(:,:),indexes(:)

      integer :: i,nq,idx_check,idx_diff
      real :: qv(3),dist_check,dist_tmp,xdist,ydist,zdist,dist_diff

      !For whatever reason, the kd-search sometimes misses edge cases if
      !nq = 2. nq = 3 fixes it
      nq = 3

      allocate(xp_data(numpart,3),index_data(numpart,2))
      allocate(distances(nq),indexes(nq))

      !Loop over particles and fill arrays
      i = 1
      part_tmp => first_particle
      do while (associated(part_tmp))

         xp_data(i,1:3) = part_tmp%xp(1:3) 
         index_data(i,1) = part_tmp%pidx
         index_data(i,2) = part_tmp%procidx

      i = i+1   
      part_tmp => part_tmp%next
      end do

      !Build the kd-tree
      tree => create_tree(xp_data) 

      !Do the search for each of the particles
      part_tmp => first_particle
      do i=1,numpart

         qv(1:3) = xp_data(i,1:3)
         call n_nearest_to(tree,qv,nq,indexes,distances)
         !call n_nearest_to_brute_force(tree,qv,nq,indexes,distances)
 
         !Go back through and assign the shortest distance and index
         !NOTE: Must use 2nd one since it finds itself as nearest neighbor
         part_tmp%dist = sqrt(distances(2))
         part_tmp%nbr_pidx = index_data(indexes(2),1)
         part_tmp%nbr_procidx = index_data(indexes(2),2)

         !if (myid==0) write(*,*) 'dist = ',part_tmp%dist
         !if (myid==0) write(*,*) 'index = ',part_tmp%pidx

      part_tmp => part_tmp%next
      end do


      !Do a check on the first particle in the array that a brute force
      !search finds the same particle and distance
      dist_check = 1.0e6
      part_tmp => first_particle%next
      do while (associated(part_tmp))

            xdist = first_particle%xp(1) - part_tmp%xp(1)
            ydist = first_particle%xp(2) - part_tmp%xp(2)
            zdist = first_particle%xp(3) - part_tmp%xp(3)

            dist_tmp = sqrt(xdist**2 + ydist**2 + zdist**2)

         if (dist_tmp .lt. dist_check) then
            dist_check = dist_tmp
            idx_check = part_tmp%pidx
          end if

      part_tmp => part_tmp%next
      end do

      idx_diff = idx_check-first_particle%nbr_pidx
      dist_diff = abs(dist_check - first_particle%dist)
      if (idx_diff .ne. 0) then
         write(*,'(a30,3i)') 'WARNING in neighbor,idx:',
     +   idx_check,first_particle%nbr_pidx,myid
      end if
      if (dist_diff .gt. 1.0e-8) then
         write(*,'(a30,2e15.6,i)') 'WARNING in neighbor,dist:',
     +   dist_check,first_particle%dist,myid
      end if

      deallocate(xp_data,index_data)


      end subroutine particle_neighbor_search_kd
      subroutine particle_neighbor_search_brute
      use particles
      use pars
      implicit none 

      type(particle), pointer :: part_ref,part_query
      real :: dist_tmp,xdist,ydist,zdist,distance
      integer :: nbr_pidx,nbr_procidx,ip

      
      !Loop over the points that need to perform search
      part_query => first_particle
      do while (associated(part_query))


         distance = 1.0e6
         ip = 1
         !Loop over points to be searched
         part_ref => first_particle
         do while (associated(part_ref))

            xdist = part_query%xp(1) - part_ref%xp(1)
            ydist = part_query%xp(2) - part_ref%xp(2)
            zdist = part_query%xp(3) - part_ref%xp(3)

            dist_tmp = sqrt(xdist**2 + ydist**2 + zdist**2)


            !Record distance and ID of nearest particle
            if (dist_tmp .lt. part_query%dist) then
               distance = dist_tmp
               nbr_pidx = part_ref%pidx
               nbr_procidx = part_ref%procidx
            end if

         ip = ip + 1
         part_ref => part_ref%next
         end do

         part_query%dist = distance
         part_query%nbr_pidx = nbr_pidx
         part_query%nbr_procidx = nbr_procidx
         

      part_query => part_query%next
      end do


      end subroutine particle_neighbor_search_brute
      subroutine destroy_particle
      use particles
      implicit none

      type(particle), pointer :: tmp

      !Is it the first and last in the list?
      if (associated(part,first_particle) .AND. 
     +    (.NOT. associated(part%next)) ) then
          nullify(first_particle)
          deallocate(part)
      else
        if (associated(part,first_particle)) then !Is it the first particle?
           first_particle => part%next
           part => first_particle
           deallocate(part%prev)
        elseif (.NOT. associated(part%next)) then !Is it the last particle?
           nullify(part%prev%next)
           deallocate(part)
        else
           tmp => part
           part => part%next
           tmp%prev%next => tmp%next
           tmp%next%prev => tmp%prev
           deallocate(tmp)
        end if
      end if
   
      end subroutine destroy_particle
      subroutine particle_exchange
      use pars
      use particles
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'

      type(particle), pointer :: tmp
      integer :: idx,psum,csum
      integer :: ir,itr,itop,itl,il,ibl,ib,ibr
      integer :: istatus(mpi_status_size),ierr
      integer :: status_array(mpi_status_size,16),req(16)
      type(particle), allocatable :: rbuf_s(:),trbuf_s(:)
      type(particle), allocatable :: tbuf_s(:),tlbuf_s(:)
      type(particle), allocatable :: lbuf_s(:),blbuf_s(:)
      type(particle), allocatable :: bbuf_s(:),brbuf_s(:)
      type(particle), allocatable :: rbuf_r(:),trbuf_r(:)
      type(particle), allocatable :: tbuf_r(:),tlbuf_r(:)
      type(particle), allocatable :: lbuf_r(:),blbuf_r(:)
      type(particle), allocatable :: bbuf_r(:),brbuf_r(:)
      type(particle), allocatable :: totalbuf(:)
      
      !Zero out the counters for how many particles to send each dir.
      pr_s=0;ptr_s=0;pt_s=0;ptl_s=0;pl_s=0;pbl_s=0;pb_s=0;pbr_s=0
      
      !As soon as the location is updated, must check to see if it left the proc:
      !May be a better way of doing this, but it seems most reasonable:
      part => first_particle
      do while (associated(part))     

         !First get numbers being sent to all sides:
         if (part%xp(2) .GT. ymax) then 
            if (part%xp(1) .GT. xmax) then !top right
               ptr_s = ptr_s + 1
            elseif (part%xp(1) .LT. xmin) then !bottom right
               pbr_s = pbr_s + 1
            else  !right
               pr_s = pr_s + 1
            end if
         elseif (part%xp(2) .LT. ymin) then
            if (part%xp(1) .GT. xmax) then !top left
               ptl_s = ptl_s + 1
            else if (part%xp(1) .LT. xmin) then !bottom left
               pbl_s = pbl_s + 1
            else  !left
               pl_s = pl_s + 1
            end if
         elseif ( (part%xp(1) .GT. xmax) .AND.
     +           (part%xp(2) .LT. ymax) .AND.
     +           (part%xp(2) .GT. ymin) ) then !top
            pt_s = pt_s + 1
         elseif ( (part%xp(1) .LT. xmin) .AND.
     +           (part%xp(2) .LT. ymax) .AND.
     +           (part%xp(2) .GT. ymin) ) then !bottom
            pb_s = pb_s + 1
         end if
         
         part => part%next
      end do
      
      !Now allocate the send buffers based on these counts:
      allocate(rbuf_s(pr_s),trbuf_s(ptr_s),tbuf_s(pt_s),tlbuf_s(ptl_s))
      allocate(lbuf_s(pl_s),blbuf_s(pbl_s),bbuf_s(pb_s),brbuf_s(pbr_s))

      !Now loop back through the particles and fill the buffers:
      !NOTE: If it finds one, add it to buffer and REMOVE from list
      ir=1;itr=1;itop=1;itl=1;il=1;ibl=1;ib=1;ibr=1

      part => first_particle
      do while (associated(part))
         
         if (part%xp(2) .GT. ymax) then 
            if (part%xp(1) .GT. xmax) then !top right
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',trproc,'TOP RIGHT'
               trbuf_s(itr) = part
               call destroy_particle
               itr = itr + 1 
            elseif (part%xp(1) .LT. xmin) then !bottom right
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',brproc,'BOTTOM RIGHT'
               brbuf_s(ibr) = part
               call destroy_particle
               ibr = ibr + 1
            else   !right
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',rproc,'RIGHT'
               rbuf_s(ir) = part
               call destroy_particle
               ir = ir + 1
            end if
         elseif (part%xp(2) .LT. ymin) then
            if (part%xp(1) .GT. xmax) then !top left
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',tlproc,'TOP LEFT'
               tlbuf_s(itl) = part
               call destroy_particle
               itl = itl + 1
            else if (part%xp(1) .LT. xmin) then !bottom left
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',blproc,'BOTTOM LEFT'
               blbuf_s(ibl) = part
               call destroy_particle
               ibl = ibl + 1
            else  !left
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',lproc,'LEFT'
               lbuf_s(il) = part
               call destroy_particle
               il = il + 1
            end if
         elseif ( (part%xp(1) .GT. xmax) .AND.
     +           (part%xp(2) .LT. ymax) .AND.
     +           (part%xp(2) .GT. ymin) ) then !top
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',tproc,'TOP'
            tbuf_s(itop) = part
            call destroy_particle
            itop = itop + 1
         elseif ( (part%xp(1) .LT. xmin) .AND.
     +           (part%xp(2) .LT. ymax) .AND.
     +           (part%xp(2) .GT. ymin) ) then !bottom
!               write(*,*) 'Proc',myid,'about to send part to 
!     +proc',bproc,'BOTTOM'
            bbuf_s(ib) = part
            call destroy_particle
            ib = ib + 1 
         else
         part => part%next
         end if 
         
      end do

      !Now everyone exchanges the counts with all neighbors:
      !Left/right:
      call MPI_Sendrecv(pr_s,1,mpi_integer,rproc,3,
     +        pl_r,1,mpi_integer,lproc,3,mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pl_s,1,mpi_integer,lproc,4,
     +        pr_r,1,mpi_integer,rproc,4,mpi_comm_world,istatus,ierr)

      !Top/bottom:
      call MPI_Sendrecv(pt_s,1,mpi_integer,tproc,5,
     +        pb_r,1,mpi_integer,bproc,5,mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pb_s,1,mpi_integer,bproc,6,
     +        pt_r,1,mpi_integer,tproc,6,mpi_comm_world,istatus,ierr)

      !Top right/bottom left:
      call MPI_Sendrecv(ptr_s,1,mpi_integer,trproc,7,
     +        pbl_r,1,mpi_integer,blproc,7,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pbl_s,1,mpi_integer,blproc,8,
     +        ptr_r,1,mpi_integer,trproc,8,
     +        mpi_comm_world,istatus,ierr)

       !Top left/bottom right:
      call MPI_Sendrecv(ptl_s,1,mpi_integer,tlproc,9,
     +        pbr_r,1,mpi_integer,brproc,9,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pbr_s,1,mpi_integer,brproc,10,
     +         ptl_r,1,mpi_integer,tlproc,10,
     +         mpi_comm_world,istatus,ierr)

      !Now everyone has the number of particles arriving from every neighbor
      !If the count is greater than zero, exchange:

      !Allocate room to receive from each side
      allocate(rbuf_r(pr_r),trbuf_r(ptr_r),tbuf_r(pt_r),tlbuf_r(ptl_r))
      allocate(lbuf_r(pl_r),blbuf_r(pbl_r),bbuf_r(pb_r),brbuf_r(pbr_r))
     
      !Send to right:
      if (pr_s .GT. 0) then
      call mpi_isend(rbuf_s,pr_s,particletype,rproc,11,
     +               mpi_comm_world,req(1),ierr)
      else
      req(1) = mpi_request_null
      end if

      !Receive from left:
      if (pl_r .GT. 0) then
      call mpi_irecv(lbuf_r,pl_r,particletype,lproc,11,
     +               mpi_comm_world,req(2),ierr)
      else
      req(2) = mpi_request_null
      end if

      !Send to left:
      if (pl_s .GT. 0) then
      call mpi_isend(lbuf_s,pl_s,particletype,lproc,12,
     +               mpi_comm_world,req(3),ierr)
      else
      req(3) = mpi_request_null
      end if

      !Receive from right:
      if (pr_r .GT. 0) then
      call mpi_irecv(rbuf_r,pr_r,particletype,rproc,12,
     +               mpi_comm_world,req(4),ierr)
      else
      req(4) = mpi_request_null
      end if

      !Send to top:
      if (pt_s .GT. 0) then
      call mpi_isend(tbuf_s,pt_s,particletype,tproc,13,
     +                mpi_comm_world,req(5),ierr)
      else
      req(5) = mpi_request_null
      end if
      
      !Receive from bottom:
      if (pb_r .GT. 0) then
      call mpi_irecv(bbuf_r,pb_r,particletype,bproc,13,
     +                mpi_comm_world,req(6),ierr)
      else
      req(6) = mpi_request_null
      end if

      !Send to bottom:
      if (pb_s .GT. 0) then
      call mpi_isend(bbuf_s,pb_s,particletype,bproc,14,
     +                mpi_comm_world,req(7),ierr)
      else
      req(7) = mpi_request_null
      end if
      
      !Recieve from top:
      if (pt_r .GT. 0) then
      call mpi_irecv(tbuf_r,pt_r,particletype,tproc,14,
     +                mpi_comm_world,req(8),ierr)
      else
      req(8) = mpi_request_null
      end if

      !Send to top right:
      if (ptr_s .GT. 0) then
      call mpi_isend(trbuf_s,ptr_s,particletype,trproc,15,
     +                mpi_comm_world,req(9),ierr)
      else
      req(9) = mpi_request_null
      end if
     
      !Receive from bottom left:
      if (pbl_r .GT. 0) then
      call mpi_irecv(blbuf_r,pbl_r,particletype,blproc,15,
     +                mpi_comm_world,req(10),ierr)
      else 
      req(10) = mpi_request_null
      end if
    
      !Send to bottom left:
      if (pbl_s .GT. 0) then
      call mpi_isend(blbuf_s,pbl_s,particletype,blproc,16,
     +                mpi_comm_world,req(11),ierr)
      else
      req(11) = mpi_request_null
      end if
     
      !Receive from top right:
      if (ptr_r .GT. 0) then
      call mpi_irecv(trbuf_r,ptr_r,particletype,trproc,16,
     +                mpi_comm_world,req(12),ierr)
      else 
      req(12) = mpi_request_null
      end if

      !Send to top left:
      if (ptl_s .GT. 0) then
      call mpi_isend(tlbuf_s,ptl_s,particletype,tlproc,17,
     +                mpi_comm_world,req(13),ierr)
      else 
      req(13) = mpi_request_null
      end if
    
      !Receive from bottom right:
      if (pbr_r .GT. 0) then
      call mpi_irecv(brbuf_r,pbr_r,particletype,brproc,17,
     +                mpi_comm_world,req(14),ierr)
      else 
      req(14) = mpi_request_null
      end if
  
      !Send to bottom right:
      if (pbr_s .GT. 0) then
      call mpi_isend(brbuf_s,pbr_s,particletype,brproc,18,
     +                mpi_comm_world,req(15),ierr)
      else
      req(15) = mpi_request_null
      end if
  
      !Receive from top left:
      if (ptl_r .GT. 0) then
      call mpi_irecv(tlbuf_r,ptl_r,particletype,tlproc,18,
     +                mpi_comm_world,req(16),ierr)
      else
      req(16) = mpi_request_null
      end if

      call mpi_waitall(16,req,status_array,ierr)

      !Now add incoming particles to linked list:
      !NOTE: add them to beginning since it's easiest to access (first_particle)

      !Form one large buffer to loop through and add:
      psum = pr_r+ptr_r+pt_r+ptl_r+pl_r+pbl_r+pb_r+pbr_r
      csum = 0
      allocate(totalbuf(psum))
      if (pr_r .GT. 0) then 
         totalbuf(1:pr_r) = rbuf_r(1:pr_r)
         csum = csum + pr_r 
      end if
      if (ptr_r .GT. 0) then 
         totalbuf(csum+1:csum+ptr_r) = trbuf_r(1:ptr_r)
         csum = csum + ptr_r
      end if
      if (pt_r .GT. 0) then 
         totalbuf(csum+1:csum+pt_r) = tbuf_r(1:pt_r)
         csum = csum + pt_r
      end if
      if (ptl_r .GT. 0) then 
         totalbuf(csum+1:csum+ptl_r) = tlbuf_r(1:ptl_r)
         csum = csum + ptl_r
      end if
      if (pl_r .GT. 0) then 
         totalbuf(csum+1:csum+pl_r) = lbuf_r(1:pl_r)
         csum = csum + pl_r
      end if
      if (pbl_r .GT. 0) then 
         totalbuf(csum+1:csum+pbl_r) = blbuf_r(1:pbl_r)
         csum = csum + pbl_r
      end if
      if (pb_r .GT. 0) then 
         totalbuf(csum+1:csum+pb_r) = bbuf_r(1:pb_r)
         csum = csum + pb_r
      end if
      if (pbr_r .GT. 0) then 
         totalbuf(csum+1:csum+pbr_r) = brbuf_r(1:pbr_r)
         csum = csum + pbr_r
      end if

      do idx = 1,psum
        if (.NOT. associated(first_particle)) then
           allocate(first_particle)
           first_particle = totalbuf(idx)
           nullify(first_particle%next,first_particle%prev)
        else
           allocate(first_particle%prev)
           tmp => first_particle%prev
           tmp = totalbuf(idx)
           tmp%next => first_particle
           nullify(tmp%prev)
           first_particle => tmp
           nullify(tmp)
        end if
      end do  
      
      deallocate(rbuf_s,trbuf_s,tbuf_s,tlbuf_s)
      deallocate(lbuf_s,blbuf_s,bbuf_s,brbuf_s)
      deallocate(rbuf_r,trbuf_r,tbuf_r,tlbuf_r)
      deallocate(lbuf_r,blbuf_r,bbuf_r,brbuf_r)
      deallocate(totalbuf)

      end subroutine particle_exchange

      subroutine fill_ext 
      use pars
      use particles
      use fields
      use con_stats
      use con_data
      implicit none
      include 'mpif.h'

      integer :: istatus(mpi_status_size),ierr
      integer :: ix,iy,iz
      !preceding letter: r=right,l=left,t=top,b=bot.
      !_s: buf of things to send TO r,l,t,b
      !_r: buf of things to recv FROM r,l,t,b 
      real :: tbuf_s(nnz+2,iye-iys+1,2,5),tbuf_r(nnz+2,iye-iys+1,3,5)
      real :: bbuf_s(nnz+2,iye-iys+1,3,5),bbuf_r(nnz+2,iye-iys+1,2,5)
      real :: rbuf_s(nnz+2,2,mxe-mxs+1,5),rbuf_r(nnz+2,3,mxe-mxs+1,5)
      real :: lbuf_s(nnz+2,3,mxe-mxs+1,5),lbuf_r(nnz+2,2,mxe-mxs+1,5)

      !Corners:
      real :: trbuf_s(nnz+2,2,2,5),trbuf_r(nnz+2,3,3,5) 
      real :: brbuf_s(nnz+2,2,3,5),brbuf_r(nnz+2,3,2,5)
      real :: blbuf_s(nnz+2,3,3,5),blbuf_r(nnz+2,2,2,5)
      real :: tlbuf_s(nnz+2,3,2,5),tlbuf_r(nnz+2,2,3,5)

      !MPI send counts:
      integer :: rc_s,rc_r,trc_s,trc_r,tc_s,tc_r,tlc_s,tlc_r
      integer :: lc_s,lc_r,blc_s,blc_r,bc_s,bc_r,brc_s,brc_r

      !Debugging:
      real :: xv,yv,zv
 
      !To update the particle ODE in time, need the interpolated
      !velocity field
      !This requires filling uext,vext,wext from nearby procs
      uext = 0.0
      vext = 0.0
      wext = 0.0
      Text = 0.0
      T2ext = 0.0   ! Brian 7/14/14

      !FOR CHECKING INTERPOLATION: create an artificial velocity field
      !to interpolate from:
!      do iz=izs-1,ize+1
!        do iy=iys,iye
!          do ix=1,nnx
!            xv = dx*(ix-1)
!            yv = dy*(iy-1)
!            zv = z(iz)
!            u(ix,iy,iz) = 10.0 
!            v(ix,iy,iz) = 0.0 
!            w(ix,iy,iz) = 0.0 
!            u(ix,iy,iz) = 0.3*xv**6 + 0.5*xv*yv**3*zv + 3.0*xv**3*zv**2+
!     +                      0.7*yv**2*zv**3 + 0.5*zv**3 + 0.1*yv**2
!          end do
!        end do
!      end do

      !First fill the center, since this is just u,v,w on that proc:

      !In the column setup, need to tranpose u,v,w first into u_t,v_t,w_t:
      call xtoz_trans(u(1:nnx,iys:iye,izs-1:ize+1),u_t,nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(v(1:nnx,iys:iye,izs-1:ize+1),v_t,nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(w(1:nnx,iys:iye,izs-1:ize+1),w_t,nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(t(1:nnx,iys:iye,1,izs-1:ize+1),T_t,nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(t(1:nnx,iys:iye,2,izs-1:ize+1),T2_t,nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs) 

      uext(0:nnz+1,iys:iye,mxs:mxe) = u_t(0:nnz+1,iys:iye,mxs:mxe)
      vext(0:nnz+1,iys:iye,mxs:mxe) = v_t(0:nnz+1,iys:iye,mxs:mxe)
      wext(0:nnz+1,iys:iye,mxs:mxe) = w_t(0:nnz+1,iys:iye,mxs:mxe)
      Text(0:nnz+1,iys:iye,mxs:mxe) = T_t(0:nnz+1,iys:iye,mxs:mxe)
      T2ext(0:nnz+1,iys:iye,mxs:mxe) = T2_t(0:nnz+1,iys:iye,mxs:mxe)


      !Recall that SR assign_nbrs assigned rproc,lproc, etc.

      !Going to call 6 sendrecv calls - one for each proc. nbr.:
      
      !Fill the send buffers:
      
      !I know these are redundant, but so I can keep them straight...
      tc_s = 5*(nnz+2)*2*(iye-iys+1)
      tc_r = 5*(nnz+2)*3*(iye-iys+1)
      trc_s = 5*(nnz+2)*2*2
      trc_r = 5*(nnz+2)*3*3
      rc_s = 5*(nnz+2)*(mxe-mxs+1)*2 
      rc_r = 5*(nnz+2)*(mxe-mxs+1)*3 
      tlc_s = 5*(nnz+2)*3*2 
      tlc_r = 5*(nnz+2)*2*3 
      bc_s = 5*(nnz+2)*3*(iye-iys+1) 
      bc_r = 5*(nnz+2)*2*(iye-iys+1)
      blc_s = 5*(nnz+2)*3*3
      blc_r = 5*(nnz+2)*2*2
      lc_s = 5*(nnz+2)*(mxe-mxs+1)*3 
      lc_r = 5*(nnz+2)*(mxe-mxs+1)*2
      brc_s = 5*(nnz+2)*2*3 
      brc_r = 5*(nnz+2)*3*2 
     
      !First u:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,1) = u_t(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,1) = u_t(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,1) = u_t(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,1) = u_t(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,1) = u_t(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,1) = u_t(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,1) = u_t(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,1) = u_t(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !v:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,2) = v_t(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,2) = v_t(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,2) = v_t(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,2) = v_t(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,2) = v_t(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,2) = v_t(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,2) = v_t(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,2) = v_t(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !w:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,3) = w_t(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,3) = w_t(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,3) = w_t(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,3) = w_t(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,3) = w_t(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,3) = w_t(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,3) = w_t(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,3) = w_t(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !T:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,4) = T_t(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,4) = T_t(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,4) = T_t(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,4) = T_t(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,4) = T_t(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,4) = T_t(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,4) = T_t(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,4) = T_t(0:nnz+1,iys:iys+2,mxe-1:mxe)
     
      !T2:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,5) =
     +      T2_t(0:nnz+1,iys:iye,mxe-1:mxe) !Brian 7/14/14
      trbuf_s(1:nnz+2,1:2,1:2,5) = 
     +      T2_t(0:nnz+1,iye-1:iye,mxe-1:mxe) !Brian 7/14/14
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,5) = 
     +      T2_t(0:nnz+1,iye-1:iye,mxs:mxe) !Brian 7/14/14
      brbuf_s(1:nnz+2,1:2,1:3,5) = 
     +      T2_t(0:nnz+1,iye-1:iye,mxs:mxs+2) !Brian 7/14/14
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,5) = 
     +      T2_t(0:nnz+1,iys:iye,mxs:mxs+2) !Brian 7/14/14
      blbuf_s(1:nnz+2,1:3,1:3,5) = 
     +      T2_t(0:nnz+1,iys:iys+2,mxs:mxs+2) !Brian 7/14/14
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,5) = 
     +      T2_t(0:nnz+1,iys:iys+2,mxs:mxe) !Brian 7/14/14
      tlbuf_s(1:nnz+2,1:3,1:2,5) = 
     +      T2_t(0:nnz+1,iys:iys+2,mxe-1:mxe) !Brian 7/14/14
     
      !Zero out recieve buffers 
      rbuf_r=0.0;trbuf_r=0.0;tbuf_r=0.0;tlbuf_r=0.0;lbuf_r=0.0
      blbuf_r=0.0;bbuf_r=0.0;brbuf_r=0.0

      !Left/right:
      call MPI_Sendrecv(rbuf_s,rc_s,mpi_real8,rproc,3,
     +        lbuf_r,lc_r,mpi_real8,lproc,3,mpi_comm_world,istatus,ierr)

      call mpi_barrier(mpi_comm_world,ierr)
      call MPI_Sendrecv(lbuf_s,lc_s,mpi_real8,lproc,4,
     +        rbuf_r,rc_r,mpi_real8,rproc,4,mpi_comm_world,istatus,ierr)

      !Top/bottom:
      call MPI_Sendrecv(tbuf_s,tc_s,mpi_real8,tproc,5,
     +        bbuf_r,bc_r,mpi_real8,bproc,5,mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(bbuf_s,bc_s,mpi_real8,bproc,6,
     +        tbuf_r,tc_r,mpi_real8,tproc,6,mpi_comm_world,istatus,ierr)

      !Top right/bottom left:
      call MPI_Sendrecv(trbuf_s,trc_s,mpi_real8,trproc,7,
     +        blbuf_r,blc_r,mpi_real8,blproc,7,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(blbuf_s,blc_s,mpi_real8,blproc,8,
     +        trbuf_r,trc_r,mpi_real8,trproc,8,
     +        mpi_comm_world,istatus,ierr)

       !Top left/bottom right:
      call MPI_Sendrecv(tlbuf_s,tlc_s,mpi_real8,tlproc,9,
     +        brbuf_r,brc_r,mpi_real8,brproc,9,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(brbuf_s,brc_s,mpi_real8,brproc,10,
     +         tlbuf_r,tlc_r,mpi_real8,tlproc,10,
     +         mpi_comm_world,istatus,ierr)

      !Now fill the ext arrays with the recieved buffers:
      uext(0:nnz+1,iys:iye,mxe+1:mxe+3) = 
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,1)
      uext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,1)
      uext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,1)
      uext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,1)
      uext(0:nnz+1,iys:iye,mxs-2:mxs-1) = 
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,1)
      uext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,1)
      uext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,1)
      uext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,1)
 
      vext(0:nnz+1,iys:iye,mxe+1:mxe+3) = 
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,2)
      vext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,2)
      vext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,2)
      vext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,2)
      vext(0:nnz+1,iys:iye,mxs-2:mxs-1) = 
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,2)
      vext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,2)
      vext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,2)
      vext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,2)

      wext(0:nnz+1,iys:iye,mxe+1:mxe+3) = 
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,3)
      wext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,3)
      wext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,3)
      wext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,3)
      wext(0:nnz+1,iys:iye,mxs-2:mxs-1) = 
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,3)
      wext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,3)
      wext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,3)
      wext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,3)

      Text(0:nnz+1,iys:iye,mxe+1:mxe+3) = 
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,4)
      Text(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,4)
      Text(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,4)
      Text(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,4)
      Text(0:nnz+1,iys:iye,mxs-2:mxs-1) = 
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,4)
      Text(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,4)
      Text(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,4)
      Text(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,4)


      T2ext(0:nnz+1,iys:iye,mxe+1:mxe+3) = 
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,5) !Brian 7/14/14
      T2ext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = 
     +     trbuf_r(1:nnz+2,1:3,1:3,5)!Brian 7/14/14
      T2ext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,5)!Brian 7/14/14
      T2ext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = 
     +     brbuf_r(1:nnz+2,1:3,1:2,5)!Brian 7/14/14
      T2ext(0:nnz+1,iys:iye,mxs-2:mxs-1) = 
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,5)!Brian 7/14/14
      T2ext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = 
     +     blbuf_r(1:nnz+2,1:2,1:2,5)!Brian 7/14/14
      T2ext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,5)!Brian 7/14/14
      T2ext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = 
     +     tlbuf_r(1:nnz+2,1:2,1:3,5)!Brian 7/14/14

      end subroutine fill_ext
      subroutine fill_extSFS
c       This subroutine calculte the extented fields for SFS velocity
c       calculations

      use pars
      use particles
      use fields
      use con_stats
      use con_data
      implicit none
      include 'mpif.h'

      integer :: istatus(mpi_status_size),ierr
      integer :: ix,iy,iz
      real :: sigm_st(0:nnz+1,iys:iye,mxs:mxe)
      real :: sigm_sdxt(0:nnz+1,iys:iye,mxs:mxe)
      real :: sigm_sdyt(0:nnz+1,iys:iye,mxs:mxe)
      real :: sigm_sdzt(0:nnz+1,iys:iye,mxs:mxe)
      real :: vis_st(0:nnz+1,iys:iye,mxs:mxe)
      !preceding letter: r=right,l=left,t=top,b=bot.
      !_s: buf of things to send TO r,l,t,b
      !_r: buf of things to recv FROM r,l,t,b
      real :: tbuf_s(nnz+2,iye-iys+1,2,5),tbuf_r(nnz+2,iye-iys+1,3,5)
      real :: bbuf_s(nnz+2,iye-iys+1,3,5),bbuf_r(nnz+2,iye-iys+1,2,5)
      real :: rbuf_s(nnz+2,2,mxe-mxs+1,5),rbuf_r(nnz+2,3,mxe-mxs+1,5)
      real :: lbuf_s(nnz+2,3,mxe-mxs+1,5),lbuf_r(nnz+2,2,mxe-mxs+1,5)

      !Corners:
      real :: trbuf_s(nnz+2,2,2,5),trbuf_r(nnz+2,3,3,5)
      real :: brbuf_s(nnz+2,2,3,5),brbuf_r(nnz+2,3,2,5)
      real :: blbuf_s(nnz+2,3,3,5),blbuf_r(nnz+2,2,2,5)
      real :: tlbuf_s(nnz+2,3,2,5),tlbuf_r(nnz+2,2,3,5)


      !MPI send counts:
      integer :: rc_s,rc_r,trc_s,trc_r,tc_s,tc_r,tlc_s,tlc_r
      integer :: lc_s,lc_r,blc_s,blc_r,bc_s,bc_r,brc_s,brc_r
      sigm_sext = 0.0
      sigm_sdxext = 0.0
      sigm_sdyext = 0.0
      sigm_sdzext = 0.0
      vis_sext = 0.0
      vis_st = 0.0
      sigm_st = 0.0
      sigm_sdxt = 0.0
      sigm_sdyt = 0.0
      sigm_sdzt = 0.0

              !In the column setup, need to tranpose u,v,w first into
              !u_t,v_t,w_t:
      call xtoz_trans(sigm_s(1:nnx,iys:iye,izs-1:ize+1),sigm_st,nnx,
     +                nnz,mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(sigm_sdx(1:nnx,iys:iye,izs-1:ize+1),sigm_sdxt,nnx,
     +                nnz,mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(sigm_sdy(1:nnx,iys:iye,izs-1:ize+1),sigm_sdyt,nnx,
     +                nnz,mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(sigm_sdz(1:nnx,iys:iye,izs-1:ize+1),sigm_sdzt,nnx,
     +                nnz,mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)

      call xtoz_trans(vis_ss(1:nnx,iys:iye,izs-1:ize+1),vis_st,nnx,
     +                nnz,mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)

      sigm_sext(0:nnz+1,iys:iye,mxs:mxe)
     +                              = sigm_st(0:nnz+1,iys:iye,mxs:mxe)
      sigm_sdxext(0:nnz+1,iys:iye,mxs:mxe)
     +                              = sigm_sdxt(0:nnz+1,iys:iye,mxs:mxe)
      sigm_sdyext(0:nnz+1,iys:iye,mxs:mxe)
     +                              = sigm_sdyt(0:nnz+1,iys:iye,mxs:mxe)
      sigm_sdzext(0:nnz+1,iys:iye,mxs:mxe)
     +                              = sigm_sdzt(0:nnz+1,iys:iye,mxs:mxe)
      vis_sext(0:nnz+1,iys:iye,mxs:mxe)
     +                              = vis_st(0:nnz+1,iys:iye,mxs:mxe)



      !Fill the send buffers:
      tc_s = 5*(nnz+2)*2*(iye-iys+1) 
      tc_r = 5*(nnz+2)*3*(iye-iys+1)
      trc_s = 5*(nnz+2)*2*2 
      trc_r = 5*(nnz+2)*3*3 
      rc_s = 5*(nnz+2)*(mxe-mxs+1)*2 
      rc_r = 5*(nnz+2)*(mxe-mxs+1)*3 
      tlc_s = 5*(nnz+2)*3*2 
      tlc_r = 5*(nnz+2)*2*3 
      bc_s = 5*(nnz+2)*3*(iye-iys+1) 
      bc_r = 5*(nnz+2)*2*(iye-iys+1)
      blc_s = 5*(nnz+2)*3*3 
      blc_r = 5*(nnz+2)*2*2 
      lc_s = 5*(nnz+2)*(mxe-mxs+1)*3 
      lc_r = 5*(nnz+2)*(mxe-mxs+1)*2 
      brc_s = 5*(nnz+2)*2*3
      brc_r = 5*(nnz+2)*3*2

      !First sigm_s:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,1)
     +                          = sigm_st(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,1)
     +                          = sigm_st(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,1)
     +                          = sigm_st(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,1)
     +                          = sigm_st(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,1)
     +                          = sigm_st(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,1)
     +                          = sigm_st(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,1)
     +                          = sigm_st(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,1)
     +                          = sigm_st(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !sigm_sdx:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,2)
     +                          = sigm_sdxt(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,2)
     +                          = sigm_sdxt(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,2)
     +                          = sigm_sdxt(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,2)
     +                          = sigm_sdxt(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,2)
     +                          = sigm_sdxt(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,2)
     +                          = sigm_sdxt(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,2)
     +                          = sigm_sdxt(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,2)
     +                          = sigm_sdxt(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !sigm_sdy:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,3)
     +                         = sigm_sdyt(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,3)
     +                         = sigm_sdyt(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,3)
     +                         = sigm_sdyt(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,3)
     +                         = sigm_sdyt(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,3)
     +                         = sigm_sdyt(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,3)
     +                         = sigm_sdyt(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,3)
     +                         = sigm_sdyt(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,3)
     +                         = sigm_sdyt(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !sigm_sdz
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,4)
     +                          = sigm_sdzt(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,4)
     +                          = sigm_sdzt(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,4)
     +                          = sigm_sdzt(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,4)
     +                          = sigm_sdzt(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,4)
     +                          = sigm_sdzt(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,4)
     +                          = sigm_sdzt(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,4)
     +                          = sigm_sdzt(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,4)
     +                          = sigm_sdzt(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !vis_s
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,5)
     +                          = vis_st(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,5)
     +                          = vis_st(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,5)
     +                          = vis_st(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,5)
     +                          = vis_st(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,5)
     +                          = vis_st(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,5)
     +                          = vis_st(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,5)
     +                          = vis_st(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,5)
     +                          = vis_st(0:nnz+1,iys:iys+2,mxe-1:mxe)


      !Zero out recieve buffers
      rbuf_r=0.0;trbuf_r=0.0;tbuf_r=0.0;tlbuf_r=0.0;lbuf_r=0.0
      blbuf_r=0.0;bbuf_r=0.0;brbuf_r=0.0

      !Left/right:
      call MPI_Sendrecv(rbuf_s,rc_s,mpi_real8,rproc,3,
     +        lbuf_r,lc_r,mpi_real8,lproc,3,mpi_comm_world,istatus,ierr)

      call mpi_barrier(mpi_comm_world,ierr)
      call MPI_Sendrecv(lbuf_s,lc_s,mpi_real8,lproc,4,
     +        rbuf_r,rc_r,mpi_real8,rproc,4,mpi_comm_world,istatus,ierr)

      !Top/bottom:
      call MPI_Sendrecv(tbuf_s,tc_s,mpi_real8,tproc,5,
     +        bbuf_r,bc_r,mpi_real8,bproc,5,mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(bbuf_s,bc_s,mpi_real8,bproc,6,
     +        tbuf_r,tc_r,mpi_real8,tproc,6,mpi_comm_world,istatus,ierr)

      !Top right/bottom left:
      call MPI_Sendrecv(trbuf_s,trc_s,mpi_real8,trproc,7,
     +        blbuf_r,blc_r,mpi_real8,blproc,7,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(blbuf_s,blc_s,mpi_real8,blproc,8,
     +        trbuf_r,trc_r,mpi_real8,trproc,8,
     +        mpi_comm_world,istatus,ierr)

       !Top left/bottom right:
      call MPI_Sendrecv(tlbuf_s,tlc_s,mpi_real8,tlproc,9,
     +        brbuf_r,brc_r,mpi_real8,brproc,9,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(brbuf_s,brc_s,mpi_real8,brproc,10,
     +         tlbuf_r,tlc_r,mpi_real8,tlproc,10,
     +         mpi_comm_world,istatus,ierr)

c     Now fill the ext arrays with the recieved buffers:
      sigm_sext(0:nnz+1,iys:iye,mxe+1:mxe+3) =
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,1)
      sigm_sext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) =
     +     trbuf_r(1:nnz+2,1:3,1:3,1)
      sigm_sext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,1)
      sigm_sext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) =
     +     brbuf_r(1:nnz+2,1:3,1:2,1)
      sigm_sext(0:nnz+1,iys:iye,mxs-2:mxs-1) =
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,1)
      sigm_sext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) =
     +     blbuf_r(1:nnz+2,1:2,1:2,1)
      sigm_sext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,1)
      sigm_sext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) =
     +     tlbuf_r(1:nnz+2,1:2,1:3,1)

      sigm_sdxext(0:nnz+1,iys:iye,mxe+1:mxe+3) =
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,2)
      sigm_sdxext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) =
     +     trbuf_r(1:nnz+2,1:3,1:3,2)
      sigm_sdxext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,2)
      sigm_sdxext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) =
     +     brbuf_r(1:nnz+2,1:3,1:2,2)
      sigm_sdxext(0:nnz+1,iys:iye,mxs-2:mxs-1) =
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,2)
      sigm_sdxext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) =
     +     blbuf_r(1:nnz+2,1:2,1:2,2)
      sigm_sdxext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,2)
      sigm_sdxext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) =
     +     tlbuf_r(1:nnz+2,1:2,1:3,2)

      sigm_sdyext(0:nnz+1,iys:iye,mxe+1:mxe+3) =
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,3)
      sigm_sdyext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) =
     +     trbuf_r(1:nnz+2,1:3,1:3,3)
      sigm_sdyext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,3)
      sigm_sdyext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) =
     +     brbuf_r(1:nnz+2,1:3,1:2,3)
      sigm_sdyext(0:nnz+1,iys:iye,mxs-2:mxs-1) =
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,3)
      sigm_sdyext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) =
     +     blbuf_r(1:nnz+2,1:2,1:2,3)
      sigm_sdyext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,3)
      sigm_sdyext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) =
     +     tlbuf_r(1:nnz+2,1:2,1:3,3)

      sigm_sdzext(0:nnz+1,iys:iye,mxe+1:mxe+3) =
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,4)
      sigm_sdzext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) =
     +     trbuf_r(1:nnz+2,1:3,1:3,4)
      sigm_sdzext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,4)
      sigm_sdzext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) =
     +     brbuf_r(1:nnz+2,1:3,1:2,4)
      sigm_sdzext(0:nnz+1,iys:iye,mxs-2:mxs-1) =
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,4)
      sigm_sdzext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) =
     +     blbuf_r(1:nnz+2,1:2,1:2,4)
      sigm_sdzext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,4)
      sigm_sdzext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) =
     +     tlbuf_r(1:nnz+2,1:2,1:3,4)


      vis_sext(0:nnz+1,iys:iye,mxe+1:mxe+3) =
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,5)
      vis_sext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) =
     +     trbuf_r(1:nnz+2,1:3,1:3,5)
      vis_sext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,5)
      vis_sext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) =
     +     brbuf_r(1:nnz+2,1:3,1:2,5)
      vis_sext(0:nnz+1,iys:iye,mxs-2:mxs-1) =
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,5)
      vis_sext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) =
     +     blbuf_r(1:nnz+2,1:2,1:2,5)
      vis_sext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,5)
      vis_sext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) =
     +     tlbuf_r(1:nnz+2,1:2,1:3,5)

      end subroutine fill_extSFS

      subroutine uf_interp
      use pars
      use particles
      use fields
      use con_stats
      use con_data
      implicit none
       
      integer :: ix,iy,izuv,izw,iz,i,k,j
      integer :: first,last
      real :: xkval,xjval,pj,dxvec(2)
      integer :: ijpts(2,6),kuvpts(6),kwpts(6)
      real :: wt(4,6)
      real :: ran2
      
      !Following Orlando's code, get the "leftmost" node
      !This is just the minimum (i,j,k) on the volume 

      ijpts(1,3) = floor(part%xp(1)/dx) + 1 
      ijpts(2,3) = floor(part%xp(2)/dy) + 1
     
      !Fill in the neighbors:
      ijpts(1,2) = ijpts(1,3)-1
      ijpts(1,1) = ijpts(1,2)-1
      ijpts(1,4) = ijpts(1,3)+1
      ijpts(1,5) = ijpts(1,4)+1
      ijpts(1,6) = ijpts(1,5)+1

      ijpts(2,2) = ijpts(2,3)-1
      ijpts(2,1) = ijpts(2,2)-1
      ijpts(2,4) = ijpts(2,3)+1
      ijpts(2,5) = ijpts(2,4)+1
      ijpts(2,6) = ijpts(2,5)+1
     
      !Finding the k-lhnode is different since grid may be stretched
      !AND since (u,v) and w stored differently
      !Will get a k-index for (u,v) and one for w
      
      !Do (u,v) loop first:
      do iz=0,nnz+1
         if (zz(iz) .GT. part%xp(3)) go to 1199
      end do
1199  continue

      kuvpts(3) = iz-1
      !Then fill in the rest:
      kuvpts(4) = kuvpts(3)+1
      kuvpts(5) = kuvpts(4)+1
      kuvpts(6) = kuvpts(5)+1
      kuvpts(2) = kuvpts(3)-1
      kuvpts(1) = kuvpts(2)-1

      !Do again for w:
      do iz = 0,nnz+1 
         if (z(iz) .GT. part%xp(3)) go to 1120
      end do
1120  continue
   
      kwpts(3) = iz-1
      !Then fill in the rest:
      kwpts(4) = kwpts(3)+1
      kwpts(5) = kwpts(4)+1
      kwpts(6) = kwpts(5)+1
      kwpts(2) = kwpts(3)-1
      kwpts(1) = kwpts(2)-1


!---------
!      !As an aside, use ijpts,kwpts to update the particle numbers at each cell:
!
!      !call inputs  are the nodes associated with the volume that each part. lies in
!      !The +1 on kwpts is since stats are stored at uv (center) points in z
!      call particle_stats(ijpts(1,3),ijpts(2,3),kwpts(3)+1)
!   
!---------

      !Fill in the weights:
      !First for x and y since they are periodic:
      wt(1:4,1:6) = 0.0
      dxvec(1) = dx
      dxvec(2) = dy
      do iz = 1,2
      do j = 1,6
         xjval = dxvec(iz)*(ijpts(iz,j)-1)
         pj = 1.0
         do k = 1,6
            xkval = dxvec(iz)*(ijpts(iz,k)-1)
            if (j .NE. k) then
                  pj = pj*(part%xp(iz)-xkval)/(xjval-xkval)
            end if
         end do
         wt(iz,j) = pj
       end do
       end do
      
     
       !Now compute weights in z-dir
       !There are 2 sections: weights at (u,v) nodes (kuvpts) 
       !And weights computed at w nodes (kwpts)

       !Compute weights at kuvpts
       !Must check to see how close we are to a top/bot boundary
       if (kuvpts(3) == 1) then
          first = 3
          last = 4
          !Set these equal to 1 so uext(-1) won't be accessed
          !Note: the value doesn't matter since weight will be 0
          kuvpts(1) = 1
          kuvpts(2) = 1
       elseif (kuvpts(3) == 0) then
          first = 4
          last = 5
          kuvpts(1) = 1
          kuvpts(2) = 1
          kuvpts(3) = 1
       elseif (kuvpts(3) .LT. 0) then 
          first = 0
          last = 0
       elseif (kuvpts(3) == 2) then 
          first = 2
          last = 5
       !Between top cell center and the domain boundary
       elseif (kuvpts(3) == nnz) then
          first = 2
          last = 3
          kuvpts(4) = nnz
          kuvpts(5) = nnz
          kuvpts(6) = nnz
       elseif (kuvpts(3) .GT. nnz) then
          first = 0
          last = 0
       !Between 2nd to last and last cell center at top
       elseif (kuvpts(3) == nnz-1) then
          first = 3
          last = 4
          kuvpts(5) = nnz
          kuvpts(6) = nnz
       elseif (kuvpts(3) == nnz-2) then
          first = 2
          last = 5
       else
          first = 1
          last = 6
       end if

       !Recall that wt has been set to zero, so
       !weights will be zero if (first,last) isn't (1,6)
       do j = first,last
           xjval = zz(kuvpts(j))
           pj = 1.0
           do k = first,last
              xkval = zz(kuvpts(k))
              if (j .NE. k) then
                 pj = pj*(part%xp(3)-xkval)/(xjval-xkval)
              end if
           end do
           wt(3,j) = pj
      end do

       !Now compute weights at kwpts
       !Again must check to see how close we are to a top/bot boundary
       if (kwpts(3) == 0) then
          first = 3
          last = 4
          kwpts(1) = 1
          kwpts(2) = 1
       elseif (kwpts(3) .LT. 0) then 
          first = 0
          last = 0
          kwpts(1) = 1
          kwpts(2) = 1
          kwpts(3) = 1
       elseif (kwpts(3) == 1) then 
          first = 2
          last = 5
          kwpts(1) = 1
       elseif (kwpts(3) == nnz-1) then
          first = 3
          last = 4
          kwpts(5) = nnz
          kwpts(6) = nnz
       elseif (kwpts(3) .GE. nnz) then
          first = 0
          last = 0
          kwpts(3) = nnz
          kwpts(4) = nnz
          kwpts(5) = nnz
          kwpts(6) = nnz
       elseif (kwpts(3) == nnz-2) then
          first = 2
          last = 5
          kwpts(6) = nnz
       else
          first = 1
          last = 6
       end if

       !Recall that wt has been set to zero, so
       !weights will be zero if (first,last) isn't (1,6)
       do j = first,last
           xjval = z(kwpts(j))
           pj = 1.0
           do k = first,last
              xkval = z(kwpts(k))
              if (j .NE. k) then
                 pj = pj*(part%xp(3)-xkval)/(xjval-xkval)
              end if
           end do
           wt(4,j) = pj
      end do

      !Now we have the weights - compute the velocity at xp:
        part%uf(1:3) = 0.0
        part%Tf = 0.0
        part%qinf = 0.0
        do k = 1,6
        do j = 1,6
        do i = 1,6
            ix = ijpts(1,i)
            iy = ijpts(2,j)
            izuv = kuvpts(k)
            izw = kwpts(k)

            if (ix .gt. mxe+3) write(*,'(i,7e15.6)') myid,
     + part%xp(1),part%xp(2),part%xp(3),xmin,xmax,ymin,ymax
            part%uf(1) = part%uf(1)+uext(izuv,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(3,k) 
            part%uf(2) = part%uf(2)+vext(izuv,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(3,k) 
            part%uf(3) = part%uf(3)+wext(izw,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(4,k) 
            part%Tf = part%Tf+Text(izuv,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(3,k)
            part%qinf = part%qinf+T2ext(izuv,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(3,k) 
         end do
         end do 
         end do


      end subroutine uf_interp 

      subroutine uf_interp_lin
      use pars
      use particles
      use fields
      use con_stats
      use con_data
      implicit none
       
      integer :: ix,iy,izuv,izw,iz,i,k,j
      integer :: ipt,jpt,kpt,kwpt
      real :: wtx,wty,wtz,wtzw,wtt,wttw
      real :: xv,yv,zv,zwv
      

      ipt = floor(part%xp(1)/dx) + 1 
      jpt = floor(part%xp(2)/dy) + 1
     
     
      !Finding the k-lhnode is different since grid may be stretched
      !AND since (u,v) and w stored differently
      !Will get a k-index for (u,v) and one for w
      
      !Do (u,v) loop first:
      do iz=0,nnz+1
         if (zz(iz) .GT. part%xp(3)) go to 1199
      end do
1199  continue

      kpt = iz-1

      !Do again for w:
      do iz = 0,nnz+1 
         if (z(iz) .GT. part%xp(3)) go to 1120
      end do
1120  continue
   
      kwpt = iz-1


      part%uf(1:3) = 0.0
      part%Tf = 0.0
      part%qinf = 0.0
      do i=0,1
      do j=0,1
      do k=0,1

      xv = dx*(i+ipt-1)
      yv = dy*(j+jpt-1)
      zv = zz(k+kpt)
      zwv = z(k+kwpt)

      wtx = (1.0 - abs(part%xp(1)-xv)/dx)
      wty = (1.0 - abs(part%xp(2)-yv)/dy)
      wtz = (1.0 - abs(part%xp(3)-zv)/dzu(kpt+1))
      wtzw = (1.0 - abs(part%xp(3)-zwv)/dzw(kwpt+1))

      ix = ipt+i
      iy = jpt+i
      izuv = kpt+k
      izw = kwpt+k


      part%uf(1) = part%uf(1) + uext(izuv,iy,ix)*
     +      wtx*wty*wtz
      part%uf(2) = part%uf(2) + vext(izuv,iy,ix)*
     +      wtx*wty*wtz
      part%uf(3) = part%uf(3) + wext(izw,iy,ix)*
     +      wtx*wty*wtzw

      part%Tf = part%Tf + Text(izuv,iy,ix)*
     +      wtx*wty*wtz
      part%qinf = part%qinf + T2ext(izuv,iy,ix)*
     +      wtx*wty*wtz


      end do
      end do
      end do


      end subroutine uf_interp_lin

      subroutine sigm_interp(sigm_sdxp,sigm_sdyp,sigm_sdzp,
     +                          vis_sp,iz_part)
      use pars
      use particles
      use fields
      use con_stats
      use con_data
      implicit none

      integer :: ix,iy,izuv,izw,iz,i,k,j
      integer :: first,last
      real :: xkval,xjval,pj,dxvec(2)
      integer :: ijpts(2,6),kuvpts(6),kwpts(6),iz_part
      real :: wt(4,6)
      real :: ran2
      real :: sigm_sdxp,sigm_sdyp,sigm_sdzp,vis_sp
      !Following Orlando's code, get the "leftmost" node
      !This is just the minimum (i,j,k) on the volume

      ijpts(1,3) = floor(part%xp(1)/dx) + 1
      ijpts(2,3) = floor(part%xp(2)/dy) + 1

      !Fill in the neighbors:
      ijpts(1,2) = ijpts(1,3)-1
      ijpts(1,1) = ijpts(1,2)-1
      ijpts(1,4) = ijpts(1,3)+1
      ijpts(1,5) = ijpts(1,4)+1
      ijpts(1,6) = ijpts(1,5)+1

      ijpts(2,2) = ijpts(2,3)-1
      ijpts(2,1) = ijpts(2,2)-1
      ijpts(2,4) = ijpts(2,3)+1
      ijpts(2,5) = ijpts(2,4)+1
      ijpts(2,6) = ijpts(2,5)+1

      !Finding the k-lhnode is different since grid may be stretched
      !AND since (u,v) and w stored differently
      !Will get a k-index for (u,v) and one for w

      !Do (u,v) loop first:
      do iz=0,nnz+1
         if (zz(iz) .GT. part%xp(3)) go to 1199
      end do
1199  continue
c      iz_part = iz-1
      kuvpts(3) = iz-1
      !Then fill in the rest:
      kuvpts(4) = kuvpts(3)+1
      kuvpts(5) = kuvpts(4)+1
      kuvpts(6) = kuvpts(5)+1
      kuvpts(2) = kuvpts(3)-1
      kuvpts(1) = kuvpts(2)-1

      !Do again for w:
      do iz = 0,nnz+1
         if (z(iz) .GT. part%xp(3)) go to 1120
      end do
1120  continue
      iz_part = iz-1

      kwpts(3) = iz-1
      !Then fill in the rest:
      kwpts(4) = kwpts(3)+1
      kwpts(5) = kwpts(4)+1
      kwpts(6) = kwpts(5)+1
      kwpts(2) = kwpts(3)-1
      kwpts(1) = kwpts(2)-1

      !Fill in the weights:
      !First for x and y since they are periodic:
      wt(1:4,1:6) = 0.0
      dxvec(1) = dx
      dxvec(2) = dy
      do iz = 1,2
      do j = 1,6
         xjval = dxvec(iz)*(ijpts(iz,j)-1)
         pj = 1.0
         do k = 1,6
            xkval = dxvec(iz)*(ijpts(iz,k)-1)
            if (j .NE. k) then
                  pj = pj*(part%xp(iz)-xkval)/(xjval-xkval)
            end if
         end do
         wt(iz,j) = pj
       end do
       end do
       !Compute weights at kuvpts
       !Must check to see how close we are to a top/bot boundary
       if (kuvpts(3) == 1) then
          first = 3
          last = 4
          !Set these equal to 1 so uext(-1) won't be accessed
          !Note: the value doesn't matter since weight will be 0
          kuvpts(1) = 1
          kuvpts(2) = 1
       elseif (kuvpts(3) == 0) then
          first = 4
          last = 5
          kuvpts(1) = 1
          kuvpts(2) = 1
          kuvpts(3) = 1
       elseif (kuvpts(3) .LT. 0) then
          first = 0
          last = 0
       elseif (kuvpts(3) == 2) then
          first = 2
          last = 5
          kuvpts(1) = 1
          kuvpts(6) = 1
       !Between top cell center and the domain boundary
       elseif (kuvpts(3) == nnz) then
          first = 2
          last = 3
          kuvpts(4) = nnz
          kuvpts(5) = nnz
          kuvpts(6) = nnz
       elseif (kuvpts(3) .GT. nnz) then
          first = 0
          last = 0
       !Between 2nd to last and last cell center at top
       elseif (kuvpts(3) == nnz-1) then
          first = 3
          last = 4
          kuvpts(5) = nnz
          kuvpts(6) = nnz
       elseif (kuvpts(3) == nnz-2) then
          first = 2
          last = 5
          kuvpts(1) = nnz
          kuvpts(6) = nnz
       else
          first = 1
          last = 6
       end if

       !Recall that wt has been set to zero, so
       !weights will be zero if (first,last) isn't (1,6)
       do j = first,last
           xjval = zz(kuvpts(j))
           pj = 1.0
           do k = first,last
              xkval = zz(kuvpts(k))
              if (j .NE. k) then
                 pj = pj*(part%xp(3)-xkval)/(xjval-xkval)
              end if
           end do
           wt(3,j) = pj
      end do
       !Now compute weights at kwpts
       !Again must check to see how close we are to a top/bot boundary
       if (kwpts(3) == 0) then
          first = 3
          last = 4
          kwpts(1) = 1
          kwpts(2) = 1
       elseif (kwpts(3) .LT. 0) then
          first = 0
          last = 0
          kwpts(1) = 1
          kwpts(2) = 1
          kwpts(3) = 1
       elseif (kwpts(3) == 1) then
          first = 2
          last = 5
          kwpts(1) = 1
       elseif (kwpts(3) == nnz-1) then
          first = 3
          last = 4
          kwpts(5) = nnz
          kwpts(6) = nnz
       elseif (kwpts(3) .GE. nnz) then
          first = 0
          last = 0
          kwpts(3) = nnz
          kwpts(4) = nnz
          kwpts(5) = nnz
          kwpts(6) = nnz
       elseif (kwpts(3) == nnz-2) then
          first = 2
          last = 5
          kwpts(6) = nnz
       else
          first = 1
          last = 6
       end if
       !Recall that wt has been set to zero, so
       !weights will be zero if (first,last) isn't (1,6)
       do j = first,last
           xjval = z(kwpts(j))
           pj = 1.0
           do k = first,last
              xkval = z(kwpts(k))
              if (j .NE. k) then
                 pj = pj*(part%xp(3)-xkval)/(xjval-xkval)
              end if
           end do
           wt(4,j) = pj
      end do
      !Now we have the weights - compute the sigma_s and derivatives at
      !xp:
c        part%uf(1:3) = 0.0
c        part%Tf = 0.0
c        part%qinf = 0.0 ! Brian 7/14/14
        part%sigm_s = 0.0
        sigm_sdxp = 0.0
        sigm_sdyp = 0.0
        sigm_sdzp = 0.0
        vis_sp = 0.0
        do k = 1,6
        do j = 1,6
        do i = 1,6
            ix = ijpts(1,i)
            iy = ijpts(2,j)
            izuv = kuvpts(k)
            izw = kwpts(k)
            part%sigm_s = part%sigm_s+sigm_sext(izw,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(4,k)
            sigm_sdxp = sigm_sdxp+ sigm_sdxext(izw,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(4,k)
            sigm_sdyp = sigm_sdyp+sigm_sdyext(izw,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(4,k)
            sigm_sdzp = sigm_sdzp+sigm_sdzext(izuv,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(3,k)

            vis_sp =vis_sp +vis_sext(izw,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(4,k)

         if (isnan(vis_sp)) then
         write(*,*) 'P:',vis_sext(izw,iy,ix),wt(1,i),wt(2,j),wt(4,k),
     +   part%xp(3),i,j,k,kwpts(3)
         end if


         end do
         end do
         end do


      end subroutine sigm_interp

      subroutine particle_init
      use particles
      use pars
      use con_data
      implicit none
      include 'mpif.h' 
      integer :: values(8)
      integer :: idx,ierr
      integer*8 :: mult
      real :: xv,yv,zv,ran2,deltaz
      real :: maxx,maxy,maxz
      real :: xp_init(3),rad_init,m_s,Os
      real :: S,M,prob
      integer*8 :: num_a,num_c,totnum_a,totnum_c

      !Create the seed for the random number generator:
      call date_and_time(VALUES=values)
      iseed = -(myid+values(8)+values(7)+values(6))

      write(*,*) 'iseed=',myid,iseed

      !Initialize ngidx, the particle global index
      ngidx = 1
  
      !For the channel case, set the total number of particles:
      deltaz = zmax-zmin

      numpart = tnumpart/numprocs
      if (myid == 0) then
      numpart = numpart + MOD(tnumpart,numprocs)
      endif

      !Lognormal distribution parameters
      mult_factor = 500
      pdf_factor = 0.02*real(mult_factor)
      prob = pdf_factor/(1 + pdf_factor)

      !Adjust the multiplicity so that the total number of particles
      !isn't altered:
      mult_c = mult_init/mult_factor
      mult_a = mult_init/(1.0-prob)*(1 - prob/real(mult_factor))

      !Initialize the linked list of particles:
      nullify(part,first_particle)
      
      !Now initialize all particles with a random location on that processor
      maxx=0.0
      maxy=0.0
      maxz=0.0
      num_a = 0
      num_c = 0
      do idx=1,numpart
      xv = ran2(iseed)*(xmax-xmin) + xmin
      yv = ran2(iseed)*(ymax-ymin) + ymin
      zv = ran2(iseed)*(zi-zw1) + zw1
      xp_init = (/xv,yv,zv/)

      !Get the initial droplet size based on some distribution:

      !From params file -- all particles identical
      !rad_init = radius_init
      !m_s = Sal/1000.0*2.0/3.0*pi2*radius_init**3*rhow

      !call exponential_dist(rad_init)  !From Shima et al. (2009) test case
      !m_s = Sal/1000.0*4.0/3.0*pi*radius_init**3*rhow

      !Generate a double lognormal from individual coarse & accumulation modes
      !call double_lognormal_dist(rad_init,m_s,Os,mult)  !First attempt!

      !Generate
      if (ran2(iseed) .gt. prob) then   !It's accumulation mode
         S = 0.5
         M = -1.95
         Os = 0.5
         mult = mult_a

         !With these parameters, get m_s and rad_init from distribution
         call lognormal_dist(rad_init,m_s,Os,M,S)
         num_a = num_a + 1

       else  !It's coarse mode

         S = 0.45
         M = 0.0
         Os = 1.0
         mult = mult_c

         !With these parameters, get m_s and rad_init from distribution
         call lognormal_dist(rad_init,m_s,Os,M,S)
         num_c = num_c + 1

       end if


      call create_particle(xp_init,vp_init,
     +         Tp_init,m_s,Os,mult,rad_init,idx) 
      end do


      partTsrc = 0.0
      partTsrc_t = 0.0
      partHsrc = 0.0
      partHsrc_t = 0.0
      partTEsrc = 0.0
      partTEsrc_t = 0.0


      end subroutine particle_init

      subroutine exponential_dist(rad_init)
      use particles
      use con_data
      implicit none
      include 'mpif.h'

      real :: R0,X0,rad_init,vol_init
      real :: ran2

      !Exponential distribution of droplet size to match Shima et al.
      !(2009) test case
      R0 = 30.531e-6
      X0 = pi2*2.0/3.0*R0**3

      vol_init = -X0*log(1.0-ran2(iseed))
      rad_init = (vol_init/(pi2*2.0/3.0))**(1.0/3.0)

      end subroutine exponential_dist

      subroutine lognormal_dist(rad_init,m_s,Os,M,S)
      use particles
      use pars
      use con_data
      implicit none
      include 'mpif.h'

      real, intent(inout) :: rad_init,m_s,Os
      real :: ran2,cdf_func_single
      real :: M,S
      real :: d1,d2,err,dhalf,ftest,CDF
      real :: daerosol
      real :: a,b,c,d,p,q
      integer :: iter,prob

      !Use bisection to get the radius based on the CDF contained in function cdf_func
      d1 = 100.0
      d2 = 0.01
      err = 1.0
      iter = 0
      CDF = ran2(iseed)

      do while (err > 1.0e-10)

         dhalf = 0.5*(d1+d2)
         ftest = cdf_func_single(dhalf,CDF,M,S)
         if (ftest < 0.0) then
            d2 = dhalf
         else 
            d1 = dhalf
         end if

         err = abs(cdf_func_single(dhalf,CDF,M,S))
         iter = iter + 1

         if (iter .gt. 1000) then
            d1 = 1.0e-6
            write(*,'(a12)') 'CDF convergence error'
            exit
         end if


      end do
      
      daerosol = d1*1.0e-6  !Don't forget to convert micron to m
      m_s = 2.0/3.0*pi2*(daerosol/2.0)**3*rhos

      !Now have the dry aerosol diameter and mass, must rehydrate it using Kohler
      !theory to get the proper initial condition
      !Rehydrate to 95% RH

      !The equation for the equilibrium radius is a cubic -- take positive root:
      
      a = log(0.95)  !RH = 95%
      b = -2.0*Mw*Gam/Ru/rhow/Tp_init
      c = 0.0
      d = Ion*Os*m_s*(Mw/Ms)/(2.0/3.0*pi2*rhow)

      p = -b/3.0/a
      q = p**3 - 3.0*a*d/6.0/a**2

      rad_init = (q + (q**2 - p**6)**(0.5))**(1./3.) + (q -
     +         (q**2-p**6)**(0.5))**(1./3.) + p

      end subroutine lognormal_dist
      subroutine double_lognormal_dist(rad_init,m_s,Os,mult)
      use particles
      use pars
      use con_data
      implicit none
      include 'mpif.h'

      real, intent(inout) :: rad_init,m_s,Os
      real :: ran2,cdf_func
      real :: M_a,S_a,M_c,S_c,totarea
      real :: d1,d2,err,dhalf,ftest,CDF
      real :: daerosol
      real :: a,b,c,d,p,q
      integer :: iter,prob
      integer*8 :: mult

      !Set the parameters of the two lognormals:
      M_a = -1.95
      S_a = 0.5
      M_c = 0.0
      S_c = 0.45

      !The ratio of the coarse mode PDF to accumulation mode
      !NOTE: Making the factor inflated by 100 but decreasing the
      !multiplicity by same amount to get better coarse mode resolution

      totarea = 1.0 + pdf_factor  !The total area of combined PDFs, for normalization purposes 

      !Use bisection to get the radius based on the CDF contained in function cdf_func
      d1 = 100.0
      d2 = 0.01
      err = 1.0
      iter = 0
      CDF = ran2(iseed)

      do while (err > 1.0e-10)

         dhalf = 0.5*(d1+d2)
         ftest = cdf_func(dhalf,CDF,pdf_factor,M_a,S_a,M_c,S_c)
         if (ftest < 0.0) then
            d2 = dhalf
         else 
            d1 = dhalf
         end if

         err = abs(cdf_func(dhalf,CDF,pdf_factor,M_a,S_a,M_c,S_c))
         iter = iter + 1

         if (iter .gt. 1000) then
            d1 = 1.0e-6
            write(*,'(a12)') 'CDF convergence error'
            exit
         end if


      end do
      
      daerosol = d1*1.0e-6  !Don't forget to convert micron to m
      m_s = 2.0/3.0*pi2*(daerosol/2.0)**3*rhos

      !Actually set Os (to set kappa of kappa-Kohler theory) based on
      !size -- accumulation and coarse modes likely have different
      !compositions

      if (daerosol .lt. 1.0e-6) then
         Os = 0.5
         mult = mult_init
      elseif (daerosol .gt. 1.0e-6) then
         Os = 1.0
         mult = mult_init
      else
         if (ran2(iseed) .lt. 0.5) then
            prob = 0
         else
            prob = 1
         end if
        Os = real(prob)*0.5 + real((1-prob))*1.0
        mult = prob*mult_init + (1-prob)*mult_init
      end if
         

      !Now have the dry aerosol diameter and mass, must rehydrate it using Kohler
      !theory to get the proper initial condition
      !Rehydrate to 95% RH

      !The equation for the equilibrium radius is a cubic -- take positive root:
      
      a = log(0.95)  !RH = 95%
      b = -2.0*Mw*Gam/Ru/rhow/Tp_init
      c = 0.0
      d = Ion*Os*m_s*(Mw/Ms)/(2.0/3.0*pi2*rhow)

      p = -b/3.0/a
      q = p**3 - 3.0*a*d/6.0/a**2

      rad_init = (q + (q**2 - p**6)**(0.5))**(1./3.) + (q -
     +         (q**2-p**6)**(0.5))**(1./3.) + p

      end subroutine double_lognormal_dist


      subroutine particle_setup

      use particles
      use pars
      implicit none 
      include 'mpif.h'

      integer :: blcts(3),types(3)
      integer :: ierr
      real :: pi
      integer(kind=MPI_ADDRESS_KIND) :: extent,lb
      integer(kind=MPI_ADDRESS_KIND) :: extent2,lb2,displs(3)
      integer :: num_reals,num_integers,num_longs

      !First set up the neighbors for the interpolation stage:
      call assign_nbrs

      !Also assign the x,y,z max and mins to track particles leaving
      call set_bounds


      call set_binsdata(bins_rad,histbins+2,1.0e-8,1.0e-3)


      if (myid==0) write(*,*) 'Init Particle radius = ',radius_init
      if (myid==0) write(*,*) 'Init Particle Salinity = ',Sal     
      if (myid==0) write(*,*) 'Particle gravity = ',part_grav
      if (myid==0) write(*,*) 'nuf = ',nuf
      if (myid==0) write(*,*) 'CpaCpp = ',CpaCpp
      if (myid==0) write(*,*) 'Pra = ',Pra

      !Initialize the linked list of particles:
      nullify(part,first_particle)

      !Set up MPI datatypes for sending particle information
      !MUST UPDATE IF THINGS ARE ADDED/REMOVED FROM PARTICLE STRUCTURE

      num_reals = 6*3+13
      num_integers = 4
      num_longs = 3
      
      blcts(1:3) = (/num_integers,num_reals,num_longs/)
      displs(1) = 0
      types(1) = mpi_integer
      call mpi_type_get_extent(mpi_integer,lb,extent,ierr)
      
      !Displace num_integers*mpi_integer
      displs(2) = extent*num_integers
      types(2) = mpi_real8
      call mpi_type_get_extent(mpi_real8,lb,extent,ierr)
      !Displace num_reals*size of mpi_real8
      displs(3) = displs(2) + extent*num_reals
      types(3) = mpi_integer8

      !Now define the type:
      call mpi_type_create_struct(3,blcts,
     +            displs,types,particletype,ierr)


       call mpi_type_get_true_extent(particletype,lb2,extent2,ierr)
       call mpi_type_get_extent(particletype,lb2,extent,ierr)
       if (extent .NE. sizeof(part) ) then
          if (myid==0) then
          write(*,*) 'WARNING: extent of particletype not equal
     +                  to sizeof(part):'
          write(*,*) 'sizeof(part) = ', sizeof(part)
!          write(*,*) 'sizeof(part%pidx) = ', sizeof(part%pidx)
          write(*,*) 'mpi_type_get_true_extent(particletype) = ',extent2
          write(*,*) 'mpi_type_get_extent(particletype) = ',extent
          end if
       end if
      
      !Need to compute any padding which may exist in particle struct:
      pad_diff = extent-extent2 
      if (myid==0) then
      write(*,*) 'mpi_get_extent = ',extent
      write(*,*) 'mpi_get_true_extent = ',extent2
      write(*,*) 'sizeof(part) = ',sizeof(part)
      write(*,*) 'DIFF = ',pad_diff
      end if
      if (pad_diff .LT. 0) then
        write(*,*) 'WARNING: mpi_get_extent - mpi_get_true_extent LT 0!'
        call mpi_finalize(ierr)
        stop
      end if
      
      if (myid==0) then
      write(*,*) 'huge(tnumpart) = ',huge(tnumpart)
      write(*,*) 'huge(part%pidx) = ',huge(part%pidx)
      write(*,*) 'huge(part%mult) = ',huge(part%mult)
      end if


      call mpi_type_commit(particletype,ierr)

      end subroutine particle_setup
      subroutine assign_nbrs
        use pars
        use particles
        include 'mpif.h'
      !Figure out which processors lie to all sides: 
      !NOTE: For this updated case, where particles lie in columns not 
      !aligning with the velocity, there will be no MPI_PROC_NULL since
      !x and y are BOTH periodic
     
      !On right boundary:
      if ( mod(myid+1,ncpu_s) == 0 ) then
         !On the top:
         if ( myid .GE. ncpu_s*(ncpu_z-1) ) then
            rproc = myid-ncpu_s+1
            trproc = 0 
            tproc = ncpu_s-1 
            tlproc = ncpu_s-2 
            lproc = myid-1
            blproc = myid-ncpu_s-1
            bproc = myid-ncpu_s
            brproc = myid-ncpu_s - ncpu_s+1
         !On the bottom:
         elseif ( myid .LT. ncpu_s ) then
            rproc = myid-ncpu_s+1
            trproc = myid+1
            tproc = myid+ncpu_s
            tlproc = myid+ncpu_s-1
            lproc = myid-1
            blproc = myid+ncpu_s*(ncpu_z-1)-1 
            bproc = myid+ncpu_s*(ncpu_z-1) 
            brproc = ncpu_s*(ncpu_z-1) 
         !In the middle of right side:
         else 
            rproc = myid-ncpu_s+1
            trproc = myid+1
            tproc = myid+ncpu_s
            tlproc = myid+ncpu_s-1
            lproc = myid-1
            blproc = myid-ncpu_s-1
            bproc = myid-ncpu_s
            brproc = myid-ncpu_s - ncpu_s+1
         end if 

      !On the left boundary:
      elseif ( mod(myid,ncpu_s) == 0) then
         !On the top:
         if ( myid .GE. ncpu_s*(ncpu_z-1) ) then
            rproc = myid+1
            trproc = 1 
            tproc = 0 
            tlproc = ncpu_s-1
            lproc = myid+ncpu_s-1
            blproc = myid-1
            bproc = myid-ncpu_s
            brproc = myid-ncpu_s+1
         !On the bottom:
         elseif ( myid .LT. ncpu_s ) then
            rproc = myid+1
            trproc = myid+ncpu_s+1
            tproc = myid+ncpu_s
            tlproc = myid+ncpu_s+ncpu_s-1
            lproc = myid+ncpu_s-1
            blproc = numprocs-1 
            bproc = ncpu_s*(ncpu_z-1) 
            brproc = ncpu_s*(ncpu_z-1)+1 
         !In the middle of left side:
         else
            rproc = myid+1
            trproc = myid+ncpu_s+1
            tproc = myid+ncpu_s
            tlproc = myid+ncpu_s + ncpu_s-1
            lproc = myid+ncpu_s-1
            blproc = myid-1
            bproc = myid-ncpu_s
            brproc = myid-ncpu_s+1
         end if
      !On the top boundary
      elseif ( myid .GE. ncpu_s*(ncpu_z-1) ) then
         !Only check if in the middle:
         if ( .NOT. ( mod(myid,ncpu_s) == 0) ) then
            if ( .NOT. (mod(myid+1,ncpu_s) == 0) ) then
               rproc = myid+1
               trproc = myid-(ncpu_s*(ncpu_z-1))+1 
               tproc = myid-(ncpu_s*(ncpu_z-1)) 
               tlproc = myid-(ncpu_s*(ncpu_z-1))-1 
               lproc = myid-1
               blproc = myid-ncpu_s-1
               bproc = myid-ncpu_s
               brproc = myid-ncpu_s+1
            end if
         end if 
      !On the bottom boundary
      elseif ( myid .LT. ncpu_s) then
         if ( .NOT. ( mod(myid,ncpu_s) == 0) ) then
            if ( .NOT. (mod(myid+1,ncpu_s) == 0) ) then
               rproc = myid+1
               trproc = myid+ncpu_s+1
               tproc = myid+ncpu_s
               tlproc = myid+ncpu_s-1
               lproc = myid-1
               blproc = myid+ncpu_s*(ncpu_z-1)-1
               bproc = myid+ncpu_s*(ncpu_z-1) 
               brproc = myid+ncpu_s*(ncpu_z-1)+1 
            end if
         end if
      !Everywhere else:
      else 
         rproc = myid+1
         trproc = myid+ncpu_s+1
         tproc = myid+ncpu_s
         tlproc = myid+ncpu_s-1
         lproc = myid-1
         blproc = myid-ncpu_s-1
         bproc = myid-ncpu_s
         brproc = myid-ncpu_s+1
      end if

      return
      end subroutine assign_nbrs
      subroutine set_bounds  
        use pars
        use particles
        use con_data
        use con_stats
        implicit none
        include 'mpif.h'

      !Each processor must figure out at what ymin,ymax,zmin,zmax a particle leaves
      ymin = dy*(iys-1)
      ymax = dy*(iye)
      zmin = z(izs-1)
      zmax = z(ize)  
      xmin = dx*(mxs-1)
      xmax = dx*(mxe)

      end subroutine set_bounds
      subroutine save_particles
      use particles
      use pars
      implicit none
      include 'mpif.h'

      integer :: istatus(mpi_status_size), ierr, fh
      integer(kind=mpi_offset_kind) :: zoffset,offset
      integer :: pnum_vec(numprocs)
      integer :: iproc,i
      type(particle) :: writebuf(numpart),tmp

      !Do this with mpi_write_at_all
      !Need to figure out the displacements - need numpart from each proc
      call mpi_allgather(numpart,1,mpi_integer,pnum_vec,1,mpi_integer,
     +                   mpi_comm_world,ierr)

      !Package all the particles into writebuf:
      i = 1
      part => first_particle
      do while (associated(part))
      writebuf(i) = part
      !write(*,'a5,3e15.6') 'xp:',part%xp(1:3)
      part => part%next
      i = i + 1
      end do

      !Now only write to the file if you actually have particles
      !EXCEPTION: proc 0, which needs to write tnumpart regardless
      call mpi_file_open(mpi_comm_world, path_sav_part,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null,fh,ierr)

      zoffset = 0
      !Write tnumpart first:
      if (myid==0) then
      call mpi_file_write_at(fh,zoffset,tnumpart,1,mpi_integer,
     +                       istatus,ierr)
      write(*,*) 'wrote tnumpart = ',tnumpart
      end if

      zoffset = zoffset + 4
     
      !Now compute the offset (in bytes!):
      offset = zoffset 
      do iproc = 0,myid-1
         offset = offset + pnum_vec(iproc+1)*(sizeof(tmp)-pad_diff) 
      end do

      !Now everyone else write, ONLY if numpart > 0
      if (numpart .GT. 0) then
      call mpi_file_write_at(fh,offset,writebuf,numpart,
     +                       particletype,istatus,ierr)
!      write(*,*) 'wrote writebuf(1)%pidx = ',writebuf(1)%pidx
!      write(*,'a30,3e15.6') 'wrote writebuf(1)%vp =',writebuf(1)%vp(1:3)
!      write(*,'a30,3e15.6') 'wrote writebuf(1)%xp =',writebuf(1)%xp(1:3)
!      write(*,*) 'wrote writebuf(2)%pidx = ',writebuf(2)%pidx
!      write(*,'a30,3e15.6') 'wrote writebuf(2)%vp =',writebuf(2)%vp(1:3)
!      write(*,'a30,3e15.6') 'wrote writebuf(2)%xp =',writebuf(2)%xp(1:3)
      end if

      call mpi_file_close(fh,ierr)

      write(*,*) 'proc',myid,'wrote numpart = ',numpart

      if (myid==0) write(*,7000) path_sav_part
 7000 format(' PARTICLE DATA IS WRITTEN IN FILE  ',a80)

      end subroutine save_particles
      subroutine read_part_res
      use pars
      use particles
      implicit none
      include 'mpif.h'

      integer :: istatus(mpi_status_size), ierr, fh
      integer(kind=mpi_offset_kind) :: zoffset,offset
      integer :: myp,totalp 
      integer :: iproc,i,pidxmax
      type(particle), allocatable :: readbuf(:)

      if (myid==0) write(*,7000) path_part 
 7000 format(' READING PARTICLE DATA FROM  ',a80)


      call mpi_file_open(mpi_comm_world,path_part,
     +                   mpi_mode_rdonly,
     +                   mpi_info_null,fh,ierr)


      !Read in the total number of particles:
      offset = 0
      call mpi_file_read_at_all(fh,offset,tnumpart,1,
     +                      mpi_integer,istatus,ierr)
      if (myid==0) write(*,*) 'read tnumpart = ',tnumpart
    
      offset = 4
      allocate(readbuf(tnumpart))
      call mpi_file_read_at_all(fh,offset,readbuf,tnumpart,
     +                       particletype,istatus,ierr)

      do i = 1,tnumpart
        !Now - does it lie within this proc's bounds?
        if (readbuf(i)%xp(2) .GT. ymin .AND.
     +       readbuf(i)%xp(2) .LT. ymax .AND.
     +       readbuf(i)%xp(1) .GT. xmin .AND.
     +       readbuf(i)%xp(1) .LT. xmax) then 
            if (.NOT. associated(first_particle)) then
               allocate(first_particle)
               first_particle = readbuf(i)
               nullify(first_particle%prev,first_particle%next)
               part => first_particle
            else
               allocate(part%next)
               part%next = readbuf(i)
               part%next%prev => part
               part => part%next
               nullify(part%next)
            end if

        end if
      end do

      deallocate(readbuf)

      call mpi_file_close(fh,ierr)
      
      !Now just check how many each processor obtained:
      !At the same time, figure out max(pidx) and set ngidx 
      !to one plus this value:
      pidxmax = 0
      part => first_particle
      myp = 0
      do while (associated(part))
         myp = myp+1
         if (part%pidx .gt. pidxmax) pidxmax = part%pidx
         part => part%next
      end do

      !Set ngidx (the index for creating new particles) to 1+pidmax:
      ngidx = pidxmax + 1

      numpart = myp
     
      call mpi_allreduce(myp,totalp,1,mpi_integer,mpi_sum,
     +                   mpi_comm_world,ierr)

      write(*,*) 'proc',myid,'read in numpart:',myp
      if (myid==0) write(*,*) 'total number of particles read:',totalp

      end subroutine read_part_res
      subroutine particle_stats(ipt,jpt,kpt)
      use particles
      use pars
      implicit none
      integer :: i,ipt,jpt,kpt
      real :: rhop,pi


      pi   = 4.0*atan(1.0)

      rhop = (part%m_s+4.0/3.0*pi*part%radius**3*rhow)/
     +                           (4.0/3.0*pi*part%radius**3)


      !Takes in ipt,jpt,kpt as the node to the "bottom left" of the particle
      !(i.e. the node in the negative direction for x,y,z)
      !and computes quantities needed to get particle statistics

      partcount_t(kpt,jpt,ipt) = partcount_t(kpt,jpt,ipt) + 1.0
      
      !Get su mean, mean-squared of particle velocities at each level
      upwp_t(kpt,jpt,ipt) = upwp_t(kpt,jpt,ipt) + part%vp(1)*part%vp(3)
      do i = 1,3
      vpsum_t(kpt,jpt,ipt,i) = vpsum_t(kpt,jpt,ipt,i) + part%vp(i)
      vpsqrsum_t(kpt,jpt,ipt,i)=vpsqrsum_t(kpt,jpt,ipt,i)+part%vp(i)**2
      end do

      Tpsum_t(kpt,jpt,ipt) = Tpsum_t(kpt,jpt,ipt) + part%Tp
      Tpsqrsum_t(kpt,jpt,ipt) = Tpsqrsum_t(kpt,jpt,ipt) + part%Tp**2

      Tfsum_t(kpt,jpt,ipt) = Tfsum_t(kpt,jpt,ipt) + part%Tf

      qfsum_t(kpt,jpt,ipt) = qfsum_t(kpt,jpt,ipt) + part%qinf

      wpTpsum_t(kpt,jpt,ipt) = wpTpsum_t(kpt,jpt,ipt) + 
     +                          part%Tp*part%vp(3)


      radsum_t(kpt,jpt,ipt) = radsum_t(kpt,jpt,ipt) + part%radius 

      rad2sum_t(kpt,jpt,ipt) = rad2sum_t(kpt,jpt,ipt) + part%radius**2  

      multcount_t(kpt,jpt,ipt) = multcount_t(kpt,jpt,ipt) + 
     +               real(part%mult)

      mwsum_t(kpt,jpt,ipt) = mwsum_t(kpt,jpt,ipt) + 
     +               real(part%mult)*(rhow*4.0/3.0*pi*part%radius**3)

      qstarsum_t(kpt,jpt,ipt) = qstarsum_t(kpt,jpt,ipt) + part%qstar



      end subroutine particle_stats
      subroutine particle_coupling_update
      use particles
      use pars
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'
      real :: wtx,wty,wtz,wtt,dV
      real :: rhop,taup_i,partmass
      real :: xv,yv,zv
      real :: ctbuf_s(nnz+2,1:iye-iys+2,6),cbbuf_r(nnz+2,1:iye-iys+2,6)
      real :: crbuf_s(nnz+2,1:mxe-mxs+1,6),clbuf_r(nnz+2,1:mxe-mxs+1,6)
      integer :: i,j,k,ncount,ipt,jpt,kpt,kwpt
      integer :: istatus(mpi_status_size),ierr
      integer :: ix,iy,iz
      real :: lref,Cvv,Cva,Cpv

      partsrc_t = 0.0
      partTsrc_t = 0.0
      partHsrc_t = 0.0
      partTEsrc_t = 0.0


      part => first_particle
      do while (associated(part))

      !First, as done in uf_interp, must find the "leftmost" node
      !of volume where particle belongs:
      !(must repeat since now particle locations have been updated)
     
      ipt = floor(part%xp(1)/dx) + 1
      jpt = floor(part%xp(2)/dy) + 1
      do iz=0,nnz+1
         if (zz(iz) .GT. part%xp(3)) go to 1299
      end do
1299  continue
      kpt = iz-1

      !Do again for w:
      do iz = 0,nnz+1
         if (z(iz) .GT. part%xp(3)) go to 1120
      end do
1120  continue

      kwpt = iz


!---------
      !As an aside, use ipt,jpt,kwpt to update the particle numbers at each cell:

      !call inputs are the nodes associated with the volume that each part. lies in
      call particle_stats(ipt,jpt,kwpt)
  
!---------

      !Add contribution to each of the 8 surrounding nodes:
      do i=0,1
      do j=0,1
      do k=0,1

      xv = dx*(i+ipt-1)
      yv = dy*(j+jpt-1)
      zv = zz(k+kpt)

      dV = dx*dy*dzu(kpt+1)

      wtx = (1.0 - abs(part%xp(1)-xv)/dx)
      wty = (1.0 - abs(part%xp(2)-yv)/dy)
      wtz = (1.0 - abs(part%xp(3)-zv)/dzu(kpt+1))
      wtt = wtx*wty*wtz

      rhop = (part%m_s+pi2*2.0/3.0*part%radius**3*rhow)/
     +                           (pi2*2.0/3.0*part%radius**3)
      partmass = rhop*2.0/3.0*pi2*(part%radius)**3
      taup_i = 18.0*rhoa*nuf/rhop/(2.0*part%radius)**2 !Brian 8/4/14

      !write(*,'a6,3i') 'node:',ipt+i,jpt+j,kpt+k
      !write(*,'a6,5e15.6') 'Vwt: ',dV,wtx,wty,wtz,wtt

      ix = ipt+i
!      !Account for periodicity:
!      if (ix == nnx+1) ix = 1
      iy = jpt+j
      iz = kpt+k

      if (ix .gt. mxe+1) write(*,*) 'proc',myid,'has ix = ',ix
      if (ix .lt. mxs) write(*,*) 'proc',myid,'has ix = ',ix
      if (iy .gt. iye+1) write(*,*) 'proc',myid,'has iy = ',iy
      if (iy .lt. iys) write(*,*) 'proc',myid,'has iy = ',iy
      if (iz .gt. nnz+1) write(*,*) 'proc',myid,'has iz = ',iz
      if (iz .lt. 0) then
          write(*,*) 'proc',myid,'has iz = ',iz
          write(*,*) 'DHR:',part%radius,part%xp(3),part%mult,part%vp(3)
      end if

      !Recall to subtract g since momentum is extracted form
      !fluid only through drag term - NOT the gravity term as well
      if (icouple == 1) then
      !drag momentum coupling
      partsrc_t(iz,iy,ix,1:3) =
     +     partsrc_t(iz,iy,ix,1:3) +
     +     partmass/rhoa*(part%vrhs(1:3)-part_grav(1:3))*wtt/dV*
     +     real(part%mult)

      !vapor momentum coupling
      partsrc_t(iz,iy,ix,1:3) =
     +     partsrc_t(iz,iy,ix,1:3) +
     +     rhow/rhoa*pi2*2*part%radius**2*
     +     part%radrhs*part%vp(1:3)*wtt/dV*
     +     real(part%mult)
      endif

      if (iTcouple == 1) then
      partTsrc_t(iz,iy,ix) =
     +     partTsrc_t(iz,iy,ix) +
     +     (part%Tprhs_s*6.0*rhow/rhop/CpaCpp/taup_i*
     +     (pi2/2.0)*part%radius*nuf)*wtt/dV*
     +     real(part%mult)
      endif

      if (iHcouple == 1) then
      partHsrc_t(iz,iy,ix) =     
     +     partHsrc_t(iz,iy,ix) +  !minus comes later
     +     rhow/rhoa*pi2*2*part%radius**2*part%radrhs*wtt/dV*
     +     real(part%mult)


       Cvv = 1463.0 !J/degkg
       Cpv = 1952.0 !from E. Toolbox -> 1864.0 !J/degkg
       !Cva = 717.0 !J/degkg
       partTEsrc_t(iz,iy,ix) = 
     +     partTEsrc_t(iz,iy,ix) +  !minus comes later
     +         rhow/rhoa*pi2*2*part%radius**2*
     +         part%radrhs*Cpv/Cpa*part%Tp*wtt/dV*real(part%mult) -
     +         rhow/rhoa*pi2*2*part%radius**2*
     +         part%radrhs*Cpv/Cpa*part%Tf*wtt/dV*real(part%mult)


      endif


      end do
      end do
      end do

      part => part%next
      end do


      end subroutine particle_coupling_update
      subroutine particle_coupling_exchange
      use particles
      use pars
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'
      real :: ctbuf_s(nnz+2,1:iye-iys+2,6),cbbuf_r(nnz+2,1:iye-iys+2,6)
      real :: crbuf_s(nnz+2,1:mxe-mxs+1,6),clbuf_r(nnz+2,1:mxe-mxs+1,6)
      integer :: istatus(mpi_status_size),ierr,ncount


      !Now, partsrc and partTsrc have halos on each processor - give
      !these to the rightful owner:
      crbuf_s=0.0;ctbuf_s=0.0
      clbuf_r=0.0;cbbuf_r=0.0

      !First send top: 
      !get send buffer ready:
      ctbuf_s(1:nnz+2,1:iye-iys+2,1:3)=partsrc_t(0:nnz+1,
     +                                iys:iye+1,mxe+1,1:3)
      ctbuf_s(1:nnz+2,1:iye-iys+2,4)=partTsrc_t(0:nnz+1,
     +                                iys:iye+1,mxe+1)

      ctbuf_s(1:nnz+2,1:iye-iys+2,5)=partHsrc_t(0:nnz+1,
     +                                iys:iye+1,mxe+1)

       ctbuf_s(1:nnz+2,1:iye-iys+2,6)=partTEsrc_t(0:nnz+1,
     +                                iys:iye+1,mxe+1)

      ncount = 6*(nnz+2)*(iye-iys+2)
      call mpi_sendrecv(ctbuf_s,ncount,mpi_real8,tproc,1,
     +     cbbuf_r,ncount,mpi_real8,bproc,1,mpi_comm_world,istatus,ierr)

      !Now just add the contents of the receive buffer into the entire 
      !iys column of this proc:

      partsrc_t(0:nnz+1,iys:iye+1,mxs,1:3) =
     +   partsrc_t(0:nnz+1,iys:iye+1,mxs,1:3) +
     +    cbbuf_r(1:nnz+2,1:iye-iys+2,1:3)
      partTsrc_t(0:nnz+1,iys:iye+1,mxs) =
     +   partTsrc_t(0:nnz+1,iys:iye+1,mxs) +
     +     cbbuf_r(1:nnz+2,1:iye-iys+2,4)

      partHsrc_t(0:nnz+1,iys:iye+1,mxs) =
     +   partHsrc_t(0:nnz+1,iys:iye+1,mxs) +
     +     cbbuf_r(1:nnz+2,1:iye-iys+2,5)


      partTEsrc_t(0:nnz+1,iys:iye+1,mxs) =
     +   partTEsrc_t(0:nnz+1,iys:iye+1,mxs) +
     +   cbbuf_r(1:nnz+2,1:iye-iys+2,6)

      !Now get the right send buffer ready:
      crbuf_s(1:nnz+2,1:mxe-mxs+1,1:3)=
     +        partsrc_t(0:nnz+1,iye+1,mxs:mxe,1:3)
      crbuf_s(1:nnz+2,1:mxe-mxs+1,4)=
     +        partTsrc_t(0:nnz+1,iye+1,mxs:mxe)
      crbuf_s(1:nnz+2,1:mxe-mxs+1,5)=
     +        partHsrc_t(0:nnz+1,iye+1,mxs:mxe)

      crbuf_s(1:nnz+2,1:mxe-mxs+1,6)=
     +        partTEsrc_t(0:nnz+1,iye+1,mxs:mxe)

      !Now send to right:
      ncount = 6*(nnz+2)*(mxe-mxs+1)
      call mpi_sendrecv(crbuf_s,ncount,mpi_real8,rproc,2,
     +     clbuf_r,ncount,mpi_real8,lproc,2,mpi_comm_world,istatus,ierr)

      !And again add the contents to the top/bottom rows of partsrc:
      partsrc_t(0:nnz+1,iys,mxs:mxe,1:3) =
     +           partsrc_t(0:nnz+1,iys,mxs:mxe,1:3) +
     +                clbuf_r(1:nnz+2,1:mxe-mxs+1,1:3)

      partTsrc_t(0:nnz+1,iys,mxs:mxe) =
     +           partTsrc_t(0:nnz+1,iys,mxs:mxe) +
     +                 clbuf_r(1:nnz+2,1:mxe-mxs+1,4)
      partHsrc_t(0:nnz+1,iys,mxs:mxe) =
     +           partHsrc_t(0:nnz+1,iys,mxs:mxe) +
     +                 clbuf_r(1:nnz+2,1:mxe-mxs+1,5)
      partTEsrc_t(0:nnz+1,iys,mxs:mxe) =
     +           partTEsrc_t(0:nnz+1,iys,mxs:mxe) +
     +                 clbuf_r(1:nnz+2,1:mxe-mxs+1,6)


      end subroutine particle_coupling_exchange
      subroutine write_histograms
      use pars
      use fields
      use con_data
      use con_stats
      use particles
      implicit none

      include 'mpif.h'
      integer :: i,j
      integer :: ierr
      logical :: isfirst

      integer :: iblnk
      integer :: num_entries

      real :: sumbuf_rad(histbins+2)
      real :: binval_rad(histbins+2)

      character :: fformat*10, num*3


c -------------- record bins for histogram file -----------
        binval_rad = bins_rad(:)


c -------------- Collect all histogram from different processors
        num_entries = (histbins+2)
        call mpi_reduce(hist_rad,sumbuf_rad,num_entries,
     +   mpi_real8,mpi_sum,0,mpi_comm_world,ierr)

        if (myid==0) then


c ---------------- save data 

        write(num,'(i3.3)') histbins+2
        fformat = '('//num//'e15.6)'

        write(nrad,fformat),sumbuf_rad(1:histbins+2)

        end if

      end subroutine write_histograms
      subroutine set_binsdata(binsdata,sizea,rmin,rmax)
      use pars
      use fields
      use con_data
      use con_stats
      use particles
      implicit none

      integer :: i,nbin,ibin,nbinnew
      real :: dh

      integer,intent(in) :: sizea
      real,intent(inout) :: binsdata(sizea)

      real :: rmin,rmax,rmin10,rmax10

        nbin = histbins !From Module Particle 

        rmin10 = log10(rmin)
        rmax10 = log10(rmax)

c       Calculate size of interval
        dh = (rmax10-rmin10)/(nbin-1)

c       ===== update x-axis for each bin ===== 
        binsdata(1) = rmin10-dh
        do i = 1,histbins+1
          binsdata(i+1)= dh+binsdata(i)
        end do


      end subroutine set_binsdata
      subroutine add_histogram(binsdata,histdata,sizea,val1,mult)

      use pars
      use fields
      use con_data
      use con_stats
      use particles
      implicit none

      integer :: i,nbin,nbinnew,ibin
      integer*8, intent(in) :: mult
      real :: dh
      real,intent(in) :: val1

      integer,intent(in) :: sizea
      real,intent(in) :: binsdata(sizea)
      real,intent(inout) :: histdata(sizea)

      real :: rmin,rmax,logval1

        rmin = binsdata(2)
        nbin = histbins !From Module Particle 
        rmax = binsdata(nbin+1)

c       Calculate size of interval
        dh = (rmax-rmin)/(nbin-1)

        logval1 = log10(val1)
        if (logval1 .gt. rmax+0.5*dh) then
                ibin = nbin + 2
        elseif (logval1 .lt. rmin-0.5*dh) then
                ibin = 1
        else
                ibin = (floor((logval1-(rmin-0.5*dh))/dh)+1)+1
        end if


c       Add the current event to the histogram
        histdata(ibin) = histdata(ibin) + real(mult)

      end subroutine add_histogram
      subroutine radius_histogram
      use particles
      implicit none

      hist_rad = 0.0
      part => first_particle
      do while (associated(part))
         !if (part%xp(3) .lt. 60.0) then
         call add_histogram(bins_rad,hist_rad,
     +    histbins+2,part%radius,part%mult)
         !end if
      part => part%next
      end do

      end subroutine radius_histogram

      subroutine read_input_file
      use pars
      use particles
      use con_data
      implicit none

      character(48) :: label
      character(180) :: params_dir
      namelist /step_params/ iti,itmax,imean,ihst,itape,
     +                       itstr,it_his,it_viz,i_viz,itn

      namelist /grid_params/ ncpu_s, Uo, Ttop, Tbot,
     +         qstar, tsfcc, ugcont, vgcont,
     +         zi, zl, xl, yl, zw1, dpdx,
     +         cfl,dt_new,surf_RH,psurf,surf_temp,
     +         rad_kappa,rad_Fct,rad_Fcb

      namelist /path_names/ path_seed,path_part,path_res,
     +         path_sav,path_his,path_viz_xy,path_viz_xz,
     +         path_viz_yz,path_stuf,path_ran,path_histog

      namelist /flags/ ismlt,ifree,isfc,iradup,
     +         iupwnd,ibuoy,ifilt,itcut,isubs,ibrcl,iocean,
     +         method,idebug,iz_space,ivis0,ifix_dt,new_vis,iDNS,
     +         icouple,iTcouple,iHcouple,ievap,ifields,ilin,
     +         ineighbor,icoalesce,ipart_method,
     +         ipartdiff,isfs,iexner,ilongwave


      namelist /constants/ rhoa, nuf, Cpa, Pra, Sc,
     +         tnumpart,mult_init,rhow,part_grav,rhos,
     +         Cpp,Mw,Ms,Ru,Gam,Ion,Sal,Rd,
     +         radius_init,
     +         grav, t00,fcor,zo,
     +         vp_init,Tp_init


      !params.in contains namelists to read
      !open(12, file="./params.in", status="old")
      call get_command_argument(1,params_dir)
      open(12,file=params_dir,status="old")

      read(12,nml=step_params)
      if (myid==0) print step_params

      read(12,nml=flags)
      if (myid==0) print flags

      read(12,nml=grid_params)
      if (myid==0) print grid_params

      read(12,nml=path_names)
      if (myid==0) print path_names

      read(12,nml=constants)
      if (myid==0) print constants
      CpaCpp = Cpa/Cpp



      end subroutine read_input_file

      function ran2(idum)
      integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real :: ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     +     IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     +     IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER :: idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/,iv/NTAB*0/,iy/0/

      if (idum .le. 0) then
          idum=max(-idum,1)
          idum2 = idum
          do j = NTAB+8,1,-1
             k=idum/IQ1
             idum=IA1*(idum-k*IQ1)-k*IR1
             if (idum .lt. 0) idum=idum+IM1
             if (j .le. NTAB) iv(j) = idum
          end do
          iy=iv(1)
      end if
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum .lt. 0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2 .lt. 0) idum2=idum2+IM2
      j = 1+iy/NDIV
      iy = iv(j) - idum2
      iv(j) = idum
      if (iy .lt. 1) iy = iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      end function ran2
      function gasdev(idum)
       integer :: idum
       real :: gasdev
       integer :: iset
       real :: fac,gset,rsq,v1,v2,ran2
       save iset,gset
       data iset/0/

       if (idum .lt. 0) iset=0
       if (iset .eq. 0) then
 1000    v1 = 2.*ran2(idum)-1.0
         v2 = 2.*ran2(idum)-1.0
         rsq = v1**2+v2**2
         if ( (rsq .ge. 1) .or. (rsq .eq. 0)) goto 1000
         fac = sqrt(-2.0*log(rsq)/rsq)
         gset = v1*fac
         gasdev = v2*fac
         iset = 1
       else
         gasdev = gset
         iset = 0
       end if
       return
      end
      function mod_Magnus(T)
      implicit none

      !Take in T in Kelvin and return saturation vapor pressure using
      !function of Alduchov and Eskridge, 1996
      real,intent(in) :: T  
      real :: mod_Magnus

      
      mod_Magnus = 610.94 *exp((17.6257*(T-273.15))/
     +  (243.04+(T-273.15)))
       
      
      end function mod_Magnus
      function cdf_func_single(d,CDF,M,S)
      implicit none
      real :: d,CDF,M,S,cdf_func_single

      cdf_func_single = 1.0/2.0*
     + (1.0 + erf((log(d) - M)/S/sqrt(2.0)))  - CDF

      end function cdf_func_single
      function cdf_func(d,CDF,factor,Ma,Sa,Mc,Sc)
      implicit none
      real :: d,CDF,factor,Ma,Sa,Mc,Sc,cdf_func

      cdf_func = 1.0/(1.0+factor)/2.0*
     + ( 1.0 + erf((log(d) - Ma)/Sa/sqrt(2.0)) +
     + factor*(1 + erf((log(d) - Mc)/Sc/sqrt(2.0))) ) - CDF

      end function cdf_func
