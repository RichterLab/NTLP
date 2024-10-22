module netcdf_io
implicit none


integer :: ncid,ncid_histog,ncid_viz
integer :: time_dimid,zu_vid,zu_dimid,s_dimid
integer :: time_histog_dimid
integer :: histbins_dimid
integer :: time_vid,dt_vid,time_histog_vid
integer :: utau_vid,uwsfc_vid
integer :: Tsfc_vid,qsfc_vid,wtsfc_vid,wqsfc_vid
integer :: tnumpart_vid,tnum_destroy_vid,tot_reintro_vid
integer :: tnumdrop_vid,tnumaerosol_vid
integer :: Swall_vid
integer :: meanRH_vid,varRH_vid
integer :: radavg_vid,radmsqr_vid
integer :: radavg_center_vid,radmsqr_center_vid
integer :: tdenum_vid,tactnum_vid,tnumimpos_vid,tnum100_vid
integer :: zw_vid,zw_dimid
integer :: uxym_vid,vxym_vid,wxym_vid,txym_vid,RHxym_vid,tempxym_vid,exym_vid,RHmsqr_vid
integer :: ups_vid,vps_vid,wps_vid,tps_vid
integer :: wtle_vid,wtsb_vid
integer :: uwle_vid,uwsb_vid
integer :: vwle_vid,vwsb_vid
integer :: t_diss_vid,t_rprod_vid,t_sprod_vid,t_tran_vid,t_buoy_vid
integer :: zconc_vid
integer :: pflux_vid,pfluxdiff_vid
integer :: pmassflux_vid,penegflux_vid
integer :: vp1mean_vid,vp2mean_vid,vp3mean_vid
integer :: vp1msqr_vid,vp2msqr_vid,vp3msqr_vid
integer :: uf1mean_vid,uf2mean_vid,uf3mean_vid
integer :: uf1msqr_vid,uf2msqr_vid,uf3msqr_vid
integer :: m1src_vid,m2src_vid,m3src_vid
integer :: Tpsrc_vid,TEpsrc_vid,Hpsrc_vid
integer :: Tpmean_vid,Tpmsqr_vid
integer :: Tfmean_vid,qfmean_vid
integer :: radmean_vid,rad2mean_vid
integer :: qstarm_vid
integer :: Nc_vid,ql_vid,radsrc_vid
integer :: rho_base_vid,p_base_vid,T_base_vid,theta_base_vid
integer :: radbins_vid,resbins_vid
integer :: radhist_vid,reshist_vid,radhistdeath_vid
integer :: actresbins_vid,actreshist_vid
integer :: acttodeathbins_vid,acttodeathhist_vid
integer :: numactbins_vid,numacthist_vid

!Viz:
integer :: time_viz_dimid,viz_nx_dimid,viz_ny_dimid,viz_nzu_dimid,viz_nzw_dimid
integer :: time_viz_vid
integer :: u_yz_vid,t_yz_vid,q_yz_vid,w_yz_vid,v_yz_vid
integer :: u_xy_vid,t_xy_vid,q_xy_vid,w_xy_vid,v_xy_vid
integer :: u_xz_vid,t_xz_vid,q_xz_vid,w_xz_vid,v_xz_vid
integer :: xgrid_vid,ygrid_vid,zugrid_vid,zwgrid_vid
integer :: LWP_vid,surf_precip_vid

character(len=80) :: path_netcdf_his,path_netcdf_histog,path_netcdf_viz

CONTAINS

subroutine netcdf_check(status)
      use netcdf
      implicit none
      integer, intent ( in) :: status

      if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
      end if

end subroutine netcdf_check

subroutine netcdf_init
      use netcdf
      use pars
      use con_data
      implicit none

      integer :: dimids(1),dimids_zu(2),dimids_zw(2),dimids_zu_s(3),dimids_zw_s(3)
      integer :: dimids_zu_only(1),dimids_zw_only(1)

      path_netcdf_his = trim(adjustl(path_seed))//"history.nc"

      call netcdf_check( nf90_create(path_netcdf_his,nf90_clobber,ncid))

      call netcdf_check( nf90_def_dim(ncid, "time",NF90_UNLIMITED, time_dimid) )

      call netcdf_check( nf90_def_dim(ncid, "zu",nnz, zu_dimid) )
      call netcdf_check( nf90_def_dim(ncid, "zw",nnz+1, zw_dimid) )

      call netcdf_check( nf90_def_dim(ncid,"nscl",nscl,s_dimid) )

      dimids = (/ time_dimid /)
      dimids_zu = (/ zu_dimid, time_dimid/)
      dimids_zw = (/ zw_dimid, time_dimid/)
      dimids_zu_only = (/ zu_dimid/)
      dimids_zw_only = (/ zw_dimid/)
      dimids_zu_s = (/ zu_dimid, s_dimid, time_dimid/)
      dimids_zw_s = (/ zw_dimid, s_dimid, time_dimid/)

!!! Single quantities
      call netcdf_check( nf90_def_var(ncid,"time",NF90_REAL,dimids,time_vid) )
      call netcdf_check( nf90_put_att(ncid,time_vid,"title","Simulation time") )

      call netcdf_check( nf90_def_var(ncid, "dt", NF90_REAL, dimids,dt_vid) )
      call netcdf_check( nf90_put_att(ncid,dt_vid,"title","Model time step") )

      call netcdf_check( nf90_def_var(ncid, "utau", NF90_REAL, dimids,utau_vid) )
      call netcdf_check( nf90_put_att(ncid,utau_vid,"title","Friction velocity ustar") )

      call netcdf_check( nf90_def_var(ncid, "uwsfc", NF90_REAL, dimids,uwsfc_vid) )
      call netcdf_check( nf90_put_att(ncid,uwsfc_vid,"title","Lower surface stress") )

      call netcdf_check( nf90_def_var(ncid, "tnumpart", NF90_REAL, dimids,tnumpart_vid) )
      call netcdf_check( nf90_put_att(ncid,tnumpart_vid,"title","Total number of particles") )

      call netcdf_check( nf90_def_var(ncid, "tnumdrop", NF90_REAL, dimids,tnumdrop_vid) )
      call netcdf_check( nf90_put_att(ncid,tnumdrop_vid,"title","Total number of droplets defined by r>rc") )

      call netcdf_check( nf90_def_var(ncid, "tnumaerosol", NF90_REAL, dimids,tnumaerosol_vid) )
      call netcdf_check( nf90_put_att(ncid,tnumaerosol_vid,"title","Total number of aerosols defined by r<rc") )

      call netcdf_check( nf90_def_var(ncid, "tnum_destroy", NF90_REAL, dimids,tnum_destroy_vid) )
      call netcdf_check( nf90_put_att(ncid,tnum_destroy_vid,"title","Particles killed since last write") )

      call netcdf_check( nf90_def_var(ncid, "tdenum", NF90_REAL, dimids,tdenum_vid) )
      call netcdf_check( nf90_put_att(ncid,tdenum_vid,"title","Number of particles deactivated") )

      call netcdf_check( nf90_def_var(ncid, "tactnum", NF90_REAL, dimids,tactnum_vid) )
      call netcdf_check( nf90_put_att(ncid,tactnum_vid,"title","Number of particles activated") )

      call netcdf_check( nf90_def_var(ncid, "tnum100", NF90_REAL, dimids,tnum100_vid) )
      call netcdf_check( nf90_put_att(ncid,tnum100_vid,"title","Number of particles which did not converge Gauss-Newton solver since last history file write") )

      call netcdf_check( nf90_def_var(ncid, "tnumimpos", NF90_REAL, dimids,tnumimpos_vid) )
      call netcdf_check( nf90_put_att(ncid,tnumimpos_vid,"title","Number of particles which did not converge BOTH implicit solvers since last history file write") )

      call netcdf_check( nf90_def_var(ncid, "tot_reintro", NF90_REAL, dimids,tot_reintro_vid) )
      call netcdf_check( nf90_put_att(ncid,tot_reintro_vid,"title","Particles introduced this time step") )

      call netcdf_check( nf90_def_var(ncid, "Tsfc", NF90_REAL, dimids,Tsfc_vid) )
      call netcdf_check( nf90_put_att(ncid,Tsfc_vid,"title","Surface temperature") )

      call netcdf_check( nf90_def_var(ncid, "qsfc", NF90_REAL, dimids,qsfc_vid) )
      call netcdf_check( nf90_put_att(ncid,qsfc_vid,"title","Surface specific humidity qv") )

      call netcdf_check( nf90_def_var(ncid, "wtsfc", NF90_REAL, dimids,wtsfc_vid) )
      call netcdf_check( nf90_put_att(ncid,wtsfc_vid,"title","Surface sensible heat flux <w'T'>") )

      call netcdf_check( nf90_def_var(ncid, "wqsfc", NF90_REAL, dimids,wqsfc_vid) )
      call netcdf_check( nf90_put_att(ncid,wqsfc_vid,"title","Surface vapor flux <w'q'>") )

      call netcdf_check( nf90_def_var(ncid, "Swall", NF90_REAL, dimids,Swall_vid) )
      call netcdf_check( nf90_put_att(ncid,Swall_vid,"title","Humidity sink for Pi Chamber") )

      call netcdf_check( nf90_def_var(ncid, "meanRH", NF90_REAL, dimids,meanRH_vid) )
      call netcdf_check( nf90_put_att(ncid,meanRH_vid,"title","Volume average RH") )

      call netcdf_check( nf90_def_var(ncid, "varRH", NF90_REAL, dimids,varRH_vid) )
      call netcdf_check( nf90_put_att(ncid,varRH_vid,"title","Volume-averaged RH variance") )

      call netcdf_check( nf90_def_var(ncid, "radavg", NF90_REAL, dimids,radavg_vid) )
      call netcdf_check( nf90_put_att(ncid,radavg_vid,"title","Mean radius of activated droplets") )

      call netcdf_check( nf90_def_var(ncid, "radmsqr", NF90_REAL, dimids,radmsqr_vid) )
      call netcdf_check( nf90_put_att(ncid,radmsqr_vid,"title","Mean-squared radius of activated droplets") )

      call netcdf_check( nf90_def_var(ncid, "radavg_center", NF90_REAL, dimids,radavg_center_vid) )
      call netcdf_check( nf90_put_att(ncid,radavg_center_vid,"title","Mean radius of activated droplets only in domain interior") )

      call netcdf_check( nf90_def_var(ncid, "radmsqr_center", NF90_REAL, dimids,radmsqr_center_vid) )
      call netcdf_check( nf90_put_att(ncid,radmsqr_center_vid,"title","Mean-squared radius of activated droplets only in domain interior") )


!!! Profiles
      call netcdf_check( nf90_def_var(ncid, "zu", NF90_REAL, dimids_zu_only,zu_vid) )
      call netcdf_check( nf90_put_att(ncid,zu_vid,"title","z levels at u-points") )

      call netcdf_check( nf90_def_var(ncid, "zw", NF90_REAL, dimids_zw_only,zw_vid) )
      call netcdf_check( nf90_put_att(ncid,zw_vid,"title","z levels at w-points") )

      call netcdf_check( nf90_def_var(ncid,"uxym",NF90_REAL, dimids_zu,uxym_vid) )
      call netcdf_check( nf90_put_att(ncid,uxym_vid,"title","Horiz. avg. u vel") )

      call netcdf_check( nf90_def_var(ncid,"vxym",NF90_REAL, dimids_zu,vxym_vid) )
      call netcdf_check( nf90_put_att(ncid,vxym_vid,"title","Horiz. avg. v vel") )

      call netcdf_check( nf90_def_var(ncid,"wxym",NF90_REAL, dimids_zw,wxym_vid) )
      call netcdf_check( nf90_put_att(ncid,wxym_vid,"title","Horiz. avg. w vel") )

      call netcdf_check( nf90_def_var(ncid,"txym",NF90_REAL, dimids_zu_s,txym_vid) )
      call netcdf_check( nf90_put_att(ncid,txym_vid,"title","Horiz. avg. scalars, potential temp, specific humidity,...") )

      call netcdf_check( nf90_def_var(ncid,"exym",NF90_REAL, dimids_zw,exym_vid) )
      call netcdf_check( nf90_put_att(ncid,exym_vid,"title","Horiz. avg. subgrid energy e") )

      call netcdf_check( nf90_def_var(ncid,"RHxym",NF90_REAL, dimids_zu,RHxym_vid) )
      call netcdf_check( nf90_put_att(ncid,RHxym_vid,"title","Horiz. avg. relative humidity") )

      call netcdf_check( nf90_def_var(ncid,"RHmsqr",NF90_REAL, dimids_zu,RHmsqr_vid) )
      call netcdf_check( nf90_put_att(ncid,RHmsqr_vid,"title","Horiz. <RH^2>") )

      call netcdf_check( nf90_def_var(ncid,"tempxym",NF90_REAL, dimids_zu,tempxym_vid) )
      call netcdf_check( nf90_put_att(ncid,tempxym_vid,"title","Horiz. avg. temperature") )

      call netcdf_check( nf90_def_var(ncid,"ups",NF90_REAL, dimids_zu,ups_vid) )
      call netcdf_check( nf90_put_att(ncid,ups_vid,"title","Fluctuating velocity <u'^2>") )

      call netcdf_check( nf90_def_var(ncid,"vps",NF90_REAL, dimids_zu,vps_vid) )
      call netcdf_check( nf90_put_att(ncid,vps_vid,"title","Fluctuating velocity <v'^2>") )

      call netcdf_check( nf90_def_var(ncid,"wps",NF90_REAL, dimids_zw,wps_vid) )
      call netcdf_check( nf90_put_att(ncid,wps_vid,"title","Fluctuating velocity <w'^2>") )

      call netcdf_check( nf90_def_var(ncid,"tps",NF90_REAL, dimids_zu_s,tps_vid) )
      call netcdf_check( nf90_put_att(ncid,tps_vid,"title","Fluctuating scalars <t'^2>") )

      call netcdf_check( nf90_def_var(ncid,"uwle",NF90_REAL, dimids_zw,uwle_vid) )
      call netcdf_check( nf90_put_att(ncid,uwle_vid,"title","Resolved <u'w'>") )

      call netcdf_check( nf90_def_var(ncid,"uwsb",NF90_REAL, dimids_zw,uwsb_vid) )
      call netcdf_check( nf90_put_att(ncid,uwsb_vid,"title","Subgrid <u'w'>") )
     
      call netcdf_check( nf90_def_var(ncid,"vwle",NF90_REAL, dimids_zw,vwle_vid) )
      call netcdf_check( nf90_put_att(ncid,vwle_vid,"title","Resolved <v'w'>") )

      call netcdf_check( nf90_def_var(ncid,"vwsb",NF90_REAL, dimids_zw,vwsb_vid) )
      call netcdf_check( nf90_put_att(ncid,vwsb_vid,"title","Subgrid <v'w'>") )

      call netcdf_check( nf90_def_var(ncid,"wtle",NF90_REAL, dimids_zw_s,wtle_vid) )
      call netcdf_check( nf90_put_att(ncid,wtle_vid,"title","Resolved <w't'>") )

      call netcdf_check( nf90_def_var(ncid,"wtsb",NF90_REAL, dimids_zw_s,wtsb_vid) )
      call netcdf_check( nf90_put_att(ncid,wtsb_vid,"title","Subgrid <w't'>") )

      call netcdf_check( nf90_def_var(ncid,"t_rprod",NF90_REAL, dimids_zw,t_rprod_vid) )
      call netcdf_check( nf90_put_att(ncid,t_rprod_vid,"title","TKE resolved production") )

      call netcdf_check( nf90_def_var(ncid,"t_sprod",NF90_REAL, dimids_zw,t_sprod_vid) )
      call netcdf_check( nf90_put_att(ncid,t_sprod_vid,"title","TKE subgrid production") )

      call netcdf_check( nf90_def_var(ncid,"t_tran",NF90_REAL, dimids_zw,t_tran_vid) )
      call netcdf_check( nf90_put_att(ncid,t_tran_vid,"title","TKE transport") )

      call netcdf_check( nf90_def_var(ncid,"t_buoy",NF90_REAL, dimids_zw,t_buoy_vid) )
      call netcdf_check( nf90_put_att(ncid,t_buoy_vid,"title","TKE buoyancy production") )

      call netcdf_check( nf90_def_var(ncid,"t_diss",NF90_REAL, dimids_zw,t_diss_vid) )
      call netcdf_check( nf90_put_att(ncid,t_diss_vid,"title","TKE dissipation") )

      call netcdf_check( nf90_def_var(ncid,"zconc",NF90_REAL, dimids_zu,zconc_vid) )
      call netcdf_check( nf90_put_att(ncid,zconc_vid,"title","Computational droplet number concentration") )

      call netcdf_check( nf90_def_var(ncid,"pflux",NF90_REAL, dimids_zw,pflux_vid) )
      call netcdf_check( nf90_put_att(ncid,pflux_vid,"title","Flux of particles through zw points due to gravity/advection") )

      call netcdf_check( nf90_def_var(ncid,"pfluxdiff",NF90_REAL, dimids_zw,pfluxdiff_vid) )
      call netcdf_check( nf90_put_att(ncid,pfluxdiff_vid,"title","Flux of particles through zw points due to SFS diffusion") )

      call netcdf_check( nf90_def_var(ncid,"pmassflux",NF90_REAL, dimids_zw,pmassflux_vid) )
      call netcdf_check( nf90_put_att(ncid,pmassflux_vid,"title","Flux of droplet mass gained/lost across zw points") )

      call netcdf_check( nf90_def_var(ncid,"penegflux",NF90_REAL, dimids_zw,penegflux_vid) )
      call netcdf_check( nf90_put_att(ncid,penegflux_vid,"title","Flux of droplet energy gained/lost across zw points") )

      call netcdf_check( nf90_def_var(ncid,"vp1mean",NF90_REAL, dimids_zu,vp1mean_vid) )
      call netcdf_check( nf90_put_att(ncid,vp1mean_vid,"title","Horiz. avg. particle velocity u") )

      call netcdf_check( nf90_def_var(ncid,"vp2mean",NF90_REAL, dimids_zu,vp2mean_vid) )
      call netcdf_check( nf90_put_att(ncid,vp2mean_vid,"title","Horiz. avg. particle velocity v") )

      call netcdf_check( nf90_def_var(ncid,"vp3mean",NF90_REAL, dimids_zu,vp3mean_vid) )
      call netcdf_check( nf90_put_att(ncid,vp3mean_vid,"title","Horiz. avg. particle velocity w") )

      call netcdf_check( nf90_def_var(ncid,"vp1msqr",NF90_REAL, dimids_zu,vp1msqr_vid) )
      call netcdf_check( nf90_put_att(ncid,vp1msqr_vid,"title","Mean squared particle velocity <u'^2>") )

      call netcdf_check( nf90_def_var(ncid,"vp2msqr",NF90_REAL, dimids_zu,vp2msqr_vid) )
      call netcdf_check( nf90_put_att(ncid,vp2msqr_vid,"title","Mean squared particle velocity <v'^2>") )

      call netcdf_check( nf90_def_var(ncid,"vp3msqr",NF90_REAL, dimids_zu,vp3msqr_vid) )
      call netcdf_check( nf90_put_att(ncid,vp3msqr_vid,"title","Mean squared particle velocity <w'^2>") )

      call netcdf_check( nf90_def_var(ncid,"uf1mean",NF90_REAL, dimids_zu,uf1mean_vid) )
      call netcdf_check( nf90_put_att(ncid,uf1mean_vid,"title","Horiz. avg. velocity u seen by particle") )

      call netcdf_check( nf90_def_var(ncid,"uf2mean",NF90_REAL, dimids_zu,uf2mean_vid) )
      call netcdf_check( nf90_put_att(ncid,uf2mean_vid,"title","Horiz. avg. velocity v seen by particle") )

      call netcdf_check( nf90_def_var(ncid,"uf3mean",NF90_REAL, dimids_zu,uf3mean_vid) )
      call netcdf_check( nf90_put_att(ncid,uf3mean_vid,"title","Horiz. avg. velocity w seen by particle") )

      call netcdf_check( nf90_def_var(ncid,"uf1msqr",NF90_REAL, dimids_zu,uf1msqr_vid) )
      call netcdf_check( nf90_put_att(ncid,uf1msqr_vid,"title","Mean squared u-velocity seen by particle <uf'^2>") )

      call netcdf_check( nf90_def_var(ncid,"uf2msqr",NF90_REAL, dimids_zu,uf2msqr_vid) )
      call netcdf_check( nf90_put_att(ncid,uf2msqr_vid,"title","Mean squared v-velocity seen by particle <vf'^2>") )

      call netcdf_check( nf90_def_var(ncid,"uf3msqr",NF90_REAL, dimids_zu,uf3msqr_vid) )
      call netcdf_check( nf90_put_att(ncid,uf3msqr_vid,"title","Mean squared w-velocity seen by particle <wf'^2>") )

      call netcdf_check( nf90_def_var(ncid,"m1src",NF90_REAL, dimids_zu,m1src_vid) )
      call netcdf_check( nf90_put_att(ncid,m1src_vid,"title","Particle x-momentum source") )

      call netcdf_check( nf90_def_var(ncid,"m2src",NF90_REAL, dimids_zu,m2src_vid) )
      call netcdf_check( nf90_put_att(ncid,m2src_vid,"title","Particle y-momentum source") )

      call netcdf_check( nf90_def_var(ncid,"m3src",NF90_REAL, dimids_zu,m3src_vid) )
      call netcdf_check( nf90_put_att(ncid,m3src_vid,"title","Particle z-momentum source") )

      call netcdf_check( nf90_def_var(ncid,"Tpsrc",NF90_REAL, dimids_zu,Tpsrc_vid) )
      call netcdf_check( nf90_put_att(ncid,Tpsrc_vid,"title","Particle temp. source (sensible)") )

      call netcdf_check( nf90_def_var(ncid,"TEpsrc",NF90_REAL, dimids_zu,TEpsrc_vid) )
      call netcdf_check( nf90_put_att(ncid,TEpsrc_vid,"title","Particle temp. source (latent)") )

      call netcdf_check( nf90_def_var(ncid,"Hpsrc",NF90_REAL, dimids_zu,Hpsrc_vid) )
      call netcdf_check( nf90_put_att(ncid,Hpsrc_vid,"title","Particle vapor source") )

      call netcdf_check( nf90_def_var(ncid,"Tpmean",NF90_REAL, dimids_zu,Tpmean_vid) )
      call netcdf_check( nf90_put_att(ncid,Tpmean_vid,"title","Horiz. avg. particle temp") )

      call netcdf_check( nf90_def_var(ncid,"Tpmsqr",NF90_REAL, dimids_zu,Tpmsqr_vid) )
      call netcdf_check( nf90_put_att(ncid,Tpmsqr_vid,"title","Mean squared particle temp <Tp^2>") )

      call netcdf_check( nf90_def_var(ncid,"Tfmean",NF90_REAL, dimids_zu,Tfmean_vid) )
      call netcdf_check( nf90_put_att(ncid,Tfmean_vid,"title","Horiz. avg. fluid temp at particle") )

      call netcdf_check( nf90_def_var(ncid,"qfmean",NF90_REAL, dimids_zu,qfmean_vid) )
      call netcdf_check( nf90_put_att(ncid,qfmean_vid,"title","Horiz. avg. fluid qv at particle") )

      call netcdf_check( nf90_def_var(ncid,"radmean",NF90_REAL, dimids_zu,radmean_vid) )
      call netcdf_check( nf90_put_att(ncid,radmean_vid,"title","Horiz. avg. particle radius") )

      call netcdf_check( nf90_def_var(ncid,"rad2mean",NF90_REAL, dimids_zu,rad2mean_vid) )
      call netcdf_check( nf90_put_att(ncid,rad2mean_vid,"title","Mean squared particle radius <rp^2>") )

      call netcdf_check( nf90_def_var(ncid,"qstarm",NF90_REAL, dimids_zu,qstarm_vid) )
      call netcdf_check( nf90_put_att(ncid,qstarm_vid,"title","Horiz. avg. qstar") )

      call netcdf_check( nf90_def_var(ncid,"Nc",NF90_REAL, dimids_zu,Nc_vid) )
      call netcdf_check( nf90_put_att(ncid,Nc_vid,"title","Horiz. avg. total number concentration") )

      call netcdf_check( nf90_def_var(ncid,"ql",NF90_REAL, dimids_zu,ql_vid) )
      call netcdf_check( nf90_put_att(ncid,ql_vid,"title","Horiz. avg. liquid droplet mixing ratio") )

      call netcdf_check( nf90_def_var(ncid,"radsrc",NF90_REAL, dimids_zu,radsrc_vid) )
      call netcdf_check( nf90_put_att(ncid,radsrc_vid,"title","Radiation tendency from the longwave model") )

      call netcdf_check( nf90_def_var(ncid,"rho_base",NF90_REAL, dimids_zu,rho_base_vid) )
      call netcdf_check( nf90_put_att(ncid,rho_base_vid,"title","Base state density used in Boussinesq") )

      call netcdf_check( nf90_def_var(ncid,"p_base",NF90_REAL, dimids_zu,p_base_vid) )
      call netcdf_check( nf90_put_att(ncid,p_base_vid,"title","Base state pressure used in Boussinesq") )

      call netcdf_check( nf90_def_var(ncid,"T_base",NF90_REAL, dimids_zu,T_base_vid) )
      call netcdf_check( nf90_put_att(ncid,T_base_vid,"title","Base state temperature used in Boussinesq") )

      call netcdf_check( nf90_def_var(ncid,"theta_base",NF90_REAL, dimids_zu,theta_base_vid) )
      call netcdf_check( nf90_put_att(ncid,theta_base_vid,"title","Base state potential temperature used in Boussinesq") )

      call netcdf_check( nf90_enddef(ncid) )

      his_counter = 1

end subroutine netcdf_init

subroutine netcdf_res
!https://www.unidata.ucar.edu/software/netcdf/docs-fortran/f90-use-of-the-netcdf-library.html#f90-writing-data-in-an-existing-netcdf-dataset
      use netcdf
      use con_data
      use pars
      implicit none

      integer :: dimids(1),dimids_zu(2),dimids_zw(2),dimids_zu_s(3),dimids_zw_s(3)
      integer :: dimids_zu_only(1),dimids_zw_only(1)
      integer :: Ntime

      path_netcdf_his = trim(adjustl(path_seed))//"history.nc"

      call netcdf_check( nf90_open(path_netcdf_his,NF90_WRITE,ncid) )

      call netcdf_check (nf90_inq_dimid(ncid,'time',time_dimid) )
      call netcdf_check (nf90_inq_dimid(ncid,'zu',zu_dimid) )
      call netcdf_check (nf90_inq_dimid(ncid,'zw',zw_dimid) )
      call netcdf_check (nf90_inq_dimid(ncid,'nscl',s_dimid) )

      dimids = (/ time_dimid /)
      dimids_zu = (/ zu_dimid, time_dimid/)
      dimids_zw = (/ zw_dimid, time_dimid/)
      dimids_zu_only = (/ zu_dimid/)
      dimids_zw_only = (/ zw_dimid/)
      dimids_zu_s = (/ zu_dimid, s_dimid, time_dimid/)
      dimids_zw_s = (/ zw_dimid, s_dimid, time_dimid/)

      if (myid==0) write(*,*) 'Restart his_counter for netCDF = ',his_counter

!!! Single quantities
      call netcdf_check( nf90_inq_varid(ncid,"time",time_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"dt",dt_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"utau",utau_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"uwsfc",uwsfc_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"tnumpart",tnumpart_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"tnumdrop",tnumdrop_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"tnumaerosol",tnumaerosol_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"tnum_destroy",tnum_destroy_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"tdenum",tdenum_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"tactnum",tactnum_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"tnum100",tnum100_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"tnumimpos",tnumimpos_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"tot_reintro",tot_reintro_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"Tsfc",Tsfc_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"qsfc",qsfc_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"wtsfc",wtsfc_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"wqsfc",wqsfc_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"Swall",Swall_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"meanRH",meanRH_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"varRH",varRH_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"radavg",radavg_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"radmsqr",radmsqr_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"radavg_center",radavg_center_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"radmsqr_center",radmsqr_center_vid) )


!!! Profiles
      call netcdf_check( nf90_inq_varid(ncid,"zu",zu_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"zw",zw_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"uxym",uxym_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vxym",vxym_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"wxym",wxym_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"txym",txym_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"exym",exym_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"RHxym",RHxym_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"RHmsqr",RHmsqr_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"tempxym",tempxym_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"ups",ups_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vps",vps_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"wps",wps_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"tps",tps_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"uwle",uwle_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"uwsb",uwsb_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vwle",vwle_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vwsb",vwsb_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"wtle",wtle_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"wtsb",wtsb_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"t_rprod",t_rprod_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"t_sprod",t_sprod_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"t_tran",t_tran_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"t_buoy",t_buoy_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"t_diss",t_diss_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"zconc",zconc_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"pflux",pflux_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"pfluxdiff",pfluxdiff_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"pmassflux",pmassflux_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"penegflux",penegflux_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vp1mean",vp1mean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vp2mean",vp2mean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vp3mean",vp3mean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vp1msqr",vp1msqr_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vp2msqr",vp2msqr_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vp3msqr",vp3msqr_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"uf1mean",uf1mean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"uf2mean",uf2mean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"uf3mean",uf3mean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"uf1msqr",uf1msqr_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"uf2msqr",uf2msqr_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"uf3msqr",uf3msqr_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"m1src",m1src_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"m2src",m2src_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"m3src",m3src_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"Tpsrc",Tpsrc_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"TEpsrc",TEpsrc_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"Hpsrc",Hpsrc_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"Tpmean",Tpmean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"Tpmsqr",Tpmsqr_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"Tfmean",Tfmean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"qfmean",qfmean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"radmean",radmean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"rad2mean",rad2mean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"qstarm",qstarm_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"Nc",Nc_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"ql",ql_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"radsrc",radsrc_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"rho_base",rho_base_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"p_base",p_base_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"T_base",T_base_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"theta_base",theta_base_vid) )


end subroutine netcdf_res

subroutine write_his_netcdf
      use netcdf
      use pars
      use fields
      use con_data
      use con_stats
      use particles

      implicit none
      real :: tmp(0:nnz),tmp_s(0:nnz,nscl)

      write(*,'(a35,i)') 'WRITE HISTORY FILE, HIS_COUNTER = ',his_counter

      call netcdf_check( nf90_put_var(ncid, time_vid, real(time),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, dt_vid, real(dt),start=(/his_counter/)) )

      !Only write the grid on the first call
      if (his_counter.eq.1) then
      call netcdf_check( nf90_put_var(ncid, zu_vid, real(zz(1:nnz)),start=(/1/)) )
      call netcdf_check( nf90_put_var(ncid, zw_vid, real(z(0:nnz)),start=(/1/)) )
      end if

      call netcdf_check( nf90_put_var(ncid, utau_vid, real(utau),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, uwsfc_vid, real(uwsfc),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, tnumpart_vid, real(tnumpart),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, tnumdrop_vid, real(tnumdrop),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, tnumaerosol_vid, real(tnumaerosol),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, tnum_destroy_vid, real(tnum_destroy_accum),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, tdenum_vid, real(tdenum),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, tactnum_vid, real(tactnum),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, tnum100_vid, real(tnum100),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, tnumimpos_vid, real(tnumimpos),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, tot_reintro_vid, real(tot_reintro),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, Tsfc_vid, real(tsfcc(1)),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, qsfc_vid, real(tsfcc(2)),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, wtsfc_vid, real(wtsfc(1)),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, wqsfc_vid, real(wtsfc(2)),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, Swall_vid, real(Swall),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, meanRH_vid, real(meanRH),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, varRH_vid, real(varRH),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, radavg_vid, real(radavg),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, radmsqr_vid, real(radmsqr),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, radavg_center_vid, real(radavg_center),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, radmsqr_center_vid, real(radmsqr_center),start=(/his_counter/)) )

      call netcdf_check( nf90_put_var(ncid,uxym_vid,real(uxym(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vxym_vid,real(vxym(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,wxym_vid,real(wxym(0:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,txym_vid,real(txym(1:nnz,1:nscl)),start=(/1,1,his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,exym_vid,real(e_mn(0:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,RHxym_vid,real(RHxym(1:nnz)),start=(/1,his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,RHmsqr_vid,real(RHmsqr(1:nnz)),start=(/1,his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,tempxym_vid,real(tempxym(1:nnz)),start=(/1,his_counter/)) )

      call netcdf_check( nf90_put_var(ncid,ups_vid,real(ups(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vps_vid,real(vps(1:nnz)),start=(/1, his_counter/)) )

      tmp(0) = 0.0
      tmp(1:nnz) = wps(1:nnz)
      call netcdf_check( nf90_put_var(ncid,wps_vid,real(tmp),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,tps_vid,real(tps(1:nnz,1:nscl)),start=(/1,1,his_counter/)) )

      tmp(0) = 0
      tmp(1:nnz) = uwle(1:nnz)
      call netcdf_check( nf90_put_var(ncid,uwle_vid,real(tmp),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,uwsb_vid,real(uwsb(0:nnz)),start=(/1, his_counter/)) )

      tmp(0) = 0
      tmp(1:nnz) = vwle(1:nnz)
      call netcdf_check( nf90_put_var(ncid,vwle_vid,real(tmp),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vwsb_vid,real(vwsb(0:nnz)),start=(/1, his_counter/)) )

      tmp_s(0,1:nscl) = 0
      tmp_s(1:nnz,1:nscl) = wtle(1:nnz,1:nscl)
      call netcdf_check( nf90_put_var(ncid,wtle_vid,real(tmp_s),start=(/1,1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,wtsb_vid,real(wtsb(0:nnz,1:nscl)),start=(/1,1, his_counter/)) )


      tmp(0) = 0  !Need to change this
      tmp(1:nnz) = t_rprod(1:nnz)
      call netcdf_check( nf90_put_var(ncid,t_rprod_vid,real(tmp),start=(/1, his_counter/)) )

      tmp(0) = 0  !Need to change this
      tmp(1:nnz) = t_sprod(1:nnz)
      call netcdf_check( nf90_put_var(ncid,t_sprod_vid,real(tmp),start=(/1, his_counter/)) )

      tmp(0) = 0  !Need to change this
      tmp(1:nnz) = t_tran(1:nnz)
      call netcdf_check( nf90_put_var(ncid,t_tran_vid,real(tmp),start=(/1, his_counter/)) )

      tmp(0) = 0  !Need to change this
      tmp(1:nnz) = t_buoy(1:nnz)
      call netcdf_check( nf90_put_var(ncid,t_buoy_vid,real(tmp),start=(/1, his_counter/)) )

      tmp(0) = 0  !Need to change this
      tmp(1:nnz) = t_diss(1:nnz)
      call netcdf_check( nf90_put_var(ncid,t_diss_vid,real(tmp),start=(/1, his_counter/)) )

      call netcdf_check( nf90_put_var(ncid,zconc_vid,real(zconc(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,pflux_vid,real(pflux(0:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,pfluxdiff_vid,real(pfluxdiff(0:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,pmassflux_vid,real(pmassflux(0:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,penegflux_vid,real(penegflux(0:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vp1mean_vid,real(vp1mean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vp2mean_vid,real(vp2mean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vp3mean_vid,real(vp3mean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vp1msqr_vid,real(vp1msqr(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vp2msqr_vid,real(vp2msqr(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vp3msqr_vid,real(vp3msqr(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,uf1mean_vid,real(uf1mean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,uf2mean_vid,real(uf2mean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,uf3mean_vid,real(uf3mean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,uf1msqr_vid,real(uf1msqr(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,uf2msqr_vid,real(uf2msqr(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,uf3msqr_vid,real(uf3msqr(1:nnz)),start=(/1, his_counter/)) )

      call netcdf_check( nf90_put_var(ncid,m1src_vid,real(m1src(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,m2src_vid,real(m2src(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,m3src_vid,real(m3src(1:nnz)),start=(/1, his_counter/)) )

      call netcdf_check( nf90_put_var(ncid,Tpsrc_vid,real(Tpsrc(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,TEpsrc_vid,real(TEpsrc(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,Hpsrc_vid,real(Hpsrc(1:nnz)),start=(/1, his_counter/)) )

      call netcdf_check( nf90_put_var(ncid,Tpmean_vid,real(Tpmean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,Tpmsqr_vid,real(Tpmsqr(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,Tfmean_vid,real(Tfmean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,qfmean_vid,real(qfmean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,radmean_vid,real(radmean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,rad2mean_vid,real(rad2mean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,qstarm_vid,real(qstarm(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,Nc_vid,real(Nc(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,ql_vid,real(ql(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,radsrc_vid,real(radsrc(1:nnz)),start=(/1, his_counter/)) )

      call netcdf_check( nf90_put_var(ncid,rho_base_vid,real(rho_base(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,p_base_vid,real(p_base(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,T_base_vid,real(T_base(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,theta_base_vid,real(theta_base(1:nnz)),start=(/1, his_counter/)) )

      !Reset some cumulative quantities

      num100 = 0
      numimpos = 0

      tnum_destroy_accum = 0
      his_counter = his_counter + 1

end subroutine write_his_netcdf

subroutine netcdf_init_histog
      use netcdf
      use pars
      use con_data
      use particles
      implicit none

      integer :: dimids(1),dimids_bins(1),dimids_t_bins(2)

      path_netcdf_histog = trim(adjustl(path_seed))//"histograms.nc"

      call netcdf_check( nf90_create(path_netcdf_histog,nf90_clobber,ncid_histog))

      call netcdf_check( nf90_def_dim(ncid_histog, "time",NF90_UNLIMITED, time_histog_dimid) )

      call netcdf_check( nf90_def_dim(ncid_histog,"histbins",histbins+2,histbins_dimid) )

      dimids = (/ time_histog_dimid /)
      dimids_bins = (/ histbins_dimid /)
      dimids_t_bins = (/histbins_dimid, time_histog_dimid /)


!!! Single quantities
      call netcdf_check( nf90_def_var(ncid_histog,"time",NF90_REAL,dimids,time_histog_vid) )
      call netcdf_check( nf90_put_att(ncid_histog,time_histog_vid,"title","Simulation time") )

!! Store the bin definitions
      call netcdf_check( nf90_def_var(ncid_histog,"radbins",NF90_DOUBLE,dimids_bins,radbins_vid) )
      call netcdf_check( nf90_put_att(ncid_histog,radbins_vid,"title","Bin centers of particle radii histogram") )

      call netcdf_check( nf90_def_var(ncid_histog,"resbins",NF90_DOUBLE,dimids_bins,resbins_vid) )
      call netcdf_check( nf90_put_att(ncid_histog,resbins_vid,"title","Bin centers of particle residence time histogram") )

      call netcdf_check( nf90_def_var(ncid_histog,"actresbins",NF90_DOUBLE,dimids_bins,actresbins_vid) )
      call netcdf_check( nf90_put_att(ncid_histog,actresbins_vid,"title","Bin centers of particle activated residence time histogram") )

      call netcdf_check( nf90_def_var(ncid_histog,"acttodeathbins",NF90_DOUBLE,dimids_bins,acttodeathbins_vid) )
      call netcdf_check( nf90_put_att(ncid_histog,acttodeathbins_vid,"title","Bin centers of particle activated residence time till death histogram") )

      call netcdf_check( nf90_def_var(ncid_histog,"numactbins",NF90_DOUBLE,dimids_bins,numactbins_vid) )
      call netcdf_check( nf90_put_att(ncid_histog,numactbins_vid,"title","Bin centers of histogram of number of activations in lifetime") )


!!! Histograms
      call netcdf_check( nf90_def_var(ncid_histog, "radhist", NF90_REAL, dimids_t_bins,radhist_vid) )
      call netcdf_check( nf90_put_att(ncid_histog,radhist_vid,"title","Histogram of particle radius at snapshot") )

      call netcdf_check( nf90_def_var(ncid_histog, "radhistdeath", NF90_REAL, dimids_t_bins,radhistdeath_vid) )
      call netcdf_check( nf90_put_att(ncid_histog,radhistdeath_vid,"title","Histogram of particle radius of those which died over past ihst") )

      call netcdf_check( nf90_def_var(ncid_histog, "reshist", NF90_REAL, dimids_t_bins,reshist_vid) )
      call netcdf_check( nf90_put_att(ncid_histog,reshist_vid,"title","Histogram of residence time of particles which died over past ihst") )

      call netcdf_check( nf90_def_var(ncid_histog, "actreshist", NF90_REAL, dimids_t_bins,actreshist_vid) )
      call netcdf_check( nf90_put_att(ncid_histog,actreshist_vid,"title","Histogram of residence time of particles which deactivated OR died over past ihst") )

      call netcdf_check( nf90_def_var(ncid_histog, "acttodeathhist", NF90_REAL, dimids_t_bins,acttodeathhist_vid) )
      call netcdf_check( nf90_put_att(ncid_histog,acttodeathhist_vid,"title","Histogram of residence time of particles which deactivated OR died over past ihst") )

      call netcdf_check( nf90_def_var(ncid_histog, "numacthist", NF90_REAL, dimids_t_bins,numacthist_vid) )
      call netcdf_check( nf90_put_att(ncid_histog,numacthist_vid,"title","Histogram of number of activations of particles which died over past ihst") )

      call netcdf_check( nf90_enddef(ncid_histog) )

      histog_counter = 1

end subroutine netcdf_init_histog

subroutine netcdf_res_histog
      use netcdf
      use pars
      use con_data
      use particles
      implicit none

      integer :: dimids(1),dimids_bins(1),dimids_t_bins(2)
      integer :: Ntime


      path_netcdf_histog = trim(adjustl(path_seed))//"histograms.nc"

      call netcdf_check( nf90_open(path_netcdf_histog,NF90_WRITE,ncid_histog) )

      call netcdf_check (nf90_inq_dimid(ncid_histog,'time',time_histog_dimid) )
      call netcdf_check (nf90_inq_dimid(ncid_histog,'histbins',histbins_dimid) )

      if (myid==0) write(*,*) 'Restart histog_counter for netCDF = ',histog_counter

      dimids = (/ time_histog_dimid /)
      dimids_bins = (/ histbins_dimid /)
      dimids_t_bins = (/histbins_dimid, time_histog_dimid /)


!!! Single quantities
      call netcdf_check( nf90_inq_varid(ncid_histog,"time",time_histog_vid) )

!! Store the bin definitions
      call netcdf_check( nf90_inq_varid(ncid_histog,"radbins",radbins_vid) )
      call netcdf_check( nf90_inq_varid(ncid_histog,"resbins",resbins_vid) )
      call netcdf_check( nf90_inq_varid(ncid_histog,"actresbins",actresbins_vid) )
      call netcdf_check( nf90_inq_varid(ncid_histog,"acttodeathbins",acttodeathbins_vid) )
      call netcdf_check( nf90_inq_varid(ncid_histog,"numactbins",numactbins_vid) )

!!! Histograms
      call netcdf_check( nf90_inq_varid(ncid_histog, "radhist",radhist_vid) )
      call netcdf_check( nf90_inq_varid(ncid_histog, "radhistdeath",radhistdeath_vid) )
      call netcdf_check( nf90_inq_varid(ncid_histog, "reshist",reshist_vid) )
      call netcdf_check( nf90_inq_varid(ncid_histog, "actreshist",actreshist_vid) )
      call netcdf_check( nf90_inq_varid(ncid_histog, "acttodeathhist",acttodeathhist_vid) )
      call netcdf_check( nf90_inq_varid(ncid_histog, "numacthist",numacthist_vid) )

end subroutine netcdf_res_histog

subroutine write_histog_netcdf
      use netcdf
      use pars
      use fields
      use con_data
      use con_stats
      use particles

      implicit none

      call netcdf_check( nf90_put_var(ncid_histog, time_histog_vid, real(time),start=(/histog_counter/)) )

      !Store bin definitions only once
      if (histog_counter == 1) then

      call netcdf_check( nf90_put_var(ncid_histog, radbins_vid, bins_rad,start=(/ 1 /)) )
      call netcdf_check( nf90_put_var(ncid_histog, resbins_vid, bins_res,start=(/ 1 /)) )
      call netcdf_check( nf90_put_var(ncid_histog, actresbins_vid, bins_actres,start=(/ 1 /)) )
      call netcdf_check( nf90_put_var(ncid_histog, acttodeathbins_vid, bins_acttodeath,start=(/ 1 /)) )
      call netcdf_check( nf90_put_var(ncid_histog, numactbins_vid, bins_numact,start=(/ 1 /)) )

      end if

      call netcdf_check( nf90_put_var(ncid_histog, radhist_vid, real(hist_rad),start=(/1,histog_counter/)) )
      call netcdf_check( nf90_put_var(ncid_histog, radhistdeath_vid, real(hist_raddeath),start=(/1,histog_counter/)) )
      call netcdf_check( nf90_put_var(ncid_histog, reshist_vid, real(hist_res),start=(/1,histog_counter/)) )
      call netcdf_check( nf90_put_var(ncid_histog, actreshist_vid, real(hist_actres),start=(/1,histog_counter/)) )
      call netcdf_check( nf90_put_var(ncid_histog, acttodeathhist_vid, real(hist_acttodeath),start=(/1,histog_counter/)) )
      call netcdf_check( nf90_put_var(ncid_histog, numacthist_vid, real(hist_numact),start=(/1,histog_counter/)) )

      histog_counter = histog_counter + 1

end subroutine write_histog_netcdf

subroutine netcdf_init_viz
      use netcdf
      use pars
      use con_data
      use particles
      implicit none
      include 'mpif.h'

      integer :: dimids(1),dimids_xy(3),dimids_yz(3),dimids_xz(3)
      integer :: dimids_xgrid(1),dimids_ygrid(1),dimids_zugrid(1),dimids_zwgrid(1)
      integer :: ierr

      path_netcdf_viz = trim(adjustl(path_seed))//"viz.nc"

      call netcdf_check( nf90_create(path_netcdf_viz,nf90_clobber,ncid_viz))

      call netcdf_check( nf90_def_dim(ncid_viz, "time",NF90_UNLIMITED, time_viz_dimid) )

      call netcdf_check( nf90_def_dim(ncid_viz,"nx",nnx,viz_nx_dimid) )
      call netcdf_check( nf90_def_dim(ncid_viz,"ny",nny,viz_ny_dimid) )
      call netcdf_check( nf90_def_dim(ncid_viz,"nzu",nnz,viz_nzu_dimid) )
      call netcdf_check( nf90_def_dim(ncid_viz,"nzw",nnz+1,viz_nzw_dimid) )

      dimids = (/ time_viz_dimid /)
      dimids_xy = (/viz_nx_dimid, viz_ny_dimid, time_viz_dimid /)
      dimids_yz = (/viz_ny_dimid, viz_nzu_dimid, time_viz_dimid /)
      dimids_xz = (/viz_nx_dimid, viz_nzu_dimid, time_viz_dimid /)

      dimids_xgrid = (/viz_nx_dimid /)
      dimids_ygrid = (/viz_ny_dimid /)
      dimids_zugrid = (/viz_nzu_dimid /)
      dimids_zwgrid = (/viz_nzw_dimid /)

!!! Single quantities
      call netcdf_check( nf90_def_var(ncid_viz,"time",NF90_REAL,dimids,time_viz_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,time_viz_vid,"title","Simulation time") )


!! Would need to store grid values

      call netcdf_check( nf90_def_var(ncid_viz, "x", NF90_REAL, dimids_xgrid,xgrid_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,xgrid_vid,"title","x locations") )

      call netcdf_check( nf90_def_var(ncid_viz, "y", NF90_REAL, dimids_ygrid,ygrid_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,ygrid_vid,"title","x locations") )

      call netcdf_check( nf90_def_var(ncid_viz, "zu", NF90_REAL, dimids_zugrid,zugrid_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,zugrid_vid,"title","z locations at u,v points") )

      call netcdf_check( nf90_def_var(ncid_viz, "zw", NF90_REAL, dimids_zwgrid,zwgrid_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,zwgrid_vid,"title","z locations at w points") )


!!! Slices

      !yz slices
      call netcdf_check( nf90_def_var(ncid_viz, "u_yz", NF90_REAL, dimids_yz,u_yz_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,u_yz_vid,"title","yz slice of u-velocity") )

      call netcdf_check( nf90_def_var(ncid_viz, "t_yz", NF90_REAL, dimids_yz,t_yz_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,t_yz_vid,"title","yz slice of temperature") )

      call netcdf_check( nf90_def_var(ncid_viz, "q_yz", NF90_REAL, dimids_yz,q_yz_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,q_yz_vid,"title","yz slice of water vapor mixing ratio") )

      call netcdf_check( nf90_def_var(ncid_viz, "w_yz", NF90_REAL, dimids_yz,w_yz_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,w_yz_vid,"title","yz slice of w-velocity") )

      call netcdf_check( nf90_def_var(ncid_viz, "v_yz", NF90_REAL, dimids_yz,v_yz_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,v_yz_vid,"title","yz slice of v-velocity") )

      !xy slices
      call netcdf_check( nf90_def_var(ncid_viz, "u_xy", NF90_REAL, dimids_xy,u_xy_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,u_xy_vid,"title","xy slice of u-velocity") )

      call netcdf_check( nf90_def_var(ncid_viz, "v_xy", NF90_REAL, dimids_xy,v_xy_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,v_xy_vid,"title","xy slice of v-velocity") )

      call netcdf_check( nf90_def_var(ncid_viz, "w_xy", NF90_REAL, dimids_xy,w_xy_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,w_xy_vid,"title","xy slice of w-velocity") )

      call netcdf_check( nf90_def_var(ncid_viz, "t_xy", NF90_REAL, dimids_xy,t_xy_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,t_xy_vid,"title","xy slice of temperature") )

      call netcdf_check( nf90_def_var(ncid_viz, "q_xy", NF90_REAL, dimids_xy,q_xy_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,q_xy_vid,"title","xy slice of water vapor mixing ratio") )

      call netcdf_check( nf90_def_var(ncid_viz, "LWP", NF90_REAL, dimids_xy,LWP_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,LWP_vid,"title","Liquid water path [kg/m^2]") )

      call netcdf_check( nf90_def_var(ncid_viz, "surf_precip", NF90_REAL, dimids_xy,surf_precip_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,surf_precip_vid,"title","Surface precipitation rate [kg/m^2/s]") )

      !xz slices
      call netcdf_check( nf90_def_var(ncid_viz, "u_xz", NF90_REAL, dimids_xz,u_xz_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,u_xz_vid,"title","xz slice of u-velocity") )

      call netcdf_check( nf90_def_var(ncid_viz, "v_xz", NF90_REAL, dimids_xz,v_xz_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,v_xz_vid,"title","xz slice of v-velocity") )

      call netcdf_check( nf90_def_var(ncid_viz, "w_xz", NF90_REAL, dimids_xz,w_xz_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,w_xz_vid,"title","xz slice of w-velocity") )

      call netcdf_check( nf90_def_var(ncid_viz, "t_xz", NF90_REAL, dimids_xz,t_xz_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,t_xz_vid,"title","xz slice of temperature") )

      call netcdf_check( nf90_def_var(ncid_viz, "q_xz", NF90_REAL, dimids_xz,q_xz_vid) )
      call netcdf_check( nf90_put_att(ncid_viz,q_xz_vid,"title","xz slice of water vapor mixing ratio") )

      call netcdf_check( nf90_enddef(ncid_viz) )

      viz_counter = 1

end subroutine netcdf_init_viz

subroutine netcdf_res_viz
      use netcdf
      use pars
      use con_data
      use particles
      implicit none
      include 'mpif.h'

      integer :: dimids(1),dimids_xy(3),dimids_yz(3),dimids_xz(3)
      integer :: dimids_xgrid(1),dimids_ygrid(1),dimids_zugrid(1),dimids_zwgrid(1)
      integer :: ierr, Ntime


      path_netcdf_viz = trim(adjustl(path_seed))//"viz.nc"

      call netcdf_check( nf90_open(path_netcdf_viz,NF90_WRITE,ncid_viz) )

      call netcdf_check (nf90_inq_dimid(ncid_viz,'time',time_viz_dimid) )
      call netcdf_check (nf90_inq_dimid(ncid_viz,'nx',viz_nx_dimid) )
      call netcdf_check (nf90_inq_dimid(ncid_viz,'ny',viz_ny_dimid) )
      call netcdf_check (nf90_inq_dimid(ncid_viz,'nzu',viz_nzu_dimid) )
      call netcdf_check (nf90_inq_dimid(ncid_viz,'nzw',viz_nzw_dimid) )


      if (myid==0) write(*,*) 'Restart viz_counter for netCDF = ',viz_counter


      dimids = (/ time_viz_dimid /)
      dimids_xy = (/viz_nx_dimid, viz_ny_dimid, time_viz_dimid /)
      dimids_yz = (/viz_ny_dimid, viz_nzu_dimid, time_viz_dimid /)
      dimids_xz = (/viz_nx_dimid, viz_nzu_dimid, time_viz_dimid /)

      dimids_xgrid = (/viz_nx_dimid /)
      dimids_ygrid = (/viz_ny_dimid /)
      dimids_zugrid = (/viz_nzu_dimid /)
      dimids_zwgrid = (/viz_nzw_dimid /)

!!! Single quantities
      call netcdf_check( nf90_inq_varid(ncid_viz,"time",time_viz_vid) )

!!store grid values

      call netcdf_check( nf90_inq_varid(ncid_viz, "x",xgrid_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "y",ygrid_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "zu",zugrid_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "zw",zwgrid_vid) )

!!! Slices

      !yz slices
      call netcdf_check( nf90_inq_varid(ncid_viz, "u_yz",u_yz_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "t_yz",t_yz_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "q_yz",q_yz_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "w_yz",w_yz_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "v_yz",v_yz_vid) )

      !xy slices
      call netcdf_check( nf90_inq_varid(ncid_viz, "u_xy",u_xy_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "v_xy",v_xy_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "w_xy",w_xy_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "t_xy",t_xy_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "q_xy",q_xy_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "LWP",LWP_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "surf_precip",surf_precip_vid) )

      !xz slices
      call netcdf_check( nf90_inq_varid(ncid_viz, "u_xz",u_xz_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "v_xz",v_xz_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "w_xz",w_xz_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "t_xz",t_xz_vid) )
      call netcdf_check( nf90_inq_varid(ncid_viz, "q_xz",q_xz_vid) )

end subroutine netcdf_res_viz

subroutine write_viz_netcdf
      use netcdf
      use pars
      use fields
      use con_data
      use con_stats
      use particles
      implicit none
      include 'mpif.h'

      integer :: yzslice,xyslice,xzslice
      integer :: ix,iy
      real :: tmpyz(nny,nnz),tmpxy(nnx,nny),tmpxz(nnx,nnz)
      real :: xvec(nnx),yvec(nny)
 

      if (myid==0) write(*,'(a35,i)') 'WRITE VIZ FILE, VIZ_COUNTER = ',viz_counter

      !Fill the grid on the first time through
      if (viz_counter == 1) then

      if (myid==0) then

      do ix=1,nnx
         xvec(ix) = dx*(ix-1)
      end do
      do iy=1,nny
         yvec(iy) = dy*(iy-1)
      end do

      call netcdf_check( nf90_put_var(ncid_viz, xgrid_vid, real(xvec),start=(/1/)) )
      call netcdf_check( nf90_put_var(ncid_viz, ygrid_vid, real(yvec),start=(/1/)) )
      call netcdf_check( nf90_put_var(ncid_viz, zugrid_vid, real(zz(1:nnz)),start=(/1/)) )
      call netcdf_check( nf90_put_var(ncid_viz, zwgrid_vid, real(z(0:nnz)),start=(/1/)) )
          
      end if
      end if

      !Write the time
      if (myid==0) then
      call netcdf_check( nf90_put_var(ncid_viz, time_viz_vid, real(time),start=(/viz_counter/)) )
      end if

      !Cycle through yz slices
      tmpyz = 0.0

      !Define the slices that should be sent
      yzslice = nnx/2  !Choose the nz value for the xy slice

    
      call fill_yz_slice(yzslice,u(1:nnx,iys:iye,izs:ize),tmpyz)
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, u_yz_vid,real(tmpyz(1:nny,1:nnz)),start=(/1,1,viz_counter/)) ) 

      call fill_yz_slice(yzslice,t(1:nnx,iys:iye,1,izs:ize),tmpyz)
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, t_yz_vid,real(tmpyz(1:nny,1:nnz)),start=(/1,1,viz_counter/)) ) 

      call fill_yz_slice(yzslice,t(1:nnx,iys:iye,2,izs:ize),tmpyz)
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, q_yz_vid,real(tmpyz(1:nny,1:nnz)),start=(/1,1,viz_counter/)) ) 

      call fill_yz_slice(yzslice,w(1:nnx,iys:iye,izs:ize),tmpyz)
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, w_yz_vid,real(tmpyz(1:nny,1:nnz)),start=(/1,1,viz_counter/)) ) 

      call fill_yz_slice(yzslice,v(1:nnx,iys:iye,izs:ize),tmpyz)
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, v_yz_vid,real(tmpyz(1:nny,1:nnz)),start=(/1,1,viz_counter/)) ) 


      !Now do xy slices
      tmpxy = 0.0
      xyslice = nnz/2  !Choose the nz value for the xy slice

      call fill_xy_slice(xyslice,u(1:nnx,iys:iye,izs:ize),tmpxy)
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, u_xy_vid,real(tmpxy(1:nnx,1:nny)),start=(/1,1,viz_counter/)) ) 

      call fill_xy_slice(xyslice,v(1:nnx,iys:iye,izs:ize),tmpxy)
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, v_xy_vid,real(tmpxy(1:nnx,1:nny)),start=(/1,1,viz_counter/)) ) 

      call fill_xy_slice(xyslice,w(1:nnx,iys:iye,izs:ize),tmpxy)
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, w_xy_vid,real(tmpxy(1:nnx,1:nny)),start=(/1,1,viz_counter/)) ) 

      call fill_xy_slice(xyslice,t(1:nnx,iys:iye,1,izs:ize),tmpxy)
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, t_xy_vid,real(tmpxy(1:nnx,1:nny)),start=(/1,1,viz_counter/)) ) 

      call fill_xy_slice(xyslice,t(1:nnx,iys:iye,2,izs:ize),tmpxy)
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, q_xy_vid,real(tmpxy(1:nnx,1:nny)),start=(/1,1,viz_counter/)) ) 

      !Now do xz slices
      tmpxz = 0.0
      xzslice = nny/2  !Choose the nz value for the xy slice

      call fill_xz_slice(xzslice,u(1:nnx,iys:iye,izs:ize),tmpxz)
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, u_xz_vid,real(tmpxz(1:nnx,1:nnz)),start=(/1,1,viz_counter/)) ) 

      call fill_xz_slice(xzslice,v(1:nnx,iys:iye,izs:ize),tmpxz)
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, v_xz_vid,real(tmpxz(1:nnx,1:nnz)),start=(/1,1,viz_counter/)) ) 

      call fill_xz_slice(xzslice,w(1:nnx,iys:iye,izs:ize),tmpxz)
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, w_xz_vid,real(tmpxz(1:nnx,1:nnz)),start=(/1,1,viz_counter/)) ) 

      call fill_xz_slice(xzslice,t(1:nnx,iys:iye,1,izs:ize),tmpxz)
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, t_xz_vid,real(tmpxz(1:nnx,1:nnz)),start=(/1,1,viz_counter/)) ) 

      call fill_xz_slice(xzslice,t(1:nnx,iys:iye,2,izs:ize),tmpxz)
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, q_xz_vid,real(tmpxz(1:nnx,1:nnz)),start=(/1,1,viz_counter/)) ) 


      !Now do the 2D particle stats
      call fill_xy_part_slice(LWP(mxs:mxe,iys:iye),tmpxy)
      tmpxy = tmpxy/dx/dy
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, LWP_vid,real(tmpxy(1:nnx,1:nny)),start=(/1,1,viz_counter/)) ) 

      call fill_xy_part_slice(surf_precip(mxs:mxe,iys:iye),tmpxy)
      tmpxy = tmpxy/dx/dy/viz_t_elapsed
      if (myid==0)  call netcdf_check( nf90_put_var(ncid_viz, surf_precip_vid,real(tmpxy(1:nnx,1:nny)),start=(/1,1,viz_counter/)) ) 


      surf_precip = 0.0
      viz_t_elapsed = 0.0
      viz_counter = viz_counter + 1


end subroutine write_viz_netcdf


subroutine close_his_netcdf
      use netcdf
      use pars
      implicit none

      if (myid==0) call netcdf_check( nf90_close(ncid) )

end subroutine close_his_netcdf

subroutine close_histog_netcdf
      use netcdf
      use pars
      implicit none

      if (myid==0) call netcdf_check( nf90_close(ncid_histog) )

end subroutine close_histog_netcdf

subroutine close_viz_netcdf
      use netcdf
      use pars
      implicit none

      if (myid==0) call netcdf_check( nf90_close(ncid_viz) )

end subroutine close_viz_netcdf

subroutine open_his_netcdf
      use netcdf
      use pars
      implicit none

      if (myid==0) call netcdf_check( nf90_open(path_netcdf_his,NF90_WRITE,ncid) )

end subroutine open_his_netcdf

subroutine open_histog_netcdf
      use netcdf
      use pars
      implicit none

      if (myid==0) call netcdf_check( nf90_open(path_netcdf_histog,NF90_WRITE,ncid_histog) )

end subroutine open_histog_netcdf

subroutine open_viz_netcdf
      use netcdf
      use pars
      implicit none
      include 'mpif.h'

      if (myid==0) call netcdf_check( nf90_open(path_netcdf_viz,NF90_WRITE,ncid_viz) )

end subroutine open_viz_netcdf

subroutine fill_yz_slice(yzslice,dat,tmp)
!Does all the MPI communication to fill a slice on proc 0 that can be written
      use pars
      use fields
      use particles
      implicit none
      include 'mpif.h'
      integer :: yzslice
      real,intent(in) :: dat(1:nnx,iys:iye,izs:ize)
      real,intent(out) :: tmp(nny,nnz)
      integer :: ierr,ip,mynx,myny,mynz,iystmp,iyetmp,izstmp,izetmp,istatus(MPI_STATUS_SIZE)
      integer :: sbuf_limits(4),rbuf_limits(4)
      real,allocatable :: sbuf_data(:,:),rbuf_data(:,:)

      !Now start the process of collecting the slices on myid==0

!      if (xyslice .le. ize .and. xyslice .ge. izs) then

      !For the yz slice, all processors will have to send to root:
      if (myid==0) then

         tmp(iys:iye,izs:ize) = dat(yzslice,iys:iye,izs:ize)

         do ip = 1,numprocs-1

            call mpi_recv(rbuf_limits,4,mpi_integer,ip,2,mpi_comm_world,istatus,ierr)

            myny = (rbuf_limits(2)-rbuf_limits(1))+1
            mynz = (rbuf_limits(4)-rbuf_limits(3))+1

            iystmp = rbuf_limits(1)
            iyetmp = rbuf_limits(2)
            izstmp = rbuf_limits(3)
            izetmp = rbuf_limits(4)

            allocate(rbuf_data(myny,mynz))

            call mpi_recv(rbuf_data,myny*mynz,mpi_real8,ip,2,mpi_comm_world,istatus,ierr)

            tmp(iystmp:iyetmp,izstmp:izetmp) = rbuf_data
            
            deallocate(rbuf_data)
         end do
         

      else  !All other processors

         mynx = nnx
         myny = (iye-iys)+1
         mynz = (ize-izs)+1

         sbuf_limits = (/iys,iye,izs,ize/)
         call mpi_send(sbuf_limits,4,mpi_integer,0,2,mpi_comm_world,ierr)


         allocate(sbuf_data(myny,mynz))

         sbuf_data = dat(yzslice,iys:iye,izs:ize)
         call mpi_send(sbuf_data,mynz*myny,mpi_real8,0,2,mpi_comm_world,ierr)
         
         deallocate(sbuf_data)

      end if

     !Now tmp contains the whole slice on processor 0

end subroutine fill_yz_slice

subroutine fill_xy_slice(xyslice,dat,tmp)
!Does all the MPI communication to fill a slice on proc 0 that can be written
      use pars
      use fields
      use particles
      implicit none
      include 'mpif.h'
      integer :: xyslice
      real,intent(in) :: dat(1:nnx,iys:iye,izs:ize)
      real,intent(out) :: tmp(nnx,nny)
      integer :: ierr,ip,mynx,myny,mynz,iystmp,iyetmp,izstmp,izetmp,istatus(MPI_STATUS_SIZE)
      integer :: sbuf_limits(4),rbuf_limits(4)
      real,allocatable :: sbuf_data(:,:),rbuf_data(:,:)

      !Now start the process of collecting the slices on myid==0

      !For the xy slice, not all processors participate
      if (myid==0) then

         if (xyslice .le. ize .and. xyslice .ge. izs) then
            tmp(1:nnx,iys:iye) = dat(1:nnx,iys:iye,xyslice)
         end if

         do ip = 1,numprocs-1

            call mpi_recv(rbuf_limits,4,mpi_integer,ip,2,mpi_comm_world,istatus,ierr)

            if (maxval(rbuf_limits) .ge. 0) then  !This is how it checks whether to accept data based on slice location

            mynx = nnx
            myny = (rbuf_limits(2)-rbuf_limits(1))+1
            mynz = (rbuf_limits(4)-rbuf_limits(3))+1

            iystmp = rbuf_limits(1)
            iyetmp = rbuf_limits(2)
            izstmp = rbuf_limits(3)
            izetmp = rbuf_limits(4)

            allocate(rbuf_data(mynx,myny))

            call mpi_recv(rbuf_data,mynx*myny,mpi_real8,ip,2,mpi_comm_world,istatus,ierr)

            tmp(1:nnx,iystmp:iyetmp) = rbuf_data
            
            deallocate(rbuf_data)

            end if

         end do
         

      else  !All other processors

         mynx = nnx
         myny = (iye-iys)+1
         mynz = (ize-izs)+1

         if (xyslice .le. ize .and. xyslice .ge. izs) then
            sbuf_limits = (/iys,iye,izs,ize/)
         else
            sbuf_limits = (/-1,-1,-1,-1/)
         end if

         call mpi_send(sbuf_limits,4,mpi_integer,0,2,mpi_comm_world,ierr)

         if (xyslice .le. ize .and. xyslice .ge. izs) then  !Only send to root if you have this slice

         allocate(sbuf_data(mynx,myny))

         sbuf_data = dat(1:nnx,iys:iye,xyslice)
         call mpi_send(sbuf_data,mynx*myny,mpi_real8,0,2,mpi_comm_world,ierr)
         
         deallocate(sbuf_data)

         end if

      end if

     !Now tmp contains the whole slice on processor 0

end subroutine fill_xy_slice

subroutine fill_xy_part_slice(dat,tmp)
!Does all the MPI communication to fill a slice on proc 0 that can be written
!This is specifically for particle stats decomposed across the XY plane
      use pars
      use fields
      use particles
      implicit none
      include 'mpif.h'
      real,intent(in) :: dat(mxs:mxe,iys:iye)
      real,intent(out) :: tmp(nnx,nny)
      integer :: ierr,ip,mynx,myny,mynz,iystmp,iyetmp,ixstmp,ixetmp,istatus(MPI_STATUS_SIZE)
      integer :: sbuf_limits(4),rbuf_limits(4)
      real,allocatable :: sbuf_data(:,:),rbuf_data(:,:)

      !Now start the process of collecting the slices on myid==0

      !For the xy slice, not all processors participate
      if (myid==0) then

	 !Proc 0 puts its data in the full array
         tmp(mxs:mxe,iys:iye) = dat(mxs:mxe,iys:iye)

         do ip = 1,numprocs-1

            call mpi_recv(rbuf_limits,4,mpi_integer,ip,2,mpi_comm_world,istatus,ierr)

            mynx = (rbuf_limits(2)-rbuf_limits(1))+1
            myny = (rbuf_limits(4)-rbuf_limits(3))+1

            ixstmp = rbuf_limits(1)
            ixetmp = rbuf_limits(2)
            iystmp = rbuf_limits(3)
            iyetmp = rbuf_limits(4)

            allocate(rbuf_data(mynx,myny))

            call mpi_recv(rbuf_data,mynx*myny,mpi_real8,ip,2,mpi_comm_world,istatus,ierr)

            tmp(ixstmp:ixetmp,iystmp:iyetmp) = rbuf_data
            
            deallocate(rbuf_data)

         end do
         

      else  !All other processors

         mynx = (mxe-mxs)+1
         myny = (iye-iys)+1

         sbuf_limits = (/mxs,mxe,iys,iye/)

         call mpi_send(sbuf_limits,4,mpi_integer,0,2,mpi_comm_world,ierr)

         allocate(sbuf_data(mynx,myny))

         sbuf_data = dat(mxs:mxe,iys:iye)
         call mpi_send(sbuf_data,mynx*myny,mpi_real8,0,2,mpi_comm_world,ierr)
         
         deallocate(sbuf_data)

      end if

     !Now tmp contains the whole slice on processor 0

end subroutine fill_xy_part_slice

subroutine fill_xz_slice(xzslice,dat,tmp)
!Does all the MPI communication to fill a slice on proc 0 that can be written
      use pars
      use fields
      use particles
      implicit none
      include 'mpif.h'
      integer :: xzslice
      real,intent(in) :: dat(1:nnx,iys:iye,izs:ize)
      real,intent(out) :: tmp(nnx,nnz)
      integer :: ierr,ip,mynx,myny,mynz,iystmp,iyetmp,izstmp,izetmp,istatus(MPI_STATUS_SIZE)
      integer :: sbuf_limits(4),rbuf_limits(4)
      real,allocatable :: sbuf_data(:,:),rbuf_data(:,:)

      !Now start the process of collecting the slices on myid==0

      !For the xy slice, not all processors participate
      if (myid==0) then

         if (xzslice .le. iye .and. xzslice .ge. iys) then
            tmp(1:nnx,izs:ize) = dat(1:nnx,xzslice,izs:ize)
         end if

         do ip = 1,numprocs-1

            call mpi_recv(rbuf_limits,4,mpi_integer,ip,2,mpi_comm_world,istatus,ierr)

            if (maxval(rbuf_limits) .ge. 0) then  !This is how it checks whether to accept data based on slice location

            mynx = nnx
            myny = (rbuf_limits(2)-rbuf_limits(1))+1
            mynz = (rbuf_limits(4)-rbuf_limits(3))+1

            iystmp = rbuf_limits(1)
            iyetmp = rbuf_limits(2)
            izstmp = rbuf_limits(3)
            izetmp = rbuf_limits(4)

            allocate(rbuf_data(mynx,mynz))

            call mpi_recv(rbuf_data,mynx*mynz,mpi_real8,ip,2,mpi_comm_world,istatus,ierr)

            tmp(1:nnx,izstmp:izetmp) = rbuf_data
            
            deallocate(rbuf_data)

            end if

         end do
         

      else  !All other processors

         mynx = nnx
         myny = (iye-iys)+1
         mynz = (ize-izs)+1

         if (xzslice .le. iye .and. xzslice .ge. iys) then
            sbuf_limits = (/iys,iye,izs,ize/)
         else
            sbuf_limits = (/-1,-1,-1,-1/)
         end if

         call mpi_send(sbuf_limits,4,mpi_integer,0,2,mpi_comm_world,ierr)

         if (xzslice .le. iye .and. xzslice .ge. iys) then  !Only send to root if you have this slice

         allocate(sbuf_data(mynx,mynz))

         sbuf_data = dat(1:nnx,xzslice,izs:ize)
         call mpi_send(sbuf_data,mynx*mynz,mpi_real8,0,2,mpi_comm_world,ierr)
         
         deallocate(sbuf_data)

         end if

      end if

     !Now tmp contains the whole slice on processor 0

end subroutine fill_xz_slice



subroutine write_histograms
use particles
use pars
use fields
use con_data
use con_stats
implicit none

include 'mpif.h'
integer :: i,j
integer :: ierr
logical :: isfirst

integer :: iblnk
integer :: num_entries

real :: sumbuf_rad(histbins+2)
real :: sumbuf_raddeath(histbins+2)
real :: sumbuf_res(histbins+2)
real :: sumbuf_actres(histbins+2)
real :: sumbuf_acttodeath(histbins+2)
real :: sumbuf_numact(histbins+2)

character :: fformat*10, num*3



! -------------- Collect all histogram from different processors

  num_entries = (histbins+2)

  !Radius PDF
  sumbuf_rad = hist_rad
  call mpi_reduce(sumbuf_rad,hist_rad,num_entries,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)

  !Radius PDF at death
  sumbuf_raddeath = hist_raddeath
  call mpi_reduce(sumbuf_raddeath,hist_raddeath,num_entries,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)

  !Residence time PDF
  sumbuf_res = hist_res
  call mpi_reduce(sumbuf_res,hist_res,num_entries,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)

  !Activation residence time PDF
  sumbuf_actres = hist_actres
  call mpi_reduce(sumbuf_actres,hist_actres,num_entries,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)

  !Activation residence time until death PDF
  sumbuf_acttodeath = hist_acttodeath
  call mpi_reduce(sumbuf_acttodeath,hist_acttodeath,num_entries,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)

  !# of activations PDF
  sumbuf_numact = hist_numact
  call mpi_reduce(sumbuf_numact,hist_numact,num_entries,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)


  if (myid==0) then

! ---------------- save data 

  if (inetcdf .eq. 1) then

      call write_histog_netcdf

  else

  write(num,'(i3.3)') histbins+2
  fformat = '('//num//'e15.6)'

  write(nrad,fformat),hist_rad(1:histbins+2)
  write(nrad,fformat),hist_raddeath(1:histbins+2)
  write(nres,fformat),hist_res(1:histbins+2)
  write(nactres,fformat),hist_actres(1:histbins+2)

  end if

  end if

  hist_res = 0.0  !Reset this one every time it writes
  hist_actres = 0.0  !Reset this one every time it writes
  hist_acttodeath = 0.0  !Reset this one every time it writes
  hist_numact = 0.0  !Reset this one every time it writes
  hist_raddeath = 0.0 !Reset this one every time it writes

end subroutine write_histograms

end module netcdf_io
