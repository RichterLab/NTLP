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
  real, allocatable :: ufsum_t(:,:,:,:),ufsqrsum_t(:,:,:,:)
  real, allocatable :: upwp_t(:,:,:),upwp(:,:,:)
  real, allocatable :: partcount(:,:,:),partsrc(:,:,:,:)
  real, allocatable :: vpsum(:,:,:,:),vpsqrsum(:,:,:,:)
  real, allocatable :: ufsum(:,:,:,:),ufsqrsum(:,:,:,:)
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

  !--- SFS velocity calculation ---------
  real, allocatable :: sigm_s(:,:,:),sigm_sdx(:,:,:),sigm_sdy(:,:,:)
  real, allocatable :: sigm_sdz(:,:,:),sigm_sext(:,:,:)
  real, allocatable :: sigm_sdxext(:,:,:),sigm_sdyext(:,:,:)
  real, allocatable :: sigm_sdzext(:,:,:)
  real, allocatable :: vis_ss(:,:,:),vis_sext(:,:,:)

  integer :: particletype,pad_diff
  integer :: numpart,tnumpart,ngidx
  integer :: numdrop,tnumdrop
  integer :: numaerosol,tnumaerosol
  integer :: iseed
  integer :: num100=0, num1000=0, numimpos=0
  integer :: tnum100, tnum1000, tnumimpos
  integer :: denum, actnum, tdenum, tactnum
  integer :: num_destroy=0,tnum_destroy=0
  integer :: tot_reintro=0

  real :: Rep_avg,part_grav(3)
  real :: radavg,radmin,radmax,radmsqr,tempmin,tempmax,qmin,qmax
  real :: radavg_center,radmsqr_center
  real :: vp_init(3),Tp_init,radius_init,radius_std,kappas_init,kappas_std
  real :: pdf_factor,pdf_prob
  integer*8 :: mult_init,mult_factor,mult_a,mult_c

  real,parameter :: Cvv=1463.0
  real,parameter :: Cpv = 1952.0
  real,parameter :: Cva = 717.04

  real :: avgres=0,tavgres=0

  integer, parameter :: histbins = 512
  real :: hist_rad(histbins+2)
  real :: hist_raddeath(histbins+2)
  real :: bins_rad(histbins+2)

  real :: hist_res(histbins+2)
  real :: bins_res(histbins+2)

  real :: hist_actres(histbins+2)
  real :: bins_actres(histbins+2)

  real :: hist_acttodeath(histbins+2)
  real :: bins_acttodeath(histbins+2)

  real :: hist_numact(histbins+2)
  real :: bins_numact(histbins+2)

  !REMEMBER: IF ADDING ANYTHING, MUST UPDATE MPI DATATYPE!
  type :: particle
    integer :: pidx,procidx,nbr_pidx,nbr_procidx
    real :: vp(3),xp(3),uf(3),xrhs(3),vrhs(3),Tp,Tprhs_s
    real :: Tprhs_L,Tf,radius,radrhs,qinf,qstar,dist
    real :: res,m_s,kappa_s,rc,actres,numact
    real :: u_sub(3),sigm_s
    integer*8 :: mult
    type(particle), pointer :: prev,next
  end type particle

  type(particle), pointer :: part,first_particle

CONTAINS

  subroutine fill_ext
    use pars
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
    T2ext = 0.0

    !First fill the center, since this is just u,v,w on that proc:

    !In the column setup, need to tranpose u,v,w first into u_t,v_t,w_t:
    call xtoz_trans(u(1:nnx,iys:iye,izs-1:ize+1),u_t,nnx,nnz, &
    mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e, &
    myid,ncpu_s,numprocs)
    call xtoz_trans(v(1:nnx,iys:iye,izs-1:ize+1),v_t,nnx,nnz, &
    mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e, &
    myid,ncpu_s,numprocs)
    call xtoz_trans(w(1:nnx,iys:iye,izs-1:ize+1),w_t,nnx,nnz, &
    mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e, &
    myid,ncpu_s,numprocs)
    call xtoz_trans(t(1:nnx,iys:iye,1,izs-1:ize+1),T_t,nnx,nnz, &
    mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e, &
    myid,ncpu_s,numprocs)
    call xtoz_trans(t(1:nnx,iys:iye,2,izs-1:ize+1),T2_t,nnx,nnz, &
    mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e, &
    myid,ncpu_s,numprocs)

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
    tbuf_s(1:nnz+2,1:iye-iys+1,1:2,5) = T2_t(0:nnz+1,iys:iye,mxe-1:mxe)
    trbuf_s(1:nnz+2,1:2,1:2,5) = T2_t(0:nnz+1,iye-1:iye,mxe-1:mxe)
    rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,5) = T2_t(0:nnz+1,iye-1:iye,mxs:mxe)
    brbuf_s(1:nnz+2,1:2,1:3,5) = T2_t(0:nnz+1,iye-1:iye,mxs:mxs+2)
    bbuf_s(1:nnz+2,1:iye-iys+1,1:3,5) = T2_t(0:nnz+1,iys:iye,mxs:mxs+2)
    blbuf_s(1:nnz+2,1:3,1:3,5) = T2_t(0:nnz+1,iys:iys+2,mxs:mxs+2)
    lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,5) = T2_t(0:nnz+1,iys:iys+2,mxs:mxe)
    tlbuf_s(1:nnz+2,1:3,1:2,5) = T2_t(0:nnz+1,iys:iys+2,mxe-1:mxe)

    !Zero out recieve buffers
    rbuf_r=0.0;trbuf_r=0.0;tbuf_r=0.0;tlbuf_r=0.0;lbuf_r=0.0
    blbuf_r=0.0;bbuf_r=0.0;brbuf_r=0.0

    !Left/right:
    call MPI_Sendrecv(rbuf_s,rc_s,mpi_real8,rproc,3, &
    lbuf_r,lc_r,mpi_real8,lproc,3,mpi_comm_world,istatus,ierr)

    call mpi_barrier(mpi_comm_world,ierr)
    call MPI_Sendrecv(lbuf_s,lc_s,mpi_real8,lproc,4, &
    rbuf_r,rc_r,mpi_real8,rproc,4,mpi_comm_world,istatus,ierr)

    !Top/bottom:
    call MPI_Sendrecv(tbuf_s,tc_s,mpi_real8,tproc,5, &
    bbuf_r,bc_r,mpi_real8,bproc,5,mpi_comm_world,istatus,ierr)

    call MPI_Sendrecv(bbuf_s,bc_s,mpi_real8,bproc,6, &
    tbuf_r,tc_r,mpi_real8,tproc,6,mpi_comm_world,istatus,ierr)

    !Top right/bottom left:
    call MPI_Sendrecv(trbuf_s,trc_s,mpi_real8,trproc,7, &
    blbuf_r,blc_r,mpi_real8,blproc,7, &
    mpi_comm_world,istatus,ierr)

    call MPI_Sendrecv(blbuf_s,blc_s,mpi_real8,blproc,8, &
    trbuf_r,trc_r,mpi_real8,trproc,8, &
    mpi_comm_world,istatus,ierr)

    !Top left/bottom right:
    call MPI_Sendrecv(tlbuf_s,tlc_s,mpi_real8,tlproc,9, &
    brbuf_r,brc_r,mpi_real8,brproc,9, &
    mpi_comm_world,istatus,ierr)

    call MPI_Sendrecv(brbuf_s,brc_s,mpi_real8,brproc,10, &
    tlbuf_r,tlc_r,mpi_real8,tlproc,10, &
    mpi_comm_world,istatus,ierr)

    !Now fill the ext arrays with the recieved buffers:
    uext(0:nnz+1,iys:iye,mxe+1:mxe+3) = tbuf_r(1:nnz+2,1:iye-iys+1,1:3,1)
    uext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,1)
    uext(0:nnz+1,iye+1:iye+3,mxs:mxe) = rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,1)
    uext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,1)
    uext(0:nnz+1,iys:iye,mxs-2:mxs-1) = bbuf_r(1:nnz+2,1:iye-iys+1,1:2,1)
    uext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,1)
    uext(0:nnz+1,iys-2:iys-1,mxs:mxe) = lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,1)
    uext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,1)

    vext(0:nnz+1,iys:iye,mxe+1:mxe+3) = tbuf_r(1:nnz+2,1:iye-iys+1,1:3,2)
    vext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,2)
    vext(0:nnz+1,iye+1:iye+3,mxs:mxe) = rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,2)
    vext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,2)
    vext(0:nnz+1,iys:iye,mxs-2:mxs-1) = bbuf_r(1:nnz+2,1:iye-iys+1,1:2,2)
    vext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,2)
    vext(0:nnz+1,iys-2:iys-1,mxs:mxe) = lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,2)
    vext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,2)

    wext(0:nnz+1,iys:iye,mxe+1:mxe+3) = tbuf_r(1:nnz+2,1:iye-iys+1,1:3,3)
    wext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,3)
    wext(0:nnz+1,iye+1:iye+3,mxs:mxe) = rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,3)
    wext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,3)
    wext(0:nnz+1,iys:iye,mxs-2:mxs-1) = bbuf_r(1:nnz+2,1:iye-iys+1,1:2,3)
    wext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,3)
    wext(0:nnz+1,iys-2:iys-1,mxs:mxe) = lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,3)
    wext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,3)

    Text(0:nnz+1,iys:iye,mxe+1:mxe+3) = tbuf_r(1:nnz+2,1:iye-iys+1,1:3,4)
    Text(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,4)
    Text(0:nnz+1,iye+1:iye+3,mxs:mxe) = rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,4)
    Text(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,4)
    Text(0:nnz+1,iys:iye,mxs-2:mxs-1) = bbuf_r(1:nnz+2,1:iye-iys+1,1:2,4)
    Text(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,4)
    Text(0:nnz+1,iys-2:iys-1,mxs:mxe) = lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,4)
    Text(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,4)


    T2ext(0:nnz+1,iys:iye,mxe+1:mxe+3) = tbuf_r(1:nnz+2,1:iye-iys+1,1:3,5)
    T2ext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,5)
    T2ext(0:nnz+1,iye+1:iye+3,mxs:mxe) = rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,5)
    T2ext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,5)
    T2ext(0:nnz+1,iys:iye,mxs-2:mxs-1) = bbuf_r(1:nnz+2,1:iye-iys+1,1:2,5)
    T2ext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,5)
    T2ext(0:nnz+1,iys-2:iys-1,mxs:mxe) = lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,5)
    T2ext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,5)


  end subroutine fill_ext
  subroutine fill_extSFS
    !       This subroutine calculte the extented fields for SFS velocity

    use pars
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

    call xtoz_trans(sigm_s(1:nnx,iys:iye,izs-1:ize+1),sigm_st,nnx, &
    nnz,mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e, &
    myid,ncpu_s,numprocs)
    call xtoz_trans(sigm_sdx(1:nnx,iys:iye,izs-1:ize+1),sigm_sdxt,nnx, &
    nnz,mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e, &
    myid,ncpu_s,numprocs)
    call xtoz_trans(sigm_sdy(1:nnx,iys:iye,izs-1:ize+1),sigm_sdyt,nnx, &
    nnz,mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e, &
    myid,ncpu_s,numprocs)
    call xtoz_trans(sigm_sdz(1:nnx,iys:iye,izs-1:ize+1),sigm_sdzt,nnx, &
    nnz,mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e, &
    myid,ncpu_s,numprocs)

    call xtoz_trans(vis_ss(1:nnx,iys:iye,izs-1:ize+1),vis_st,nnx, &
    nnz,mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e, &
    myid,ncpu_s,numprocs)

    sigm_sext(0:nnz+1,iys:iye,mxs:mxe) = sigm_st(0:nnz+1,iys:iye,mxs:mxe)
    sigm_sdxext(0:nnz+1,iys:iye,mxs:mxe) = sigm_sdxt(0:nnz+1,iys:iye,mxs:mxe)
    sigm_sdyext(0:nnz+1,iys:iye,mxs:mxe) = sigm_sdyt(0:nnz+1,iys:iye,mxs:mxe)
    sigm_sdzext(0:nnz+1,iys:iye,mxs:mxe) = sigm_sdzt(0:nnz+1,iys:iye,mxs:mxe)
    vis_sext(0:nnz+1,iys:iye,mxs:mxe) = vis_st(0:nnz+1,iys:iye,mxs:mxe)



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
    tbuf_s(1:nnz+2,1:iye-iys+1,1:2,1) = sigm_st(0:nnz+1,iys:iye,mxe-1:mxe)
    trbuf_s(1:nnz+2,1:2,1:2,1) = sigm_st(0:nnz+1,iye-1:iye,mxe-1:mxe)
    rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,1) = sigm_st(0:nnz+1,iye-1:iye,mxs:mxe)
    brbuf_s(1:nnz+2,1:2,1:3,1) = sigm_st(0:nnz+1,iye-1:iye,mxs:mxs+2)
    bbuf_s(1:nnz+2,1:iye-iys+1,1:3,1) = sigm_st(0:nnz+1,iys:iye,mxs:mxs+2)
    blbuf_s(1:nnz+2,1:3,1:3,1) = sigm_st(0:nnz+1,iys:iys+2,mxs:mxs+2)
    lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,1) = sigm_st(0:nnz+1,iys:iys+2,mxs:mxe)
    tlbuf_s(1:nnz+2,1:3,1:2,1) = sigm_st(0:nnz+1,iys:iys+2,mxe-1:mxe)

    !sigm_sdx:
    tbuf_s(1:nnz+2,1:iye-iys+1,1:2,2) = sigm_sdxt(0:nnz+1,iys:iye,mxe-1:mxe)
    trbuf_s(1:nnz+2,1:2,1:2,2) = sigm_sdxt(0:nnz+1,iye-1:iye,mxe-1:mxe)
    rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,2) = sigm_sdxt(0:nnz+1,iye-1:iye,mxs:mxe)
    brbuf_s(1:nnz+2,1:2,1:3,2) = sigm_sdxt(0:nnz+1,iye-1:iye,mxs:mxs+2)
    bbuf_s(1:nnz+2,1:iye-iys+1,1:3,2) = sigm_sdxt(0:nnz+1,iys:iye,mxs:mxs+2)
    blbuf_s(1:nnz+2,1:3,1:3,2) = sigm_sdxt(0:nnz+1,iys:iys+2,mxs:mxs+2)
    lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,2) = sigm_sdxt(0:nnz+1,iys:iys+2,mxs:mxe)
    tlbuf_s(1:nnz+2,1:3,1:2,2) = sigm_sdxt(0:nnz+1,iys:iys+2,mxe-1:mxe)

    !sigm_sdy:
    tbuf_s(1:nnz+2,1:iye-iys+1,1:2,3) = sigm_sdyt(0:nnz+1,iys:iye,mxe-1:mxe)
    trbuf_s(1:nnz+2,1:2,1:2,3) = sigm_sdyt(0:nnz+1,iye-1:iye,mxe-1:mxe)
    rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,3) = sigm_sdyt(0:nnz+1,iye-1:iye,mxs:mxe)
    brbuf_s(1:nnz+2,1:2,1:3,3) = sigm_sdyt(0:nnz+1,iye-1:iye,mxs:mxs+2)
    bbuf_s(1:nnz+2,1:iye-iys+1,1:3,3) = sigm_sdyt(0:nnz+1,iys:iye,mxs:mxs+2)
    blbuf_s(1:nnz+2,1:3,1:3,3) = sigm_sdyt(0:nnz+1,iys:iys+2,mxs:mxs+2)
    lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,3) = sigm_sdyt(0:nnz+1,iys:iys+2,mxs:mxe)
    tlbuf_s(1:nnz+2,1:3,1:2,3) = sigm_sdyt(0:nnz+1,iys:iys+2,mxe-1:mxe)

    !sigm_sdz
    tbuf_s(1:nnz+2,1:iye-iys+1,1:2,4) = sigm_sdzt(0:nnz+1,iys:iye,mxe-1:mxe)
    trbuf_s(1:nnz+2,1:2,1:2,4) = sigm_sdzt(0:nnz+1,iye-1:iye,mxe-1:mxe)
    rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,4) = sigm_sdzt(0:nnz+1,iye-1:iye,mxs:mxe)
    brbuf_s(1:nnz+2,1:2,1:3,4) = sigm_sdzt(0:nnz+1,iye-1:iye,mxs:mxs+2)
    bbuf_s(1:nnz+2,1:iye-iys+1,1:3,4) = sigm_sdzt(0:nnz+1,iys:iye,mxs:mxs+2)
    blbuf_s(1:nnz+2,1:3,1:3,4) = sigm_sdzt(0:nnz+1,iys:iys+2,mxs:mxs+2)
    lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,4) = sigm_sdzt(0:nnz+1,iys:iys+2,mxs:mxe)
    tlbuf_s(1:nnz+2,1:3,1:2,4) = sigm_sdzt(0:nnz+1,iys:iys+2,mxe-1:mxe)

    !vis_s
    tbuf_s(1:nnz+2,1:iye-iys+1,1:2,5) = vis_st(0:nnz+1,iys:iye,mxe-1:mxe)
    trbuf_s(1:nnz+2,1:2,1:2,5) = vis_st(0:nnz+1,iye-1:iye,mxe-1:mxe)
    rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,5) = vis_st(0:nnz+1,iye-1:iye,mxs:mxe)
    brbuf_s(1:nnz+2,1:2,1:3,5) = vis_st(0:nnz+1,iye-1:iye,mxs:mxs+2)
    bbuf_s(1:nnz+2,1:iye-iys+1,1:3,5) = vis_st(0:nnz+1,iys:iye,mxs:mxs+2)
    blbuf_s(1:nnz+2,1:3,1:3,5)  = vis_st(0:nnz+1,iys:iys+2,mxs:mxs+2)
    lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,5) = vis_st(0:nnz+1,iys:iys+2,mxs:mxe)
    tlbuf_s(1:nnz+2,1:3,1:2,5) = vis_st(0:nnz+1,iys:iys+2,mxe-1:mxe)


    !Zero out recieve buffers
    rbuf_r=0.0;trbuf_r=0.0;tbuf_r=0.0;tlbuf_r=0.0;lbuf_r=0.0
    blbuf_r=0.0;bbuf_r=0.0;brbuf_r=0.0

    !Left/right:
    call MPI_Sendrecv(rbuf_s,rc_s,mpi_real8,rproc,3, &
    lbuf_r,lc_r,mpi_real8,lproc,3,mpi_comm_world,istatus,ierr)

    call MPI_Sendrecv(lbuf_s,lc_s,mpi_real8,lproc,4, &
    rbuf_r,rc_r,mpi_real8,rproc,4,mpi_comm_world,istatus,ierr)

    !Top/bottom:
    call MPI_Sendrecv(tbuf_s,tc_s,mpi_real8,tproc,5, &
    bbuf_r,bc_r,mpi_real8,bproc,5,mpi_comm_world,istatus,ierr)

    call MPI_Sendrecv(bbuf_s,bc_s,mpi_real8,bproc,6, &
    tbuf_r,tc_r,mpi_real8,tproc,6,mpi_comm_world,istatus,ierr)

    !Top right/bottom left:
    call MPI_Sendrecv(trbuf_s,trc_s,mpi_real8,trproc,7, &
    blbuf_r,blc_r,mpi_real8,blproc,7, &
    mpi_comm_world,istatus,ierr)

    call MPI_Sendrecv(blbuf_s,blc_s,mpi_real8,blproc,8, &
    trbuf_r,trc_r,mpi_real8,trproc,8, &
    mpi_comm_world,istatus,ierr)

    !Top left/bottom right:
    call MPI_Sendrecv(tlbuf_s,tlc_s,mpi_real8,tlproc,9, &
    brbuf_r,brc_r,mpi_real8,brproc,9, &
    mpi_comm_world,istatus,ierr)

    call MPI_Sendrecv(brbuf_s,brc_s,mpi_real8,brproc,10, &
    tlbuf_r,tlc_r,mpi_real8,tlproc,10, &
    mpi_comm_world,istatus,ierr)

    !Now fill the ext arrays with the recieved buffers:
    sigm_sext(0:nnz+1,iys:iye,mxe+1:mxe+3) = tbuf_r(1:nnz+2,1:iye-iys+1,1:3,1)
    sigm_sext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,1)
    sigm_sext(0:nnz+1,iye+1:iye+3,mxs:mxe) = rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,1)
    sigm_sext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,1)
    sigm_sext(0:nnz+1,iys:iye,mxs-2:mxs-1) = bbuf_r(1:nnz+2,1:iye-iys+1,1:2,1)
    sigm_sext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,1)
    sigm_sext(0:nnz+1,iys-2:iys-1,mxs:mxe) = lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,1)
    sigm_sext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,1)

    sigm_sdxext(0:nnz+1,iys:iye,mxe+1:mxe+3) = tbuf_r(1:nnz+2,1:iye-iys+1,1:3,2)
    sigm_sdxext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,2)
    sigm_sdxext(0:nnz+1,iye+1:iye+3,mxs:mxe) = rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,2)
    sigm_sdxext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,2)
    sigm_sdxext(0:nnz+1,iys:iye,mxs-2:mxs-1) = bbuf_r(1:nnz+2,1:iye-iys+1,1:2,2)
    sigm_sdxext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,2)
    sigm_sdxext(0:nnz+1,iys-2:iys-1,mxs:mxe) = lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,2)
    sigm_sdxext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,2)

    sigm_sdyext(0:nnz+1,iys:iye,mxe+1:mxe+3) = tbuf_r(1:nnz+2,1:iye-iys+1,1:3,3)
    sigm_sdyext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,3)
    sigm_sdyext(0:nnz+1,iye+1:iye+3,mxs:mxe) = rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,3)
    sigm_sdyext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,3)
    sigm_sdyext(0:nnz+1,iys:iye,mxs-2:mxs-1) = bbuf_r(1:nnz+2,1:iye-iys+1,1:2,3)
    sigm_sdyext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,3)
    sigm_sdyext(0:nnz+1,iys-2:iys-1,mxs:mxe) = lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,3)
    sigm_sdyext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,3)

    sigm_sdzext(0:nnz+1,iys:iye,mxe+1:mxe+3) = tbuf_r(1:nnz+2,1:iye-iys+1,1:3,4)
    sigm_sdzext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,4)
    sigm_sdzext(0:nnz+1,iye+1:iye+3,mxs:mxe) = rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,4)
    sigm_sdzext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,4)
    sigm_sdzext(0:nnz+1,iys:iye,mxs-2:mxs-1) = bbuf_r(1:nnz+2,1:iye-iys+1,1:2,4)
    sigm_sdzext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,4)
    sigm_sdzext(0:nnz+1,iys-2:iys-1,mxs:mxe) = lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,4)
    sigm_sdzext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,4)


    vis_sext(0:nnz+1,iys:iye,mxe+1:mxe+3) = tbuf_r(1:nnz+2,1:iye-iys+1,1:3,5)
    vis_sext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,5)
    vis_sext(0:nnz+1,iye+1:iye+3,mxs:mxe) = rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,5)
    vis_sext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,5)
    vis_sext(0:nnz+1,iys:iye,mxs-2:mxs-1) = bbuf_r(1:nnz+2,1:iye-iys+1,1:2,5)
    vis_sext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,5)
    vis_sext(0:nnz+1,iys-2:iys-1,mxs:mxe) = lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,5)
    vis_sext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,5)

  end subroutine fill_extSFS

  subroutine uf_interp
  use pars
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
      
  !get the "leftmost" node
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
  kuvpts(3) = minloc(zz,1,mask=(zz.gt.part%xp(3))) - 2
  !Then fill in the rest:
  kuvpts(4) = kuvpts(3)+1
  kuvpts(5) = kuvpts(4)+1
  kuvpts(6) = kuvpts(5)+1
  kuvpts(2) = kuvpts(3)-1
  kuvpts(1) = kuvpts(2)-1


  kwpts(3) = minloc(z,1,mask=(z.gt.part%xp(3))) - 2
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

        part%uf(1) = part%uf(1)+uext(izuv,iy,ix)*wt(1,i)*wt(2,j)*wt(3,k) 
        part%uf(2) = part%uf(2)+vext(izuv,iy,ix)*wt(1,i)*wt(2,j)*wt(3,k) 
        part%uf(3) = part%uf(3)+wext(izw,iy,ix)*wt(1,i)*wt(2,j)*wt(4,k) 
        part%Tf = part%Tf+Text(izuv,iy,ix)*wt(1,i)*wt(2,j)*wt(3,k)
        part%qinf = part%qinf+T2ext(izuv,iy,ix)*wt(1,i)*wt(2,j)*wt(3,k) 
     end do
     end do 
     end do

  end subroutine uf_interp 

  subroutine uf_interp_lin
  use pars
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
  
  kpt = minloc(zz,1,mask=(zz.gt.part%xp(3))) - 2
  kwpt = minloc(z,1,mask=(z.gt.part%xp(3))) - 2


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
     iy = jpt+j
     izuv = kpt+k
     izw = kwpt+k


     part%uf(1) = part%uf(1) + uext(izuv,iy,ix)*wtx*wty*wtz
     part%uf(2) = part%uf(2) + vext(izuv,iy,ix)*wtx*wty*wtz
     part%uf(3) = part%uf(3) + wext(izw,iy,ix)*wtx*wty*wtzw

     part%Tf = part%Tf + Text(izuv,iy,ix)*wtx*wty*wtz
     part%qinf = part%qinf + T2ext(izuv,iy,ix)*wtx*wty*wtz

  end do

  !Since the ghost points don't follow the filling of "ext", a quick
  !fix is to use 0th order interpolation if between last (uv) point
  !and wall
  if (kpt .eq. nnz) then
     part%uf(1) = uext(kpt,iy,ix)
     part%uf(2) = vext(kpt,iy,ix)
     part%Tf = Text(kpt,iy,ix)
     part%qinf = T2ext(kpt,iy,ix)
  end if
  if (kpt .eq. 0) then
     part%uf(1) = uext(1,iy,ix)
     part%uf(2) = vext(1,iy,ix)
     part%Tf = Text(1,iy,ix)
     part%qinf = T2ext(1,iy,ix)
  end if
  end do
  end do


  end subroutine uf_interp_lin

  subroutine sigm_interp(sigm_sdxp,sigm_sdyp,sigm_sdzp,vis_sp,iz_part)
  use pars
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
  !get the "leftmost" node
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

  !Do (u,v) first
  kuvpts(3) = minloc(zz,1,mask=(zz.gt.part%xp(3))) - 2
  !Then fill in the rest:
  kuvpts(4) = kuvpts(3)+1
  kuvpts(5) = kuvpts(4)+1
  kuvpts(6) = kuvpts(5)+1
  kuvpts(2) = kuvpts(3)-1
  kuvpts(1) = kuvpts(2)-1

  kwpts(3) = minloc(z,1,mask=(z.gt.part%xp(3))) - 2
  iz_part = kwpts(3)
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
   !Recall that wt has been set to zero, so weights will be zero if (first,last) isn't (1,6)
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

  !Now we have the weights - compute the sigma_s and derivatives at xp:
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
      part%sigm_s = part%sigm_s+sigm_sext(izw,iy,ix)*wt(1,i)*wt(2,j)*wt(4,k)
      sigm_sdxp = sigm_sdxp+ sigm_sdxext(izw,iy,ix)*wt(1,i)*wt(2,j)*wt(4,k)
      sigm_sdyp = sigm_sdyp+sigm_sdyext(izw,iy,ix)*wt(1,i)*wt(2,j)*wt(4,k)
      sigm_sdzp = sigm_sdzp+sigm_sdzext(izuv,iy,ix)*wt(1,i)*wt(2,j)*wt(3,k)

      vis_sp =vis_sp +vis_sext(izw,iy,ix)*wt(1,i)*wt(2,j)*wt(4,k)

  if (isnan(vis_sp)) then
        write(*,*) 'WARNING:',vis_sext(izw,iy,ix),wt(1,i),wt(2,j),wt(4,k),part%xp(3),i,j,k,kwpts(3)
  end if


  end do
  end do
  end do


  end subroutine sigm_interp

  subroutine particle_coupling_exchange
      use pars
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'
      real :: ctbuf_s(nnz+2,1:iye-iys+2,6),cbbuf_r(nnz+2,1:iye-iys+2,6)
      real :: crbuf_s(nnz+2,1:mxe-mxs+1,6),clbuf_r(nnz+2,1:mxe-mxs+1,6)
      integer :: istatus(mpi_status_size),ierr,ncount


      !Now, partsrc and partTsrc have halos on each processor - give these to the rightful owner:
      crbuf_s=0.0;ctbuf_s=0.0
      clbuf_r=0.0;cbbuf_r=0.0

      !First send top: 
      !get send buffer ready:
      ctbuf_s(1:nnz+2,1:iye-iys+2,1:3)=partsrc_t(0:nnz+1,iys:iye+1,mxe+1,1:3)
      ctbuf_s(1:nnz+2,1:iye-iys+2,4)=partTsrc_t(0:nnz+1,iys:iye+1,mxe+1)
      ctbuf_s(1:nnz+2,1:iye-iys+2,5)=partHsrc_t(0:nnz+1,iys:iye+1,mxe+1)
      ctbuf_s(1:nnz+2,1:iye-iys+2,6)=partTEsrc_t(0:nnz+1,iys:iye+1,mxe+1)

      ncount = 6*(nnz+2)*(iye-iys+2)
      call mpi_sendrecv(ctbuf_s,ncount,mpi_real8,tproc,1, &
           cbbuf_r,ncount,mpi_real8,bproc,1,mpi_comm_world,istatus,ierr)

      !Now just add the contents of the receive buffer into the entire iys column of this proc:

      partsrc_t(0:nnz+1,iys:iye+1,mxs,1:3) = partsrc_t(0:nnz+1,iys:iye+1,mxs,1:3) + cbbuf_r(1:nnz+2,1:iye-iys+2,1:3)
      partTsrc_t(0:nnz+1,iys:iye+1,mxs) = partTsrc_t(0:nnz+1,iys:iye+1,mxs) + cbbuf_r(1:nnz+2,1:iye-iys+2,4)
      partHsrc_t(0:nnz+1,iys:iye+1,mxs) = partHsrc_t(0:nnz+1,iys:iye+1,mxs) + cbbuf_r(1:nnz+2,1:iye-iys+2,5)
      partTEsrc_t(0:nnz+1,iys:iye+1,mxs) = partTEsrc_t(0:nnz+1,iys:iye+1,mxs) + cbbuf_r(1:nnz+2,1:iye-iys+2,6)

      !Now get the right send buffer ready:
      crbuf_s(1:nnz+2,1:mxe-mxs+1,1:3)=partsrc_t(0:nnz+1,iye+1,mxs:mxe,1:3)
      crbuf_s(1:nnz+2,1:mxe-mxs+1,4)=partTsrc_t(0:nnz+1,iye+1,mxs:mxe)
      crbuf_s(1:nnz+2,1:mxe-mxs+1,5)=partHsrc_t(0:nnz+1,iye+1,mxs:mxe)
      crbuf_s(1:nnz+2,1:mxe-mxs+1,6)=partTEsrc_t(0:nnz+1,iye+1,mxs:mxe)

      !Now send to right:
      ncount = 6*(nnz+2)*(mxe-mxs+1)
      call mpi_sendrecv(crbuf_s,ncount,mpi_real8,rproc,2, &
           clbuf_r,ncount,mpi_real8,lproc,2,mpi_comm_world,istatus,ierr)

      !And again add the contents to the top/bottom rows of partsrc:
      partsrc_t(0:nnz+1,iys,mxs:mxe,1:3) = partsrc_t(0:nnz+1,iys,mxs:mxe,1:3) + clbuf_r(1:nnz+2,1:mxe-mxs+1,1:3)

      partTsrc_t(0:nnz+1,iys,mxs:mxe) = partTsrc_t(0:nnz+1,iys,mxs:mxe) + clbuf_r(1:nnz+2,1:mxe-mxs+1,4)
      partHsrc_t(0:nnz+1,iys,mxs:mxe) = partHsrc_t(0:nnz+1,iys,mxs:mxe) + clbuf_r(1:nnz+2,1:mxe-mxs+1,5)
      partTEsrc_t(0:nnz+1,iys,mxs:mxe) = partTEsrc_t(0:nnz+1,iys,mxs:mxe) + clbuf_r(1:nnz+2,1:mxe-mxs+1,6)


  end subroutine particle_coupling_exchange

  subroutine particle_coupling_update
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
  kpt = minloc(zz,1,mask=(zz.gt.part%xp(3))) - 2
  kwpt = minloc(z,1,mask=(z.gt.part%xp(3))) - 1


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

     rhop = (part%m_s+pi2*2.0/3.0*part%radius**3*rhow)/(pi2*2.0/3.0*part%radius**3)
     partmass = rhop*2.0/3.0*pi2*(part%radius)**3
     taup_i = 18.0*rhoa*nuf/rhop/(2.0*part%radius)**2 !Brian 8/4/14

     ix = ipt+i
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
      partsrc_t(iz,iy,ix,1:3) = &
          partsrc_t(iz,iy,ix,1:3) - partmass/rhoa*(part%vrhs(1:3)-part_grav(1:3))*wtt/dV*real(part%mult)

      !vapor momentum coupling
      partsrc_t(iz,iy,ix,1:3) = &
          partsrc_t(iz,iy,ix,1:3) - rhow/rhoa*pi2*2*part%radius**2*part%radrhs*part%vp(1:3)*wtt/dV*real(part%mult)
     endif

     if (iTcouple == 1) then
      partTsrc_t(iz,iy,ix) = &
          partTsrc_t(iz,iy,ix) - (part%Tprhs_s*6.0*rhow/rhop/CpaCpp/taup_i*(pi2/2.0)*part%radius*nuf)*wtt/dV*real(part%mult)
     endif

     if (iHcouple == 1) then
      partHsrc_t(iz,iy,ix) = &
          partHsrc_t(iz,iy,ix) - rhow/rhoa*pi2*2*part%radius**2*part%radrhs*wtt/dV*real(part%mult)


      partTEsrc_t(iz,iy,ix) = &
          partTEsrc_t(iz,iy,ix) - rhow/rhoa*pi2*2*part%radius**2*part%radrhs*Cpv/Cpa*part%Tp*wtt/dV*real(part%mult) + &
              rhow/rhoa*pi2*2*part%radius**2*part%radrhs*Cpv/Cpa*part%Tf*wtt/dV*real(part%mult)

     endif


     end do
     end do
     end do

     part => part%next
  end do

  end subroutine particle_coupling_update

  subroutine assign_nbrs
        use pars
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

  subroutine particle_exchange
      use pars
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
         elseif ( (part%xp(1) .GT. xmax) .AND. &
                  (part%xp(2) .LT. ymax) .AND. &
                  (part%xp(2) .GT. ymin) ) then !top
            pt_s = pt_s + 1
         elseif ( (part%xp(1) .LT. xmin) .AND. &
                  (part%xp(2) .LT. ymax) .AND. &
                  (part%xp(2) .GT. ymin) ) then !bottom
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
               trbuf_s(itr) = part
               call destroy_particle
               itr = itr + 1 
            elseif (part%xp(1) .LT. xmin) then !bottom right
               brbuf_s(ibr) = part
               call destroy_particle
               ibr = ibr + 1
            else   !right
               rbuf_s(ir) = part
               call destroy_particle
               ir = ir + 1
            end if
         elseif (part%xp(2) .LT. ymin) then
            if (part%xp(1) .GT. xmax) then !top left
               tlbuf_s(itl) = part
               call destroy_particle
               itl = itl + 1
            else if (part%xp(1) .LT. xmin) then !bottom left
               blbuf_s(ibl) = part
               call destroy_particle
               ibl = ibl + 1
            else  !left
               lbuf_s(il) = part
               call destroy_particle
               il = il + 1
            end if
         elseif ( (part%xp(1) .GT. xmax) .AND. &
                  (part%xp(2) .LT. ymax) .AND. &
                  (part%xp(2) .GT. ymin) ) then !top
            tbuf_s(itop) = part
            call destroy_particle
            itop = itop + 1
         elseif ( (part%xp(1) .LT. xmin) .AND. &
                  (part%xp(2) .LT. ymax) .AND. &
                  (part%xp(2) .GT. ymin) ) then !bottom
            bbuf_s(ib) = part
            call destroy_particle
            ib = ib + 1 
         else
         part => part%next
         end if 
         
      end do

      !Now everyone exchanges the counts with all neighbors:
      !Left/right:
      call MPI_Sendrecv(pr_s,1,mpi_integer,rproc,3, &
             pl_r,1,mpi_integer,lproc,3,mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pl_s,1,mpi_integer,lproc,4, &
             pr_r,1,mpi_integer,rproc,4,mpi_comm_world,istatus,ierr)

      !Top/bottom:
      call MPI_Sendrecv(pt_s,1,mpi_integer,tproc,5, &
             pb_r,1,mpi_integer,bproc,5,mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pb_s,1,mpi_integer,bproc,6, &
             pt_r,1,mpi_integer,tproc,6,mpi_comm_world,istatus,ierr)

      !Top right/bottom left:
      call MPI_Sendrecv(ptr_s,1,mpi_integer,trproc,7, &
             pbl_r,1,mpi_integer,blproc,7,mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pbl_s,1,mpi_integer,blproc,8, &
             ptr_r,1,mpi_integer,trproc,8,mpi_comm_world,istatus,ierr)

       !Top left/bottom right:
      call MPI_Sendrecv(ptl_s,1,mpi_integer,tlproc,9, &
             pbr_r,1,mpi_integer,brproc,9,mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pbr_s,1,mpi_integer,brproc,10, &
              ptl_r,1,mpi_integer,tlproc,10,mpi_comm_world,istatus,ierr)

      !Now everyone has the number of particles arriving from every neighbor
      !If the count is greater than zero, exchange:

      !Allocate room to receive from each side
      allocate(rbuf_r(pr_r),trbuf_r(ptr_r),tbuf_r(pt_r),tlbuf_r(ptl_r))
      allocate(lbuf_r(pl_r),blbuf_r(pbl_r),bbuf_r(pb_r),brbuf_r(pbr_r))
     
      !Send to right:
      if (pr_s .GT. 0) then
      call mpi_isend(rbuf_s,pr_s,particletype,rproc,11,mpi_comm_world,req(1),ierr)
      else
      req(1) = mpi_request_null
      end if

      !Receive from left:
      if (pl_r .GT. 0) then
      call mpi_irecv(lbuf_r,pl_r,particletype,lproc,11,mpi_comm_world,req(2),ierr)
      else
      req(2) = mpi_request_null
      end if

      !Send to left:
      if (pl_s .GT. 0) then
      call mpi_isend(lbuf_s,pl_s,particletype,lproc,12,mpi_comm_world,req(3),ierr)
      else
      req(3) = mpi_request_null
      end if

      !Receive from right:
      if (pr_r .GT. 0) then
      call mpi_irecv(rbuf_r,pr_r,particletype,rproc,12,mpi_comm_world,req(4),ierr)
      else
      req(4) = mpi_request_null
      end if

      !Send to top:
      if (pt_s .GT. 0) then
      call mpi_isend(tbuf_s,pt_s,particletype,tproc,13,mpi_comm_world,req(5),ierr)
      else
      req(5) = mpi_request_null
      end if
      
      !Receive from bottom:
      if (pb_r .GT. 0) then
      call mpi_irecv(bbuf_r,pb_r,particletype,bproc,13,mpi_comm_world,req(6),ierr)
      else
      req(6) = mpi_request_null
      end if

      !Send to bottom:
      if (pb_s .GT. 0) then
      call mpi_isend(bbuf_s,pb_s,particletype,bproc,14,mpi_comm_world,req(7),ierr)
      else
      req(7) = mpi_request_null
      end if
      
      !Recieve from top:
      if (pt_r .GT. 0) then
      call mpi_irecv(tbuf_r,pt_r,particletype,tproc,14,mpi_comm_world,req(8),ierr)
      else
      req(8) = mpi_request_null
      end if

      !Send to top right:
      if (ptr_s .GT. 0) then
      call mpi_isend(trbuf_s,ptr_s,particletype,trproc,15,mpi_comm_world,req(9),ierr)
      else
      req(9) = mpi_request_null
      end if
     
      !Receive from bottom left:
      if (pbl_r .GT. 0) then
      call mpi_irecv(blbuf_r,pbl_r,particletype,blproc,15,mpi_comm_world,req(10),ierr)
      else 
      req(10) = mpi_request_null
      end if
    
      !Send to bottom left:
      if (pbl_s .GT. 0) then
      call mpi_isend(blbuf_s,pbl_s,particletype,blproc,16,mpi_comm_world,req(11),ierr)
      else
      req(11) = mpi_request_null
      end if
     
      !Receive from top right:
      if (ptr_r .GT. 0) then
      call mpi_irecv(trbuf_r,ptr_r,particletype,trproc,16,mpi_comm_world,req(12),ierr)
      else 
      req(12) = mpi_request_null
      end if

      !Send to top left:
      if (ptl_s .GT. 0) then
      call mpi_isend(tlbuf_s,ptl_s,particletype,tlproc,17,mpi_comm_world,req(13),ierr)
      else 
      req(13) = mpi_request_null
      end if
    
      !Receive from bottom right:
      if (pbr_r .GT. 0) then
      call mpi_irecv(brbuf_r,pbr_r,particletype,brproc,17,mpi_comm_world,req(14),ierr)
      else 
      req(14) = mpi_request_null
      end if
  
      !Send to bottom right:
      if (pbr_s .GT. 0) then
      call mpi_isend(brbuf_s,pbr_s,particletype,brproc,18,mpi_comm_world,req(15),ierr)
      else
      req(15) = mpi_request_null
      end if
  
      !Receive from top left:
      if (ptl_r .GT. 0) then
      call mpi_irecv(tlbuf_r,ptl_r,particletype,tlproc,18,mpi_comm_world,req(16),ierr)
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

  subroutine set_bounds  
        use pars
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

  subroutine particle_init
      use pars
      use con_data
      implicit none
      include 'mpif.h' 
      integer :: values(8)
      integer :: idx,ierr

      !Create the seed for the random number generator:
      call date_and_time(VALUES=values)
      iseed = -(myid+values(8)+values(7)+values(6))


      numpart = tnumpart/numprocs
      if (myid == 0) then
      numpart = numpart + MOD(tnumpart,numprocs)
      endif


      !Initialize ngidx, the particle global index for this processor
      ngidx = 1

      !Initialize the linked list of particles:
      nullify(part,first_particle)
      

      do idx=1,numpart

         ! Call new_particle, which creates a new particle basd on some strategy dictated by inewpart
         call new_particle(idx,myid)


      ngidx = ngidx + 1
      end do


      partTsrc = 0.0
      partTsrc_t = 0.0
      partHsrc = 0.0
      partHsrc_t = 0.0
      partTEsrc = 0.0
      partTEsrc_t = 0.0


  end subroutine particle_init

  subroutine particle_setup

      use pars
      implicit none 
      include 'mpif.h'

      integer :: blcts(3),types(3)
      integer :: ierr
      real :: pi
      integer(kind=MPI_ADDRESS_KIND) :: extent,lb
      integer(kind=MPI_ADDRESS_KIND) :: extent2,lb2,displs(3)
      integer :: num_reals,num_integers,num_longs
      character*4 :: myid_char
      character*80 :: traj_file

      !First set up the neighbors for the interpolation stage:
      call assign_nbrs

      !Also assign the x,y,z max and mins to track particles leaving
      call set_bounds

      !Set up the path for the trajectory files
      if (itrajout) then
      write(myid_char,'(i4.4)') myid
      path_traj = trim(adjustl(path_seed))//"particle_traj/"//myid_char//".dat"
      ntraj = 128
      open(ntraj,file=path_traj,form='formatted',status='replace')
      end if


      !Lognormal distribution parameters  -- Must be called even on restart!
      mult_factor = 25
      pdf_factor = 4.8081e-04*real(mult_factor)
      pdf_prob = pdf_factor/(1 + pdf_factor)

      !Adjust the multiplicity so that the total number of particles isn't altered:
      !mult_c = mult_init/mult_factor
      !mult_a = mult_init/(1.0-pdf_prob)*(1 - pdf_prob/real(mult_factor))

      mult_a = mult_init*(1+4.8081e-04*mult_factor)/(1+4.8081e-04)

      mult_c = mult_init*(1+4.8081e-04*mult_factor)/(mult_factor*(1+4.8081e-04))      


      !set_binsdata does logarithmic binning!
      !Radius histogram
      call set_binsdata(bins_rad,histbins+2,1.0e-8,1.0e-3)

      !Residence time histogram
      call set_binsdata(bins_res,histbins+2,1.0e-1,1.0e4)

      !Activated time histogram
      call set_binsdata(bins_actres,histbins+2,1.0e-1,1.0e4)

      !Activated time until death histogram
      call set_binsdata(bins_acttodeath,histbins+2,1.0e-1,1.0e4)

      !Num activations histogram
      call set_binsdata_integer(bins_numact,histbins+2,0.0)


      !Initialize the linked list of particles:
      nullify(part,first_particle)

      !Set up MPI datatypes for sending particle information
      !MUST UPDATE IF THINGS ARE ADDED/REMOVED FROM PARTICLE STRUCTURE

      num_reals = 6*3+16
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
      call mpi_type_create_struct(3,blcts,displs,types,particletype,ierr)

       call mpi_type_get_true_extent(particletype,lb2,extent2,ierr)
       call mpi_type_get_extent(particletype,lb2,extent,ierr)
       if (extent .NE. sizeof(part) ) then
          if (myid==0) then
          write(*,*) 'WARNING: extent of particletype not equalto sizeof(part):'
          write(*,*) 'sizeof(part) = ', sizeof(part)
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

  subroutine save_particles
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
      call mpi_allgather(numpart,1,mpi_integer,pnum_vec,1,mpi_integer,mpi_comm_world,ierr)

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
      call mpi_file_open(mpi_comm_world, path_sav_part, &
                        mpi_mode_create+mpi_mode_rdwr, &
                        mpi_info_null,fh,ierr)

      zoffset = 0
      !Write tnumpart first:
      if (myid==0) then
      call mpi_file_write_at(fh,zoffset,tnumpart,1,mpi_integer,istatus,ierr)
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
      call mpi_file_write_at(fh,offset,writebuf,numpart,particletype,istatus,ierr)
      end if

      call mpi_file_close(fh,ierr)

      write(*,*) 'proc',myid,'wrote numpart = ',numpart

      if (myid==0) write(*,7000) path_sav_part
 7000 format(' PARTICLE DATA IS WRITTEN IN FILE  ',a80)

  end subroutine save_particles

  subroutine read_part_res
      use pars
      implicit none
      include 'mpif.h'

      integer :: istatus(mpi_status_size), ierr, fh
      integer(kind=mpi_offset_kind) :: zoffset,offset
      integer :: myp,totalp 
      integer :: iproc,i,pidxmax,numloop,partloop,readpart
      type(particle), allocatable :: readbuf(:)

      if (myid==0) write(*,7000) path_part 
 7000 format(' READING PARTICLE DATA FROM  ',a80)


      call mpi_file_open(mpi_comm_world,path_part,  &
                        mpi_mode_rdonly, &
                        mpi_info_null,fh,ierr)


      !Read in the total number of particles:
      offset = 0
      call mpi_file_read_at_all(fh,offset,tnumpart,1,mpi_integer,istatus,ierr)
      if (myid==0) write(*,*) 'read tnumpart = ',tnumpart
    
      offset = 4

      !For many particles (above ~10 million), I can't read them all
      !into the readbuf at the same time - must break up into chunks.
      !Arbitrarily choose 5 million particles at a time (~840 MB)

      !numloop will be 1 if tnumpart < 5 million, increasing from there
      numloop = floor(tnumpart/5e6)+1

      do partloop = 1,numloop

      if (partloop == numloop) then
         readpart = tnumpart - (numloop-1)*5e6
      else
         readpart = 5e6
      end if

      allocate(readbuf(readpart))

      call mpi_file_read_at_all(fh,offset,readbuf,readpart,particletype,istatus,ierr)

      do i = 1,readpart
        !Now - does it lie within this proc's bounds?
        if (readbuf(i)%xp(2) .GT. ymin .AND. &
            readbuf(i)%xp(2) .LT. ymax .AND. &
            readbuf(i)%xp(1) .GT. xmin .AND. &
            readbuf(i)%xp(1) .LT. xmax) then 
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

      offset = offset + sizeof(part)*readpart

      end do

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
     
      call mpi_allreduce(myp,totalp,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)

      write(*,*) 'proc',myid,'read in numpart:',myp
      if (myid==0) write(*,*) 'total number of particles read:',totalp

  end subroutine read_part_res

  subroutine particle_reintro(it)
      use pars
      use con_data
      use con_stats
      use fields
      implicit none
      include 'mpif.h'

      integer :: it,it_delay
      integer :: ierr,randproc,np,my_reintro
      real :: totdrops,t_reint


      if (inewpart.eq.4) then

      !!Sea spray, given by Andreas SSGF 98
      call andreas_dist_num(totdrops) ! drops per m^2 per s based on Andreas 98      

      t_reint = 1800.0  !Time after which to start injecting
      it_delay = 200    !Num of time steps between injection events (sometimes too few are produced and if it's < numprocs then it gets rounded to zero)
      
      
      if (time .gt. t_reint) then
         my_reintro = xl*yl*dt*real(it_delay)*totdrops/(numprocs*mult_init)
         tot_reintro = 0
      else
         my_reintro = 0
         tot_reintro = 0
      endif
      
      elseif (inewpart.eq.2) then

      !!Pi Chamber, given by constant injection rate (nprime)
      it_delay = 20

      my_reintro = nprime*(1./60.)*(10.**6.)*dt*4/numprocs*real(it_delay) !4m^3 (vol chamber)
      tot_reintro = 0

      else

      my_reintro = 0
      tot_reintro = 0

      end if

      if (mod(it, it_delay)==0) then

      tot_reintro = my_reintro*numprocs


      if (myid==0) write(*,*) 'time,tot_reintro:',time,tot_reintro

      do np=1,my_reintro

         call new_particle(np,myid)

         !Update this processor's global ID for each one created:
         ngidx = ngidx + 1

      end do
      end if

      !Now update the total number of particles
      numpart = 0
      part => first_particle
      do while (associated(part))
      numpart = numpart + 1
      part => part%next
      end do

      call mpi_allreduce(numpart,tnumpart,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)


  end subroutine particle_reintro

  subroutine create_particle(xp,vp,Tp,m_s,kappa_s,mult,rad_init,idx,procidx)
      use pars
      implicit none

      real :: xp(3),vp(3),Tp,qinfp,rad_init,pi,m_s,kappa_s
      integer :: idx,procidx
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
      part%procidx = procidx
      part%nbr_pidx = -1
      part%nbr_procidx = -1
      part%mult = mult
      part%res = 0.0
      part%actres = 0.0
      part%m_s = m_s
      part%kappa_s = kappa_s
      part%dist = 0.0
      part%u_sub(1:3) = 0.0
      part%sigm_s = 0.0
      part%numact = 0.0

      
  end subroutine create_particle

  subroutine new_particle(idx,procidx)
  use pars
  use con_data
  implicit none

  real :: xv,yv,zv,ran2,m_s
  real :: kappas_dinit,radius_dinit
  real :: xp_init(3)
  integer :: idx,procidx

  !C-FOG and FATIMA parameters: lognormal of accumulation + lognormal of coarse, with extra "resolution" on the coarse mode
  real :: S,M,kappa_s,rad_init
  real :: num_a,num_c,totnum_a,totnum_c
  integer*8 :: mult

  if (inewpart.eq.1) then  !Simple: properties as in params.in, randomly located in domain

      xv = ran2(iseed)*(xmax-xmin) + xmin
      yv = ran2(iseed)*(ymax-ymin) + ymin
      zv = ran2(iseed)*zl
      xp_init = (/xv,yv,zv/) 

      m_s = radius_init**3*pi2*2.0/3.0*rhow*Sal  !Using the salinity specified in params.in

      call create_particle(xp_init,vp_init,Tp_init,m_s,kappas_init,mult_init,radius_init,ngidx,procidx)




   elseif (inewpart.eq.2) then  !Same as above, but with NORMAL distribution of radius and kappa given by radius_std and kappas_std
                              !Use for Pi Chamber

      xv = ran2(iseed)*(xmax-xmin) + xmin
      yv = ran2(iseed)*(ymax-ymin) + ymin
      !zv = ran2(iseed)*(zmax-zmin) + zmin
      zv = zl/2.0   !Midplane of the Pi Chamber domain
      xp_init = (/xv,yv,zv/) 

      ! Set distribution for initial radius
      radius_dinit = abs(radius_std*sqrt(-2*log(ran2(iseed)))*cos(pi2*ran2(iseed)) + radius_init)

      ! Set distribution for kappa_s
      kappas_dinit = abs(kappas_std * sqrt(-2*log(ran2(iseed)))*cos(pi2*ran2(iseed)) + kappas_init)

      m_s = radius_dinit**3*pi2*2.0/3.0*rhow*Sal  !Using the salinity specified in params.in

      call create_particle(xp_init,vp_init,Tp_init,m_s,kappas_dinit,mult_init,radius_dinit,ngidx,procidx)




    elseif (inewpart.eq.3) then !Special for C-FOG: the scheme used in Richter et al., BLM, 2021

      xv = ran2(iseed)*(xmax-xmin) + xmin
      yv = ran2(iseed)*(ymax-ymin) + ymin
      zv = ran2(iseed)*(zi-zw1) + zw1
      xp_init = (/xv,yv,zv/) 


      !Generate
      if (ran2(iseed) .gt. pdf_prob) then   !It's accumulation mode
         S = 0.5
         M = -1.95
         kappa_s = 0.6
         mult = mult_a

         !With these parameters, get m_s and rad_init from distribution
         call lognormal_dist(rad_init,m_s,kappa_s,M,S)
         num_a = num_a + 1

      else  !It's coarse mode

         S = 0.45
         M = 0.0
         kappa_s = 1.2
         mult = mult_c

         !With these parameters, get m_s and rad_init from distribution
         call lognormal_dist(rad_init,m_s,kappa_s,M,S)
         num_c = num_c + 1

      end if      

      !Force the output particle to be a coarse mode particle
      if (idx==1 .and. procidx==0) then
         !Force particle log output to be a coarse mode in fog layer -- "giant mode"
         S = 0.45
         M = 0.0
         kappa_s = 1.2
         mult = mult_c

         !With these parameters, get m_s and rad_init from distribution
         call lognormal_dist(rad_init,m_s,kappa_s,M,S)
         xp_init(3) = 10.0
      end if


      call create_particle(xp_init,vp_init,Tp_init,m_s,kappa_s,mult,rad_init,idx,procidx)



   elseif (inewpart.eq.4) then  !Sea spray, from Andreas 98

      
         m_s = rad_init**3*pi2*2.0/3.0*rhow*Sal  !Using the salinity specified in params.in
         
         Tp_init = tsfcc(1)

         xp_init(1) = ran2(iseed)*(xmax-xmin) + xmin
         xp_init(2) = ran2(iseed)*(ymax-ymin) + ymin
         xp_init(3) = ran2(iseed)*8.0  !Distributing between 0 and 8 meters (like a sig. wave height)

         call andreas_dist(rad_init)

         m_s = rad_init**3*pi2*2.0/3.0*rhow*Sal  !Using the salinity specified in params.in

         vp_init(3) = ran2(iseed)*4.0

         call create_particle(xp_init,vp_init,Tp_init,m_s,kappas_init,mult_init,rad_init,ngidx,procidx) 

   elseif (inewpart.eq.5) then !Special for Sc: still working on it
      
      xv = ran2(iseed)*(xmax-xmin) + xmin
      yv = ran2(iseed)*(ymax-ymin) + ymin
      zv = ran2(iseed)*zl
      xp_init = (/xv,yv,zv/)
      
      
      !Generate
      S = 0.45
      M = 0.0
      kappa_s = 1.2
      mult = mult_init
      
      !With these parameters, get m_s and rad_init from distribution
      call lognormal_dist(rad_init,m_s,kappa_s,M,S)
       
      call create_particle(xp_init,vp_init,Tp_init,m_s,kappa_s,mult,rad_init,idx,procidx)
      
   elseif (inewpart.eq.6) then !FATIMA IOP-5
      
      xv = ran2(iseed)*(xmax-xmin) + xmin
      yv = ran2(iseed)*(ymax-ymin) + ymin
      zv = ran2(iseed)*(zi-zw1) + zw1
      xp_init = (/xv,yv,zv/)


      !Generate
      if (ran2(iseed) .gt. pdf_prob) then   !It's accumulation mode
         S = 0.2403
         M = -1.7570
         kappa_s = 0.6
         mult = mult_a

         !With these parameters, get m_s and rad_init from distribution
         call lognormal_dist(rad_init,m_s,kappa_s,M,S)
         num_a = num_a + 1

      else  !It's coarse mode

         S = 0.2997
         M = -0.1930
         kappa_s = 1.2
         mult = mult_c

         !With these parameters, get m_s and rad_init from distribution
         call lognormal_dist(rad_init,m_s,kappa_s,M,S)
         num_c = num_c + 1

      end if

      !Force the output particle to be a coarse mode particle
      if (idx==1 .and. procidx==0) then
         !Force particle log output to be a coarse mode in fog layer -- "giant mode"
         S = 0.2997
         M = -0.1930
         kappa_s = 1.2
         mult = mult_c

         !With these parameters, get m_s and rad_init from distribution
         call lognormal_dist(rad_init,m_s,kappa_s,M,S)
         xp_init(3) = 10.0
      end if


      call create_particle(xp_init,vp_init,Tp_init,m_s,kappa_s,mult,rad_init,idx,procidx)

   end if

  end subroutine new_particle

  subroutine lognormal_dist(rad_init,m_s,kappa_s,M,S)
  use pars
  use con_data
  implicit none
  include 'mpif.h'

  real, intent(inout) :: rad_init,m_s,kappa_s
  real :: ran2,cdf_func_single
  real :: M,S
  real :: d1,d2,err,dhalf,ftest,CDF
  real :: daerosol
  real :: a(4), rtr(3), rti(3)
  integer :: iter, k

  !Use bisection to get the radius based on the CDF contained in
  !function cdf_func
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
        d1 = 0.1
        write(*,'(a12)') 'CDF convergence error'
        exit
     end if


  end do

  daerosol = d1*1.0e-6  !Don't forget to convert micron to m
  m_s = 2.0/3.0*pi2*(daerosol/2.0)**3*rhos

  !Now have the dry aerosol diameter and mass, must rehydrate it using Kohler
  !theory to get the proper initial condition
  !Rehydrate to 95% RH

  !The equation for the equilibrium radius is a cubic:
   a(4) = log(0.95)  !RH = 95%
   a(3) = -2.0*Mw*Gam/Ru/rhow/Tp_init
   a(2) = 0.0
   a(1) = kappa_s*m_s / (2.0/3.0*pi2*rhos)
   ! calculate all roots of the cubic equation (real and complex)
   !The method is to construct an upper Hessenberg matrix
   !whose eigenvalues are the desired roots, and then use the
   !routines balanc and hqr . The real and imaginary parts of the
   !roots are returned in rtr(1:m) and rti(1:m) , respectively.
   call eigen_roots(a,3,rtr,rti)
   !now select the real and positive one
   do k=1,3
     if(rti(k).eq.0.) then
        if(rtr(k).gt.0.)then
           rad_init = rtr(k)
        endif
     endif
   enddo

  end subroutine lognormal_dist

  subroutine andreas_dist(rad_init)
  use pars
  use fields
  use con_data
  implicit none
  include 'mpif.h'

  real, intent(inout) :: rad_init
  real :: ran2,cdf_func,prob
  real :: M_a,S_a,M_c,S_c,totarea
  real :: daerosol,totdrops
  real :: c1,c2,c3,u14,cdn10,a
  real :: a1,a2,r80,r0,dh,binsdata(100)
  real :: dFssum(101)
  real :: dFsdr80(100),dFmsdr0(100),dr80_dr0
  integer :: iter,i,nbin,num_create,j

  ! implementing the Andreas 1998 sea spray generation function
  !Set the parameters of the two lognormals:
  real :: rmin,rmax,rmin10,rmax10

    rmin10 = log10(2e-06)
    rmax10 = log10(500e-06)
    ! Testing to see if u10 is the problem
    !u10 = 15.0
!   Calculate size of interval
    dh = 0.0242

!   ===== update x-axis for each bin =====
    binsdata(1) = rmin10
    do i = 1,99
      binsdata(i+1)= dh+binsdata(i)
    end do

    binsdata = 1e6*(10**binsdata)

    ! Now figure out how many particles are produced in each of
    ! those bins. This is a function of u10

    c1 = 10*smithssgf(u10,10.0,0.4)
    c2 = c1*(37.5**1.8)
    c3 = c2*(100**5.2)
    dFssum(1) = 0
    do i = 1,100
       r0 = binsdata(i)
       !calculate radius at 80%RH
       r80 = 0.518*r0**0.976;

      !now apply these to Eqs. 3.5 in Andreas 1998
       if (r80 .lt. 10) then
         dFsdr80(i) =  smithssgf(u10,r80,0.4)
       elseif ((r80 .ge. 10) .and. (r80 .lt. 37.5)) then
         dFsdr80(i) = c1/r80
       elseif ((r80 .ge. 37.5) .and. (r80 .lt. 100)) then
         dFSdr80(i) = c2*r80**-2.8
       elseif ((r80 .ge. 100) .and. (r80 .lt. 250)) then
         dFsdr80(i) = c3*r80**-8
       endif

      !apply eq. 3.8 from Andreas 1998
      dFmsdr0(i) = 3.5*dFsdr80(i)*0.506*r0**-0.024

      if (i==1) then
        dFssum(i+1) = dFmsdr0(i)
      elseif (i .gt. 1) then
        dFssum(i+1) = dFssum(i) + dFmsdr0(i)
      endif

    end do


    totdrops = dFssum(101)
    a = totdrops*ran2(iseed)

    do i = 1,100
       if ((a .gt. dFssum(i)) .and. (a .le. dFssum(i+1))) then
          rad_init = binsdata(i) + ran2(iseed)*(binsdata(i+1)-binsdata(i))
       endif
    enddo

    rad_init = 1e-6*rad_init

  end subroutine andreas_dist

  subroutine andreas_dist_num(totdrops)
  use pars
  use fields
  use con_data
  implicit none
  include 'mpif.h'

  real, intent(inout) :: totdrops
  real :: ran2,cdf_func,prob
  real :: M_a,S_a,M_c,S_c,totarea
  real :: daerosol
  real :: c1,c2,c3,u14,cdn10,a
  real :: a1,a2,r80,r0,dh,binsdata(100)
  real :: dFssum(101)
  real :: dFsdr80(100),dFmsdr0(100),dr80_dr0
  integer :: iter,i,nbin,num_create,j

  !implementing the Andreas 1998 sea spray generation function
  !Set the parameters of the two lognormals:
  real :: rmin,rmax,rmin10,rmax10

    rmin10 = log10(2e-06)
    rmax10 = log10(500e-06)
    ! Testing to see if u10 is the problem
    !u10 = 15.0
!   Calculate size of interval
    dh = 0.0242

!   ===== update x-axis for each bin =====
    binsdata(1) = rmin10
    do i = 1,99
      binsdata(i+1)= dh+binsdata(i)
    end do

    binsdata = 1e6*(10**binsdata)

    ! Now figure out how many particles are produced in each of
    ! those bins. This is a function of u10

    c1 = 10*smithssgf(u10,10.0,0.4)
    c2 = c1*(37.5**1.8)
    c3 = c2*(100**5.2)
    dFssum(1) = 0
    do i = 1,100
     r0 = binsdata(i)
     ! calculate radius at 80%RH
     r80 = 0.518*r0**0.976;

  ! now apply these to Eqs. 3.5 in Andreas 1998
     if (r80 .lt. 10) then
       dFsdr80(i) =  smithssgf(u10,r80,0.4)
     elseif ((r80 .ge. 10) .and. (r80 .lt. 37.5)) then
       dFsdr80(i) = c1/r80
     elseif ((r80 .ge. 37.5) .and. (r80 .lt. 100)) then
       dFSdr80(i) = c2*r80**-2.8
     elseif ((r80 .ge. 100) .and. (r80 .lt. 250)) then
       dFsdr80(i) = c3*r80**-8
     endif

    ! apply eq. 3.8 from Andreas 1998
    dFmsdr0(i) = 3.5*dFsdr80(i)*0.506*r0**-0.024

    if (i==1) then
      dFssum(i+1) = dFmsdr0(i)
    elseif (i .gt. 1) then
      dFssum(i+1) = dFssum(i) + dFmsdr0(i)
    endif

  end do


  totdrops = dFssum(101)

  end subroutine andreas_dist_num


  subroutine particle_bcs_nonperiodic
  use con_stats
  use pars
  implicit none
  real :: top,bot
  integer :: idx,procidx,idx_old,procidx_old

  real :: xv,yv,zv,ran2,m_s
  real :: kappas_dinit,radius_dinit
  real :: xp_init(3)

  !C-FOG and FATIMA parameters: lognormal of accumulation + lognormal of coarse, with extra "resolution" on the coarse mode
  real :: S,M,kappa_s,rad_init
  real :: num_a,num_c,totnum_a,totnum_c
  integer*8 :: mult


  !Assumes domain goes from [0,xl),[0,yl),[0,zl]
  !Also maintain the number of particles on each proc

  part => first_particle
  do while (associated(part))

    !perfectly elastic collisions on top, bottom walls
    !i.e. location is reflected, w-velocity is negated

    top = z(nnz)
    !bot = 0.0 + part%radius
    !bot = zw1
    bot = 0.0

    if (part%xp(3) .GT. top) then
       part%xp(3) = top - (part%xp(3)-top)
       part%vp(3) = -part%vp(3)
       part => part%next
    elseif (part%xp(3) .LT. bot) then
       !part%xp(3) = bot + (bot-part%xp(3))
       !part%vp(3) = -part%vp(3)
       !part => part%next

       idx_old = part%pidx
       procidx_old = part%procidx

       !Before destroying it, put its residence time in histogram
       call add_histogram(bins_res,hist_res,histbins+2,part%res,part%mult)

       !Also record this in the "activation till death" residence time
       if (part%radius .gt. part%rc) then
         call add_histogram(bins_acttodeath,hist_acttodeath,histbins+2,part%actres,part%mult)
       end if

       !Also record the number of activations
       call add_histogram_integer(bins_numact,hist_numact,histbins+2,part%numact)

       !Also record the size of the dead droplet
       call add_histogram(bins_rad,hist_raddeath,histbins+2,part%radius,part%mult)

       if (ireintro.eq.1 .and. inewpart.eq.6) then  !FATIMA can reintroduce particle of the same type

            xv = ran2(iseed)*(xmax-xmin) + xmin
            yv = ran2(iseed)*(ymax-ymin) + ymin
            zv = ran2(iseed)*(zi-zw1) + zw1
            xp_init = (/xv,yv,zv/)

            ! acummulation mode
            if (part%mult .eq. mult_a) then

                S = 0.2403
                M = -1.7570
                kappa_s = 0.6
                mult = mult_a

            ! coarse mode
            elseif (part%mult .eq. mult_c) then

                S = 0.2997
                M = -0.1930
                kappa_s = 1.2
                mult = mult_c

            end if

            ! destroy old particle before creating new one
            call destroy_particle
            num_destroy = num_destroy + 1

            !With these parameters, get m_s and rad_init from distribution
            call lognormal_dist(rad_init,m_s,kappa_s,M,S)

            call create_particle(xp_init,vp_init,Tp_init,m_s,kappa_s,mult,rad_init,idx_old,procidx_old)

       else

            call destroy_particle
            num_destroy = num_destroy + 1

       end if
       
       if (icase.eq.5 .or. icase.eq.3) then
          call new_particle(idx_old,procidx_old)
       end if


    else
       part => part%next
    end if

  end do

  end subroutine particle_bcs_nonperiodic

  subroutine particle_bcs_periodic
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

  subroutine particle_update_rk3(it,istage)
      use pars
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
      ufsum_t = 0.0
      upwp_t = 0.0
      vpsqrsum_t = 0.0
      ufsqrsum_t = 0.0
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
             part%Tf = part%Tf*(psurf/(psurf-part%xp(3)*rhoa*grav))**(-Rd/Cpa)
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
         Eff_S = part%kappa_s*part%m_s*rhow/rhos/(Volp*rhop-part%m_s)
         estar = einf*exp(Mw*Lv/Ru*(1.0/part%Tf-1.0/part%Tp)+Eff_C-Eff_S)
         part%qstar = Mw/Ru*estar/part%Tp/rhoa

  
         xtmp(1:3) = part%xp(1:3) + dt*zetas(istage)*part%xrhs(1:3)
         vtmp(1:3) = part%vp(1:3) + dt*zetas(istage)*part%vrhs(1:3) 
         Tptmp = part%Tp + dt*zetas(istage)*part%Tprhs_s
         Tptmp = Tptmp + dt*zetas(istage)*part%Tprhs_L
         radiustmp = part%radius + dt*zetas(istage)*part%radrhs

         part%xrhs(1:3) = part%vp(1:3)
         part%vrhs(1:3) = corrfac*taup_i*(part%uf(1:3)-part%vp(1:3)) + part_grav(1:3)

         if (ievap .EQ. 1) then      
            part%radrhs = Shp/9.0/Sc*rhop/rhow*part%radius*taup_i*(part%qinf-part%qstar) !assumes qinf=rhov/rhoa rather than rhov/rhom
         else
            part%radrhs = 0.0
         end if


         part%Tprhs_s = -Nup/3.0/Pra*CpaCpp*rhop/rhow*taup_i*(part%Tp-part%Tf)
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


      !Now that particles are in their updated position, 
      !compute their contribution to the momentum coupling:
      !t_s = mpi_wtime()
      call particle_coupling_update
      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time cpl: ', t_f - t_s

      call particle_coupling_exchange

      call particle_stats

 
      !Finally, now that coupling and statistics arrays are filled, 
      !Transpose them back to align with the velocities:
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,1), &
                     partsrc(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,2), &
                     partsrc(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,3), &
                     partsrc(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(partTsrc_t(0:nnz+1,iys:iye,mxs:mxe), &
                     partTsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(partHsrc_t(0:nnz+1,iys:iye,mxs:mxe), &
                     partHsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(partTEsrc_t(0:nnz+1,iys:iye,mxs:mxe), &
                     partTEsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(mwsum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     mwsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      !Try only calling these when the history data is being written:
      if(mtrans  .and. istage .eq. 3) then
      call ztox_trans(upwp_t(0:nnz+1,iys:iye,mxs:mxe), &
                     upwp(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,1), &
                     vpsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,2), &
                     vpsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,3), &
                     vpsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,1), &
                     vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,2), &
                     vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,3), &
                     vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(ufsum_t(0:nnz+1,iys:iye,mxs:mxe,1), &
                     ufsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(ufsum_t(0:nnz+1,iys:iye,mxs:mxe,2), &
                     ufsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(ufsum_t(0:nnz+1,iys:iye,mxs:mxe,3), &
                     ufsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(ufsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,1), &
                     ufsqrsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(ufsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,2), &
                     ufsqrsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(ufsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,3), &
                     ufsqrsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(Tpsum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     Tpsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(Tpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     Tpsqrsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(Tfsum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     Tfsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(qfsum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     qfsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(wpTpsum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     wpTpsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(partcount_t(0:nnz+1,iys:iye,mxs:mxe), &
                     partcount(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(radsum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     radsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(rad2sum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     rad2sum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs) 

      call ztox_trans(multcount_t(0:nnz+1,iys:iye,mxs:mxe), &
                     multcount(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)


      call ztox_trans(qstarsum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     qstarsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

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
      call mpi_allreduce(numpart,tnumpart,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
      !Compute average particle Reynolds number
      call mpi_allreduce(myRep_avg,Rep_avg,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

      Rep_avg = Rep_avg/tnumpart

      call mpi_allreduce(mylwc_sum,lwc,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

      call mpi_allreduce(myphiw_sum,phiw,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

      call mpi_allreduce(myphiv_sum,phiv,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

      
      phiw = phiw/xl/yl/zl/rhoa
      phiv = phiv/xl/yl/zl

      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time mpi_allreduce: ', t_f - t_s

  end subroutine particle_update_rk3

  subroutine particle_update_BE(it)
      use pars
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'

      integer :: ierr,it,fluxloc,fluxloci
      real :: tmpbuf(9),tmpbuf_rec(9)
      integer :: intbuf(10),intbuf_rec(10)
      integer :: numdrop_center,tnumdrop_center
      real :: myradavg,myradmax,myradmin,mytempmax,mytempmin
      real :: myradavg_center,myradmsqr_center
      real :: myradmsqr
      real :: myqmin,myqmax
      real :: denom,dtl,sigma
      integer :: ix,iy,iz,im,flag,mflag
      real :: Rep,diff(3),diffnorm,corrfac,myRep_avg
      real :: Nup,Shp,rhop,taup_i,estar,einf
      real :: mylwc_sum,myphiw_sum,myphiv_sum,Volp
      real :: Eff_C,Eff_S
      real :: t_s,t_f,t_s1,t_f1
      real :: rt_start(2)
      real :: rt_zeroes(2)
      real :: taup0, dt_taup0, temp_r, temp_t, guess
      real :: tmp_coeff
      real :: xp3i
      real :: mod_Magnus
      real :: rad_i,Tp_i,vp_i(3),mp_i,rhop_i



      !First fill extended velocity field for interpolation
      call fill_ext

      partcount_t = 0.0
      vpsum_t = 0.0
      ufsum_t = 0.0
      upwp_t = 0.0
      vpsqrsum_t = 0.0
      ufsqrsum_t = 0.0
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

      !partsrc_t = 0.0
      !partTsrc_t = 0.0
      !partHsrc_t = 0.0
      !partTEsrc_t = 0.0

      pflux = 0.0
      pmassflux = 0.0
      penegflux = 0.0

      denum = 0
      actnum = 0
      num100 = 0
      num1000 = 0
      numimpos = 0
      num_destroy = 0

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
             part%Tf = part%Tf*(psurf/(psurf-part%xp(3)*rhoa*grav))**(-Rd/Cpa)
         end if

        diff(1:3) = part%vp - part%uf
        diffnorm = sqrt(diff(1)**2 + diff(2)**2 + diff(3)**2)
        Volp = pi2*2.0/3.0*part%radius**3
        rhop = (part%m_s+Volp*rhow)/Volp
        taup_i = 18.0*rhoa*nuf/rhop/(2.0*part%radius)**2
        Rep = 2.0*part%radius*diffnorm/nuf
        corrfac = (1.0 + 0.15*Rep**(0.687))

        corrfac = 1.0

        xp3i = part%xp(3)   !Store this to do flux calculation

        !Store these to compute the feedback terms
        rad_i = part%radius
        Tp_i = part%Tp
        vp_i(1:3) = part%vp(1:3)
        mp_i = Volp*rhop
        rhop_i = rhop

        !implicitly calculates next velocity and position
        part%xp(1:3) = part%xp(1:3) + dt*part%vp(1:3)
        part%vp(1:3) = (part%vp(1:3)+taup_i*dt*corrfac*part%uf(1:3)+dt*part_grav(1:3))/(1+dt*corrfac*taup_i)


        ! non-dimensionalizes particle radius and temperature before
        ! iteratively solving for next radius and temperature

        taup0 = (((part%m_s)/((2./3.)*pi2*radius_init**3) + rhow)*(radius_init*2)**2)/(18*rhoa*nuf)

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

               call gauss_newton_2d(part%vp,dt_taup0,rt_start, rt_zeroes,flag)

               if (flag==1) then
               num100 = num100+1

               call LV_solver(part%vp,dt_taup0,rt_start, rt_zeroes,flag)

               end if

               if (flag == 1) num1000 = num1000 + 1

               if      (isnan(rt_zeroes(1)) &
                  .OR. (rt_zeroes(1)*part%radius<0) &
                  .OR. isnan(rt_zeroes(2)) &
                  .OR. (rt_zeroes(2)<0) &
                  .OR. (rt_zeroes(1)*part%radius>1.0e-2)) & !These last 2 are very specific to pi chamber
                  !.OR. (rt_zeroes(2)*part%Tp > Tbot(1)*1.1)  &
                  !.OR. (rt_zeroes(2)*part%Tp < Ttop(1)*0.9)) &
               then

                write(*,'(a30,14e15.6)') 'WARNING: CONVERGENCE',  &
               part%radius,part%qinf,part%Tp,part%Tf,part%xp(3), &
               part%kappa_s,part%m_s,part%vp(1),part%vp(2),part%vp(3), &
               part%res,part%sigm_s,rt_zeroes(1),rt_zeroes(2)

                numimpos = numimpos + 1  !How many have failed?
                !If they failed (should be very small number), radius,
                !temp remain unchanged
                rt_zeroes(1) = 1.0
                rt_zeroes(2) = part%Tf/part%Tp


               end if

               !Get the critical radius based on old temp
               part%rc = crit_radius(part%m_s,part%kappa_s,part%Tp) 

               !Count if activated/deactivated
               if (part%radius > part%rc .AND. part%radius*rt_zeroes(1) < part%rc) then
                   denum = denum + 1

                   !Also add activated lifetime to histogram
               call add_histogram(bins_actres,hist_actres,histbins+2,part%actres,part%mult)
                   

               elseif (part%radius < part%rc .AND. part%radius*rt_zeroes(1) > part%rc) then
                   actnum = actnum + 1
                   part%numact = part%numact + 1.0

                   !Reset the activation lifetime
                   part%actres = 0.0

               endif

               !Redimensionalize
               part%radius = rt_zeroes(1)*part%radius
               part%Tp = rt_zeroes(2)*part%Tp
        end if

         if (part%radius .gt. 1.0e-2) then
         write(*,'(a30,12e15.6)') 'WARNING: BIG DROPLET',  &
         part%radius,part%qinf,part%Tp,part%Tf,part%xp(3), &
         part%kappa_s,part%m_s,part%vp(1),part%vp(2),part%vp(3), &
         part%res,part%sigm_s
         end if

         if (part%qinf .lt. 0.0) then
         write(*,'(a30,2i,12e15.6)') 'WARNING: NEG QINF',  &
         part%pidx,part%procidx, &
         part%radius,part%qinf,part%Tp,part%Tf,part%xp(3), &
         part%kappa_s,part%m_s,part%vp(1),part%vp(2),part%vp(3), &
         part%res,part%sigm_s
         end if


         !New volume and particle density
         Volp = pi2*2.0/3.0*part%radius**3
         rhop = (part%m_s+Volp*rhow)/Volp

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
         Eff_S = part%kappa_s*part%m_s*rhow/rhos/(Volp*rhop-part%m_s)
         estar = einf*exp(Mw*Lv/Ru*(1.0/part%Tf-1.0/part%Tp)+Eff_C-Eff_S)
         part%qstar = Mw/Ru*estar/part%Tp/rhoa

        if (ievap .EQ. 1) then
            !part%radrhs = Shp/9.0/Sc*rhop/rhow*part%radius*taup_i*(part%qinf-part%qstar) !assumes qinf=rhov/rhoa rather than rhov/rhom
            part%radrhs = (part%radius-rad_i)/dt
        else

            part%radrhs = 0.0

            !Also update the temperature directly using BE:
            tmp_coeff = Nup/3.0/Pra*CpaCpp*rhop/rhow*taup_i
            part%Tp = (part%Tp + tmp_coeff*dt*part%Tf)/(1+dt*tmp_coeff)
        end if

        !part%Tprhs_s = -Nup/3.0/Pra*CpaCpp*rhop/rhow*taup_i*(part%Tp-part%Tf)
        part%Tprhs_L = 3.0*Lv/Cpp/part%radius*part%radrhs
        part%Tprhs_s = (part%Tp-Tp_i)/dt - part%Tprhs_L

        part%xrhs(1:3) = part%vp(1:3)
        !part%vrhs(1:3) = corrfac*taup_i*(part%uf(1:3)-part%vp(1:3)) + part_grav(1:3)
        part%vrhs(1:3) = (part%vp(1:3)-vp_i(1:3))/dt

        part%res = part%res + dt
        part%actres = part%actres + dt


        !Store the particle flux now that everything has been updated
        if (part%xp(3) .gt. zl) then   !This will get treated in particle_bcs_nonperiodic, but record here
           fluxloc = nnz+1
           fluxloci = minloc(z,1,mask=(z.gt.xp3i))-1
        elseif (part%xp(3) .lt. 0.0) then !This will get treated in particle_bcs_nonperiodic, but record here
           fluxloci = minloc(z,1,mask=(z.gt.xp3i))-1
           fluxloc = 0
        else

        fluxloc = minloc(z,1,mask=(z.gt.part%xp(3)))-1
        fluxloci = minloc(z,1,mask=(z.gt.xp3i))-1

        end if  !Only apply flux calc to particles in domain

        if (xp3i .lt. part%xp(3)) then !Particle moved up

        do iz=fluxloci,fluxloc-1
           pflux(iz) = pflux(iz) + part%mult
           pmassflux(iz) = pmassflux(iz) + rhop*Volp*part%mult
           penegflux(iz) = penegflux(iz) + rhop*Volp*Cpp*part%Tp*part%mult
        end do

        elseif (xp3i .gt. part%xp(3)) then !Particle moved down

        do iz=fluxloc,fluxloci-1
           pflux(iz) = pflux(iz) - part%mult
           pmassflux(iz) = pmassflux(iz) - rhop*Volp*part%mult
           penegflux(iz) = penegflux(iz) - rhop*Volp*Cpp*part%Tp*part%mult
        end do

        end if  !Up/down conditional statement


      part => part%next
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

      call particle_stats

      !Finally, now that coupling and statistics arrays are filled, 
      !Transpose them back to align with the velocities:
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,1), &
                     partsrc(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,2), &
                     partsrc(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,3), &
                     partsrc(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(partTsrc_t(0:nnz+1,iys:iye,mxs:mxe), &
                     partTsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(partHsrc_t(0:nnz+1,iys:iye,mxs:mxe), &
                     partHsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(partTEsrc_t(0:nnz+1,iys:iye,mxs:mxe), &
                     partTEsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(mwsum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     mwsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(partcount_t(0:nnz+1,iys:iye,mxs:mxe), &
                     partcount(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(multcount_t(0:nnz+1,iys:iye,mxs:mxe), &
                     multcount(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(radsum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     radsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)



      !Try only calling these when the history data is being written:
      if(mtrans) then
      call ztox_trans(upwp_t(0:nnz+1,iys:iye,mxs:mxe), &
                     upwp(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,1), &
                     vpsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,2), &
                     vpsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,3), &
                     vpsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,1), &
                     vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs) 
      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,2), &
                     vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,3), &
                     vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(ufsum_t(0:nnz+1,iys:iye,mxs:mxe,1), &
                     ufsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(ufsum_t(0:nnz+1,iys:iye,mxs:mxe,2), &
                     ufsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(ufsum_t(0:nnz+1,iys:iye,mxs:mxe,3), &
                     ufsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(ufsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,1), &
                     ufsqrsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs) 
      call ztox_trans(ufsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,2), &
                     ufsqrsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(ufsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,3), &
                     ufsqrsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(Tpsum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     Tpsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(Tpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     Tpsqrsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(Tfsum_t(0:nnz+1,iys:iye,mxs:mxe), & 
                     Tfsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(qfsum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     qfsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)
      call ztox_trans(wpTpsum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     wpTpsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(rad2sum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     rad2sum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      call ztox_trans(qstarsum_t(0:nnz+1,iys:iye,mxs:mxe), &
                     qstarsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs, &
                     mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid, &
                     ncpu_s,numprocs)

      end if

      !Get particle count:
      numpart = 0
      numdrop = 0
      numdrop_center = 0
      numaerosol = 0

      myradavg = 0.0
      myradmsqr = 0.0
      myradavg_center = 0.0
      myradmsqr_center = 0.0
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
      myradmsqr = myradmsqr + part%radius**2

     !Want to get droplet statistics only in the interior
     if (part%xp(3) .gt. 0.25*zl .AND. part%xp(3) .lt. 0.75*zl) then
        myradavg_center = myradavg_center + part%radius
        myradmsqr_center = myradmsqr_center + part%radius**2
        numdrop_center = numdrop_center + 1
     end if

      if (part%radius .gt. part%rc) then
         numdrop = numdrop + 1
      else
         numaerosol = numaerosol + 1
      end if

      if (part%radius .gt. myradmax) myradmax = part%radius
      if (part%radius .lt. myradmin) myradmin = part%radius
      if (part%Tp .gt. mytempmax) mytempmax = part%Tp
      if (part%Tp .lt. mytempmin) mytempmin = part%Tp
      if (part%qstar .gt. myqmax) myqmax = part%qinf
      if (part%qstar .lt. myqmin) myqmin = part%qinf

      part => part%next
      end do


      !Compute sums of integer quantities
      intbuf(1) = numpart
      intbuf(2) = numdrop
      intbuf(3) = numaerosol
      intbuf(4) = denum
      intbuf(5) = actnum
      intbuf(6) = num_destroy
      intbuf(7) = num100
      intbuf(8) = num1000
      intbuf(9) = numimpos
      intbuf(10) = numdrop_center

      call mpi_allreduce(intbuf,intbuf_rec,10,mpi_integer,mpi_sum,mpi_comm_world,ierr)

      tnumpart = intbuf_rec(1)
      tnumdrop = intbuf_rec(2)
      tnumaerosol = intbuf_rec(3)
      tdenum = intbuf_rec(4)
      tactnum = intbuf_rec(5)
      tnum_destroy = intbuf_rec(6)
      tnum100 = intbuf_rec(7)
      tnum1000 = intbuf_rec(8)
      tnumimpos = intbuf_rec(9)
      tnumdrop_center = intbuf_rec(10)

      
      !Compute sums of real quantities

      tmpbuf(1) = myRep_avg
      tmpbuf(2) = mylwc_sum
      tmpbuf(3) = myphiw_sum
      tmpbuf(4) = myphiv_sum
      tmpbuf(5) = myradavg
      tmpbuf(6) = myradmsqr
      tmpbuf(7) = avgres
      tmpbuf(8) = myradavg_center
      tmpbuf(9) = myradmsqr_center

      call mpi_allreduce(tmpbuf,tmpbuf_rec,9,mpi_real8,mpi_sum,mpi_comm_world,ierr)

      Rep_avg = tmpbuf_rec(1)
      lwc = tmpbuf_rec(2)
      phiw = tmpbuf_rec(3)
      phiv = tmpbuf_rec(4)
      radavg = tmpbuf_rec(5)
      radmsqr = tmpbuf_rec(6)
      tavgres = tmpbuf_rec(7)
      radavg_center = tmpbuf_rec(8)
      radmsqr_center = tmpbuf_rec(9)



      phiw = phiw/xl/yl/zl/rhoa
      phiv = phiv/xl/yl/zl
      Rep_avg = Rep_avg/tnumpart
      radavg = radavg/tnumpart
      radmsqr = radmsqr/tnumpart
      tavgres = tavgres/tnum_destroy
      radavg_center = radavg_center/tnumdrop_center
      radmsqr_center = radmsqr_center/tnumdrop_center

      !Min and max radius
      call mpi_allreduce(myradmin,radmin,1,mpi_real8,mpi_min,mpi_comm_world,ierr)
      call mpi_allreduce(myradmax,radmax,1,mpi_real8,mpi_max,mpi_comm_world,ierr)

      call mpi_allreduce(mytempmin,tempmin,1,mpi_real8,mpi_min,mpi_comm_world,ierr)
      call mpi_allreduce(mytempmax,tempmax,1,mpi_real8,mpi_max,mpi_comm_world,ierr)

      call mpi_allreduce(myqmin,qmin,1,mpi_real8,mpi_min,mpi_comm_world,ierr)
      call mpi_allreduce(myqmax,qmax,1,mpi_real8,mpi_min,mpi_comm_world,ierr)



  end subroutine particle_update_BE

  subroutine destroy_particle
      implicit none

      type(particle), pointer :: tmp

      !Is it the first and last in the list?
      if (associated(part,first_particle) .AND. (.NOT. associated(part%next)) ) then
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

  subroutine particle_stats
      use pars
      use con_stats
      use con_data
      implicit none
      integer :: i,ipt,jpt,kpt
      real :: rhop,pi

      part => first_particle
      do while (associated(part))     

      ipt = floor(part%xp(1)/dx) + 1
      jpt = floor(part%xp(2)/dy) + 1
      kpt = minloc(z,1,mask=(z.gt.part%xp(3))) - 1

      pi   = 4.0*atan(1.0)

      rhop = (part%m_s+4.0/3.0*pi*part%radius**3*rhow)/(4.0/3.0*pi*part%radius**3)

      !Takes in ipt,jpt,kpt as the node to the "bottom left" of the particle
      !(i.e. the node in the negative direction for x,y,z)
      !and computes quantities needed to get particle statistics

      partcount_t(kpt,jpt,ipt) = partcount_t(kpt,jpt,ipt) + 1.0
      
      !Get su mean, mean-squared of particle velocities at each level
      upwp_t(kpt,jpt,ipt) = upwp_t(kpt,jpt,ipt) + part%vp(1)*part%vp(3)
      do i = 1,3
      vpsum_t(kpt,jpt,ipt,i) = vpsum_t(kpt,jpt,ipt,i) + part%vp(i)
      vpsqrsum_t(kpt,jpt,ipt,i)=vpsqrsum_t(kpt,jpt,ipt,i)+part%vp(i)**2

      ufsum_t(kpt,jpt,ipt,i) = ufsum_t(kpt,jpt,ipt,i) + part%uf(i)
      ufsqrsum_t(kpt,jpt,ipt,i)=ufsqrsum_t(kpt,jpt,ipt,i)+part%uf(i)**2
      end do

      Tpsum_t(kpt,jpt,ipt) = Tpsum_t(kpt,jpt,ipt) + part%Tp
      Tpsqrsum_t(kpt,jpt,ipt) = Tpsqrsum_t(kpt,jpt,ipt) + part%Tp**2

      Tfsum_t(kpt,jpt,ipt) = Tfsum_t(kpt,jpt,ipt) + part%Tf

      qfsum_t(kpt,jpt,ipt) = qfsum_t(kpt,jpt,ipt) + part%qinf

      wpTpsum_t(kpt,jpt,ipt) = wpTpsum_t(kpt,jpt,ipt) + part%Tp*part%vp(3)


      radsum_t(kpt,jpt,ipt) = radsum_t(kpt,jpt,ipt) + part%radius 

      rad2sum_t(kpt,jpt,ipt) = rad2sum_t(kpt,jpt,ipt) + part%radius**2  

      multcount_t(kpt,jpt,ipt) = multcount_t(kpt,jpt,ipt) + real(part%mult)

      mwsum_t(kpt,jpt,ipt) = mwsum_t(kpt,jpt,ipt) + real(part%mult)*(rhow*4.0/3.0*pi*part%radius**3)

      qstarsum_t(kpt,jpt,ipt) = qstarsum_t(kpt,jpt,ipt) + part%qstar

      part => part%next
      end do

  end subroutine particle_stats

  subroutine particle_write_traj(it)
   use con_data
   use pars
   implicit none

   integer :: it

   
   part => first_particle
   do while (associated(part))
      
      if (mod(part%pidx,4000) .eq. 0) then
          write(ntraj,'(2i,14e15.6)') part%pidx,part%procidx,time,part%xp(1),part%xp(2),part%xp(3),part%vp(1),part%vp(2),part%vp(3),part%radius,part%Tp,part%Tf,part%qinf,part%qstar,part%rc,part%numact
      end if

   part => part%next
   end do

   if (it .ge. itmax) then
      close(ntraj)
   end if

  end subroutine particle_write_traj

  subroutine particle_coalesce
      use pars
      use kd_tree
      use pars
      use con_data
      implicit none 

      type(particle), pointer :: part_tmp
      type(tree_master_record), pointer :: tree

      real, allocatable :: xp_data(:,:),distances(:),rad_data(:)
      real, allocatable :: kappa_s_data(:), ms_data(:)
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
      real ::   kappa_s_j_temp, kappa_s_k_temp,  ms_j_temp,  ms_k_temp

      !For whatever reason, the kd-search sometimes misses edge cases,
      !and you should do nq+1 if you actually want nq
      nq = 11

      allocate(xp_data(numpart,3),index_data(numpart,2))
      allocate(vel_data(numpart,3))
      allocate(distances(nq),indexes(nq),ran_nq(2:nq-1))
      allocate(rad_data(numpart),mult_data(numpart),coal_data(numpart))
      allocate(destroy_data(numpart))
      allocate(kappa_s_data(numpart), ms_data(numpart))

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

         kappa_s_data(i) = part_tmp%kappa_s
         ms_data(i) = part_tmp%m_s

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
         
         veldiff = sqrt( (vel_data(i,1)-vel_data(coal_idx,1))**2 +  &
                         (vel_data(i,2)-vel_data(coal_idx,2))**2 +  &
                         (vel_data(i,3)-vel_data(coal_idx,3))**2 )

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
!         K = pi2/2.0*E*veldiff*(rad_data(i) + rad_data(coal_idx))**2

         if (ikernel.eq.0) then ! Golovin (1963) kernel
            pvol_j = pi2*2.0/3.0*rad_data(j_idx)**3.0
            pvol_k = pi2*2.0/3.0*rad_data(k_idx)**3.0
            golovin_b = 1.5e3
            K = golovin_b*(pvol_j + pvol_k)
         elseif (ikernel.eq.1) then ! Geometric kernel with collection efficiencies E=1
            E = 1.
            K = pi2/2.0*E*veldiff*(rad_data(i) + rad_data(coal_idx))**2
         elseif (ikernel.eq.2) then ! Long (1974) polynomial kernel (as in Bott (1998))
            E = long_effic(rad_data(j_idx)*1e6,rad_data(j_idx)*1e6)
            K = pi2/2.0*E*veldiff*(rad_data(i) + rad_data(coal_idx))**2
         elseif (ikernel.eq.3) then ! Geometric kernel with modified Hall collection efficiencies (as in Bott 1998)
            E = hall_effic(rad_data(i) * 1e6, rad_data(coal_idx)*1e6) ! input needs to be in [um]
            K = pi2/2.0*E*veldiff*(rad_data(i) + rad_data(coal_idx))**2
         elseif (ikernel.eq.4) then ! Geometric kernel with turbulence-enhanced Hall efficiencies (Wang and Grabowski 2009)
            !!! TO BE IMPLEMENTED SOON
         endif


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
               rad_data(k_idx) = (gam_til*rad_data(j_idx)**3 + rad_data(k_idx)**3)**(1.0/3.0)

               !Update particle k's kappa coefficient
			      !We can use the mass mixing rule for now (will be updated to volume mixing rule)
			      kappa_s_data(k_idx) = gam_til*ms_data(j_idx)*kappa_s_data(j_idx)/( gam_til*ms_data(j_idx) + ms_data(k_idx))  +   &
                ms_data(k_idx)*kappa_s_data(k_idx) / ( gam_til*ms_data(j_idx) + ms_data(k_idx))
 
               !Update particle k's solute mass
               ms_data(k_idx) = gam_til*ms_data(j_idx) + ms_data(k_idx)

            elseif (xi_j - gam_til*xi_k .eq. 0) then

               mult_tmp_j = floor(real(mult_data(k_idx))/2.0)
               mult_tmp_k = mult_data(k_idx) - floor(real(mult_data(k_idx))/2.0)
       
               mult_data(j_idx) = mult_tmp_j
               mult_data(k_idx) = mult_tmp_k

               rad_j_tmp = (gam_til*rad_data(j_idx)**3 + rad_data(k_idx)**3)**(1.0/3.0)
               rad_k_tmp = (gam_til*rad_data(j_idx)**3 + rad_data(k_idx)**3)**(1.0/3.0)

               rad_data(j_idx) = rad_j_tmp
               rad_data(k_idx) = rad_k_tmp

               !Update kappa coefficients
               kappa_s_j_temp = gam_til*ms_data(j_idx)*kappa_s_data(j_idx)/(gam_til*ms_data(j_idx) + ms_data(k_idx)) +  &
                ms_data(k_idx)*kappa_s_data(k_idx) / ( gam_til*ms_data(j_idx) + ms_data(k_idx))
   
               kappa_s_k_temp = gam_til*ms_data(j_idx)*kappa_s_data(j_idx)/(gam_til*ms_data(j_idx) + ms_data(k_idx)) +  &
                ms_data(k_idx)*kappa_s_data(k_idx) / ( gam_til*ms_data(j_idx) + ms_data(k_idx))
   
               kappa_s_data(j_idx) = kappa_s_j_temp
               kappa_s_data(k_idx) = kappa_s_k_temp
   
               !Update particle's solute masses
               ms_j_temp = gam_til*ms_data(j_idx) + ms_data(k_idx)
               ms_k_temp = gam_til*ms_data(j_idx) + ms_data(k_idx)
   
               ms_data(j_idx) = ms_j_temp
               ms_data(k_idx) = ms_k_temp

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

         !Now we also change solute mass and kappa coefficient
         part_tmp%kappa_s = kappa_s_data(i)
         part_tmp%m_s = ms_data(i)

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
      deallocate(kappa_s_data, ms_data)

  end subroutine particle_coalesce

  subroutine gauss_newton_2d(vnext,h,vec1,vec2,flag)
        implicit none

        real, intent(in) :: vnext(3), h, vec1(2)
        real, intent(out) :: vec2(2)
        integer, intent(out) :: flag
        real :: error,fv1(2),fv2(2),v1(2),v_output(3),rel,det
        real :: diff, temp1(2), temp2(2), relax, coeff, correct(2)
        real, dimension(1:2, 1:2) :: J, fancy, inv, finalJ
        integer :: iterations,neg,counts,iteration_max

        iterations = 0
        flag = 0
        error = 1.0e-8

        v1 = vec1
        fv2 = (/1., 1./) 
        coeff = 0.1
        correct = 0.0
        iteration_max = 50

        do while (iterations<iteration_max)

                iterations = iterations + 1

                call ie_vrt_nd(vnext,v1(1),v1(2),v_output,fv1,h)
                call jacob_approx_2d(vnext,v1(1),v1(2),h,J)

                fancy = matmul(transpose(J),J)

                det = fancy(1,1)*fancy(2,2)-fancy(1,2)*fancy(2,1)
                if (abs(det) .lt. 1.0e-10) then
                   flag = 1
                   EXIT
                end if

                call inverse_finder_2d(fancy,det,inv)

                finalJ = matmul(inv,transpose(J))

                correct = matmul(finalJ,fv1)
                vec2 = v1 - correct

                do while ((vec2(1)<0) .OR. (vec2(2)<0) .OR. isnan(vec2(1)))
                        counts = counts + 1
                        coeff = 0.5
                        relax = relax * coeff
                        vec2 = v1-matmul(finalJ,fv1)*relax
                        if (counts>10) EXIT
                end do

                !Excessively low or negative values are a sign that iterations have failed
                if (vec2(1)<0.01 .OR. vec2(2)<0.01 .OR. isnan(vec2(1)) .OR. isnan(vec2(2))) then
                   flag = 1
                   EXIT
                end if

                !Successful completion
                if (sqrt(dot_product(correct,correct))<error) then
                        EXIT
                end if

                relax = 1.0
                counts = 0

                v1 = vec2

        end do
      !A final check to make sure that problem droplets get flagged
      if (iterations == iteration_max) flag = 1
      if (isnan(vec2(1)) .OR. vec2(1)<0 .OR. isnan(vec2(2)) .OR. vec2(2)<0) flag = 1

  end subroutine gauss_newton_2d
  subroutine LV_solver(vnext,h,vec1,vec2,flag)
        implicit none

        real, intent(in) :: vnext(3),h, vec1(2)
        real, intent(out) :: vec2(2)
        integer, intent(out) :: flag
        real :: error,fv1(2),fv2(2),v1(2),v_output(3),rel,det
        real :: diff, lambda,lup,ldown
        real :: C(2), newC(2), gradC(2), correct(2)
        real, dimension(1:2, 1:2) :: J,I,g,invg
        integer :: iterations,neg,iterations_max

        I = reshape((/1, 0, 0, 1/),shape(I))
        iterations = 0
        flag = 0
        v1 = vec1
        fv2 = (/1., 1./)
        error = 1.0e-8

        lambda = 0.001
        lup = 2.0
        ldown = 2.0

        iterations_max = 1000

        do while (iterations<iterations_max)

        iterations = iterations + 1

        call jacob_approx_2d(vnext, v1(1), v1(2), h,J)
        call ie_vrt_nd(vnext, v1(1), v1(2),v_output,fv1,h)

        g = matmul(transpose(J),J)+lambda*I
        gradC = matmul(transpose(J),fv1)
        C = 0.5*fv1*fv1

        det = g(1,1)*g(2,2)-g(1,2)*g(2,1)
        if (abs(det) .lt. 1.0e-10) then
           flag = 1
           EXIT
        end if

        call inverse_finder_2d(g,det,invg)
        correct = matmul(invg, gradC)


        vec2 = v1 - correct

        !Successful completion
        if (sqrt(dot_product(correct,correct)) < error) then
                EXIT
        end if

        if (vec2(1) .lt. 0.0 .OR. vec2(2) .lt. 0.0 .OR. isnan(vec2(1)) .OR. isnan(vec2(2))) then
           flag = 1
           EXIT
        end if

        call ie_vrt_nd(vnext, vec2(1), vec2(2),v_output,fv2,h)
        newC = 0.5*fv2*fv2


        if (sqrt(dot_product(newC,newC))<sqrt(dot_product(C,C))) then
                v1 = vec2
                lambda = lambda/ldown
        else
                lambda = lambda*lup
        end if

        end do

        if (iterations==iterations_max) then
                flag = 1
        end if

        if (vec2(1) < 0 .OR. vec2(2) < 0) then
                flag = 1
        end if


  end subroutine LV_solver
  subroutine jacob_approx_2d(vnext, rnext, tnext, h, J)
        implicit none
        integer :: n

        real, intent(in) :: vnext(3), rnext, tnext, h
        real, intent(out), dimension(1:2, 1:2) :: J
        real :: diff = 0, v_output(3), rt_output(2),xper(2),fxper(2), ynext(2),xper2(2),fxper2(2)

        diff = 1E-12

        ynext(1) = rnext
        ynext(2) = tnext

        call ie_vrt_nd(vnext, rnext, tnext, v_output, rt_output, h)

        xper = ynext
        xper2 = ynext

        do n=1, 2
                xper(n) = xper(n) + diff
                xper2(n) = xper2(n) - diff
                call ie_vrt_nd(vnext, xper(1), xper(2),v_output,fxper,h)
                call ie_vrt_nd(vnext, xper2(1), xper2(2),v_output,fxper2,h)
                J(:, n) = (fxper-rt_output)/diff
                xper(n) = ynext(n)
                xper2(n) = ynext(n)
        end do

  end subroutine jacob_approx_2d
  subroutine inverse_finder_2d(C,det,invC)
        implicit none
        real, intent(in) :: det
        real, dimension(1:2, 1:2), intent(in) :: C
        real, dimension(1:2, 1:2), intent(out) :: invC

        invC = reshape((/C(2, 2), -C(2,1), -C(1, 2), C(1, 1)/),shape(invC))
        invC = (1./det)*invC

  end subroutine inverse_finder_2d
  subroutine ie_vrt_nd(vnext, tempr, tempt, v_output,rt_output, h)
      use pars
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'

      real, intent(in) :: vnext(3), tempr, tempt, h
      real, intent(out) :: v_output(3), rT_output(2)

      real :: esa, dnext,  m_w, rhop, Rep, taup,vprime(3), rprime, Tprime, qstr, Shp, Nup, dp, VolP
      real :: diff(3), diffnorm, Tnext, rnext, T
      real :: taup0, g(3)
      real :: mod_Magnus


        taup0 = (((part%m_s)/((2./3.)*pi2*radius_init**3) + rhow)*(radius_init*2)**2)/(18*rhoa*nuf)
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
        vprime(1:3) = (1. + 0.15 * (Rep**0.687)) * (1./taup)*diff(1:3) - g(1:3)
        vprime(1:3) = vprime(1:3) * taup0 ** 2
        !!!!!!!!!!!!!!!!

        !!! Humidity !!!
        qstr = (Mw/(Ru*Tnext*rhoa)) * esa * exp(((Lv*Mw/Ru)*((1./part%Tf) - (1./Tnext))) + ((2.*Mw*Gam)/(Ru*rhow*rnext*Tnext)) - ((part%kappa_s*part%m_s*rhow/rhos)/(Volp*rhop-part%m_s)))
        !!!!!!!!!!!!!!!!!!

        !!! Radius !!!
        Shp = 2. + 0.6 * Rep**(1./2.) * Sc**(1./3.)
        rprime = (1./9.) * (Shp/Sc) * (rhop/rhow) * (rnext/taup) * (part%qinf - qstr)
        rprime = rprime * (taup0/part%radius)
        !!!!!!!!!!!!!!!!!

        !!! Temperature !!!
        Nup = 2. + 0.6*Rep**(1./2.)*Pra**(1./3.);

        Tprime = -(1./3.)*(Nup/Pra)*CpaCpp*(rhop/rhow)*(1./taup)*(Tnext-part%Tf) + 3.*Lv*(1./(rnext*Cpp))*rprime*(part%radius/taup0)
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

      a = -(2*Mw*Gam)/(Ru*rhow*part%Tf)/LOG((Ru*part%Tf*rhoa*part%qinf)/(Mw*esa))
      c = (part%kappa_s*part%m_s)/((2.0/3.0)*pi2*rhos)/LOG((Ru*part%Tf*rhoa*part%qinf)/(Mw*esa))

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

  subroutine SFS_velocity
  !This subroutine calculate the SFS velocity for particles
  !Uses Weil et al. (2004) formulation
  use pars
  use fields
  use fftwk
  use con_data
  use con_stats
  implicit none
  include 'mpif.h'
  real :: sigm_sdxp,sigm_sdyp,sigm_sdzp,vis_sp
  real :: sigm_su,sigm_sl,us_ran,gasdev,tengz,englez_bar
  real :: engsbz_bar,sigm_w, sigm_ws
  real :: L_flt,epsn,fs,C0,a1,a2,a3,sigm_sprev,fs1
  real :: weit,weit1,weit3,weit4, T_lagr
  real :: us(3)
  real :: xp3i,Volp,rhop
  integer :: ix,iy,iz,izp1,izm1,ind,iz_part,ierr
  integer :: fluxloc,fluxloci

!       ---initialize -------
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
  pfluxdiff = 0.0

!       ------------------
!       compute sigma squre (sigm_s) based on subgrid energy field
!       -----------------       
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
       vis_ss(ix,iy,iz) = vis_s(ix,iy,1,iz)
       sigm_sdx(ix,iy,iz) =  sigm_s(ix,iy,iz)       !for xderiv
       sigm_sdy(ix,iy,iz) =  sigm_s(ix,iy,iz)       !for yderiv

!     --------------------
!     calculate z derivative of sigma_s
!     this will be at the u-point!
!     --------------------
       sigm_sdz(ix,iy,iz)=(sigm_s(ix,iy,iz)-sigm_s(ix,iy,izm1))*dzw_i(iz)
      end do
      end do

!     -------------------
!     calculate x derivatives of sigma_s
!     -------------------
     call xderivp(sigm_sdx(1,iys,iz),trigx(1,1),xk(1),nnx,iys,iye)

  end do

!     ------------------
!     calculate y derivative of sigma_s
!     ------------------

  call yd_mpi(sigm_sdy(1,iys,izs),trigx(1,2),yk(1),nnx,nny,ixs,ixe,ix_s,ix_e,iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)

!     ----------------
!    calculate extented fileds of sigm_s and its derivatives
!     ---------------

  call fill_extSFS

  !Loop over the linked list of particles:
  part => first_particle
  do while (associated(part))

    ! interpolate sigm_s and its derivative at particle location
    sigm_sprev = part%sigm_s
    call sigm_interp(sigm_sdxp,sigm_sdyp,sigm_sdzp,vis_sp,iz_part)

    part%sigm_s = abs(part%sigm_s)  !Interpolation near surface can give small negative numbers
    part%sigm_s = max(part%sigm_s,1.0e-4) !Prevent it from getting too small, makes time derivative term singular

!     -----------------
!     calculate the subgrid velocity Weil et al ,2004-isotropic turb.
!     ----------------
     l_flt = (2.25*dx*dy*dzw(iz_part+1))**(1.0/3.0)! filtered with
     epsn = (0.93/l_flt)*(3.0*part%sigm_s/2.0)**(1.5) ! tur. dis.rt
     ! TKE : resolved + Subgrid at grid center
      !---------------------------------
!      tot_eng = (engsbz(iz_part+1) + engz(iz_part+1))
     englez_bar = 0.5*(englez(iz_part)+englez(iz_part+1))
     engsbz_bar = 0.5*(engsbz(iz_part)+engsbz(iz_part+1))
     tengz = englez_bar + engsbz_bar
    !---------------------------------
    ! Calculate fs basd on w-componet of velocity
     sigm_w = 0.5*(wps(iz_part)+wps(iz_part+1))
     sigm_ws = 0.5*(engsbz(iz_part)+engsbz(iz_part+1))/3.0

!     ---------write for single droplet ---------
!       write(*,*)'sigm_w:', sigm_w,sigm_ws
!     ------------------------------------------
    if(tengz.gt.0.0)then
!      fs = engsbz(iz_part+1)/(engsbz(iz_part+1) + engz(iz_part+1))
!      fs = engsbz_bar/tengz
       fs = sigm_ws/(sigm_w + sigm_ws)
    else
       fs =0.0
    end if
    C0 = 6.0  ! Changed to 6.0 from 3.0 Indrjith 11-20-17

!     ---------Check for single-part----------
       T_lagr = 2*part%sigm_s/(C0*epsn)   ! Lagrangian time scale
!       write(*,*) 'L_time:',part%xp(3),T_lagr
!     -----------------------------------
!    ------------------
!      Calculate subgrid velocity components
!     -----------------
     us(1:3) = part%u_sub(1:3)
!    -----x component ----------------
     a1 = 0.0
     a2 =0.0
     a3 = 0.0
     a1 =(-0.5)*fs*C0*epsn*part%u_sub(1)/part%sigm_s
     a2 =0.5*(fs/part%sigm_s)*part%u_sub(1)*(part%sigm_s - sigm_sprev)/dt
     a3 = 0.5*fs*sigm_sdxp
     us_ran = sqrt(fs*C0*epsn*dt)*gasdev(iseed)
     part%u_sub(1) = (a1+a2+a3)*dt + us_ran

!     ----------------------------------
!     -----y component ----------------
    a1 = 0.0
    a2 =0.0
    a3 = 0.0
    a1 =(-0.5)*fs*C0*epsn*part%u_sub(2)/part%sigm_s
    a2 =0.5*(fs/part%sigm_s)*part%u_sub(2)*(part%sigm_s - sigm_sprev)/dt
    a3 = 0.5*fs*sigm_sdyp
    us_ran = sqrt(fs*C0*epsn*dt)*gasdev(iseed)
    part%u_sub(2) = (a1+a2+a3)*dt + us_ran

!     ----------------------------------
!     -----z component ----------------
    a1 = 0.0
    a2 =0.0
    a3 = 0.0
    a1 =(-0.5)*fs*C0*epsn*part%u_sub(3)/part%sigm_s
    a2 =0.5*(fs/part%sigm_s)*part%u_sub(3)*(part%sigm_s - sigm_sprev)/dt
    a3 = 0.5*fs*sigm_sdzp
    us_ran = sqrt(fs*C0*epsn*dt)*gasdev(iseed)
    part%u_sub(3) = (a1+a2+a3)*dt + us_ran

!     -------------------
!     Update particle location and velocity
!     ------------------
    xp3i = part%xp(3)
    do ind = 1,3
      part%xp(ind) = part%xp(ind) + part%u_sub(ind)*dt
    end do
!     ---------------------------------        

    Volp = pi2*2.0/3.0*part%radius**3
    rhop = (part%m_s+Volp*rhow)/Volp

    !Store the particle flux now that we have the new position
    if (part%xp(3) .gt. zl) then   !This will get treated in particle_bcs_nonperiodic, but record here
       fluxloc = nnz+1
       fluxloci = minloc(z,1,mask=(z.gt.xp3i))-1
    elseif (part%xp(3) .lt. 0.0) then !This will get treated in particle_bcs_nonperiodic, but record here
       fluxloci = minloc(z,1,mask=(z.gt.xp3i))-1
       fluxloc = 0
    else
       fluxloc = minloc(z,1,mask=(z.gt.part%xp(3)))-1
       fluxloci = minloc(z,1,mask=(z.gt.xp3i))-1
    end if  !Only apply flux calc to particles in domain

    if (xp3i .lt. part%xp(3)) then !Particle moved up

      do iz=fluxloci,fluxloc-1
         pfluxdiff(iz) = pfluxdiff(iz) + part%mult
         pmassflux(iz) = pmassflux(iz) + rhop*Volp*part%mult
         penegflux(iz) = penegflux(iz) + rhop*Volp*Cpp*part%Tp*part%mult
      end do

    elseif (xp3i .gt. part%xp(3)) then !Particle moved down

      do iz=fluxloc,fluxloci-1
         pfluxdiff(iz) = pfluxdiff(iz) - part%mult
         pmassflux(iz) = pmassflux(iz) - rhop*Volp*part%mult
         penegflux(iz) = penegflux(iz) - rhop*Volp*Cpp*part%Tp*part%mult
      end do

    end if  !Up/down conditional statement


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
  call mpi_allreduce(numpart,tnumpart,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)

  end subroutine SFS_velocity

  subroutine SFS_position
  !Use a more simplistic SFS treatment: Stochastic particle position rather than velocity
  !Designed to be consistent with LES subgrid eddy diffusivity
  !Does not use Weil et al. 2004 formulation at all
  use pars
  use fields
  use fftwk
  use con_data
  use con_stats
  implicit none
  include 'mpif.h'

  real :: sigm_sdxp,sigm_sdyp,sigm_sdzp,vis_sp
  real :: phim,phis,psim,psis,zeta
  real :: dadz,gasdev
  real :: xp3i
  integer :: ix,iy,iz,izp1,izm1,ind,iz_part,ierr
  integer :: fluxloc,fluxloci

    sigm_s = 0.0
    pfluxdiff = 0.0

!       ------------------
!       compute sigma squre (sigm_s) based on subgrid energy field
!       -----------------       
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

      xp3i = part%xp(3)

      !Now simply solve Langevin equation:
      part%xp(1) = part%xp(1) + gasdev(iseed)*sqrt(2.0*abs(vis_sp)*dt)
      part%xp(2) = part%xp(2) + gasdev(iseed)*sqrt(2.0*abs(vis_sp)*dt)
      part%xp(3) = part%xp(3) + gasdev(iseed)*sqrt(2.0*abs(vis_sp)*dt) + dadz*dt


      !Store the particle flux now that we have the new position
      if (part%xp(3) .gt. zl) then   !This will get treated in particle_bcs_nonperiodic, but record here
         fluxloc = nnz+1
         fluxloci = minloc(z,1,mask=(z.gt.xp3i))-1
      elseif (part%xp(3) .lt. 0.0) then !This will get treated in particle_bcs_nonperiodic, but record here
         fluxloci = minloc(z,1,mask=(z.gt.xp3i))-1
         fluxloc = 0
      else

      fluxloc = minloc(z,1,mask=(z.gt.part%xp(3)))-1
      fluxloci = minloc(z,1,mask=(z.gt.xp3i))-1

      end if  !Only apply flux calc to particles in domain

      if (xp3i .lt. part%xp(3)) then !Particle moved up

      do iz=fluxloci,fluxloc-1
         pfluxdiff(iz) = pfluxdiff(iz) + part%mult
      end do

      elseif (xp3i .gt. part%xp(3)) then !Particle moved down

      do iz=fluxloc,fluxloci-1
         pfluxdiff(iz) = pfluxdiff(iz) - part%mult
      end do

      end if  !Up/down conditional statement


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
  call mpi_allreduce(numpart,tnumpart,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)


  end subroutine SFS_position

  subroutine particle_neighbor_search_kd
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
     write(*,'(a30,3i)') 'WARNING in neighbor,idx:',idx_check,first_particle%nbr_pidx,myid
  end if
  if (dist_diff .gt. 1.0e-8) then
     write(*,'(a30,2e15.6,i)') 'WARNING in neighbor,dist:',dist_check,first_particle%dist,myid
  end if

  deallocate(xp_data,index_data)


  end subroutine particle_neighbor_search_kd

  subroutine particle_neighbor_search_brute
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


  subroutine set_binsdata(binsdata,sizea,rmin,rmax)
  use pars
  use fields
  use con_data
  use con_stats
  implicit none

  integer :: i,nbin,ibin,nbinnew
  real :: dh

  integer,intent(in) :: sizea
  real,intent(inout) :: binsdata(sizea)

  real :: rmin,rmax,rmin10,rmax10

    nbin = histbins !From Module Particle 

    rmin10 = log10(rmin)
    rmax10 = log10(rmax)

!       Calculate size of interval
    dh = (rmax10-rmin10)/(nbin-1)

!       ===== update x-axis for each bin ===== 
    binsdata(1) = rmin10-dh
    do i = 1,histbins+1
      binsdata(i+1)= dh+binsdata(i)
    end do


  end subroutine set_binsdata

  subroutine set_binsdata_integer(binsdata,sizea,rmin)
  use pars
  use fields
  use con_data
  use con_stats
  implicit none

  integer :: i,nbin,ibin,nbinnew
  real :: dh

  integer,intent(in) :: sizea
  real,intent(inout) :: binsdata(sizea)

  real :: rmin

    nbin = histbins !From Module Particle 

!       Calculate size of interval
    dh = 1.0

!       ===== update x-axis for each bin ===== 
    binsdata(1) = rmin-dh
    do i = 1,histbins+1
      binsdata(i+1)= dh+binsdata(i)
    end do


  end subroutine set_binsdata_integer

  subroutine add_histogram(binsdata,histdata,sizea,val1,mult)

  use pars
  use fields
  use con_data
  use con_stats
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

!       Calculate size of interval
    dh = (rmax-rmin)/(nbin-1)

    logval1 = log10(val1)
    if (logval1 .gt. rmax+0.5*dh) then
            ibin = nbin + 2
    elseif (logval1 .lt. rmin-0.5*dh) then
            ibin = 1
    else
            ibin = (floor((logval1-(rmin-0.5*dh))/dh)+1)+1
    end if

!       Add the current event to the histogram
    histdata(ibin) = histdata(ibin) + dble(mult)

  end subroutine add_histogram

  subroutine add_histogram_integer(binsdata,histdata,sizea,val)

  use pars
  use fields
  use con_data
  use con_stats
  implicit none

  integer :: i,nbin,nbinnew,ibin
  real :: dh
  real,intent(in) :: val

  integer,intent(in) :: sizea
  real,intent(in) :: binsdata(sizea)
  real,intent(inout) :: histdata(sizea)

  real :: rmin,rmax

    rmin = binsdata(2)
    nbin = histbins !From Module Particle 
    rmax = binsdata(nbin+1)

!       Calculate size of interval
    dh = (rmax-rmin)/(nbin-1)

    if (val .gt. rmax+0.5*dh) then
            ibin = nbin + 2
    elseif (val .lt. rmin-0.5*dh) then
            ibin = 1
    else
            ibin = (floor((val-(rmin-0.5*dh))/dh)+1)+1
    end if


!       Add the current event to the histogram
    histdata(ibin) = histdata(ibin) + 1.0

  end subroutine add_histogram_integer

  subroutine radius_histogram
  implicit none

  hist_rad = 0.0
  part => first_particle
  do while (associated(part))
     call add_histogram(bins_rad,hist_rad,histbins+2,part%radius,part%mult)
     part => part%next
  end do

  end subroutine radius_histogram

   subroutine eigen_roots(a,m,rtr,rti)
   ! USES balanc,hqr
   !Find all the roots of a polynomial with real coeﬃcients. The method is to construct an upper Hessenberg matrix
   !whose eigenvalues are the desired roots, and then use the routines balanc and hqr . The real and imaginary parts of the roots are returned in rtr(1:m) and rti(1:m) , respectively.
   
      INTEGER :: m,MAXM
      real :: a(m+1),rtr(m),rti(m)
      PARAMETER (MAXM=50)
      INTEGER :: j,k
      real :: hess(MAXM,MAXM),xr,xi
   
      if (m.gt.MAXM.or.a(m+1).eq.0.) then
         write(*,*)'bad args in zrhqr'
         stop
      end if
   
      do k=1,m
         hess(1,k)=-a(m+1-k)/a(m+1)
         do j=2,m
            hess(j,k)=0.
         enddo
         if (k.ne.m) hess(k+1,k)=1.
      enddo
   
      call balanc(hess,m,MAXM)
   
      call hqr(hess,m,MAXM,rtr,rti)
   
      do j=2,m
         xr=rtr(j)
         xi=rti(j)
         do k=j-1,1,-1
            if(rtr(k).le.xr) goto 1
            rtr(k+1)=rtr(k)
            rti(k+1)=rti(k)
         enddo
   
         k=0
1   		rtr(k+1)=xr
         rti(k+1)=xi
      enddo
      return
   end subroutine
   
   
   SUBROUTINE balanc(a,n,np)
   !Given an n by n matrix a stored in an array of physical dimensions np by np , this routine
   !replaces it by a balanced matrix with identical eigenvalues. A symmetric matrix is already
   !balanced and is unaﬀected by this procedure. The parameter RADIX should be the machine’s
   !ﬂoating-point radix.
      INTEGER n,np
      real a(np,np),RADIX,SQRDX
      PARAMETER (RADIX=2.,SQRDX=RADIX**2)
      INTEGER i,j,last
      real c,f,g,r,s
   
2   	continue
   
      last=1
      do i=1,n
         c=0.
         r=0.
         do j=1,n
            if(j.ne.i)then
               c=c+abs(a(j,i))
               r=r+abs(a(i,j))
            endif
         enddo
   
         if(c.ne.0..and.r.ne.0.)then
            g=r/RADIX
            f=1.
            s=c+r
3   			if(c.lt.g)then
               f=f*RADIX
               c=c*SQRDX
               goto 3
            endif
   
            g=r*RADIX
4   			if(c.gt.g)then
               f=f/RADIX
               c=c/SQRDX
               goto 4
            endif
   
            if((c+r)/f.lt.0.95*s)then
               last=0
               g=1./f
               do j=1,n
                  a(i,j)=a(i,j)*g
               enddo
   
               do j=1,n
                  a(j,i)=a(j,i)*f
               enddo
            endif
         endif
      enddo
   
      if(last.eq.0) goto 2
      return
   END subroutine
   
   SUBROUTINE hqr(a,n,np,wr,wi)
   !Finds all eigenvalues of an n by n upper Hessenberg matrix a that is stored in an np by np
   !array. On input a can be exactly as output from elmhes; on output it is destroyed.
   !The real and imaginary parts of the eigenvalues are returned in wr and wi , respectively.
      INTEGER n,np
      double precision a(np,np),wi(np),wr(np)
      INTEGER i,its,j,k,l,m,nn
      double precision anorm,p,q,r,s,t,u,v,w,x,y,z
   
      anorm=0.
   
      do i=1,n
         do j=max(i-1,1),n
            anorm=anorm+abs(a(i,j))
         enddo
      enddo
   
      nn=n
      t=0.
   
5   	if(nn.ge.1)then
         its=0
6   		do l=nn,2,-1
            s=abs(a(l-1,l-1))+abs(a(l,l))
            if(s.eq.0.)s=anorm
            if(abs(a(l,l-1))+s.eq.s)goto 7
         enddo
   
         l=1
   
7   		x=a(nn,nn)
   
         if(l.eq.nn)then
            wr(nn)=x+t
            wi(nn)=0.
            nn=nn-1
         else
            y=a(nn-1,nn-1)
            w=a(nn,nn-1)*a(nn-1,nn)
   
            if(l.eq.nn-1)then
               p=0.5*(y-x)
               q=p**2+w
               z=sqrt(abs(q))
               x=x+t
   
               if(q.ge.0.)then
                  z=p+sign(z,p)
                  wr(nn)=x+z
                  wr(nn-1)=wr(nn)
   
                  if(z.ne.0.)wr(nn)=x-w/z
                     wi(nn)=0.
                     wi(nn-1)=0.
                  else
                     wr(nn)=x+p
                     wr(nn-1)=wr(nn)
                     wi(nn)=z
                     wi(nn-1)=-z
                  endif
   
                  nn=nn-2
   
               else
                  if(its.eq.30) then
                     write(*,*) 'too many iterations in hqr'
                     stop
                  endif
   
                  if(its.eq.10.or.its.eq.20)then
                     t=t+x
   
                     do i=1,nn
                        a(i,i)=a(i,i)-x
                     enddo
   
                     s=abs(a(nn,nn-1))+abs(a(nn-1,nn-2))
                     x=0.75*s
                     y=x
                     w=-0.4375*s**2
                  endif
   
                  its=its+1
   
                  do m=nn-2,l,-1
                     z=a(m,m)
                     r=x-z
                     s=y-z
                     p=(r*s-w)/a(m+1,m)+a(m,m+1)
                     q=a(m+1,m+1)-z-r-s
                     r=a(m+2,m+1)
                     s=abs(p)+abs(q)+abs(r)
                     p=p/s
                     q=q/s
                     r=r/s
                     if(m.eq.l)goto 8
                     u=abs(a(m,m-1))*(abs(q)+abs(r))
                     v=abs(p)*(abs(a(m-1,m-1))+abs(z)+abs(a(m+1,m+1)))
                     if(u+v.eq.v)goto 8
                  enddo
   
8   					do i=m+2,nn
                     a(i,i-2)=0.
                     if (i.ne.m+2) a(i,i-3)=0.
                  enddo
   
                  do k=m,nn-1
                     if(k.ne.m)then
                        p=a(k,k-1)
                        q=a(k+1,k-1)
                        r=0.
   
                        if(k.ne.nn-1)r=a(k+2,k-1)
                        x=abs(p)+abs(q)+abs(r)
                        if(x.ne.0.)then
                           p=p/x
                           q=q/x
                           r=r/x
                        endif
                     endif
   
                     s=sign(sqrt(p**2+q**2+r**2),p)
   
                     if(s.ne.0.)then
                        if(k.eq.m)then
                           if(l.ne.m)a(k,k-1)=-a(k,k-1)
                        else
                           a(k,k-1)=-s*x
                        endif
   
                        p=p+s
                        x=p/s
                        y=q/s
                        z=r/s
                        q=q/p
                        r=r/p
   
                        do j=k,nn
                           p=a(k,j)+q*a(k+1,j)
                           if(k.ne.nn-1)then
                              p=p+r*a(k+2,j)
                              a(k+2,j)=a(k+2,j)-p*z
                           endif
                           a(k+1,j)=a(k+1,j)-p*y
                           a(k,j)=a(k,j)-p*x
                        enddo
   
                        do i=l,min(nn,k+3)
                           p=x*a(i,k)+y*a(i,k+1)
                           if(k.ne.nn-1)then
                              p=p+z*a(i,k+2)
                              a(i,k+2)=a(i,k+2)-p*r
                           endif
                           a(i,k+1)=a(i,k+1)-p*q
                           a(i,k)=a(i,k)-p
                        enddo
                     endif
                  enddo
   
                  goto 6
               endif
            endif
            goto 5
         endif
      return
   END subroutine


  function crit_radius(m_s,kappa_s,Tf)
    use pars
    use con_data
    implicit none

    integer :: i,maxidx(1)
    integer, parameter :: N=1000
    real :: m_s,kappa_s,Tf
    real :: radval(N),SS(N)
    real :: radstart,radend,dr
    real :: crit_radius
    real :: firstterm,secterm

    radstart = -8
    radend = -3
    dr = (radstart-radend)/N

    do i=1,N

      radval(i) = 10**(radstart - (i-1)*dr)

      firstterm = 2*Mw*Gam/Ru/rhow/radval(i)/Tf
      secterm =kappa_s*m_s/(rhos*pi2*2.0/3.0*radval(i)**3)

      SS(i) = exp( 2*Mw*Gam/Ru/rhow/radval(i)/Tf - kappa_s*m_s/(rhos*pi2*2.0/3.0*radval(i)**3))

    end do

    maxidx = maxloc(SS)

    crit_radius = radval(maxidx(1))


  end function crit_radius

  function smithssgf(u10,r80,vk)
  implicit none
  real :: u10,r80,cdn10,vk,u14
  real :: a1,a2,smithssgf


  ! added SEA93 function
  if (u10 .ge. 4 .and. u10 .lt. 11) then
     cdn10 = 1.2e-3
  elseif (u10 .ge. 11) then
     cdn10 = 1e-3*(0.49 + 0.065*u10)
  end if

  ! calculate 14 m wind speed
  u14 = u10*(1 + (sqrt(cdn10)/vk)*alog(14.0/10.0))

  a1 = 10**(0.0676*u14 + 2.43)
  a2 = 10**(0.959*sqrt(u14) - 1.476)

  !print*,'tdf',u10,u14,cdn10,a1,a2

  ! calcaulate dFs/dr80 (Eq. A1 in Andreas 1998)
  smithssgf = a1*exp(-3.1*(alog(r80/2.1))**2) + a2*exp(-3.3*(alog(r80/9.2))**2)

  end function smithssgf

! Long (1974) collection efficiencies (as in Bott 1998)
   real function long_effic(r1,r2)
      implicit none
      real :: r1, r2, effic, rbig, rsmall

      if (r1.le.r2) then
         rbig = r2
         rsmall = r1
      else
         rbig = r1
         rsmall = r2
      end if

      ! r units in [um]
      if (rbig.le.50.) then
         effic = 4.5d-4*rbig**2*(1.d0-3.d0/(max(3.d0,dble(rsmall))+1.d-2))
      else
         effic = 1.
      end if

      long_effic = effic

   end function long_effic

   ! modified Hall (1980) collection efficiencies (as in Bott 1998)
	real function hall_effic(r1,r2)
      implicit none
      dimension rat(21),r0(15), ecoll(15,21)
      integer :: k, ir, kk, iq
      real :: drop_ratio, p, q, r1, r2, rbig, rsmall
      real :: effic, r0, rat, ecoll, ek, ec

      ! selecting collector (rbig) and collected droplet (rsmall)
      ! radii in um
      if (r1 .le. r2) then
         rbig = r2
         rsmall = r1
      else
         rbig = r1
         rsmall = r2
      endif

      drop_ratio = rsmall / rbig

      !r0: collector droplet array
      !rat: droplet ratio array
      !ecoll: collection efficiencies

      data r0 /6.,8.,10.,15.,20.,25.,30.,40.,50., &
               60.,70.,100.,150.,200.,300./

      data rat /0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5, &
                0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0/

      data ecoll / &
         0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001, &
         0.001,0.001,0.001,0.001,0.001,0.003,0.003,0.003,0.004,0.005, &
         0.005,0.005,0.010,0.100,0.050,0.200,0.500,0.770,0.870,0.970, &
         0.007,0.007,0.007,0.008,0.009,0.010,0.010,0.070,0.400,0.430, &
         0.580,0.790,0.930,0.960,1.000,0.009,0.009,0.009,0.012,0.015, &
         0.010,0.020,0.280,0.600,0.640,0.750,0.910,0.970,0.980,1.000, &
         0.014,0.014,0.014,0.015,0.016,0.030,0.060,0.500,0.700,0.770, &
         0.840,0.950,0.970,1.000,1.000,0.017,0.017,0.017,0.020,0.022, &
         0.060,0.100,0.620,0.780,0.840,0.880,0.950,1.000,1.000,1.000, &
         0.030,0.030,0.024,0.022,0.032,0.062,0.200,0.680,0.830,0.870, &
         0.900,0.950,1.000,1.000,1.000,0.025,0.025,0.025,0.036,0.043, &
         0.130,0.270,0.740,0.860,0.890,0.920,1.000,1.000,1.000,1.000, &
         0.027,0.027,0.027,0.040,0.052,0.200,0.400,0.780,0.880,0.900, &
         0.940,1.000,1.000,1.000,1.000,0.030,0.030,0.030,0.047,0.064, &
         0.250,0.500,0.800,0.900,0.910,0.950,1.000,1.000,1.000,1.000, &
         0.040,0.040,0.033,0.037,0.068,0.240,0.550,0.800,0.900,0.910, &
         0.950,1.000,1.000,1.000,1.000,0.035,0.035,0.035,0.055,0.079, &
         0.290,0.580,0.800,0.900,0.910,0.950,1.000,1.000,1.000,1.000, &
         0.037,0.037,0.037,0.062,0.082,0.290,0.590,0.780,0.900,0.910, &
         0.950,1.000,1.000,1.000,1.000,0.037,0.037,0.037,0.060,0.080, &
         0.290,0.580,0.770,0.890,0.910,0.950,1.000,1.000,1.000,1.000, &
         0.037,0.037,0.037,0.041,0.075,0.250,0.540,0.760,0.880,0.920, &
         0.950,1.000,1.000,1.000,1.000,0.037,0.037,0.037,0.052,0.067, &
         0.250,0.510,0.770,0.880,0.930,0.970,1.000,1.000,1.000,1.000, &
         0.037,0.037,0.037,0.047,0.057,0.250,0.490,0.770,0.890,0.950, &
         1.000,1.000,1.000,1.000,1.000,0.036,0.036,0.036,0.042,0.048, &
         0.230,0.470,0.780,0.920,1.000,1.020,1.020,1.020,1.020,1.020, &
         0.040,0.040,0.035,0.033,0.040,0.112,0.450,0.790,1.010,1.030, &
         1.040,1.040,1.040,1.040,1.040,0.033,0.033,0.033,0.033,0.033, &
         0.119,0.470,0.950,1.300,1.700,2.300,2.300,2.300,2.300,2.300, &
         0.027,0.027,0.027,0.027,0.027,0.125,0.520,1.400,2.300,3.000, &
         4.000,4.000,4.000,4.000,4.000/

      ! two-dimensional linear interpolation of the collision efficiency
      ! query points for interpolation are drop_ratio and rbig

      ! selecting intervals
      do k=2,15
         if (rbig.le.r0(k).and.rbig.ge.r0(k-1)) then
            ir=k
         elseif (rbig.gt.r0(15)) then
            ir=16
         elseif (rbig.lt.r0(1)) then
            ir=1
         endif
      enddo

      do kk=2,21
         if (drop_ratio.le.rat(kk).and.drop_ratio.gt.rat(kk-1)) then
            iq=kk
         endif
      enddo

      ! interpolating
      if (ir.lt.16) then
         if (ir.ge.2) then
            p=(rbig-r0(ir-1))/(r0(ir)-r0(ir-1))
            q=(drop_ratio-rat(iq-1))/(rat(iq)-rat(iq-1))
            ec=(1.-p)*(1.-q)*ecoll(ir-1,iq-1)+  &
            p*(1.-q)*ecoll(ir,iq-1)+  &
            q*(1.-p)*ecoll(ir-1,iq)+  &
            (p*q*ecoll(ir,iq))
         else
            q=(drop_ratio-rat(iq-1))/(rat(iq)-rat(iq-1))
            ec=(1.-q)*ecoll(1,iq-1)+q*ecoll(1,iq)
         endif
      else
         q=(drop_ratio-rat(iq-1))/(rat(iq)-rat(iq-1))
         ek=(1.-q)*ecoll(15,iq-1)+q*ecoll(15,iq)
         ec=min(ek,1.d0)
      endif

      hall_effic = ec

   end function hall_effic

end module particles
