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

  !--- SFS velocity calculation ---------
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
  integer :: tot_reintro=0

  real :: Rep_avg,part_grav(3)
  real :: radavg,radmin,radmax,radmsqr,tempmin,tempmax,qmin,qmax
  real :: vp_init(3),Tp_init,radius_init
  real :: pdf_factor,pdf_prob
  integer*8 :: mult_init,mult_factor,mult_a,mult_c

  real :: avgres=0,tavgres=0

  integer, parameter :: histbins = 512
  real :: hist_rad(histbins+2)
  real :: bins_rad(histbins+2)

  real :: hist_res(histbins+2)
  real :: bins_res(histbins+2)

  real :: hist_actres(histbins+2)
  real :: bins_actres(histbins+2)

  real :: hist_numact(histbins+2)
  real :: bins_numact(histbins+2)

  !REMEMBER: IF ADDING ANYTHING, MUST UPDATE MPI DATATYPE!
  type :: particle
    integer :: pidx,procidx,nbr_pidx,nbr_procidx
    real :: vp(3),xp(3),uf(3),xrhs(3),vrhs(3),Tp,Tprhs_s
    real :: Tprhs_L,Tf,radius,radrhs,qinf,qstar,dist
    real :: res,m_s,Os,rc,actres,numact
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

  subroutine particle_reintro(it)
      use pars
      use con_data
      use con_stats
      use fields
      implicit none
      include 'mpif.h'

      integer :: it
      integer :: ierr,randproc,np,my_reintro
      real :: xp_init(3),ran2,Os,m_s

      my_reintro = nprime*(1./60.)*(10.**6.)*dt*4/numprocs*10.0 !4m^3 (vol chamber)

      tot_reintro = 0

      if (mod(it, 10)==0) then

      tot_reintro = my_reintro*numprocs

      if (myid==0) write(*,*) 'time,tot_reintro:',time,tot_reintro

      do np=1,my_reintro

      !Proc 0 gets a random proc ID, broadcasts it out:
      !if (myid==0) randproc = floor(ran2(iseed)*numprocs)
      !call mpi_bcast(randproc,1,mpi_integer,0,mpi_comm_world,ierr)


         xp_init(1) = ran2(iseed)*(xmax-xmin) + xmin
         xp_init(2) = ran2(iseed)*(ymax-ymin) + ymin
         xp_init(3) = ran2(iseed)*zl/2.0

         Os = 1.0
         m_s = radius_init**3*pi2*2.0/3.0*rhow*Sal  !Using the salinity specified in params.in
         

         call create_particle(xp_init,vp_init, &
              Tp_init,m_s,Os,mult_init,radius_init,2,myid) 

      end do
      end if


  end subroutine particle_reintro

  subroutine create_particle(xp,vp,Tp,m_s,Os,mult,rad_init,idx,procidx)
      use pars
      implicit none

      real :: xp(3),vp(3),Tp,qinfp,rad_init,pi,m_s,Os
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
      part%Os = Os
      part%dist = 0.0
      part%u_sub(1:3) = 0.0
      part%sigm_s = 0.0
      part%numact = 0.0

      
  end subroutine create_particle

  function mod_Magnus(T)
    implicit none

    !Take in T in Kelvin and return saturation vapor pressure using function of Alduchov and Eskridge, 1996
    real,intent(in) :: T
    real :: mod_Magnus

    mod_Magnus = 610.94 *exp((17.6257*(T-273.15))/(243.04+(T-273.15)))


  end function mod_Magnus

  function crit_radius(m_s,Os,Tf)
    use pars
    use con_data
    implicit none

    integer :: i,maxidx(1)
    integer, parameter :: N=1000
    real :: m_s,Os,Tf
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
      secterm =Ion*Os*m_s*(Mw/Ms)/(rhow*pi2*2.0/3.0*radval(i)**3)

      SS(i) = exp( 2*Mw*Gam/Ru/rhow/radval(i)/Tf - Ion*Os*m_s*(Mw/Ms)/(rhow*pi2*2.0/3.0*radval(i)**3))

    end do

    maxidx = maxloc(SS)

    crit_radius = radval(maxidx(1))


  end function crit_radius


end module particles
