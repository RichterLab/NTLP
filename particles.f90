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
      real :: radavg,radmin,radmax,tempmin,tempmax,qmin,qmax
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

      
!REMEMBER: IF ADDING ANYTHING, MUST UPDATE MPI DATATYPE!
      type :: particle
      integer :: pidx,procidx,nbr_pidx,nbr_procidx
      real :: vp(3),xp(3),uf(3),xrhs(3),vrhs(3),Tp,Tprhs_s
      real :: Tprhs_L,Tf,radius,radrhs,qinf,qstar,dist
      real :: res,m_s,Os,rc,actres
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
