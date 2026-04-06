! ----------------------------------------------------------------------
! mod_mpi.f90 — MPI communication layer
!
! Contains all parallel communication routines:
!   - Domain decomposition range setup (range, set_range)
!   - Reduction/gather wrappers (mpi_sum_xy/z, mpi_gath_root, etc.)
!   - Array transposes for FFT and pressure solver
!   - Ghost-point exchange and boundary condition broadcasts
!
! The send_*/recv_* pack/unpack helpers are private; only the high-level
! transpose routines are exposed.
! ----------------------------------------------------------------------
module mod_mpi

  use mpi
  use pars,     only: nscl, nnx, nny, nnz, &
                      ncpu_s, ncpu_z, numprocs, myid, &
                      iss, ise, iys, iye, ixs, ixe, jxs, jxe, &
                      kxs, kxe, mxs, mxe, izs, ize, jys, jye, &
                      ix_s, ix_e, jx_s, jx_e, kx_s, kx_e, &
                      mx_s, mx_e, iy_s, iy_e, jy_s, jy_e, &
                      iz_s, iz_e, is_s, is_e, &
                      nprt, l_debug, qstar
  use fields,   only: u, v, w, t, e, pbc, pbc2
  use con_data, only: utau, u10

  implicit none
  private

  ! Public API
  public :: range, set_range
  public :: mpi_sum_xy, mpi_sum_z, mpi_sum_z_s
  public :: mpi_gath_root, mpi_send_root
  public :: xtoy_trans, ytox_trans
  public :: xtoz_trans, ztox_trans
  public :: exchange
  public :: bcast_pbc, bcast_lbc

contains

  ! --------------------------------------------------------------------
  ! range — IBM load-balanced range finder
  ! Splits [n1,n2] across nprocs processors; returns [ista,iend] for irank.
  ! --------------------------------------------------------------------
  subroutine range(n1, n2, nprocs, irank, ista, iend)
    integer, intent(in)  :: n1, n2, nprocs, irank
    integer, intent(out) :: ista, iend
    integer :: iwork1, iwork2

    iwork1 = (n2 - n1 + 1) / nprocs
    iwork2 = mod(n2 - n1 + 1, nprocs)
    ista   = irank*iwork1 + n1 + min(irank, iwork2)
    iend   = ista + iwork1 - 1
    if (iwork2 > irank) iend = iend + 1
  end subroutine range

  ! --------------------------------------------------------------------
  ! set_range — build x/y/z parallel decomposition index arrays
  !
  ! Populates the ix_s/ix_e, iy_s/iy_e, iz_s/iz_e, etc. arrays in pars
  ! and sets the per-process scalars (ixs, ixe, iys, iye, ...).
  ! --------------------------------------------------------------------
  subroutine set_range
    integer :: nn, mm, ii
    integer :: lx_s, lx_e, nx_s, nx_e, ly_s, ly_e, mz_s, mz_e
    integer :: ny_s, ny_e, nxy_s, nxy_e, l2x_s, l2x_e
    integer :: ncx

    ncx = nnx/2 + 1

    write(nprt, '(a,3i5)') ' set_range: nnx/nny/nnz =', nnx, nny, nnz

    ii = -1
    do nn = 0, ncpu_z - 1
      call range(1, nnx+2,  ncpu_z, nn, lx_s,  lx_e)
      call range(1, nnx,    ncpu_z, nn, nx_s,  nx_e)
      call range(1, nny,    ncpu_z, nn, ly_s,  ly_e)
      call range(1, nnz,    ncpu_z, nn, mz_s,  mz_e)
      do mm = 0, ncpu_s - 1
        call range(1, nny,  ncpu_s, mm, ny_s,   ny_e)
        call range(1, nnx,  ncpu_s, mm, nxy_s,  nxy_e)
        call range(1, ncx,  ncpu_s, mm, l2x_s,  l2x_e)
        ii = ii + 1

        ix_s(ii) = nxy_s
        ix_e(ii) = nxy_e
        jx_s(ii) = (l2x_s - 1)*2 + 1
        jx_e(ii) = l2x_e*2
        kx_s(ii) = lx_s
        kx_e(ii) = lx_e
        mx_s(ii) = nx_s
        mx_e(ii) = nx_e

        iy_s(ii) = ny_s
        iy_e(ii) = ny_e
        jy_s(ii) = ly_s
        jy_e(ii) = ly_e

        iz_s(ii) = mz_s
        iz_e(ii) = mz_e

        is_s(ii) = (ii/ncpu_s)*ncpu_s
        is_e(ii) = is_s(ii) + ncpu_s - 1
      end do
    end do

    iys = iy_s(myid);  iye = iy_e(myid)
    jys = jy_s(myid);  jye = jy_e(myid)
    ixs = ix_s(myid);  ixe = ix_e(myid)
    jxs = jx_s(myid);  jxe = jx_e(myid)
    kxs = kx_s(myid);  kxe = kx_e(myid)
    mxs = mx_s(myid);  mxe = mx_e(myid)
    izs = iz_s(myid);  ize = iz_e(myid)
    iss = is_s(myid);  ise = is_e(myid)

    if (l_debug) then
      write(nprt, '(a,i4,/,a,/,(7i6))') &
        ' myid = ', myid, &
        ' nn    ixs    ixe    jxs    jxe    kxs    kxe', &
        (nn, ix_s(nn), ix_e(nn), jx_s(nn), jx_e(nn), &
             kx_s(nn), kx_e(nn), nn=0, numprocs-1)

      write(nprt, '(a,i3,/,a,/,(9i6))') &
        ' myid = ', myid, &
        ' nn    iys    iye    jys    jye    izs    ize    iss    ise', &
        (nn, iy_s(nn), iy_e(nn), jy_s(nn), jy_e(nn), &
             iz_s(nn), iz_e(nn), is_s(nn), is_e(nn), nn=0, numprocs-1)
    end if
  end subroutine set_range

  ! --------------------------------------------------------------------
  ! mpi_sum_xy — horizontal x-y sum over processors [iss:ise]
  ! Overwrites f(1:nsend) with the sum. No-op for single processor.
  ! --------------------------------------------------------------------
  subroutine mpi_sum_xy(f, myid_l, iss_l, ise_l, nsend)
    integer, intent(in)    :: myid_l, iss_l, ise_l, nsend
    real,    intent(inout) :: f(nsend)
    real    :: work(nsend, iss_l:ise_l)
    integer :: istatus(mpi_status_size)
    integer :: i, j, ierr

    if (iss_l == ise_l) return

    do j = 1, nsend
      work(j, myid_l) = f(j)
      f(j) = 0.0
    end do
    do i = iss_l, ise_l
      if (i /= myid_l) then
        call mpi_sendrecv(work(1, myid_l), nsend, mpi_real8, i, 1, &
                          work(1, i),      nsend, mpi_real8, i, 1, &
                          mpi_comm_world, istatus, ierr)
      end if
    end do
    do i = iss_l, ise_l
      do j = 1, nsend
        f(j) = f(j) + work(j, i)
      end do
    end do
  end subroutine mpi_sum_xy

  ! --------------------------------------------------------------------
  ! mpi_sum_z — reduce 1-D vector across all z processors
  ! iall=1: all ranks get result; iall=0: only i_root gets result.
  ! --------------------------------------------------------------------
  subroutine mpi_sum_z(f, i_root, myid_l, nsend, iall)
    integer, intent(in)    :: i_root, myid_l, nsend, iall
    real,    intent(inout) :: f(nsend)
    real    :: recv_b(nsend)
    integer :: i, ierr

    if (iall /= 1) then
      call mpi_reduce(f(1), recv_b(1), nsend, mpi_real8, mpi_sum, &
                      i_root, mpi_comm_world, ierr)
      if (myid_l == i_root) f(1:nsend) = recv_b(1:nsend)
    else
      call mpi_allreduce(f(1), recv_b(1), nsend, mpi_real8, mpi_sum, &
                         mpi_comm_world, ierr)
      f(1:nsend) = recv_b(1:nsend)
    end if
  end subroutine mpi_sum_z

  ! --------------------------------------------------------------------
  ! mpi_sum_z_s — reduce 2-D (nsend x nscl) array across z processors
  ! --------------------------------------------------------------------
  subroutine mpi_sum_z_s(f, i_root, myid_l, nsend, nscl_l, iall)
    integer, intent(in)    :: i_root, myid_l, nsend, nscl_l, iall
    real,    intent(inout) :: f(nsend, nscl_l)
    real    :: recv_b(nsend, nscl_l)
    integer :: i, iscl, ierr

    if (iall /= 1) then
      call mpi_reduce(f(1,1), recv_b(1,1), nsend*nscl_l, mpi_real8, &
                      mpi_sum, i_root, mpi_comm_world, ierr)
      if (myid_l == i_root) f = recv_b
    else
      call mpi_allreduce(f(1,1), recv_b(1,1), nsend*nscl_l, mpi_real8, &
                         mpi_sum, mpi_comm_world, ierr)
      f = recv_b
    end if
  end subroutine mpi_sum_z_s

  ! --------------------------------------------------------------------
  ! mpi_gath_root — gather vertical slabs onto root processors
  ! --------------------------------------------------------------------
  subroutine mpi_gath_root(fs, fr, iz_s_l, iz_e_l, izs_l, ize_l, &
                            nz, myid_l, np, ns)
    integer, intent(in)  :: iz_s_l(0:np-1), iz_e_l(0:np-1)
    integer, intent(in)  :: izs_l, ize_l, nz, myid_l, np, ns
    real,    intent(in)  :: fs(izs_l:ize_l)
    real,    intent(out) :: fr(1:nz)
    integer :: istatus(mpi_status_size)
    integer :: irow_r, l, ind, num, ierr

    if (np == 1) return

    irow_r = mod(myid_l, ns)
    if (myid_l > ns) then
      call mpi_send(fs(izs_l), ize_l+1-izs_l, mpi_real8, &
                    irow_r, 1, mpi_comm_world, ierr)
    else
      do l = irow_r+ns, np-1, ns
        ind = iz_s_l(l) + 1
        num = iz_e_l(l) + 1 - iz_s_l(l)
        call mpi_recv(fr(ind), num, mpi_real8, l, 1, &
                      mpi_comm_world, istatus, ierr)
      end do
    end if
  end subroutine mpi_gath_root

  ! --------------------------------------------------------------------
  ! mpi_send_root — broadcast root results to non-root processors
  ! --------------------------------------------------------------------
  subroutine mpi_send_root(fs, num, myid_l, np, ns)
    integer, intent(in)    :: num, myid_l, np, ns
    real,    intent(inout) :: fs(num)
    integer :: istatus(mpi_status_size)
    integer :: irow_r, l, ierr

    if (np == 1) return

    irow_r = mod(myid_l, ns)
    if (myid_l >= ns) then
      call mpi_recv(fs(1), num, mpi_real8, irow_r, 1, &
                    mpi_comm_world, istatus, ierr)
    else
      do l = irow_r+ns, np-1, ns
        call mpi_send(fs(1), num, mpi_real8, l, 1, &
                      mpi_comm_world, ierr)
      end do
    end if
  end subroutine mpi_send_root

  ! --------------------------------------------------------------------
  ! xtoy_trans — transpose f(nx, iys:iye, iz1:iz2)
  !                     -> g(ny, ixs:ixe, iz1:iz2)
  ! --------------------------------------------------------------------
  subroutine xtoy_trans(f, g, nx, ny, ixs_l, ixe_l, ix_s_l, ix_e_l, &
                         iys_l, iye_l, iy_s_l, iy_e_l, iz1, iz2, &
                         myid_l, ncpu_s_l, np)
    integer, intent(in) :: nx, ny, ixs_l, ixe_l, iys_l, iye_l
    integer, intent(in) :: iz1, iz2, myid_l, ncpu_s_l, np
    integer, intent(in) :: ix_s_l(0:np-1), ix_e_l(0:np-1)
    integer, intent(in) :: iy_s_l(0:np-1), iy_e_l(0:np-1)
    real, intent(inout) :: f(nx, iys_l:iye_l, iz1:iz2)
    real, intent(out)   :: g(ny, ixs_l:ixe_l, iz1:iz2)
    real    :: ft(nx*(iye_l+1-iys_l)*(iz2-iz1+1))
    real    :: gt(ny*(ixe_l+1-ixs_l)*(iz2-iz1+1))
    integer :: istatus(mpi_status_size)
    integer :: i, islab, iss_l, ise_l, jk, ik, nsend, nrecv, ierr

    jk = (iye_l - iys_l + 1)*(iz2 - iz1 + 1)
    ik = (ixe_l - ixs_l + 1)*(iz2 - iz1 + 1)

    islab = myid_l / ncpu_s_l
    iss_l = islab * ncpu_s_l
    ise_l = iss_l + ncpu_s_l - 1

    do i = iss_l, ise_l
      nsend = (ix_e_l(i) - ix_s_l(i) + 1)*jk
      nrecv = (iy_e_l(i) - iy_s_l(i) + 1)*ik
      if (i == myid_l) then
        call send_xtoy(f, gt(1), nx, ix_s_l(i), ix_e_l(i), &
                       iy_s_l(myid_l), iy_e_l(myid_l), iz1, iz2)
      else
        call send_xtoy(f, ft(1), nx, ix_s_l(i), ix_e_l(i), &
                       iy_s_l(myid_l), iy_e_l(myid_l), iz1, iz2)
        call mpi_sendrecv(ft(1), nsend, mpi_real8, i, 1, &
                          gt(1), nrecv, mpi_real8, i, 1, &
                          mpi_comm_world, istatus, ierr)
      end if
      call recv_xtoy(g, gt(1), ny, ix_s_l(myid_l), ix_e_l(myid_l), &
                     iy_s_l(i), iy_e_l(i), iz1, iz2)
    end do
  end subroutine xtoy_trans

  ! --------------------------------------------------------------------
  ! ytox_trans — inverse of xtoy_trans
  !   g(ny, ixs:ixe, iz1:iz2) -> f(nx, iys:iye, iz1:iz2)
  ! --------------------------------------------------------------------
  subroutine ytox_trans(g, f, nx, ny, ixs_l, ixe_l, ix_s_l, ix_e_l, &
                         iys_l, iye_l, iy_s_l, iy_e_l, iz1, iz2, &
                         myid_l, ncpu_s_l, np)
    integer, intent(in) :: nx, ny, ixs_l, ixe_l, iys_l, iye_l
    integer, intent(in) :: iz1, iz2, myid_l, ncpu_s_l, np
    integer, intent(in) :: ix_s_l(0:np-1), ix_e_l(0:np-1)
    integer, intent(in) :: iy_s_l(0:np-1), iy_e_l(0:np-1)
    real, intent(in)    :: g(ny, ixs_l:ixe_l, iz1:iz2)
    real, intent(out)   :: f(nx, iys_l:iye_l, iz1:iz2)
    real    :: ft(nx*(iye_l+1-iys_l)*(iz2-iz1+1))
    real    :: gt(ny*(ixe_l+1-ixs_l)*(iz2-iz1+1))
    integer :: istatus(mpi_status_size)
    integer :: i, islab, iss_l, ise_l, jk, ik, nsend, nrecv, ierr

    jk = (iye_l - iys_l + 1)*(iz2 - iz1 + 1)
    ik = (ixe_l - ixs_l + 1)*(iz2 - iz1 + 1)

    islab = myid_l / ncpu_s_l
    iss_l = islab * ncpu_s_l
    ise_l = iss_l + ncpu_s_l - 1

    do i = iss_l, ise_l
      nsend = (iy_e_l(i) - iy_s_l(i) + 1)*ik
      nrecv = (ix_e_l(i) - ix_s_l(i) + 1)*jk
      if (i == myid_l) then
        call send_ytox(g, ft(1), ny, ix_s_l(myid_l), ix_e_l(myid_l), &
                       iy_s_l(i), iy_e_l(i), iz1, iz2)
      else
        call send_ytox(g, gt(1), ny, ix_s_l(myid_l), ix_e_l(myid_l), &
                       iy_s_l(i), iy_e_l(i), iz1, iz2)
        call mpi_sendrecv(gt(1), nsend, mpi_real8, i, 1, &
                          ft(1), nrecv, mpi_real8, i, 1, &
                          mpi_comm_world, istatus, ierr)
      end if
      call recv_ytox(f, ft(1), nx, ix_s_l(i), ix_e_l(i), &
                     iy_s_l(myid_l), iy_e_l(myid_l), iz1, iz2)
    end do
  end subroutine ytox_trans

  ! --------------------------------------------------------------------
  ! xtoz_trans — transpose f(nx, iys:iye, izs-1:ize+1)
  !                     -> g(0:nz+1, iys:iye, ixs:ixe)
  ! --------------------------------------------------------------------
  subroutine xtoz_trans(f, g, nx, nz, ixs_l, ixe_l, ix_s_l, ix_e_l, &
                         iys_l, iye_l, izs_l, ize_l, iz_s_l, iz_e_l, &
                         myid_l, ncpu_s_l, np)
    integer, intent(in) :: nx, nz, ixs_l, ixe_l, iys_l, iye_l
    integer, intent(in) :: izs_l, ize_l, myid_l, ncpu_s_l, np
    integer, intent(in) :: ix_s_l(0:np-1), ix_e_l(0:np-1)
    integer, intent(in) :: iz_s_l(0:np-1), iz_e_l(0:np-1)
    real, intent(inout) :: f(nx, iys_l:iye_l, izs_l-1:ize_l+1)
    real, intent(out)   :: g(0:nz+1, iys_l:iye_l, ixs_l:ixe_l)
    real    :: ft(nx*(iye_l+1-iys_l)*(ize_l-izs_l+1))
    real    :: gt(nz*(ixe_l+1-ixs_l)*(iye_l-iys_l+1))
    integer :: istatus(mpi_status_size)
    integer :: i, iss_l, jk, ij, nsend, nrecv, ierr

    jk = (ize_l - izs_l + 1)*(iye_l - iys_l + 1)
    ij = (ixe_l - ixs_l + 1)*(iye_l - iys_l + 1)

    iss_l = myid_l - (myid_l/ncpu_s_l)*ncpu_s_l

    do i = iss_l, np-1, ncpu_s_l
      nsend = (ix_e_l(i) - ix_s_l(i) + 1)*jk
      nrecv = (iz_e_l(i) - iz_s_l(i) + 1)*ij
      if (i == myid_l) then
        call send_xtoz(f, gt(1), nx, ix_s_l(i), ix_e_l(i), &
                       iys_l, iye_l, iz_s_l(myid_l), iz_e_l(myid_l))
      else
        call send_xtoz(f, ft(1), nx, ix_s_l(i), ix_e_l(i), &
                       iys_l, iye_l, iz_s_l(myid_l), iz_e_l(myid_l))
        call mpi_sendrecv(ft(1), nsend, mpi_real8, i, 1, &
                          gt(1), nrecv, mpi_real8, i, 1, &
                          mpi_comm_world, istatus, ierr)
      end if
      call recv_xtoz(g, gt(1), nz, ix_s_l(myid_l), ix_e_l(myid_l), &
                     iys_l, iye_l, iz_s_l(i), iz_e_l(i))
    end do
  end subroutine xtoz_trans

  ! --------------------------------------------------------------------
  ! ztox_trans — inverse of xtoz_trans
  !   g(0:nz+1, iys:iye, ixs:ixe) -> f(nx, iys:iye, izs-1:ize+1)
  ! --------------------------------------------------------------------
  subroutine ztox_trans(g, f, nx, nz, ixs_l, ixe_l, ix_s_l, ix_e_l, &
                         iys_l, iye_l, izs_l, ize_l, iz_s_l, iz_e_l, &
                         myid_l, ncpu_s_l, np)
    integer, intent(in) :: nx, nz, ixs_l, ixe_l, iys_l, iye_l
    integer, intent(in) :: izs_l, ize_l, myid_l, ncpu_s_l, np
    integer, intent(in) :: ix_s_l(0:np-1), ix_e_l(0:np-1)
    integer, intent(in) :: iz_s_l(0:np-1), iz_e_l(0:np-1)
    real, intent(out)   :: f(nx, iys_l:iye_l, izs_l-1:ize_l+1)
    real, intent(inout) :: g(0:nz+1, iys_l:iye_l, ixs_l:ixe_l)
    real    :: ft(nx*(iye_l+1-iys_l)*(ize_l-izs_l+3))
    real    :: gt((nz+3)*(iye_l+1-iys_l)*(ixe_l-ixs_l+1))
    integer :: istatus(mpi_status_size)
    integer :: i, iss_l, jk, ij, nsend, nrecv, ierr

    jk = (ize_l - izs_l + 3)*(iye_l - iys_l + 1)
    ij = (ixe_l - ixs_l + 1)*(iye_l - iys_l + 1)

    iss_l = myid_l - (myid_l/ncpu_s_l)*ncpu_s_l

    do i = iss_l, np-1, ncpu_s_l
      nsend = (iz_e_l(i) - iz_s_l(i) + 3)*ij
      nrecv = (ix_e_l(i) - ix_s_l(i) + 1)*jk
      if (i == myid_l) then
        call send_ztox(g, ft(1), nz, ix_s_l(myid_l), ix_e_l(myid_l), &
                       iys_l, iye_l, iz_s_l(i), iz_e_l(i))
      else
        call send_ztox(g, gt(1), nz, ix_s_l(myid_l), ix_e_l(myid_l), &
                       iys_l, iye_l, iz_s_l(i), iz_e_l(i))
        call mpi_sendrecv(gt(1), nsend, mpi_real8, i, 1, &
                          ft(1), nrecv, mpi_real8, i, 1, &
                          mpi_comm_world, istatus, ierr)
      end if
      call recv_ztox(f, ft(1), nx, ix_s_l(i), ix_e_l(i), &
                     iys_l, iye_l, iz_s_l(myid_l), iz_e_l(myid_l))
    end do
  end subroutine ztox_trans

  ! --------------------------------------------------------------------
  ! exchange — swap ghost points in the vertical (z) direction
  ! --------------------------------------------------------------------
  subroutine exchange
    real    :: fs(nnx, iys:iye, 4+nscl)
    real    :: fr(nnx, iys:iye, 4+nscl)
    integer :: istatus(mpi_status_size)
    integer :: nb, nt, nsend, nrecv
    integer :: ix, iy, iscl, jloc, izm1, izp1, izm2, izp2, ierr

    nb = myid - ncpu_s
    nt = myid + ncpu_s
    if (iss == 0)            nb = mpi_proc_null
    if (ise == numprocs-1)   nt = mpi_proc_null

    nsend = nnx*(iye + 1 - iys)*(4 + nscl)
    nrecv = nsend

    ! Send top of myid; receive bottom from myid-ncpu_s
    do iy = iys, iye
      do ix = 1, nnx
        fs(ix,iy,1) = u(ix,iy,ize)
        fs(ix,iy,2) = v(ix,iy,ize)
        fs(ix,iy,3) = w(ix,iy,ize)
        fs(ix,iy,4) = e(ix,iy,ize)
      end do
    end do
    do iscl = 1, nscl
      jloc = 4 + iscl
      do iy = iys, iye
        do ix = 1, nnx
          fs(ix,iy,jloc) = t(ix,iy,iscl,ize)
        end do
      end do
    end do

    call mpi_sendrecv(fs(1,iys,1), nsend, mpi_real8, nt, 0, &
                      fr(1,iys,1), nrecv, mpi_real8, nb, 0, &
                      mpi_comm_world, istatus, ierr)

    if (iss /= 0) then
      izm1 = izs - 1
      do iy = iys, iye
        do ix = 1, nnx
          u(ix,iy,izm1) = fr(ix,iy,1)
          v(ix,iy,izm1) = fr(ix,iy,2)
          w(ix,iy,izm1) = fr(ix,iy,3)
          e(ix,iy,izm1) = fr(ix,iy,4)
        end do
      end do
      do iscl = 1, nscl
        jloc = 4 + iscl
        do iy = iys, iye
          do ix = 1, nnx
            t(ix,iy,iscl,izm1) = fr(ix,iy,jloc)
          end do
        end do
      end do
    end if

    ! Send bottom of myid; receive top from myid+ncpu_s
    do iy = iys, iye
      do ix = 1, nnx
        fs(ix,iy,1) = u(ix,iy,izs)
        fs(ix,iy,2) = v(ix,iy,izs)
        fs(ix,iy,3) = w(ix,iy,izs)
        fs(ix,iy,4) = e(ix,iy,izs)
      end do
    end do
    do iscl = 1, nscl
      jloc = 4 + iscl
      do iy = iys, iye
        do ix = 1, nnx
          fs(ix,iy,jloc) = t(ix,iy,iscl,izs)
        end do
      end do
    end do

    call mpi_sendrecv(fs(1,iys,1), nsend, mpi_real8, nb, 1, &
                      fr(1,iys,1), nrecv, mpi_real8, nt, 1, &
                      mpi_comm_world, istatus, ierr)

    if (ise /= numprocs-1) then
      izp1 = ize + 1
      do iy = iys, iye
        do ix = 1, nnx
          u(ix,iy,izp1) = fr(ix,iy,1)
          v(ix,iy,izp1) = fr(ix,iy,2)
          w(ix,iy,izp1) = fr(ix,iy,3)
          e(ix,iy,izp1) = fr(ix,iy,4)
        end do
      end do
      do iscl = 1, nscl
        jloc = 4 + iscl
        do iy = iys, iye
          do ix = 1, nnx
            t(ix,iy,iscl,izp1) = fr(ix,iy,jloc)
          end do
        end do
      end do
    end if

    ! Exchange extra scalar ghost points (two-level stencil)
    nsend = nnx*(iye + 1 - iys)*nscl
    nrecv = nsend

    izm1 = ize - 1
    do iscl = 1, nscl
      do iy = iys, iye
        do ix = 1, nnx
          fs(ix,iy,iscl) = t(ix,iy,iscl,izm1)
        end do
      end do
    end do

    call mpi_sendrecv(fs(1,iys,1), nsend, mpi_real8, nt, 0, &
                      fr(1,iys,1), nrecv, mpi_real8, nb, 0, &
                      mpi_comm_world, istatus, ierr)

    if (iss /= 0) then
      izm2 = izs - 2
      do iscl = 1, nscl
        do iy = iys, iye
          do ix = 1, nnx
            t(ix,iy,iscl,izm2) = fr(ix,iy,iscl)
          end do
        end do
      end do
    end if

    izp1 = izs + 1
    do iscl = 1, nscl
      do iy = iys, iye
        do ix = 1, nnx
          fs(ix,iy,iscl) = t(ix,iy,iscl,izp1)
        end do
      end do
    end do

    call mpi_sendrecv(fs(1,iys,1), nsend, mpi_real8, nb, 1, &
                      fr(1,iys,1), nrecv, mpi_real8, nt, 1, &
                      mpi_comm_world, istatus, ierr)

    if (ise /= numprocs-1) then
      izp2 = ize + 2
      do iscl = 1, nscl
        do iy = iys, iye
          do ix = 1, nnx
            t(ix,iy,iscl,izp2) = fr(ix,iy,iscl)
          end do
        end do
      end do
    end if
  end subroutine exchange

  ! --------------------------------------------------------------------
  ! bcast_pbc — send upper pressure BCs from top row to all processors
  ! --------------------------------------------------------------------
  subroutine bcast_pbc
    integer :: istatus(mpi_status_size)
    integer :: irow_r, irow_t, num, l, ierr

    if (numprocs == 1) return

    irow_r = mod(myid, ncpu_s)
    irow_t = is_s(numprocs-1) + irow_r
    num    = nnx*(iye + 1 - iys)

    if (iss /= is_s(numprocs-1)) then
      call mpi_recv(pbc(1,iys,1),  num, mpi_real8, irow_t, 1, &
                    mpi_comm_world, istatus, ierr)
      call mpi_recv(pbc2(1,iys,1), num, mpi_real8, irow_t, 1, &
                    mpi_comm_world, istatus, ierr)
    else
      do l = irow_r, irow_t-ncpu_s, ncpu_s
        call mpi_send(pbc(1,iys,1),  num, mpi_real8, l, 1, &
                      mpi_comm_world, ierr)
        call mpi_send(pbc2(1,iys,1), num, mpi_real8, l, 1, &
                      mpi_comm_world, ierr)
      end do
    end if
  end subroutine bcast_pbc

  ! --------------------------------------------------------------------
  ! bcast_lbc — synchronise lower boundary condition scalars across
  !             horizontal processors then broadcast to all
  ! --------------------------------------------------------------------
  subroutine bcast_lbc
    real    :: upars(4)
    integer :: ierr

    upars(1) = utau
    upars(2) = u10
    upars(3) = qstar(1)
    upars(4) = qstar(2)

    call mpi_sum_xy(upars, myid, iss, ise, 4)

    u10 = upars(2) / ncpu_s

    call mpi_bcast(u10,    1, mpi_real8, 0, mpi_comm_world, ierr)
    call mpi_bcast(qstar,  2, mpi_real8, 0, mpi_comm_world, ierr)
  end subroutine bcast_lbc

  ! ====================================================================
  ! Private pack/unpack helpers for transpose routines
  ! ====================================================================

  subroutine send_xtoy(f, ft, nx, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l)
    integer, intent(in) :: nx, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l
    real, intent(in)  :: f(nx, iys_l:iye_l, izs_l:ize_l)
    real, intent(out) :: ft(ixs_l:ixe_l, iys_l:iye_l, izs_l:ize_l)
    integer :: i, j, k
    do k = izs_l, ize_l
      do j = iys_l, iye_l
        do i = ixs_l, ixe_l
          ft(i,j,k) = f(i,j,k)
        end do
      end do
    end do
  end subroutine send_xtoy

  subroutine recv_xtoy(g, gt, ny, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l)
    integer, intent(in) :: ny, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l
    real, intent(out) :: g(ny, ixs_l:ixe_l, izs_l:ize_l)
    real, intent(in)  :: gt(ixs_l:ixe_l, iys_l:iye_l, izs_l:ize_l)
    integer :: i, j, k
    do k = izs_l, ize_l
      do j = iys_l, iye_l
        do i = ixs_l, ixe_l
          g(j,i,k) = gt(i,j,k)
        end do
      end do
    end do
  end subroutine recv_xtoy

  subroutine send_ytox(g, gt, ny, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l)
    integer, intent(in) :: ny, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l
    real, intent(in)  :: g(ny, ixs_l:ixe_l, izs_l:ize_l)
    real, intent(out) :: gt(iys_l:iye_l, ixs_l:ixe_l, izs_l:ize_l)
    integer :: i, j, k
    do k = izs_l, ize_l
      do i = ixs_l, ixe_l
        do j = iys_l, iye_l
          gt(j,i,k) = g(j,i,k)
        end do
      end do
    end do
  end subroutine send_ytox

  subroutine recv_ytox(f, ft, nx, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l)
    integer, intent(in) :: nx, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l
    real, intent(out) :: f(nx, iys_l:iye_l, izs_l:ize_l)
    real, intent(in)  :: ft(iys_l:iye_l, ixs_l:ixe_l, izs_l:ize_l)
    integer :: i, j, k
    do k = izs_l, ize_l
      do i = ixs_l, ixe_l
        do j = iys_l, iye_l
          f(i,j,k) = ft(j,i,k)
        end do
      end do
    end do
  end subroutine recv_ytox

  subroutine send_xtoz(f, ft, nx, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l)
    integer, intent(in) :: nx, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l
    real, intent(in)  :: f(nx, iys_l:iye_l, izs_l-1:ize_l+1)
    real, intent(out) :: ft(ixs_l:ixe_l, iys_l:iye_l, izs_l:ize_l)
    integer :: i, j, k
    do k = izs_l, ize_l
      do j = iys_l, iye_l
        do i = ixs_l, ixe_l
          ft(i,j,k) = f(i,j,k)
        end do
      end do
    end do
  end subroutine send_xtoz

  subroutine recv_xtoz(g, gt, nz, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l)
    integer, intent(in) :: nz, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l
    real, intent(out) :: g(0:nz+1, iys_l:iye_l, ixs_l:ixe_l)
    real, intent(in)  :: gt(ixs_l:ixe_l, iys_l:iye_l, izs_l:ize_l)
    integer :: i, j, k
    do k = izs_l, ize_l
      do j = iys_l, iye_l
        do i = ixs_l, ixe_l
          g(k,j,i) = gt(i,j,k)
        end do
      end do
    end do
  end subroutine recv_xtoz

  subroutine send_ztox(g, gt, nz, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l)
    integer, intent(in) :: nz, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l
    real, intent(in)  :: g(0:nz+1, iys_l:iye_l, ixs_l:ixe_l)
    real, intent(out) :: gt(izs_l-1:ize_l+1, iys_l:iye_l, ixs_l:ixe_l)
    integer :: i, j, k
    do j = iys_l, iye_l
      do i = ixs_l, ixe_l
        do k = izs_l-1, ize_l+1
          gt(k,j,i) = g(k,j,i)
        end do
      end do
    end do
  end subroutine send_ztox

  subroutine recv_ztox(f, ft, nx, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l)
    integer, intent(in) :: nx, ixs_l, ixe_l, iys_l, iye_l, izs_l, ize_l
    real, intent(out) :: f(nx, iys_l:iye_l, izs_l-1:ize_l+1)
    real, intent(in)  :: ft(izs_l-1:ize_l+1, iys_l:iye_l, ixs_l:ixe_l)
    integer :: i, j, k
    do i = ixs_l, ixe_l
      do j = iys_l, iye_l
        do k = izs_l-1, ize_l+1
          f(i,j,k) = ft(k,j,i)
        end do
      end do
    end do
  end subroutine recv_ztox

end module mod_mpi
