! ----------------------------------------------------------------------
! mod_fft — spectral derivative and transform routines
!
! Public interface
!   xderivp      : x-derivatives via real FFT (FFTPACK rfftf/rfftb)
!   yderiv       : y-derivatives via real FFT, serial (FFTPACK rfftf/rfftb)
!   fft2d_mpi    : 2-D forward/inverse FFT with MPI slab transposes
!   yd_mpi       : y-derivatives via real FFT with MPI transposes
!   dealias      : 2/3-rule spectral dealiasing of all prognostic fields
!   sharp        : sharp spectral cutoff filter on a 2-D FFT coefficient array
!
! Notes
!   FFTPACK routines (rfftf, rfftb, cfftf, cfftb) are called as external
!   F77 routines from fft.f; they carry implicit interfaces here.
!   Dummy arguments whose names would shadow host-associated module variables
!   (trigx, trigc from fftwk) carry an _l suffix.
! ----------------------------------------------------------------------
module mod_fft
  use mod_mpi,  only: xtoy_trans, ytox_trans
  use pars,     only: nscl,   &
                      nnx, nnxp2, nny, nnyp2, &
                      jxs, jxe, jx_s, jx_e,  &
                      iys, iye, iy_s, iy_e,  &
                      izs, ize,               &
                      myid, ncpu_s, numprocs
  use fields,   only: u, v, w, t, e
  use fftwk,    only: trigx, trigc
  implicit none
  private
  public :: xderivp, yderiv, fft2d_mpi, yd_mpi, dealias, sharp

contains

  ! --------------------------------------------------------------------
  subroutine xderivp(ax, trigx_l, xk, nnx_l, iys_l, iye_l)
  ! Get multiple x-derivatives using FFTPACK.
  ! Storage order: a0, (a1,b1), (a2,b2), ..., an
  ! Assumes xk(1) = 0 and wavenumbers normalized by number of points.
  ! --------------------------------------------------------------------
    real,    intent(inout) :: ax(nnx_l, iys_l:iye_l)
    real,    intent(in)    :: trigx_l(2*nnx_l+15), xk(nnx_l)
    integer, intent(in)    :: nnx_l, iys_l, iye_l
    integer :: ix, iy, ii
    real    :: temp

    do iy = iys_l, iye_l
      call rfftf(nnx_l, ax(1,iy), trigx_l)
      ii         = 1
      ax(1,iy)   = 0.0
      ax(nnx_l,iy) = 0.0
      do ix = 2, nnx_l-1, 2
        ii          = ii + 1
        temp        = ax(ix,iy)
        ax(ix,iy)   = -xk(ii)*ax(ix+1,iy)
        ax(ix+1,iy) =  xk(ii)*temp
      end do
      call rfftb(nnx_l, ax(1,iy), trigx_l)
    end do
  end subroutine xderivp

  ! --------------------------------------------------------------------
  subroutine yderiv(ay, trigy_l, yk, nnx_l, nny_l)
  ! Get multiple y-derivatives using FFTPACK (serial version).
  ! Storage order: a0, (a1,b1), (a2,b2), ...
  ! Assumes yk(1) = 0 and wavenumbers normalized by number of points.
  ! --------------------------------------------------------------------
    real,    intent(inout) :: ay(nnx_l, nny_l)
    real,    intent(in)    :: trigy_l(2*nny_l+15), yk(nny_l)
    integer, intent(in)    :: nnx_l, nny_l
    real    :: a_trans(nny_l)
    integer :: ix, iy, ii
    real    :: temp

    do ix = 1, nnx_l
      do iy = 1, nny_l
        a_trans(iy) = ay(ix,iy)
      end do
      call rfftf(nny_l, a_trans(1), trigy_l)
      ii            = 1
      a_trans(1)    = 0.0
      a_trans(nny_l) = 0.0
      do iy = 2, nny_l-1, 2
        ii              = ii + 1
        temp            = a_trans(iy)
        a_trans(iy)     = -yk(ii)*a_trans(iy+1)
        a_trans(iy+1)   =  yk(ii)*temp
      end do
      call rfftb(nny_l, a_trans(1), trigy_l)
      do iy = 1, nny_l
        ay(ix,iy) = a_trans(iy)
      end do
    end do
  end subroutine yderiv

  ! --------------------------------------------------------------------
  subroutine fft2d_mpi(ax, at, trigx_l, trigc_l, nx, ny, &
                        jxs_l, jxe_l, jx_s_l, jx_e_l,   &
                        iys_l, iye_l, iy_s_l, iy_e_l,   &
                        iz1, iz2, myid_l, ncpu, np, isgn)
  ! 2-D FFT using FFTPACK + MPI slab transposes.
  !
  ! isgn = -1  forward transform; result in ax(nx+2, iys:iye, iz1:iz2)
  ! isgn = -2  forward transform; result in at(ny,   jxs:jxe, iz1:iz2)
  ! isgn =  1  inverse transform; input  in ax(nx+2, iys:iye, iz1:iz2)
  ! isgn =  2  inverse transform; input  in at(ny,   jxs:jxe, iz1:iz2)
  ! --------------------------------------------------------------------
    real,    intent(inout) :: ax(nx+2, iys_l:iye_l, iz1:iz2)
    real,    intent(inout) :: at(ny,   jxs_l:jxe_l, iz1:iz2)
    real,    intent(in)    :: trigx_l(2*nx+15), trigc_l(4*ny+15)
    integer, intent(in)    :: nx, ny
    integer, intent(in)    :: jxs_l, jxe_l
    integer, intent(in)    :: jx_s_l(0:np-1), jx_e_l(0:np-1)
    integer, intent(in)    :: iys_l, iye_l
    integer, intent(in)    :: iy_s_l(0:np-1), iy_e_l(0:np-1)
    integer, intent(in)    :: iz1, iz2, myid_l, ncpu, np, isgn
    real    :: a2d(2,ny), a_wrk(nx)
    integer :: ix, iy, iz, nxp2
    real    :: fn

    nxp2 = nx + 2

    if (isgn < 0) then
      fn = 1.0 / (real(nx) * real(ny))

      ! 1-D FFT in x over [iys,iye] for all z
      do iz = iz1, iz2
        do iy = iys_l, iye_l
          do ix = 1, nx
            a_wrk(ix) = ax(ix,iy,iz) * fn
          end do
          call rfftf(nx, a_wrk(1), trigx_l(1))
          ax(1,iy,iz)    = a_wrk(1)
          ax(2,iy,iz)    = 0.0
          do ix = 2, nx
            ax(ix+1,iy,iz) = a_wrk(ix)
          end do
          ax(nx+2,iy,iz) = 0.0
        end do
      end do
      call xtoy_trans(ax, at, nxp2, ny, jxs_l, jxe_l, jx_s_l, jx_e_l, &
                      iys_l, iye_l, iy_s_l, iy_e_l, iz1, iz2, myid_l, ncpu, np)

      ! 1-D FFT in y over [jxs,jxe] for all z
      do iz = iz1, iz2
        do ix = jxs_l, jxe_l, 2
          do iy = 1, ny
            a2d(1,iy) = at(iy,ix,iz)
            a2d(2,iy) = at(iy,ix+1,iz)
          end do
          call cfftf(ny, a2d(1,1), trigc_l(1))
          do iy = 1, ny
            at(iy,ix,iz)   = a2d(1,iy)
            at(iy,ix+1,iz) = a2d(2,iy)
          end do
        end do
      end do

      if (isgn == -1) then
        call ytox_trans(at, ax, nxp2, ny, jxs_l, jxe_l, jx_s_l, jx_e_l, &
                        iys_l, iye_l, iy_s_l, iy_e_l, iz1, iz2, myid_l, ncpu, np)
      end if

    else   ! isgn > 0 — inverse transform

      if (isgn == 1) then
        call xtoy_trans(ax, at, nxp2, ny, jxs_l, jxe_l, jx_s_l, jx_e_l, &
                        iys_l, iye_l, iy_s_l, iy_e_l, iz1, iz2, myid_l, ncpu, np)
      end if

      ! 1-D FFT in y over [jxs,jxe] for all z
      do iz = iz1, iz2
        do ix = jxs_l, jxe_l, 2
          do iy = 1, ny
            a2d(1,iy) = at(iy,ix,iz)
            a2d(2,iy) = at(iy,ix+1,iz)
          end do
          call cfftb(ny, a2d(1,1), trigc_l(1))
          do iy = 1, ny
            at(iy,ix,iz)   = a2d(1,iy)
            at(iy,ix+1,iz) = a2d(2,iy)
          end do
        end do
      end do
      call ytox_trans(at, ax, nxp2, ny, jxs_l, jxe_l, jx_s_l, jx_e_l, &
                      iys_l, iye_l, iy_s_l, iy_e_l, iz1, iz2, myid_l, ncpu, np)

      ! 1-D FFT in x over [iys,iye] for all z
      do iz = iz1, iz2
        do iy = iys_l, iye_l
          a_wrk(1) = ax(1,iy,iz)
          do ix = 2, nx
            a_wrk(ix) = ax(ix+1,iy,iz)
          end do
          call rfftb(nx, a_wrk(1), trigx_l(1))
          do ix = 1, nx
            ax(ix,iy,iz) = a_wrk(ix)
          end do
        end do
      end do

    end if
  end subroutine fft2d_mpi

  ! --------------------------------------------------------------------
  subroutine yd_mpi(ay, trigy_l, yk, &
                    nx, ny, ixs_l, ixe_l, ix_s_l, ix_e_l, &
                    iys_l, iye_l, iy_s_l, iy_e_l, &
                    iz1, iz2, myid_l, ncpu, np)
  ! y-derivatives via real FFT + MPI slab transposes.
  ! Wavenumbers normalized by ny.
  ! --------------------------------------------------------------------
    real,    intent(inout) :: ay(nx, iys_l:iye_l, iz1:iz2)
    real,    intent(in)    :: trigy_l(2*ny+15), yk(ny)
    integer, intent(in)    :: nx, ny
    integer, intent(in)    :: ixs_l, ixe_l
    integer, intent(in)    :: ix_s_l(0:np-1), ix_e_l(0:np-1)
    integer, intent(in)    :: iys_l, iye_l
    integer, intent(in)    :: iy_s_l(0:np-1), iy_e_l(0:np-1)
    integer, intent(in)    :: iz1, iz2, myid_l, ncpu, np
    real    :: ayt(ny, ixs_l:ixe_l, iz1:iz2)
    integer :: ix, iy, iz, ii
    real    :: temp

    call xtoy_trans(ay, ayt, nx, ny, ixs_l, ixe_l, ix_s_l, ix_e_l, &
                    iys_l, iye_l, iy_s_l, iy_e_l, iz1, iz2, myid_l, ncpu, np)

    do iz = iz1, iz2
      do ix = ixs_l, ixe_l
        call rfftf(ny, ayt(1,ix,iz), trigy_l)
        ii             = 1
        ayt(1,ix,iz)   = 0.0
        ayt(ny,ix,iz)  = 0.0
        do iy = 2, ny-1, 2
          ii               = ii + 1
          temp             = ayt(iy,ix,iz)
          ayt(iy,ix,iz)    = -yk(ii)*ayt(iy+1,ix,iz)
          ayt(iy+1,ix,iz)  =  yk(ii)*temp
        end do
        call rfftb(ny, ayt(1,ix,iz), trigy_l)
      end do
    end do

    call ytox_trans(ayt, ay, nx, ny, ixs_l, ixe_l, ix_s_l, ix_e_l, &
                    iys_l, iye_l, iy_s_l, iy_e_l, iz1, iz2, myid_l, ncpu, np)
  end subroutine yd_mpi

  ! --------------------------------------------------------------------
  subroutine dealias
  ! 2/3-rule spectral dealiasing of all prognostic fields (u,v,w,e,t).
  ! Applies forward 2-D FFT, zeroes modes beyond the 2/3 cutoff,
  ! then inverse 2-D FFT.
  ! --------------------------------------------------------------------
    real :: wve (nny,    jxs:jxe, izs:ize)
    real :: wves(nnxp2, iys:iye, izs:ize)
    integer :: ix, iy, iz, iscl, ix_cut, iy_cut_l, iy_cut_u

    ix_cut   = 2*int(real(nnx)/3.) + 3
    iy_cut_l = int(real(nny)/3.) + 2
    iy_cut_u = nnyp2 - iy_cut_l

    ! u
    call fft2d_mpi(u(1,iys,izs), wve(1,jxs,izs), trigx(1,1), trigc, &
                   nnx, nny, jxs, jxe, jx_s, jx_e, iys, iye, iy_s, iy_e, &
                   izs, ize, myid, ncpu_s, numprocs, -2)
    call sharp(wve, ix_cut, iy_cut_l, iy_cut_u)
    call fft2d_mpi(u(1,iys,izs), wve(1,jxs,izs), trigx(1,1), trigc, &
                   nnx, nny, jxs, jxe, jx_s, jx_e, iys, iye, iy_s, iy_e, &
                   izs, ize, myid, ncpu_s, numprocs,  2)

    ! v
    call fft2d_mpi(v(1,iys,izs), wve(1,jxs,izs), trigx(1,1), trigc, &
                   nnx, nny, jxs, jxe, jx_s, jx_e, iys, iye, iy_s, iy_e, &
                   izs, ize, myid, ncpu_s, numprocs, -2)
    call sharp(wve, ix_cut, iy_cut_l, iy_cut_u)
    call fft2d_mpi(v(1,iys,izs), wve(1,jxs,izs), trigx(1,1), trigc, &
                   nnx, nny, jxs, jxe, jx_s, jx_e, iys, iye, iy_s, iy_e, &
                   izs, ize, myid, ncpu_s, numprocs,  2)

    ! w
    call fft2d_mpi(w(1,iys,izs), wve(1,jxs,izs), trigx(1,1), trigc, &
                   nnx, nny, jxs, jxe, jx_s, jx_e, iys, iye, iy_s, iy_e, &
                   izs, ize, myid, ncpu_s, numprocs, -2)
    call sharp(wve, ix_cut, iy_cut_l, iy_cut_u)
    call fft2d_mpi(w(1,iys,izs), wve(1,jxs,izs), trigx(1,1), trigc, &
                   nnx, nny, jxs, jxe, jx_s, jx_e, iys, iye, iy_s, iy_e, &
                   izs, ize, myid, ncpu_s, numprocs,  2)

    ! e
    call fft2d_mpi(e(1,iys,izs), wve(1,jxs,izs), trigx(1,1), trigc, &
                   nnx, nny, jxs, jxe, jx_s, jx_e, iys, iye, iy_s, iy_e, &
                   izs, ize, myid, ncpu_s, numprocs, -2)
    call sharp(wve, ix_cut, iy_cut_l, iy_cut_u)
    call fft2d_mpi(e(1,iys,izs), wve(1,jxs,izs), trigx(1,1), trigc, &
                   nnx, nny, jxs, jxe, jx_s, jx_e, iys, iye, iy_s, iy_e, &
                   izs, ize, myid, ncpu_s, numprocs,  2)

    ! scalars (t is stored with scalar index as third dimension)
    do iscl = 1, nscl
      do iz = izs, ize
        do iy = iys, iye
          do ix = 1, nnx
            wves(ix,iy,iz) = t(ix,iy,iscl,iz)
          end do
        end do
      end do
      call fft2d_mpi(wves(1,iys,izs), wve(1,jxs,izs), trigx(1,1), trigc, &
                     nnx, nny, jxs, jxe, jx_s, jx_e, iys, iye, iy_s, iy_e, &
                     izs, ize, myid, ncpu_s, numprocs, -2)
      call sharp(wve, ix_cut, iy_cut_l, iy_cut_u)
      call fft2d_mpi(wves(1,iys,izs), wve(1,jxs,izs), trigx(1,1), trigc, &
                     nnx, nny, jxs, jxe, jx_s, jx_e, iys, iye, iy_s, iy_e, &
                     izs, ize, myid, ncpu_s, numprocs,  2)
      do iz = izs, ize
        do iy = iys, iye
          do ix = 1, nnx
            t(ix,iy,iscl,iz) = wves(ix,iy,iz)
          end do
        end do
      end do
    end do
  end subroutine dealias

  ! --------------------------------------------------------------------
  subroutine sharp(wve, ix_cut, iy_cut_l, iy_cut_u)
  ! Sharp spectral cutoff filter on an array stored in 2-D FFT coefficient
  ! order (as produced by fft2d_mpi with isgn = -2).
  ! Zeroes y-modes in [iy_cut_l, iy_cut_u] on all pencils, then zeroes
  ! all modes on x-pencils with index >= ix_cut.
  ! --------------------------------------------------------------------
    real,    intent(inout) :: wve(nny, jxs:jxe, izs:ize)
    integer, intent(in)    :: ix_cut, iy_cut_l, iy_cut_u
    integer :: ix, iy, iz

    do iz = izs, ize
      do ix = jxs, jxe
        do iy = iy_cut_l, iy_cut_u
          wve(iy,ix,iz) = 0.0
        end do
      end do
    end do

    if (jxe >= ix_cut) then
      do iz = izs, ize
        do ix = max(jxs, ix_cut), jxe
          do iy = 1, nny
            wve(iy,ix,iz) = 0.0
          end do
        end do
      end do
    end if
  end subroutine sharp

end module mod_fft
