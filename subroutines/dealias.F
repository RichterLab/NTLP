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
