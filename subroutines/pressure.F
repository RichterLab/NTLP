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
