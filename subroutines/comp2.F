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
