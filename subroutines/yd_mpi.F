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
