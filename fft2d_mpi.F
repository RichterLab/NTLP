      subroutine fft2d_mpi(ax,at,trigx,trigc,nx,ny,
     +           jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           iz1,iz2,myid,ncpu,np,isgn)
c
c -------- get 2d fft using fftpack routines and parallel mpi
c          use fftpack storage a0, (a1,b1), (a2,b2),...,
c
c         isgn = -1 do forward transform, get coefficients
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
c
c         isgn = -2 do forward transform, get coefficients
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is at(ny,jxs:jxe,iz1:iz2)
c
c         isgn =  1 do inverse transform, move to physical space
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
c
c         isgn =  2 do inverse transform, move to physical space
c                   incoming array is at(ny,jxs:jxe,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
c
      real ax(nx+2,iys:iye,iz1:iz2), at(ny,jxs:jxe,iz1:iz2),
     +     trigx(2*nx+15), trigc(4*ny+15),
     +     a2d(2,ny), a_wrk(nx)
      integer jx_s(0:np-1), jx_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      nxp2 = nx + 2
      if(isgn .lt. 0) then
         fn   = 1.0/(float(nx)*float(ny))
c
c ------ 1d fft in x over [iys,iye] for all z
c
         do iz=iz1,iz2
            do iy=iys,iye
               do ix=1,nx
                  a_wrk(ix) = ax(ix,iy,iz)*fn
               enddo
               call rfftf(nx,a_wrk(1),trigx(1))
               ax(1,iy,iz) = a_wrk(1)
               ax(2,iy,iz) = 0.0
               do ix=2,nx
                  ax(ix+1,iy,iz) = a_wrk(ix)
               enddo
               ax(nx+2,iy,iz) = 0.0
            enddo
         enddo
         call xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +        iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c ------ 1d fft in y over [jxs,jxe] for all z
c
         do iz=iz1,iz2
            do ix=jxs,jxe,2
               do iy=1,ny
                  a2d(1,iy) = at(iy,ix,iz)
                  a2d(2,iy) = at(iy,ix+1,iz)
               enddo
               call cfftf(ny,a2d(1,1),trigc(1))
               do iy=1,ny
                  at(iy,ix,iz)   = a2d(1,iy)
                  at(iy,ix+1,iz) = a2d(2,iy)
               enddo
            enddo
         enddo
c
c ---- decide whether to transpose back or leave as is
c
         if(isgn .eq. -1) then
            call ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
         endif
c
      else
c
c ---- decide whether to first transpose or leave as is
c
         if(isgn .eq. 1) then
            call xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
         endif
c
c ------ 1d fft in y over [jxs,jxe] for all z
c
         do iz=iz1,iz2
            do ix=jxs,jxe,2
               do iy=1,ny
                  a2d(1,iy) = at(iy,ix,iz)
                  a2d(2,iy) = at(iy,ix+1,iz)
               enddo
               call cfftb(ny,a2d(1,1),trigc(1))
               do iy=1,ny
                  at(iy,ix,iz)   = a2d(1,iy)
                  at(iy,ix+1,iz) = a2d(2,iy)
               enddo
            enddo
         enddo
         call ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +        iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c ------  1d fft in x over [iys,iye] for all z
c
         do iz=iz1,iz2
            do iy=iys,iye
               a_wrk(1) = ax(1,iy,iz)
               do ix=2,nx
                  a_wrk(ix) = ax(ix+1,iy,iz)
               enddo
               call rfftb(nx,a_wrk(1),trigx(1))
               do ix=1,nx
                  ax(ix,iy,iz) = a_wrk(ix)
               enddo
            enddo
         enddo
      endif
c
      return
      end
