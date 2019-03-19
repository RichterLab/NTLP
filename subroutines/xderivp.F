      subroutine xderivp(ax,trigx,xk,nnx,iys,iye)
c
c -------- get multiple x derivatives using fftpack routines
c          use fftpack storage a0, (a1,b1), (a2,b2),...,an
c          assumes that first wavenumber xk(1) = 0.0
c
c          assumes that wavenumbers are normalized by number of points
c
      real xk(nnx), trigx(2*nnx+15), ax(nnx,iys:iye)
c
c     fn = 1.0/float(nnx)
      do iy=iys,iye
         call rfftf(nnx,ax(1,iy),trigx)
         ii = 1
         ax(1,iy) = 0.0
         ax(nnx,iy) = 0.0
         do ix=2,nnx-1,2
            ii          = ii + 1
            temp        = ax(ix,iy)
            ax(ix,iy)   = -xk(ii)*ax(ix+1,iy)
            ax(ix+1,iy) = xk(ii)*temp
         enddo
         call rfftb(nnx,ax(1,iy),trigx)
      enddo
c
      return
      end
