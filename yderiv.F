      subroutine yderiv(ay,trigy,yk,nnx,nny)
c
c -------- get multiple y derivatives using fftpack routines
c          use fftpack storage a_0, (a1,b1), (a2,b2), ...,
c          assumes that first wavenumber yk(1) = 0.0
c
c          assumes that wavenumbers are normalized by number of points
c
      real yk(nny), trigy(2*nny+15), ay(nnx,nny)
      real a_trans(nny)
c
c     fn = 1.0/float(nny)
      do ix=1,nnx
         do iy=1,nny
            a_trans(iy) = ay(ix,iy)
         enddo
         call rfftf(nny,a_trans(1),trigy)
         ii = 1
         a_trans(1)   = 0.0
         a_trans(nny) = 0.0
         do iy=2,nny-1,2
            ii            = ii + 1
            temp          = a_trans(iy)
            a_trans(iy)   = -yk(ii)*a_trans(iy+1)
            a_trans(iy+1) = yk(ii)*temp
         enddo
         call rfftb(nny,a_trans(1),trigy)
         do iy=1,nny
            ay(ix,iy) = a_trans(iy)
         enddo
      enddo
c
      return
      end
