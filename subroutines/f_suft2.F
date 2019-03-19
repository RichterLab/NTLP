      subroutine f_suft2(rbuf,nnx,mxs,mxe,iys,iye,nscl,
     +                  tau13m,tau23m,taut3m,t_grnd)
c
c ------ fill surface arrays on root processors
c
      real rbuf(2+2*nscl,mxs:mxe,iys:iye)
      real tau13m(nnx,iys:iye), tau23m(nnx,iys:iye),
     +     taut3m(nnx,iys:iye,nscl), t_grnd(nnx,iys:iye,nscl)
c
      do iy=iys,iye
      do ix=mxs,mxe
         tau13m(ix,iy) = rbuf(1,ix,iy)
         tau23m(ix,iy) = rbuf(2,ix,iy)
      enddo
      enddo
      do iscl=1,nscl
         do iy=iys,iye
         do ix=mxs,mxe
            taut3m(ix,iy,iscl) = rbuf(2+iscl,ix,iy)
            t_grnd(ix,iy,iscl) = rbuf(2+nscl+iscl,ix,iy)
         enddo
         enddo
      enddo
c
      return
      end
