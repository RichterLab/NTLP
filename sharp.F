      subroutine sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
c
c --------- sharp cutoff filter for field wve stored
c           in 2d-fft order
c
      use pars
      real wve(nny,jxs:jxe,izs:ize)
c
      do iz=izs,ize
         do ix=jxs,jxe
         do iy=iy_cut_l,iy_cut_u
            wve(iy,ix,iz) = 0.0
         enddo
         enddo
      enddo
c
      if(jxe .lt. ix_cut) go to 999
c
         do iz=izs,ize
            do ix=max(jxs,ix_cut),jxe
            do iy=1,nny
               wve(iy,ix,iz) = 0.0
            enddo
            enddo
         enddo
c
  999 continue
c
      return
      end
