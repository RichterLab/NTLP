      subroutine recv_xtoz(g,gt,nz,ixs,ixe,iys,iye,izs,ize)
      real g(0:nz+1,iys:iye,ixs:ixe), gt(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         g(k,j,i) = gt(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
