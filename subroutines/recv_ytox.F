      subroutine recv_ytox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
      real f(nx,iys:iye,izs:ize), ft(iys:iye,ixs:ixe,izs:ize)
c
      do k=izs,ize
      do i=ixs,ixe
      do j=iys,iye
         f(i,j,k) = ft(j,i,k)
      enddo
      enddo
      enddo
c
      return
      end
