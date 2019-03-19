      subroutine recv_ztox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
      real f(nx,iys:iye,izs-1:ize+1), ft(izs-1:ize+1,iys:iye,ixs:ixe)
c
      do i=ixs,ixe
      do j=iys,iye
      do k=izs-1,ize+1
         f(i,j,k) = ft(k,j,i)
      enddo
      enddo
      enddo
c
      return
      end
