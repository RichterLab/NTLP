      subroutine send_xtoy(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
c
c ------------- grab correct chunk of array to be sent
c
      real f(nx,iys:iye,izs:ize), ft(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         ft(i,j,k) = f(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
