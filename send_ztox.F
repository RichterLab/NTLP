      subroutine send_ztox(g,gt,nz,ixs,ixe,iys,iye,izs,ize)
c
c ------------- grab correct chunk of array to be sent,
c               account for ghost points
c
      real g(0:nz+1,iys:iye,ixs:ixe), gt(izs-1:ize+1,iys:iye,ixs:ixe)
c
      do j=iys,iye
      do i=ixs,ixe
      do k=izs-1,ize+1
         gt(k,j,i) = g(k,j,i)
      enddo
      enddo
      enddo
c
      return
      end
