      subroutine ztox_trans(g,f,nx,nz,ixs,ixe,ix_s,ix_e,
     +           iys,iye,izs,ize,iz_s,iz_e,
     +           myid,ncpu_s,numprocs)
c
c ------- transpose array g(0:nz+1,iys:iye,ixs:ixe) 
c                    ---> f(nx,iys:iye,izs-1:ize+1)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,izs-1:ize+1), g(0:nz+1,iys:iye,ixs:ixe)
      real ft(nx*(iye+1-iys)*(ize-izs+3)),
     +     gt((nz+3)*(iye+1-iys)*(ixe-ixs+1))
c
      integer ix_s(0:numprocs-1), ix_e(0:numprocs-1),
     +        iz_s(0:numprocs-1), iz_e(0:numprocs-1)
c
      jk = (ize - izs + 3)*(iye - iys + 1)
      ij = (ixe - ixs + 1)*(iye - iys + 1)
c
c ------------- get starting location
c
      iss = myid - (myid/ncpu_s)*ncpu_s
c
      do i=iss,numprocs-1,ncpu_s
         nsend = (iz_e(i) - iz_s(i) + 3)*ij
         nrecv = (ix_e(i) - ix_s(i) + 1)*jk
         if(i .eq. myid) then
            call send_ztox(g,ft(1),nz,ix_s(myid),ix_e(myid),
     +                  iys,iye,iz_s(i),iz_e(i))
         else
            call send_ztox(g,gt(1),nz,ix_s(myid),ix_e(myid),
     +                  iys,iye,iz_s(i),iz_e(i))
            call mpi_sendrecv(
     +          gt(1),nsend,mpi_real8,i,1,
     +          ft(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_ztox(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iys,iye,iz_s(myid),iz_e(myid))
      enddo
c
      return
      end
