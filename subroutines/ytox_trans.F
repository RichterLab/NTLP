      subroutine ytox_trans(g,f,nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,
     +           myid,ncpu_s,np)
c 
c ------- transpose array g(ny,ixs:ixe,iz1:iz2) ---> f(nx,iys:iye,iz1:iz2)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,iz1:iz2), 
     +     g(ny,ixs:ixe,iz1:iz2)
      real ft(nx*(iye+1-iys)*(iz2 - iz1 + 1)),
     +     gt(ny*(ixe+1-ixs)*(iz2 - iz1 + 1))
c
      integer ix_s(0:np-1), ix_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      jk = (iye - iys + 1)*(iz2 - iz1 + 1)
      ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)
c
c ----------- get cpus on slab for myid
c
      islab = myid/ncpu_s
      iss   = islab*ncpu_s
      ise   = iss + ncpu_s - 1
      do i=iss,ise
         nsend = (iy_e(i) - iy_s(i) + 1)*ik
         nrecv = (ix_e(i) - ix_s(i) + 1)*jk
         if(i .eq. myid) then
            call send_ytox(g,ft(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(i),iy_e(i),iz1,iz2)
         else
            call send_ytox(g,gt(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(i),iy_e(i),iz1,iz2)
            call mpi_sendrecv(
     +          gt(1),nsend,mpi_real8,i,1,
     +          ft(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_ytox(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)
      enddo
c
      return
      end
