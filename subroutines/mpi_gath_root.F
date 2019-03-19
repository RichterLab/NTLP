      subroutine mpi_gath_root(fs,fr,iz_s,iz_e,izs,ize,nz,myid,np,ns)
c
c ---------- gather results on root processors
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      integer iz_s(0:np-1), iz_e(0:np-1)
      real fs(izs:ize), fr(1:nz)
c
      if(np .eq. 1) go to 999
c
      irow_r = mod(myid,ns)
      if(myid .gt. ns) then
        call mpi_send(fs(izs),ize+1-izs,mpi_real8,irow_r,1,
     +       mpi_comm_world,ierr)
      else
        do l=irow_r+ns,np-1,ns
           ind = iz_s(l) + 1
           num = iz_e(l) + 1 - iz_s(l)
           call mpi_recv(fr(ind),num,mpi_real8,l,1,
     +       mpi_comm_world,istatus,ierr)
        enddo
      endif
c
  999 continue
c
      return
      end
