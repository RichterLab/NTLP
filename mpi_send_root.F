      subroutine mpi_send_root(fs,num,myid,np,ns)
c
c ---------- send root results to other processors above it
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real fs(num)
c
      if(np .eq. 1) go to 999
c
      irow_r = mod(myid,ns)
      if(myid .ge. ns) then
        call mpi_recv(fs(1),num,mpi_real8,irow_r,1,
     +       mpi_comm_world,istatus,ierr)
      else
        do l=irow_r+ns,np-1,ns
           call mpi_send(fs(1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
  999 continue
c
      return
      end
