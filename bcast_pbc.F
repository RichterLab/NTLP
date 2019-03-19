      subroutine bcast_pbc
c
c ---- send upper boundary conditions to other processors 
c      for fft solution of pressure
c
      use pars
      use fields
      include 'mpif.h'
      integer istatus(mpi_status_size),ierr
c
      if(numprocs .eq. 1) go to 999
c
      irow_r = mod(myid,ncpu_s)
      irow_t = is_s(numprocs-1) + irow_r
      num = nnx*(iye+1-iys)
c
c
      if(iss .ne. is_s(numprocs-1)) then
c
c ------ not in the top row, receive from top
c
        call mpi_recv(pbc(1,iys,1),num,mpi_real8,irow_t,1,
     +       mpi_comm_world,istatus,ierr)
      else
c
c ------ myid is in the top row, send to everyone below
c
        do l=irow_r,irow_t-ncpu_s,ncpu_s
           call mpi_send(pbc(1,iys,1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
c --------- same thing for another variable
c
      if(iss .ne. is_s(numprocs-1)) then
c
c ------ not in the top row, receive from top
c
        call mpi_recv(pbc2(1,iys,1),num,mpi_real8,irow_t,1,
     +       mpi_comm_world,istatus,ierr)
      else
c
c ------ in the top row, send to everyone below
c
        do l=irow_r,irow_t-ncpu_s,ncpu_s
           call mpi_send(pbc2(1,iys,1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
  999 continue
c
      return
      end
