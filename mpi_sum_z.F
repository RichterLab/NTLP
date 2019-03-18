      subroutine mpi_sum_z(f,i_root,myid,nsend,iall)
c
c --------- get sums on root or all processors
c           for all z for vector f(i)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real recv_b(nsend), f(nsend)
c
c -------- just root gets the result
c
      if(iall .ne. 1) then
         call mpi_reduce(f(1),recv_b(1),nsend,mpi_real8,mpi_sum,i_root,
     +                  mpi_comm_world,ierr)
         if(myid .eq. i_root) then
            do i=1,nsend
               f(i) = recv_b(i)
            enddo
         endif
      else
c
c -------- everyone gets the result
c
         call mpi_allreduce(f(1),recv_b(1),nsend,mpi_real8,mpi_sum,
     +                  mpi_comm_world,ierr)
         do i=1,nsend
            f(i) = recv_b(i)
         enddo
      endif
c
      return
      end
