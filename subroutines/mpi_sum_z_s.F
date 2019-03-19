      subroutine mpi_sum_z_s(f,i_root,myid,nsend,nscl,iall)
c
c --------- get sums on root or all processors
c           for all z for vector f(i,nscl)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real recv_b(nsend,nscl), f(nsend,nscl)
c
      if(iall .ne. 1) then
         call mpi_reduce(f(1,1),recv_b(1,1),nsend*nscl,mpi_real8,
     +        mpi_sum,i_root,mpi_comm_world,ierr)
         if(myid .eq. i_root) then
            do iscl=1,nscl
            do i=1,nsend
               f(i,iscl) = recv_b(i,iscl)
            enddo
            enddo
         endif
      else
         call mpi_allreduce(f(1,1),recv_b(1,1),nsend*nscl,mpi_real8,
     +        mpi_sum, mpi_comm_world,ierr)
         do iscl=1,nscl
         do i=1,nsend
            f(i,iscl) = recv_b(i,iscl)
         enddo
         enddo
      endif
c
      return
      end
