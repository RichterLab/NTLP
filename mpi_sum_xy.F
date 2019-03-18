      subroutine mpi_sum_xy(f,myid,iss,ise,nsend)
c
c --------- get horizontal x-y sum over a set of proccessors [iss:ise]
c           for vector f(i). f(i) is overwritten. skip if single processor
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real work(nsend,iss:ise), f(nsend)
c
      if(iss .eq. ise) go to 999
c
      do j=1,nsend
         work(j,myid) = f(j)
         f(j)         = 0.0
      enddo
      do i=iss,ise
         if(i .ne. myid) then
            call mpi_sendrecv(work(1,myid),nsend,mpi_real8,i,1,
     +               work(1,i),nsend,mpi_real8,i,1,
     +           mpi_comm_world,istatus,ierr)
         endif
      enddo
      do i=iss,ise
      do j=1,nsend
         f(j) = f(j) + work(j,i)
      enddo
      enddo
c
  999 continue
c
      return
      end
