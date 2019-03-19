      subroutine exchange
c
c ------------- exchange ghost points with mpi,
c               nb and nt are the destination and
c               source nodes. Allows for 1z per cpu
c
      use pars
      use fields
c     use fftwk
      include 'mpif.h'
c
      real fs(nnx,iys:iye,(4+nscl)),fr(nnx,iys:iye,(4+nscl))
      integer istatus(mpi_status_size)
c
      nb = myid - ncpu_s
      nt = myid + ncpu_s
c
c ------------ account for endpoints
c
      if(iss .eq. 0) then
         nb = mpi_proc_null
      endif
      if(ise .eq. numprocs-1) then
         nt = mpi_proc_null
      endif
      nsend = nnx*(iye + 1 - iys)*(4+nscl)
      nrecv = nsend
c
c --------- send top of myid, receive bottom from myid - ncpu_s
c
      do iy=iys,iye
      do ix=1,nnx
         fs(ix,iy,1) = u(ix,iy,ize)
         fs(ix,iy,2) = v(ix,iy,ize)
         fs(ix,iy,3) = w(ix,iy,ize)
         fs(ix,iy,4) = e(ix,iy,ize)
      enddo
      enddo
      do iscl=1,nscl
         jloc = 4 + iscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,jloc) = t(ix,iy,iscl,ize)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nt,0,
     +     fr(1,iys,1),nrecv,mpi_real8,nb,0,
     +     mpi_comm_world,istatus,ierr)
c
      if(iss .ne. 0) then
         izm1 = izs-1
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,izm1) = fr(ix,iy,1)
            v(ix,iy,izm1) = fr(ix,iy,2)
            w(ix,iy,izm1) = fr(ix,iy,3)
            e(ix,iy,izm1) = fr(ix,iy,4)
         enddo
         enddo
         do iscl=1,nscl
            jloc = 4 + iscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izm1) = fr(ix,iy,jloc)
            enddo
            enddo
         enddo
      endif
c
c -------- send bottom of myid, receive bottom from myid + ncpu_s
c
      do iy=iys,iye
      do ix=1,nnx
         fs(ix,iy,1) = u(ix,iy,izs)
         fs(ix,iy,2) = v(ix,iy,izs)
         fs(ix,iy,3) = w(ix,iy,izs)
         fs(ix,iy,4) = e(ix,iy,izs)
      enddo
      enddo
      do iscl=1,nscl
         jloc = 4 + iscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,jloc) = t(ix,iy,iscl,izs)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nb,1,
     +     fr(1,iys,1),nrecv,mpi_real8,nt,1,
     +     mpi_comm_world,istatus,ierr)
c
      if(ise .ne. numprocs-1) then
         izp1 = ize+1
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,izp1) = fr(ix,iy,1)
            v(ix,iy,izp1) = fr(ix,iy,2)
            w(ix,iy,izp1) = fr(ix,iy,3)
            e(ix,iy,izp1) = fr(ix,iy,4)
         enddo
         enddo
         do iscl=1,nscl
            jloc = 4 + iscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izp1) = fr(ix,iy,jloc)
            enddo
            enddo
         enddo
      endif
c
c --------------- send extra scalar points 
c
      nsend = nnx*(iye + 1 - iys)*nscl
      nrecv = nsend
c
c -------------- send top of myid, receive bottom from myid - ncpu_s
c
      izm1 = ize-1
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,iscl) = t(ix,iy,iscl,izm1)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nt,0,
     +     fr(1,iys,1),nrecv,mpi_real8,nb,0,
     +     mpi_comm_world,istatus,ierr)
c
      if(iss .ne. 0) then
         izm2 = izs-2
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izm2) = fr(ix,iy,iscl)
            enddo
            enddo
         enddo
      endif
c
c -------------- send bottom of myid, receive bottom from myid + ncpu_s
c
      izp1 = izs+1
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,iscl) = t(ix,iy,iscl,izp1)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nb,1,
     +     fr(1,iys,1),nrecv,mpi_real8,nt,1,
     +     mpi_comm_world,istatus,ierr)
c
      if(ise .ne. numprocs-1) then
         izp2 = ize+2
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izp2) = fr(ix,iy,iscl)
            enddo
            enddo
         enddo
      endif
c
      return
      end
