      subroutine comp_p
c
c --------- setup pressure solver
c
      use pars
      use fftwk
      use fields
      use con_data
      use con_stats
      include 'mpif.h'
      real fnt1(nnx,iys:iye,izs:ize)
      real fs(nnx,iys:iye,2), fr(nnx,iys:iye,2)
      integer istatus(mpi_status_size)
c
      gami = 1.0/dtgama
c
      nb = myid - ncpu_s
      nt = myid + ncpu_s
c
c ------------ Send both r3 and updated w (from comp1)
c              to processor above the current myid.
c
      if(iss .eq. 0) then
         nb = mpi_proc_null
      endif
      if(ise .eq. numprocs-1) then
         nt = mpi_proc_null
      endif
      nsend = 2*nnx*(iye + 1 - iys)
      nrecv = nsend
      do iy=iys,iye
      do ix=1,nnx
         fs(ix,iy,1) = r3(ix,iy,ize)
         fs(ix,iy,2) = w(ix,iy,ize)
      enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nt,2,
     +     fr(1,iys,1),nrecv,mpi_real8,nb,2,
     +     mpi_comm_world,istatus,ierr)
      if(iss .ne. 0) then
         do iy=iys,iye
         do ix=1,nnx
            r3(ix,iy,izs-1) = fr(ix,iy,1)
            w(ix,iy,izs-1)  = fr(ix,iy,2)
         enddo
         enddo
      endif
c
c ----------- setup general pressure calculation
c             relies on rhs from step n-1 being included 
c             in velocity-arrays already
c
      do iz=izs,ize
         izm1 = iz -1
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = u(ix,iy,iz)*gami + r1(ix,iy,iz)
         enddo
         enddo
         call xderivp(fnt1(1,iys,iz),trigx(1,1),xk(1),nnx,iys,iye)
c
         if(iz .eq. 1) then
            do iy=iys,iye
            do ix=1,nnx
                p(ix,iy,iz) = fnt1(ix,iy,iz) +  
     +                     ((w(ix,iy,iz) -wbc(ix,iy,2))*gami +
     +                       r3(ix,iy,iz))*dzw_i(iz)
            enddo
            enddo
         else if(iz .eq. nnz) then
            do iy=iys,iye
            do ix=1,nnx
                p(ix,iy,iz) = fnt1(ix,iy,iz) + 
     +                     ((wbc(ix,iy,1) - w(ix,iy,izm1))*gami -
     +                      r3(ix,iy,izm1))*dzw_i(iz)
            enddo
            enddo
         else 
            do iy=iys,iye
            do ix=1,nnx
                p(ix,iy,iz) = fnt1(ix,iy,iz) + 
     +                    ((w(ix,iy,iz)  - w(ix,iy,izm1))*gami +
     +                      r3(ix,iy,iz) - r3(ix,iy,izm1))*dzw_i(iz)
            enddo
            enddo
         endif
c
c --------- end z loop
c
      enddo
c
c ----------- check for radiation boundary condition, all processors
c
      if(ibcu .eq. 1) then
        do iy=iys,iye
        do ix=1,nnx
           ptop(ix,iy,1) = pbc(ix,iy,1)
           ptop(ix,iy,2) = pbc2(ix,iy,1)
        enddo
        enddo
      endif
c
c --------- now y contribution
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = v(ix,iy,iz)*gami + r2(ix,iy,iz)
         enddo
         enddo
      enddo
c
      call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
             p(ix,iy,iz) = p(ix,iy,iz) + fnt1(ix,iy,iz) 
         enddo
         enddo
      enddo
c
      call pressure
c
      return
      end
