      subroutine save_viz(it)
c
c --------------- save multiple (x-y), (x-z), (y-z), planes of data .
c                 modify recl in all open statements for more or less
c                 variables. 
c                 Constant - x, implies yz planes
c                 Constant - y, implies xz planes
c                 Constant - z, implies xy planes
c
c ------------- routine uses send/recv to get information in y-z planes
c
      use pars
      use fields
      use con_data
      use con_stats
      use fftwk
      use particles
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
c
      parameter(nvar_o = 6)
c
      integer ix_pick(maxnx),  iy_pick(maxny),  iz_pick(maxnz), 
     +        ix_order(maxnx), iy_order(maxny), iz_order(maxnz)
c
      integer istatus(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
c
      real(kind=4), dimension(nvar_o,nny,izs:ize) :: temp_x
      real(kind=4), dimension(nvar_o,nnx,izs:ize) :: temp_y
      real(kind=4), dimension(nvar_o,nnx,iys:iye) :: temp_z
      real, dimension(nvar_o,iys:iye,izs:ize)     :: buf_send
c
c ------------- don't touch
c
      data iviz_x,  iviz_y,  iviz_z  /0, 0, 0/
      data ionce_x, ionce_y, ionce_z, istuff /0, 0, 0, 0/
      data ix_pick, iy_pick, iz_pick /maxnx*0, maxny*0, maxnz*0/
      data ix_order, iy_order, iz_order /maxnx*0, maxny*0, maxnz*0/
      save iviz_x,  iviz_y,  iviz_z, 
     +     ix_pick, iy_pick, iz_pick, 
     +     ix_order, iy_order, iz_order,
     +     ionce_x, ionce_y, ionce_z, istuff,
     +     npln_x, npln_y, npln_z
c
c
c ----------- turn on z levels to save. Customize for your own use.
c             Set iz_pick(iz) = iz, ix_pick(ix) = ix, iy_pick(iy) = iy
c             Data is round-robin alternated in the data file for more than
c             1 plane for any particular view.
c
      iz_pick(12) = 12
      iz_pick(28) = 28
      iz_pick(64) = 64
      !iz_pick(20) = 20
      !iz_pick(45) = 45
      !iz_pick(60) = 60
c
c -------------- pick an x-z plane of data (can add more)
c
      iy_pick(nny/2) = nny/2
c     iy_pick(nny)   = nny
c
c -------------- pick a y-z plane of data (can add more)
c
      ix_pick(nnx/2) = nnx/2
c     ix_pick(nnx)   = nnx
c
c ------ find total number of z's turned on and open file once
c
      if (ionce_z .eq. 0) then
         npln_z = 0
         do k=1,nnz
            if(iz_pick(k) .eq. k) then
               npln_z = npln_z + 1
               iz_order(k) = npln_z
            endif
         enddo
         ionce_z = 1
         iviz_z =  -npln_z
         if (npln_z .ne. 0) then
            call mpi_file_open(mpi_comm_world, path_viz_xy,
     +                         mpi_mode_create+mpi_mode_rdwr,
     +                         mpi_info_null, nviz_z, ierr)
            disp = 0
            call mpi_file_set_view(nviz_z,disp,mpi_real4,mpi_real4,
     +                            'native',mpi_info_null,ierr)
         endif
      endif
c
c ------ find total number of y's turned on and open file once
c
      if (ionce_y .eq. 0) then
         npln_y = 0
         do j = 1,nny
            if(iy_pick(j) .eq. j) then
               npln_y = npln_y + 1
               iy_order(j) = npln_y
            endif
         enddo
         ionce_y = 1
         iviz_y  = -npln_y
         if (npln_y .ne. 0) then
            call mpi_file_open(mpi_comm_world, path_viz_xz,
     +                         mpi_mode_create+mpi_mode_rdwr,
     +                         mpi_info_null, nviz_y, ierr)
            disp = 0
            call mpi_file_set_view(nviz_y,disp,mpi_real4,mpi_real4,
     +                            'native',mpi_info_null,ierr)
         endif
      endif
c
c ------ find total number of x's turned on and open file once
c
      if (ionce_x .eq. 0) then
         npln_x = 0
         do i=1,nnx
            if(ix_pick(i) .eq. i) then
               npln_x = npln_x + 1
               ix_order(i) = npln_x
            endif
         enddo
         ionce_x = 1
         iviz_x  = -npln_x
         if (npln_x .ne. 0) then
            call mpi_file_open(mpi_comm_world, path_viz_yz,
     +                         mpi_mode_create+mpi_mode_rdwr,
     +                         mpi_info_null, nviz_x, ierr)
            disp = 0
            call mpi_file_set_view(nviz_x,disp,mpi_real4,mpi_real4,
     +                            'native',mpi_info_null,ierr)
         endif
      endif
c
      if(istuff .eq. 0 .and. l_root) then
         open(nviz_s,file=path_stuf)
         istuff = 1
      endif
c
c --------- write data, subtract t_ref to increase
c           resolution on 32 bit machines
c
c ---------- xy planes of data
c
      iviz_z  = iviz_z + npln_z
      nsize   = int(nvar_o,k8)*nnx*nny
      nsize2  = int(nvar_o,k8)*nnx*(iys-1)
      n_write = nvar_o*nnx*(iye+1-iys)
      do k=izs,ize
         if(iz_pick(k) .eq. k) then
            km1 = k - 1
            do j=iys,iye
            do i=1,nnx
               temp_z(1,i,j) = u(i,j,k)
               temp_z(2,i,j) = v(i,j,k)
               temp_z(3,i,j) = w(i,j,k)
               temp_z(4,i,j) = (t(i,j,1,k) - t_ref)
c
c ---------- get just the fluctuating pressure field
c
               temp_z(5,i,j) = p(i,j,k) - pxym(k)
     +                        -(e(i,j,k) + e(i,j,km1))/3.0
     +                        -0.5*((u(i,j,k) + stokes(k))**2 +
     +                               v(i,j,k)*v(i,j,k) +
     +                         0.5*(w(i,j,k)*w(i,j,k) + 
     +                              w(i,j,km1)*w(i,j,km1)))

               temp_z(6,i,j) = u(i,j,k)-uxym(k)
               !temp_z(6,i,j) = partsrc(i,j,k,1)
               !temp_z(7,i,j) = partsrc(i,j,k,2)
               !temp_z(8,i,j) = partsrc(i,j,k,3)
               !temp_z(9,i,j) = partcount(i,j,k)/dx/dy/dzw(k)
            enddo
            enddo
#if defined(SWAP)
            call byteswap(temp_z)
#endif
            offset = int((iviz_z + iz_order(k) - 1),k8)*nsize + nsize2
            call mpi_file_write_at(nviz_z,offset,temp_z,n_write,
     +                             mpi_real4,istatus,ierr)
            if (ierr .ne. 0) go to 9991
         endif
      enddo
c
c ---------- xz planes of data
c
      iviz_y = iviz_y + npln_y
      nsize  = int(nvar_o,k8)*nnx*nnz
      nsize2 = int(nvar_o,k8)*nnx*(izs-1)
      nwrite = nvar_o*nnx*(ize+1-izs)
      do j=iys,iye
         if(iy_pick(j) .eq. j) then
            do k=izs,ize
            km1 = k - 1
            do i=1,nnx
               temp_y(1,i,k) = u(i,j,k)
               temp_y(2,i,k) = v(i,j,k)
               temp_y(3,i,k) = w(i,j,k)
               temp_y(4,i,k) = (t(i,j,1,k) - t_ref)
c
c ---------- get just the fluctuating pressure field
c
               temp_y(5,i,k) =  p(i,j,k) - pxym(k)
     +                          -(e(i,j,k)+e(i,j,km1))/3.0
     +                          -0.5*((u(i,j,k)+stokes(k))**2 +
     +                               v(i,j,k)*v(i,j,k) +
     +                           0.5*(w(i,j,k)*w(i,j,k) + 
     +                                w(i,j,km1)*w(i,j,km1)))

               temp_y(6,i,k) = u(i,j,k) - uxym(k)
               !temp_y(6,i,k) = partsrc(i,j,k,1)
               !temp_y(7,i,k) = partsrc(i,j,k,2)
               !temp_y(8,i,k) = partsrc(i,j,k,3)
               !temp_y(9,i,k) = partcount(i,j,k)/dx/dy/dzw(k)
            enddo
            enddo
#if defined(SWAP)
            call byteswap(temp_y)
#endif
            offset = int((iviz_y + iy_order(j) - 1),k8)*nsize + nsize2
            call mpi_file_write_at(nviz_y,offset,temp_y,nwrite,
     +                                mpi_real4,istatus,ierr)
            if (ierr .ne. 0) goto 9992
         endif
      enddo
c
c ---------- yz planes that cut across all processors
c            just have root node on that slab write data
c
      iviz_x  = iviz_x + npln_x
      n_write = nvar_o*nny*(ize+1-izs)
      nsize   = int(nvar_o,k8)*nny*nnz
      nsize2  = int(nvar_o,k8)*nny*(izs-1)
      n_send  = nvar_o*(ize+1-izs)*(iye+1-iys)
      do i=1,nnx
         if(ix_pick(i) .eq. i) then
c
c ----------- build send buffer
c
            do k=izs,ize
            km1 = k - 1
            do j=iys,iye
               buf_send(1,j,k) = u(i,j,k)
               buf_send(2,j,k) = v(i,j,k)
               buf_send(3,j,k) = w(i,j,k)
               buf_send(4,j,k) = (t(i,j,1,k) - t_ref)
c
c ---------- get just the fluctuating pressure field
c
               buf_send(5,j,k) = p(i,j,k) - pxym(k)
     +                          -(e(i,j,k) + e(i,j,km1))/3.0
     +                          -0.5*((u(i,j,k) + stokes(k))**2 +
     +                               v(i,j,k)*v(i,j,k) +
     +                           0.5*(w(i,j,k)*w(i,j,k) + 
     +                                w(i,j,km1)*w(i,j,km1)))
              buf_send(6,j,k) = u(i,j,k)-uxym(k)
              !buf_send(6,j,k) = partsrc(i,j,k,1)
              !buf_send(7,j,k) = partsrc(i,j,k,2)
              !buf_send(8,j,k) = partsrc(i,j,k,3)
              !buf_send(9,j,k) = partcount(i,j,k)/dx/dy/dzw(k)
            enddo
            enddo
            if(myid .ne. iss) then
              call mpi_send(buf_send(1,iys,izs),n_send,
     +                      mpi_real8,iss,1,
     +                      mpi_comm_world,ierr)
            else
              do k=izs,ize
              do j=iys,iye
              do ii=1,nvar_o
                 temp_x(ii,j,k) = buf_send(ii,j,k)
              enddo
              enddo
              enddo
              do l=iss+1,ise
                 call recv_yz_var(temp_x,nvar_o,nny,
     +                            iy_s(l),iy_e(l),izs,ize,l)
              enddo
#if defined(SWAP)
              call byteswap(temp_x)
#endif
              offset = int((iviz_x + ix_order(i) - 1),k8)*nsize + nsize2
              call mpi_file_write_at(nviz_x,offset,temp_x,n_write,
     +                          mpi_real4,istatus,ierr)
              if (ierr .ne. 0) goto 9993
            endif
         endif
      enddo
c
c ------------- ascii file with facts in it that goes
c               with visualization
c
      if(l_root) then
         write(nviz_s,5000) time, amonin, zi, utau
 5000    format(4e20.8)
      endif
c
c ---- last time step close the files
c
!      if (it .eq. itmax) then
!         call mpi_file_close(nviz_z, ierr)
!         call mpi_file_close(nviz_y, ierr)
!         call mpi_file_close(nviz_x, ierr)
!         if (l_root) then
!            close(nviz_s)
!         endif
!      endif
       if (it .eq. itmax .or. mtape) then
        if(npln_z .ne. 0) then
            call mpi_file_close(nviz_z, ierr)
            ionce_z = 0
         endif
         if(npln_y .ne. 0) then
            call mpi_file_close(nviz_y, ierr)
            ionce_y = 0
         endif
         if(npln_x .ne. 0) then
            call mpi_file_close(nviz_x, ierr)
            ionce_x = 0
         endif
         if (l_root) then
            close(nviz_s)
            istuff = 0
         endif
      endif
c
      return
c --------------------------  errors in writing viz file
 9991 continue
      write(6,6000) nviz_z, iz
 6000 format(' SR. SAVE_VIS:',/,
     +       '    trouble cannot write xy viz file on unit = ',i2,/,
     +       '             at iz = ',i4)
      call mpi_finalize(ierr)
      stop
c --------------------------  errors in writing viz file
 9992 continue
      write(6,6100) nviz_y, iz, iviz_y
 6100 format(' SR. SAVE_VIS:',/,
     +       '    trouble cannot write xz viz file on unit = ',i2,/,
     +       '             at iz = ',i4,/,
     +       '            iviz_y = ',i8)
c --------------------------  errors in writing viz file
 9993 continue
      write(6,6200) nviz_x, iz, iviz_x
 6200 format(' SR. SAVE_VIS:',/,
     +       '    trouble cannot write yz viz file on unit = ',i2,/,
     +       '             at iz = ',i4,/,
     +       '            iviz_x = ',i8)
      call mpi_finalize(ierr)
      stop
      end
