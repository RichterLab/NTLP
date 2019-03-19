      subroutine save_v(it)
c
c --------------- save 3d fields
c
      use pars
      use fields
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
      logical there
c
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)                 nsize, nsize2
c
      real, allocatable, dimension(:,:,:) :: temp
      allocate(temp(nvar,nnx,iys:iye))
c
c ---- open file
c
      call mpi_file_open(mpi_comm_world, path_sav_v,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, nvel, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ---- write data
c
      nsize   = int(nvar,k8)*nnx*nny
      nsize2  = int(nvar,k8)*nnx*(iys-1)
      n_write = nvar*nnx*(iye+1-iys)
c
      do k=izs,ize
         do j = iys,iye
         do i = 1,nnx
            temp(1,i,j)    = u(i,j,k)
            temp(2,i,j)    = v(i,j,k)
            temp(3,i,j)    = w(i,j,k)
            temp(nvar,i,j) = e(i,j,k)
         enddo
         enddo
         do is = 1,nscl
            do j = iys,iye
            do i = 1,nnx
               temp(3+is,i,j) = t(i,j,is,k)
            enddo
            enddo
         enddo


#if defined(SWAP)
      call byteswap(temp)
#endif
c

         offset = int((k-1),k8)*nsize + nsize2
c         call mpi_file_write_at(nvel,offset,temp,n_write,
c     +                              mpi_real8,status,ierr)
         call mpi_file_write_at_all(nvel,offset,temp,n_write,
     +                              mpi_real8,status,ierr)
         if (ierr /= 0) goto 9991
c
      enddo

c
c ---- close file
c
      call mpi_file_close(nvel, ierr)

c
c ---- check file
c
      if (l_root) then
         inquire(file=path_sav_v,exist=there)
         if(.not.there) then
            write(6,8000) nvel,myid
            call mpi_finalize(ierr)
            stop
         endif
         write(6,7000) it,path_sav_v
      endif
c
      deallocate(temp)
c
      return
c --------------------------  errors in writing restart file
 9991 continue
      write(6,6000) nvel, iz
 6000 format(' SR. SAVE_V:',/,
     +       '    trouble cannot write restart file on unit = ',i2,/,
     +       '             at iz = ',i4)
      call mpi_finalize(ierr)
      stop
c --------------------
 7000 format(' **** DATA SET AT IT = ',I6,/,
     +       '      VELOCITY DATA IS WRITTEN IN FILE  ',a80)
 8000 format(' in SAVE_V: trouble writing file ',i5,'  myid = ',i5,
     +       ' at iz = ',i5)
      end
