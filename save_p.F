      subroutine save_p
c
c -------------- save pressure file
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
      integer(kind=k8)              :: nsize, nsize2
c
      real temp(nnx,iys:iye)
c
c ---- open file
c
      call mpi_file_open(mpi_comm_world, path_sav_p,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, npre, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(npre,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ---- write data
c
      nsize   = int(nnx,k8)*nny
      nsize2  = int(nnx,k8)*(iys -1)
      n_write = nnx*(iye+1-iys)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            temp(i,j) = p(i,j,k)
         enddo
         enddo
#if defined(SWAP)
      call byteswap(temp)
#endif
         offset = int((k-1),k8)*nsize + nsize2
         call mpi_file_write_at_all(npre,offset,temp,n_write,
     +                              mpi_real8,status,ierr)
c         call mpi_file_write_at(npre,offset,temp,n_write,
c     +                              mpi_real8,status,ierr)
      enddo
c
c ---- close file
c
      call mpi_file_close(npre, ierr)
c
c ---- check file
c
      if (l_root) then
         inquire(file=path_sav_p,exist=there)
         if(.not.there) then
            write(6,8000) path_sav_p
            call mpi_finalize(ierr)
            stop
         endif
         write(6,7000) path_sav_p
      endif
c
      return
c -------------------- process write errors
 9991 continue
      write(6,6000) npre, iz
 6000 format(' SR. SAVE_P:',/,
     +       '    trouble cannot write pressure file on unit = ',i2,/,
     +       '             at iz = ',i4)
      call mpi_finalize(ierr)
      stop
c -----------------------
 7000 format('      PRESSURE DATA IS WRITTEN IN FILE  ',a80)
 8000 format(' SR. SAVE_P: Trouble pressure file not in path =',a80)
      end
