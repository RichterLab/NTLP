      subroutine viz_output_filename(istep)
c
c ------------------- set visualization files,
c                     leaves files in scratch directory 
c
      use pars
      include 'mpif.h'
      character iblks*16

c
c --------------- build character strings for file names
c                 with time step
c
      call blnk(iblks)
      iblks(1:1) = '.'
      write(iblks(2:8),'(i7.7)') istep
      iblks(9:9) = '_'
      write(iblks(10:16),'(i7.7)') (istep + itape)
c
      iloc = index(path_seed,' ')
      path_viz_xy = path_seed(1:iloc-1)
     +         //'/viz.'//case(1:3)//iblks(1:16)//'.xy.data'
c
      path_viz_xz = path_seed(1:iloc-1)
     +         //'/viz.'//case(1:3)//iblks(1:16)//'.xz.data'
c
      path_viz_yz = path_seed(1:iloc-1)
     +         //'/viz.'//case(1:3)//iblks(1:16)//'.yz.data'
c
      path_stuf = path_seed(1:iloc-1)
     +         //'/stuff.'//case(1:3)//iblks(1:16)//'.data'
c
c     if(l_root) then
c        write(6,8001) path_viz_xy
c8001    format(' 8001: viz file = ',a80)
c        write(6,8001) path_viz_xz
c        write(6,8001) path_viz_yz
c        write(6,8001) path_stuf
c        write(6,8001) path_seed
c     endif
c
      return
      end
