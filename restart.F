
      subroutine restart
c
c ----------- get restart file from local directory
c
      use pars
      use fields
      use con_data
      use con_stats
      character*80 path_res_c
      logical there
c
c --------------------- check if file is there
c
      inquire(file=path_res,exist=there)
      if(there) then
         if(l_root) write(6,6001) path_res
      else
         if(l_root) write(6,6005) path_res
         stop
      endif
c
c ------------------ get constant file
c
      iloc = index(path_res,' ')
      path_res_c = path_res(1:iloc-1)//'.con'
      inquire(file=path_res_c,exist=there)
      if(there) then
         if(l_root) write(6,6002) path_res_c
      else
         if(l_root) write(6,6006) path_res_c
         stop
      endif
      open(nvelc,err=200,file=path_res_c,form='unformatted',
     +        status='old')
c
      call read_res
c
      return
c ---------------------------- process errors
  100 continue
      write(6,9000) path_res, nvel
      call mpi_finalize(ierr)
      stop
c -----------------------
  200 continue
      write(6,9001) path_res_c, nvelc
      call mpi_finalize(ierr)
      stop
c -----------------------
 6001 format(' SR. RESTART: FILE READ = ',A80)
 6002 format(' SR. RESTART: CONSTANT FILE READ = ',A80)
 6005 format(' 6005, SR. RESTART: cannot find restart file = ',a80)
 6006 format(' 6005, SR. RESTART: cannot find constant file = ',a80)
 9000 format(' 9000, SR. RESTART: cannot open file =',a80,/,
     +       ' to unit number = ',i2)
 9001 format(' 9001, SR. RESTART: cannot open file =',a80,/,
     +       ' to unit number = ',i2)
      end
