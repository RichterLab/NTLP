      subroutine get_units
c
      use pars
c
c -------------- unit numbers for files
c
      nvel  = 20 
      npre  = 30
      nhis1 = 40
      nvelc = 50
      nhisp = 60
      nprobe1 = 61 !SET SOMEWHERE ELSE
      nprobe2 = 62
      nprobe3 = 63
      nprobe4 = 64
      nprobe5 = 65
      nviz_z = 80
      nviz_y = 82
      nviz_x = 84
      nviz_s = 90
      nparts = 100 
c
c ------------- unit number for standard print out
c               for each mpi task
c
      nprt = 1 
c
c ------------- open unit for standard printout
c
      path_prt = case_inp(1:3)//'.mp.xxxxx.out'
      write(path_prt(8:12),'(i5.5)') myid
      open(nprt,file=path_prt,form='formatted')
c
      return
      end
