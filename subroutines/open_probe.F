

      subroutine open_probe(istep)
      use pars
      include 'mpif.h'
      character cgrid*4, iblks*16
      logical there

c
c ------------------ build character strings for ascii probe file name
c
      cgrid = '.le.'
      call blnk(iblks)
      write(iblks(1:7),'(i7.7)') istep
      iblks(8:8) = '_'
      write(iblks(9:15),'(i7.7)') (istep+itape)
      iblnk = index(path_probe,' ')
      call blnk(path_sav_probe)
      path_sav_probe = path_probe(1:iblnk-1)//'/probe'//
     +             cgrid(1:4)//case(1:3)//'.'//iblks(1:15)

c------------------- save data in directory
      if(l_root) then
      close(nprobe4)
      open(nprobe4,err=5000,file=path_sav_probe,form='formatted')
      endif

      return
 5000 continue
      write(6,6367) nprobe4, path_sav_probe
 6367 format(' 6367, SR. OPEN_PROBE:',/,
     +       '    cannot open probe history profile on unit = ',i2,/,
     +       '    path = ',a80)  
      stop
      end
