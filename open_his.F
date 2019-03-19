      subroutine open_his(istep)
c
c ------------------- open history files by root
c                     isize determined in sr. fill_cs
c
      use pars
      include 'mpif.h'
      character cgrid*4, iblks*16
      logical there
c
c --------------- build character strings for ascii history file name
c
      cgrid = '.le.'
      call blnk(iblks)
      write(iblks(1:7),'(i7.7)') istep
      iblks(8:8) = '_'
      write(iblks(9:15),'(i7.7)') (istep + itape)
      iblnk = index(path_his,' ')
      call blnk(path_sav_h)
      path_sav_h = path_his(1:iblnk-1)//'/his'//
     +         cgrid(1:4)//case(1:3)//'.'//iblks(1:15)
c
c --------------- build character strings for ieee profile history file
c                 set record counter for direct access file = 0
c
      krec = 0
      cgrid = '.le.'
      call blnk(iblks)
      write(iblks(1:7),'(i7.7)') istep
      iblks(8:8) = '_'
      write(iblks(9:15),'(i7.7)') (istep + itape)
      iblnk = index(path_his,' ')
      call blnk(path_sav_hp)
      path_sav_hp = path_his(1:iblnk-1)//'/his'//
     +         cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.ieee'
c
c ----------------- save data in directory
c
      if(l_root) then

      close(nhis1)
      open(nhis1,err=3000,file=path_sav_h,form='formatted')
c
      close(nhisp)
      open(nhisp,err=4000,file=path_sav_hp,
     +        form='unformatted',access='direct',recl=isize*j_recl,
     +        status='unknown')
      endif
c
      return
c ------------------- process errors
 3000 continue
      write(6,6301) nhis1, path_sav_h
 6301 format(' 6301, SR. OPEN_HIS:',/,
     +       '    cannot open history1 file on unit = ',i2,/,
     +       '    path = ',a80)
      stop
c-------------------
 4000 continue
      write(6,6302) nhisp, path_sav_hp
 6302 format(' 6302, SR. OPEN_HIS:',/,
     +       '    cannot open history profile file on unit = ',i2,/,
     +       '    path = ',a80)
      stop
      end
