      subroutine get_output_filenames
c
c ----------- build file names for velocity, pressure, and constants
c
      use pars
      include 'mpif.h'
      character cgrid*10, num*3
c
c --------------- build character strings for file name
c
      cgrid = '.le.'
      write(num,'(i3.3)') itn
      iblnk = index(path_sav,' ')
      call blnk(path_sav_v)
      call blnk(path_sav_p)
      call blnk(path_sav_c)
      call blnk(path_sav_part)

c     JP Probe 10/29/2018
c      call blnk(path_sav_probe)

      path_sav_v = path_sav(1:iblnk-1)//'/u'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)
      path_sav_p = path_sav(1:iblnk-1)//'/p'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)
      path_sav_c = path_sav(1:iblnk-1)//'/u'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)//'.con'
      path_sav_part = path_sav(1:iblnk-1)//'/part'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)
      path_sav_hist = path_sav(1:iblnk-1)//'/hist'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)

c     JP Probe 10/29/2018
c      path_sav_probe = path_sav(1:iblnk-1)//'/probe'//
c     +                 cgrid(1:4)//case(1:3)//num(1:3)


c
      return
      end
