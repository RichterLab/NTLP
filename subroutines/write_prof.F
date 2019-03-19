      subroutine write_prof(nhisp,krec,num,f)
      real f(num)
      real*4 f32(num)
c
c -------------- build special 32 bit arrays for profiles
c
      do i=1,num
         f32(i) = f(i)
      enddo
c
      write(nhisp,err=999,rec=krec) (f32(i),i=1,num)
c
      return
c --------------- errors
  999 continue
      write(6,9000) num,krec
 9000 format(' 9000, trouble in ',
     +       'SR. save_prof cannot write profile data ',/,
     +       ' num = ',i8, 'krec = ',i6)
      stop
      end
