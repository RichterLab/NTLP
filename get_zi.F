      subroutine get_zi(gradmax,gradout,len,itype)
c
      use pars
      real gradmax(*), gradout(*)
c
c     write(nprt,2001) myid, len
c2001 format(' 2001 in get_zi myid = ',i4,' len = ',i8)
c     write(nprt,2002) (i,gradmax(i),gradmax(i+1),i=1,len,2)
c2002 format(' i ',5x,' grad ',5x,' location ',/,
c    +      (i5,2e15.6))
c
      do i=1,len,2
         if(gradmax(i) .gt. gradout(i)) then
              gradout(i)   = gradmax(i)
              gradout(i+1) = gradmax(i+1)
         endif
      enddo
c
      return
      end
