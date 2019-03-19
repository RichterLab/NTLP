      subroutine blnk(word)
      character word*(*)
      nchar = len(word)
      do j=1,nchar
         word(j:j) = ' '
      enddo
c
      return
      end
