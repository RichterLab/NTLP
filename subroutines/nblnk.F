      subroutine nblnk(word)
      parameter (nmax=304)
      character wordt*304, word*(*)
      nchar = len(word)
      if(nchar .gt. nmax) then
         write(6,6000) nchar,nmax
 6000    format(' TROUBLE, IN SR. NBLNK : NCHAR = ',i6,
     +          ' EXCEEDS NMAX = ',i6)
         stop
      endif
      jj = 0
      do j=1,nchar
         if(word(j:j) .ne. ' ') then
            jj = jj + 1
            wordt(jj:jj) = word(j:j)
         endif
         word(j:j) = ' '
      enddo
      do j=1,jj
         word(j:j) = wordt(j:j)
      enddo
c
      return
      end
