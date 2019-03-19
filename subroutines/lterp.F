      subroutine lterp(n,zary,zpt,i,ip1,ratio)
c
c ---- linear interpolation for zpt in zary, where zary is 
c      monotonic increasing or decreasing function
c
      dimension zary(*)
      nm1 = n-1
      if(n.le.1) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
      endif
      if(zary(1) .lt. zary(2)) go to 1
                               go to 101
    1 continue
c
c **** monotonic increasing array
c
        if(zpt .lt. zary(1)) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
        else if(zpt .gt. zary(n)) then
          i = n
          ip1 = n
          ratio = 1.0
          go to 999
        endif
        do j=1,nm1
              if(zpt .ge. zary(j) .and.
     $           zpt .le. zary(j+1)) then
                 i = j
                 ip1 = j+1
                 ratio = (zpt - zary(i))/(zary(ip1) - zary(i))
                 go to 999
              endif
        enddo
c
c **** decreasing array
c
  101 continue 
        if(zpt .gt. zary(1)) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
        else if(zpt .lt. zary(n)) then
          i = n
          ip1 = n
          ratio = 1.0
          go to 999
        endif
        do j=1,nm1
              if(zpt .le. zary(j) .and.
     $           zpt .ge. zary(j+1)) then
                 i = j
                 ip1 = j+1
                 ratio = (zpt - zary(i))/(zary(ip1) - zary(I))
                 go to 999
              endif
        enddo
  999 continue
      return
      end
