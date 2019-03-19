
      subroutine probe_locator(height)
      use pars
      use fields
      use con_stats
      
c     Identifying probe in physical space
      do iz=0,maxnz
        if(z(iz) .lt. height)then
          iploc = iz
        end if
      enddo
      end
