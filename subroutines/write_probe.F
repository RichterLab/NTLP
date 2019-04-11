



      subroutine write_probe
      use pars
      use fields
      use con_data
      use con_stats
      use fftwk
      use particles
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'

!     Call a bunch of mpi_send receive communications

!     Get the 3 probe information


c      write(nprobe4,7072) u, v, w, partcount,
c     +     u, v, w, partcount, 
c     +     u, v, w, partcount
c 7072 format(12e17.8)
      end
