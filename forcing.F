      subroutine forcing
c
c ----------- update surface temperature based on a 
c             constant cooling rate
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      tsfcc(1) = t_surf_i - c_rate*time
c
      return
      end
