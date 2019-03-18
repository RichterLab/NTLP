

      !Generates a random particle given the CDF from a uniform
      !distribution
      subroutine create_random_partsize(rad_rand)
      use pars
      use fields
      use con_data
      use con_stats
      use particles
      
      real,intent(out) :: rad_rand
      real             :: ran2

      rad_rand = rpcrop(minloc(CDF,DIM=1,MASK=(CDF>=ran2(iseed))))
      end
