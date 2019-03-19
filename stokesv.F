      subroutine stokesv
c
c ----------- get stokes drift velocity for assumed wavelength stokesw
c             and wave amplitude stokesa. Changed sign for z.
c
c
      use pars
      use con_data
      use con_stats
      include 'mpif.h'
c
      if(iocean .eq. 1) then
c
c ----------- compute stokes velocity for ocean pbls
c
c        stokesw = pi2/20.0
         stokesw = pi2/76.5
c        ak      = 0.04
         ak      = 0.00
c        stokesa = 1.0
         stokesa = ak/stokesw
         sigma = sqrt(abs(grav)*stokesw)
         stokess = sigma*stokesw*stokesa**2
         do iz=1,nnzp1
            stokes(iz) = stokess*exp(2.0*stokesw*zz(iz))
         enddo
         if(l_root) then
            write(6,6000) (iz,zz(iz),stokes(iz),iz=1,nnz)
 6000       format(' iz ',10x,' zz',10x,' stokes',/,(1x,i3,2e12.4))
         endif
c
      else
c
c ----------------- set stokes velocity = 0 for atmos. pbls 
c
         do iz=1,nnzp1
            stokes(iz) = 0.0
         enddo
         stokess = 0.0
         udrift = 0.0
         vdrift = 0.0
      endif
c
      return
      end
