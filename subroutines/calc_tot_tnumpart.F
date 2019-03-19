      





      subroutine calc_tot_tnumpart(mass_fraction,tot_parts)
      use pars
      use fields
      use con_data
      use particles
      implicit none
      real, intent(in)     :: mass_fraction
      integer, intent(out) :: tot_parts
      real :: tot_numflux,vol_frac,horiz_area,domain_vol
      real :: mass_flux,vol_flux,pi
      integer :: i,j,k,l,m,n

      pi = 4.0*atan(1.0)
      !Domain horizontal area and volume
      horiz_area = xl*yl
      domain_vol = xl*yl*zl

      !Get total number flux based on unscaled source density
      !using numerical integration.
      
      n=size(Snscl)
      if(myid==0) write(*,*)'size(Snscl): ',n
      tot_numflux = Snscl(1)
      vol_flux = 4/3*pi*rpcrop(1)*Snscl(1)
      mass_flux = (4/3*pi*rpcrop(1)**3*rhow)*Snscl(1)
      do i = 2,n
        tot_numflux = tot_numflux + Snscl(i) !# particles/m^2/s
        vol_flux    = vol_flux    + 4/3*pi*rpcrop(i)**3*Snscl(i)
        mass_flux   = mass_flux   + (4/3*pi*rpcrop(i)**3*rhow)*Snscl(i)
      end do

      tot_parts = int(mass_fraction*(rhoa*domain_vol)
     +            /(mass_flux/tot_numflux))

      if (myid==0) write(*,*)'total number flux: ', tot_numflux
      if (myid==0) write(*,*)'volume flux      : ', vol_flux
      if (myid==0) write(*,*)'mass flux        : ', mass_flux
      if (myid==0) write(*,*)'domain volume    : ', domain_vol
      if (myid==0) write(*,*)'mass fraction    : ', mass_fraction
      if (myid==0) write(*,*)'total number flux: ', tot_numflux

      end
