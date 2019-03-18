      subroutine busngr(zeta,phim,phis,psim,psis)
c
c ---- Businger's version of similarity theory
c
      data pih /1.57079633/
      save pih
c
      if(zeta .lt. 0.) then
         x=(1.0 - 15.0*zeta)**0.25
         phim = 1.0/x
         psim = 2.0*alog((1.0+x)/2.0) + alog((1.0+x*x)/2.0) - 
     +          2.0*atan(x)+pih
         if(psim.gt.2.0)psim=2.0
         y = sqrt(1.0-9.0*zeta)
         phis = 0.74/y
         psis = alog((1.0+y)/2.0)*2.0
      else if(zeta .gt. 0) then
         phim = 1.0 + 4.7*zeta
         phis = 0.74 + 4.7*zeta
         psim = -4.7*zeta
         psis = -4.7*zeta
      else
         phim = 1.0
         phis = 0.74
         psim = 0.0
         psis = 0.0
      endif
      return
      end
