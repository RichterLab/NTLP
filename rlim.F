      function rlim(d1,d2,d3)
c
c ------------- Cees's kappa=1/3 scheme
c
      r = (d1-d2+1.e-100)/(d2-d3-1.e-100)
      rlim = (d2-d3)*amax1(0.,amin1(r,amin1(1./6.+1./3.*r,1.)))
c
c ------------- Cees's kappa=-1 scheme
c
c     r = (d1-d2+1.e-100)/(d2-d3-1.e-100)
c     rlim = (d2-d3)*amin1(abs(r),0.5)
c
c ------------- first order upwind
c
c     rlim = 0.0
c
c ------------- QUICK scheme
c
c     rlim = -0.25*d2 - 0.125*d3 + 0.375*d1
c
      return
      end
