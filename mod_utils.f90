! ----------------------------------------------------------------------
! mod_utils — mathematical utility functions used across the codebase
!
! Combines thermodynamic/stability helpers (formerly mod_thermo) and
! statistical/random utility functions (formerly in les.F).
!
! Thermodynamic / stability:
!   fzol          – stability functions phi/psi (Businger-Dyer)
!   mod_magnus    – saturation vapour pressure (Alduchov & Eskridge 1996)
!   exner         – Exner function T/theta
!   func_p_base   – dry-adiabatic base-state pressure profile
!   func_T_base   – dry-adiabatic base-state temperature profile
!   func_rho_base – dry-adiabatic base-state density profile
!
! Random / statistical:
!   ran2          – Numerical Recipes Park-Miller RNG
!   gasdev        – Gaussian deviate (uses ran2)
!   cdf_func_single – CDF lookup, single log-normal distribution
!   cdf_func        – CDF lookup, two-component log-normal distribution
! ----------------------------------------------------------------------
module mod_utils
  use pars, only: grav, Cpa, Rd
  implicit none
  private

  ! Thermodynamic / stability
  public :: fzol
  public :: mod_magnus, exner
  public :: func_p_base, func_T_base, func_rho_base
  ! Random / statistical
  public :: ran2, gasdev, cdf_func_single, cdf_func

contains

! ──────────────────────────────────────────────────────────────────
      subroutine fzol(zeta,phim,phis,psim,psis)
      implicit real(a-h,o-z), integer(i-n)
!        estimate the stability functions for momentum, m
!                                         and scalars,  c
!        from input of the stability parameter zeta = z/L

      data c1/5./
      data a3,b3,a4,b4/1.258,8.382,-28.862,98.9545/
      data zetam,zetas/-0.2,-1.0/
      save c1, a3, b3, a4, b4, zetam, zetas
!
      psimu(Y)  = 1.571 + 2.0*(alog(0.5*(1.0 + Y)) - atan(Y)) + &
                  alog(0.5 + 0.5*Y**2)
      psisu(Y)  = 2.0*alog(0.5 + 0.5*Y)
      psicu(Y,G)= (1.0 - G)*alog(abs(Y - 1.0)) &
                + 0.5*(G + 2.0)*alog(abs(Y**2 + Y + 1.0)) &
                - (2.0*G + 1.0) / sqrt(3.0) * &
                  atan((Y + 0.5)*2.0/sqrt(3.0))
      Xm(zol)   = (1.0 - 16.0*zol)**0.25
      Xs(zol)   = sqrt(1.0 - 16.0*zol)
      Xc(zol,f) =  abs(1.0 - f*zol)**(4.0/3.0)/(1.0 - f*zol)
!
      if(zeta.ge.0.0)       then
!                                          STABLE
      if(zeta.le.1.0) then
        phim = 1.0 + c1 * zeta
        psim = - c1 * zeta
        phis = phim
        psis = psim
                      else
!                                   use limiting form
        phim = c1 + zeta
        psim = (1.0 - c1)*(1.0 + alog(zeta) ) - zeta
        phis = phim
        psis = psim
                      endif

                            else
!                                         UNSTABLE
!                                                  momentum
       if(zeta.ge.zetam) then
         phim = 1.0 / Xm(zeta)
         psim = psimu(Xm(zeta))
                         else
!                            use convective limit for momentum
         X = (1.0 - b3/a3 * zeta)**(1.0/3.0)

         fm = a3**(-1.0/3.0)
         phim = fm / Xc(zeta,b3/a3)
         psim = psimu(Xm(zetam)) &
              + psicu(Xc(zeta,b3/a3),fm) &
              - psicu(Xc(zetam,b3/a3),fm)
                         endif

!                                         UNSTABLE scalars
       if(zeta.ge.zetas) then
         phis = 1.0/Xs(zeta)
         psis = psisu(Xs(zeta))
                         else
!                              use convective limit for scalars
         fs =   abs(a4)**(-1.0/3.0)*abs(a4)/a4
         phis = (a4 - b4*zeta)**(-1.0/3.0)
         psis = psisu(Xs(zetas)) &
              + psicu(Xc(zeta,b4/a4),fs) &
              - psicu(Xc(zetas,b4/a4),fs)
                         endif

                            endif
       return
      end subroutine fzol

  ! ──────────────────────────────────────────────────────────────────
      function mod_magnus(T)
      implicit none
      real :: mod_magnus
      !Take in T in Kelvin and return saturation vapor pressure
      !using function of Alduchov and Eskridge, 1996
      real, intent(in) :: T

      mod_magnus = 610.94 * exp((17.6257*(T-273.15))/(243.04+(T-273.15)))

      end function mod_magnus

  ! ──────────────────────────────────────────────────────────────────
      function exner(p0,p)
      implicit none
      real :: exner
      real, intent(in) :: p0, p

      !Take in the reference pressure p0 and p(z), compute exner = T/theta
      exner = (p0/p)**(-Rd/Cpa)

      end function exner

      ! The dry-adiabatic base-state profile functions:

  ! ──────────────────────────────────────────────────────────────────
      function func_p_base(p_surf,T_surf,z)
      implicit none
      real :: func_p_base
      real, intent(in) :: p_surf, T_surf, z

      func_p_base = p_surf*(1 - z*grav/Cpa/T_surf)**(Cpa/Rd)

      end function func_p_base

  ! ──────────────────────────────────────────────────────────────────
      function func_T_base(T_surf,z)
      implicit none
      real :: func_T_base
      real, intent(in) :: T_surf, z

      func_T_base = T_surf - grav/Cpa*z

      end function func_T_base

  ! ──────────────────────────────────────────────────────────────────
      function func_rho_base(p_surf,T_surf,z)
      implicit none
      real :: func_rho_base
      real, intent(in) :: p_surf, T_surf, z

      func_rho_base = p_surf*T_surf**(-Cpa/Rd)/Rd* &
               (T_surf - grav*z/Cpa)**(Cpa/Rd-1.0)

      end function func_rho_base

  ! ──────────────────────────────────────────────────────────────────
  ! Random number and statistical functions (formerly in les.F)
  ! ──────────────────────────────────────────────────────────────────

      function ran2(idum)
      integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real :: ran2, AM, EPS, RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211, &
      IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER :: idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/,iv/NTAB*0/,iy/0/

      if (idum .le. 0) then
          idum=max(-idum,1)
          idum2 = idum
          do j = NTAB+8,1,-1
             k=idum/IQ1
             idum=IA1*(idum-k*IQ1)-k*IR1
             if (idum .lt. 0) idum=idum+IM1
             if (j .le. NTAB) iv(j) = idum
          end do
          iy=iv(1)
      end if
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum .lt. 0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2 .lt. 0) idum2=idum2+IM2
      j = 1+iy/NDIV
      iy = iv(j) - idum2
      iv(j) = idum
      if (iy .lt. 1) iy = iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      end function ran2
      function gasdev(idum)
       integer :: idum
       real :: gasdev
       integer :: iset
       save iset,gset
       data iset/0/

       if (idum .lt. 0) iset=0
       if (iset .eq. 0) then
 1000    v1 = 2.*ran2(idum)-1.0
         v2 = 2.*ran2(idum)-1.0
         rsq = v1**2+v2**2
         if ( (rsq .ge. 1) .or. (rsq .eq. 0)) goto 1000
         fac = sqrt(-2.0*log(rsq)/rsq)
         gset = v1*fac
         gasdev = v2*fac
         iset = 1
       else
         gasdev = gset
         iset = 0
       end if
       return
      end function gasdev
      function cdf_func_single(d,CDF,M,S)
      implicit none
      real :: d, CDF, M, S, cdf_func_single

      cdf_func_single = 1.0/2.0* &
      (1.0 + erf((log(d) - M)/S/sqrt(2.0)))  - CDF

      end function cdf_func_single
      function cdf_func(d,CDF,factor,Ma,Sa,Mc,Sc)
      implicit none
      real :: d, CDF, factor, Ma, Sa, Mc, Sc, cdf_func

      cdf_func = 1.0/(1.0+factor)/2.0* &
      ( 1.0 + erf((log(d) - Ma)/Sa/sqrt(2.0)) + &
      factor*(1 + erf((log(d) - Mc)/Sc/sqrt(2.0))) ) - CDF

      end function cdf_func

end module mod_utils
