! ----------------------------------------------------------------------
! mod_thermo — thermodynamic and surface-layer stability helper functions
!
! Pure mathematical helpers that are called by both mod_solver and
! particles.  Keeping them here avoids a circular dependency between
! those two modules.
!
! Public interface
!   fzol          – stability functions phi/psi for momentum and scalars
!   mod_magnus    – saturation vapour pressure (Alduchov & Eskridge 1996)
!   exner         – Exner function T/theta
!   func_p_base   – dry-adiabatic base-state pressure profile
!   func_T_base   – dry-adiabatic base-state temperature profile
!   func_rho_base – dry-adiabatic base-state density profile
! ----------------------------------------------------------------------
module mod_thermo
  use pars, only: grav, Cpa, Rd
  implicit none
  private

  public :: fzol
  public :: mod_magnus, exner
  public :: func_p_base, func_T_base, func_rho_base

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

end module mod_thermo
