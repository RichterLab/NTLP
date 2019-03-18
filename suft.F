      subroutine suft(it)
c
c ---------- iterate for zeta = z/L using bisection method
c            either businger or large functions can be specified
c
c            isfc = 0, specified surface heat flux
c                 = 1, specified surface temperature
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
      real buf(3+nscl)
c
      parameter (iter_mo = 30, zeta_min = -6.0, zeta_max = 3.0)
c
c ---------- limiting value for wind
c
      free = 0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)
c
c ---- save old utau
c
      utausv = utau
      utau2  = utau*utau
c
      iz   = 1
      izp1 = iz + 1
      izm1 = iz - 1
c
      buf(1)  = 0.0
      buf(2)  = 0.0
      buf(3)  = 0.0
      tol     = 0.01
      do iy=iys,iye
      do ix=1,nnx
         buf(1) = buf(1) + u(ix,iy,iz)
         buf(2) = buf(2) + v(ix,iy,iz)
         wind(ix,iy) = sqrt((u(ix,iy,iz)+ugal)**2
     +                    +v(ix,iy,iz)*v(ix,iy,iz))
         buf(3) = buf(3) + wind(ix,iy)
      enddo
      enddo
      do iscl=1,nscl
         buf(3+iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            buf(3+iscl) = buf(3+iscl) + t(ix,iy,iscl,iz)
         enddo
         enddo
      enddo
c
c -------- get x-y slab sums
c
      call mpi_sum_xy(buf,myid,iss,ise,(3+nscl))
      u1xy  = buf(1)*fnxy + ugal
      v1xy  = buf(2)*fnxy
      windm = buf(3)*fnxy
      do iscl=1,nscl
         t1xy(iscl) = buf(3+iscl)*fnxy
      enddo
      vsfc  = sqrt(u1xy*u1xy+v1xy*v1xy)
      windm = amax1(windm,ufree)
      vsfc  = amax1(vsfc,ufree)
c
c ---------- limits for zeta
c
      zeta_mn = zeta_min
      zeta_mx = zeta_max
      if(isfc(1) .eq. 0) then
         f_con = z1*batag*vk*qstar(1)/((windm*vk)**3)
      else
         d_theta = vk74in*(tsfcc(1) - t1xy(1))
         f_con   = z1*batag*vk*d_theta/((windm*vk)**2)
      endif
c
c --------- iteration for zeta
c
      do iter=1,iter_mo
         zeta_a = 0.5*(zeta_mn + zeta_mx)
         if(ismlt .eq. 1) then
             call busngr(zeta_a,phim,phis,psim,psis)
         else
             call fzol(zeta_a,phim,phis,psim,psis)
         endif
         u_fac = (zody - psim)
         if(isfc(1) .eq. 0) then
            f_new =  zeta_a + f_con*u_fac**3
         else
            t_fac = 1.0/(zody - psis)
            f_new =  zeta_a + f_con*u_fac*u_fac*t_fac
         endif
         if(f_new .lt. 0.0) then
            zeta_mn = zeta_a
         else
            zeta_mx = zeta_a
         endif
c
c ----------- iteration details
c
c        utau      = windm*vk/(zody-psim)
c        write(nprt,1000) iter, zeta_a, utau, phim, psim
c1000    format(' 1000 iter = ',i5,' zeta = ',e15.6,' u_* = ',e15.6,
c    +          ' phim = ',e15.6,' psim = ',e15.6)
      enddo
c
c --------- check if neutral surface layer
c
      if (ibuoy.eq.0 .or.istbl.eq. 0.) then
          amonin    = 1000.
          zeta      = 0.0
          utau      = windm*vk/zody
          thstar(1) = 0.0
          t10xy(1)  = 0.0
          tsfcc(1)  = t1xy(1)
          dnom = zody
      else
         utau = windm*vk/(zody-psim)
         dnom = (zody-psis)*vk74in
         if(isfc(1) .eq. 0) then
            thstar(1) = -qstar(1)/utau
            tsfcc(1)  = t1xy(1)-thstar(1)*dnom
            t10xy(1)  = thstar(1)*dnom
         else
            thstar(1) = (t1xy(1) - tsfcc(1))/dnom
            t10xy(1)  = thstar(1)*dnom
            qstar(1)  = -utau*thstar(1)
         endif
         amonin = -utau**3/(batagk*qstar(1))
         zeta   = z1/amonin
      endif
c
c ------- surface details, for debug
c
c     write(nprt,2000) windm, utau, qstar(1), tsfcc(1), amonin, zeta,
c    +              z1, batag, vk, batagk, zo
 2000 format(' 2000 suft ',/,
     +       '    windm = ',e15.6,' utau = ',e15.6,' qstar = ',e15.6,/,
     +       '    tsfcc = ',e15.6,' MO L = ',e15.6,' z1/L = ',e15.6,/,
     +       '    z1 = ',e15.6,' batag = ',e15.6,' vk = ',e15.6,/,
     +       '    batagk = ',e15.6,' zo = ',e15.6)
c
      if (utau.gt.10.0) then
         write(6,9000)
         write(6,9200) utau,windm
         go to 9999
      endif
      if (t10xy(1).gt.0. .and. qstar(1) .gt. 0.) then
         write(6,9000)
         write(6,9300) u1xy,v1xy,t1xy(1),
     +                 tsfcc(1),amonin,utau,it
         go to 9999
      endif

!Now that temperature and momentum have been done, apply BCs on remaining scalars:

      iscl = 2
      if(isfc(iscl) .eq. 0) then
         thstar(iscl) = -qstar(iscl)/utau
         tsfcc(iscl)  = t1xy(iscl)-thstar(iscl)*dnom
         t10xy(iscl)  = thstar(iscl)*dnom
      else
!----------------------------------------------------------------------------------------------
         tsfcc(2) = 1.0* Mw/Ru/tsfcc(1)*610.94*
     +               EXP(17.6257*(tsfcc(1)-273.15)/(tsfcc(1)-
     +               273.15+243.04))/rhoa
!this was inserted to set the bottom boundarylayer q value to saturated RH  -- assumes 100% RH!
!----------------------------------------------------------------------------------------------	
         thstar(iscl) = (t1xy(iscl) - tsfcc(iscl))/dnom
         t10xy(iscl)  = thstar(iscl)*dnom
         qstar(iscl)  = -utau*thstar(iscl)
      endif


c---------John Park Sept 2017 Added David's Eulerian Conc.
c---------Computing BCs for spray concentration
      iscl = 3

      if (isfc(iscl) .eq. 0) then
         alpha = ws/vk/utau
         thstar(iscl) = -qstar(iscl)/utau
         tsfcc(iscl) = (t1xy(iscl) + qstar(iscl)/ws -
     +   qstar(iscl)/ws*(z1/zo)**(-alpha))*(z1/zo)**alpha
         t10xy(iscl) = thstar(iscl)*dnom
      else

         alpha = ws/vk/utau
         qstar(iscl) = ws*(t1xy(iscl) - tsfcc(iscl)*(z1/zo)**(-alpha))/
     +   ((z1/zo)**(-alpha) - 1.0)
         thstar(iscl) = -qstar(iscl)/utau
         t10xy(iscl)  = thstar(iscl)*dnom
      end if




c ---------- examples of two other scalars
c
c     c
c     c **** get flux of b scalar, specified surface value
c     c
c           dnom      = (zody-psis)*vk74in
c           thstar(2) = (t1xy(2)-tsfcc(2))/dnom
c           qstar(2)  = -thstar(2)*utau
c           t10xy(2)  = thstar(2)*dnom
c           aut3m(2)  =  qstar(2)
c
c **** get surface value of c scalar, specified surface flux
c
c     dnom      = (zody-psis)*vk74in
c     thstar(2) = -qstar(2)/utau
c     tsfcc(2)  = t1xy(2) - dnom*thstar(2)
c     t10xy(2)  = thstar(2)*dnom
c     aut3m(2)  = qstar(2)
c
      zol = zeta
      hol = zol*zi/z1
c
c ---- note roundoff problem in angles if close to multiples of pi
c
      tep = u1xy/vsfc
      if(tep.gt.1.)  tep = 1.0
      if(tep.lt.-1.) tep = -1.0
      thta      = acos(tep)
      utau2     = utau*utau
      au13m     = -utau2*cos(thta)
      au23m     = -utau2*sin(thta)*sign(1.,v1xy)
      do iscl = 1,nscl
      aut3m(iscl)  =  qstar(iscl)
      end do
c
      return
c
c -------- iteration did not converge
c
 9999 continue
 9000 format(' Trouble in SR. suft')
 9200 format(' Stop because utau = ',e15.6,' windm = ',e15.6)
 9300 format(' ** CHECK SFC U = ',e15.6,' V=',e15.6,' T,TS = ',2e15.6,
     +       ' L =',e15.6,' U_* = ',e15.6,' AT IT = ',i5)
      call mpi_finalize(ierr)
      stop
      end
