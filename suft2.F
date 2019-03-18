      subroutine suft2(u_level1,it)
c
      use pars
      use fields
      use con_data
      use con_stats
c
      real u_level1(nnx,iys:iye,2+nscl)
c
c     u_level1(.,.,1) = u
c     u_level1(.,.,2) = v
c     u_level1(.,.,3) = theta
c     u_level1(.,.,4) = more scalars
c
      tol = 0.01
      ufree=0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)
      zeta_mn = -6.0
      zeta_mn_i = 1.0/zeta_mn
      iz   = 1
c     izm1 = iz - 1
c     izp1 = iz + 1
c
c      write(nprt,3131) myid, utau, zody, vk74in, batagk,
c    +               u_level1(jxs,1,3), u_level1(jxe,1,3),
c    +               u_level1(jxs,1,1), u_level1(jxe,1,1),
c    +               u_level1(jxs,1,2), u_level1(jxe,1,2)
c3131  format(' in suft2 myid = ',i4,/,
c    +        ' utau = ',e15.6,' zody = ',e15.6,/,
c    +        ' vk74in = ',e15.6,' batagk = ',e15.6,/,
c    +        ' t(jxs) = ',e15.6,' t(jxe) = ',e15.6,/,
c    +        ' u(jxs) = ',e15.6,' u(jxe) = ',e15.6,/,
c    +        ' v(jxs) = ',e15.6,' v(jxe) = ',e15.6)
c
      do iy=iys,iye
      do ix=mxs,mxe
c
c ----------------- first guess for utau
c
      utau = .001
c
      t10xy(1) = -qstar(1)/utau*zody*vk74in
      tsfcc(1) = u_level1(ix,iy,3) - t10xy(1)
      vsfc2    = u_level1(ix,iy,1)**2 + u_level1(ix,iy,2)**2
      vsfc     = sqrt(vsfc2)
      windm    = ufree+vsfc
      utausv   = utau
      utau2    = utau*utau
      amonin   = -utau2*utau/(batagk*qstar(1))
      if(amonin.eq.0.) then
            write(6,5050) ix,iy,it,utau,amonin
 5050       format(' 5050, sr. suft2, trouble at ',/,
     +             ' ix = ',i6,'iy = ',i6,' it = ',i6,' utau = ',e15.6,
     +             ' amonin = ',e15.6)
            stop
      endif
c
c ---- for unstable, free convection pbl
c
      iter = 0
 100  continue
c
c ----------------- limit the min (-l/z) change to accmmodate stable flow
c
      zeta_i = amin1(amonin/z1,zeta_mn_i)
      zeta_a = 1.0/zeta_i
c
      if(ismlt .eq. 1) then
          call busngr(zeta_a,phim,phis,psim,psis)
      else
          call fzol(zeta_a,phim,phis,psim,psis)
      endif
      utau     = windm*vk/(zody-psim)
      thstar(1)=-qstar(1)/utau
      amonold  = amonin
      amonin   = utau*utau/(batagk*thstar(1))
      diff     = abs(amonin - amonold)
c      write(nprt,5656)iter,psim,utau,zeta,amonin,dmonin,diff
c 5656 format(' iter=',i4,' phm=',e10.3,' utau=',e10.3,
c     1      ' zeta=',e10.3,' l=',e10.3,' diff = ',e12.4)
      iter = iter+1
      if(iter.gt.10)go to 1000
      if(diff.gt.abs(tol*amonin)) go to 100
 1000 continue
c
 2000 continue
c
      if (utau.gt.10.) then
         write(6,232)utau,windm
  232    format(' stop because utau=',e15.6,' windm=',e15.6)
         stop 9999
      endif
      t10xy(1) = -qstar(1)/utau*vk74in*(zody-psis)
      t_grnd(ix,iy,1) = u_level1(ix,iy,3) - t10xy(1)
c
      zol = zeta_a
      hol = zol*zi/z1
      tep = u_level1(ix,iy,1)/windm
      if(tep.gt.1.)  tep = 1.0
      if(tep.lt.-1.) tep = -1.0
      thta  = acos(tep)
      utau2 = utau*utau
c     au13m=-utau2*cos(thta)
c     au23m=-utau2*sin(thta)*sign(1.,u_level1(ix,iy,2))
c     aut3m(1)= qstar(1)
c
      tau13m(ix,iy)   = -utau2*cos(thta)
      tau23m(ix,iy)   = -utau2*sin(thta)*sign(1.,u_level1(ix,iy,2))
      taut3m(ix,iy,1) = qstar(1)
c
c **** get surface value of c scalar, specified surface flux
c
c     dnom      = (zody-psis)*vk74in
c     thstar(2) = -qstar(2)/utau
c     tsfcc(2)  = u_level1(ix,iy,4) - dnom*thstar(2)
c     t_grnd(ix,iy,2)  = u_level1(ix,iy,4) - dnom*thstar(2)
c     t10xy(2)  = thstar(2)*dnom
c     taut3m(ix,iy,2)  = qstar(2)
c
c
c ------- end of x-y loops
c
      enddo
      enddo
c
      return
      end
