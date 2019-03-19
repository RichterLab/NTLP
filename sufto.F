      subroutine sufto(it)
c
      use pars
      use fields
      use con_data
      use con_stats
      real buf(3+nscl)
c
c ------- version of similarity theory adpated for ocean flows
c      option to use businger or large version of similarity theory
c
      iz    = 1
      izm1  = iz - 1
      izp1  = iz + 1
      z1_a  = abs(z1)
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
      t10xy(1)=-qstar(1)/utau*zody*vk74in
c
c ---- check for temperature boundary condition
c
      if(isfc(iscl) .eq. 0 ) then
         tsfcc(1)=t1xy(1)-t10xy(1)
      endif
c     vsfc=sqrt(u1xy*u1xy+v1xy*v1xy)
c     if(windm.le.0.01)windm=0.01
c     if(vsfc .le.0.01)vsfc =0.01
c
c ----------- input surface wind stress (tau = 0.0184n/m*m)
c             density rho = 1000kg/m^3
c          
c     utau = 4.29e-03
c     utau = 6.10e-03
      utau = 7.00e-03
c
c **** save old utau
      utausv = utau
      utau2  = utau*utau
c      if (ibuoy.eq.0 .or. qstar(1) .eq. 0.) then
        if (ibuoy.eq.0 .or. istbl .eq. 0.) then
          amonin    = 1000.
          zeta      = 0.
          thstar(1) = 0.0
          t10xy(1)  = 0.0
      else
          amonin = -utau2*utau/(batagk*qstar(1))
          zeta   = z1_a/amonin
      endif
      if (t10xy(1).lt.0. .and. qstar(1) .lt. 0.) then
         write(6,1234)u1xy,v1xy,t1xy(1),tsfcc(1),amonin,utau,it
 1234    format(' ** check sfc u=',e12.3,' v=',e12.3,' t,ts=',2f10.3,
     +     ' l=',e12.3,' u*=',e12.3,' at it=',i5)
         go to 9999
      endif
c
c -------- for stable,neutral and unstable pbl get drift velocity
c
      if(ismlt .eq. 1) then
          call busngr(zeta,phim,phis,psim,psis)
      else
          call fzol(zeta,phim,phis,psim,psis)
      endif
      udrift = windm + stokes(1) - stokess + utau*(zody-psim)*vkin
      vdrift = 0.0
      dnom      = (zody-psis)*vk74in
      if (isfc(iscl).eq.1) then
         thstar(1) = (t1xy(1) - tsfcc(1))/dnom
         t10xy(1)  = thstar(1)*dnom
         qstar(1)  = - utau*thstar(1)
      else
         thstar(1)  = -qstar(1)/utau
         tsfcc(1)   = t1xy(1)-thstar(1)*dnom
         t10xy(1)   = thstar(1)*dnom
      endif
      zol = zeta
      hol = zol*zi/z1
c
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
c     c
c     c **** get surface value of c scalar, specified surface flux
c     c
c           dnom      = (zody-psis)*vk74in
c           thstar(3) = -qstar(3)/utau
c           tsfcc(3)  = t1xy(3) - dnom*thstar(3)
c           t10xy(3)  = thstar(3)*dnom
c           aut3m(3)  = qstar(3)
c
c **** note roundoff problem in angles are close to multiples of pi
c     tep=u1xy/vsfc
c     if(tep.gt.1.)tep=1.
c     if(tep.lt.-1.)tep=-1.
c     thta=acos(tep)
      utau2 = utau*utau
c     au13m = -utau2*cos(thta)
c     au23m = -utau2*sin(thta)*sign(1.,v1xy)
      au13m = utau2
      au23m = 0.0
      aut3m(1)= qstar(1)
c
      return
c
c --------- trouble in sl routine
c
 9999 continue
c
      write(nprt,9000)
 9000 format(' Trouble in SR. sufto')
      call mpi_finalize(ierr)
      stop
      end
