      subroutine surfvis(it)
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      real xkvis(nnx,iys:iye), alwk(nnx,iys:iye)
c
      real send(3), buf(3)
c
      xksurf = 0.0
      viscon = 0.0
      vise   = 0.0
c
c ----------- only root process(es) compute 
c
      if(iss .eq. 0) then

c     ck = 0.1
c     csmag = 0.18
c     xkmax  = dzdz/dt/5.
      iz   = 1
      izm1 = iz - 1
      izp1 = iz + 1
c     xkmax  = dzu(izp1)*dzu(izp1)/(5.0*dt)
      dz_i = dzu_i(izp1)
      if(iocean .eq. 1) then
         call sufto(it)
      else
         call suft(it)
      endif
      if(qstar(1) .eq. 0.) then
         zeta = 0.0
      else
         zeta = abs(z(1))/amonin
      endif
      if(ismlt .eq. 1) then
          call busngr(zeta,phim,phis,psim,psis)
      else
          call fzol(zeta,phim,phis,psim,psis)
      endif
      viscon = vk*abs(z(1))/(utau*phim)
      vise   = utau*vk*abs(z(1))/phim
c
c ---- get special value at z1 to match with surface layer
c
      uws = 0.0
      vws = 0.0
      do iy=iys,iye
      do ix=1,nnx
         uws = uws + 0.5*(u(ix,iy,iz)-u_mn(iz) + 
     +         u(ix,iy,izp1) - u_mn(izp1))*(w(ix,iy,iz)-w_mn(iz))
         vws = vws + 0.5*(v(ix,iy,iz)-v_mn(iz) + 
     +         v(ix,iy,izp1) - v_mn(izp1))*(w(ix,iy,iz)-w_mn(iz))
      enddo
      enddo
      uws = uws*fnxy
      vws = vws*fnxy !Brian corrected from fxy 1/15/15
c
c ---- get average fluctuating eddy viscsoity
c
      do iy=iys,iye
      do ix=1,nnx
         e(ix,iy,iz)=amax1(e(ix,iy,iz),sml_eg)
      enddo
      enddo
      dslk = amin1(dsl,vk*abs(z(iz))/csmag)
c     stabmin = 1.e-12
c     almin = 0.0001*dsl
      almin = almin_c*dsl_z(iz)
      do iy=iys,iye
      do ix=1,nnx
         alwk(ix,iy)=dslk
c
c --------no stability corrected length scales when
c         new eddy viscosity is on
c
c         stab=batag*(t(ix,iy,1,izp1)-t(ix,iy,1,iz))*dz_i
c         if(stab.gt.stabmin) then
c           als = stab_c*sqrt(e(ix,iy,iz)/stab)
c           alwk(ix,iy) = amin1(dslk,als)
c         endif
c         alwk(ix,iy)=amax1(almin,alwk(ix,iy))
         xkvis(ix,iy)=ck*alwk(ix,iy)*sqrt(e(ix,iy,iz))*dfac(1)
c        xkvis(ix,iy)=amin1(xkvis(ix,iy),xkmax)
      enddo
      enddo
c
c ---- get average viscosity
c
      xkavg = 0.0
      do iy=iys,iye
      do ix=1,nnx
         xkavg = xkavg + xkvis(ix,iy)
      enddo
      enddo
      xkavg = xkavg*fnxy
c
      buf(1) = uws
      buf(2) = vws
      buf(3) = xkavg
      call mpi_sum_xy(buf,myid,iss,ise,3)
c
      uws   = buf(1)
      vws   = buf(2)
      xkavg = buf(3)
c
      xkz1 = vise - sqrt(uws**2 + vws**2)*viscon
      xksurf =  xkz1 - xkavg
      xksurf = amax1(xksurf,0.0)
      xksurf = amin1(xksurf,vise)
c     if(l_root) write(6,6000) dfac(1), xkavg, xkz1, vise, xksurf
 6000 format(' dfac = ',e12.4,' xkavg = ',e12.4,' xkz1 = ',e12.4,/,
     +       ' vise = ',e12.4,' xksurf = ',e12.4)
c
      endif

c
c ---------- broadcast values to other processes
c
      send(1) = xksurf
      send(2) = viscon
      send(3) = vise
c
      call mpi_bcast(send,3,mpi_real8,
     +              i_root,mpi_comm_world,ierr)
c
      xksurf = send(1)
      viscon = send(2)
      vise   = send(3)
c
      return
      end
