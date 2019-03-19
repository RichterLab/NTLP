      subroutine upper
c
c ---- set boundary condition on upper boundary iz=nnz
c      option for special radiation boundary condition
c                 index f(.,.,1)  indicates upper. 
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      iz   = nnz
      izm1 = iz - 1
      izm2 = iz - 2
      izp1 = iz + 1
      izp2 = iz + 2
c
      if(ibcu .eq. 0) then
c
c --------- boundary conditions are gradient conditions
c
c            dudzbc = 0.0
c            dvdzbc = 0.0
c            dtdzbc = dtdzf
c            wbc    = 0.0
c            ebc    = 0.0
c
        do iy=iys,iye
        do ix=1,nnx
           wbc(ix,iy,1) = 0.0
           ebc(ix,iy,1) = 0.0
           ubc(ix,iy,1) = u(ix,iy,iz)
           vbc(ix,iy,1) = v(ix,iy,iz)
           pbc(ix,iy,1) = 0.0
           pbc2(ix,iy,1)= 0.0
        enddo
        enddo
        do iscl=1,nscl
c
c ---------- first get average scalar gradient
c
           dtdzf(iscl) = 0.0
           do iy=iys,iye
           do ix=1,nnx
              dtdzf(iscl) = dtdzf(iscl) + (t(ix,iy,iscl,nnz) -
     +                      t(ix,iy,iscl,nnz-1))*dzu_i(nnz)
           enddo
           enddo
           dtdzf(iscl) = dtdzf(iscl)*fnxy
        enddo
c
        call mpi_sum_xy(dtdzf,myid,iss,ise,nscl)
c
        do iscl=1,nscl
           do iy=iys,iye
           do ix=1,nnx
              tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) + 
     +                            dtdzf(iscl)*dzu(nnzp1)
           enddo
           enddo
        enddo
      else if(ibcu .eq. 1) then
c
c ------------- special if iradup boundary condition
c               get estimate of w from continuity and 
c               linearized relation for pressure
c
      xmeanp = 0.0
      grad_ug = ug(nnz) - ug((nnz-1))
      do iy=iys,iye
      do ix=1,nnx
         wbc(ix,iy,1) = w(ix,iy,izm1)-
     +                  (ux(ix,iy,iz)+vy(ix,iy,iz))*dzw(iz)
         pbc(ix,iy,1) = .5*(w(ix,iy,izm1)+wbc(ix,iy,1))
         ebc(ix,iy,1) = 0.0
         ubc(ix,iy,1) = u(ix,iy,iz) + grad_ug
         vbc(ix,iy,1) = v(ix,iy,iz)
         pbc2(ix,iy,1)=0.5*(u(ix,iy,iz)**2 + v(ix,iy,iz)**2) +
     +              0.25*(w(ix,iy,izm1)**2 + wbc(ix,iy,1)**2)
         xmeanp = xmeanp + pbc2(ix,iy,1)
      enddo
      enddo
      call mpi_sum_xy(xmeanp,myid,iss,ise,1)
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) + 
     +                          dtdzf(iscl)*dzu(nnzp1)
         enddo
         enddo
      enddo
      xmeanp = xmeanp*fnxy
      do iy=iys,iye
      do ix=1,nnx
         pbc2(ix,iy,1) = pbc2(ix,iy,1) - xmeanp
      enddo
      enddo


      else if(ibcu .eq. 2) then
c
c ----------Constant flux at the top (fix RH and T)
c
       do iy=iys,iye
        do ix=1,nnx
           wbc(ix,iy,1) = 0.0
           ebc(ix,iy,1) = 0.0
           ubc(ix,iy,1) = u(ix,iy,iz)
           vbc(ix,iy,1) = v(ix,iy,iz)
           pbc(ix,iy,1) = 0.0
           pbc2(ix,iy,1)= 0.0
        enddo
       enddo
c
c ---------- first set the values for Temperature and RH
c
              tbc_top(1) = 310.0 
              T_real = tbc_top(1)*(1.0-9.81*z(iz)/(Cpa*tsfcc(1)))    ! real temperature
       
             tbc_top(2)=  0.05*Mw/Ru/T_real*610.94*
     +               EXP(17.6257*(T_real-273.15)/(T_real-
     +               273.15+243.04))/rhoa  ! RH set to 5%  Indrajith 08/17

         do iscl=1,nscl
           do iy=iys,iye
           do ix=1,nnx
              tbc(ix,iy,iscl,1) = tbc_top(iscl)
c              if(iy.eq.32) then
c                write(*,*) 'top-boundary',iscl,tbc(ix,iy,iscl,1)
c              end if
           enddo
           enddo
         enddo

c
c ---------- end if block
c
      endif
c
      do iy=iys,iye
      do ix=1,nnx
         w(ix,iy,iz)   = wbc(ix,iy,1)
         e(ix,iy,iz)   = ebc(ix,iy,1)
         r3(ix,iy,iz)  = 0.0
         r5(ix,iy,iz)  = 0.0
         u(ix,iy,izp1) = ubc(ix,iy,1)
         v(ix,iy,izp1) = vbc(ix,iy,1)
c ------------- note w and e nnz+1 values are not needed
         w(ix,iy,izp1) = wbc(ix,iy,1)
         e(ix,iy,izp1) = ebc(ix,iy,1)
         r3(ix,iy,izp1)= 0.0
         r5(ix,iy,izp1)= 0.0
c
c ---------- set derivatives at top of box (wx,wy not needed)
c            ux,uy,vx,vy are used in e production, but neglect
c            at top of box becuase of bc
c
         wx(ix,iy,izp1) = 0.0
         wy(ix,iy,izp1) = 0.0
         ux(ix,iy,izp1) = 0.0
         uy(ix,iy,izp1) = 0.0
         vx(ix,iy,izp1) = 0.0
         vy(ix,iy,izp1) = 0.0
c        ux(ix,iy,izp1) = ubc(ix,iy,1)
c        uy(ix,iy,izp1) = ubc(ix,iy,1)
c        vx(ix,iy,izp1) = vbc(ix,iy,1)
c        vy(ix,iy,izp1) = vbc(ix,iy,1)
      enddo
      enddo
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,iscl,izp1) = tbc(ix,iy,iscl,1)
            t(ix,iy,iscl,izp2) = tbc(ix,iy,iscl,1)
         enddo
         enddo
      enddo
c
      return
      end
