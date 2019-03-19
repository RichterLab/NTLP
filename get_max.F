      subroutine get_max
c
c --------- routine computes max velocities as sweep through
c           the velocity field 
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      real u_send(5), u_recv(5)
c
      dx_i = 1.0/dx
      dy_i = 1.0/dy
c
      u_temp = 0.0
      v_temp = 0.0
      w_temp = 0.0
      vis_temp = 0.0
      vis_temp = 0.0
      do iz=izs,ize
c
        u_xy = 0.0
        v_xy = 0.0
        w_xy = 0.0
        vis_xym = 0.0
        vis_xys = 0.0
        vis_xy = 0.0
        do iy=iys,iye
        do ix=1,nnx
           u_xy = amax1(u_xy,abs(u(ix,iy,iz)+stokes(iz)))
           v_xy = amax1(v_xy,abs(v(ix,iy,iz)))
           w_xy = amax1(w_xy,abs(w(ix,iy,iz)))
           vis_xym = amax1(vis_xym,vis_m(ix,iy,iz))
           vis_xys = amax1(vis_xys,vis_s(ix,iy,1,iz))
           vis_xy = amax1(vis_xys,vis_xym)
        enddo
        enddo
        u_xy   = u_xy*dx_i
        v_xy   = v_xy*dy_i
        wsav   = w_xy
        w_xy   = w_xy/abs(dzw(iz))
        vis_xy = vis_xy/amin1(dx,dy,dzw(iz))**2
c
        u_temp = amax1(u_xy,u_temp)
        v_temp = amax1(v_xy,v_temp)
        w_temp = amax1(w_xy,w_temp)
        vis_temp = amax1(vis_xy,vis_temp)
c
c       if(iz .le. 15) then
c         write(6,6000) iz, wmax
c6000     format(' in get_dt iz = ',i3,' wmax = ',e15.6)
c       endif
c
      enddo
      u_send(1) = u_temp
      u_send(2) = v_temp
      u_send(3) = w_temp
      u_send(4) = wsav
      u_send(5) = vis_temp 
c
      call mpi_allreduce(u_send,u_recv,5,mpi_real8,
     +     mpi_max,mpi_comm_world,ierror)
c
      umax = u_recv(1)
      vmax = u_recv(2)
      wmax = u_recv(3)
      wabs = u_recv(4)
      vismax = u_recv(5)
c
      return
      end
