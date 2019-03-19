      subroutine lower_dns

      !Make lower BC by setting w = 0, u,v equal to mirror of interior points
      !Also set the scalar boundary condition based on either Neumann or Dirichlet

      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      implicit none

      integer :: iz,izm1,iscl,ix,iy
      real :: dz_i,dtdz,flux

      iz   = 1
      izm1 = iz - 1
      dz_i = dzu_i(1)

c ------------ initialize u, v, w, t and derivatives at izm1
c
      do iy=iys,iye
      do ix=1,nnx
         u(ix,iy,izm1)  = -2.0*Uo-u(ix,iy,iz)
         v(ix,iy,izm1)  = -v(ix,iy,iz) 
         w(ix,iy,izm1)  = 0.0 
         r3(ix,iy,izm1) = 0.0
         e(ix,iy,izm1)  = 0.0 
         ux(ix,iy,izm1) = 0.0
         uy(ix,iy,izm1) = 0.0
         vx(ix,iy,izm1) = 0.0
         vy(ix,iy,izm1) = 0.0
         wx(ix,iy,izm1) = 0.0 
         wy(ix,iy,izm1) = 0.0 
      enddo
      enddo


!NOTE: sign convention is that wtsfc is flux upwards INTO domain (not just in vertical direction)
! Set the scalar boundary condition based on isfc
c            isfc = 0, specified surface heat flux (through qstar)
c                 = 1, specified surface temperature (Tbot)

      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            if (isfc(iscl)==1) then
            t(ix,iy,iscl,izm1) = 2.0*Tbot(iscl) - t(ix,iy,iscl,iz) 
            end if

            if (isfc(iscl)==0) then 
            flux = dzu(0)*wtsfc(iscl)/vis_s(ix,iy,iscl,izm1)
            t(ix,iy,iscl,izm1) = t(ix,iy,iscl,iz)     
     +       + flux  
            end if
         enddo
         enddo
      enddo

      end subroutine lower_dns
      subroutine upper_dns

      !Make upper BC by setting w = 0, u,v equal to mirror of interior points
      !Also set the scalar boundary condition based on either Neumann or Dirichlet

      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      implicit none

      integer :: iz,izp1,iscl,ix,iy
      real :: dz_i,dtdz,flux

      iz   = nnz
      izp1 = iz + 1

c
c ------------ initialize u, v, w, t and derivatives at izp1
c
      do iy=iys,iye
      do ix=1,nnx
         u(ix,iy,izp1)  = 2.0*Uo-u(ix,iy,iz)
         v(ix,iy,izp1)  = -v(ix,iy,iz) 
         w(ix,iy,izp1)  = 0.0 
         r3(ix,iy,izp1) = 0.0
         e(ix,iy,izp1)  = 0.0 
         ux(ix,iy,izp1) = 0.0
         uy(ix,iy,izp1) = 0.0
         vx(ix,iy,izp1) = 0.0
         vy(ix,iy,izp1) = 0.0
         wx(ix,iy,izp1) = 0.0 
         wy(ix,iy,izp1) = 0.0 
      enddo
      enddo

! Set the scalar boundary condition based on isfc
c            isfc = 0, specified surface heat flux (through qstar)
c                 = 1, specified surface temperature (Ttop)

!NOTE: sign convention is that wtsfc is flux upwards INTO domain (not just in vertical direction)

      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            if (isfc(iscl)==1) then
            t(ix,iy,iscl,izp1) = 2.0*Ttop(iscl)-t(ix,iy,iscl,iz) 
            end if
            if (isfc(iscl)==0) then 
            flux = dzu(izp1)*wtsfc(iscl)/vis_s(ix,iy,iscl,izp1)
            t(ix,iy,iscl,izp1) = t(ix,iy,iscl,iz)      
     +       - flux 
            end if
         enddo
         enddo
      enddo

      end subroutine upper_dns
      subroutine lower_free(it)
c
c --------------- setup lower boundary condition for free
c                 convection where each processor applies
c                 log-law at several (ix,iy) for iz = 1
c
c                 index f(.,.,2)  indicates lower
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      real u_level1(nnx,iys:iye,2+nscl), buf(2+2*nscl)
      real sbuf(2+2*nscl,mxs:mxe,iys:iye)
      real rbuf((2+2*nscl)*nnx*(iye+1-iys))
c
c -------------- broadcast level 1 data everywhere
c
      if(iss .eq. 0) then
         do iy=iys,iye
         do ix=1,nnx
            u_level1(ix,iy,1) = u(ix,iy,1)
            u_level1(ix,iy,2) = v(ix,iy,1)
         enddo
         enddo
         do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            u_level1(ix,iy,2+iscl) = t(ix,iy,iscl,1)
         enddo
         enddo
         enddo
      endif
      num = nnx*(iye + 1 - iys)*(2+nscl)
c
c ------ send all of root data to other processors
c
      call mpi_send_root(u_level1(1,iys,1),
     +             num,myid,numprocs,ncpu_s)
c
c --------- every task gets their own fluxes and surface scalars
c
      call suft2(u_level1,it)
c
c --------- send surface scalars and momentum fluxes
c           back to root(s)
c
      if(numprocs .eq. 1) go to 999
c
      do iy=iys,iye
      do ix=mxs,mxe
         sbuf(1,ix,iy)  = tau13m(ix,iy)
         sbuf(2,ix,iy)  = tau23m(ix,iy)
      enddo
      enddo
      do iscl=1,nscl
      do iy=iys,iye
      do ix=mxs,mxe
         sbuf(2+iscl,ix,iy)      = taut3m(ix,iy,iscl)
         sbuf(2+nscl+iscl,ix,iy) = t_grnd(ix,iy,iscl)
      enddo
      enddo
      enddo
c
      irow_r = mod(myid,ncpu_s)
      if(myid .ge. ncpu_s) then
        num = (2+2*nscl)*(mxe+1-mxs)*(iye+1-iys)
        call mpi_send(sbuf(1,mxs,iys),num,mpi_real8,irow_r,1,
     +       mpi_comm_world,ierr)
      else
        do l=irow_r+ncpu_s,numprocs-1,ncpu_s
           num = (2+2*nscl)*(mx_e(l)+1-mx_s(l))*(iye+1-iys)
           call mpi_recv(rbuf(1),num,mpi_real8,l,1,
     +          mpi_comm_world,istatus,ierr)
c          call f_suft2(rbuf,maxnx,maxny,mx_s(l),mx_e(l),iys,iye,nscl,
           call f_suft2(rbuf,nnx,mx_s(l),mx_e(l),iys,iye,nscl,
     +                  tau13m,tau23m,taut3m,t_grnd)
        enddo
      endif
c
  999 continue
c
c ------------ only for root row = 0
c              get sums of surface conditions
c              and set surface boundary conditions
c
      if(iss .eq. 0) then
c
         buf(1) = 0.0
         buf(2) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            buf(1) = buf(1) + tau13m(ix,iy)
            buf(2) = buf(2) + tau23m(ix,iy)
         enddo
         enddo
         do iscl=1,nscl
            buf(2+iscl)      = 0.
            buf(2+nscl+iscl) = 0.
            do iy=iys,iye
            do ix=1,nnx
               buf(2+iscl)      = buf(2+iscl) + taut3m(ix,iy,iscl)
               buf(2+nscl+iscl) = buf(2+nscl+iscl) + t_grnd(ix,iy,iscl)
            enddo
            enddo
         enddo
c
         call mpi_sum_xy(buf,myid,iss,ise,2+2*nscl)
         uwsfc = buf(1)*fnxy
         vwsfc = buf(2)*fnxy
         do iscl=1,nscl
            wtsfc(iscl) = buf(2+iscl)*fnxy
            tsfcc(iscl) = buf(2+nscl+iscl)*fnxy
         enddo
c
         iz   = 1
         izm1 = iz - 1
         dz_i = dzu_i(iz)
c
         do iy=iys,iye
         do ix=1,nnx
            ebc(ix,iy,2)=amax1(e(ix,iy,iz),sml_eg)
            wbc(ix,iy,2)= 0.0
            pbc(ix,iy,2) = 0.0
            pbc2(ix,iy,2) = 0.0
         enddo
         enddo
c
         do iy=iys,iye
         do ix=1,nnx
            dudz     = 2.*u(ix,iy,iz)*dz_i
            dvdz     = 2.*v(ix,iy,iz)*dz_i
            ubc(ix,iy,2) = u(ix,iy,iz) - dudz*dzu(iz)
            vbc(ix,iy,2) = v(ix,iy,iz) - dvdz*dzu(iz)
         enddo
         enddo
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               dtdz     = 2.*(t(ix,iy,iscl,iz)-tsfcc(iscl))*dz_i
               tbc(ix,iy,iscl,2) = t(ix,iy,iscl,iz) - dtdz*dzu(iz)
            enddo
            enddo
         enddo
c
c ------------ initialize u, v, w, t and derivatives at izm1
c
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,izm1)  = ubc(ix,iy,2)
            v(ix,iy,izm1)  = vbc(ix,iy,2)
            w(ix,iy,izm1)  = wbc(ix,iy,2)
            r3(ix,iy,izm1) =  0.0
            e(ix,iy,izm1)  = ebc(ix,iy,2)
            ux(ix,iy,izm1) = 0.0
            uy(ix,iy,izm1) = 0.0
            vx(ix,iy,izm1) = 0.0
            vy(ix,iy,izm1) = 0.0
            wx(ix,iy,izm1) = wbc(ix,iy,2)
            wy(ix,iy,izm1) = wbc(ix,iy,2)
         enddo
         enddo
c
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izm1) = tbc(ix,iy,iscl,2)
            enddo
            enddo
         enddo
c
c ----- end of if block for root row
c
      endif
c
 7999 continue
c
      return
      end
