      subroutine comp1(istep,it)
c
c ----- 3-order runge-kutta time stepping and monotone scalar fluxes in x,y,z.
c       designed to use mpi in x & y directions.
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      parameter(js = 6, ns = 3, nstat = js + ns*nscl+1)
      real stat(0:nnz,nstat)
c
c ------ temp arrays to hold rhs from step n-1 and
c        field variables from step n 
c
      real urhs(nnx,iys:iye,izs:ize), 
     +     vrhs(nnx,iys:iye,izs:ize),
     +     wrhs(nnx,iys:iye,izs:ize),
     +     erhs(nnx,iys:iye,izs:ize),
     +     trhs(nnx,iys:iye,nscl,izs:ize)
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            urhs(ix,iy,iz) = u(ix,iy,iz) + dtzeta*r1(ix,iy,iz)
            vrhs(ix,iy,iz) = v(ix,iy,iz) + dtzeta*r2(ix,iy,iz)
            wrhs(ix,iy,iz) = w(ix,iy,iz) + dtzeta*r3(ix,iy,iz)
            erhs(ix,iy,iz) = e(ix,iy,iz) + dtzeta*r5(ix,iy,iz)
         enddo
         enddo
      enddo
      do iz=izs,ize
         do l=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            trhs(ix,iy,l,iz) = t(ix,iy,l,iz) + dtzeta*r4(ix,iy,l,iz)
         enddo
         enddo
         enddo
      enddo
c
c --------- get viscosity and rhs of (e,u,v,w)-equations
c           at next step
c
      if (iDNS .eq. 1) then
         call dns_vis
         call rhs_uvw_DNS(istep)
      else
         call tke_vis(istep)   ! WARNING: this has not been checked for iscl=2
         call rhs_uvw(istep)   ! WARNING: this has not been checked for iscl=2
      end if
c
c -------- evaluate rhs of scalar equations
c
      do l=1,nscl
         if (iDNS .eq. 1) then
            call rhs_scl_dns(istep,l,it) 
         else
            call rhs_scl(istep,l,it) 
         end if
      enddo
c
c ---------- gather stat sums on root processor
c            using mpi_reduction over all processors
c
      if(istep .eq. 1) then
c
        do j=1,nstat
        do iz=0,nnz
           stat(iz,j) = 0.0
        enddo
        enddo
        do iz=izs,ize
           stat(iz,1) = uwsb(iz)
           stat(iz,2) = vwsb(iz)
           stat(iz,3) = wwsb(iz)
!           stat(iz,4) = tr_tau(iz)
           stat(iz,5) = triz(iz)
           stat(iz,6) = shrz(iz)
           !stat(iz,7) = t_diss(iz)
        enddo
        m1 = js
        m2 = js + nscl
        m3 = js + 2*nscl
        do l=1,nscl
           do iz=izs,ize
              stat(iz,m1+l) = utsb(iz,l)
              stat(iz,m2+l) = vtsb(iz,l)
              stat(iz,m3+l) = wtsb(iz,l)
           enddo
        enddo

        !Populate the iz=0 locations if processor lies at the bottom -- only for certain quantities
        if (iss .eq. 0) then
           stat(0,1) = uwsb(0)
           stat(0,2) = vwsb(0)
           do l=1,nscl
              stat(0,m3+l) = wtsb(0,l)
           enddo
        end if

        call mpi_sum_z(stat(0,1),i_root,myid,nstat*(nnz+1),1)
        do iz=0,nnz
           uwsb(iz)   = stat(iz,1)
           vwsb(iz)   = stat(iz,2)
        enddo
        do iz=1,nnz
           wwsb(iz)   = stat(iz,3)
!           tr_tau(iz) = stat(iz,4)
           triz(iz)   = stat(iz,5)
           shrz(iz)   = stat(iz,6)
           !t_diss(iz) = stat(iz,7)
        enddo
        do l=1,nscl
           do iz=1,nnz
              utsb(iz,l) = stat(iz,m1+l)
              vtsb(iz,l) = stat(iz,m2+l)
           enddo
           do iz=0,nnz
              wtsb(iz,l) = stat(iz,m3+l)
           enddo
        enddo
        do iz=1,nnz
           buyz(iz) = batag*wtsb(iz,1)
        enddo

c
c -------- end if block
c
      endif
c
c ------- save old rhs in field variables for RK-advancement
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = urhs(ix,iy,iz)
            v(ix,iy,iz) = vrhs(ix,iy,iz)
            w(ix,iy,iz) = wrhs(ix,iy,iz)
            e(ix,iy,iz) = erhs(ix,iy,iz)
         enddo
         enddo
      enddo
      do iz=izs,ize
         do l=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,l,iz) = trhs(ix,iy,l,iz)
         enddo
         enddo
         enddo
      enddo
c
      return
      end
