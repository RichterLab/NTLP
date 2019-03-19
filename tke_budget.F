      subroutine tke_budget
c
c -------- get terms in resolved scale tke budget
c          as in gabls writeup at w-points
c          at istage = 1. 
c
      use pars
      use particles
      use fields
      use con_data
      use con_stats
c
      real :: stat(1:nnz,5)
      real :: s11s,s22s,s33s,s12s,s13s,s23s,wz,wzp
      real :: s13,s23,s33
      real :: ufluc,ufluc_t,ufluc_b,vfluc,vfluc_t,vfluc_b,wfluc
      real :: ffluc1,ffluc2,ffluc3
      real :: ffluc1p,ffluc2p,ffluc3p
      real :: weit,weit1
      integer :: iz,i,j,izp1,izm1
c
c -------- stat(.,1) = tke transport  = wq
c          stat(.,2) = pressure transport  = wp
c          stat(.,3) = tke dissipation
c          stat(.,4) = tke dissipation
c          stat(.,5) = particle force correlation
c
      do iz=1,nnz
         stat(iz,1) = 0.0
         stat(iz,2) = 0.0
         stat(iz,3) = 0.0
         stat(iz,4) = 0.0
         stat(iz,5) = 0.0
      enddo

!Compute DNS dissipation, since there is no subgrid dissipation now:
      do iz=izs,ize
         izp1 = iz + 1
         izm1 = iz - 1
         weit = dzw(iz)/(dzw(iz) + dzw(izp1))
         weit1 = 1.0 - weit
!
! ---- get fluctuating strain rates:
!      here, sij = 1/2*(duidxj + dujdxi)
!      then t_diss = 2*nu*<sij sij>
! ---- NOTE: these are computed at the w-locations!  (not u,v locations)
!
         t_diss(iz) = 0.0
         do j=iys,iye
         do i=1,nnx

            !Things for dissipation - these are computed at w-locations since
            !there is no z-derivative after
            s11s = weit1*ux(i,j,iz)**2 + weit*ux(i,j,izp1)**2
            s22s = weit1*vy(i,j,iz)**2 + weit*vy(i,j,izp1)**2
            wz  = (w(i,j,iz)-w(i,j,izm1))*dzw_i(iz)
            wzp = (w(i,j,izp1)-w(i,j,iz))*dzw_i(izp1)
            s33s = weit*wzp**2 + weit1*wz**2
            s12s = weit1*(0.5*(uy(i,j,iz) + vx(i,j,iz)))**2 +
     +            weit*(0.5*(uy(i,j,izp1) + vx(i,j,izp1)))**2
            s13s = (0.5*((u(i,j,izp1) - u(i,j,iz) +
     +            uxym(iz) - uxym(izp1))*dzu_i(izp1) +
     +            wx(i,j,iz)))**2
            s23s = (0.5*((v(i,j,izp1) - v(i,j,iz) +
     +          vxym(iz) - vxym(izp1))*dzu_i(izp1) +
     +          wy(i,j,iz)))**2

         stat(iz,3) = stat(iz,3) + 2.0*vis_m(i,j,iz)*(
     +               s11s+s22s+s33s+2.0*(s12s+s13s+s23s))

            !Things for viscous transport - these are computed at u-locations since
            !a z-derivative is done to the average, which lands t_tau on w-locations
            
            ufluc_t   = u(i,j,izp1) - uxym(izp1)
            ufluc   = u(i,j,iz) - uxym(iz)
            ufluc_b   = u(i,j,izm1) - uxym(izm1)
            vfluc_t   = v(i,j,izp1) - vxym(izp1)
            vfluc   = v(i,j,iz) - vxym(iz)
            vfluc_b   = v(i,j,izm1) - vxym(izm1)
            wfluc = 0.5*( (w(i,j,iz)-wxym(iz)) 
     +                  + (w(i,j,izm1)-wxym(izm1)) )

            uz_t = (ufluc_t-ufluc)*dzu_i(izp1)
            uz_b = (ufluc-ufluc_b)*dzu_i(iz)
            vz_t = (vfluc_t-vfluc)*dzu_i(izp1)
            vz_b = (vfluc-vfluc_b)*dzu_i(iz)
            
            uz = 0.5*(uz_t+uz_b)
            vz = 0.5*(vz_t+vz_b)

            s13 = 0.5*(uz + 0.5*(wx(i,j,iz)+wx(i,j,izm1)))
            s23 = 0.5*(vz + 0.5*(wy(i,j,iz)+wy(i,j,izm1)))
            s33 = wz

         !Note: just uses vis_m(1,1,iz) since it's equal everywhere:
         stat(iz,4) = stat(iz,4) + 2.0*vis_m(i,j,iz)*(
     +               ufluc*s13 + vfluc*s23 + wfluc*s33)
                    
         !Finally get the particle force correlation term:
         ffluc1 = partsrc(i,j,iz,1)-m1src(iz)
         ffluc2 = partsrc(i,j,iz,2)-m2src(iz)
         ffluc3 = partsrc(i,j,iz,3)-m3src(iz)
         if (iz==nnz) then
         ffluc1p = 0.0
         ffluc2p = 0.0
         ffluc3p = 0.0
         else
         ffluc1p = partsrc(i,j,izp1,1)-m1src(izp1)
         ffluc2p = partsrc(i,j,izp1,2)-m2src(izp1)
         ffluc3p = partsrc(i,j,izp1,3)-m3src(izp1)
         end if
         stat(iz,5) = stat(iz,5) +
     +       weit*(ufluc_t*ffluc1p)+weit1*(ufluc*ffluc1)+
     +       weit*(vfluc_t*ffluc2p) + weit1*(vfluc*ffluc2) + 
     +       (w(i,j,iz)-wxym(iz))*(weit*ffluc3p+weit1*ffluc3)
         enddo
         enddo
         stat(iz,3) = stat(iz,3)*fnxy
         stat(iz,4) = stat(iz,4)*fnxy
         stat(iz,5) = stat(iz,5)*fnxy
      enddo
c
c --------------- get transport terms as vertical arrays
c
      do iz=izs,ize
c
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
c
c --------- get estimate of turbulent transport term
c
            ufluc   = u(ix,iy,iz) - uxym(iz)
            vfluc   = v(ix,iy,iz) - vxym(iz)
            wfluc   = w(ix,iy,iz) - wxym(iz)
            wfluc_l = w(ix,iy,izm1) - wxym(izm1)
            stat(iz,1)  = stat(iz,1) + 0.25*(wfluc + wfluc_l)*
     +             (ufluc**2 + vfluc**2 + 0.5*(wfluc_l**2 + wfluc**2))
c
c --------- get estimate of pressure transport term
c
            pfluc = p(ix,iy,iz) - pxym(iz)
     +           -0.5*((u(ix,iy,iz)+stokes(iz))**2 + 
     +                 v(ix,iy,iz)*v(ix,iy,iz) + 
     +      0.5*(w(ix,iy,iz)*w(ix,iy,iz)+w(ix,iy,izm1)*w(ix,iy,izm1)))
            stat(iz,2) = stat(iz,2) + pfluc*0.5*(wfluc_l + wfluc)
         enddo
         enddo
         stat(iz,1) = stat(iz,1)*fnxy
         stat(iz,2) = stat(iz,2)*fnxy
      enddo
      call mpi_sum_z(stat(1,1),i_root,myid,nnz*5,1)
c
c ------ we have all terms on all processors for all z, add them up
c        treat tr_tau at bottom special, tr_tau = 0
c
!      tr_tau(0) = 0.0
      do iz=1,nnz
c
         izp1 = iz + 1
         izm1 = iz - 1
         if(iz .eq. nnz) then
            t_tau(iz) = 0.0
            t_wp(iz)  = 0.0
            t_wq(iz)  = 0.0
         else
            t_tau(iz) = (stat(izp1,4) - stat(iz,4))*dzu_i(izp1) 
            t_wq(iz)  = -(stat(izp1,1) - stat(iz,1))*dzu_i(izp1)
            t_wp(iz)  = -(stat(izp1,2) - stat(iz,2))*dzu_i(izp1)
         endif
         dudz = (uxym(izp1) - uxym(iz))*dzu_i(izp1)
         dvdz = (vxym(izp1) - vxym(iz))*dzu_i(izp1)
c
c ------------- gather all the budget terms
c
         t_tran(iz)  = t_wq(iz) + t_wp(iz) + t_tau(iz)
         t_rprod(iz) = -(dudz*uwle(iz) + dvdz*vwle(iz))
         !Old t_sprod had subgrid stuff
         !t_sprod(iz) =  (dudz*uwsb(iz) + dvdz*vwsb(iz))
         !Now make t_sprod the spray tke term to reduce new variables:
         t_sprod(iz) = -stat(iz,5)
         t_buoy(iz)  =  batag*wtle(iz,1)
         t_diss(iz) = stat(iz,3)
c
      enddo
c
      return
      end
