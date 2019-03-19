      subroutine Tvar_budget(iscl)
      use pars
      use particles
      use fields
      use con_data
      use con_stats
      use fftwk
      implicit none

      real :: stat(1:nnz,4)
      real :: weit,weit1
      real :: tx_tmp(nnx,iys:iye), ty(nnx,iys:iye,izs-1:ize+1)
      real :: tx(nnx,iys:iye,izs-1:ize+1)
      real :: trans(izs:ize)
      real :: gradTp(3),Tfluc,dTpdz1,dTpdz,dTdz
      real :: Tflucp1,Tflucm1,qfluc,qflucp1,Tmean
      integer :: iz,i,j,izp1,izm1,iscl

c -------- stat(.,1) = Transport: -del.[U<T'2> + <u T'2> - alpha*grad(T'2)]
c          stat(.,2) = Dissipation: -2*alpha <grad(T').grad(T')> 
c          stat(.,3) = Particle: <T' Q'>

      !Need the y gradient of temp:
      do iz=izs-1,ize+1
      do j=iys,iye
      do i=1,nnx
         ty(i,j,iz)  = t(i,j,iscl,iz)
         tx_tmp(i,j)  = t(i,j,iscl,iz)
      enddo
      enddo
      call xderivp(tx_tmp(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      tx(1:nnx,iys:iye,iz) = tx_tmp(1:nnx,iys:iye)
      enddo
      call yd_mpi(ty(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)

      stat = 0.0
      Tv_part2 = 0.0
      do iz=izs,ize
         izp1 = iz + 1
         izm1 = iz - 1
         weit = dzw(iz)/(dzw(iz) + dzw(izp1))
         weit1 = 1.0 - weit

       do j=iys,iye
       do i=1,nnx

         if (iz==1)  then
         Tmean = 2.0*Tbot(iscl) - txym(iz,iscl)
         Tflucm1 = t(i,j,iscl,izm1)-Tmean
         Tflucp1 = t(i,j,iscl,izp1)-txym(izp1,iscl)
         elseif (iz==nnz) then
         Tmean = 2.0*Ttop(iscl) - txym(iz,iscl)
         Tflucp1 = t(i,j,iscl,izp1)-Tmean
         Tflucm1 = t(i,j,iscl,izm1)-txym(izm1,iscl)
         else
         Tflucp1 = t(i,j,iscl,izp1)-txym(izp1,iscl)
         Tflucm1 = t(i,j,iscl,izm1)-txym(izm1,iscl)
         end if
         Tfluc = t(i,j,iscl,iz)-txym(iz,iscl)

         !First dissipation: 
         !Note that gradients of total T and T' in x,y directions are equal since d<T>/dx = d<T>/dy = 0
         gradTp(1) = weit1*tx(i,j,iz) + weit*tx(i,j,izp1)
         gradTp(2) = weit1*ty(i,j,iz) + weit*ty(i,j,izp1)

         !Now get dT'/dz:
         gradTp(3) = (Tflucp1 - Tfluc)*dzu_i(izp1)
         
         stat(iz,2) = stat(iz,2)  - 2.0*vis_s(i,j,iscl,iz)*
     +                (gradTp(1)**2+gradTp(2)**2+gradTp(3)**2)

         !Next transport

         !Store the transport sum at u-locations since z-derivative at the end

         stat(iz,1) = stat(iz,1) + w(i,j,iz)*Tfluc**2

         dTpdz1 = (Tflucp1**2-Tfluc**2)*dzu_i(izp1)
         dTpdz = (Tfluc**2-Tflucm1**2)*dzu_i(iz)
         stat(iz,1) = stat(iz,1) - vis_s(i,j,iscl,iz)*0.5*(dTpdz1+dTpdz)

          !Particle source:
          if (iscl == 1) then
          if (iTcouple ==1) then
            qfluc = partTsrc(i,j,iz)-Tpsrc(iz)
            if (iz == nnz) then
            qflucp1 = 0.0
            else
            qflucp1 = partTsrc(i,j,izp1)-Tpsrc(izp1)
            end if
            stat(iz,3) = stat(iz,3) + 
     +         weit*(qflucp1*Tflucp1) + weit1*(qfluc*Tfluc)
          endif
          endif


          if (iscl == 2) then
          if (iHcouple ==1) then
            qfluc = partHsrc(i,j,iz)-Hpsrc(iz)
            if (iz == nnz) then
            qflucp1 = 0.0
            else
            qflucp1 = partHsrc(i,j,izp1)-Hpsrc(izp1)
            end if
            stat(iz,3) = stat(iz,3) + 
     +         weit*(qflucp1*Tflucp1) + weit1*(qfluc*Tfluc)
          endif
          endif

          if (iscl .eq. 2 .and. iHcouple .eq. 1 ) then
              qfluc = partTEsrc(i,j,iz)-TEpsrc(iz)
               if (iz == nnz) then
               qflucp1 = 0.0
               else
              qflucp1 = partTEsrc(i,j,izp1)-TEpsrc(izp1)
              end if
              Tv_part2(iz) = Tv_part2(iz) + 
     +          weit*(qflucp1*Tflucp1) + weit1*(qfluc*Tfluc)
  
         endif
     
       end do
       end do

         stat(iz,1) = stat(iz,1)*fnxy
         stat(iz,2) = stat(iz,2)*fnxy
         stat(iz,3) = stat(iz,3)*fnxy
         if (iscl .eq. 2 .and. iHcouple .eq. 1 ) then
         Tv_part2(iz) = Tv_part2(iz)*fnxy
         endif
      end do
         

      call mpi_sum_z(stat(1,1),i_root,myid,nnz*3,1)
      if (iscl .eq. 2 .and. iHcouple .eq. 1 ) then
      call mpi_sum_z(Tv_part2(1),i_root,myid,nnz,1)
      endif


c
c ------ we have all terms on all processors for all z, add them up
c
      do iz=1,nnz
         izp1 = iz + 1
         izm1 = iz - 1
         if(iz .eq. nnz) then
            Tv_tran(iz,iscl) = 0.0
         else
            Tv_tran(iz,iscl) = -(stat(izp1,1)-stat(iz,1))*dzu_i(izp1)
         endif
         dTdz = (txym(izp1,iscl)-txym(iz,iscl))*dzu_i(izp1)
c
c ------------- gather all the budget terms
c
         Tv_prod(iz,iscl) = -2.0*wtle(iz,iscl)*dTdz
         Tv_diss(iz,iscl) = stat(iz,2)
         Tv_part1(iz,iscl) = -stat(iz,3)
         !Tv_part2 is already gathered
      enddo
c
      return
      end
