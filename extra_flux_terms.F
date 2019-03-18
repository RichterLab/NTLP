      subroutine extra_flux_terms
      use pars
      use particles
      use fields
      use con_data
      use con_stats
      use fftwk
      implicit none

      real :: stat(1:nnz,3)
      real :: weit,weit1
      real :: Tflucp1,Tflucm1,Tfluc,Tmean
      real :: qflucp1,qflucm1,qfluc,qmean,dqpdz
      real :: Sq,Sqp1
      real :: gradTp(3),dTpdz1,dTpdz,dTdz
      integer :: iz,i,j,izp1,izm1,iscl

!     Compute the "extra" enthalpy budget terms, located at w locations:
c -------- stat(.,1) = < q' T' w' > 
c          stat(.,2) = < T' Sq > where Sq is q source
c          stat(.,3) = Dv*< T' dq'/dz > 

       

      stat = 0.0
      do iz=izs,ize
         izp1 = iz + 1
         izm1 = iz - 1
         weit = dzw(iz)/(dzw(iz) + dzw(izp1))
         weit1 = 1.0 - weit

       do j=iys,iye
       do i=1,nnx

         if (iz==1)  then
           Tmean = 2.0*Tbot(1) - txym(iz,1)
           Tflucm1 = t(i,j,1,izm1)-Tmean
           Tflucp1 = t(i,j,1,izp1)-txym(izp1,1)
         elseif (iz==nnz) then
           Tmean = 2.0*Ttop(1) - txym(iz,1)
           Tflucp1 = t(i,j,1,izp1)-Tmean
           Tflucm1 = t(i,j,1,izm1)-txym(izm1,1)
         else
           Tflucp1 = t(i,j,1,izp1)-txym(izp1,1)
           Tflucm1 = t(i,j,1,izm1)-txym(izm1,1)
         end if
         Tfluc = t(i,j,1,iz)-txym(iz,1)

         if (iz==1)  then
           qmean = 2.0*Tbot(2) - txym(iz,2)
           qflucm1 = t(i,j,2,izm1)-qmean
           qflucp1 = t(i,j,2,izp1)-txym(izp1,2)
           Sqp1 = -partHsrc(i,j,izp1)
         elseif (iz==nnz) then
           qmean = 2.0*Ttop(2) - txym(iz,2)
           qflucp1 = t(i,j,2,izp1)-qmean
           qflucm1 = t(i,j,2,izm1)-txym(izm1,2)
           Sqp1 = 0.0
         else
           qflucp1 = t(i,j,2,izp1)-txym(izp1,2)
           qflucm1 = t(i,j,2,izm1)-txym(izm1,2)
         end if
         qfluc = t(i,j,2,iz)-txym(iz,2)
         Sq = -partHsrc(i,j,iz)

         !Get dq'/dz at w-locations:
         dqpdz = (qflucp1 - qfluc)*dzu_i(izp1)

         !Then can get source #3:
         stat(iz,3) = stat(iz,3) + vis_s(i,j,2,iz)*
     +   (weit1*Tfluc + weit*Tflucp1)*dqpdz


         !Now source #1, which must be evaluated at w-location:

         stat(iz,1) = stat(iz,1) + w(i,j,iz)*
     +   (weit1*Tfluc + weit*Tflucp1)*
     +   (weit1*qfluc + weit*qflucp1)


         !Now source #2:
         stat(iz,2) = stat(iz,2) + 
     +   (weit1*Tfluc + weit*Tflucp1)*
     +   (weit1*Sq + weit*Sqp1)

     
       end do
       end do

         stat(iz,1) = stat(iz,1)*fnxy
         stat(iz,2) = stat(iz,2)*fnxy
         stat(iz,3) = stat(iz,3)*fnxy
      end do
         

      call mpi_sum_z(stat(1,1),i_root,myid,nnz*3,1)


c
c ------ we have all terms on all processors for all z, add them up
c
      do iz=1,nnz
c
c ------------- gather all the budget terms
c
         trip(iz) = stat(iz,1)
         TpSq(iz) = stat(iz,2)
         Tpdqp(iz) = stat(iz,3)
      enddo
c
      return
      end
