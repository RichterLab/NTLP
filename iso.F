      subroutine iso(it)
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      real sfk(1:nnz)



c
c ---- get isotropy factor and scale it to match at the matching
c      height. uses boundary conditions from lower and upper. 
c
      do iz=1,nnz
c        dfac(iz) = 1.0
         dfac(iz) = 0.0
         sfk(iz)  = 0.0
      enddo
      do iz=izs,ize
         dfac(iz) = 1.0
      enddo
c
c ------ set nmatch equal to fraction of initial zi in sr. random
c
c     nmatch = izi/2
c     nmatch = 16
      nmatch = 48
      do i=0,numprocs-1,ncpu_s
         if(nmatch .ge. iz_s(i) .and.
     +      nmatch .le. iz_e(i)) myid_newvis = i
      enddo
c
      do iz=izs,min(ize,nmatch)
         izp1 = iz + 1
         izm1 = iz - 1
         weit = dzw(iz)/(dzw(iz) + dzw(izp1))
         weit1 = 1.0 - weit
c
c ---- get fluctuating strains
c
         do j=iys,iye
         do i=1,nnx
            s11 = weit1*ux(i,j,iz)**2 + weit*ux(i,j,izp1)**2
            s22 = weit1*vy(i,j,iz)**2 + weit*vy(i,j,izp1)**2
            wz  = (w(i,j,iz)-w(i,j,izm1))*dzw_i(iz)
            wzp = (w(i,j,izp1)-w(i,j,iz))*dzw_i(izp1)
            s33 = weit*wzp**2 + weit1*wz**2
            s12 = weit1*(uy(i,j,iz) + vx(i,j,iz))**2 +
     +            weit*(uy(i,j,izp1) + vx(i,j,izp1))**2
            s13 = (((u(i,j,izp1) - u(i,j,iz) +
     +            u_mn(iz) - u_mn(izp1))*dzu_i(izp1) +
     +            wx(i,j,iz)))**2
            s23 = (((v(i,j,izp1) - v(i,j,iz) +
     +          v_mn(iz) - v_mn(izp1))*dzu_i(izp1) +
     +          wy(i,j,iz)))**2
            sfk(iz) = sfk(iz) + 2.0*(s11 + s22 + s33) +
     +                       s12 + s13 + s23
         enddo
         enddo
         sfk(iz) = sfk(iz)*fnxy
      enddo
      call mpi_sum_z(sfk,i_root,myid,nnz,1)
c
      do iz=izs,min(ize,nmatch)
         izp1 = iz + 1
         izm1 = iz - 1
c
         sfk(iz) = sqrt(sfk(iz))
         smk = sqrt((u_mn(izp1)-u_mn(iz))**2 +
     +              (v_mn(izp1)-v_mn(iz))**2)*abs(dzu_i(izp1))
         if(sfk(iz) .le. 0. .and. smk .le. 0.) then
           dfac(iz) = 1.0
         else
           dfac(iz) = sfk(iz)/(sfk(iz) + smk)
         endif
c     if(l_root) write(6,6001) iz,sfk(iz),smk,dfac(iz)
 6001 format(' iz = ',i3,' sfk = ',e15.6,
     +       ' smk = ',e15.6,' dfac = ',e15.6)
      enddo

c
c
c ---- rescale ratio to give unity at match height
c      and if nested grid match value at upper boundary
c      of coarser grid
c
      if(myid .eq. myid_newvis) then
         dfacm = dfac(nmatch)
      endif
c

      call mpi_bcast(dfacm,1,mpi_real8,
     +              myid_newvis,mpi_comm_world,ierr)
c
      do iz=izs,min(ize,nmatch)
         dfac(iz) = dfac(iz)/dfacm
         dfac(iz) = amax1(dfac(iz), 0.1)
         dfac(iz) = amin1(dfac(iz), 1.0)
      enddo
c
c --------- gather dfac on all processes for printing and use in tke_vis
c           use reduce and divide by number of slab cpus
c

      call mpi_sum_z(dfac,i_root,myid,nnz,1)
      fncpu_s = 1.0/float(ncpu_s)
      do iz=1,nnz
         dfac(iz) = dfac(iz)*fncpu_s
      enddo
c
c     call mpi_gath_root(dfac(izs),dfac(1),iz_s,iz_e,izs,ize,nnz,
c    +                   myid,numprocs,ncpu_s)
c
c     if(l_root) write(6,6000) nmatch,ivis,(iz,dfac(iz),iz=1,nnz)
 6000 format(' in sr. iso, nmatch = ',i3,/,
     +       ' ivis = ',i3,'iz',5x,'dfac',/,(i3,1x,e15.6))
c     write(nprt,3001) (iz,dfac(iz),iz=1,nnz)
 3001 format(' iz ',5x,' dfac ',/,(i5,e15.6))
      return
      end
