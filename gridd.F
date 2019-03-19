
      subroutine gridd
c
c ----------- allocate space and pass arrays using modules
c
      use pars
      use fields
      use particles 
      use fftwk
      use con_data
      use con_stats
c
c ------------- establish association between pointers
c               and data structures
c
      call fill_cc
      call fill_cs

      if (myid==0) write(6,5001) isize
 5001 format(' size of stats array = ',i8)

c
c
c ---------------- debug for arrays
c
      big = -99.0e+300
c
c ---------------- setup grid
c
      nnx = maxnx
      nny = maxny
      nnz = maxnz
c     izs = 1
c     ize = nnz
c
c
c ----------- make sure problem and cpu's match
c
      maxp   = numprocs-1
      ncpu_z = numprocs/ncpu_s
      if(mod(numprocs,ncpu_s) .ne. 0 .or.
     +   ncpu_z .gt. nnz) then
         go to 999
      endif
      if(l_root) write(6, 1100) ncpu_s, ncpu_z, numprocs,
     +                          maxp
      write(nprt,1100) ncpu_s, ncpu_z, numprocs, maxp
 1100 format(' Number of x-y slab cpus = ',i5,/,
     +       ' Number of z-level cpus  = ',i5,/,
     +       ' Total number of cpus    = ',i5,/,
     +       ' Max-p for index arrays  = ',i5)
c
c ---------------- allocate arrays for (i,j,k)-indexing on
c                  each processor (see set_range)
c
      allocate(ix_s(0:maxp), ix_e(0:maxp),
     +         jx_s(0:maxp), jx_e(0:maxp),
     +         kx_s(0:maxp), kx_e(0:maxp),
     +         mx_s(0:maxp), mx_e(0:maxp),
     +         iy_s(0:maxp), iy_e(0:maxp),
     +         jy_s(0:maxp), jy_e(0:maxp),
     +         is_s(0:maxp), is_e(0:maxp),
     +         iz_s(0:maxp), iz_e(0:maxp))
c
c ---------------- setup array sizes and variable dimensions
c
      nxy   = nnx*nny
      ncx   = nnx/2 + 1
      ncy   = nny/2 + 1
      nnxp1 = nnx + 1
      nnyp1 = nny + 1
      nnxp2 = nnx + 2
      nnyp2 = nny + 2
      nnzp1 = nnz + 1
      nnzm1 = nnz - 1
      ivis = ivis0
      fnxy  = 1.0/float(nnx*nny)
c
      write(nprt,7001) nnx,nny,nnz
 7001 format(' 7001 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)
c
      call set_range
c
      write(nprt,7002) nnx,nny,nnz
 7002 format(' 7002 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)
c
      num_y = iye + 1 - iys


c
c ------------- allocate solution arrays
c               account for nnxp2 for fields but not in rhs
c               and possible monotone for scalars
c
      allocate(u(nnxp2,iys:iye,izs-1:ize+1), 
     +         v(nnxp2,iys:iye,izs-1:ize+1), 
     +         w(nnxp2,iys:iye,izs-1:ize+1), 
     +         t(nnxp2,iys:iye,nscl,izs-2:ize+2), 
     +         e(nnxp2,iys:iye,izs-1:ize+1), 
     +         r1(nnx,iys:iye,izs-1:ize+1),
     +         r2(nnx,iys:iye,izs-1:ize+1),
     +         r3(nnx,iys:iye,izs-1:ize+1),
     +         r4(nnx,iys:iye,nscl,izs-1:ize+1),
     +         r5(nnx,iys:iye,izs-1:ize+1))

c   ------Allocate arrays for SFS stress and its derivatives
c        Indrajith 07/17
      allocate(sigm_s(nnxp2,iys:iye,izs-1:ize+1),
     +         sigm_sdx(nnxp2,iys:iye,izs-1:ize+1),
     +         sigm_sdy(nnxp2,iys:iye,izs-1:ize+1),
     +         sigm_sdz(nnxp2,iys:iye,izs-1:ize+1),
     +         vis_ss(nnxp2,iys:iye,izs-1:ize+1))

c
c ------------- allocate extended arrays for interpolation of
c               particle/spray location
c
      if (ispray==1) then
      allocate(uext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3), 
     +         vext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3), 
     +         wext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
c     +         Text(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3))
     +         Text(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),!Brian 7/14/14 
     +         T2ext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3))  !Brian 7/14/14 

c    ---allocate extented arrys for interpolation of SFS-stress 
c       and its derivatives, Indrajith 07/17      
      allocate (sigm_sext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         sigm_sdxext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         sigm_sdyext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         sigm_sdzext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         vis_sext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3))

c     Transposed velocities to do the uf interpolation:
      allocate(u_t(0:nnz+1,iys:iye,mxs:mxe),
     +         v_t(0:nnz+1,iys:iye,mxs:mxe),
     +         w_t(0:nnz+1,iys:iye,mxs:mxe), 
     +         T_t(0:nnz+1,iys:iye,mxs:mxe),!Brian 7/14/14 
     +         T2_t(0:nnz+1,iys:iye,mxs:mxe)) !Brian 7/14/14 
      end if
c      Keep track of particle counts at each node (its surrounding volume)
c      NOTE: allocate even if ispray == 0, since it's in xy_stats
      allocate(partcount(nnx,iys:iye,izs-1:ize+1))
      allocate(partcount_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(partflux(nnx,iys:iye,izs-1:ize+1))
      allocate(partflux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(partcountU(nnx,iys:iye,izs-1:ize+1))
      allocate(partcountU_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pcountrange(nnx,iys:iye,izs-1:ize+1))
      allocate(pcountrange_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(Spartflux(nnx,iys:iye,izs-1:ize+1))
      allocate(Spartflux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(tpartflux(nnx,iys:iye,izs-1:ize+1))
      allocate(tpartflux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(tupartflux(nnx,iys:iye,izs-1:ize+1))
      allocate(tupartflux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(tdpartflux(nnx,iys:iye,izs-1:ize+1))
      allocate(tdpartflux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(xytpartflux(nnx,iys:iye,izs-1:ize+1))
      allocate(xytpartflux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(xyutpartflux(nnx,iys:iye,izs-1:ize+1))
      allocate(xyutpartflux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(xydtpartflux(nnx,iys:iye,izs-1:ize+1))
      allocate(xydtpartflux_t(0:nnz+1,iys:iye,mxs:mxe))

      allocate(pp1flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp1flux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pp2flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp2flux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pp3flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp3flux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pp4flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp4flux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pp5flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp5flux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pp6flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp6flux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pp7flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp7flux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pp8flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp8flux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pp9flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp9flux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pp10flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp10flux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pp11flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp11flux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pp12flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp12flux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pp13flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp13flux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pp14flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp14flux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pp15flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp15flux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(pp16flux(nnx,iys:iye,izs-1:ize+1))
      allocate(pp16flux_t(0:nnz+1,iys:iye,mxs:mxe))


      allocate(upwp_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(upwp(nnx,iys:iye,izs-1:ize+1))
      allocate(vpsum(nnx,iys:iye,izs-1:ize+1,1:3))
      allocate(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,1:3))
      allocate(vpsqrsum(nnx,iys:iye,izs-1:ize+1,1:3))
      allocate(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,1:3))
      allocate(Tpsum(nnx,iys:iye,izs-1:ize+1))
      allocate(Tpsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(Tpsqrsum(nnx,iys:iye,izs-1:ize+1))
      allocate(Tpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(Tfsum(nnx,iys:iye,izs-1:ize+1))
      allocate(Tfsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(qfsum(nnx,iys:iye,izs-1:ize+1))
      allocate(qfsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(wpTpsum(nnx,iys:iye,izs-1:ize+1))
      allocate(wpTpsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(partsrc(nnx,iys:iye,izs-1:ize+1,1:3))
      allocate(partsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1,1:3))
      allocate(partTsrc(nnx,iys:iye,izs-1:ize+1))
      allocate(partTsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1))
      allocate(radsum(nnx,iys:iye,izs-1:ize+1))
      allocate(radsum_t(0:nnz+1,iys:iye,mxs:mxe)) 
      allocate(rad2sum(nnx,iys:iye,izs-1:ize+1)) 
      allocate(rad2sum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(mpsum(nnx,iys:iye,izs-1:ize+1)) 
      allocate(mpsum_t(0:nnz+1,iys:iye,mxs:mxe)) 
      allocate(mwsum(nnx,iys:iye,izs-1:ize+1)) 
      allocate(mwsum_t(0:nnz+1,iys:iye,mxs:mxe)) 
      allocate(partHsrc(nnx,iys:iye,izs-1:ize+1))
      allocate(partHsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1))
      allocate(partTEsrc(nnx,iys:iye,izs-1:ize+1))
      allocate(partTEsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1))
      allocate(qstarsum(nnx,iys:iye,izs-1:ize+1)) 
      allocate(qstarsum_t(0:nnz+1,iys:iye,mxs:mxe)) 


c ------------- allocate space for boundary condition arrays
c               on top and bottom of domain
c
      allocate(ubc(nnx,iys:iye,2),
     +         vbc(nnx,iys:iye,2),
     +         wbc(nnx,iys:iye,2),
     +         tbc(nnx,iys:iye,nscl,2),
     +         ebc(nnx,iys:iye,2),
     +         pbc(nnx,iys:iye,2),
     +         pbc2(nnx,iys:iye,2))
c
c ------------ allocate space for wind and surface arrays
c
      allocate(wind(nnx,iys:iye), 
     +         tau13m(nnx,iys:iye), 
     +         tau23m(nnx,iys:iye), 
     +         taut3m(nnx,iys:iye,nscl), 
     +         t_grnd(nnx,iys:iye,nscl))
c
c ------------------- allocate space for derivative arrays
c
      allocate(ux(nnx,iys:iye,izs-1:ize+1),
     +         uy(nnx,iys:iye,izs-1:ize+1),
     +         vx(nnx,iys:iye,izs-1:ize+1),
     +         vy(nnx,iys:iye,izs-1:ize+1),
     +         wx(nnx,iys:iye,izs-1:ize+1),
     +         wy(nnx,iys:iye,izs-1:ize+1))

c
c ------------- allocate space for subgrid velocity Indrajith 06/16/2017

c      allocate(u_sub(nnx,iys:iye,izs-1:ize+1),
c     +          v_sub(nnx,iys:iye,izs-1:ize+1),
c     +          w_sub(nnx,iys:iye,izs-1:ize+1))
c ------------- allocate space for pressure, pressure bcs
c
      allocate(p(nnxp2,iys:iye,izs-1:ize+1),
     +         ptop(nnxp2,iys:iye,2))
c
c ------------- allocate space for viscosity and diffusivity
c
      allocate(vis_m(nnx,iys:iye,izs-1:ize+1),
     +         vis_s(nnx,iys:iye,nscl,izs-1:ize+1))
c
c ------------- allocate space for fft trig factors
c
      nq_trig = max(nnx,nny)
      allocate(trigx(2*nq_trig+15,2),
     +         trigc(4*nq_trig+15))
c
      return
  999 continue
c
      if(l_root) write(6,1000) numprocs, ncpu_s, mmz
      write(nprt,1000) numprocs, ncpu_s, nnz
 1000 format(' Gridd Trouble number of processors and grid',
     +          ' partitioning do not match!',/,
     +          ' Total num of cpus   = ',i5,
     +          ' Num cpu on x-y slab = ',i5,/,
     +          ' Num of z-levels     = ',i5)
      call mpi_finalize(ierr)
      end
