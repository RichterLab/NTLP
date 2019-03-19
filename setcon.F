      subroutine setcon
c
      use pars
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      external get_zi
c
c ----------------- get machine type, can erad datadir also
c
      open(99,file='mach.file')
      read(99,9000) imach
 9000 format(i1)
c
      close(99)
c
c ----------- initialize fft
c
      call rffti(nnx,trigx(1,1))
      call rffti(nny,trigx(1,2))
      call cffti(nny,trigc(1))
c
c ----------- start step for history files
c
      it_his_nxt = it_his
      it_viz_nxt = it_viz
      it_pro_nxt = it_probe
c
c ---------------- set min value of e
c
      if(iocean .eq. 1) then
c
         smal_e = 0.0
         smal_e = 1.0e-12
c        smal_e = 6.0e-03
      else
         smal_e = 1.0e-09
c        smal_e = 0.0
      endif
c
c ---------------------- set constants in eddy viscosity model
c
      ck       = 0.1
      ceps     = 0.93
      csmag    = sqrt(ck*sqrt(ck/ceps))
      stab_c   = 0.76
c
c ----------------- set stability constant
c
      stabmin = 1.0e-12
c
c ---------------- minimum dsl length constant
c 
      almin_c = 0.0001
c
c -------------- initialize grid restart flag
c
      igrdr = 1
c
c -------------- create mpi operation to find max and location
c                using local gradient method
c
      call mpi_op_create(get_zi,.true.,ziloc,ierror)
c
c ------------------- define coefficients for 3-order runge-kutta
c                     time stepping scheme, borrowed from Spalart,
c                     Moser and Rogers, J. Comp. Physics 3/21/90
c                     Note this is a simplier version since all terms
c                     are lumped in the non-linear terms.
c                     cfl number is for an entire runge-kutta step
c                     in this case three stages. cfl = max(u)*dt/dx
c
      zetas(1) = 0.0
      zetas(2) = -17.0/60.0
      zetas(3) = -5.0/12.0
      gama(1)  = 8.0/15.0
      gama(2)  = 5.0/12.0
      gama(3)  = 3.0/4.0
      etas(1)  = -1.0
      etas(2)  = -1.0 + 8.0/15.0
      etas(3)  = -1.0 + 2.0/3.0
c
c ----------- a full step, at the new time
c
      etas(4)  =  0.0
c
      cfl = 0.40
      !cfl = 0.2  !For phi0.25_rho1000_r100_Mun_Tcpl
c     cfl = 0.50
c
      return
c
      end
