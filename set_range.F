      subroutine set_range
c
c ---- build special x,y,z-ranges. dimensioned for 0:numprocs-1
c      indexed with myid
c
c      [ixs:ixe] x-range designated for a specific cpu
c      [iys:iye] y-range designated for a specific cpu
c      [izs:ize] z-range designated for a specific cpu
c      for example, myid==0 will have ixs==1 and ixs==22 for
c      a maxnx==128 and ncpu_slab=6 
c
c
c
c      [ix_s:ix_e] x-range for computing y-derivatives nx-pts/ncpu_s
c                  in xtoy and ytox tranposes
c
c      [jx_s:jx_e] x-range for computing 2d fft (nx+2)-pts/ncpu_s
c                  must be even in each x-interval for complex fft in y
c
c      [kx_s:kx_e] x-range for pressure solver transpose (nx+2)-pts/ncpu_z
c                  nx+2 fourier coefficients for xtoz and ztox transposes
c
c      [mx_s:mx_e] x-range split across z cpus as nx-pts/ncpu_z
c                  for use in surface layer routines 
c
c      [is_s:is_e] starting and ending processor id's for a
c                  particular z-level
c
c      [iy_s:iy_e] y-range for computing y-derivatives ny-pts/ncpu_s
c                  in xtoy and ytox tranposes
c
c      [jy_s:jy_e] y-range for use in xtoz and ztox transposes 
c                  in pressure solution
c
c      [iz_s:iz_e] z-range for a particular vertical slab
c
c
      use pars
c
      write(nprt,7002) nnx,nny,nnz
 7002 format(' 7002 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)
c
      ii = -1
      do nn=0,ncpu_z-1
         call range(1,nnx+2,ncpu_z,nn,lx_s,lx_e)
         call range(1,nnx,ncpu_z,nn,nx_s,nx_e)
         call range(1,nny,ncpu_z,nn,ly_s,ly_e)
         call range(1,nnz,ncpu_z,nn,mz_s,mz_e)
         do mm=0,ncpu_s-1
            call range(1,nny,ncpu_s,mm,ny_s,ny_e)
            call range(1,nnx,ncpu_s,mm,nxy_s,nxy_e)
            call range(1,ncx,ncpu_s,mm,l2x_s,l2x_e)
            ii       = ii + 1
c
            ix_s(ii) = nxy_s
            ix_e(ii) = nxy_e
            jx_s(ii) = (l2x_s - 1)*2 + 1
            jx_e(ii) = l2x_e*2
            kx_s(ii) = lx_s
            kx_e(ii) = lx_e
            mx_s(ii) = nx_s
            mx_e(ii) = nx_e
c
            iy_s(ii) = ny_s
            iy_e(ii) = ny_e
            jy_s(ii) = ly_s
            jy_e(ii) = ly_e
c
            iz_s(ii) = mz_s
            iz_e(ii) = mz_e
c
            is_s(ii) = (ii/ncpu_s)*ncpu_s
            is_e(ii) = is_s(ii) + ncpu_s - 1
         enddo
      enddo
c
      iys = iy_s(myid)
      iye = iy_e(myid)
      jys = jy_s(myid)
      jye = jy_e(myid)
      ixs = ix_s(myid)
      ixe = ix_e(myid)
      jxs = jx_s(myid)
      jxe = jx_e(myid)
      kxs = kx_s(myid)
      kxe = kx_e(myid)
      mxs = mx_s(myid)
      mxe = mx_e(myid)
      izs = iz_s(myid)
      ize = iz_e(myid)
c
c ----------- get starting and  ending processor id's on each
c             vertical slab
c
      iss = is_s(myid)
      ise = is_e(myid)
c
c ------------ debug ranges 
c
      if(l_debug) then
         write(nprt,1200) myid, (nn, ix_s(nn), ix_e(nn), jx_s(nn),
     +                     jx_e(nn), kx_s(nn), kx_e(nn),
     +                     nn = 0,numprocs-1)
 1200    format(' myid =  ',i4,/,
     +       ' nn',5x,' ixs ',5x,' ixe ',5x,' jxs ',5x,' jxe '
     +       ,5x,' kxs ',5x,' kxe',/,(7i6))
c
         write(nprt,1213) myid, (nn, iy_s(nn), iy_e(nn),
     +                   jy_s(nn), jy_e(nn),
     +                   iz_s(nn), iz_e(nn), is_s(nn), is_e(nn),
     +                   nn=0,numprocs-1)

 1213    format(' myid = ',i3,/,
     +       ' nn ',3x,' iys ',5x,' iye ',5x,
     +       ' jys ',5x,' jye ',5x,
     +       ' izs ',5x,' ize',5x,' iss ',5x,' ise ',/,
     +       (9i6))
      endif
c
      return
      end
