      subroutine vgrid(z1,zi,zl,nnz,z,l_root,ldebug)
c
      real z(0:nnz+1)
      logical l_root, l_debug
c
c ----------------- build grid up to zi first
c
      z_frst = z1
      z_cntr = zi*0.5
      n_pbl  = nnz/2
c     n_pbl  = (5*nnz)/8
      z_fac1 = z_cntr/z_frst
      z_fac2 = 1.0/float(n_pbl/2)
      z_fac  = 1.1
      knt = 0
      tol = 0.00001
   10 continue
        knt = knt + 1
        z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
        test   = abs(1.0 - z_facn/z_fac)
        if(knt .gt. 50) then
            if(l_root) write(6,9000) z_fac, z_facn, knt
 9000       format(' Cannot find stretching factor',/,
     +             ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 10
      if(l_root) write(6,9100) z_fac, z_cntr, z1, knt
 9100 format(' Stretching factor = ',e15.6,/,
     +       ' Match point       = ',e15.6,/,
     +       ' First z           = ',e15.6,/,
     +       ' Number of iters   = ',i4)
      z(1) = z_frst
      do iz=2,n_pbl/2-1
         z(iz) = z_frst*(z_fac**(float(iz)) - 1.0)/(z_fac - 1.0)
      enddo
      z(n_pbl/2) = z_cntr
      do iz=1,n_pbl/2 - 1
         z(n_pbl-iz) = zi - z(iz)
      enddo
      z(n_pbl) = zi
      z(0)   = 0.0
c
      if(l_root) write(6,5300) n_pbl
 5300 format(' n_pbl = ',i4)
c
c -------------- build grid from zi on up
c
      z_frst = z1
      z_cntr = zl - zi
      n_top  = nnz - n_pbl
      z_fac1 = z_cntr/z_frst
      z_fac2 = 1.0/float(n_top)
      z_fac  = 1.1
      knt = 0
      tol = 0.00001
   20 continue
        knt = knt + 1
        z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
        test   = abs(1.0 - z_facn/z_fac)
        if(knt .gt. 50) then
            if(l_root) write(6,8000) z_fac, z_facn, knt
 8000       format(' Cannot find stretching factor',/,
     +             ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 20
      if(l_root) write(6,8100) z_fac, z_cntr, z1, knt
 8100 format(' Stretching factor = ',e15.6,/,
     +       ' Match point       = ',e15.6,/,
     +       ' First z           = ',e15.6,/,
     +       ' Number of iters   = ',i4)
c
      z(n_pbl+1) = zi + z_frst
      do iz=n_pbl+2,nnz-1
         z(iz) = zi + z_frst*
     +           (z_fac**(float(iz-n_pbl)) - 1.0)/(z_fac - 1.0)
      enddo
      z(nnz) = zl
      z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))
c     if(l_root) write(6,5600) (iz,z(iz),iz=0,nnz+1)
 5600 format(' 5600 in vgrid ',/,
     +       ' iz ',5x,' zw ',/,(i3,e15.6))
c
c     write(1,2000)
c2000 format('#k ',/,
c    +       '#lw 0.5 ',/,
c    +       '#m 1',/,
c    +       '#x 0 100 50',/,
c    +       '#y -50 2100 500')
c     x1 = 30.0
c     x2 = 80.0
c     do iz=0,nnz+1
c        write(1,1000) x1,z(iz)
c1000    format('#k ',/,
c    +          (2e15.6))
c        write(1,1100) x2,z(iz)
c1100    format(2e15.6)
c     enddo
c
      return
      end
