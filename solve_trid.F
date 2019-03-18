      subroutine solve_trid(pt, ptop)
c 
c --------- tridiagonal solver. odd order for ptop, ptop2
c           because of 2d-fft
c
      use pars
      use con_data
      use con_stats
c
      real ptop(nny,jxs:jxe,1:2)
      real pt(0:nnz+1,jxs:jxe,jys:jye)
      real aa(nnz,jxs:jxe),bb(nnz,jxs:jxe),
     +     dd(nnz,jxs:jxe),rh(nnz,jxs:jxe)
      real fac_u(nnz), fac_l(nnz), fac_a(nnz)
c
      do iz=1,nnz
         fac_u(iz) = 1.0/(dzw(iz)*dzu(iz+1))
         fac_l(iz) = 1.0/(dzw(iz)*dzu(iz))
         fac_a(iz) = fac_l(iz) + fac_u(iz)
      enddo
c
      do kp=jys,jye    
         do lp=jxs,jxe
         do iz=2,nnz-1
            bb(iz,lp)  = fac_l(iz)
            aa(iz,lp)  = fac_u(iz)
            dd(iz,lp)  = -xks(lp,kp) - fac_a(iz)
            rh(iz,lp)  = pt(iz,lp,kp)
         enddo
         enddo
c
c --------------- lower boundary, fill exterior pressure (not used)
c
         do lp=jxs,jxe
            bb(1,lp)  = 1.0
            aa(1,lp)  = fac_u(1)
            dd(1,lp)  = -xks(lp,kp) - fac_u(1)
            rh(1,lp)  = pt(1,lp,kp)
            pt(0,lp,kp) = 0.0
         enddo
c
c --------------- upper boundary, fill exterior pressure (not used)
c
         if(ibcu .eq. 1) then
            do lp=jxs,jxe
              bb(nnz,lp) = 0.0
              aa(nnz,lp) = 0.0
              dd(nnz,lp) = 1.0
              rh(nnz,lp) = ptop(kp,lp,1)*wavexy(lp,kp) + ptop(kp,lp,2)
              pt(nnz+1,lp,kp) = 0.0
            enddo
         else
            do lp=jxs,jxe
               bb(nnz,lp) = fac_l(nnz)
               aa(nnz,lp) = 1.0
               dd(nnz,lp) = -xks(lp,kp) - fac_l(nnz)
               rh(nnz,lp) = pt(nnz,lp,kp)
               pt(nnz+1,lp,kp) = 0.0
            enddo
         endif
c
c ---------------- special situation for zeroth mode
c                  makes mean pressure = 0
c
         if(kp .eq. 1 .and. jxs .eq. 1) then
           do iz=1,nnz
              dd(iz,1) = 1.0
              rh(iz,1) = 0.0
              aa(iz,1) = 0.0
              bb(iz,1) = 0.0
              dd(iz,2) = 1.0
              rh(iz,2) = 0.0
              aa(iz,2) = 0.0
              bb(iz,2) = 0.0
           enddo
         endif
c
c --------------- solve system
c
         call tridv(bb,dd,aa,rh,nnz,jxs,jxe)
         do lp=jxs,jxe
         do iz=1,nnz
            pt(iz,lp,kp) = rh(iz,lp)
         enddo
         enddo
      enddo
c
      return
      end
