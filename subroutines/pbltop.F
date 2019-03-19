      subroutine pbltop(itop)
c
c ---------- get estimate of pbl top
c
c            method = 0, min of wt flux
c                        (good for buoyancy cases)
c            method = 1, uw flux less than critical value
c                        (good for ekman cases)
c            method = 2, running t average exceeds criterion
c                        (good for neutral cases with capping
c                         inversions)
c            method = 3, maximum gradient in temperature field
c                        (good for finding local zi see jas paper)
c                        with minimum search height (sr. setup)
c
c ------------ if method uses average statistics then only root
c              process need find zi
c
      use pars
      use fields
      use con_data
      use con_stats
      real trun(maxnz)
      include 'mpif.h'
      real gradloc(2,nnx,nny), gradmax(2,nnx,nny)
      external get_zi
c
      if(method .le. 2 .and. l_root) then
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = 1.0
      if (method .le. 0 .or. method .gt. 2) then
         itop=1
         wttot=wtle(1,1)+wtsb(1,1)
         wtmin=wttot*sgn
         do iz=2,nnz
            wttot=(wtle(iz,1)+wtsb(iz,1))*sgn
            if (wttot.le.wtmin) then
               itop=iz
               wtmin=wttot
            endif
         enddo
         zi=z(itop)
      else if (method .eq. 1) then
         itop = 1
         crit = 0.05
         uwsf = utau*utau
         do iz=1,nnzm1
               uwtot = (uwle(iz) + uwsb(iz))**2 +
     $                 (vwle(iz) + vwsb(iz))**2
               uwtot = sqrt(uwtot)
               if(uwtot/uwsf .gt. crit) then
                  itop=iz
               endif
         enddo
         zi=z(itop)
      else if (method .eq. 2) then
         trun(1) = txym(1,1)
         do iz=2,nnz
             weight = z(iz-1)/z(iz)
             trun(iz) = trun(iz-1)*weight + (1.0-weight)*txym(iz,1)
         enddo
         itop = 1
         tcrit = 0.25
         if(iocean .eq. 1) tcrit = 0.1
         do iz=2,nnz
                if(txym(iz,1) .gt. (trun(iz) + tcrit)) then
                  itop = iz
                  go to 320
                endif
         enddo
  320    continue
         zi=z(itop)
      endif
      do iy=1,nny
      do ix=1,nnx
         gradmax(2,ix,iy) = zi
      enddo
      enddo
c
c ----------- use gradient method, every process computes
c
      elseif(method .eq. 3) then
c
c ---------------- get local zi from gradient in temperaure field
c
c     dz_i = dzu_i(izs+1)
c     do iy=1,nny
c     do ix=1,nnx
c        gradloc(1,ix,iy) = (t(ix,iy,1,izs+1) - t(ix,iy,1,izs))*dz_i
c        gradloc(2,ix,iy) = z(izs)
c     enddo
c     enddo
c
c ------- similar to zeroing the stat array in sr. mean_stat
c
      do iy=1,nny
      do ix=1,nnx
         gradloc(1,ix,iy) = 0.0
         gradloc(2,ix,iy) = z(iz_min)
      enddo
      enddo
c
c ------------- now all z in this process
c
      if(iz_min .le. ize) then
      do iz=max(izs,iz_min),ize
         izp1 = iz + 1
         do iy=iys,iye
         do ix=1,nnx
            grad = (t(ix,iy,1,izp1) - t(ix,iy,1,iz))*dzu_i(izp1)
            if(grad .gt. gradloc(1,ix,iy)) then
               gradloc(1,ix,iy) = grad
               gradloc(2,ix,iy) = z(iz)
            endif
         enddo
         enddo
      enddo
      endif
c
c     call mpi_reduce(gradloc,gradmax,2*nnx*nny,mpi_real8,ziloc,
c    +                i_root,mpi_comm_world,ierror)
c
c ----------- alternate version using already defined function in mpi
c             passes 2 real8 variables
c
      call mpi_reduce(gradloc,gradmax,nnx*nny,mpi_2double_precision,
     +                mpi_maxloc,i_root,mpi_comm_world,ierror)
c
c ------------ get average on root process
c
      if(l_root) then
         zi_avg = 0.0
         do iy=1,nny
         do ix=1,nnx
            zi_avg = zi_avg + gradmax(2,ix,iy)
         enddo
         enddo
         zi = zi_avg*fnxy
c        itop = nint(zi/dz)
      endif
c
      endif
c
c -------- send average zi everywhere
c
      call mpi_bcast(zi,1,mpi_real8,
     +              i_root,mpi_comm_world,ierr)
c
      if(iocean .ne. 1) then
         do iz=1,nnz
            if(zi .ge. z(iz) .and.
     +         zi .lt. z(iz+1)) itop = iz
         enddo
      else
         do iz=1,nnz
            if(zi .le. z(iz) .and.
     +         zi .gt. z(iz+1)) itop = iz
         enddo
      endif
c
c     if(l_root) write(6,7001) myid,zi,itop
 7001 format(' 7001 in pbltop myid = ',i4,' zi = ',e15.6,
     +       ' itop = ',i3)
c
      return
      end
