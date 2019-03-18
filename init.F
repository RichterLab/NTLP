      subroutine init
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
c
      real :: RHT,RHB
      character(len=50) :: wname

      debug_var = 1
      pi   = 4.0*atan(1.0)
      pi2  = 2.0*pi
      bfac = 1.0
      if(ibuoy.eq.0) bfac = 0.
c
c -------------------- case specific data
c
      if(iocean .eq. 1) then
         t00     = 283.
         t00b    = 5000.0
         cp      = 4.20e03
         gcp     = grav/cp
         batag   = bfac*grav/t00b
c        fcor    = 0.0
         fcor    = 1.39e-04
         fcor_h  = 0.0
         ugcont  = 0.0
         vgcont  = 0.0
c        wtsfc(1)=0.00
c        wtsfc(1)=4.96e-07
         wtsfc(1)=1.190476e-06
         qstar(1)=wtsfc(1)
c        dtdzf(1)=0.000
         dtdzf(1)=0.2548
         dtjump  = 0.
         divgls  = 0.
         zo      = 0.0001
         zi      = -5.
c        izi     = (55*nnz)/100
c        izi     = nnz
         izi     = 55
         xl      = 50.
         yl      = 50.
         zl      = -20.
c
c ---------- if stretched grid specify location of first point
c
         zw1 = -0.5
      else
         gcp     = grav/Cpa
         batag   = bfac*grav/t00
         fcor_h  = 0.0

         !!!If doing DNS, BCs are set somewhat differently than when doing LES
         if (iDNS .eq. 1) then

         !Get BCs of temperature in Celsius for use in specific expressions
         topC=(Ttop(1)-273.15)
         botC=(Tbot(1)-273.15)  
         
         !Normalize the RH given in input file
         RHT = Ttop(2)/100.0 
         RHB = Tbot(2)/100.0

         !Convert RH given in input file into specific humidity
         Ttop(2) =  RHT*Mw/Ru/Ttop(1)*610.94*
     +               EXP(17.6257*topC/(topC+243.04))/rhoa
         Tbot(2) =  RHB*Mw/Ru/Tbot(1)*610.94*
     +               EXP(17.6257*botC/(botC+243.04))/rhoa

         wtsfc(1) = qstar(1)
         wtsfc(2) = qstar(2)

         else !(doing LES)
c
         !qstar specified in input file for isfc = 0
        
         wtsfc(1) = qstar(1)
         wtsfc(2) = qstar(2)
         wtsfc(3) = qstar(3)

         !tsfcc specified in input file for isfc = 1

         !Assuming tsfcc(2) is giving RH:

!         botC=(tsfcc(1)-273.15)  
!         RHB = tsfcc(2)/100.0
!
!         tsfcc(2) =  RHB*Mw/Ru/tsfcc(1)*610.94*
!     +               EXP(17.6257*botC/(botC+243.04))/rhoa
!

         end if
c
c
c
         dtjump  = 0.0
         divgls  = 0.0

  
c
c ---------- if stretched grid specify location of first point
c
         !Particle channel:
         izi     = (60*nnz)/100
      endif
c
      time  = 0.0
c 
c ---------- outermost coarse grid  indicies are bounds of grid
c
      izlow = 1
      izup  = nnz
      dz    = zl/nnz
      dzg   = abs(dz)
      if(l_root) write(6,4040) zl,nnz,dzg
c
c --------------- generate z grids for particular mesh from
c                 iz = 0,1,...,nnz+1; this allows indexing
c                 to array elements z(0), etc.
c
      zwstrt = 0.0
c
c ------------ if uniform vertical spacing then
c
      if(iz_space .eq. 0) then
c
c ------------ build z grid for w points
c
         do iz=0,nnz+1
            z(iz) = dz*float(iz) + zwstrt
         enddo
      else
        call vgrid(zw1,zi,zl,nnz,z(0),l_root,l_debug)
        !call vgrid_channel(zw1,zi,zl,nnz,z(0),l_root,l_debug)
      endif
c
      call get_dz
c
      if(l_root) then
         write(6,8002) zwstrt
         write(6,8003) (iz,z(iz),zz(iz),iz=0,nnz+1)
      endif
c
      nnzm1 = nnz-1
      dx    = xl/nnx
      dy    = yl/nny
      fnxy  = 1./float(nxy)
      dzdz  = dzw(1)*dzw(1)
      z1    = zz(1)
c
      c23  = 2.0/3.0
      dsl  = (dx*1.5*dy*1.5*abs(dzw(1)))**(1./3.)
      dslg = dsl
      cs   = 0.2
c
      vk     = 0.4
      batagk = batag*vk
      vkin   = 1./vk
      ttmean = 0.
      zody   = alog(abs(z1/zo))
      write(nprt, 9901) z1,zo,zody
 9901 format(' 9901 z1 = ',e15.6,' zo = ',e15.6,/,
     +       ' zody = ',e15.6)
      zodyin = 1./zody
      wstar  = abs(batag*zi*wtsfc(1))**(1./3.)
      if(ismlt .eq. 1) then
c
c ---- set constants for businger similarity functions
c
         vk74   = vk*0.74
         vk74in = 0.74/vk
         zody74 = zody*0.74
      else 
c
c ---- set constants for large similarity functions
c
        vk74    = vk
        vk74in  = 1.0/vk
        zody74  = zody
      endif
      ugal   = 0.0
c      ugal   = ugcont*0.5
c     ugcont = ugcont - ugal
      cdbtm  = vk*vk/zody/zody
      if(iocean .eq. 1) then
c ----------- set surface friction velocity here and in sr. sufto
c        utau = 4.29e-03
         utau = 7.00e-03
      else
         ufree = 0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)
c
c ---- note : new estimate for utau !!!
c
         utau  = vk*(ufree+ugcont)/zody
c        utau  = vk*(ufree)/zody
      endif
      utau2    = utau*utau
      if(ibuoy .eq. 0 .or. qstar(1) .eq. 0.) then
        amonin = 1000.0
      else
        amonin = -utau2*utau/(batagk*qstar(1))
      endif
      hol   = abs(zi)/amonin
      zol   = abs(z1)/amonin
      uwsfc = -utau*utau
      vwsfc = -utau*utau
c
      if(l_root) then
         write(6,80)
         write(6,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +         ,cdbtm,ugcont
      endif
c
      if(l_debug) then
         write(nprt,80)
         write(nprt,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +         ,cdbtm,ugcont
      endif
c

c     Virtual probe placement/locator 6/27/2018 JP
      probe_loc = -1
      probe_id  = -1
      probe_z_locations=(/ 10, 300, 500 /)
      call probe_locator(probe_z_locations(1))
      call prob_proc(nprobe1,iploc)
      call probe_locator(probe_z_locations(2))
      call prob_proc(nprobe2,iploc)
      call probe_locator(probe_z_locations(3))
      call prob_proc(nprobe3,iploc)
      if(myid==0)then
        write(*,*)'first flux plane locagtion z: ',0.5*zl/nnz
      end if

      write(wname,"(A,I0.3,A)")'part',myid,'.out'
      open(nparts+myid,file=wname,form='unformatted',
     +     status='replace') 

      return
c ------------------------
   2  format(10x,' WT =',e12.4,',  U* =',e12.4,',  L =',e12.4,/,
     +       10x,' DTDZ FREE =',e12.4,',  ZODY=',e12.4,/,10x,
     +       ' ZO(BTM) =',e12.4,',  CDBTM=',e12.4,
     +       ',  UG = ',e12.4)
  80  format(///,' ***** SCRATCH RUN ***** ',//)
 4040 format(' zl = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 4043 format(' znest = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 8002 format(' zwstrt = ',e12.4)
 8003 format(' iz ',5x,' zw',5x,' zu ',5x,/,(i3,2e12.4))
      end
