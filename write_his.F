      subroutine write_his(iloc)
c
c ----- write history file with global parameters
c       write tsfcc specially to preserve digits!
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
c
      divgmax = 0.0
      do iz=1,nnz
         divgmax = amax1(divgmax, divz(iz))
      enddo
c
      ziavg = zi
      holtop = hol
      wt_min = wtsb(iloc,1)
      wt_le  = wtle(iloc,1)
      krec = krec + 1
      mid = nnz/4
      write(nhis1,6000) time,dt,utau,ziavg,amonin,holtop,
     +         (tsfcc(1)-t_ref),uwsfc,vwsfc,divgmax, wt_min, wt_le,
     +         ucfl, vcfl, wcfl, wtsfc(1),
     +         ups(mid),vps(mid),wps(mid),tps(mid,1),
     +         uwle(mid),uwsb(mid),uw_tot(mid),
     +         vwle(mid),vwsb(mid),vw_tot(mid),
     +         wtle(mid,1),wtsb(mid,1),wt_tot(mid,1),
     +         englez(mid),eavg(mid), wabs,
     +         maxval(partcount(1:nnx,iys:iye,1)),
     +         Rep_avg,phip,phiw,phiv,
c    +         tps(mid,2), tps(mid,3),
c    +         wtle(mid,2),wtsb(mid,2),wt_tot(mid,2),
c    +         wtle(mid,3),wtsb(mid,3),wt_tot(mid,3)

c     Eddy Covariance method writing stuff out on ascii
     +         t(64,64,3,2),w(64,64,2),partcount(64,64,2)
     
c    Instantaneous statistics at some height in the boundary layer
 

c       format(number of quantities,e17.8 output style)
c 6000 format(5e17.8)
 6000 format(40e17.8)
c
c -------------- write profile information
c
      call write_prof(nhisp,krec,isize,c_s%wwsb)
c
      return
      end
