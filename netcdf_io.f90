module netcdf_io
implicit none

integer :: ncid
integer :: time_dimid,zu_varid,zu_dimid,s_dimid
integer :: time_varid,dt_varid
integer :: zw_varid,zw_dimid
integer :: uxym_varid,vxym_varid,wxym_varid,txym_varid
integer :: ups_varid,vps_varid,wps_varid,tps_varid
integer :: uwle_varid,uwsb_varid
integer :: dimids(1),dimids_zu(2),dimids_zw(2),dimids_zu_s(3),dimids_zw_s(3)
integer :: his_counter
character(len=80) :: path_netcdf_his

CONTAINS

subroutine netcdf_check(status)
      use netcdf
      implicit none
      integer, intent ( in) :: status

      if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
      end if

end subroutine netcdf_check

subroutine netcdf_init
      use netcdf
      use pars
      implicit none

      path_netcdf_his = "/scratch365/drichte2/tutorial/case1/netcdf_test.nc"

      call netcdf_check( nf90_create(path_netcdf_his,nf90_clobber,ncid))

      call netcdf_check( nf90_def_dim(ncid, "time",NF90_UNLIMITED, time_dimid) )

      call netcdf_check( nf90_def_dim(ncid, "zu",nnz, zu_dimid) )
      call netcdf_check( nf90_def_dim(ncid, "zw",nnz+1, zw_dimid) )

      call netcdf_check( nf90_def_dim(ncid,"nscl",nscl,s_dimid) )

      dimids = (/ time_dimid /)
      dimids_zu = (/ zu_dimid, time_dimid/)
      dimids_zw = (/ zw_dimid, time_dimid/)
      dimids_zu_s = (/ zu_dimid, s_dimid, time_dimid/)
      dimids_zw_s = (/ zw_dimid, s_dimid, time_dimid/)

      call netcdf_check( nf90_def_var(ncid,"time",NF90_REAL,dimids,time_varid) )
      call netcdf_check( nf90_put_att(ncid,time_varid,"title","Simulation time") )

      call netcdf_check( nf90_def_var(ncid, "dt", NF90_REAL, dimids,dt_varid) )
      call netcdf_check( nf90_put_att(ncid,dt_varid,"title","Model time step") )

      call netcdf_check( nf90_def_var(ncid, "zu", NF90_REAL, dimids_zu,zu_varid) )
      call netcdf_check( nf90_put_att(ncid,zu_varid,"title","z levels at u-points") )

      call netcdf_check( nf90_def_var(ncid, "zw", NF90_REAL, dimids_zw,zw_varid) )
      call netcdf_check( nf90_put_att(ncid,zw_varid,"title","z levels at w-points") )

      call netcdf_check( nf90_def_var(ncid,"uxym",NF90_REAL, dimids_zu,uxym_varid) )
      call netcdf_check( nf90_put_att(ncid,uxym_varid,"title","Horiz. avg. u vel") )

      call netcdf_check( nf90_def_var(ncid,"vxym",NF90_REAL, dimids_zu,vxym_varid) )
      call netcdf_check( nf90_put_att(ncid,vxym_varid,"title","Horiz. avg. v vel") )

      call netcdf_check( nf90_def_var(ncid,"wxym",NF90_REAL, dimids_zw,wxym_varid) )
      call netcdf_check( nf90_put_att(ncid,wxym_varid,"title","Horiz. avg. w vel") )

      call netcdf_check( nf90_def_var(ncid,"txym",NF90_REAL, dimids_zu_s,txym_varid) )
      call netcdf_check( nf90_put_att(ncid,txym_varid,"title","Horiz. avg. scalars") )

      call netcdf_check( nf90_def_var(ncid,"ups",NF90_REAL, dimids_zu,ups_varid) )
      call netcdf_check( nf90_put_att(ncid,ups_varid,"title","Fluctuating velocity <u'^2>") )

      call netcdf_check( nf90_def_var(ncid,"vps",NF90_REAL, dimids_zu,vps_varid) )
      call netcdf_check( nf90_put_att(ncid,vps_varid,"title","Fluctuating velocity <v'^2>") )

      call netcdf_check( nf90_def_var(ncid,"wps",NF90_REAL, dimids_zw,wps_varid) )
      call netcdf_check( nf90_put_att(ncid,wps_varid,"title","Fluctuating velocity <w'^2>") )

      call netcdf_check( nf90_def_var(ncid,"tps",NF90_REAL, dimids_zu_s,tps_varid) )
      call netcdf_check( nf90_put_att(ncid,tps_varid,"title","Fluctuating scalars <t'^2>") )

      call netcdf_check( nf90_def_var(ncid,"uwle",NF90_REAL, dimids_zw,uwle_varid) )
      call netcdf_check( nf90_put_att(ncid,uwle_varid,"title","Resolved <u'w'>") )

      call netcdf_check( nf90_def_var(ncid,"uwsb",NF90_REAL, dimids_zw,uwsb_varid) )
      call netcdf_check( nf90_put_att(ncid,uwsb_varid,"title","Subgrid <u'w'>") )
     

      call netcdf_check( nf90_enddef(ncid) )

      his_counter = 1

end subroutine netcdf_init

subroutine netcdf_restart
      use netcdf
      use pars
      implicit none

       

      path_netcdf_his = "/scratch365/drichte2/tutorial/case1/netcdf_test.nc"
      call netcdf_check( nf90_open(path_netcdf_his,NF90_WRITE,ncid) )

      !Could get the length of the time dimension to initialize his_counter
      !call netcdf_check( nf90_inq_dimlen(ncid,

      !... but realized that I would need to somehow populate all of the dimids
      !and varids! Can I loop through all of the vars in the file? How to do
      !this efficienetly?


end subroutine netcdf_restart

subroutine write_his_netcdf
      use netcdf
      use pars
      use fields
      use con_data
      use con_stats
      use particles

      implicit none
      real :: tmp(0:nnz)

      call netcdf_check( nf90_put_var(ncid, time_varid, real(time),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, dt_varid, real(dt),start=(/his_counter/)) )

      call netcdf_check( nf90_put_var(ncid, zu_varid, real(zz(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, zw_varid, real(z(0:nnz)),start=(/1, his_counter/)) )

      call netcdf_check( nf90_put_var(ncid,uxym_varid,real(uxym(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vxym_varid,real(vxym(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,wxym_varid,real(wxym(0:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,txym_varid,real(txym(1:nnz,1:nscl)),start=(/1,1,his_counter/)) )

      call netcdf_check( nf90_put_var(ncid,ups_varid,real(ups(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vps_varid,real(vps(1:nnz)),start=(/1, his_counter/)) )

      tmp(0) = 0.0
      tmp(1:nnz) = wps(1:nnz)
      call netcdf_check( nf90_put_var(ncid,wps_varid,real(tmp),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,tps_varid,real(tps(1:nnz,1:nscl)),start=(/1,1,his_counter/)) )

      tmp(0) = 0
      tmp(1:nnz) = uwle(1:nnz)
      call netcdf_check( nf90_put_var(ncid,uwle_varid,real(tmp),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,uwsb_varid,real(uwsb(0:nnz)),start=(/1, his_counter/)) )


      his_counter = his_counter + 1

end subroutine write_his_netcdf
subroutine close_his_netcdf
      use netcdf
      use pars
      implicit none

      call netcdf_check( nf90_close(ncid) )

end subroutine close_his_netcdf
subroutine open_his_netcdf
      use netcdf
      use pars
      implicit none

      call netcdf_check( nf90_open(path_netcdf_his,NF90_WRITE,ncid) )

end subroutine open_his_netcdf

end module netcdf_io
