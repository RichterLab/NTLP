module tec_io
implicit none


integer :: plt_counter
character(len=80) :: path_plt

CONTAINS


subroutine plt_fields
  use particles
  use fields
  use pars
  use con_stats
  use con_data
  implicit none
  include 'mpif.h'

  integer :: i,j,k
  integer :: myny,mynz
  integer :: rankarray(numprocs)
  integer :: teci,isdouble,filetype,debug,fileformat
  integer :: tecini142,teczne142,tecdat142,tecend142,tecmpiinit142,tecznemap142,tecijkptn142
  character*9 :: output_name
  character*160 :: data_names
  character*1 :: nullchr
  integer :: sharevar(8),null(8),nodeloc(8),passive(8)

  real :: wplt(1:nnx,iys:iye,izs:ize),x(1:nnx,iys:iye,izs:ize),y(1:nnx,iys:iye,izs:ize),zplt(1:nnx,iys:iye,izs:ize)
  real*4 :: toplt(1:nnx,iys:iye,izs:ize)

  path_plt = trim(adjustl(trim(adjustl(path_seed))//"dummy.plt"))

  do k=izs,ize
  do j=iys,iye
  do i=1,nnx

     if (k==1) then
        wplt(i,j,k) = 0.5*w(i,j,k)
     else
        wplt(i,j,k) = 0.5*(w(i,j,k)+w(i,j,k-1))
     end if

  end do
  end do
  end do

  myny = iye-iys+1
  mynz = ize-izs+1

  do i=0,numprocs-1
     rankarray(i+1) = i
  end do

  !Fill x,y,z:
  do k=izs,ize
     do j=iys,iye
        do i=1,nnx
           x(i,j,k) = dx*(i-1)
           y(i,j,k) = dy*(j-1)
           zplt(i,j,k) = zz(k)
        end do
     end do
  end do


  output_name = 'dummy.plt'
  data_names = 'x y z u v w th q'

  nullchr = char(0)
  null = 0
  sharevar = 0
  passive = 0
  nodeloc = 1
  isdouble = 0
  fileformat = 1
  filetype = 0
  debug = 1

  teci = tecini142('veldata'//nullchr,&
       data_names//nullchr,&
       !output_name//nullchr,&
       path_plt//nullchr,&
       trim((path_seed)//nullchr,&
       fileformat,&
       filetype,&
       debug,&
       isdouble)

  teci = tecmpiinit142(mpi_comm_world,0)

  teci = teczne142('FLOW'//nullchr,&
       0,&              !0 = ordered data
       nnx,&             !num elements in first direction
       nny,&             !num elements in second direction
       nnz,&             !num elements in third direction
       0,&              !ICellMax (don't use)
       0,&              !JCellMax (don't use)
       0,&              !KCellMax (don't use)
       null,&           !time stamp
       0,&              !strandid
       0,&              !parent zone
       1,&              !0=point, 1=block
       0,&              !numfaceconnections
       0,&              !faceneighbormode
       1,&              !totalnumfacenodes
       1,&              !numconnectedboundaryfaces
       1,&              !totalnumboundaryconnections
       passive,&        !passivevarlist
       nodeloc,&        !valuelocation, null means centered
       sharevar,&       !sharevarzone, share variable from zone
       0)

  
  teci = tecznemap142(numprocs,rankarray)

  teci = tecijkptn142(myid+1,&
       1,&
       iys,&
       izs,&
       nnx,&
       iye,&
       ize)

  toplt = x(1:nnx,iys:iye,izs:ize)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,iys:iye,izs:ize),isdouble)

  toplt = y(1:nnx,iys:iye,izs:ize)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,iys:iye,izs:ize),isdouble)

  toplt = zplt(1:nnx,iys:iye,izs:ize)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,iys:iye,izs:ize),isdouble)

  toplt = u(1:nnx,iys:iye,izs:ize)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,iys:iye,izs:ize),isdouble)

  toplt = v(1:nnx,iys:iye,izs:ize)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,iys:iye,izs:ize),isdouble)

  toplt = wplt(1:nnx,iys:iye,izs:ize)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,iys:iye,izs:ize),isdouble)

  toplt = t(1:nnx,iys:iye,1,izs:ize)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,iys:iye,izs:ize),isdouble)

  toplt = t(1:nnx,iys:iye,2,izs:ize)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,iys:iye,izs:ize),isdouble)


  !Now get the particle locations to put into the *.plt file:

!  passive = 1
!  passive(1:7) = 0
!
!  teci = teczne142('PART'//nullchr,&
!       0,&              !0 = ordered data
!       np,&             !num elements in first direction
!       1,&              !num elements in second direction
!       1,&              !num elements in third direction
!       0,&              !ICellMax (don't use)
!       0,&              !JCellMax (don't use)
!       0,&              !KCellMax (don't use)
!       null,&           !time stamp
!       0,&              !strandid
!       0,&              !parent zone
!       1,&              !0=point, 1=block
!       0,&              !numfaceconnections
!       0,&              !faceneighbormode
!       0,&              !totalnumfacenodes
!       0,&              !numconnectedboundaryfaces
!       0,&              !totalnumboundaryconnections
!       passive,&        !passivevarlist
!       nodeloc,&        !valuelocation, null means centered
!       sharevar,&       !sharevarzone, share variable from zone
!       0)  
!
!  teci = tecdat142(np,xp(1:np),isdouble)
!  teci = tecdat142(np,yp(1:np),isdouble)
!  teci = tecdat142(np,zp(1:np),isdouble)
!  teci = tecdat142(np,upart(1:np),isdouble)
!  teci = tecdat142(np,vpart(1:np),isdouble)
!  teci = tecdat142(np,wpart(1:np),isdouble)
!  teci = tecdat142(np,rad(1:np),isdouble)

  teci = tecend142()  


end subroutine plt_fields

end module tec_io
