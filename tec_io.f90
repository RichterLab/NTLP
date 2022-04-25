module tec_io
implicit none


integer :: plt_counter
character(len=:),allocatable :: path_plt

CONTAINS


subroutine plt_fields
  use particles
  use fields
  use pars
  use con_stats
  use con_data
  implicit none
  include 'mpif.h'

  integer :: i,j,k,kmin,kmax,jmin,jmax
  integer :: myny,mynz
  integer :: rankarray(numprocs)
  integer :: teci,isdouble,filetype,debug,fileformat
  integer :: tecini142,teczne142,tecdat142,tecend142,tecmpiinit142,tecznemap142,tecijkptn142
  character*9 :: output_name
  character*160 :: data_names
  character*1 :: nullchr
  integer :: sharevar(8),null(8),nodeloc(8),passive(8)

  real,allocatable :: uplt(:,:,:),vplt(:,:,:),wplt(:,:,:),wtmp(:,:,:),tplt(:,:,:),qplt(:,:,:),xplt(:,:,:),yplt(:,:,:),zplt(:,:,:)
  real*4,allocatable :: toplt(:,:,:)

  path_plt = trim(trim(path_seed)//"tec_viz.szplt")

  kmin = izs
  kmax = min(ize+1,nnz)
  jmin = max(iys-1,1)
  jmax = iye

  allocate(uplt(1:nnx,jmin:jmax,kmin:kmax),vplt(1:nnx,jmin:jmax,kmin:kmax),wplt(1:nnx,jmin:jmax,kmin:kmax),wtmp(1:nnx,jmin:jmax,kmin:kmax),tplt(1:nnx,jmin:jmax,kmin:kmax),qplt(1:nnx,jmin:jmax,kmin:kmax),xplt(1:nnx,jmin:jmax,kmin:kmax),yplt(1:nnx,jmin:jmax,kmin:kmax),zplt(1:nnx,jmin:jmax,kmin:kmax),toplt(1:nnx,jmin:jmax,kmin:kmax))

  call fill_plt_fields(jmin,jmax,kmin,kmax,uplt,vplt,wtmp,tplt,qplt)

  do k=kmin,kmax
  do j=jmin,jmax
  do i=1,nnx

     !if (k==1) then
     !   wplt(i,j,k) = 0.5*wtmp(i,j,k)
     !else
        wplt(i,j,k) = 0.5*(wtmp(i,j,k)+wtmp(i,j,k-1))
     !end if

  end do
  end do
  end do

  myny = jmax-jmin+1
  mynz = kmax-kmin+1

  do i=0,numprocs-1
     rankarray(i+1) = i
  end do

  !Fill x,y,z:
  do k=kmin,kmax
     do j=jmin,jmax
        do i=1,nnx
           xplt(i,j,k) = dx*(i-1)
           yplt(i,j,k) = dy*(j-1)
           zplt(i,j,k) = zz(k)
        end do
     end do
  end do


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
       path_plt//nullchr,&
       trim(path_seed)//nullchr,&
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
       jmin,&
       kmin,&
       nnx,&
       jmax,&
       kmax)

  toplt = xplt(1:nnx,jmin:jmax,kmin:kmax)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,jmin:jmax,kmin:kmax),isdouble)

  toplt = yplt(1:nnx,jmin:jmax,kmin:kmax)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,jmin:jmax,kmin:kmax),isdouble)

  toplt = zplt(1:nnx,jmin:jmax,kmin:kmax)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,jmin:jmax,kmin:kmax),isdouble)

  toplt = uplt(1:nnx,jmin:jmax,kmin:kmax)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,jmin:jmax,kmin:kmax),isdouble)

  toplt = vplt(1:nnx,jmin:jmax,kmin:kmax)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,jmin:jmax,kmin:kmax),isdouble)

  toplt = wplt(1:nnx,jmin:jmax,kmin:kmax)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,jmin:jmax,kmin:kmax),isdouble)

  toplt = tplt(1:nnx,jmin:jmax,kmin:kmax)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,jmin:jmax,kmin:kmax),isdouble)

  toplt = qplt(1:nnx,jmin:jmax,kmin:kmax)
  teci = tecdat142(nnx*myny*mynz,toplt(1:nnx,jmin:jmax,kmin:kmax),isdouble)


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

  deallocate(uplt,vplt,wplt,wtmp,tplt,qplt,xplt,yplt,zplt,toplt)

end subroutine plt_fields
subroutine fill_plt_fields(jmin,jmax,kmin,kmax,uplt,vplt,wtmp,tplt,qplt)
use pars
use fields
implicit none
include 'mpif.h'

real, intent(inout) :: uplt(1:nnx,jmin:jmax,kmin:kmax),vplt(1:nnx,jmin:jmax,kmin:kmax),wtmp(1:nnx,jmin:jmax,kmin:kmax),tplt(1:nnx,jmin:jmax,kmin:kmax),qplt(1:nnx,jmin:jmax,kmin:kmax)
integer, intent(in) :: jmin,jmax,kmin,kmax
integer :: yfore,yback,nsend,nrecv
integer :: ix,iy,iz,jloc,iscl
integer :: ierr,istatus(mpi_status_size)
real :: sendbuf(1:nnx,kmin:kmax,6),recvbuf(1:nnx,kmin:kmax,6)

  !First fill everything except y-halo
  do iz=kmin,kmax
  do iy=iys,iye
  do ix=1,nnx
     uplt(ix,iy,iz) = u(ix,iy,iz)
     vplt(ix,iy,iz) = v(ix,iy,iz)
     wtmp(ix,iy,iz) = w(ix,iy,iz)
     tplt(ix,iy,iz) = t(ix,iy,1,iz)
     qplt(ix,iy,iz) = t(ix,iy,2,iz)
  end do
  end do
  end do

  !To fill the halos in the y direction, need to figure out the proc numbers ahead and behind myid

  if ( mod(myid+1,ncpu_s) == 0 ) then  !At max y
     !yfore = myid-ncpu_s+1
     yfore = mpi_proc_null
     yback = myid-1
  elseif (mod(myid,ncpu_s) == 0) then  !At min y
     yfore = myid+1
     !yback = myid+ncpu_s-1
     yback = mpi_proc_null
  else
     yfore = myid+1
     yback = myid-1
  endif

  nsend = nnx*(kmax-kmin+1)*6
  nrecv = nsend

  !Send to right, receive from left -- that's all that tecplot needs
  do iz=kmin,kmax
  do ix=1,nnx
     sendbuf(ix,iz,1) = u(ix,iye,iz)
     sendbuf(ix,iz,2) = v(ix,iye,iz)
     sendbuf(ix,iz,3) = w(ix,iye,iz)
     sendbuf(ix,iz,4) = e(ix,iye,iz)
     sendbuf(ix,iz,5) = t(ix,iye,1,iz)
     sendbuf(ix,iz,6) = t(ix,iye,2,iz)
  end do
  end do

    
  call mpi_sendrecv(sendbuf(1,kmin,1),nsend,mpi_real8,yfore,0, &
                    recvbuf(1,kmin,1),nrecv,mpi_real8,yback,0, &
                    mpi_comm_world,istatus,ierr)

  if (yback .ne. mpi_proc_null) then

  do iz=kmin,kmax
  do ix=1,nnx
     uplt(ix,iys-1,iz) = recvbuf(ix,iz,1)
     vplt(ix,iys-1,iz) = recvbuf(ix,iz,2)
     wtmp(ix,iys-1,iz) = recvbuf(ix,iz,3)
     tplt(ix,iys-1,iz) = recvbuf(ix,iz,5)
     qplt(ix,iys-1,iz) = recvbuf(ix,iz,6)
  end do
  end do

  end if
  

!  !Send to left, receive from right
!  do iz=izs,ize
!  do ix=1,nnx
!     sendbuf(ix,iz,1) = u(ix,iys,iz)
!     sendbuf(ix,iz,2) = v(ix,iys,iz)
!     sendbuf(ix,iz,3) = w(ix,iys,iz)
!     sendbuf(ix,iz,4) = e(ix,iys,iz)
!  end do
!  end do
!
!  do iscl=1,nscl
!     jloc = 4+iscl
!     do iz=izs,ize
!     do ix=1,nnx
!        sendbuf(ix,iz,jloc) = t(ix,iys,iscl,iz)
!     end do
!     end do
!  end do
!    
!  call mpi_sendrecv(sendbuf(1,izs,1),nsend,mpi_real8,yback,0, &
!                    recvbuf(1,izs,1),nrecv,mpi_real8,yfore,0, &
!                    mpi_comm_world,istatus,ierr)
!
!  do iz=izs,ize
!  do ix=1,nnx
!     uplt(ix,iye+1,iz) = recvbuf(ix,iz,1)
!     vplt(ix,iye+1,iz) = recvbuf(ix,iz,2)
!     wtmp(ix,iye+1,iz) = recvbuf(ix,iz,3)
!     tplt(ix,iye+1,iz) = recvbuf(ix,iz,5)
!     qplt(ix,iye+1,iz) = recvbuf(ix,iz,6)
!  end do
!  end do
  


end subroutine fill_plt_fields

end module tec_io
