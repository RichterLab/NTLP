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

  integer :: i,j,k,kmin,kmax,jmin,jmax,ip
  integer :: myny,mynz,mynp,npsize,nps,npe
  integer :: rankarray(numprocs),nparray(numprocs)
  integer :: ierr
  integer :: teci,isdouble,filetype,debug,fileformat
  integer :: tecini142,teczne142,tecdat142,tecend142,tecmpiinit142,tecznemap142,tecijkptn142
  character*9 :: output_name
  character*160 :: data_names
  character*1 :: nullchr
  integer :: sharevar(10),null(10),nodeloc(10),passive(10)

  real,allocatable :: uplt(:,:,:),vplt(:,:,:),wplt(:,:,:),wtmp(:,:,:),tplt(:,:,:),qplt(:,:,:),xplt(:,:,:),yplt(:,:,:),zplt(:,:,:)
  real*4,allocatable :: toplt(:,:,:),xpart(:),ypart(:),zpart(:),upart(:),vpart(:),wpart(:),tpart(:),rpart(:)

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


  data_names = 'x y z u v w th q tp rp'

  nullchr = char(0)
  null = 0
  sharevar = 0
  passive = 1
  passive(1:8) = 0
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

  passive = 0
  passive(7:8) = 1

  mynp = numpart  !From the last time particle number was calculated

  if (myid .eq. 0) then
     npsize=mynp
  else
     npsize=mynp+1
  end if

  !allocate(xpart(npsize),ypart(npsize),zpart(npsize),upart(npsize),vpart(npsize),wpart(npsize),tpart(npsize),rpart(npsize))
  allocate(xpart(tnumpart),ypart(tnumpart),zpart(tnumpart),upart(tnumpart),vpart(tnumpart),wpart(tnumpart),tpart(tnumpart),rpart(tnumpart))


  !Need the beginning and end index of particle number for each proc
  call mpi_allgather(mynp,1,mpi_integer,nparray,1,mpi_integer,mpi_comm_world,ierr)

  if (myid .eq. 0) then
     nps=1
     npe = mynp
  else
    nps = 0
    do ip = 1,myid
       nps = nps + nparray(ip) 
    end do
    npe = nps + mynp
  end if

  !call fill_plt_part(npsize,xpart,ypart,zpart,upart,vpart,wpart,tpart,rpart)

  call fill_plt_part_ROOT(nparray,xpart,ypart,zpart,upart,vpart,wpart,tpart,rpart)

  if (myid .eq. 0) then


  teci = teczne142('PART'//nullchr,&
       0,&              !0 = ordered data
       tnumpart,&       !num elements in first direction
       1,&              !num elements in second direction
       1,&              !num elements in third direction
       0,&              !ICellMax (don't use)
       0,&              !JCellMax (don't use)
       0,&              !KCellMax (don't use)
       null,&           !time stamp
       0,&              !strandid
       0,&              !parent zone
       1,&              !0=point, 1=block
       0,&              !numfaceconnections
       0,&              !faceneighbormode
       0,&              !totalnumfacenodes
       0,&              !numconnectedboundaryfaces
       0,&              !totalnumboundaryconnections
       passive,&        !passivevarlist
       nodeloc,&        !valuelocation, null means centered
       sharevar,&       !sharevarzone, share variable from zone
       0)  

  teci = tecznemap142(1,0)

  !teci = tecijkptn142(myid+1,&
  !     nps,&
  !     1,&
  !     1,&
  !     npe,&
  !     1,&
  !     1)

  teci = tecdat142(tnumpart,xpart,isdouble)
  teci = tecdat142(tnumpart,ypart,isdouble)
  teci = tecdat142(tnumpart,zpart,isdouble)
  teci = tecdat142(tnumpart,upart,isdouble)
  teci = tecdat142(tnumpart,vpart,isdouble)
  teci = tecdat142(tnumpart,wpart,isdouble)
  teci = tecdat142(tnumpart,tpart,isdouble)
  teci = tecdat142(tnumpart,rpart,isdouble)

  end if

  teci = tecend142()  

  deallocate(uplt,vplt,wplt,wtmp,tplt,qplt,xplt,yplt,zplt,toplt)
  deallocate(xpart,ypart,zpart,upart,vpart,wpart,tpart,rpart)

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

end subroutine fill_plt_fields
subroutine fill_plt_part(npsize,xpart,ypart,zpart,upart,vpart,wpart,tpart,rpart)
  use pars
  use particles
  implicit none
  include 'mpif.h'

  integer :: ip,npsize,nsend,nrecv
  integer :: procnext,procbefore
  integer :: ierr,istatus(mpi_status_size)
  real :: sendbuf(8),recvbuf(8)
  real*4, intent(inout) :: xpart(npsize),ypart(npsize),zpart(npsize),upart(npsize),vpart(npsize),wpart(npsize),tpart(npsize),rpart(npsize)

  !Proc 0 has no "halo", sends its last one to proc 1, etc.
  if (myid .eq. 0) then
     ip = 1
  else
     ip = 2
  end if

  part => first_particle
  do while (associated(part))

    xpart(ip) = part%xp(1)
    ypart(ip) = part%xp(2)
    zpart(ip) = part%xp(3)
    upart(ip) = part%vp(1)
    vpart(ip) = part%vp(2)
    wpart(ip) = part%vp(3)
    tpart(ip) = part%Tp
    rpart(ip) = part%radius
     
    ip = ip+1
  part => part%next
  end do

  !Need to get "halo" of particles by sending to procnext
  procnext = myid+1
  procbefore = myid-1

  if (myid==0) then
     procbefore = mpi_proc_null
  end if
  if (myid==numprocs-1) then
     procnext = mpi_proc_null
  end if

  sendbuf(1) = xpart(npsize)
  sendbuf(2) = ypart(npsize)
  sendbuf(3) = zpart(npsize)
  sendbuf(4) = upart(npsize)
  sendbuf(5) = vpart(npsize)
  sendbuf(6) = wpart(npsize)
  sendbuf(7) = tpart(npsize)
  sendbuf(8) = rpart(npsize)

  nsend = 8
  nrecv = nsend
   
  !Send to procnext, receive from procbefore
  call mpi_sendrecv(sendbuf(1),nsend,mpi_real8,procnext,0, &
                    recvbuf(1),nrecv,mpi_real8,procbefore,0, &
                    mpi_comm_world,istatus,ierr)

  if (myid .ne. 0) then
     xpart(1) = recvbuf(1)
     ypart(1) = recvbuf(2)
     zpart(1) = recvbuf(3)
     upart(1) = recvbuf(4)
     vpart(1) = recvbuf(5)
     wpart(1) = recvbuf(6)
     tpart(1) = recvbuf(7)
     rpart(1) = recvbuf(8)
   end if
   

end subroutine fill_plt_part
subroutine fill_plt_part_ROOT(nparray,xpart,ypart,zpart,upart,vpart,wpart,tpart,rpart)
  use pars
  use particles
  implicit none
  include 'mpif.h'

  integer :: i,ip,nparray(numprocs),npsarray(numprocs)
  integer :: ierr
  real*4, intent(inout) :: xpart(:),ypart(:),zpart(:),upart(:),vpart(:),wpart(:),tpart(:),rpart(:)
  real*4, allocatable :: myxpart(:),myypart(:),myzpart(:),myupart(:),myvpart(:),mywpart(:),mytpart(:),myrpart(:)

  allocate(myxpart(numpart),myypart(numpart),myzpart(numpart),myupart(numpart),myvpart(numpart),mywpart(numpart),mytpart(numpart),myrpart(numpart))

  !Everyone collects into "my" arrays
  ip = 1
  part => first_particle
  do while (associated(part))

    myxpart(ip) = part%xp(1)
    myypart(ip) = part%xp(2)
    myzpart(ip) = part%xp(3)
    myupart(ip) = part%vp(1)
    myvpart(ip) = part%vp(2)
    mywpart(ip) = part%vp(3)
    mytpart(ip) = part%Tp
    myrpart(ip) = part%radius

    ip = ip+1
  part => part%next
  end do

!  !Start collecting, starting with root
!  if (myid .eq. 0) then
!     xpart(1:numpart) = myxpart(1:numpart)
!     ypart(1:numpart) = myypart(1:numpart)
!     zpart(1:numpart) = myzpart(1:numpart)
!     upart(1:numpart) = myupart(1:numpart)
!     vpart(1:numpart) = myvpart(1:numpart)
!     wpart(1:numpart) = mywpart(1:numpart)
!     tpart(1:numpart) = mytpart(1:numpart)
!     rpart(1:numpart) = myrpart(1:numpart)
!  end if

  ip = 1
  do i=1,numprocs

     if (i==1) then
        npsarray(i) = 0
     else
        npsarray(i) = npsarray(i-1)+nparray(i-1)
     end if

  end do

  call mpi_gatherv(myxpart,numpart,mpi_real4,xpart,nparray,npsarray,mpi_real4,0,mpi_comm_world,ierr)
  call mpi_gatherv(myypart,numpart,mpi_real4,ypart,nparray,npsarray,mpi_real4,0,mpi_comm_world,ierr)
  call mpi_gatherv(myzpart,numpart,mpi_real4,zpart,nparray,npsarray,mpi_real4,0,mpi_comm_world,ierr)
  call mpi_gatherv(myupart,numpart,mpi_real4,upart,nparray,npsarray,mpi_real4,0,mpi_comm_world,ierr)
  call mpi_gatherv(myvpart,numpart,mpi_real4,vpart,nparray,npsarray,mpi_real4,0,mpi_comm_world,ierr)
  call mpi_gatherv(mywpart,numpart,mpi_real4,wpart,nparray,npsarray,mpi_real4,0,mpi_comm_world,ierr)
  call mpi_gatherv(mytpart,numpart,mpi_real4,tpart,nparray,npsarray,mpi_real4,0,mpi_comm_world,ierr)
  call mpi_gatherv(myrpart,numpart,mpi_real4,rpart,nparray,npsarray,mpi_real4,0,mpi_comm_world,ierr)

  deallocate(myxpart,myypart,myzpart,myupart,myvpart,mywpart,mytpart,myrpart)
  

end subroutine fill_plt_part_ROOT


end module tec_io
