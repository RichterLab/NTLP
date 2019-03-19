
c     This is Monahan's SSGF
c     Currently hard coded into doing 12m/s u10
      function SSGF0(rp80)
      real,intent(in) :: rp80
      real :: SSGF0

      SSGF0 = 1.373*(12**3.41)*rp80**(-3)*(1+0.057*r**1.05)
     +        *10**(1.19*exp(-((0.380-log10(rp80)/0.650))**2))

      end function

      function SSGF1(rp80)
      real,intent(in) :: rp80
      real :: SSGF1
      SSGF1 = -log(rp80**3)
      end function

c     Using a gaussian distribution
      function SSGF2(rp80)
      real,intent(in) :: rp80
      real :: SSGF2
      real :: sigma, pi
      pi = 4.0*atan(1.0)

      sigma=1e-6

      SSGF2 = 1/sqrt(2*pi*sigma**2)*exp(-(rp80-5e-6)**2/(2*sigma**2))
      end function

      !Calculates the truncated Cumulative Density Function
      !Must make sample radius size large enough, if not
      !then radius max/min might not be at close enough
      !locations to consider
      subroutine calc_truncated_CDF
      use pars
      use fields
      use con_data
      use particles
      real, allocatable :: rpinit(:)
      !Sn,Snscl,rpcrop,CDF are global values
      integer :: irad_max,irad_min
      integer :: m,n,i,j
      real    :: SSGF
      m           = 3000 !sample radius size
      allocate(rpinit(m+2))
      rpinit(1)   = 0.0
c      Tianze's original rpinit designation
c      rpinit(2:m) = (/((5e-7*i),i=1,m)/)

      !my code for rpinit
      rpinit(2:m+2)=(/((radius_min+(radius_max-radius_min)*(i-1)/m),i=1,
     + m+2)/)
      allocate(Sn(m+2))
      !Non-scaled SSGF up to this point
      do i=1,m+2
      !Select an SSGF based on iSSGF (Control this!)
        if(iSSGF==0) then
          SSGF=SSGF0(rpinit(i)/1e-6) !Monahan's relation!
        else if (iSSGF==1) then
          SSGF=SSGF1(rpinit(i))
        else if (iSSGF==2) then
          SSGF=SSGF2(rpinit(i))
        end if
        Sn(i) = SSGF
      end do
      if (myid==0) write(*,*)'iSSGF: '    ,iSSGF
      if (myid==0) write(*,*)'rpinit(2): ',rpinit(2)
      if (myid==0) write(*,*)'Sn(2): '    ,Sn(2)

      !Find index for useful range
      irad_max=minloc(rpinit,DIM=1,MASK=(rpinit>=radius_max))
      irad_min=maxloc(rpinit,DIM=1,MASK=(rpinit<=radius_min))

      !Now scale the SSGF
      n=size(Sn(irad_min:irad_max))
      allocate(Snscl(n),rpcrop(n),CDF(n))
      rpcrop=rpinit(irad_min:irad_max)
      Snscl=Sn(irad_min:irad_max)/sum(Sn(irad_min:irad_max))
      CDF(1) = Snscl(1)

      !Print out the scaled quantities:
c      if(myid==0) write(*,*)'rpcrop(2): ',rpcrop(2)
      if(myid==0) write(*,*)'Sn: ',Sn
c      if(myid==0) write(*,*)'sum of Sn crop:',
c     + sum(Sn(irad_min:irad_max))

      do j=2,n
        CDF(j)=CDF(j-1)+Snscl(j)
      end do
      end
