module particles
      integer :: rproc,trproc,tproc,tlproc,lproc,blproc,bproc,brproc
      integer :: pr_r,pl_r,pt_r,pb_r,ptr_r,ptl_r,pbl_r,pbr_r
      integer :: pr_s,pl_s,pt_s,pb_s,ptr_s,ptl_s,pbl_s,pbr_s
      real :: ymin,ymax,zmin,zmax,xmax,xmin
      real, allocatable :: uext(:,:,:), vext(:,:,:), wext(:,:,:)
      real, allocatable :: u_t(:,:,:), v_t(:,:,:), w_t(:,:,:)
      real, allocatable :: Text(:,:,:),T_t(:,:,:)
      real, allocatable :: T2ext(:,:,:),T2_t(:,:,:)
      real, allocatable :: partTsrc(:,:,:),partTsrc_t(:,:,:)
      real, allocatable :: partHsrc(:,:,:),partHsrc_t(:,:,:) 
      real, allocatable :: partTEsrc(:,:,:),partTEsrc_t(:,:,:) 
      real, allocatable :: partcount_t(:,:,:),partsrc_t(:,:,:,:)
      real, allocatable :: vpsum_t(:,:,:,:),vpsqrsum_t(:,:,:,:)
      real, allocatable :: upwp_t(:,:,:),upwp(:,:,:)
      real, allocatable :: partcount(:,:,:),partsrc(:,:,:,:)
      real, allocatable :: vpsum(:,:,:,:),vpsqrsum(:,:,:,:)
      real, allocatable :: Tpsum(:,:,:),Tpsum_t(:,:,:)
      real, allocatable :: Tpsqrsum(:,:,:),Tpsqrsum_t(:,:,:)
      real, allocatable :: Tfsum(:,:,:),Tfsum_t(:,:,:)
      real, allocatable :: qfsum(:,:,:),qfsum_t(:,:,:)
      real, allocatable :: wpTpsum(:,:,:),wpTpsum_t(:,:,:)
      real, allocatable :: radsum(:,:,:),radsum_t(:,:,:) 
      real, allocatable :: rad2sum(:,:,:),rad2sum_t(:,:,:) 
      real, allocatable :: multcount(:,:,:),multcount_t(:,:,:) 
      real, allocatable :: mwsum(:,:,:),mwsum_t(:,:,:)
      real, allocatable :: qstarsum(:,:,:),qstarsum_t(:,:,:) 

!--- SFS velocity calculation ---------
      real, allocatable :: sigm_s(:,:,:),sigm_sdx(:,:,:),sigm_sdy(:,:,:)
      real, allocatable :: sigm_sdz(:,:,:),sigm_sext(:,:,:)
      real, allocatable :: sigm_sdxext(:,:,:),sigm_sdyext(:,:,:)
      real, allocatable :: sigm_sdzext(:,:,:)
      real, allocatable :: vis_ss(:,:,:),vis_sext(:,:,:) 

      integer :: particletype,pad_diff
      integer :: numpart,tnumpart,ngidx
      integer :: iseed
      integer :: num100=0, num1000=0, numimpos=0
      integer :: tnum100, tnum1000, tnumimpos
      integer :: denum, actnum, tdenum, tactnum
      integer :: num_destroy=0,tnum_destroy=0
      integer :: tot_reintro=0

      real :: Rep_avg,part_grav(3)
      real :: radavg,radmin,radmax,tempmin,tempmax,qmin,qmax
      real :: vp_init(3),Tp_init,radius_init
      real :: pdf_factor,pdf_prob
      integer*8 :: mult_init,mult_factor,mult_a,mult_c

      real :: avgres=0,tavgres=0

      integer, parameter :: histbins = 512
      real :: hist_rad(histbins+2)
      real :: bins_rad(histbins+2)      

      real :: hist_res(histbins+2)
      real :: bins_res(histbins+2)      

      
!REMEMBER: IF ADDING ANYTHING, MUST UPDATE MPI DATATYPE!
      type :: particle
      integer :: pidx,procidx,nbr_pidx,nbr_procidx
      real :: vp(3),xp(3),uf(3),xrhs(3),vrhs(3),Tp,Tprhs_s
      real :: Tprhs_L,Tf,radius,radrhs,qinf,qstar,dist
      real :: res,m_s,Os,rc
      real :: u_sub(3),sigm_s
      integer*8 :: mult
      type(particle), pointer :: prev,next
      end type particle

      type(particle), pointer :: part,first_particle

contains 

function mod_Magnus(T)
      implicit none

      !Take in T in Kelvin and return saturation vapor pressure using function of Alduchov and Eskridge, 1996
      real,intent(in) :: T  
      real :: mod_Magnus

      mod_Magnus = 610.94 *exp((17.6257*(T-273.15))/(243.04+(T-273.15)))
       
      
end function mod_Magnus

function crit_radius(m_s,Os,Tf)
      use pars
      use con_data
      implicit none

      integer :: i,maxidx(1)
      integer, parameter :: N=1000
      real :: m_s,Os,Tf
      real :: radval(N),SS(N)
      real :: radstart,radend,dr
      real :: crit_radius
      real :: firstterm,secterm

      radstart = -8
      radend = -3
      dr = (radstart-radend)/N
       
      do i=1,N

         radval(i) = 10**(radstart - (i-1)*dr)

         firstterm = 2*Mw*Gam/Ru/rhow/radval(i)/Tf
         secterm =Ion*Os*m_s*(Mw/Ms)/(rhow*pi2*2.0/3.0*radval(i)**3)

         SS(i) = exp( 2*Mw*Gam/Ru/rhow/radval(i)/Tf - Ion*Os*m_s*(Mw/Ms)/(rhow*pi2*2.0/3.0*radval(i)**3))

      end do

      maxidx = maxloc(SS)
      
      crit_radius = radval(maxidx(1))
      

end function crit_radius


end module particles
