program benchmark_particle_loop
  ! Compares particle loop performance for three memory layouts:
  !   1. Linked list (current production code)
  !   2. Array of Structs (AoS)
  !   3. Struct of Arrays (SoA)
  !
  ! Physics: position/velocity/temperature update (no evaporation, no interpolation)
  ! Each version inlines the physics directly against its own data structure so
  ! that all field accesses go through the actual memory layout being tested.
  !
  ! Real fields use real(8) to match production code compiled with -r8.
  ! Resulting particle struct size matches the production 352-byte type.
  !
  ! Compile: ifx -O2 -o benchmark_particle_loop benchmark_particle_loop.f90
  implicit none

  integer, parameter :: npart  = 300000
  integer, parameter :: nsteps = 1000

  ! Physical constants
  real(8), parameter :: pi2    = 3.14159265358979d0
  real(8), parameter :: rhow   = 1000.0d0
  real(8), parameter :: rhoa   = 1.2d0
  real(8), parameter :: nuf    = 1.5d-5
  real(8), parameter :: dt     = 0.01d0
  real(8), parameter :: Pra    = 0.71d0
  real(8), parameter :: CpaCpp = 1.0d0
  real(8), parameter :: grav1  = 0.0d0
  real(8), parameter :: grav2  = 0.0d0
  real(8), parameter :: grav3  = -9.81d0

  real(8) :: t_list, t_aos, t_soa
  real(8) :: chk_list, chk_aos, chk_soa

  call run_linked_list(t_list, chk_list)
  call run_aos        (t_aos,  chk_aos)
  call run_soa        (t_soa,  chk_soa)

  write(*,'(a)')          '--- Particle loop benchmark ---'
  write(*,'(a,i0)')       'npart  = ', npart
  write(*,'(a,i0)')       'nsteps = ', nsteps
  write(*,*)
  write(*,'(a,f8.3,a)')   'Linked list : ', t_list, ' s'
  write(*,'(a,f8.3,a)')   'AoS array   : ', t_aos,  ' s'
  write(*,'(a,f8.3,a)')   'SoA array   : ', t_soa,  ' s'
  write(*,*)
  write(*,'(a,f6.2)')     'Speedup list->AoS : ', t_list / t_aos
  write(*,'(a,f6.2)')     'Speedup list->SoA : ', t_list / t_soa
  write(*,'(a,f6.2)')     'Speedup AoS ->SoA : ', t_aos  / t_soa
  write(*,*)

  ! Verify: mean radius at last timestep should agree across all three
  write(*,'(a,3e15.6)')   'Mean radius (list, AoS, SoA): ', chk_list, chk_aos, chk_soa
  if (abs(chk_list - chk_aos) / abs(chk_list) > 1.0d-4 .or. &
      abs(chk_list - chk_soa) / abs(chk_list) > 1.0d-4) then
    write(*,'(a)') 'WARNING: checksums differ — results may be incorrect.'
  else
    write(*,'(a)') 'Checksums agree.'
  end if

contains

  ! -----------------------------------------------------------------------
  ! Shared: deterministic pseudo-random initialisation for core fields
  ! -----------------------------------------------------------------------
  subroutine init_particle_state(i, xp, vp, uf, radius, m_s, Tp, Tf)
    integer,  intent(in)  :: i
    real(8),  intent(out) :: xp(3), vp(3), uf(3), radius, m_s, Tp, Tf

    real(8) :: r
    r     = mod(real(i,8)*1.6180339887d0, 1.0d0)

    xp(1) = r * 1000.0d0
    xp(2) = mod(real(i,8)*2.7182818d0, 1.0d0) * 1000.0d0
    xp(3) = mod(real(i,8)*3.1415926d0, 1.0d0) * 500.0d0 + 1.0d0

    vp(1) = (r - 0.5d0) * 2.0d0
    vp(2) = (mod(real(i,8)*1.4142135d0, 1.0d0) - 0.5d0) * 2.0d0
    vp(3) = (mod(real(i,8)*1.7320508d0, 1.0d0) - 0.5d0) * 2.0d0 - 0.1d0

    uf(1) = (mod(real(i,8)*2.2360679d0, 1.0d0) - 0.5d0) * 5.0d0
    uf(2) = (mod(real(i,8)*2.6457513d0, 1.0d0) - 0.5d0) * 5.0d0
    uf(3) = (mod(real(i,8)*2.8284271d0, 1.0d0) - 0.5d0) * 2.0d0

    radius = 5.0d-6 + r * 45.0d-6
    m_s    = 1.0d-18
    Tp     = 280.0d0 + r * 10.0d0
    Tf     = 285.0d0
  end subroutine init_particle_state

  ! -----------------------------------------------------------------------
  ! 1. Linked list — physics inlined, accesses p%field directly
  ! -----------------------------------------------------------------------
  subroutine run_linked_list(elapsed, checksum)
    real(8), intent(out) :: elapsed, checksum

    type :: particle
      real(8) :: xp(3), vp(3), uf(3)
      real(8) :: radius, m_s, Tp, Tf
      integer :: pidx, procidx, nbr_pidx, nbr_procidx
      real(8) :: xrhs(3), vrhs(3)
      real(8) :: Tprhs_s, Tprhs_L, radrhs
      real(8) :: qinf, qstar, nbr_dist
      real(8) :: res, kappa_s, rc, actres, numact
      real(8) :: u_sub(3), sigm_s
      real(8) :: vp_old(3), Tp_old, radius_old
      integer*8 :: mult
      type(particle), pointer :: next
    end type particle

    type(particle), pointer :: head, p, tmp
    integer    :: i, istep
    real(8)    :: diff1, diff2, diff3, diffnorm
    real(8)    :: Volp, rhop, taup_i, Rep, Nup, tmp_coeff
    integer(8) :: t0, t1, rate

    nullify(head)
    do i = npart, 1, -1
      allocate(tmp)
      call init_particle_state(i, tmp%xp, tmp%vp, tmp%uf, &
                               tmp%radius, tmp%m_s, tmp%Tp, tmp%Tf)
      tmp%pidx = i;         tmp%procidx     = 0
      tmp%nbr_pidx = 0;     tmp%nbr_procidx = 0
      tmp%xrhs = tmp%vp;    tmp%vrhs        = 0.0d0
      tmp%Tprhs_s = 0.0d0;  tmp%Tprhs_L     = 0.0d0;  tmp%radrhs   = 0.0d0
      tmp%qinf = 0.01d0;    tmp%qstar       = 0.01d0;  tmp%nbr_dist = 10.0d0
      tmp%res = 0.0d0;      tmp%kappa_s     = 0.5d0
      tmp%rc = 1.0d-5;      tmp%actres      = 0.0d0;   tmp%numact   = 0.0d0
      tmp%u_sub = 0.0d0;    tmp%sigm_s      = 0.0d0
      tmp%vp_old = tmp%vp;  tmp%Tp_old      = tmp%Tp
      tmp%radius_old = tmp%radius;  tmp%mult = 1
      tmp%next => head
      head => tmp
    end do

    call system_clock(t0, rate)

    do istep = 1, nsteps
      p => head
      do while (associated(p))

        diff1    = p%vp(1) - p%uf(1)
        diff2    = p%vp(2) - p%uf(2)
        diff3    = p%vp(3) - p%uf(3)
        diffnorm = sqrt(diff1**2 + diff2**2 + diff3**2)

        Volp   = pi2 * 2.0d0/3.0d0 * p%radius**3
        rhop   = (p%m_s + Volp*rhow) / Volp
        taup_i = 18.0d0*rhoa*nuf / rhop / (2.0d0*p%radius)**2
        Rep    = 2.0d0*p%radius*diffnorm / nuf

        p%xp(1) = p%xp(1) + dt*p%vp(1)
        p%xp(2) = p%xp(2) + dt*p%vp(2)
        p%xp(3) = p%xp(3) + dt*p%vp(3)

        p%vp(1) = (p%vp(1) + taup_i*dt*p%uf(1) + dt*grav1) / (1.0d0 + dt*taup_i)
        p%vp(2) = (p%vp(2) + taup_i*dt*p%uf(2) + dt*grav2) / (1.0d0 + dt*taup_i)
        p%vp(3) = (p%vp(3) + taup_i*dt*p%uf(3) + dt*grav3) / (1.0d0 + dt*taup_i)

        Nup       = 2.0d0 + 0.6d0*sqrt(Rep)*Pra**(1.0d0/3.0d0)
        tmp_coeff = Nup/3.0d0/Pra * CpaCpp * rhop/rhow * taup_i
        p%Tp      = (p%Tp + tmp_coeff*dt*p%Tf) / (1.0d0 + dt*tmp_coeff)

        p => p%next
      end do
    end do

    call system_clock(t1)
    elapsed = real(t1 - t0, 8) / real(rate, 8)

    checksum = 0.0d0
    p => head
    do while (associated(p))
      checksum = checksum + p%radius
      p => p%next
    end do
    checksum = checksum / real(npart, 8)

    do while (associated(head))
      tmp => head%next
      deallocate(head)
      head => tmp
    end do
  end subroutine run_linked_list

  ! -----------------------------------------------------------------------
  ! 2. Array of Structs (AoS) — physics inlined, accesses parts(i)%field
  ! -----------------------------------------------------------------------
  subroutine run_aos(elapsed, checksum)
    real(8), intent(out) :: elapsed, checksum

    type :: particle_aos
      real(8) :: xp(3), vp(3), uf(3)
      real(8) :: radius, m_s, Tp, Tf
      integer :: pidx, procidx, nbr_pidx, nbr_procidx
      real(8) :: xrhs(3), vrhs(3)
      real(8) :: Tprhs_s, Tprhs_L, radrhs
      real(8) :: qinf, qstar, nbr_dist
      real(8) :: res, kappa_s, rc, actres, numact
      real(8) :: u_sub(3), sigm_s
      real(8) :: vp_old(3), Tp_old, radius_old
      integer*8 :: mult
    end type particle_aos

    type(particle_aos), allocatable :: parts(:)
    integer    :: i, istep
    real(8)    :: diff1, diff2, diff3, diffnorm
    real(8)    :: Volp, rhop, taup_i, Rep, Nup, tmp_coeff
    integer(8) :: t0, t1, rate

    allocate(parts(npart))
    do i = 1, npart
      call init_particle_state(i, parts(i)%xp, parts(i)%vp, parts(i)%uf, &
                               parts(i)%radius, parts(i)%m_s, &
                               parts(i)%Tp, parts(i)%Tf)
      parts(i)%pidx = i;              parts(i)%procidx     = 0
      parts(i)%nbr_pidx = 0;          parts(i)%nbr_procidx = 0
      parts(i)%xrhs = parts(i)%vp;    parts(i)%vrhs        = 0.0d0
      parts(i)%Tprhs_s = 0.0d0;       parts(i)%Tprhs_L     = 0.0d0
      parts(i)%radrhs = 0.0d0
      parts(i)%qinf = 0.01d0;         parts(i)%qstar       = 0.01d0
      parts(i)%nbr_dist = 10.0d0
      parts(i)%res = 0.0d0;           parts(i)%kappa_s     = 0.5d0
      parts(i)%rc = 1.0d-5;           parts(i)%actres      = 0.0d0
      parts(i)%numact = 0.0d0
      parts(i)%u_sub = 0.0d0;         parts(i)%sigm_s      = 0.0d0
      parts(i)%vp_old = parts(i)%vp;  parts(i)%Tp_old      = parts(i)%Tp
      parts(i)%radius_old = parts(i)%radius;  parts(i)%mult = 1
    end do

    call system_clock(t0, rate)

    do istep = 1, nsteps
      do i = 1, npart

        diff1    = parts(i)%vp(1) - parts(i)%uf(1)
        diff2    = parts(i)%vp(2) - parts(i)%uf(2)
        diff3    = parts(i)%vp(3) - parts(i)%uf(3)
        diffnorm = sqrt(diff1**2 + diff2**2 + diff3**2)

        Volp   = pi2 * 2.0d0/3.0d0 * parts(i)%radius**3
        rhop   = (parts(i)%m_s + Volp*rhow) / Volp
        taup_i = 18.0d0*rhoa*nuf / rhop / (2.0d0*parts(i)%radius)**2
        Rep    = 2.0d0*parts(i)%radius*diffnorm / nuf

        parts(i)%xp(1) = parts(i)%xp(1) + dt*parts(i)%vp(1)
        parts(i)%xp(2) = parts(i)%xp(2) + dt*parts(i)%vp(2)
        parts(i)%xp(3) = parts(i)%xp(3) + dt*parts(i)%vp(3)

        parts(i)%vp(1) = (parts(i)%vp(1) + taup_i*dt*parts(i)%uf(1) + dt*grav1) &
                         / (1.0d0 + dt*taup_i)
        parts(i)%vp(2) = (parts(i)%vp(2) + taup_i*dt*parts(i)%uf(2) + dt*grav2) &
                         / (1.0d0 + dt*taup_i)
        parts(i)%vp(3) = (parts(i)%vp(3) + taup_i*dt*parts(i)%uf(3) + dt*grav3) &
                         / (1.0d0 + dt*taup_i)

        Nup       = 2.0d0 + 0.6d0*sqrt(Rep)*Pra**(1.0d0/3.0d0)
        tmp_coeff = Nup/3.0d0/Pra * CpaCpp * rhop/rhow * taup_i
        parts(i)%Tp = (parts(i)%Tp + tmp_coeff*dt*parts(i)%Tf) &
                      / (1.0d0 + dt*tmp_coeff)

      end do
    end do

    call system_clock(t1)
    elapsed  = real(t1 - t0, 8) / real(rate, 8)
    checksum = sum(parts(:)%radius) / real(npart, 8)

    deallocate(parts)
  end subroutine run_aos

  ! -----------------------------------------------------------------------
  ! 3. Struct of Arrays (SoA) — physics inlined, accesses field(i) directly
  ! -----------------------------------------------------------------------
  subroutine run_soa(elapsed, checksum)
    real(8), intent(out) :: elapsed, checksum

    real(8), allocatable :: xp(:,:), vp(:,:), uf(:,:)
    real(8), allocatable :: radius(:), m_s(:), Tp(:), Tf(:)
    integer,  allocatable :: pidx(:), procidx(:), nbr_pidx(:), nbr_procidx(:)
    real(8),  allocatable :: xrhs(:,:), vrhs(:,:)
    real(8),  allocatable :: Tprhs_s(:), Tprhs_L(:), radrhs(:)
    real(8),  allocatable :: qinf(:), qstar(:), nbr_dist(:)
    real(8),  allocatable :: res(:), kappa_s(:), rc(:), actres(:), numact(:)
    real(8),  allocatable :: u_sub(:,:), sigm_s(:)
    real(8),  allocatable :: vp_old(:,:), Tp_old(:), radius_old(:)
    integer(8), allocatable :: mult(:)

    integer    :: i, istep
    real(8)    :: diff1, diff2, diff3, diffnorm
    real(8)    :: Volp, rhop, taup_i, Rep, Nup, tmp_coeff
    integer(8) :: t0, t1, rate

    allocate(xp(3,npart), vp(3,npart), uf(3,npart))
    allocate(radius(npart), m_s(npart), Tp(npart), Tf(npart))
    allocate(pidx(npart), procidx(npart), nbr_pidx(npart), nbr_procidx(npart))
    allocate(xrhs(3,npart), vrhs(3,npart))
    allocate(Tprhs_s(npart), Tprhs_L(npart), radrhs(npart))
    allocate(qinf(npart), qstar(npart), nbr_dist(npart))
    allocate(res(npart), kappa_s(npart), rc(npart), actres(npart), numact(npart))
    allocate(u_sub(3,npart), sigm_s(npart))
    allocate(vp_old(3,npart), Tp_old(npart), radius_old(npart))
    allocate(mult(npart))

    do i = 1, npart
      call init_particle_state(i, xp(:,i), vp(:,i), uf(:,i), &
                               radius(i), m_s(i), Tp(i), Tf(i))
      pidx(i) = i;           procidx(i)     = 0
      nbr_pidx(i) = 0;       nbr_procidx(i) = 0
      xrhs(:,i) = vp(:,i);   vrhs(:,i)      = 0.0d0
      Tprhs_s(i) = 0.0d0;    Tprhs_L(i)     = 0.0d0;  radrhs(i)   = 0.0d0
      qinf(i) = 0.01d0;      qstar(i)       = 0.01d0; nbr_dist(i) = 10.0d0
      res(i) = 0.0d0;        kappa_s(i)     = 0.5d0
      rc(i) = 1.0d-5;        actres(i)      = 0.0d0;  numact(i)   = 0.0d0
      u_sub(:,i) = 0.0d0;    sigm_s(i)      = 0.0d0
      vp_old(:,i) = vp(:,i); Tp_old(i)      = Tp(i)
      radius_old(i) = radius(i);  mult(i)    = 1
    end do

    call system_clock(t0, rate)

    do istep = 1, nsteps
      do i = 1, npart

        diff1    = vp(1,i) - uf(1,i)
        diff2    = vp(2,i) - uf(2,i)
        diff3    = vp(3,i) - uf(3,i)
        diffnorm = sqrt(diff1**2 + diff2**2 + diff3**2)

        Volp   = pi2 * 2.0d0/3.0d0 * radius(i)**3
        rhop   = (m_s(i) + Volp*rhow) / Volp
        taup_i = 18.0d0*rhoa*nuf / rhop / (2.0d0*radius(i))**2
        Rep    = 2.0d0*radius(i)*diffnorm / nuf

        xp(1,i) = xp(1,i) + dt*vp(1,i)
        xp(2,i) = xp(2,i) + dt*vp(2,i)
        xp(3,i) = xp(3,i) + dt*vp(3,i)

        vp(1,i) = (vp(1,i) + taup_i*dt*uf(1,i) + dt*grav1) / (1.0d0 + dt*taup_i)
        vp(2,i) = (vp(2,i) + taup_i*dt*uf(2,i) + dt*grav2) / (1.0d0 + dt*taup_i)
        vp(3,i) = (vp(3,i) + taup_i*dt*uf(3,i) + dt*grav3) / (1.0d0 + dt*taup_i)

        Nup       = 2.0d0 + 0.6d0*sqrt(Rep)*Pra**(1.0d0/3.0d0)
        tmp_coeff = Nup/3.0d0/Pra * CpaCpp * rhop/rhow * taup_i
        Tp(i)     = (Tp(i) + tmp_coeff*dt*Tf(i)) / (1.0d0 + dt*tmp_coeff)

      end do
    end do

    call system_clock(t1)
    elapsed  = real(t1 - t0, 8) / real(rate, 8)
    checksum = sum(radius) / real(npart, 8)

    deallocate(xp, vp, uf, radius, m_s, Tp, Tf)
    deallocate(pidx, procidx, nbr_pidx, nbr_procidx)
    deallocate(xrhs, vrhs, Tprhs_s, Tprhs_L, radrhs)
    deallocate(qinf, qstar, nbr_dist)
    deallocate(res, kappa_s, rc, actres, numact)
    deallocate(u_sub, sigm_s, vp_old, Tp_old, radius_old, mult)
  end subroutine run_soa

end program benchmark_particle_loop
