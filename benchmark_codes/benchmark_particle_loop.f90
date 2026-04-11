program benchmark_particle_loop
  ! Compares particle loop performance for three memory layouts:
  !   1. Linked list (current production code)
  !   2. Array of Structs (AoS)
  !   3. Struct of Arrays (SoA)
  !
  ! Physics: position/velocity/temperature update (no evaporation, no interpolation)
  ! Compile: ifx -O2 -o benchmark_particle_loop benchmark_particle_loop.f90
  implicit none

  integer, parameter :: npart  = 300000
  integer, parameter :: nsteps = 1000

  ! Physical constants
  real, parameter :: pi2    = 3.14159265358979
  real, parameter :: rhow   = 1000.0
  real, parameter :: rhoa   = 1.2
  real, parameter :: nuf    = 1.5e-5
  real, parameter :: dt     = 0.01
  real, parameter :: Pra    = 0.71
  real, parameter :: CpaCpp = 1.0
  real, parameter :: grav1  = 0.0
  real, parameter :: grav2  = 0.0
  real, parameter :: grav3  = -9.81

  real(8) :: t_list, t_aos, t_soa
  real    :: chk_list, chk_aos, chk_soa

  call run_linked_list(t_list, chk_list)
  call run_aos        (t_aos,  chk_aos)
  call run_soa        (t_soa,  chk_soa)

  write(*,'(a)')         '--- Particle loop benchmark ---'
  write(*,'(a,i0)')      'npart  = ', npart
  write(*,'(a,i0)')      'nsteps = ', nsteps
  write(*,*)
  write(*,'(a,f8.3,a)')  'Linked list : ', t_list, ' s'
  write(*,'(a,f8.3,a)')  'AoS array   : ', t_aos,  ' s'
  write(*,'(a,f8.3,a)')  'SoA array   : ', t_soa,  ' s'
  write(*,*)
  write(*,'(a,f6.2)')    'Speedup list->AoS : ', t_list / t_aos
  write(*,'(a,f6.2)')    'Speedup list->SoA : ', t_list / t_soa
  write(*,'(a,f6.2)')    'Speedup AoS ->SoA : ', t_aos  / t_soa
  write(*,*)

  ! Verify all three produce the same result (sum of all radii after nsteps)
  write(*,'(a,3e15.6)')  'Checksums (list, AoS, SoA): ', chk_list, chk_aos, chk_soa
  if (abs(chk_list - chk_aos) / abs(chk_list) > 1.0e-4 .or. &
      abs(chk_list - chk_soa) / abs(chk_list) > 1.0e-4) then
    write(*,'(a)') 'WARNING: checksums differ — results may be incorrect.'
  else
    write(*,'(a)') 'Checksums agree.'
  end if

contains

  ! -----------------------------------------------------------------------
  ! Shared: compute one particle update step given its state variables.
  ! Returns updated xp, vp, Tp (all intent(inout)).
  ! -----------------------------------------------------------------------
  subroutine update_particle(xp, vp, uf, radius, m_s, Tp, Tf, &
                              Volp, rhop, taup_i)
    real, intent(inout) :: xp(3), vp(3), Tp
    real, intent(in)    :: uf(3), radius, m_s, Tf
    real, intent(out)   :: Volp, rhop, taup_i

    real :: diff(3), diffnorm, Rep, Nup, tmp_coeff

    diff     = vp - uf
    diffnorm = sqrt(diff(1)**2 + diff(2)**2 + diff(3)**2)

    Volp   = pi2 * 2.0/3.0 * radius**3
    rhop   = (m_s + Volp*rhow) / Volp
    taup_i = 18.0*rhoa*nuf / rhop / (2.0*radius)**2
    Rep    = 2.0*radius*diffnorm / nuf

    ! Position and velocity (backward Euler)
    xp(1:3) = xp(1:3) + dt*vp(1:3)
    vp(1) = (vp(1) + taup_i*dt*uf(1) + dt*grav1) / (1.0 + dt*taup_i)
    vp(2) = (vp(2) + taup_i*dt*uf(2) + dt*grav2) / (1.0 + dt*taup_i)
    vp(3) = (vp(3) + taup_i*dt*uf(3) + dt*grav3) / (1.0 + dt*taup_i)

    ! Temperature (Ranz-Marshall Nusselt, backward Euler)
    Nup       = 2.0 + 0.6*sqrt(Rep)*Pra**(1.0/3.0)
    tmp_coeff = Nup/3.0/Pra * CpaCpp * rhop/rhow * taup_i
    Tp        = (Tp + tmp_coeff*dt*Tf) / (1.0 + dt*tmp_coeff)
  end subroutine update_particle

  ! -----------------------------------------------------------------------
  ! 1. Linked list
  ! -----------------------------------------------------------------------
  subroutine run_linked_list(elapsed, checksum)
    real(8), intent(out) :: elapsed
    real,    intent(out) :: checksum

    type :: particle
      real    :: xp(3), vp(3), uf(3)
      real    :: radius, m_s, Tp, Tf
      type(particle), pointer :: next
    end type particle

    type(particle), pointer :: head, p, tmp
    integer :: i, istep
    real    :: Volp, rhop, taup_i
    integer(8) :: t0, t1, rate

    ! Build linked list
    nullify(head)
    do i = npart, 1, -1
      allocate(tmp)
      call init_particle_state(i, tmp%xp, tmp%vp, tmp%uf, &
                               tmp%radius, tmp%m_s, tmp%Tp, tmp%Tf)
      tmp%next => head
      head     => tmp
    end do

    call system_clock(t0, rate)

    do istep = 1, nsteps
      p => head
      do while (associated(p))
        call update_particle(p%xp, p%vp, p%uf, p%radius, p%m_s, &
                             p%Tp, p%Tf, Volp, rhop, taup_i)
        p => p%next
      end do
    end do

    call system_clock(t1)
    elapsed = real(t1 - t0, 8) / real(rate, 8)

    ! Checksum
    checksum = 0.0
    p => head
    do while (associated(p))
      checksum = checksum + p%radius
      p => p%next
    end do

    ! Cleanup
    do while (associated(head))
      tmp => head%next
      deallocate(head)
      head => tmp
    end do
  end subroutine run_linked_list

  ! -----------------------------------------------------------------------
  ! 2. Array of Structs (AoS)
  ! -----------------------------------------------------------------------
  subroutine run_aos(elapsed, checksum)
    real(8), intent(out) :: elapsed
    real,    intent(out) :: checksum

    type :: particle_aos
      real :: xp(3), vp(3), uf(3)
      real :: radius, m_s, Tp, Tf
    end type particle_aos

    type(particle_aos), allocatable :: parts(:)
    integer   :: i, istep
    real      :: Volp, rhop, taup_i
    integer(8) :: t0, t1, rate

    allocate(parts(npart))
    do i = 1, npart
      call init_particle_state(i, parts(i)%xp, parts(i)%vp, parts(i)%uf, &
                               parts(i)%radius, parts(i)%m_s, &
                               parts(i)%Tp, parts(i)%Tf)
    end do

    call system_clock(t0, rate)

    do istep = 1, nsteps
      do i = 1, npart
        call update_particle(parts(i)%xp, parts(i)%vp, parts(i)%uf, &
                             parts(i)%radius, parts(i)%m_s, &
                             parts(i)%Tp, parts(i)%Tf, Volp, rhop, taup_i)
      end do
    end do

    call system_clock(t1)
    elapsed  = real(t1 - t0, 8) / real(rate, 8)
    checksum = sum(parts(:)%radius)

    deallocate(parts)
  end subroutine run_aos

  ! -----------------------------------------------------------------------
  ! 3. Struct of Arrays (SoA)
  ! -----------------------------------------------------------------------
  subroutine run_soa(elapsed, checksum)
    real(8), intent(out) :: elapsed
    real,    intent(out) :: checksum

    real, allocatable :: xp(:,:), vp(:,:), uf(:,:)
    real, allocatable :: radius(:), m_s(:), Tp(:), Tf(:)

    integer   :: i, istep
    real      :: xp_i(3), vp_i(3), uf_i(3), radius_i, m_s_i, Tp_i, Tf_i
    real      :: Volp, rhop, taup_i
    integer(8) :: t0, t1, rate

    allocate(xp(3,npart), vp(3,npart), uf(3,npart))
    allocate(radius(npart), m_s(npart), Tp(npart), Tf(npart))

    do i = 1, npart
      call init_particle_state(i, xp(:,i), vp(:,i), uf(:,i), &
                               radius(i), m_s(i), Tp(i), Tf(i))
    end do

    call system_clock(t0, rate)

    do istep = 1, nsteps
      do i = 1, npart
        xp_i = xp(:,i); vp_i = vp(:,i); uf_i = uf(:,i)
        radius_i = radius(i); m_s_i = m_s(i); Tp_i = Tp(i); Tf_i = Tf(i)
        call update_particle(xp_i, vp_i, uf_i, radius_i, m_s_i, &
                             Tp_i, Tf_i, Volp, rhop, taup_i)
        xp(:,i) = xp_i; vp(:,i) = vp_i; Tp(i) = Tp_i
      end do
    end do

    call system_clock(t1)
    elapsed  = real(t1 - t0, 8) / real(rate, 8)
    checksum = sum(radius)

    deallocate(xp, vp, uf, radius, m_s, Tp, Tf)
  end subroutine run_soa

  ! -----------------------------------------------------------------------
  ! Shared initializer: deterministic pseudo-random state for particle i
  ! -----------------------------------------------------------------------
  subroutine init_particle_state(i, xp, vp, uf, radius, m_s, Tp, Tf)
    integer, intent(in)  :: i
    real,    intent(out) :: xp(3), vp(3), uf(3), radius, m_s, Tp, Tf

    real :: r
    r       = mod(real(i)*1.6180339887, 1.0)   ! quasi-random in [0,1)

    xp(1)   = r * 1000.0
    xp(2)   = mod(real(i)*2.7182818, 1.0) * 1000.0
    xp(3)   = mod(real(i)*3.1415926, 1.0) * 500.0 + 1.0

    vp(1)   = (r - 0.5) * 2.0
    vp(2)   = (mod(real(i)*1.4142135, 1.0) - 0.5) * 2.0
    vp(3)   = (mod(real(i)*1.7320508, 1.0) - 0.5) * 2.0 - 0.1

    uf(1)   = (mod(real(i)*2.2360679, 1.0) - 0.5) * 5.0
    uf(2)   = (mod(real(i)*2.6457513, 1.0) - 0.5) * 5.0
    uf(3)   = (mod(real(i)*2.8284271, 1.0) - 0.5) * 2.0

    radius  = 5.0e-6 + r * 45.0e-6      ! 5–50 µm
    m_s     = 1.0e-18                    ! fixed solute mass
    Tp      = 280.0 + r * 10.0          ! 280–290 K
    Tf      = 285.0
  end subroutine init_particle_state

end program benchmark_particle_loop
