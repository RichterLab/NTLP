! Originally developed for use in the Diamessis group at Cornell by Greg
! Thomsen.  Lightly adapted for use in NTLP, mostly to remove the explicit
! use of a 64-bit floating point type which is set via compiler flags.

module stopwatch
    ! Timing module.  Contains routines to compute time spent in a section of
    ! code by marking a reference point (via stopwatch_tick()) and computing the
    ! elapsed time (via stopwatch_get_elapsed()).  Prior to use, the module must
    ! be initialized (via stopwatch_initialize()).

    use, intrinsic :: iso_c_binding, only:    C_DOUBLE

    use precision, only:                      int_64

    implicit none
    save

    ! Rate at which the clock advances on this system.  This is necessary to
    ! convert clock ticks into actual time.
    !
    ! NOTE: We explicitly use 64-bit integers to ensure the highest precision
    !       timer available is used by SYSTEM_CLOCK().
    !
    ! NOTE: We explicitly initialize this to 0 so we get NaN's if the module
    !       isn't initialized properly (via initialize_stopwatches()).  That
    !       should be noticeable to someone without adding overhead to check
    !       for initialization.
    integer(kind=int_64) :: clock_rate = 0_int_64

    real(kind=C_DOUBLE)  :: timeGlobal, timeLocal

    integer(kind=int_64) :: clock_start

contains

    subroutine stopwatch_initialize()
        ! Initializes the Stopwatches module so that its routines may be called.
        ! This *must* be called before any other subroutines or functions from
        ! this module otherwise invalid results will be computed.

        implicit none

        ! Determine how finely we can time things, though only do it once
        ! should someone accidentally call this routine multiple times.
        if (clock_rate == 0) then
            call SYSTEM_CLOCK( COUNT_RATE=clock_rate )
        end if

    end subroutine stopwatch_initialize

    subroutine stopwatch_tick( watch )
        ! Ticks the stopwatch and records the current time.  This can be used as
        ! a reference point in time so as to compute elapsed time via
        ! stopwatch_get_elapsed().  The Stopwatches module must be initialized
        ! via initialize_stopwatches() before this function may be called.

        use precision, only:    int_64

        implicit none

        integer(kind=int_64), intent(inout) :: watch

        call SYSTEM_CLOCK( COUNT=watch )

    end subroutine stopwatch_tick

    function stopwatch_get_elapsed( before ) result( elapsed )
        ! Computes the elapsed time since the supplied stopwatch tick.  The
        ! Stopwatches module must be initialized via initialize_stopwatches()
        ! before this function may be called.

        use precision, only:    int_64

        implicit none

        integer(kind=int_64) :: before
        real                 :: elapsed

        integer(kind=int_64) :: now

        ! Get the current clock value and compute the elapsed time based on the
        ! previously recorded rate.
        call SYSTEM_CLOCK( COUNT=now )

        elapsed = real( (now - before) ) / clock_rate

    end function stopwatch_get_elapsed

    subroutine seconds_since_epoch( timestamp )
        ! Computes the fractional seconds since midnight on January 1st, 1970 in
        ! UTC.

        implicit none

        ! Fractional seconds since midnight on January 1st, 1970 in UTC.
        real, intent(out)       :: timestamp

        ! Used to compute timestamp.
        integer, dimension(1:8) :: current_time_array
        integer                 :: todays_day
        integer                 :: epoch_day

        ! Get the current time, in component form.
        call DATE_AND_TIME( VALUES=current_time_array )

        ! Compute the number of days from some time in the distant past to the
        ! Unix Epoch and today.
        call julian_day( 1970, 1, 1, epoch_day )
        call julian_day( current_time_array(1), current_time_array(2), &
                         current_time_array(3), todays_day )

        ! Compute the current time, in fractional epoch seconds, in UTC.  Note
        ! that the 4th component of the current time is an offset to be added
        ! to convert from UTC to the local time, so we have to subtract it to
        ! go back to UTC.
        timestamp = (todays_day - epoch_day) * 86400.0 + &
                    current_time_array(5) * 3600.0 + &
                    (current_time_array(6) - current_time_array(4)) * 60.0 + &
                    current_time_array(7) + &
                    current_time_array(8) / 1000.0
        return

    end subroutine seconds_since_epoch

    subroutine julian_day( year, month, day, jday )
        ! Computes the Julian day for the supplyed year, month, and day of
        ! month.  The computed day is relative to some date in the distant past,
        ! such that January 1st, 1970 is Julian day 2,440,588.
        !
        ! Adapted from the equation given in Communications of the ACM, Volume
        ! 11, Issue 10 of October 1968, page 657.  See "Letters to the editor: a
        ! machine algorithm for processing calendar dates" by Henry F. Fliegel
        ! and Thomas C. Van Flandern.

        implicit none

        integer, intent(in)  :: year
        integer, intent(in)  :: month
        integer, intent(in)  :: day
        integer, intent(out) :: jday

        ! Magic!
        jday = day - 32075 + 1461 * (year + 4800 + (month - 14) / 12) / 4 +  &
               367 * (month - 2 - ((month - 14) / 12) * 12) / 12 - &
               3 * ((year + 4900 + (month - 14) / 12) / 100) / 4

        return

    end subroutine julian_day

end module stopwatch

module measure

    ! simple interface to measure the wall time of individual phases (sections)
    ! of code.  multiple phases can be active at once, so long as the sequence
    ! of phases started is the reverse of the sequence of phases ended.
    ! measurements are made using the stopwatch module and have the same
    ! precision that it provides.  additionally, measurements are OpenMP-aware
    ! so individual OpenMP threads may have their execution measured.
    !
    ! this module's interface is designed so that the measurement process
    ! is roughly:
    !
    !   1. initialize the module's state
    !   2. create the phases to measure
    !   3. start and end phases throughout the code
    !   4. query the measurements collected, or dump them to a file unit
    !   5. deallocate the resources associated with the module
    !
    ! NOTE: this module does not provide any phase identifiers as it is intended
    !       to serve as the foundation for application-specific measurement
    !       modules.
    !
    ! while the module's measurements are OpenMP-aware, there are many
    ! gotchas that make measuring OpenMP code more difficult than one would
    ! expect.  specifically, module users should be aware of the following:
    !
    !   * initialize_measurements() *must* be called from outside of all
    !     OpenMP parallel regions.  this is due to the vagaries of how OpenMP
    !     exposes the number of threads and the fact that the specification
    !     does not attempt to make this possible given all the ways they can
    !     be spawned.  see below for more details.
    !
    !     all of the module's public routines may be called from within parallel
    !     regions, provided the following restrictions are adhered to.
    !
    !   * it is assumed that OpenMP will be used to divide the same type of
    !     work amongst its threads in a simple "distribute this loop's
    !     iterations across available threads" manner and not to distribute
    !     different tasks across threads.  this design decision was made based
    !     on how hybrid MPI/OpenMPI typically work.
    !
    !     this assumption remains unless significant additional complexity is
    !     added to the module to track task-based phases, or by forcing users to
    !     be very careful about which phases are used by which threads, neither
    !     of which were needed for the initial implementation.
    !
    !   * phases are blind to any OpenMP parallel regions that are started
    !     and completed within them.  the measurements made for such phases
    !     will have the correct elapsed time, but will not attribute run-time
    !     to each of the threads that participated in the parallel region.
    !
    !     start and end phases within an OpenMP parallel region to measure
    !     their execution across threads.
    !
    !   * phases must start and complete in the same execution type of the
    !     program.  that is, if a phase starts in a sequential region then it
    !     must end in a sequential region, just as a phase starting in a
    !     parallel region must end in a parallel region.  phases that do
    !     not respect this requirement will generate invalid measurements
    !     across the known threads.
    !
    !   * this implementation does *NOT* support nested parallel regions!
    !     the reason for this is that uniquely identifying how many total
    !     threads can run (so as to size our measurements array) and correctly
    !     mapping a nested thread identifier, which is not globally unique, to
    !     something that is globally unique, is non-trivial.  we check for
    !     nested parallelism when we initialize things and bail there so as to
    !     avoid silently producing invalid measurements.
    !
    !     as a result of not supporting nested regions, the following points
    !     are worth noting:
    !
    !       - getting the current thread identifier requires care!  you cannot
    !         just use omp_get_thread_num() as that will tell you the identifier
    !         in the current region, which will normally be 1 if you are inside
    !         a nested parallel region.  in this case, one needs to use max( 1,
    !         omp_get_ancestor_thread_num( 1 ) + 1 ) to get a globally unique
    !         thread identifier.
    !
    !       - initialize_measurements() *must* be called from outside all
    !         parallel regions, otherwise it will stop execution.
    !         initialization requires assessing the maximum number of threads
    !         that can run concurrently, which is hard when we're inside a
    !         parallel region.
    !
    !     if for some reason this functionality is required, then the following
    !     needs to be done:
    !
    !       - empirically determine how many threads can run concurrently
    !         with multiple nested regions.  this could recursively call
    !         a function that returns the largest thread number across a
    !         parallel region and take the product across each nested
    !         invocation.  recursion would be terminated when only one
    !         thread was encountered, or possibly when a sequence of one
    !         thread is seen.  it does not appear that there is an OpenMP
    !         routine that will tell us the current product of the number
    !         of threads at each active level.
    !
    !         each instance of omp_get_max_threads() would need to be
    !         replaced.
    !
    !       - start_phase() and end_phase() need to be updated to correctly
    !         identify the globally unique thread identifier currently
    !         executing.  this is necessary to avoid a race condition when
    !         accessing entries in the measurements array.  if this does
    !         not work, threads will stop execution when they detect that
    !         the phase ending does not match the currently active phase.
    !         this is because two threads enter a region (and one is
    !         lost/never recorded in measurements) and another thread
    !         leaves a region.
    !
    !         this likelihood that this race is encountered is maximized
    !         by starting/ending a phase in an OpenMP-parallelized loop.
    !
    !         getting a globally unique thread identifier will require
    !         tracking the maximum number of threads at each nesting level
    !         and then mapping the current thread by the product of each of
    !         current active nested region's maximum number of threads.
    !         this is *hard* to do in a general case (if the application
    !         can change the maximum number of threads at any time) but
    !         seems feasible if they're static throughout the lifetime of
    !         the measurements.
    !

    use Int64Map_h, only: Int64Map
    use IntStack_h, only: IntStack
    use RealMap_h, only:  RealMap

    implicit none

    ! everything in the module is private by default.  the external interface is
    ! marked as appropriate.
    private

    integer, public :: REDUCTION_TYPE_MAX = 0
    integer, public :: REDUCTION_TYPE_SUM = 1

    ! internal type used to track phase measurements across OpenMP threads.
    type :: Measurement

        ! stack of phase identifiers representing the current sequence of active
        ! phases.  the phase identifiers stored are indices in the phase states
        ! map.
        type(IntStack) :: active_phases

        ! maps from phase name to the most recent start time and its cumulative
        ! duration for all previously completed phases.  start times are in
        ! seconds since the Epoch.
        !
        ! NOTE: we must use 64-bit integers as this is matches the interface
        !       provided by the stopwatch module.
        !
        type(Int64Map) :: phase_start
        type(RealMap)  :: phase_duration

    end type Measurement

    ! one per OpenMP thread to track its journey through the application's
    ! phases.  private since this is the module's internal state.
    !
    ! NOTE: each Measurement object has its own phase stack so individual
    !       threads can have different call trees.  while this isn't
    !       anticipated, it's really not that much more complicated to be
    !       flexible up front.
    !
    type(Measurement), allocatable :: measurements(:)

    ! the number of entries in measurements.  this will be the number of
    ! threads executing when the module is initialized.
    !
    ! NOTE: this is stored as a module variable so as to provide more
    !       flexibility to non-initialization routines.  since nested parallel
    !       regions report thread numbers and counts relative to the size of
    !       that region, we set this once at initialization so all other
    !       module routines can be called without restriction.
    !
    integer :: number_measurements

    ! declare the external interface.
    public :: &
         create_phase, &
         dump_measurements, &
         end_phase, &
         get_duration, &
         get_durations, &
         initialize_measurements, &
         is_initialized, &
         print_duration, &
         reset_duration, &
         shutdown_measurements, &
         start_phase

contains

    ! initializes the module's state so that measurements can be performed.
    !
    ! NOTE: this must be called before any of the other routines provided
    !       by this module.  it also should only be called at most once
    !       before shutdown_measurements() is called.  once the module has
    !       been shutdown, it may be initialized again.
    !
    subroutine initialize_measurements()

        !$ use omp_lib, only:    omp_get_level, omp_get_max_threads, omp_get_max_active_levels, &
        !$                       omp_get_nested

        use stopwatch, only:    stopwatch_initialize

        implicit none

        integer :: thread_number
        logical :: nesting_enabled_flag
        integer :: parallel_region_level

        ! let the user know their code is broken if they attempt to initialize
        ! multiple times.
        if( is_initialized() ) then
            deallocate( measurements )

            write( 6, "(A)" ) "WARNING: The measurement module has been initialized multiple times.  Fix this."
            return
        end if

        ! detect when nesting is enabled and blow up.
        nesting_enabled_flag = .false.
        !$ nesting_enabled_flag = (omp_get_max_active_levels() > 1)
        if( nesting_enabled_flag .eqv. .true. ) then

            ! double check whether the deprecated interface confirms that
            ! nesting is enabled.  some compilers incorrectly report a
            ! maximum level larger than one, but .false. when using
            ! omp_get_nested().
            !
            ! NOTE: this is a workaround to support both GNU gfortran 7.5, which
            !       does not correctly report the maximum nesting level, and
            !       newer versions of Intel's oneAPI, which complain that
            !       omp_get_nested() is deprecated in newer versions of OpenMP
            !       (presumably in 5.0, info message #269).  sadly, we can't use
            !       the _OPENMP preprocessor symbol to pick a single conditional
            !       as gfortran 7.x through at least 13.2 all report supporting
            !       OpenMP 4.5 (reported as 201511).
            !
            !       falling back to the deprecated interface ensures correct
            !       behavior in legacy compilers while avoiding the warning
            !       from Intel's compilers, without having to include
            !       preprocessor checks against the compiler and its version.
            !
            nesting_enabled_flag    = .false.
            !$ nesting_enabled_flag = omp_get_nested()
            if( nesting_enabled_flag .eqv. .true. ) then
                write( 6, "(A)" ) "Nested OpenMP parallel regions are enabled!  The measurement " // &
                     "module does not support this!"
                stop
            end if
        end if

        ! refuse to initialize when we're in a parallel region.  we have no easy
        ! way to identify the maximum number of threads in the application if we
        ! are.
        parallel_region_level = 0
        !$ parallel_region_level = omp_get_level()
        if( parallel_region_level .ne. 0 ) then
            write( 6, "(A)" ) "The measurement module must be initialized outside of all " // &
                        "parallel regions!"
            stop
        end if

        ! we track time spent in individual phases on a per-thread basis.
        number_measurements = 1
        !$ number_measurements = omp_get_max_threads()

        allocate( measurements(number_measurements) )

        do thread_number = 1, number_measurements
            call measurements(thread_number)%active_phases%initialize()
        end do

        ! setup our fine-grained timer.  this is used to track the start and end
        ! of each phase we're interested in.
        call stopwatch_initialize()

    end subroutine initialize_measurements

    ! returns the status of the module.  if .false. is returned, then a call to
    ! initialize_measurements() must be made before using any of the module's
    ! interfaces.
    function is_initialized() result( initialized_flag )

        implicit none

        logical :: initialized_flag

        initialized_flag = allocated( measurements )

    end function is_initialized

    ! shuts down the module and releases the measurement resources allocated
    ! when it was initialized.  additional calls to the module's routines,
    ! except another call to initialize_measurements() to re-initialize the
    ! module, are forbidden.
    subroutine shutdown_measurements()

        implicit none

        ! release our measurement resources.  take care to not explode if this
        ! is accidentally called twice without an intervening call to
        ! initialize_measurements().
        if( is_initialized() ) then
            deallocate( measurements )
        end if

    end subroutine shutdown_measurements

    ! writes a single phase's duration to the supplied unit.  care is taken to
    ! format the duration and the percentage relative to the parent phase's
    ! duration.  output is of the form:
    !
    !    <description><duration>s (<duration>/<total_duration>%)
    !
    subroutine print_duration( file_unit, description, duration, total_duration )

        use precision, only:    dp_equals

        implicit none

        integer, intent(in)          :: file_unit
        character(len=*), intent(in) :: description
        real, intent(in)             :: duration
        real, intent(in)             :: total_duration

        real                         :: duration_percentage

        ! format descriptor of the form:
        !
        !     <description><duration>s (<duration>/<total_duration>%)
        !
        ! NOTE: we format our durations with F9.1 so we maintain aligned output
        !       for reasonably long durations.  30 days in seconds is 2.59M
        !       seconds, which requires 7 digits for the integral seconds, one
        !       for the decimal point, and one more for the 10th of a second
        !       remainder.  we use the 'F' edit specifier so we get a leading
        !       zero for sub-second durations.
        !
        !       we format our percentages with F5.1 so we maintain aligned
        !       outputs for each section.  we don't report fine-grained
        !       percentages so we only need three digits for the integral
        !       precision (in the range of [0, 100]), one for the decimal point,
        !       and one for the 10th of a percentage point fractional
        !       percentage.  we use the 'F' edit specifier so we get a leading
        !       zero for sub-1% percentages.
        !
        character(len=*), parameter :: descriptor_timing = "(A,F9.1,A,F9.1,A,F5.1,A)"

        ! hard-coded strings for joining the durations and the computed
        ! percentage.
        character(len=*), parameter :: of                = "s of "
        character(len=*), parameter :: open_parenthesis  = "s ("
        character(len=*), parameter :: close_parenthesis = "%)"

        ! guard against a division by zero.
        if( dp_equals( total_duration, 0.0 ) ) then
            duration_percentage = 0.0
        else
            duration_percentage = duration / total_duration * 100.0
        end if

        write( file_unit, descriptor_timing ) description, &
             duration, of, total_duration, open_parenthesis, &
             duration_percentage, close_parenthesis

    end subroutine print_duration

    ! dumps the raw measurements structure to the supplied unit for debugging
    ! purposes.  measurements are dumped as time spent in each section on a
    ! per-thread basis.
    !
    ! NOTE: this is intended for debugging purposes and should not be used
    !       in "production" code.  use report_measurements() instead when
    !       timings need to be conveyed to end users.
    !
    subroutine dump_measurements( file_unit )

        implicit none

        integer, intent(in)           :: file_unit

        integer                       :: number_phases
        integer                       :: phase_number
        integer                       :: thread_number

        real, allocatable             :: phase_durations(:)
        character(len=:), allocatable :: phase_name

        if( .not. is_initialized() ) then
            write( 6, "(A)" ) "The measurement module has not been initialized before " // &
                 "calling dump_measurements()."
            stop
        end if

        ! we report each phase's duration for each of the threads that executed
        ! it.
        allocate( phase_durations(number_measurements) )

        number_phases = measurements(1)%phase_start%getSize()

        ! remind the user how many threads this execution has.
        write( file_unit, "(A,I0)" ) "Number of threads: ", number_measurements

        ! collect the phase names, and accumulated per-thread durations, and
        ! write them to the supplied unit.
        !
        ! XXX: this is a hack!  assumes that phase identifiers start at 1 and go
        !      to N.
        do phase_number = 1, number_phases
            do thread_number = 1, number_measurements
                phase_durations(thread_number) = measurements(thread_number)%phase_duration%getValue( phase_number )
            end do

            phase_name = measurements(1)%phase_duration%getKey( phase_number )

            !
            ! NOTE: we specify F9.5 to ensure we have a leading zero on our
            !       durations and also provide some level of detail that is
            !       useful for debugging.
            !
            write( file_unit, "(A,A,A,*(F9.5:', '))" ) "'", trim( phase_name ), "' = ", phase_durations
        end do

    end subroutine dump_measurements

    ! registers a phase name and returns the corresponding identifier.
    ! registration is done for all threads' internal state.  the returned phase
    ! identifier is only valid until measurements are reset via
    ! shutdown_measurements().
    function create_phase( phase_name ) result( phase_id )

        !
        ! NOTE: this relies on the assumption that unique ArrayMap objects that
        !       are updated in the same sequence will be identical.
        !

        use precision, only:            int_64

        implicit none

        character(len=*), intent(in) :: phase_name
        integer                      :: phase_id

        integer                      :: thread_number

        integer, allocatable         :: phase_ids(:)

        if( .not. is_initialized() ) then
            write( 6, "(A)" ) "The measurement module has not been initialized before " // &
                 "calling create_phase()."
            stop
        end if

        allocate( phase_ids(number_measurements) )

        ! ensure we don't register a duplicate phase name.
        if( measurements(1)%phase_start%find( phase_name ) .ne. -1 ) then
            write( 6, "(A)" ) "'" // trim( phase_name ) // "' has already been registered!"
            stop
        end if

        ! initialize each thread's state.
        do thread_number = 1, number_measurements

            !
            ! NOTE: we should not initialize the phase stack since phases may be
            !       created on the fly.
            !

            ! add this phase into each of the underlying maps.
            call measurements(thread_number)%phase_start%addEntry( phase_name, int( 0, kind=int_64 ) )
            call measurements(thread_number)%phase_duration%addEntry( phase_name, 0.0 )

            ! ensure that each of the maps used the same identifier.  this
            ! module assumes they will be the same to simplify its logic
            ! and we need to fix the underlying container if that is not the
            ! case.
            phase_ids(thread_number) = measurements(thread_number)%phase_start%find( phase_name )
            if( phase_ids(thread_number) .ne. measurements(thread_number)%phase_duration%find( phase_name ) ) then
                write( 6, "(A)" ) "Phase identifiers are inconsistent between maps!  " // &
                     "Fix the IntMap and RealMap implementations!"
                stop
            end if
        end do

        ! ensure that each thread used the same identifier for this phase.
        if( .not. all( phase_ids .eq. phase_ids(1) ) ) then
            write( 6, "(A)" ) "Phase identifiers are inconsistent between threads!  Fix the ArrayMap implementation."
            stop
        end if

        ! all of the phase identifiers are consistent, so we can return on that
        ! works for all threads.
        phase_id = phase_ids(1)

    end function create_phase

    ! computes the maximum duration spent by any one thread for a given phase.
    !
    ! this is effectively the wall-time elapsed while executing a phase.  since
    ! measurements are made on a per-thread basis, the caller is responsible
    ! for specifying the type of reduction to perform across threads.  when
    ! omitted, the maximum duration is returned.
    !
    ! NOTE: this is assumes that threads divide a region of work instead of
    !       doing disparate work.  if that assumption is false, and threads
    !       measure disparate phases, then this returns an underestimate of
    !       the time spent as we don't have enough information to compute the
    !       effective wall-time.
    !
    function get_duration( phase_id, reduction_type ) result( phase_duration )

        implicit none

        integer, intent(in)           :: phase_id
        integer, intent(in), optional :: reduction_type

        integer                       :: actual_reduction_type
        integer                       :: thread_number

        real                          :: phase_duration

        if( .not. is_initialized() ) then
            write( 6, "(A)" ) "The measurement module has not been initialized before " // &
                 "calling get_duration()."
            stop
        end if

        call validate_phase_id( phase_id, "get_duration()" )

        if( present( reduction_type ) ) then
            actual_reduction_type = reduction_type
        else
            ! we default to taking the maximum duration across threads as this
            ! corresponds to the most common case of measuring parallel sections
            ! where work is divided evenly across threads and all measurements
            ! are roughly the same.  as a result, this simplifies the call sites
            ! since they can omit this argument from the function call.
            actual_reduction_type = REDUCTION_TYPE_MAX
        end if

        ! take the maximum wall-time across each of the threads so we report
        ! the amount of time each section took, from the user's perspective.
        !
        ! NOTE: this assumes that threads work on the same thing, rather than
        !       disparate things.  if this is not the case, then the reported
        !       duration is smaller than the actual time observed.
        !
        !       for example, a parallel loop partitioned amongst threads should
        !       take the maximum (the slowest) across threads for its duration,
        !       while distinct tasks (i.e. using OpenMP sections) will not
        !       reflect the observed duration if the tasks are not concurrent.
        !
        ! NOTE: there is the potential for a data race here!  we do not lock
        !       access to the measurements to minimize access overhead.  we
        !       do not expect to query phase durations while measurements are
        !       being performed.
        !
        phase_duration = 0.0
        do thread_number = 1, size( measurements )
            if( actual_reduction_type == REDUCTION_TYPE_MAX ) then
                phase_duration = max( phase_duration, &
                                      measurements(thread_number)%phase_duration%getValue( phase_id ) )
            else if( actual_reduction_type == REDUCTION_TYPE_SUM ) then
                phase_duration = (phase_duration + &
                                  measurements(thread_number)%phase_duration%getValue( phase_id ))
            end if
        end do

    end function get_duration

    ! returns the current durations measured by each thread for a given phase.
    !
    ! NOTE: this is not thread safe!  this may return out of date measurements
    !       if called when other threads are making measurements.
    !
    function get_durations( phase_id ) result( phase_durations )

        implicit none

        integer, intent(in) :: phase_id

        integer             :: thread_number
        real, allocatable   :: phase_durations(:)

        if( .not. is_initialized() ) then
            write( 6, "(A)" ) "The measurement module has not been initialized before " // &
                 "calling get_durations()."
            stop
        end if

        allocate( phase_durations(number_measurements) )

        !
        ! NOTE: there is the potential for a data race here!  we do not lock
        !       access to the measurements to minimize access overhead.  we
        !       do not expect measurements to be taken in one thread while
        !       another thread queries the durations.
        !
        do thread_number = 1, number_measurements
            phase_durations(thread_number) = measurements(thread_number)%phase_duration%getValue( phase_id )
        end do

    end function get_durations

    ! starts a phase that continues until a matching call to end_phase() is
    ! made.  phases identifiers are supplied so the start of a new phase is a
    ! cheap operation.
    subroutine start_phase( phase_id )

        !$ use omp_lib, only:    omp_get_ancestor_thread_num

        use precision, only:     int_64
        use stopwatch, only:     stopwatch_tick

        implicit none

        integer, intent(in)  :: phase_id

        integer              :: thread_number

        integer(kind=int_64) :: start_time

        call validate_phase_id( phase_id, "start_phase()" )

        thread_number = 1
        !$ thread_number = max( 1, omp_get_ancestor_thread_num( 1 ) + 1 )

        ! mark the start of a new phase for this thread.
        call measurements(thread_number)%active_phases%push( phase_id )

        ! update this phase's start time.
        start_time = measurements(thread_number)%phase_start%getValue( phase_id )
        call stopwatch_tick( start_time )
        call measurements(thread_number)%phase_start%setValue( phase_id, start_time )

#if ENABLE_TRACING
        block
            ! this block of code is useful to have fine-grained phase tracking
            ! when debugging or analyzing code performance.  that said, it adds
            ! non-trivial overhead to this routine so it is only enabled when
            ! explicitly requested.

            use print_h, only:    i_ioUnit

            write( i_ioUnit, "(A,A,I0,A,I0,A)" ) trim( measurements(thread_number)%phase_start%getKey( phase_id ) ), &
                 " started at ", start_time, " [thread #", thread_number, "]"
        end block
#endif

    end subroutine start_phase

    ! ends a phase started by a matching call to start_phase() and accumulates
    ! the elapsed time into the phase's duration.
    !
    ! NOTE: this routine stops execution if the phase specified does not match
    !       the currently active phase, as set by the most recent call to
    !       start_phase().
    !
    ! NOTE: ending a phase in a different OpenMP region than when start_phase()
    !       was called can either corrupt existing measurements or misattribute
    !       measurements to the wrong thread.  said another way, if start_phase()
    !       is called from a sequential portion of the code then end_phase()
    !       must be as well, and similarly for starting a phase in a parallel
    !       region.
    !
    subroutine end_phase( phase_id )

        !$ use omp_lib, only:            omp_get_ancestor_thread_num

        use precision, only:             int_64
        use stopwatch, only:             stopwatch_get_elapsed

        implicit none

        integer, intent(in)           :: phase_id

        integer                       :: thread_number
        integer                       :: active_phase_id

        integer(kind=int_64)          :: phase_start
        real                          :: phase_duration

        character(len=:), allocatable :: active_phase_name
        character(len=:), allocatable :: current_phase_name

        ! variables for printing each thread's active phases.
        integer                       :: this_thread_number, stack_level
        integer                       :: this_stack_depth
        integer                       :: this_phase_id
        character(len=:), allocatable :: this_phase_name

        if( .not. is_initialized() ) then
            write( 6, "(A)" ) "The measurement module has not been initialized before " // &
                 "calling end_phase()."
            stop
        end if

        call validate_phase_id( phase_id, "end_phase()" )

        thread_number = 1
        !$ thread_number = max( 1, omp_get_ancestor_thread_num( 1 ) + 1 )

        ! end the active phase.
        active_phase_id = measurements(thread_number)%active_phases%pop()

        ! ensure that we just ended the phase we expected to.  if not, bail so
        ! someone can fix the code.  if we blindly continue the timings will be
        ! invalid.
        if( active_phase_id .ne. phase_id ) then
            ! attempt to report the active phases for each of the threads so we
            ! can aide the user in identifying the root cause.  this may be
            ! jumbled if multiple threads get here at the same time (e.g. via
            ! incorrect identification of thread numbers) but it is better
            ! than nothing when tracking down a non-deterministic problem.
            do this_thread_number = 1, number_measurements
                this_stack_depth = measurements(this_thread_number)%active_phases%depth()

                write( *, "(A,I0,A,I0,A,/)" ) "Thread #", thread_number, "'s phase stack (", &
                     this_stack_depth, " active phases):"

                do stack_level = 1, measurements(this_thread_number)%active_phases%depth()
                    this_phase_id   = measurements(this_thread_number)%active_phases%get_element( stack_level )
                    this_phase_name = measurements(this_thread_number)%phase_start%getKey( this_phase_id )

                    write( *, "(A,I0,A,A,A,I0,A)" ) "    ", stack_level, ".  ", this_phase_name, &
                         " (", this_phase_id, ")"
                end do
                write( *, "(/)" )
            end do

            ! now report the mismatch and abort execution.
            current_phase_name = measurements(thread_number)%phase_start%getKey( phase_id )
            active_phase_name  = measurements(thread_number)%phase_start%getKey( active_phase_id )

            write( 6, "(A,A,A,A,A)" ) "Attempted to end phase '", &
                 trim( measurements(thread_number)%phase_start%getKey( phase_id ) ), &
                 "' but was in phase '", &
                 trim( measurements(thread_number)%phase_start%getKey( active_phase_id ) ), &
                 "'!"
            stop
        end if

        ! accumulate the elapsed time in this phase and add it to the running
        ! sum of its past durations.
        phase_start    = measurements(thread_number)%phase_start%getValue( phase_id )
        phase_duration = (measurements(thread_number)%phase_duration%getValue( phase_id ) + &
                          stopwatch_get_elapsed( phase_start ))
        call measurements(thread_number)%phase_duration%setValue( phase_id, phase_duration )

#if ENABLE_TRACING
        block

            real(kind=dp) :: this_phase_duration

            this_phase_duration = stopwatch_get_elapsed( phase_start )

            write( 6, "(A,A,F5.2,A,I0,A)" ) trim( measurements(thread_number)%phase_start%getKey( phase_id ) ), &
                 " took ", this_phase_duration, " seconds. [thread #", thread_number, "]"
        end block
#endif

    end subroutine end_phase

    ! clears the accumulated duration for a given phase so a call to
    ! get_duration()/get_durations() returns zero until the phase is measured
    ! again.  this is useful for ad hoc measurements where a section of code is
    ! timed and reported for multiple parameter configurations.
    subroutine reset_duration( phase_id )

        implicit none

        integer, intent(in) :: phase_id

        integer             :: thread_number

        if( .not. is_initialized() ) then
            write( 6, "(A)" ) "The measurement module has not been initialized before " // &
                 "calling end_phase()."
            stop
        end if

        call validate_phase_id( phase_id, "reset_duration()" )

        do thread_number = 1, number_measurements
            call measurements(thread_number)%phase_duration%setValue( phase_id, 0.0 )
        end do

    end subroutine reset_duration

    ! validates a phase identifier to ensure it is known by the module.  aborts
    ! execution if the specified identifier is unknown.
    !
    ! NOTE: this is for internal use only!
    !
    subroutine validate_phase_id( phase_id, routine_name )

        implicit none

        integer, intent(in)          :: phase_id
        character(len=*), intent(in) :: routine_name
        integer                      :: number_phases

        if( .not. is_initialized() ) then
            write( 6, "(A)" ) "The measurement module has not been initialized before " // &
                 "calling validate_phase_id()."
            stop
        end if

        number_phases = measurements(1)%phase_start%getSize()

        if( (phase_id < 1) .or. (phase_id > number_phases) ) then
            write( 6, "(A,I0,A,A,A,I0,A)" ) "Invalid phase identifier ", &
                 phase_id, " supplied to ", routine_name, ".  Must be in the range of [1, ", &
                 number_phases, "]."
            stop
        end if

    end subroutine validate_phase_id

end module measure

module profiling

    ! Specialized set of phases to measure NTLP's execution.  This extends
    ! the measure interface and provides NTLP-specific phase identifiers
    ! corresponding to its major components.  Additionally, a report method is
    ! provided that acquires the elapsed durations for said phases and pretty
    ! prints their timings in a hierarchical manner.
    !
    ! See the measure module for details on starting and ending phases, as well
    ! as other lower-level routines to interact with/control the measurements.

    ! Bring all of the measurement module's public interface
    ! (e.g. start_phase()/end_phase()) into scope.  We don't use this directly,
    ! but provide it to this module's users as a convenience.
    use measure

    ! Named identifiers, one per phase we measure during execution.  These are
    ! public so they're accessible to users of this module.  The purpose and
    ! descriptions for each phase are next to the identifier, below.
    integer, public :: &
         ! Time spent computing flow derivatives.
         measurement_id_derivatives, &

         ! Time spent computing the eddy viscosity and boundary conditions.
         measurement_id_eddy_viscosity_and_bcs, &

         ! Time spent in each of the three steps of the flow solve (comp_1(),
         ! comp_p(), and comp_2(), respectively).
         measurement_id_flow_solve_1, &
         measurement_id_flow_solve_p, &
         measurement_id_flow_solve_2, &

         ! Time spent performing humidity control.
         measurement_id_humidity, &

         ! Time spent writing outputs.
         measurement_id_io_histograms, &
         measurement_id_io_history, &
         measurement_id_io_particles, &
         measurement_id_io_traj, &
         measurement_id_io_pressure, &
         measurement_id_io_tecio, &
         measurement_id_io_viz, &

         ! Time spent advecting particles with the computed flow.
         measurement_id_particle_solver, &
         measurement_id_particle_reintro, &
         measurement_id_particle_diff, &
         measurement_id_particle_coalesce, &
         measurement_id_particle_fill_ext, &
         measurement_id_particle_loop, &
         measurement_id_particle_bcs, &
         measurement_id_particle_exchange, &
         measurement_id_particle_coupling, &
         measurement_id_particle_stats, &
         measurement_id_particle_coupling_exchange, &

         ! Time spent setting the simulation up.
         measurement_id_setup, &

         ! Time spent in the simulation.  This captures everything while MPI is
         ! initialized and available.
         measurement_id_solver, &

         ! Time spent in the main timestepping loop.
         measurement_id_timestepping_loop

    public :: &
         initialize_profiling, &
         report_profile, &
         shutdown_profiling

contains

    ! Initializes the base measure module and registers each of the phases known
    ! to NTLP's execution.
    !
    ! NOTE: This must be called before any of the other routines provided by
    !       this module.  It also should only be called at most once before
    !       shutdown_measurements() is called.  Once the module has been
    !       shutdown, it may be initialized again.
    !
    subroutine initialize_profiling()

        use measure, only:    create_phase, &
                              initialize_base_measurements => initialize_measurements

        implicit none

        call initialize_base_measurements()

        measurement_id_derivatives                         = create_phase( "calculating derivatives" )
        measurement_id_eddy_viscosity_and_bcs              = create_phase( "eddy viscosity and BCs" )
        measurement_id_flow_solve_1                        = create_phase( "flow solve comp1" )
        measurement_id_flow_solve_2                        = create_phase( "flow solve comp2" )
        measurement_id_flow_solve_p                        = create_phase( "flow solve comp_p" )
        measurement_id_humidity                            = create_phase( "humidity" )
        measurement_id_io_histograms                       = create_phase( "I/O - histograms" )
        measurement_id_io_history                          = create_phase( "I/O - history" )
        measurement_id_io_particles                        = create_phase( "I/O - particles" )
        measurement_id_io_traj                             = create_phase( "I/O - trajectories" )
        measurement_id_io_pressure                         = create_phase( "I/O - pressure field" )
        measurement_id_io_tecio                            = create_phase( "I/O - TecIO" )
        measurement_id_io_viz                              = create_phase( "I/O - viz" )
        measurement_id_particle_solver                     = create_phase( "particle_solver" )
        measurement_id_particle_reintro                    = create_phase( "particle_reintro" )
        measurement_id_particle_diff                       = create_phase( "particle_diff" )
        measurement_id_particle_coalesce                   = create_phase( "particle_coalesce" )
        measurement_id_particle_fill_ext                   = create_phase( "particle_fill_ext" )
        measurement_id_particle_loop                       = create_phase( "particle_loop" )
        measurement_id_particle_bcs                        = create_phase( "particle_bcs" )
        measurement_id_particle_exchange                   = create_phase( "particle_exchange" )
        measurement_id_particle_coupling                   = create_phase( "particle_coupling" )
        measurement_id_particle_coupling_exchange          = create_phase( "particle_coupling_exchange" )
        measurement_id_particle_stats                      = create_phase( "particle_stats" )
        measurement_id_setup                               = create_phase( "setup" )
        measurement_id_solver                              = create_phase( "solver" )
        measurement_id_timestepping_loop                   = create_phase( "solver time stepping loop" )

    end subroutine initialize_profiling

    subroutine report_profile( file_unit )

        use measure, only:    is_initialized, &
                              get_duration, &
                              print_duration

        implicit none

        integer, intent(in)  :: file_unit

        ! Measured durations.
        real                 :: duration_derivatives, &
                                duration_eddy_viscosity_and_bcs, &
                                duration_flow_solve_1, &
                                duration_flow_solve_2, &
                                duration_flow_solve_p, &
                                duration_humidity, &
                                duration_io_histograms, &
                                duration_io_history, &
                                duration_io_particles, &
                                duration_io_traj, &
                                duration_io_pressure, &
                                duration_io_tecio, &
                                duration_io_viz, &
                                duration_particle_solver, &
                                duration_particle_reintro, &
                                duration_particle_diff, &
                                duration_particle_coalesce, &
                                duration_particle_fill_ext, &
                                duration_particle_loop, &
                                duration_particle_bcs, &
                                duration_particle_exchange, &
                                duration_particle_coupling, &
                                duration_particle_coupling_exchange, &
                                duration_particle_stats, &
                                duration_setup, &
                                duration_solver, &
                                duration_timestepping_loop

        ! Computed durations.
        real                 :: total_duration, &
                                io_duration, &
                                particles_duration, &
                                flow_duration

        character, parameter :: newline = new_line( "a" )

        ! We can't report anything if we're not initialized.
        if( .not. is_initialized() ) then
            write( file_unit, "(A)" ) "The measurements module is not initialized.  Cannot report measurements!"

            return
        end if

        ! Get each phase's duration so we can report it and its percentage
        ! relative to the total duration.
        duration_derivatives                         = get_duration( measurement_id_derivatives )
        duration_eddy_viscosity_and_bcs              = get_duration( measurement_id_eddy_viscosity_and_bcs )
        duration_flow_solve_1                        = get_duration( measurement_id_flow_solve_1 )
        duration_flow_solve_2                        = get_duration( measurement_id_flow_solve_2 )
        duration_flow_solve_p                        = get_duration( measurement_id_flow_solve_p )
        duration_humidity                            = get_duration( measurement_id_humidity )
        duration_io_histograms                       = get_duration( measurement_id_io_histograms )
        duration_io_history                          = get_duration( measurement_id_io_history )
        duration_io_particles                        = get_duration( measurement_id_io_particles )
        duration_io_traj                             = get_duration( measurement_id_io_traj )
        duration_io_pressure                         = get_duration( measurement_id_io_pressure )
        duration_io_tecio                            = get_duration( measurement_id_io_tecio )
        duration_io_viz                              = get_duration( measurement_id_io_viz )
        duration_particle_solver                     = get_duration( measurement_id_particle_solver )
        duration_particle_reintro                    = get_duration( measurement_id_particle_reintro )
        duration_particle_diff                       = get_duration( measurement_id_particle_diff )
        duration_particle_coalesce                   = get_duration( measurement_id_particle_coalesce )
        duration_particle_fill_ext                   = get_duration( measurement_id_particle_fill_ext )
        duration_particle_loop                       = get_duration( measurement_id_particle_loop )
        duration_particle_bcs                        = get_duration( measurement_id_particle_bcs )
        duration_particle_exchange                   = get_duration( measurement_id_particle_exchange )
        duration_particle_coupling                   = get_duration( measurement_id_particle_coupling )
        duration_particle_coupling_exchange          = get_duration( measurement_id_particle_coupling_exchange )
        duration_particle_stats                      = get_duration( measurement_id_particle_stats )
        duration_setup                               = get_duration( measurement_id_setup )
        duration_solver                              = get_duration( measurement_id_solver )
        duration_timestepping_loop                   = get_duration( measurement_id_timestepping_loop )

        ! Sum the measured durations into higher-level phases.

        !
        ! NOTE: This is intended to capture any setup, pre-processing, the main
        !       solve and any post-processing.  If NTLP is only the main solve
        !       then there is no need to compute total_duration.
        !
        total_duration = ( &
             duration_solver &
             )

        ! Compute the aggregate time spent performing I/O, regardless of the
        ! file type or format.
        io_duration = ( &
             duration_io_histograms + &
             duration_io_history + &
             duration_io_particles + &
             duration_io_traj + &
             duration_io_pressure + &
             duration_io_tecio + &
             duration_io_viz &
             )

        particles_duration = ( &
             duration_particle_solver + &
             duration_particle_reintro + &
             duration_particle_diff + &
             duration_particle_coalesce &
             )

        flow_duration = ( &
             duration_flow_solve_1 + &
             duration_flow_solve_p + &
             duration_flow_solve_2 + &
             duration_eddy_viscosity_and_bcs + &
             duration_humidity &
             )

        write( file_unit, "(A,2A)" ) "Measurements report:", newline

        write( file_unit, "(A,2A)" )    "  Program run-time:", newline
        call print_duration( file_unit, "      Total:                         ", &
             duration_solver, total_duration )

        write( file_unit, "(A)" ) ""

        call print_duration( file_unit, "      Setup:                         ", &
             duration_setup, duration_solver )

        write( file_unit, "(A)" ) ""

        call print_duration( file_unit, "      Flow solver:                   ", &
             flow_duration, total_duration )

        call print_duration( file_unit, "          get_derv:                      ", &
             duration_derivatives, flow_duration )
        call print_duration( file_unit, "          comp1:                         ", &
             duration_flow_solve_1, flow_duration )
        call print_duration( file_unit, "          comp_p:                        ", &
             duration_flow_solve_p, flow_duration )
        call print_duration( file_unit, "          comp2:                         ", &
             duration_flow_solve_2, flow_duration )
        call print_duration( file_unit, "          Eddy viscosity/BCs:            ", &
             duration_eddy_viscosity_and_bcs, flow_duration )
        call print_duration( file_unit, "          Humidity control:              ", &
             duration_humidity, flow_duration )

        write( file_unit, "(A)" ) ""

        call print_duration( file_unit, "      Particles:                     ", &
             particles_duration, duration_solver )
        call print_duration( file_unit, "          particle_reintro:              ", &
             duration_particle_reintro, particles_duration )
        call print_duration( file_unit, "          particle_diff:                 ", &
             duration_particle_diff, particles_duration )
        call print_duration( file_unit, "          particle_coalesce:             ", &
             duration_particle_coalesce, particles_duration )
        call print_duration( file_unit, "          particle_solver:               ", &
             duration_particle_solver, particles_duration )
        call print_duration( file_unit, "               particle_fill_ext:              ", &
             duration_particle_fill_ext, duration_particle_solver )
        call print_duration( file_unit, "               particle_loop:                  ", &
             duration_particle_loop, duration_particle_solver )
        call print_duration( file_unit, "               particle_bcs:                   ", &
             duration_particle_bcs, duration_particle_solver )
        call print_duration( file_unit, "               particle_exchange:              ", &
             duration_particle_exchange, duration_particle_solver )
        call print_duration( file_unit, "               particle_coupling:              ", &
             duration_particle_coupling, duration_particle_solver )
        call print_duration( file_unit, "               particle_coupling_exchange:     ", &
             duration_particle_coupling_exchange, duration_particle_solver )
        call print_duration( file_unit, "               particle_stats:                 ", &
             duration_particle_stats, duration_particle_solver )

        write( file_unit, "(A)" ) ""

        call print_duration( file_unit, "      I/O:                           ", &
             io_duration, duration_solver )
        call print_duration( file_unit, "          Histograms:                    ", &
             duration_io_histograms, io_duration )
        call print_duration( file_unit, "          History:                       ", &
             duration_io_history, io_duration )
        call print_duration( file_unit, "          Particles:                     ", &
             duration_io_particles, io_duration )
        call print_duration( file_unit, "          Trajectories:                  ", &
             duration_io_traj, io_duration )
        call print_duration( file_unit, "          Pressure field:                ", &
             duration_io_pressure, io_duration )
        call print_duration( file_unit, "          TecIO:                         ", &
             duration_io_tecio, io_duration )
        call print_duration( file_unit, "          Viz:                           ", &
             duration_io_viz, io_duration )

        write( file_unit, "(A)" ) ""

    end subroutine report_profile

    ! Shuts down the module and releases the measurement resources allocated
    ! when it was initialized.  Additional calls to the module's routines,
    ! except another call to initialize_measurements() to re-initialize the
    ! module, are forbidden.  Each of the phase identifiers exposed by the
    ! module are set to zero to mark them as invalid.
    subroutine shutdown_profiling()

        use measure, only:    shutdown_base_measurements => shutdown_measurements

        implicit none

        call shutdown_base_measurements()

        ! Clear each of the phase identifiers so they can't accidentally be used
        ! after shutdown.
        measurement_id_derivatives                  = 0
        measurement_id_eddy_viscosity_and_bcs       = 0
        measurement_id_flow_solve_1                 = 0
        measurement_id_flow_solve_2                 = 0
        measurement_id_flow_solve_p                 = 0
        measurement_id_humidity                     = 0
        measurement_id_io_histograms                = 0
        measurement_id_io_history                   = 0
        measurement_id_io_particles                 = 0
        measurement_id_io_traj                      = 0
        measurement_id_io_pressure                  = 0
        measurement_id_io_tecio                     = 0
        measurement_id_io_viz                       = 0
        measurement_id_particle_solver              = 0
        measurement_id_particle_reintro             = 0
        measurement_id_particle_diff                = 0
        measurement_id_particle_coalesce            = 0
        measurement_id_particle_fill_ext            = 0
        measurement_id_particle_loop                = 0
        measurement_id_particle_bcs                 = 0
        measurement_id_particle_exchange            = 0
        measurement_id_particle_coupling            = 0
        measurement_id_particle_coupling_exchange   = 0
        measurement_id_particle_stats               = 0
        measurement_id_setup                        = 0
        measurement_id_solver                       = 0
        measurement_id_timestepping_loop            = 0

    end subroutine shutdown_profiling

end module profiling
