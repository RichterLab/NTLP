module precision

    implicit none

    public

    ! Fixed length for "small" strings used in data structures.
    integer, parameter :: str = 64

    ! Kind representing a 64-bit signed integer.
    integer, parameter :: int_64 = selected_int_kind( 18 )

contains

    ! Compares two floating point values to see if they are close enough to be
    ! considered equal.  Close enough means that their difference is either
    ! within the module's absolute tolerance (dp_epsilon) or within a tight
    ! relative tolerance scaled by the maximum absolute value of both values.
    !
    ! This function is elemental so it can be applied element-by-element to
    ! arbitrary dimensional data (i.e. both scalars and arrays of arbitrary
    ! rank).
    !
    ! NOTE: Since this is a combined absolute and relative tolerance check, it
    !       is *NOT* equivalent to abs( lhs - rhs ) < tolerance.  Be careful.
    !
    elemental function dp_equals( lhs, rhs ) result( equality_flag )

        implicit none

        real, intent(in) :: lhs, rhs
        logical          :: equality_flag

        ! Tight tolerance for relative comparisons.  This is scaled against the
        ! larger, absolute magnitude of the two input values to set a looser
        ! tolerance than the module's absolute, while only being permissive of
        ! suitably close values when the inputs have large magnitudes.
        !
        ! NOTE: This was chosen arbitrarily so that regression tests pass but do
        !       not permit very large differences in values.  A more principled
        !       choice should be made at some point in the future and address
        !       the call sites that break because of it.
        !
        real, parameter  :: RELATIVE_TOLERANCE = 1e-15

        ! Local variables.
        real             :: difference

        difference = abs( lhs - rhs )

        ! Both values are considered "equal" if their absolute difference is
        ! smaller than the larger of:
        !
        !   1. The absolute tolerance
        !   2. The maximum absolute value of either, scaled by the relative
        !      tolerance
        !
        ! #1 is typically used when inputs are smallish and less than one, while
        ! #2 covers the other cases where the values become "large".
        equality_flag = difference <= max( epsilon( 1.0 ), &
                                           RELATIVE_TOLERANCE * max( abs( lhs ), &
                                                                     abs( rhs ) ) )

    end function dp_equals

end module precision

module Map_h

    use precision, only:    str

    implicit none

    private

    type, public :: Map
        private

        character(str), allocatable :: c_keys(:)

    contains

        final :: map_finalize

        procedure, pass, public :: &
             expand   => map_expand, &
             find     => map_find, &
             addKey   => map_addKey, &
             getKey   => map_getKey, &
             getSize  => map_getSize

    end type Map

contains

    subroutine map_finalize( this )

        implicit none

        type(Map) :: this

        if( allocated( this%c_keys ) ) deallocate( this%c_keys )

    end subroutine map_finalize

    subroutine map_expand( this )

        implicit none

        class(Map), intent(inout)   :: this

        ! Local variables.
        character(str), allocatable :: c_temp(:)
        integer                     :: i_size

        i_size = this%getSize()

        if( i_size == 0 ) then
            allocate( this%c_keys(1) )
        else
            call move_alloc( from=this%c_keys, to=c_temp )
            allocate( this%c_keys(i_size+1) )

            this%c_keys(1:i_size) = c_temp
        end if

    end subroutine map_expand

    subroutine map_addKey( this, c_key )

        implicit none

        class(Map), intent(inout)    :: this
        character(len=*), intent(in) :: c_key

        call this%expand()
        this%c_keys(size( this%c_keys )) = c_key

    end subroutine map_addKey

    function map_getKey( this, i_index ) result( c_key )

        implicit none

        class(Map), intent(in)        :: this
        integer, intent(in)           :: i_index
        character(len=:), allocatable :: c_key

        ! XXX: deal with unallocated c_keys
        c_key = this%c_keys(i_index)

    end function map_getKey

    pure function map_getSize( this ) result( i_size )

        implicit none

        class(Map), intent(in) :: this
        integer                :: i_size

        if( .not. allocated( this%c_keys ) ) then
            i_size = 0
        else
            i_size = size( this%c_keys )
        end if

    end function map_getSize

    ! Returns the index in the ordered map for a key of interest.  Returns -1 if
    ! the key was not found.
    pure function map_find( this, c_key ) result( i_index )

        implicit none

        class(Map), intent(in)       :: this
        character(len=*), intent(in) :: c_key
        integer                      :: i_index

        ! Local variables.
        integer                      :: i_current, i_last

        i_last = this%getSize()

        do i_current = 1, i_last
            if( c_key .eq. this%c_keys(i_current) ) then
                i_index = i_current
                return
            end if
        end do

        ! Indicate the key wasn't found.
        i_index = -1

    end function map_find

end module Map_h

module Int64map_h

    use Map_h, only:        Map
    use precision, only:    int_64

    implicit none

    private

    type, public, extends(Map) :: Int64map

        private

        integer(kind=int_64), allocatable :: i_values(:)

    contains

        final :: int64map_finalize

        procedure, pass, public :: &
             expand   => int64map_expand, &
             addEntry => int64map_addEntry

        procedure, pass, private ::  &
             int64map_getValueAtIndex, &
             int64map_getValueAtKey, &
             int64map_setValueAtIndex, &
             int64map_setValueAtKey

        generic, public :: &
             getValue => int64map_getValueAtIndex, &
             int64map_getValueAtKey
        generic, public :: &
             setValue => int64map_setValueAtIndex, &
             int64map_setValueAtKey

    end type Int64map

contains

    subroutine int64map_finalize( this )

        implicit none

        type(Int64map) :: this

        if( allocated( this%i_values ) ) deallocate( this%i_values )

    end subroutine int64map_finalize

    subroutine int64map_expand( this )

        implicit none

        class(Int64map), intent(inout)    :: this

        ! Local variables.
        integer(kind=int_64), allocatable :: i_temp(:)
        integer                           :: i_size

        i_size = this%getSize()

        if( i_size == 0 ) then
            allocate( this%i_values(1) )
        else
            call move_alloc( from=this%i_values, to=i_temp )
            allocate( this%i_values(i_size+1) )

            this%i_values(1:i_size) = i_temp
        end if

        call this%map%expand()
    end subroutine int64map_expand

    ! Add the key and value if it doesn't already exist in the map, otherwise
    ! ignore them.
    subroutine int64map_addEntry( this, c_key, i_value, l_keyFound )

        implicit none

        class(Int64map), intent(inout)   :: this
        character(len=*), intent(in)     :: c_key
        integer(kind=int_64), intent(in) :: i_value
        logical, optional, intent(out)   :: l_keyFound

        if( this%find( c_key ) > 0 ) then
            if( present( l_keyFound ) ) l_keyFound = .true.
            return
        else
            if( present( l_keyFound ) ) l_keyFound = .false.
            call this%addKey( c_key )
            this%i_values(this%getSize()) = i_value
        end if

    end subroutine int64map_addEntry

    pure function int64map_getValueAtIndex( this, i_index ) result( i_value )

        implicit none

        class(Int64map), intent(in) :: this
        integer, intent(in)         :: i_index
        integer(kind=int_64)        :: i_value

        i_value = this%i_values(i_index)

    end function int64map_getValueAtIndex

    ! Get the value of the supplied key, or return huge(1) if it doesn't.
    function int64map_getValueAtKey( this, c_key, l_keyFound ) result( i_value )

        implicit none

        class(Int64map), intent(in)    :: this
        character(len=*), intent(in)   :: c_key
        logical, optional, intent(out) :: l_keyFound
        integer(kind=int_64)           :: i_value

        ! Local variables.
        integer                        :: i_index

        i_index = this%find( c_key )
        if( i_index > 1 ) then
            if( present( l_keyFound ) ) l_keyFound = .true.
            i_value = this%i_values(i_index)
        else
            if( present( l_keyFound ) ) l_keyFound = .false.
            i_value = huge(1)
        end if

    end function int64map_getValueAtKey

    subroutine int64map_setValueAtIndex( this, i_index, i_value )

        implicit none

        class(Int64map), intent(inout)   :: this
        integer, intent(in)              :: i_index
        integer(kind=int_64), intent(in) :: i_value

        this%i_values(i_index) = i_value

    end subroutine int64map_setValueAtIndex

    ! Set the value of the supplied key.  Nothing is done if the key doesn't
    ! already exist.
    subroutine int64map_setValueAtKey( this, c_key, i_value, l_keyFound )

        implicit none

        class(Int64map), intent(inout)   :: this
        character(len=*), intent(in)     :: c_key
        integer(kind=int_64), intent(in) :: i_value
        logical, optional, intent(out)   :: l_keyFound

        ! Local variables.
        integer                          :: i_index

        i_index = this%find( c_key )
        if( i_index > 0 ) then
            if( present( l_keyFound ) ) l_keyFound = .true.
            this%i_values(i_index) = i_value
        else
            if( present( l_keyFound ) ) l_keyFound = .false.
        end if

    end subroutine int64map_setValueAtKey

end module Int64map_h

module IntStack_h

    implicit none

    public

    ! number of stack elements to allocate on first use.
    !
    ! NOTE: we don't need a large stack based on expected usage.
    !
    integer, parameter :: DEFAULT_NUMBER_ELEMENTS = 8

    type, public :: IntStack

        private

        integer :: current_depth
        integer :: maximum_depth

        integer, allocatable :: elements(:)

    contains

        procedure, pass, public :: &
             depth,       &
             get_element, &
             initialize,  &
             pop,         &
             push

    end type IntStack

    contains

        function depth( this ) result( stack_depth )

            implicit none

            class(IntStack), intent(in) :: this
            integer                     :: stack_depth

            stack_depth = this%current_depth

        end function depth

        function get_element( this, stack_position ) result( element )

            implicit none

            class(IntStack), intent(in) :: this
            integer, intent(in)         :: stack_position
            integer                     :: element

            if( stack_position > this%depth() ) then
                write( 6, "(A,I0,A,I0,A)" ) "Attempted to access stack element ", &
                     stack_position, " when the maximum element is ", &
                     this%depth(), "."
                stop
            end if

            element = this%elements(stack_position)

        end function get_element

        subroutine initialize( this, requested_depth )

            implicit none

            class(IntStack), intent(out)  :: this
            integer, optional, intent(in) :: requested_depth

            this%current_depth = 0

            if( present( requested_depth ) ) then
                this%maximum_depth = requested_depth
            else
                this%maximum_depth = DEFAULT_NUMBER_ELEMENTS
            end if

            allocate( this%elements(this%maximum_depth) )

            ! initialize everything to zero to improve debugging.
            !
            ! NOTE: don't do this if the stack is used in performance critical
            !       code.
            !
            this%elements = 0

        end subroutine initialize

        function pop( this ) result( element )

            implicit none

            class(IntStack), intent(inout) :: this
            integer                        :: element

            if( this%depth() == 0 ) then
                write( 6, "(A)" ) "Cannot pop an empty stack!"
            end if

            element = this%elements(this%depth())

            ! make it clear that this entry has been popped off.  this has
            ! no bearing on the externally visible behavior but is useful
            ! during debugging.
            this%elements(this%current_depth) = -1

            this%current_depth = this%current_depth - 1

        end function pop

        subroutine push( this, element, element_position )

            implicit none

            class(IntStack), intent(inout) :: this
            integer, intent(in)            :: element
            integer, optional, intent(out) :: element_position

            integer, allocatable           :: original_elements(:)

            ! grow the internal array if we're out of space.
            if( this%current_depth == this%maximum_depth ) then
                original_elements = this%elements

                deallocate( this%elements )
                allocate( this%elements(this%maximum_depth * 2) )

                this%elements(1:this%maximum_depth)                        = original_elements
                this%elements(this%maximum_depth+1:this%maximum_depth * 2) = 0

                this%maximum_depth = this%maximum_depth * 2
            end if

            this%current_depth                = this%current_depth + 1
            this%elements(this%current_depth) = element

            if( present( element_position ) ) then
                element_position = this%current_depth
            end if

        end subroutine push

end module IntStack_h

module RealMap_h

    use map_h, only:        Map

    implicit none

    private

    type, public, extends( Map ) :: RealMap

        private

        real, allocatable :: dp_values(:)

    contains

        final :: realMap_finalize

        procedure, pass, public :: &
             expand   => realMap_expand, &
             addEntry => realMap_addEntry

        procedure, pass, private :: &
             realMap_getValueAtIndex, &
             realMap_getValueAtKey, &
             realMap_setValueAtIndex, &
             realMap_setValueAtKey

        generic, public :: &
             getValue => realMap_getValueAtIndex, &
             realMap_getValueAtKey
        generic, public :: &
             setValue => realMap_setValueAtIndex, &
             realMap_setValueAtKey

    end type RealMap

contains

    subroutine realmap_finalize( this )

        implicit none

        type(RealMap) :: this

        if( allocated( this%dp_values ) ) deallocate( this%dp_values )

    end subroutine realmap_finalize

    subroutine realmap_expand( this )

        implicit none

        class(RealMap), intent(inout) :: this

        ! Local variables.
        real, allocatable             :: dp_temp(:)
        integer                       :: i_size

        i_size = this%getSize()

        if( i_size == 0 ) then
            allocate( this%dp_values(1) )
        else
            call move_alloc( from=this%dp_values, to=dp_temp )
            allocate( this%dp_values(i_size+1) )

            this%dp_values(1:i_size) = dp_temp
        end if

        call this%map%expand()

    end subroutine realmap_expand

    ! Add the key and value if it doesn't already exist in the map, otherwise
    ! ignore them.
    subroutine realmap_addEntry( this, c_key, dp_value, l_keyFound )

        implicit none

        class(RealMap), intent(inout)  :: this
        character(len=*), intent(in)   :: c_key
        real, intent(in)               :: dp_value
        logical, optional, intent(out) :: l_keyFound

        if( this%find( c_key ) > 0 ) then
            if( present( l_keyFound ) ) l_keyFound = .true.
            return
        else
            call this%addKey( c_key )
            this%dp_values( this%getSize() ) = dp_value
        end if

    end subroutine realmap_addEntry

    pure function realmap_getValueAtIndex( this, i_index ) result( dp_value )

        implicit none

        class(RealMap), intent(in) :: this
        integer, intent(in)        :: i_index
        real                       :: dp_value

        dp_value = this%dp_values(i_index)

    end function realmap_getValueAtIndex

    ! Get the value of the supplied key, or return positive infinity if it
    ! doesn't.
    function realmap_getValueAtKey( this, c_key, l_keyFound ) result( dp_value )

        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf

        implicit none

        class(RealMap), intent(in)     :: this
        character(len=*), intent(in)   :: c_key
        logical, optional, intent(out) :: l_keyFound
        real                           :: dp_value

        ! Local variables.
        integer                        :: i_index

        i_index  = this%find(c_key)
        if(i_index > 0) then
            if(present(l_keyFound)) l_keyFound = .true.
            dp_value = this%dp_values(i_index)
        else
            if(present(l_keyFound)) l_keyFound = .false.
            dp_value = ieee_value( 0.0d0, ieee_positive_inf )
        end if

    end function realmap_getValueAtKey

    subroutine realmap_setValueAtIndex( this, i_index, dp_value )

        implicit none

        class(RealMap), intent(inout) :: this
        integer, intent(in)           :: i_index
        real, intent(in)              :: dp_value

        this%dp_values(i_index) = dp_value

    end subroutine realmap_setValueAtIndex

    ! Set the value of the supplied key.  Nothing is done if the key doesn't
    ! already exist.
    subroutine realmap_setValueAtKey( this, c_key, dp_value, l_keyFound )

        implicit none

        class(RealMap), intent(inout)  :: this
        character(len=*), intent(in)   :: c_key
        real, intent(in)               :: dp_value
        logical, optional, intent(out) :: l_keyFound

        ! Local variables.
        integer                        :: i_index

        i_index = this%find( c_key )
        if( i_index > 0 ) then
            if( present( l_keyFound ) ) l_keyFound = .true.
            this%dp_values(i_index) = dp_value
        else
            if( present( l_keyFound ) ) l_keyFound = .false.
        end if

    end subroutine realmap_setValueAtKey

end module RealMap_h
