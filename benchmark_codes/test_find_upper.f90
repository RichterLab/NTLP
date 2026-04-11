program test_find_upper
  ! Checks that find_upper(arr, val) == minloc(arr, 1, mask=(arr .gt. val))
  ! for a variety of inputs.  Exits with a non-zero status on any mismatch.
  implicit none

  integer, parameter :: n = 256
  real               :: arr(n)
  real               :: val
  integer            :: i, expected, got, failures

  failures = 0

  ! ── Test 1: uniform grid ─────────────────────────────────────────
  do i = 1, n
    arr(i) = real(i)          ! 1.0, 2.0, ..., 256.0
  end do
  call run_sweep("uniform grid", arr, n, failures)

  ! ── Test 2: non-uniform grid (like a stretched vertical grid) ────
  do i = 1, n
    arr(i) = (real(i)/real(n))**1.5 * 3000.0   ! clustered near bottom
  end do
  call run_sweep("non-uniform grid", arr, n, failures)

  ! ── Test 3: specific edge cases ──────────────────────────────────

  ! Value exactly equal to an element — mask is .gt. so should NOT match
  arr = 0.0
  do i = 1, n
    arr(i) = real(i)
  end do
  val = 5.0
  expected = minloc(arr, 1, mask=(arr .gt. val))
  got      = find_upper(arr, val)
  call check("val equals element", expected, got, failures)

  ! Value below all elements
  val = 0.0
  expected = minloc(arr, 1, mask=(arr .gt. val))
  got      = find_upper(arr, val)
  call check("val below all", expected, got, failures)

  ! Value between first and second element
  val = 1.5
  expected = minloc(arr, 1, mask=(arr .gt. val))
  got      = find_upper(arr, val)
  call check("val between 1st and 2nd", expected, got, failures)

  ! Value equal to last element — no element is greater
  val = real(n)
  expected = minloc(arr, 1, mask=(arr .gt. val))
  got      = find_upper(arr, val)
  call check("val equals last element", expected, got, failures)

  ! ── Summary ──────────────────────────────────────────────────────
  if (failures == 0) then
    write(*,*) "All tests passed."
  else
    write(*,'(a,i0,a)') "FAILED: ", failures, " mismatches."
    stop 1
  end if

contains

  pure function find_upper(arr, val) result(idx)
    real, intent(in) :: arr(:)
    real, intent(in) :: val
    integer          :: idx
    integer          :: lo, hi, mid

    lo = 1
    hi = size(arr) + 1
    do while (lo < hi)
      mid = lo + (hi - lo) / 2
      if (arr(mid) <= val) then
        lo = mid + 1
      else
        hi = mid
      end if
    end do
    idx = lo
  end function find_upper

  subroutine check(label, expected, got, failures)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: expected, got
    integer,          intent(inout) :: failures
    if (expected /= got) then
      write(*,'(a,a,a,i0,a,i0)') "MISMATCH [", trim(label), &
        "]: minloc=", expected, "  find_upper=", got
      failures = failures + 1
    end if
  end subroutine check

  subroutine run_sweep(label, arr, n, failures)
    ! Tests every value from just below arr(1) to just above arr(n),
    ! plus exact hits on every element.
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: n
    real,             intent(in)    :: arr(n)
    integer,          intent(inout) :: failures

    real    :: val
    integer :: i, expected, got

    ! Sweep between every adjacent pair, plus exact element values
    do i = 1, n
      ! Exact hit on element i
      ! Skip i==n: minloc result is processor-dependent when mask is all-false;
      ! find_upper returns n+1 as a sentinel.  Tested separately below.
      if (i == n) cycle
      val      = arr(i)
      expected = minloc(arr, 1, mask=(arr .gt. val))
      got      = find_upper(arr, val)
      if (expected /= got) then
        write(*,'(a,a,a,i0,a,f12.4,a,i0,a,i0)') &
          "MISMATCH [", trim(label), "] i=", i, " val=", val, &
          " minloc=", expected, " find_upper=", got
        failures = failures + 1
      end if

      ! Midpoint between element i and i+1
      if (i < n) then
        val      = (arr(i) + arr(i+1)) / 2.0
        expected = minloc(arr, 1, mask=(arr .gt. val))
        got      = find_upper(arr, val)
        if (expected /= got) then
          write(*,'(a,a,a,i0,a,f12.4,a,i0,a,i0)') &
            "MISMATCH [", trim(label), "] mid i=", i, " val=", val, &
            " minloc=", expected, " find_upper=", got
          failures = failures + 1
        end if
      end if
    end do

  end subroutine run_sweep

end program test_find_upper
