SUBROUTINE Calcul_time_spent(instant_init, instant_final, time_s)
    ! Note: instant_init/final follow DATE_AND_TIME format:
    ! (1)Year, (2)Month, (3)Day, (4)UTC Diff, (5)Hour, (6)Min, (7)Sec, (8)Ms

    Implicit none
    integer, dimension(8), INTENT(IN)  :: instant_init
    integer, dimension(8), INTENT(IN)  :: instant_final
    integer, dimension(4), INTENT(OUT) :: time_s
    
    integer :: d1, d2, month, year
    integer, dimension(12) :: days_in_month
    
    ! 1. Calculate absolute days from a reference (simplified)
    ! To handle month/year changes correctly, we calculate 'total days'
    d1 = instant_init(3)
    do month = 1, instant_init(2) - 1
        d1 = d1 + get_days_in_month(month, instant_init(1))
    end do
    
    d2 = instant_final(3)
    do month = 1, instant_final(2) - 1
        d2 = d2 + get_days_in_month(month, instant_final(1))
    end do
    
    ! Handle year difference for days
    if (instant_final(1) > instant_init(1)) then
        do year = instant_init(1), instant_final(1) - 1
            if (is_leap(year)) then
                d2 = d2 + 366
            else
                d2 = d2 + 365
            end if
        end do
    end if

    ! 2. Calculate raw differences
    time_s(4) = instant_final(7) - instant_init(7) ! Seconds
    time_s(3) = instant_final(6) - instant_init(6) ! Minutes
    time_s(2) = instant_final(5) - instant_init(5) ! Hours
    time_s(1) = d2 - d1                            ! Days

    ! 3. Adjust for negatives (The "Borrowing" Logic)
    if (time_s(4) < 0) then
        time_s(4) = time_s(4) + 60
        time_s(3) = time_s(3) - 1
    end if

    if (time_s(3) < 0) then
        time_s(3) = time_s(3) + 60
        time_s(2) = time_s(2) - 1
    end if

    if (time_s(2) < 0) then
        time_s(2) = time_s(2) + 24
        time_s(1) = time_s(1) - 1
    end if

CONTAINS

    function is_leap(y) result(res)
        integer, intent(in) :: y
        logical :: res
        res = (mod(y,4) == 0 .and. mod(y,100) /= 0) .or. (mod(y,400) == 0)
    end function

    function get_days_in_month(m, y) result(d)
        integer, intent(in) :: m, y
        integer :: d
        integer, dimension(12) :: m_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        d = m_days(m)
        if (m == 2 .and. is_leap(y)) d = 29
    end function

END SUBROUTINE Calcul_time_spent