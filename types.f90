module types

    implicit none

    integer, parameter :: dp = selected_real_kind(12)
    integer, parameter :: max_star = 300

    type starstruct
        real :: x,y
        real :: v
    end type


end module
