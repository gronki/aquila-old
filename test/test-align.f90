program test_align

    use imutil

    implicit none

    type(starstruct), allocatable :: stars0(:), stars(:)
    integer :: nstar0,nstar,i
    real :: szx,szy
    real :: mx(6,2)

    read *,szx,szy
    read *,nstar0
    allocate(stars0(nstar0))
    do i = 1, nstar0
        read *, stars0(i) % v, stars0(i) % x, stars0(i) % y
    end do

    read *,nstar
    allocate(stars(nstar))
    do i = 1, nstar
        read *, stars(i) % v, stars(i) % x, stars(i) % y
    end do


    call align_xyr(stars0,stars,mx)


    deallocate(stars0)
    deallocate(stars)


end program test_align
