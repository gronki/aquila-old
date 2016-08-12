module utils

    use types

contains

    subroutine coordgen(xx,yy)
        real, intent(out) :: xx(:,:), yy(size(xx,1),size(xx,2))
        integer :: x,y
        do y = 1, size(xx,2)
            do x = 1, size(xx,1)
                xx(x,y) = x
                yy(x,y) = y
            end do
        end do
    end subroutine

    subroutine grow_mask(mask,r)
        logical, intent(inout) :: mask(:,:)
        logical :: mask0(size(mask,1),size(mask,2))
        integer, intent(in) :: r
        integer :: x, y, x0, x1, y0, y1, nx, ny
        nx = size(mask,1)
        ny = size(mask,2)

        mask0 = mask

        do y = 1,ny
            do x = 1,nx
                if ( mask0(x,y) ) then
                    x0 = x - r
                    x1 = x + r
                    y0 = y - r
                    y1 = y + r
                    if ( x0 .lt. 1 ) x0 = 1
                    if ( x1 .gt. nx ) x1 = nx
                    if ( y0 .lt. 1 ) y0 = 1
                    if ( y1 .gt. ny ) y1 = ny
                    mask(x0:x1,y0:y1) = .true.
                end if
            end do
        end do

    end subroutine

    subroutine fill_mask(mask,master_mask)
        logical, intent(in) :: master_mask(:,:)
        logical, intent(inout) :: mask(size(master_mask,1),size(master_mask,2))
        logical :: mask0(size(master_mask,1),size(master_mask,2))

        do
            ! write (0,*) 'growing mask:', sum(l2i(mask))
            mask0 = mask
            call grow_mask(mask,1)
            mask = mask .and. master_mask
            if ( all(mask0.eqv.mask) ) exit
        end do
    end subroutine

    subroutine qblur(im,imout,r)
        real, intent(in) :: im(:,:)
        real, intent(out) :: imout(size(im,1),size(im,2))
        integer, intent(in) :: r
        integer :: x, y, x0, x1, y0, y1, nx, ny
        nx = size(im,1)
        ny = size(im,2)

        do y = 1,ny
            do x = 1,nx
                x0 = x - r
                x1 = x + r
                y0 = y - r
                y1 = y + r
                if ( x0 .lt. 1 ) x0 = 1
                if ( x1 .gt. nx ) x1 = nx
                if ( y0 .lt. 1 ) y0 = 1
                if ( y1 .gt. ny ) y1 = ny
                imout(x,y) = sum(im(x0:x1,y0:y1))   &
                        / ( (abs(x1-x0)+1) * (abs(y1-y0)+1) )
            end do
        end do

    end subroutine

    elemental function l2i(l) result(i)
        logical, intent(in) :: l
        integer :: i
        if (l) then
            i = 1
        else
            i = 0
        end if
    end function



end module
