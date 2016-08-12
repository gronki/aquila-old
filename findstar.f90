module findstar

    use types
    use utils

    implicit none

contains


    subroutine aqfindstar(im,list)
        real, intent(in) :: im(:,:)
        logical :: mask(size(im,1),size(im,2)),  master_mask(size(im,1),size(im,2))
        ! real :: xax(size(im,1),size(im,2)), yax(size(im,1),size(im,2))
        type(starstruct), intent(out) :: list(:)
        integer :: i, ix,iy,nx,ny, xymax(2)
        real :: smax, smean, sthr
        real(dp) :: m_tx, m_ty, m_t
        real :: t

        nx = size(im,1)
        ny = size(im,2)

        ! calculate the threshold
        smax = maxval(im)
        smean = sum(im) * 1.0 / size(im)
        t = 0.1
        sthr =  t * smax + (1-t) * smean
        ! we consider only pixels brighter than half the maximum
        master_mask = (im .ge. sthr)

        list % x = 0
        list % y = 0
        list % v = -1
        do i = 1, size(list)
            ! if none left, exit
            if ( .not.any(master_mask) ) exit
            ! localize maximum pixel within good pixels
            xymax = maxloc(im, mask=master_mask)
            ! set this pixel to true in child mask
            mask = .false.
            mask(xymax(1),xymax(2)) = .true.
            ! make the child mask fill the entire blob
            call fill_mask(mask,master_mask)
            ! subtract the child mask from the major one
            master_mask = master_mask .and. (.not. mask)
            ! extend the child mask by one
            call grow_mask(mask,1)
            m_tx = 0
            m_ty = 0
            m_t = 0
            do iy = 1,ny
                do ix = 1,nx
                    if ( .not. mask(ix,iy) ) cycle
                    m_tx = m_tx + im(ix,iy) * 1.0 * ix
                    m_ty = m_ty + im(ix,iy) * 1.0 * iy
                    m_t = m_t + im(ix,iy)
                end do
            end do
            list(i) % v = m_t
            list(i) % x = m_tx / m_t - (nx+1)/2.0
            list(i) % y = m_ty / m_t - (ny+1)/2.0
            write (0,'(A,I0,A,2F14.1)') 'found star ', i,' at: ', list(i) % x, list(i) % y
        end do

    end subroutine


end module findstar
