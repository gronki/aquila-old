module imutil

    integer, parameter :: dp = selected_real_kind(12)
    integer, parameter :: max_star = 300

    type starstruct
        real :: x,y
        real :: v
    end type



contains

    subroutine aqfindstar(im,list)
        real, intent(in) :: im(:,:)
        logical :: mask(size(im,1),size(im,2)),  master_mask(size(im,1),size(im,2))
        ! real :: xax(size(im,1),size(im,2)), yax(size(im,1),size(im,2))
        type(starstruct), intent(out) :: list(:)
        integer :: i, ix,iy,nx,ny, xymax(2)
        real :: smax, smean, sth
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
            write (0,'(A,2F14.1)') 'found star at: ', list(i) % x, list(i) % y
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

    subroutine apply_transform(xy,mx,xyout)
        type(starstruct), intent(in) :: xy(:)
        real, intent(in) :: mx(6,2)
        ! (n+1)(n+2)/2
        type(starstruct), intent(out) :: xyout(size(xy))
        xyout % x = mx(1,1) + mx(2,1) * (xy%x) + mx(3,1) * (xy%y)   &
            + mx(4,1) * (xy%x)**2 + mx(5,1) * (xy%x)*(xy%y) + mx(6,1) * (xy%y)**2
        xyout % x = mx(1,2) + mx(2,2) * (xy%x) + mx(3,2) * (xy%y)   &
            + mx(4,2) * (xy%x)**2 + mx(5,2) * (xy%x)*(xy%y) + mx(6,2) * (xy%y)**2
    end subroutine

    subroutine align_xyr(xy,xy0,mx)
        type(starstruct), intent(in) :: xy(:)
        type(starstruct), intent(in) :: xy0(:)
        type(starstruct) :: xy1(size(xy))
        real, intent(out) ::mx(6,2)
        integer :: nstar0,nstar,ii,i
        integer, parameter :: i_x = 1, i_y = 2, i_r = 3
        real(dp) :: y0, y0_dv(3), y0n_dv(3), v0(3)
        real(dp) :: lam, lambdas(100), y_dlam, len0, r0

        forall (i = 1:size(lambdas))
            lambdas(i) = 1.1 * (i-1.0)/(size(lambdas)-1.0) - 0.1
        end forall

        nstar = size(xy)
        nstar0 = size(xy0)
        len0 = sqrt(sum((xy0%x)**2 + (xy0%y)**2)/nstar0)
        r0 = len0

        k0 = 0.33 * sqrt(len0**2 / nstar0)

        print *, 'k0 = ', k0
        print *, 'r0 = ', r0

        v0 = 0

        do ii = 1,4

            call comp_ydv(v0,y0,y0_dv)
            y0n_dv = y0_dv / sqrt(sum(y0_dv**2))

            lam = -0.02

            call minimi(lam,real(r0/10))

            v0 = v0 + y0n_dv * lam

            if ( lam .lt. 1e-5 ) exit

        print '(3F10.4)',v0
            print *,'#############################'
        end do

        print '(3F10.4)',v0 * (/ 1d0, 1d0, 1/r0 /)

        call v2mx(v0,mx)

    contains

        subroutine minimi(x,scale)
            real(dp), intent(inout) :: x
            real, intent(in) :: scale
            integer :: i, ii
            real(dp) :: dx
            integer, parameter :: n1 = 13
            real(dp) :: v(3), y, y_dv(3), y_dx
            logical :: success
            success = .false.

            dx = scale

            do ii = 1,8

                do i = 1, n1
                    x = x + dx

                    v = v0 + y0n_dv * x
                    call comp_ydv(v,y,y_dv)
                    y_dx = sum(y_dv * y0n_dv)
                    print '(I5,F20.8,2Es14.6)', i, x, y, y_dx


                    if (y_dx .lt. 0) then
                        x = x - dx
                        dx = dx / n1
                        success = .true.
                        exit
                    end if
                end do
            end do


        end subroutine

        subroutine comp_ydv(v,y,y_dv)
            real(dp), intent(in) :: v(3)
            real(dp), intent(out) :: y,y_dv(3)
            real(dp) :: aa, y_dx1, y_dy1, x1_dr, y1_dr
            integer :: i0,i
            y = 0
            y_dv = 0
            xy1 % x = v(i_x) + cos(v(i_r)/r0)*(xy%x) - sin(v(i_r)/r0)*(xy%y)
            xy1 % y = v(i_y) + sin(v(i_r)/r0)*(xy%x) + cos(v(i_r)/r0)*(xy%y)
            xy1 % v = xy % v
            do i0=1,nstar0
                do i=1,nstar
                    aa = sqrt( &
                        (xy1(i)%x - xy0(i0)%x)**2 + (xy1(i)%y - xy0(i0)%y)**2 &
                        + k0**2 )
                    bb = sqrt(xy0(i0)%v * xy(i)%v ) * k0
                    y = y + bb / aa
                    y_dx1 = - (xy1(i)%x-xy0(i0)%x) * bb / aa**3
                    y_dy1 = - (xy1(i)%y-xy0(i0)%y) * bb / aa**3
                    x1_dr = - (sin(v(i_r)/r0) * (xy(i)%x) + cos(v(i_r)/r0) * (xy(i)%y)) / r0
                    y1_dr =   (cos(v(i_r)/r0) * (xy(i)%x) - sin(v(i_r)/r0) * (xy(i)%y)) / r0
                    y_dv(i_x) = y_dv(i_x) + y_dx1
                    y_dv(i_y) = y_dv(i_y) + y_dy1
                    y_dv(i_r) = y_dv(i_r) + y_dx1 * x1_dr + y_dy1 * y1_dr
                end do
            end do

        end subroutine

        subroutine v2mx(v,mx)
            real(dp), intent(in) :: v(3)
            real, intent(out) :: mx(6,2)
            mx = 0
            mx(1,1) = v(1)
            mx(2,1) = cos(v(i_r))
            mx(3,1) = -sin(v(i_r))
            mx(1,2) = v(2)
            mx(2,2) = sin(v(i_r))
            mx(3,2) = cos(v(i_r))
        end subroutine

    end subroutine

    elemental function mexha(x,y,sg) result(yf)
        real, intent(in) :: x,y,sg
        real :: yf,k
        real(dp), parameter :: pi = 4 * atan(1d0)
        k = (x**2 + y**2) / (2 * sg**2)
        yf =  (1 - k)  / (pi * sg**4 )  * exp( -k )
    end function

    subroutine wavelet_kernel(k,sg)
        real, intent(in) :: sg
        real, intent(out) :: k(:,:)
        integer :: szx,szy,i,j
        real :: x,y
        szx = size(k,1)
        szy = size(k,2)
        do j = 1,szy
            do i = 1,szx
                x = i - 0.5*(szx+1)
                y = j - 0.5*(szy+1)
                k(i,j) = mexha(x,y,sg)
            end do
        end do
    end subroutine

    subroutine convol_dumb(img,krn,imgout)
        real, intent(in) :: img(:,:), krn(:,:)
        real, intent(out) :: imgout(size(img,1),size(img,2))
        integer :: x,y,szx,szy, rx, ry
        integer ::  x1, y1, x2, y2
        integer ::  x1k, y1k, x2k, y2k
        szx = size(img,1)
        szy = size(img,2)
        kx = size(krn,1)
        ky = size(krn,2)
        if ( mod(kx,2) .eq. 0 .or. mod(ky,2) .eq. 0 ) then
            stop "Kernel must have uneven dimension"
        end if
        rx = floor(kx * 0.5)
        ry = floor(ky * 0.5)
        do y=1,szy
            do x=1,szx
                x1 = max( 1, x-rx )
                x1k = 1 + x1 - (x-rx)
                x2 = min( szx, x+rx )
                x2k = kx + x2 - (x+rx)
                y1 = max( 1, y-ry )
                y1k = 1 + y1 - (y-ry)
                y2 = min( szy, y+ry )
                y2k = ky + y2 - (y+ry)
                imgout(x,y) = sum( krn(x1k:x2k,y1k:y2k) * img(x1:x2,y1:y2) )
            end do
        end do

    end subroutine
    subroutine convol_dumb_trim(img,krn,imgout)
        real, intent(in) :: img(:,:), krn(:,:)
        real, intent(out) :: imgout(size(img,1),size(img,2))
        integer :: x,y,szx,szy, rx, ry
        integer ::  x1, y1, x2, y2
        integer ::  x1k, y1k, x2k, y2k
        szx = size(img,1)
        szy = size(img,2)
        kx = size(krn,1)
        ky = size(krn,2)
        if ( mod(kx,2) .eq. 0 .or. mod(ky,2) .eq. 0 ) then
            stop "Kernel must have uneven dimension"
        end if
        imgout = 0
        rx = floor(kx * 0.5)
        ry = floor(ky * 0.5)
        do y=1+ry,szy-ry
            do x=1+rx,szx-rx
                imgout(x,y) = sum( krn * img((x-rx):(x+rx),(y-ry):(y+ry)) )
            end do
        end do

    end subroutine


end module
