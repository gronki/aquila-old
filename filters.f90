module filters

    use types
    use utils

    implicit none

contains

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
        integer :: x,y,szx,szy, rx, ry, kx, ky
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
        integer :: x,y,szx,szy, rx, ry, kx, ky
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
