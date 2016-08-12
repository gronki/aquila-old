module align

    use types
    use utils

    implicit none

contains


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


end module align
