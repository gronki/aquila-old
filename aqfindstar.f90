program program_aqfindstar

    use imutil
    use iso_fortran_env

    implicit none

    real(real32), allocatable :: img(:,:), img2(:,:), data(:,:,:)
    type(starstruct) :: list(2**16)
    integer :: i

    integer :: status, bsize, sz(3), ndim
    logical :: anyf
    ! real(real32), allocatable :: img(:,:)
    character(len=50) :: errstr, comment
    character(len=1200) :: fn
    character(len=1200) :: buf
    real :: krn(25,25), thr
    read (*,'(A)') fn
    ! print *,trim(fn)
    status = 0
    bsize = 1
    call ftopen(33,trim(fn),0,bsize,status)

    if ( status .ne. 0 ) then
        call ftgerr(status,errstr)
        write (0,*) 'error: ', trim(errstr)
        stop -1
    end if

    call ftgidm(33,ndim,status)
    if ( ndim .eq. 2 ) then
        call ftgisz(33,2,sz(1:2),status)
        sz(3) = 1
    else if ( ndim .eq. 3 ) then
        call ftgisz(33,3,sz,status)
    else
        stop "Unknown image dimension"
    end if

    ! allocate a generic 3D array
    allocate( data(sz(1),sz(2),sz(3)) )
    ! read the data from file
    call ftgpve(33,1,1,product(sz),0,data,anyf,status)
    call ftclos(33,status)

    ! print *,sz
    allocate( img(sz(1),sz(2)), img2(sz(1),sz(2)) )



    call wavelet_kernel(krn,3.0)
    call convol_dumb_trim(img,krn,img2 )
    thr = sum(img2**2) / ( product(sz) )
    ! thr = maxval(img2) * 0.025
    where (img2 .lt. thr)  img2 = 0
    call aqfindstar(img2,list)
    deallocate(img2)

    write (6,'(2I10)') sz(1:2)
    write (6,'(I10)') size(list)

    do i = 1,size(list)
        if ( list(i) % v .lt. 0 ) exit
        write (6,'(3F15.2)') list(i) % v, list(i) % x, list(i) % y
    end do
    deallocate(img)
end program
