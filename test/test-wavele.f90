program test_wavele

    use iso_fortran_env
    use imutil

    implicit none

    integer :: status, bsize, naxes(2)
    logical :: anyf
    real(real32), allocatable :: img(:,:),img2(:,:)
    character(len=50) :: errstr, comment
    real :: krn(33,33), thr

    integer :: i,j
    character(len=80) :: row

    status = 0

    call ftopen(33,'fit/al1.fit',0,bsize,status)

    ! call ftgkyj(33,'NAXIS1',naxes(1),comment,status)
    ! call ftgkyj(33,'NAXIS2',naxes(2),comment,status)

    call ftgisz(33,2,naxes,status)
    allocate( img(naxes(1),naxes(2)), img2(naxes(1),naxes(2)) )

    call ftgpve(33,1,1,naxes(1)*naxes(2),0,img,anyf,status)
    call ftclos(33,status)

    print *, 'wymiary obrazka:', naxes

    call wavelet_kernel(krn,3.0)
    call convol_dumb_trim(img,krn,img2 )
    thr = sum(img2**2) / ( product(naxes) )
    ! thr = maxval(img2) * 0.025
    where (img2 .lt. thr)  img2 = 0

    call ftinit(33,'temp/wavele.fit',bsize,status)
    if ( status .ne. 0 ) then
        call ftgerr(status,errstr)
        print *, 'Error: ', errstr
    end if
    call ftphpr(33,.true.,-32,2,naxes,0,1,.true.,status)
    call ftppre(33,1,1,naxes(1)*naxes(2),img2,status)
    call ftclos(33, status)

    deallocate(img,img2)


end program
