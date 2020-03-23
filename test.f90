! make 2an at any angle you want
! test orientation of rotLight

program main
    use jcor
    implicit none
    integer,parameter :: ax(3)=(/1,0,0/)
    integer,parameter :: ay(3)=(/0,1,0/)
    integer,parameter :: az(3)=(/0,0,1/)
    real*8,parameter:: lx=1.0/sqrt(2.0),la=1.d0,lb=sqrt(3.d0)
    real*8  :: rtmp
    character*20 :: r2c
    type(atom):: aAtom,oAtom,p1,p2,p3,p4

    oAtom%v=(/0.d0,0.d0,0.d0/)
    p1%v=(/lb,la,1.d0/)
    p2%v=(/-lb,la,2.d0/)
    p3%v=(/-lb,-la,3.d0/)
    p4%v=(/lb,-la,4.d0/)

    print *, getPhi(oAtom,p1), getPhi(oAtom,p2), getPhi(oAtom,p3), getPhi(oAtom,p4) 
    print *, getPhi(oAtom,oAtom)

end program main
