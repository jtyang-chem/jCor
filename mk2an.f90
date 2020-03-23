! make 2an at any angle you want
program main
    use jcor
    use vCommon
    implicit none
    real*8,parameter::  vx(3)=(/1.d0,0.d0,0.d0/), vy(3)=(/0.d0,1.d0,0.d0/), vz(3)=(/0.d0,0.d0,1.d0/)
    character*200::fout
    integer,parameter :: NAn=3, AtomAn=6 
    integer,parameter :: NOutCor=NAn*AtomAn
    type(atom):: p0(AtomAn),p1(AtomAn),p2(AtomAn),p3(AtomAn),triP(NOutCor)
    type(atom):: a1,a2,a3
    call readin()
    totcor= setOriginAtom(totcor,totcor(6))
    totcor= rotateVA(totcor,vx,90.d0)
    ! save origin cor
    p0=totcor

    p1= rotateVA(p0,vy,50.d0)
    p1= rotateVA(p1,vz,30.d0)
    p1= transV(p1,(/0.d0,3.2d0,0.d0/))

    p2= rotateVA(p0,vy,50.d0)
    p2= rotateVA(p2,vz,150.d0)

    triP(1:AtomAn)=p2
    triP(AtomAn+1:2*AtomAn)=p1
    triP(2*AtomAn+1:3*AtomAn)=transV(p2,(/-6.3d0,0.0d0,0.d0/))

    fout="triP.xyz"
    call saveCorXyz(triP,fout)
end program main

