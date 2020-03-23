! atom tools
module jcor
    implicit none
    public
    type atom
        character*5 :: resn,atyp,chai
        integer :: iato,icha,ires
        real*8 :: v(3)
    end type atom
contains
    ! operate

    ! translation with v
    function transV(corIn,v)
        implicit  none
        real*8,intent(in):: v(3)
        type(atom),intent(in):: corIn(:)
        type(atom):: transV(size(corIn))
        integer i
        do i=1,size(corIn)
            transV(i)%v=corIn(i)%v+v
            transV(i)%atyp=corIn(i)%atyp
        enddo
    end function transV

    ! set an atom as origin
    function setOriginAtom(corIn,a)
        implicit none
        type(atom),intent(in):: corIn(:),a
        type(atom):: setOriginAtom(size(corIn))
        setOriginAtom=transV(corIn,-1.d0*a%v)
    end function setOriginAtom

    ! rotate by axis/vector and angle, angle is in degree
    function rotateVA(corIn,axis,angle)
        implicit none
        type(atom),intent(in):: corIn(:)
        real*8,intent(in) :: angle, axis(3)
        type(atom):: rotateVA(size(corIn))
        real*8 :: M(3,3)
        real*8 :: rA,pi,sinA,cosA,x,y,z,radAngle
        integer :: i
        pi=acos(-1.d0)
        radAngle=angle/180.d0*pi
        cosA=dcos(radAngle)
        sinA=dsin(radAngle)
        rA=1-cosA
        x=axis(1)
        y=axis(2)
        z=axis(3)
        ! rotate matrix
        M(1,1)=cosA+rA*x*x
        M(1,2)=rA*x*y-sinA*z
        M(1,3)=rA*x*z+sinA*y

        M(2,1)=rA*y*x+sinA*z
        M(2,2)=cosA+rA*y*y
        M(2,3)=rA*y*z-sinA*x

        M(3,1)=rA*z*x-sinA*y
        M(3,2)=rA*z*y+sinA*x
        M(3,3)=cosA+rA*z*z
        !

        ! rotate
        do i=1,size(corIn)
            rotateVA(i)%v = matmul(M,corIn(i)%v)
            rotateVA(i)%atyp = corIn(i)%atyp
        enddo

    end function rotateVA
        



    ! measure

    ! get azimuth angle, project a1, a2 to xOy plane, a1 set as origin point
    function getPhi(a1,a2)
    implicit none
    type(atom),intent(in):: a1,a2
    type(atom):: a1p,a2p,ex
    real*8 :: pi,getPhi
    pi=dacos(-1.d0)

    ! set a1=#origin
    a1p%v=0

    ! trans and set z=0
    a2p%v=a2%v-a1%v
    if (dabs(a2p%v(1)*a2p%v(2)).lt.0.00001d0)  then
        print *, "getPhi warning: Close Atoms!"
    endif

    a2p%v(3)=0

    ex%v=(/1.d0,0.d0,0.d0/)

    getPhi=getAngle(a2p,a1p,ex)
    if (a2p%v(2).lt.0.d0) getPhi=360.d0-getPhi
    end function getPhi

    ! calculate angleTheta, angle(zAxis,N,C)
    function getTheta(a1,a2)
        implicit none
        real*8 ::  getTheta
        type(atom):: a1,a2
        real*8,parameter:: PI=dcos(-1.d0)
        integer :: i,j,k
        type(atom) :: atmp,a1p
        a1p%v(1:2)=a1%v(1:2)
        a1p%v(3)=a1%v(3)+1
        getTheta= getAngle(a1p,a1,a2)

    end  function getTheta

    function adist(a1,a2)
        type(atom) a1,a2
        real*8 adist
        adist=sqrt(&
            &(a1%v(1)-a2%v(1))*(a1%v(1)-a2%v(1))+&
            &(a1%v(2)-a2%v(2))*(a1%v(2)-a2%v(2))+&
            &(a1%v(3)-a2%v(3))*(a1%v(3)-a2%v(3)))
    end function adist
    ! get angle, in rad
    function aangl(a1,a2,a3)
        ! a2 is the tip 
        type(atom) a1,a2,a3,x0
        real*8 v1(3),v2(3),aangl
        real*8,parameter :: pi=dacos(-1.d0)
        x0%v=0
        v1=a2%v-a1%v
        v2=a2%v-a3%v
        aangl=dacos((v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))/adist(a1,a2)/&
            &adist(a2,a3))
    end function aangl
    
    ! get angle in degree
    function getAngle(a1,a2,a3)
        ! a2 is the tip 
        type(atom) a1,a2,a3,x0
        real*8 v1(3),v2(3),getAngle
        real*8,parameter :: pi=dacos(-1.d0)
        getAngle=aangl(a1,a2,a3)/pi*180.d0
    end function getAngle


    function adihe(a1,a2,a3,a4)
        type(atom) a1,a2,a3,a4,x1,x2,x0
        real*8 v21(3),v23(3),v32(3),v34(3),adihe,n1(3),n2(3)
        v21=a1%v-a2%v
        v23=a3%v-a2%v
        v34=a4%v-a3%v
        v32=-1*v23
        ! cross product
        n1=(/v23(2)*v21(3)-v23(3)*v21(2),&
            &v23(3)*v21(1)-v23(1)*v21(3),&
            &v23(1)*v21(2)-v23(2)*v21(1)/)
        n2=(/v34(2)*v32(3)-v34(3)*v32(2),&
            &     v34(3)*v32(1)-v34(1)*v32(3),&
            &     v34(1)*v32(2)-v34(2)*v32(1)/)
        x0%v=0
        x1%v=n1
        x2%v=n2
        adihe=aangl(x1,x0,x2)
    end function adihe
    ! read data
    subroutine rl_cor(a,line)
        character*200 :: line
        integer :: i,j,k
        type(atom) :: a
        read(line,999) a%iato,a%ires,a%resn,a%atyp,a%v(1),a%v(2),a%v(3),&
            &a%chai,a%icha
        999   format(2I10,2x,a5,5x,a5,3x,3f20.10,2x,a5,5x,I10)
    end subroutine rl_cor
    subroutine rl_xyz(a,line)
        character*200 :: line,line2
        integer :: i,j,k
        type(atom) :: a
        read(line,'(a4,a194)')a%atyp(1:4),line2(1:194)
        read(line2,*) a%v
    end subroutine rl_xyz
    ! show data
    subroutine showAtom(a)
        type(atom):: a
        write(*,'(a4,1x,3f8.2)')a%atyp(1:4),a%v
    end subroutine showAtom

    ! save data in xyz format
    subroutine saveCorXyz(corIn,fout)
        implicit none
        integer i,NAtom
        type(atom):: corIn(:)
        character*200::  fout
        NAtom=size(corIn)
        open(36,file=fout)
        write(36,*)NAtom
        write(36,*)
        do i=1,NAtom    
            write(36,'(a4,1x,3f15.5)') corIn(i)%atyp(1:4),corIn(i)%v
        enddo
        close(36)
    end subroutine saveCorXyz
end module jcor
