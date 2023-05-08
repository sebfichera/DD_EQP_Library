! 2D SQUARE ELEMENT CLASS V 1.0
!
! BY:
! SEBASTIANO FICHERA,
! sebastiano.fichera@polito.it
! GREGORIO MARIGGIO',
! gregorio.mariggio@polito.it
!
! File containing the 2D square element Class. The Class contains all the
! necessary methods to set the element coordinates in a 2D global reference
! system, to map the element into a standard [-1,1] parent reference system,
! to manage the discontinuities using the DD EQP library and to perform
! the integration of all subdomains.
! The file also contains the module 'parameters' to define parameters needed by
! the element Class and the module 'interfaces' defining the interaces of the
! functions from the DD EQP library.
!
! USAGE
! Declare a type(Quad) variable in the main program and call the method EvalQuad
! or EvalQuadFromFile, passing the Quad variable defined above. Note that if the
! EvalQuadFromFile is used, it is necessary to pass to the method also the input
! parameters read from the file (element coordinates, discontinuity coefficients,
! number of discontinuities and subdomain to integrate).
! See the Library example file main.f90 for usage hints.
!
! Copyright (C) 2022 Equivalent-Polynomials.net
!
! This file is part of DD EQP Library.
!
! DD EQP Library is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DD EQP Library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DD EQP Library.  If not, see <http://www.gnu.org/licenses/>.

! ------------------------ PARAMETERS MODULE ------------------------------
module parameters

    implicit none
    integer, parameter :: rk = selected_real_kind(p=15,r=307)
    real(kind=rk), parameter :: Pi = 4 * atan (1.0_8)
    real(kind=rk), parameter :: tol = 1E-4

end module parameters

! ------------------------ INTERFACES BLOCK MODULE ------------------------
module interfaces

    interface

    ! Functions from DD EQP library
    double precision function HEqPol(x,y,z,eqcv,etype)
        implicit none
        integer :: etype
        double precision :: x,y,z,eqcv(:)
    end function HEqPol

    subroutine Heqpol_coefficients(a,b,c,d,x,y,z,eqcv,etype)
        implicit none
        integer :: etype
        double precision :: a,b,c,d,x(:),y(:),z(:),eqcv(:)
    end subroutine Heqpol_coefficients

    subroutine DD_Heqpol_coefficients(a,b,c,eqcv,etype,part)
        implicit none
        integer :: etype
        character :: part
        double precision :: a(2),b(2),c(2),eqcv(:)
    end subroutine DD_Heqpol_coefficients

    subroutine DiscontMapping(a,b,c,d,x,y,z,etype)
        implicit none
        integer :: etype
        double precision :: a,b,c,d,L1,L2,L3,L4
        double precision,allocatable :: xi(:),eta(:),zeta(:)
        double precision :: x(:),y(:),z(:)
    end subroutine DiscontMapping

    end interface

end module interfaces

! ------------------------ QUAD CLASS -------------------------------------
module class_Quad

    use parameters
    use interfaces

    implicit none
    private
    public :: Quad, SetGaussPts, SetCoords, SetDiscont, CheckSelection, NodeScheme, PrintRes
    public :: HowManyDiscont, RearrangingCoefficients, SetCoordsFromFile, SetNumDscnt, SetDiscontFromFile
    public :: detJ, EvalN, Mapping, Integrate, EvalQuad, EvalQuadFromFile

    integer, parameter :: etype=21,pbdim=2,nodes=4,theSize=3,ndscnt=2
    real(kind=rk), parameter :: xi(nodes)=reshape((/-1.d0, 1.d0, 1.d0, -1.d0/),shape(xi))
    real(kind=rk), parameter :: eta(nodes)=reshape((/-1.d0, -1.d0, 1.d0, 1.d0/),shape(eta))
    real(kind=rk), parameter :: gp(theSize)=reshape((/-0.77459667,0.,0.77459667/),shape(gp))
    real(kind=rk), parameter :: gw(theSize)=reshape((/0.55555555,0.88888889,0.55555555/),shape(gw))

    type Quad
        integer :: dscnt
        real(kind=rk) :: agbl(ndscnt),bgbl(ndscnt),cgbl(ndscnt),apar(ndscnt),bpar(ndscnt),cpar(ndscnt)
        real(kind=rk) :: w(theSize*pbdim),eqcv(theSize*pbdim),t(theSize*pbdim),out_vec(theSize*pbdim)
        real(kind=rk) :: x(nodes),y(nodes)
        real(kind=rk) :: gbl_gpts(theSize*theSize,pbdim)
        real(kind=rk) :: gpts(theSize*theSize,pbdim)
        real(kind=rk) :: gwts(theSize*theSize,pbdim)
        real(kind=rk) :: N(theSize*theSize,nodes)
    end type Quad

    contains

    subroutine NodeScheme()
        implicit none
        write(*,*)
        write(*,*)'Insert element coordinates following the scheme shown below:'
        write(*,*)
        write(*,*)'4-------------3'
        write(*,*)'|             |'
        write(*,*)'|             |'
        write(*,*)'|             |'
        write(*,*)'|             |'
        write(*,*)'|             |'
        write(*,*)'1-------------2'
        write(*,*)
    end subroutine NodeScheme

    subroutine SetGaussPts(this)
        implicit none
        type(Quad), intent(out) :: this
        integer :: i,j,cnt=1
        do i=1,theSize
            do j=1,theSize
                this%gpts(cnt,1)=gp(i)
                this%gpts(cnt,2)=gp(j)
                this%gwts(cnt,1)=gw(i)
                this%gwts(cnt,2)=gw(j)
                cnt=cnt+1
            end do
        end do
    end subroutine SetGaussPts

    subroutine HowManyDiscont(this_w)
        implicit none
        type(Quad), intent(out) :: this_w
        integer :: selection
        write(*,*)
        write(*,'(1X,A)',advance='no') 'How many discontinuities exist within the same element? '
        read(*,*) selection
        call CheckSelection(selection,1,2)
        this_w%dscnt=selection
    end subroutine HowManyDiscont

    subroutine SetDiscont(this_r,this_w,i)
        implicit none
        type(Quad), intent(in) :: this_r
        type(Quad), intent(out) :: this_w
        integer :: i
        real(kind=rk) :: a,b,c,d,z(size(this_r%x))
        write(*,'(1X,A)',advance='no') 'Discontinuity no. '
        write(*,'(I1)',advance='no') i
        write(*,'(1X,A)',advance='no') 'of '
        write(*,'(I1)',advance='no') this_r%dscnt
        write(*,*) ': '
        write(*,'(1X,A)',advance='no') 'a,b,c : '
        read(*,*) a,b,c
        do while ((a == 0) .and. (b == 0))
            write(*,'(1X,A)',advance='no') 'WARNING: "a" and "b" cannot be both equal '&
                                            & // 'to 0 at the same time. Please choose '&
                                            & // 'a,b,c again : '
            read (*,*) a,b,c
        end do
        if (((abs(a).lt.tol).and.(abs(a).gt.0)).or.&
           &((abs(b).lt.tol).and.(abs(b).gt.0)).or.&
           &((abs(c).lt.tol).and.(abs(c).gt.0))) then
            write(*,'(1X,A)',advance='no') 'WARNING: some coefficients are too close to 0. '&
                                            & // 'They will be set to '
            write(*,'(F7.4)',advance='no') tol
            write(*,'(1X,A)',advance='no') ' (max tolerance)'
            write(*,*)
            if ((abs(a).lt.tol).and.(abs(a).gt.0)) a = tol
            if ((abs(b).lt.tol).and.(abs(b).gt.0)) b = tol
            if ((abs(c).lt.tol).and.(abs(c).gt.0)) c = tol
        end if
        this_w%agbl(i)=a
        this_w%bgbl(i)=b
        this_w%cgbl(i)=c
        call DiscontMapping(a,b,c,d,this_r%x,this_r%y,z,etype)
        this_w%apar(i)=a
        this_w%bpar(i)=b
        this_w%cpar(i)=c
    end subroutine SetDiscont

    subroutine SetDiscontFromFile(this_r,this_w,i,coefficients)
        implicit none
        type(Quad), intent(in) :: this_r
        type(Quad), intent(out) :: this_w
        integer :: i
        real(kind=rk) :: a,b,c,d,z(size(this_r%x)),coefficients(2,3)
        this_w%agbl(i)=coefficients(i,1)
        this_w%bgbl(i)=coefficients(i,2)
        this_w%cgbl(i)=coefficients(i,3)
        call DiscontMapping(coefficients(i,1),coefficients(i,2),coefficients(i,3),d,this_r%x,this_r%y,z,etype)
        this_w%apar(i)=coefficients(i,1)
        this_w%bpar(i)=coefficients(i,2)
        this_w%cpar(i)=coefficients(i,3)
    end subroutine SetDiscontFromFile

    subroutine RearrangingCoefficients(this_r,this_w)
        implicit none
        type(Quad), intent(in) :: this_r
        type(Quad), intent(out) :: this_w
        real(kind=rk) atemp,btemp,ctemp
        if ((-this_r%apar(1)/this_r%bpar(1)).gt.(-this_r%apar(2)/this_r%bpar(2))) then
            atemp=this_r%apar(1) ;  btemp=this_r%bpar(1) ;  ctemp=this_r%cpar(1)
            this_w%apar(1)=this_r%apar(2) ;  this_w%bpar(1)=this_r%bpar(2) ;  this_w%cpar(1)=this_r%cpar(2)
            this_w%apar(2)=atemp ;  this_w%bpar(2)=btemp ;  this_w%cpar(2)=ctemp
        end if
        if (this_r%bpar(1).lt.0) then
            this_w%apar(1)=-this_r%apar(1)
            this_w%bpar(1)=-this_r%bpar(1)
            this_w%cpar(1)=-this_r%cpar(1)
        end if
        if (this_r%bpar(2).gt.0) then
            this_w%apar(2)=-this_r%apar(2)
            this_w%bpar(2)=-this_r%bpar(2)
            this_w%cpar(2)=-this_r%cpar(2)
        end if
    end subroutine RearrangingCoefficients

    subroutine SetCoords(this)
        implicit none
        type(Quad), intent(out) :: this
        integer :: i
        real(kind=rk) :: x(nodes),y(nodes)
        do i=1,nodes
            write(*,'(a)',advance='no')'x('
            write(*,'(i1)',advance='no')i
            write(*,'(a)',advance='no')'), y('
            write(*,'(i1)',advance='no')i
            write(*,'(a)',advance='no')'): '
            read(*,*)x(i),y(i)
        end do
        this%x=x
        this%y=y
    end subroutine SetCoords

    subroutine SetCoordsFromFile(this,coords)
        implicit none
        type(Quad), intent(out) :: this
        integer :: i
        real(kind=rk) :: x(nodes),y(nodes),coords(4,2)
        do i=1,nodes
            x(i)=coords(i,1)
            y(i)=coords(i,2)
        end do
        this%x=x
        this%y=y
    end subroutine SetCoordsFromFile

    subroutine SetNumDscnt(this,dscnt)
        implicit none
        type(Quad), intent(out) :: this
        integer :: dscnt
        this%dscnt=dscnt
    end subroutine SetNumDscnt

    subroutine CheckSelection(selection,lower,upper)
        implicit none
        integer :: selection,lower,upper
        do while((selection.lt.lower).or.(selection.gt.upper))
            write(*,'(1X,A)',advance='no') 'Invalid number. Please choose one integer'&
            & // 'from ', lower, ' to ', upper, ' : '
            read (*,*) selection
        end do
    end subroutine CheckSelection

    subroutine Integrate(this_r,this_w)
        implicit none
        type(Quad), intent(in) :: this_r
        type(Quad), intent(out) :: this_w
        integer :: i
        real(kind=rk) :: t
        do i=1,size(this_r%eqcv)
            t=0; this_w%w=0; this_w%w(i)=1
            call GaussQuadrature(this_r,t)
            this_w%out_vec(i)=t                               ! Store n-th integration result
        end do
    end subroutine Integrate

    double precision function detJ(this_r,xi,eta)
        implicit none
        type(Quad), intent(in) :: this_r
        real(kind=rk) :: xi,eta
        detJ=abs((-((1.d0-eta)*this_r%x(1))/4.d0+((1.d0-eta)*this_r%x(2))/4.d0+((1.d0+eta)*this_r%x(3))&
            &/4.d0-((1.d0+eta)*this_r%x(4))/4.d0)*(-((1.d0-xi)*this_r%y(1))/4.d0-((1.d0+xi)*this_r%y(2)&
            &)/4.d0+((1.d0+xi)*this_r%y(3))/4.d0+((1.d0-xi)*this_r%y(4))/4.d0)-(-(this_r%x(1)*(1.d0-xi)&
            &)/4.d0-(this_r%x(2)*(1.d0+xi))/4.d0+(this_r%x(3)*(1.d0+xi))/4.d0+(this_r%x(4)*(1.d0-xi))/4&
            &.d0)*(-((1.d0-eta)*this_r%y(1))/4.d0+((1.d0-eta)*this_r%y(2))/4.d0+((1.d0+eta)*this_r%y(3)&
            &)/4.d0-((1.d0+eta)*this_r%y(4))/4.d0))
        end function detJ

    subroutine EvalN(this_r,this_w)
        implicit none
        type(Quad), intent(in) :: this_r
        type(Quad), intent(out) :: this_w
        integer :: i
        do i=1,size(this_r%gpts(:,1))
            this_w%N(i,1)=0.25*(1.d0-this_r%gpts(i,1))*(1.d0-this_r%gpts(i,2))
            this_w%N(i,2)=0.25*(1.d0+this_r%gpts(i,1))*(1.d0-this_r%gpts(i,2))
            this_w%N(i,3)=0.25*(1.d0+this_r%gpts(i,1))*(1.d0+this_r%gpts(i,2))
            this_w%N(i,4)=0.25*(1.d0-this_r%gpts(i,1))*(1.d0+this_r%gpts(i,2))
        end do
    end subroutine EvalN

    subroutine Mapping(this_r,this_w)
        implicit none
        type(Quad), intent(in) :: this_r
        type(Quad), intent(out) :: this_w
        integer :: i,j
        call EvalN(this_r,this_w)
        do i=1,size(this_r%gpts(:,1))
            do j=1,nodes
                this_w%gbl_gpts(i,1)=this_r%gbl_gpts(i,1)+this_r%x(j)*this_r%N(i,j)
                this_w%gbl_gpts(i,2)=this_r%gbl_gpts(i,2)+this_r%y(j)*this_r%N(i,j)
            end do
        end do
    end subroutine Mapping

    subroutine GaussQuadrature(this,res)
        implicit none
        type(Quad), intent(in) :: this
        integer :: i
        real(kind=rk) :: res
        do i=1,theSize*theSize
            res=res+&                                                           ! Summation
            &HEqPol(this%gpts(i,1),this%gpts(i,2),0.d0,this%eqcv,etype)&        ! Equivalent polynomial
            &*this%gwts(i,1)*this%gwts(i,2)&                                    ! Gauss weights
            &*HEqPol(this%gbl_gpts(i,1),this%gbl_gpts(i,2),0.d0,this%w,etype)&  ! n-th monomial basis function
            &*detJ(this,this%gpts(i,1),this%gpts(i,2))                          ! Jacobian
        end do
    end subroutine GaussQuadrature

    subroutine PrintRes(this,part,filename)
        implicit none
        type(Quad), intent(in) :: this
        integer :: i, ppos
        logical :: itsopen
        character :: part
        character*132 :: filename
        integer, parameter :: fileid=44
        ppos = scan(trim(filename),'.', BACK= .true.)
        if ( ppos > 0 ) filename = filename(1:ppos)//'out'
        inquire(unit=fileid, opened=itsopen)
        if (.not.itsopen) then
            open(fileid,file=filename,action='write')
            write(fileid,*)'****** EQP INTEGRATION RESULTS ******'
            write(*,*)
            write(*,*)'****** EQP INTEGRATION RESULTS ******'
        end if
        if (this%dscnt==1) then
            write(fileid,*)
            write(fileid,*)'Integration of part ' // part // ': '
            write(*,*)
            write(*,*)'Integration of part ' // part // ': '
            do i=1,size(this%out_vec)
                write(fileid,'(I2,5X,F15.7)') i,this%out_vec(i)
                write(*,'(I2,5X,F15.7)') i,this%out_vec(i)
            end do
        else
            write(fileid,*)
            write(fileid,*)'Integration of part ' // part // ': '
            write(*,*)
            write(*,*)'Integration of part ' // part // ': '
            do i=1,size(this%out_vec)
                write(fileid,'(I2,5X,F15.7)') i,this%out_vec(i)
                write(*,'(I2,5X,F15.7)') i,this%out_vec(i)
            end do
        end if
    end subroutine PrintRes

    subroutine EvalQuad(this_r,this_w)
        implicit none
        type(Quad), intent(out) :: this_r
        type(Quad), intent(out) :: this_w
        integer :: i
        character :: portion(4)=(/'A','B','C','D'/)
        character :: part
        character*132 :: filename='output.out'
        real(kind=rk) :: d,z(size(this_r%x))
        call NodeScheme()
        call SetCoords(this_w)
        call HowManyDiscont(this_w)
        do i=1,this_r%dscnt
            call SetDiscont(this_r,this_w,i)
        end do
        if (this_r%dscnt==2)    call RearrangingCoefficients(this_r,this_w)
        call SetGaussPts(this_w)
        call Mapping(this_r,this_w)
        if (this_r%dscnt==1) then
            part='A'
            call Heqpol_coefficients(this_r%apar(1),this_r%bpar(1),this_r%cpar(1),d,this_r%x,&
                                    &this_r%y,z,this_r%eqcv,etype)
            call Integrate(this_r,this_w)
            call PrintRes(this_r,part,filename)
            part='B'
            call Heqpol_coefficients(-this_r%apar(1),-this_r%bpar(1),-this_r%cpar(1),-d,this_r%x,&
                                    &this_r%y,z,this_r%eqcv,etype)
            call Integrate(this_r,this_w)
            call PrintRes(this_r,part,filename)
        else
            do i=1,size(portion)
                call DD_Heqpol_coefficients(this_r%apar,this_r%bpar,this_r%cpar,this_r%eqcv,etype,portion(i))
                call Integrate(this_r,this_w)
                call PrintRes(this_r,portion(i),filename)
            end do
        end if
        write(*,*)
        write(*,*)'Input any key to quit the program: '
        read(*,*)part
    end subroutine EvalQuad

    subroutine EvalQuadFromFile(this_r,this_w,coords,coefficients,dscnt,part,filename)
        implicit none
        type(Quad), intent(out) :: this_r
        type(Quad), intent(out) :: this_w
        integer :: i,dscnt
        character :: portion(4)=(/'A','B','C','D'/)
        character*3 :: part
        character*132 :: filename
        real(kind=rk) :: d,z(size(this_r%x)),coords(4,2),coefficients(2,3)
        call SetCoordsFromFile(this_w,coords)
        call SetNumDscnt(this_w,dscnt)
        do i=1,this_r%dscnt
            call SetDiscontFromFile(this_r,this_w,i,coefficients)
        end do
        if (this_r%dscnt==2)    call RearrangingCoefficients(this_r,this_w)
        call SetGaussPts(this_w)
        call Mapping(this_r,this_w)
        if (this_r%dscnt==1) then
            part='A'
            call Heqpol_coefficients(this_r%apar(1),this_r%bpar(1),this_r%cpar(1),d,this_r%x,&
                                    &this_r%y,z,this_r%eqcv,etype)
            call Integrate(this_r,this_w)
            call PrintRes(this_r,part,filename)
            part='B'
            call Heqpol_coefficients(-this_r%apar(1),-this_r%bpar(1),-this_r%cpar(1),-d,this_r%x,&
                                    &this_r%y,z,this_r%eqcv,etype)
            call Integrate(this_r,this_w)
            call PrintRes(this_r,part,filename)
        else
        select case (part)
            case ('all')
                do i=1,size(portion)
                    call DD_Heqpol_coefficients(this_r%apar,this_r%bpar,this_r%cpar,this_r%eqcv,&
                                                &etype,portion(i))
                    call Integrate(this_r,this_w)
                    call PrintRes(this_r,portion(i),filename)
                end do
            case default
                call DD_Heqpol_coefficients(this_r%apar,this_r%bpar,this_r%cpar,this_r%eqcv,etype,part)
                call Integrate(this_r,this_w)
                call PrintRes(this_r,part,filename)
        end select
        end if
        write(*,*)
        write(*,*)'Input any key to quit the program: '
        read(*,*)part
    end subroutine EvalQuadFromFile

end module class_Quad
