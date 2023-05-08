! DD EQUIVALENT POLYNOMIAL LIBRARY V 1.0
!
! BY:
! GIULIO VENTURA,
! giulio.ventura@polito.it
! MAURO CORRADO,
! mauro.corrado@polito.it
! GREGORIO MARIGGIO',
! gregorio.mariggio@polito.it
! SEBASTIANO FICHERA,
! sebastiano.fichera@polito.it
!
! www.equivalent-polynomials.net
!
! WHEN USING THIS LIBRARY PLEASE ALWAYS CITE
!
! Ventura G., On the elimination of quadrature subcells for discontinuous
! functions in the extended finite-element method, International Journal
! for Numerical Methods in Engineering 66 (2006) 761-795.
!
! Ventura G., Benvenuti E., Equivalent polynomials for quadrature in
! heaviside function enriched elements, International Journal for
! Numerical Methods in Engineering 102 (2015) 688-710.
!
! G. Mariggio', S. Fichera, M. Corrado, G. Ventura, Eqp - a 2d/3d library
! for integration of polynomials times step function, SoftwareX 12 (2020)
! 100636. doi:https://doi.org/10.1016/j.softx.2020.100636
!
! USAGE
! Call first DD_Heqpol_coefficients for a given element and discontinuity lines.
! Then, at each Gauss point of coords x,y,z, evaluate HEqPol, passing the eqcv
! vector previously computed.
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

double precision function HEqPol(x,y,z,eqcv,etype)

    ! On input
    !
    ! x,y coordinate of the evaluation point in the PARENT element reference system
    ! x,y in [-1,1] for quad
    !
    ! eqcv vector of polynomial coefficients as computed by DDHeqpol_coefficients
    !
    ! etype the element type. It is:
    ! 21  linear quadrilateral, required length for eqcv is 6
    !
    ! On output
    ! the value of the equivalent polynomial at point x,y

    implicit none
    integer etype,l
    double precision x,y,z
    double precision eqcv(:)
    double precision v(size(eqcv))

    l = size(eqcv)
    select case (etype)
    case (20)
        if (l.ne.1) stop 'ERROR: HEqPol -> Invalid lenght for vector eqcv'
            v=(/ 1.d0 /)
    case (21)
        if (l.ne.6) stop 'ERROR: HEqPol -> Invalid lenght for vector eqcv'
            v=(/ 1.d0, x, x**2, y, x*y, y**2 /)
    case default
        stop 'ERROR: HEqPol -> Invalid etype specified'
    case (30)
        if (l.ne.1) stop 'ERROR: HEqPol -> Invalid lenght for vector eqcv'
            v=(/ 1.D0 /)

    case (31)
        if (l.ne.23) stop 'ERROR: HEqPol -> Invalid lenght for vector eqcv'
            v=(/ 1.D0, x, x**2, y, x*y, x**2*y, y**2, x*y**2, x**2*y**2, z, x*z, x**2*z, y*z, x*y*z, &
            x**2*y*z, y**2*z, x*y**2*z, z**2, x*z**2, x**2*z**2, y*z**2, x*y*z**2, y**2*z**2 /)

    end select

    HEqPol=dot_product(v,eqcv)

end function HEqPol

subroutine Heqpol_coefficients(a,b,c,d,x,y,z,eqcv,etype)

    ! On input: a,b,c,d, x,y,z, etype
    ! On output: eqcv
    !
    ! This sub computed the vector of equivalent polynomial coefficients for given
    ! H discontinuity plane coefficients a,b,c,d in the parent element domain
    ! for 3D elements the discontinuity plane has equation a xi + b eta + c zeta + d = 0
    ! for 2D elements the discontinuity plane has equation a xi + b eta + c = 0
    !
    ! On output the vector eqcv of polinomial coefficients w.r.t. the base
    ! defined in the function Hex_HeqPol
    !
    ! etype the element type. It is:
    ! 20  linear triangle, required length for eqcv is 1
    ! 21  linear quadrilateral, required length for eqcv is 6
    ! 30  linear tetrahedron, required length for eqcv is 1
    ! 31  linear hexahedron, required length for eqcv is 23

    implicit none

    interface
        subroutine DiscontMapping(a,b,c,d,x,y,z,etype)
            implicit none
            integer :: etype
            double precision :: a,b,c,d,L1,L2,L3,L4
            double precision,allocatable :: xi(:),eta(:),zeta(:)
            double precision :: x(:),y(:),z(:)
        end subroutine DiscontMapping
    end interface

    integer :: etype
    double precision :: a,b,c,d,x(:),y(:),z(:),eqcv(:)
    double precision :: BV(size(eqcv))
    double precision, parameter :: tol = 1E-4 ! tolerance for vanishing plane coefficients
    double precision, parameter :: Pi = 4 * atan (1.0_8)

    double precision, parameter :: InvA20(1,1) = reshape( (/ 2.d0 /), shape(InvA20))

    double precision, parameter :: InvA21(6,6) = reshape( &
    (/ 7.d0/8.d0, 0.d0, -15.d0/16.d0, 0.d0, 0.d0, -15.d0/16.d0, &
    0.d0, 3.d0/4.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    -15.d0/16.d0, 0.d0, 45.d0/16.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 3.d0/4.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, 9.d0/4.d0, 0.d0, &
    -15.d0/16.d0, 0.d0, 0.d0, 0.d0, 0.d0, 45.d0/16.d0 &
    /),shape(InvA21))

    double precision, parameter :: InvA22(6,6) = reshape( &
    (/ 4.d0/Pi,0.d0,-6.d0/Pi,0.d0,0.d0,-6.d0/Pi, &
    0.d0,4.d0/Pi,0.d0,0.d0,0.d0,0.d0, &
    -6.d0/Pi,0.d0,18.d0/Pi,0.d0,0.d0,6.d0/Pi, &
    0.d0,0.d0,0.d0,4.d0/Pi,0.d0,0.d0, &
    0.d0,0.d0,0.d0,0.d0,24.d0/Pi,0.d0, &
    -6.d0/Pi,0.d0,6.d0/Pi,0.d0,0.d0,18.d0/Pi &
    /),shape(InvA22))

    double precision, parameter :: InvA30(1,1) = reshape( (/ 6.d0 /), shape(InvA30))

    double precision, parameter :: InvA31(23,23) = reshape( &
    (/ 151.d0/128.d0, 0.d0, -(105.d0/64.d0), 0.d0, 0.d0, 0.d0, -(105.d0/64.d0), 0.d0, 225.d0/128.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, -(105.d0/64.d0), 0.d0, 225.d0/128.d0, 0.d0, 0.d0, 225.d0/128.d0, 0.d0, 21.d0/16.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, &
    0.d0, -(105.d0/64.d0), 0.d0, 315.d0/64.d0, 0.d0, 0.d0, 0.d0, 225.d0/128.d0, 0.d0, -(675.d0/128.d0), 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 225.d0/128.d0, 0.d0, -(675.d0/128.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 21.d0/16.d0, &
    0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(45.d0/32.d0), 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 81.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, -(135.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, -(45.d0/32.d0), 0.d0, 135.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(105.d0/64.d0), 0.d0, 225.d0/128.d0, 0.d0, 0.d0, 0.d0, &
    315.d0/64.d0, 0.d0, -(675.d0/128.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 225.d0/128.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, -(675.d0/128.d0), 0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 135.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 225.d0/128.d0, 0.d0, -(675.d0/128.d0), 0.d0, 0.d0, &
    0.d0, -(675.d0/128.d0), 0.d0, 2025.d0/128.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 21.d0/16.d0, 0.d0, -(45.d0/32.d0), 0.d0, 0.d0, &
    0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    81.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(135.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, -(45.d0/32.d0), 0.d0, 135.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 81.d0/32.d0, 0.d0, -(135.d0/32.d0), 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 27.d0/8.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, -(135.d0/32.d0), 0.d0, 405.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 135.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(135.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 405.d0/32.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, -(105.d0/64.d0), 0.d0, 225.d0/128.d0, 0.d0, 0.d0, 0.d0, 225.d0/128.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 315.d0/64.d0, 0.d0, -(675.d0/128.d0), 0.d0, &
    0.d0, -(675.d0/128.d0), 0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 135.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 225.d0/128.d0, 0.d0, -(675.d0/128.d0), 0.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(675.d0/128.d0), 0.d0, 2025.d0/128.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    135.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(135.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 405.d0/32.d0, 0.d0, 225.d0/128.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(675.d0/128.d0), &
    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(675.d0/128.d0), 0.d0, 0.d0, 0.d0, 0.d0, 2025.d0/128.d0 &
    /),shape(InvA31))

    integer i
    double precision a2,a3,a4,a5,a6,b2,b3,b4,b5,b6,c2,c3,c4,c5,c6,d2,d3,d4,d5,d6
    double precision t,am,bm,cm
    double precision abs01,abs02,abs03,abs04,abs05,abs06,abs07,abs08,abs09,abs10,abs11,abs12
    double precision abs13,abs14,abs15,abs16

    ! plane coefficient normalization
    select case (etype)

    case (20,21)
    !call DiscontMapping(a,b,c,d,x,y,z,etype) discontinuity is already mapped by element_Class
    t=sqrt(a**2+b**2)
    a=a/t; b=b/t; c=c/t
    am=abs(a)
    bm=abs(b)

    case (30,31)
    call DiscontMapping(a,b,c,d,x,y,z,etype)
    t=sqrt(a**2+b**2+c**2)
    a=a/t; b=b/t; c=c/t; d=d/t
    am=abs(a)
    bm=abs(b)
    cm=abs(c)

    case default
    stop 'Invalid etype'

    end select

    select case (etype)

    case (20)
    a2=a**2; b2=b**2
    if ( am<bm*tol ) then
        ! a=0, b<>0
        ! b
        BV(1) = (b2-(2*b+c)*abs(c)+(b+c)*abs(b+c))/(4.d0*b2)
    else if ( bm<am*tol ) then
        ! a<>0, b=0
        ! a
        BV(1) = (a2-(2*a+c)*abs(c)+(a+c)*abs(a+c))/(4.d0*a2)
    else if (abs(a-b)<tol) then
        ! a<>0, b<>0, a=b
        ! a=b
        BV(1) = (a2+c*abs(c)+(a-c)*abs(a+c))/(4.d0*a2)
    else
        ! a<>0, b<>0
        ! a b
        BV(1) = ((a-b)*c*abs(c)+b*(a+c)*abs(a+c)+a*((a-b)*b-(b+c)*abs(b+c)))/(4.d0*a*(a-b)*b)
    end if
    eqcv=matmul(InvA20,BV)

    case (21)
    a2=a**2; b2=b**2; c2=c**2
    a3=a**3; b3=b**3; c3=c**3
    if ( am<bm*tol ) then
        ! a=0, b<>0
        ! b
        abs01=abs(b-c)
        abs02=abs(b+c)
        include 'quad_ph_6_b.fi'
    else if ( bm<am*tol ) then
        ! a<>0, b=0
        ! a
        abs01=abs(a-c)
        abs02=abs(a+c)
        include 'quad_ph_6_a.fi'
    else
        ! a<>0, b<>0
        ! a b
        abs01=abs(a-b-c)
        abs02=abs(a+b-c)
        abs03=abs(a-b+c)
        abs04=abs(a+b+c)
        include 'quad_ph_6_ab.fi'
    end if
    eqcv=matmul(InvA21,BV)

    case (30)
    a2=a**2; b2=b**2; c2=c**2; d2=d**2
    a3=a**3; b3=b**3; c3=c**3; d3=d**3
    if ( am<bm*tol .and. am<cm*tol .and. min(bm,cm)/max(bm,cm)>tol ) then
        ! a=0, b<>0, c<>0
        ! b c
        if (abs(b-c)>tol) then
        ! b<>c
        BV(1)=((b-c)*d*(c*d+b*(3*c+d))*abs(d)+c2*(b + d)**2*abs(b+d)+ &
                b2*((b-c)*c2-(c+d)**2*abs(c+d)))/(12.d0*b2*(b-c)*c2)
        else
        ! b=c
        BV(1) = (b3 + d*(3*b + 2*d)*abs(d) + (b2 - b*d - 2*d2)*Abs(b + d))/(12.*b3)
        end if
    else if ( bm<am*tol .and. bm<cm*tol .and. min(am,cm)/max(am,cm)>tol ) then
        ! a<>0, b=0, c<>0
        ! a c
        if (abs(a-c)>tol) then
        ! a<>c
        BV(1)=((a-c)*d*(c*d+a*(3*c+d))*abs(d)+c2*(a+d)**2*abs(a+d)+ &
                a2*((a-c)*c2-(c+d)**2*abs(c+d)))/(12.d0*a2*(a - c)*c2)
        else
        ! a=c
        BV(1) = (a3 + d*(3*a + 2*d)*abs(d) + (a2 - a*d - 2*d2)*abs(a + d))/(12.*a3)
        end if
    else if ( cm<am*tol .and. cm<bm*tol .and. min(am,bm)/max(am,bm)>tol ) then
        ! a<>0, b<>0, c=0
        ! a b
        if (abs(a-b)>tol) then
        ! a<>b
        BV(1)=((a-b)*d*(b*d+a*(3*b+d))*abs(d)+b2*(a+d)**2*abs(a+d)+ &
                a2*((a-b)*b2-(b+d)**2*abs(b+d)))/(12.d0*a2*(a-b)*b2)
        else
        ! a=b
        BV(1) = (a3 + d*(3*a + 2*d)*abs(d) + (a2 - a*d - 2*d2)*abs(a + d))/(12.*a3)
        end if
    else if ( am<cm*tol .and. bm<cm*tol ) then
        ! a=0, b=0, c<>0
        ! c
        BV(1)=(c3-(3*c2+3*c*d+d2)*abs(d)+(c+d)**2*abs(c+d))/(12.d0*c3)
    else if ( bm<am*tol .and. cm<am*tol ) then
        ! a<>0, b=0, c=0
        ! a
        BV(1)=(a3-(3*a2+3*a*d+d2)*abs(d)+(a+d)**2*abs(a+d))/(12.d0*a3)
    else if ( am<bm*tol .and. cm<bm*tol ) then
        ! a=0, b<>0, c=0
        ! b
        BV(1)=(b3-(3*b2+3*b*d+d2)*abs(d)+(b+d)**2*abs(b+d))/(12.d0*b3)
    else
        ! a<>0, b<>0, c<>0
        ! a b c
        if (abs(a-b)>tol .and. abs(a-c)>tol .and. abs(b-c)>tol) then
        ! a<>b, a<>c, b<>c
        BV(1)=((-a+b)*(a-c)*(b-c)*d2*abs(d)+b*(b-c)*c*(a+d)**2*abs(a+d)+a*(c*(-a+c)*(b+d)**2*abs(b+d)+ &
                (a-b)*b*((a-c)*(b-c)*c+(c+d)**2*abs(c+d))))/(12.d0*a*(a-b)*b*(a-c)*(b-c)*c)
        else if (abs(a-b)<tol .and. abs(a-c)>tol .and. abs(b-c)>tol) then
        ! a=b, a<>c, b<>c
        BV(1)=(-((a - c)**2*d2*abs(d)) + c*(a + d)*(a2 + c*d - 2*a*(c + d))*abs(a + d) + &
                a2*((a - c)**2*c + (c + d)**2*abs(c + d)))/(12.*a2*(a - c)**2*c)
        else if (abs(a-b)>tol .and. abs(a-c)<tol .and. abs(b-c)>tol) then
        ! a<>b, a=c, b<>c
        BV(1)=(-((a - b)**2*d2*abs(d)) + b*(a + d)*(a2 + b*d - 2*a*(b + d))*abs(a + d) + &
                a2*((a - b)**2*b + (b + d)**2*abs(b + d)))/(12.*a2*(a - b)**2*b)
        else if (abs(a-b)>tol .and. abs(a-c)>tol .and. abs(b-c)<tol) then
        ! a<>b, a<>c, b=c
        BV(1)=(-((a - b)**2*d2*abs(d)) + b2*(a + d)**2*abs(a + d) + &
                a*((a - b)**2*b2 + (b + d)*(b*(b - 2*d) + a*(-2*b + d))*abs(b + d)))/ &
                (12.*a*(a - b)**2*b2)
        end if
    end if
    eqcv=matmul(InvA30,BV)

    case (31)

    a2=a**2; b2=b**2; c2=c**2; d2=d**2
    a3=a**3; b3=b**3; c3=c**3; d3=d**3
    a4=a**4; b4=b**4; c4=c**4; d4=d**4
    a5=a**5; b5=b**5; c5=c**5; d5=d**5
    a6=a**6; b6=b**6; c6=c**6; d6=d**6
    if ( am<bm*tol .and. am<cm*tol .and. min(bm,cm)/max(bm,cm)>tol ) then
        ! a=0, b<>0, c<>0
        ! b c
        abs01=abs(b+c+d)
        abs02=abs(b+c-d)
        abs03=abs(b-c+d)
        abs04=abs(b-c-d)
        abs05=abs(-b+c+d)
        abs06=abs(-b+c-d)
        abs07=abs(-b-c+d)
        abs08=abs(-b-c-d)
        include 'hex_ph_23_bc.fi'
    else if ( bm<am*tol .and. bm<cm*tol .and. min(am,cm)/max(am,cm)>tol ) then
        ! a<>0, b=0, c<>0
        ! a c
        abs01=abs(a+c+d)
        abs02=abs(a+c-d)
        abs03=abs(a-c+d)
        abs04=abs(a-c-d)
        abs09=abs(-a+c+d)
        abs10=abs(-a+c-d)
        abs11=abs(-a-c+d)
        abs12=abs(-a-c-d)
        include 'hex_ph_23_ac.fi'
    else if ( cm<am*tol .and. cm<bm*tol .and. min(am,bm)/max(am,bm)>tol ) then
        ! a<>0, b<>0, c=0
        ! a b
        abs01=abs(a+b+d)
        abs02=abs(a+b-d)
        abs05=abs(a-b+d)
        abs06=abs(a-b-d)
        abs09=abs(-a+b+d)
        abs10=abs(-a+b-d)
        abs13=abs(-a-b+d)
        abs14=abs(-a-b-d)
        include 'hex_ph_23_ab.fi'
    else if ( am<cm*tol .and. bm<cm*tol ) then
        ! a=0, b=0, c<>0
        ! c
        abs01=abs(c+d)
        abs02=abs(c-d)
        abs03=abs(-c+d)
        abs04=abs(-c-d)
        include 'hex_ph_23_c.fi'
    else if ( bm<am*tol .and. cm<am*tol ) then
        ! a<>0, b=0, c=0
        ! a
        abs01=abs(a+d)
        abs02=abs(a-d)
        abs09=abs(-a+d)
        abs10=abs(-a-d)
        include 'hex_ph_23_a.fi'
    else if ( am<bm*tol .and. cm<bm*tol ) then
        ! a=0, b<>0, c=0
        ! b
        abs01=abs(b+d)
        abs02=abs(b-d)
        abs05=abs(-b+d)
        abs06=abs(-b-d)
        include 'hex_ph_23_b.fi'
    else
        ! a<>0, b<>0, c<>0
        ! a b c
        abs01=abs(a+b+c+d)
        abs02=abs(a+b+c-d)
        abs03=abs(a+b-c+d)
        abs04=abs(a+b-c-d)
        abs05=abs(a-b+c+d)
        abs06=abs(a-b+c-d)
        abs07=abs(a-b-c+d)
        abs08=abs(a-b-c-d)
        abs09=abs(-a+b+c+d)
        abs10=abs(-a+b+c-d)
        abs11=abs(-a+b-c+d)
        abs12=abs(-a+b-c-d)
        abs13=abs(-a-b+c+d)
        abs14=abs(-a-b+c-d)
        abs15=abs(-a-b-c+d)
        abs16=abs(-a-b-c-d)
        include 'hex_ph_23_abc.fi'
    end if
    eqcv=matmul(InvA31,BV)

    case default
    stop 'Invalid etype specified'

    end select

end subroutine Heqpol_coefficients

subroutine DD_Heqpol_coefficients(a,b,c,eqcv,etype,part)

    ! On input: a1,b1,c1, a2,b2,c2, etype
    ! On output: eqcv
    !
    ! This sub computed the vector of equivalent polynomial coefficients for given
    ! double discontinuity plane coefficients a1,b1,c1 a2,b2,c2 in the parent element domain
    ! for 2D elements the discontinuity plane has equation a xi + b eta + c = 0
    !
    ! On output the vector eqcv of polinomial coefficients w.r.t. the base
    ! defined in the function Hex_HeqPol
    !
    ! etype the element type. It is:
    ! 21  linear quadrilateral, required length for eqcv is 6

    implicit none
    integer :: etype
    character :: part
    logical :: s_x_out
    double precision :: a(2),b(2),c(2),s_vec(2),eqcv(:)
    double precision :: a1,b1,c1,a2,b2,c2,s
    double precision :: BVa(size(eqcv)),BVb(size(eqcv)),BVc(size(eqcv)),BVd(size(eqcv))
    real, parameter :: tol = 1E-4 ! tolerance for vanishing plane coefficients
    double precision, parameter :: InvA21(6,6) = reshape( &
    (/ 7.d0/8.d0, 0.d0, -15.d0/16.d0, 0.d0, 0.d0, -15.d0/16.d0, &
    0.d0, 3.d0/4.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
    -15.d0/16.d0, 0.d0, 45.d0/16.d0, 0.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 3.d0/4.d0, 0.d0, 0.d0, &
    0.d0, 0.d0, 0.d0, 0.d0, 9.d0/4.d0, 0.d0, &
    -15.d0/16.d0, 0.d0, 0.d0, 0.d0, 0.d0, 45.d0/16.d0 &
    /),shape(InvA21))
    integer i
    double precision t1,t2,am1,am2,bm1,bm2
    call evalS(a,b,c,s_vec) ! evaluating the intersection point position
    if (s_vec(1).ge.1) then
      s_x_out=.true. ; s=1
    else if (s_vec(1).le.-1) then
      s_x_out=.true. ; s=-1
    else
      s_x_out=.false. ; s=s_vec(1)
    end if
    ! set the function variables
    a1=a(1) ; b1=b(1) ; c1=c(1)
    a2=a(2) ; b2=b(2) ; c2=c(2)
    ! plane coefficient normalization
    select case (etype)
        case (21)
        t1=sqrt(a1**2+b1**2)
        t2=sqrt(a2**2+b2**2)
        am1=abs(a1/t1)
        am2=abs(a2/t2)
        bm1=abs(b1/t1)
        bm2=abs(b2/t2)
        case default
        stop 'ERROR: DD_Heqpol_coefficients -> Invalid etype'
    end select
    select case (etype)
        case (21)
        if (am1<bm1*tol.and.&       ! a1==0, b1<>0
            &am2>tol.and.&          ! a2<>0
            &bm2>tol&               ! b2<>0
            &) then
            ! a1==0, b1<>0 && a2<>0, b2<>0
            ! b1 a2 b2
            if (.not.s_x_out) then
              include 'dd_quad_b1_a2b2.fi'
            else
              include 'dd_quad_b1_a2b2_s_out.fi'
            end if
        else if (bm1<am1*tol.and.&      ! b1==0, a1<>0
                &am2>tol.and.&          ! a2<>0
                &bm2>tol&               ! b2<>0
                &) then
            ! a1<>0, b1==0 && a2<>0, b2<>0
            ! a1 a2 b2
            if (.not.s_x_out) then
              include 'dd_quad_a1_a2b2.fi'
            else
              include 'dd_quad_a1_a2b2_s_out.fi'
            end if
        else if (am2<bm2*tol.and.&      ! a2==0, b2<>0
                &am1>tol.and.&          ! a1<>0
                &bm1>tol&               ! b1<>0
                &) then
            ! a1<>0, b1<>0 && a2==0, b2<>0
            ! a1 b1 b2
            if (.not.s_x_out) then
              include 'dd_quad_a1b1_b2.fi'
            else
              include 'dd_quad_a1b1_b2_s_out.fi'
            end if
        else if (am1<bm1*tol.and.&      ! a1==0, b1<>0
                &am2<bm2*tol&           ! a2==0, b2<>0
                &) then
            ! a1==0, b1<>0 && a2==0, b2<>0
            ! b1 b2
            include 'dd_quad_b1_b2.fi'
        else if (bm1<am1*tol.and.&      ! a1<>0, b1==0
                &am2<bm2*tol&           ! a2==0, b2<>0
                &) then
            ! a1<>0, b1==0 && a2==0, b2<>0
            ! a1 b2
            if (.not.s_x_out) then
              include 'dd_quad_a1_b2.fi'
            else
              include 'dd_quad_a1_b2_s_out.fi'
            end if
        else if (bm2<am2*tol.and.&      ! a2<>0, b2==0
                &am1>tol.and.&          ! a1<>0
                &bm1>tol&               ! b1<>0
                &) then
            ! a1<>0, b1<>0 && a2<>0, b2==0
            ! a1 b1 a2
            if (.not.s_x_out) then
              include 'dd_quad_a1b1_a2.fi'
            else
              include 'dd_quad_a1b1_a2_s_out.fi'
            end if
        else if (am1<bm1*tol.and.&      ! a1==0, b1<>0
                &bm2<am2*tol&           ! a2<>0, b2==0
                &) then
            ! a1==0, b1<>0 && a2<>0, b2==0
            ! b1 a2
            if (.not.s_x_out) then
              include 'dd_quad_b1_a2.fi'
            else
              include 'dd_quad_b1_a2_s_out.fi'
            end if
        else if (bm1<am1*tol.and.&      ! a1<>0, b1==0
                &bm2<am2*tol&           ! a2<>0, b2==0
                &) then
            ! a1<>0, b1==0 && a2<>0, b2==0
            ! a1 a2
            include 'dd_quad_a1_a2.fi'
        else
            ! a1<>0, b1<>0 && a2<>0, b2<>0
            ! a1 b1 a2 b2
            if (.not.s_x_out) then
              include 'dd_quad_a1b1_a2b2.fi'
            else
              include 'dd_quad_a1b1_a2b2_s_out.fi'
            end if
        end if
        select case (part)
            case('A')
                eqcv=matmul(InvA21,BVa)
            case('B')
                eqcv=matmul(InvA21,BVb)
            case('C')
                eqcv=matmul(InvA21,BVc)
            case('D')
                eqcv=matmul(InvA21,BVd)
            case default
                stop 'ERROR: DD_Heqpol_coefficients -> Invalid part specified'
        end select
    case default
        stop 'ERROR: DD_Heqpol_coefficients -> Invalid etype specified'
    end select

end subroutine DD_Heqpol_coefficients

module distance

    implicit none
    interface DistanceFromLine
        module procedure DistanceFromLine_2D
        module procedure DistanceFromLine_3D
    end interface
    contains

    ! Distance between a line and a point
    double precision function DistanceFromLine_2D(x,y,a,b,c)

    ! On input
    !
    ! x,y coordinates of a generic point
    ! a,b,c coefficients of a generic line
    !
    ! On output
    ! Distance between the line and the point in a 2D space

    implicit none
    double precision :: x,y,a,b,c

    DistanceFromLine_2D = (a*x+b*y+c)/sqrt(a**2.d0+b**2.d0)

    end function DistanceFromLine_2D

    ! Distance between a plane and a point in a 3D space
    double precision function DistanceFromLine_3D(x,y,z,a,b,c,d)

    ! On input
    !
    ! x,y,z coordinates of a generic point
    ! a,b,c,d coefficients of a generic plane
    !
    ! On output
    ! Distance between the plane and the point in a 3D space

    implicit none
    double precision :: x,y,z,a,b,c,d

    DistanceFromLine_3D = (a*x+b*y+c*z+d)/sqrt(a**2.d0+b**2.d0+c**2.d0)

    end function DistanceFromLine_3D

end module distance

subroutine DiscontMapping(a,b,c,d,x,y,z,etype)

    ! On input
    !
    ! x,y,z vectors containing the element nodes in the GLOBAL coordinate system
    ! a,b,c,d discontinuity coefficients in the GLOBAL coordinate system
    !
    ! etype the element type. It is:
    ! 20  linear triangle
    ! 21  linear quadrilateral
    ! 30  linear tetrahedron
    ! 31  linear hexahedron, required length for eqcv is 23
    !
    ! On output
    ! a,b,c,d discontinuity coefficients in the PARENT coordinate system

    use distance

    implicit none
    integer :: etype
    double precision :: a,b,c,d,L1,L2,L3,L4
    double precision,allocatable :: xi(:),eta(:),zeta(:)
    double precision :: x(:),y(:),z(:)

    if (etype .lt. 30) then
    allocate (zeta(1))
    end if
    select case (etype)
        case(20)
            allocate(xi(3),eta(3))
            d=0.d0; L4=0.d0; zeta=0.d0; z=0.d0
            ! Initialization of the parent system
            xi      = (/0.d0, 1.d0, 0.d0/)
            eta     = (/0.d0, 0.d0, 1.d0/)
            ! Evaluating levelset values L1,L2,L3 for the global system
            L1 = DistanceFromLine(x(1),y(1),a,b,c)
            L2 = DistanceFromLine(x(2),y(2),a,b,c)
            L3 = DistanceFromLine(x(3),y(3),a,b,c)
            ! Compute a,b,c coefficients for the parent system
            a = -((L2*eta(1) - L3*eta(1) - L1*eta(2) + L3*eta(2) + L1*eta(3) &
            - L2*eta(3))/(-(xi(2)*eta(1)) + xi(3)*eta(1) + xi(1)*eta(2) - xi(3)*eta(2) - xi(1)*eta(3) + xi(2)*eta(3)))

            b = -((L2*xi(1) - L3*xi(1) - L1*xi(2) + L3*xi(2) + L1*xi(3) &
            - L2*xi(3))/(xi(2)*eta(1) - xi(3)*eta(1) - xi(1)*eta(2) + xi(3)*eta(2) + xi(1)*eta(3) - xi(2)*eta(3)))

            c = -((-(L3*xi(2)*eta(1)) + L2*xi(3)*eta(1) + L3*xi(1)*eta(2) &
            - L1*xi(3)*eta(2) - L2*xi(1)*eta(3) + L1*xi(2)*eta(3))/(xi(2)*eta(1) - xi(3)*eta(1) - xi(1)*eta(2) + xi(3)*eta(2) + xi(1)*eta(3) - xi(2)*eta(3)))
        case(21)
            allocate(xi(4),eta(4))
            d=0; L4=0
            ! Initialization of the parent system
            xi      = (/-1.d0,  1.d0, 1.d0, -1.d0/)
            eta     = (/-1.d0, -1.d0, 1.d0,  1.d0/)
            ! Evaluating levelset values L1,L2,L3 for the global system
            L1 = DistanceFromLine(x(1),y(1),a,b,c)
            L2 = DistanceFromLine(x(2),y(2),a,b,c)
            L3 = DistanceFromLine(x(3),y(3),a,b,c)
            ! Compute a,b,c coefficients for the parent system
            a = -((L2*eta(1) - L3*eta(1) - L1*eta(2) + L3*eta(2) + L1*eta(3) &
            - L2*eta(3))/(-(xi(2)*eta(1)) + xi(3)*eta(1) + xi(1)*eta(2) - xi(3)*eta(2) - xi(1)*eta(3) + xi(2)*eta(3)))

            b = -((L2*xi(1) - L3*xi(1) - L1*xi(2) + L3*xi(2) + L1*xi(3) &
            - L2*xi(3))/(xi(2)*eta(1) - xi(3)*eta(1) - xi(1)*eta(2) + xi(3)*eta(2) + xi(1)*eta(3) - xi(2)*eta(3)))

            c = -((-(L3*xi(2)*eta(1)) + L2*xi(3)*eta(1) + L3*xi(1)*eta(2) &
            - L1*xi(3)*eta(2) - L2*xi(1)*eta(3) + L1*xi(2)*eta(3))/(xi(2)*eta(1) - xi(3)*eta(1) - xi(1)*eta(2) + xi(3)*eta(2) + xi(1)*eta(3) - xi(2)*eta(3)))
    end select

end subroutine DiscontMapping

subroutine evalS(a,b,c,s)

  ! On input: a1,b1,c1, a2,b2,c2
  !
  ! This fun evaluates the position of the discontinuities intersection point S
  !
  ! On output the position of S w.r.t. the X-axis and Y-axis

  implicit none
  double precision :: a(2),b(2),c(2),s(2),am1,bm1,t1,am2,bm2,t2
  double precision, parameter :: tol=1e-4

  t1=sqrt(a(1)**2+b(1)**2)
  t2=sqrt(a(2)**2+b(2)**2)
  am1=abs(a(1)/t1)
  am2=abs(a(2)/t2)
  bm1=abs(b(1)/t1)
  bm2=abs(b(2)/t2)

  if (am1<bm1*tol.and.&       ! a1==0, b1<>0
      &am2>tol.and.&          ! a2<>0
      &bm2>tol&               ! b2<>0
      &) then
      ! a1==0, b1<>0 && a2<>0, b2<>0
      ! b1 a2 b2
      s(1)=-((-b(2)*c(1)+b(1)*c(2))/(a(2)*b(1)))
      s(2)=-(c(1)/b(1))
  else if (bm1<am1*tol.and.&      ! b1==0, a1<>0
          &am2>tol.and.&          ! a2<>0
          &bm2>tol&               ! b2<>0
          &) then
      ! a1<>0, b1==0 && a2<>0, b2<>0
      ! a1 a2 b2
      s(1)=-(c(1)/a(1))
      s(2)=-((-a(2)*c(1)+a(1)*c(2))/(a(1)*b(2)))
  else if (am2<bm2*tol.and.&      ! a2==0, b2<>0
          &am1>tol.and.&          ! a1<>0
          &bm1>tol&               ! b1<>0
          &) then
      ! a1<>0, b1<>0 && a2==0, b2<>0
      ! a1 b1 b2
      s(1)=-((b(2)*c(1)-b(1)*c(2))/(a(1)*b(2)))
      s(2)=-(c(2)/b(2))
  else if (am1<bm1*tol.and.&      ! a1==0, b1<>0
          &am2<bm2*tol&           ! a2==0, b2<>0
          &) then
      ! a1==0, b1<>0 && a2==0, b2<>0
      ! b1 b2
      s(1)=-99999
      s(2)=-99999
  else if (bm1<am1*tol.and.&      ! a1<>0, b1==0
          &am2<bm2*tol&           ! a2==0, b2<>0
          &) then
      ! a1<>0, b1==0 && a2==0, b2<>0
      ! a1 b2
      s(1)=-(c(1)/a(1))
      s(2)=-(c(2)/b(2))
  else if (bm2<am2*tol.and.&      ! a2<>0, b2==0
          &am1>tol.and.&          ! a1<>0
          &bm1>tol&               ! b1<>0
          &) then
      ! a1<>0, b1<>0 && a2<>0, b2==0
      ! a1 b1 a2
      s(1)=-(c(2)/a(2))
      s(2)=-((a(2)*c(1)-a(1)*c(2))/(a(2)*b(1)))
  else if (am1<bm1*tol.and.&      ! a1==0, b1<>0
          &bm2<am2*tol&           ! a2<>0, b2==0
          &) then
      ! a1==0, b1<>0 && a2<>0, b2==0
      ! b1 a2
      s(1)=-(c(2)/a(2))
      s(2)=-(c(1)/b(1))
  else if (bm1<am1*tol.and.&      ! a1<>0, b1==0
          &bm2<am2*tol&           ! a2<>0, b2==0
          &) then
      ! a1<>0, b1==0 && a2<>0, b2==0
      ! a1 a2
      s(1)=-99999
      s(2)=-99999
  else
      ! a1<>0, b1<>0 && a2<>0, b2<>0
      ! a1 b1 a2 b2
      s(1)=-((-b(2)*c(1)+b(1)*c(2))/(a(2)*b(1)-a(1)*b(2)))
      s(2)=-((a(2)*c(1)-a(1)*c(2))/(a(2)*b(1)-a(1)*b(2)))
  end if

end subroutine evalS
