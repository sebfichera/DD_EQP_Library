! DD EQUIVALENT POLYNOMIAL LIBRARY APPLICATION V 1.0
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
! WHEN USING THIS LIBRARY PLEASE ALWAYS CITE:
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
! IN THIS FILE EQUIVALENT POLYNOMIALS ARE USED TO COMPUTE THE VOLUMES
! AND MOMENTS OF INERTIA OF ALL SUBDOMAINS OF A 2D QUADRANGULAR ELEMENT
! CUT BY 2 DISCONTINUITIES.
!
! See the Library Application Note and the above papers for details
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

! ------------------------ MAIN PROGRAM -----------------------------------
program DD_EQP_TEST

    use i_functions
    use class_Quad

    implicit none
    character :: inputfromfile
    integer, parameter :: fileid = 33
    character*132 :: filename
    type(Quad) :: quadri                                        ! Declare a variable of type Quad
    write(*,'(1X,A)',advance='no')'Input from file? (y/n) '
    read(*,*)inputfromfile
    write(*,*)
    if (inputfromfile=='y') then
        call SetFileName(filename)
        call ReadInput(filename,fileid)
        if (etype==21) call EvalQuadFromFile(quadri,quadri,coords,coefficients,dscnt,part,filename)
    else
        call EvalQuad(quadri,quadri)
    end if

end program DD_EQP_TEST
