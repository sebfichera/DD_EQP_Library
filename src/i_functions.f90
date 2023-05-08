! ------------------------ I_FUNCTIONS MODULE -------------------------------
!
! BY:
! SEBASTIANO FICHERA,
! sebastiano.fichera@polito.it
!
! Module containing the input functions to read data from file.
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

module i_functions

    integer :: etype,dscnt
    character*3 :: part
    double precision :: coords(4,2)=0.d0,coefficients(2,3)=0.d0

    contains

    subroutine ReadInput(filename,fileid)
        implicit none
        character*132 :: str,filename
        integer :: i,fileid
        open(fileid,file=filename,action='read')
        do
            call readext(str,fileid)
            select case (str)
                case ('$ElementType')
                    call readext(str,fileid)
                    read(str,*) etype
                case ('$Coords')
                    select case (etype)
                        case(20)
                           do i=1,3
                                call readext(str,fileid)
                                read(str,*) coords(i,1),coords(i,2)
                           end do
                        case(21)
                           do i=1,4
                                call readext(str,fileid)
                                read(str,*) coords(i,1),coords(i,2)
                           end do
                    end select
                case ('$NumOfDiscont')
                    call readext(str,fileid)
                    read(str,*) dscnt
                case ('$DiscontCoefficients')
                    do i=1,dscnt
                        call readext(str,fileid)
                        read(str,*) coefficients(i,1),coefficients(i,2),coefficients(i,3)
                    end do
                case ('$ElementPart')
                    call readext(str,fileid)
                    read(str,*) part
                case ('EOF')
                    exit
                case default
                    write(*,*)
                    write(*,*) 'FATAL ERROR'
                    write(*,*) 'Invalid carachter in ' // filename // ' or missing "$..." identifier.'
                    write(*,*) 'Last character read: ' // str
                    read(*,*)
                    stop
            end select
        end do
        close(fileid)
    end subroutine ReadInput

    subroutine readext(str,fileid)
        implicit none
        character*132 :: str
        integer :: fileid
        do
        read(fileid,'(A132)',END=1) str
        str = trim(str)                         ! delete after text blanks
        if (index(str(1:2),'\\').eq.0) exit     ! delete comment lines, starting with: \\
        end do
        return
        1 str='EOF'
    end subroutine readext

    subroutine SetFileName(filename)
        implicit none
        character*132 :: filename
        write(*,'(1X,A)',advance='no') 'Enter the input file name: '
        read(*,*) filename
    end subroutine SetFileName

end module i_functions
