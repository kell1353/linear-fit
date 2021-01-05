!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! The subroutines in this module are for reading the users input file
!! and checking if the file exist then writing the data to arrays. (first subroutine) 
!!
!! Then writing the results of num_protons', 'num_neutrons', 'exp_BE', 'exp_error', 
!! 'theoretical_BE, theoretical_error' into results.dat for analysis in Jupyter. 
!! (second subroutine)
!!
!! Then writing the results of 'num_protons', 'pos_neutron_stable_isotopes', 'pos_neutron_drip'
!! into results_advanced.dat for analysis in Jupyter. (third subroutine)
!!
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_exp_data
!! write_predictions
!! write_nuclear_lines
!!
!!----------------------------------------------------------------------
module read_write

use types
use nuclear_model, only : semi_empirical_mass, semi_empirical_error, most_bound_nuclei, &
                        neutron_drip_position

implicit none

private
public :: read_exp_data, write_predictions, write_nuclear_lines

contains

!-----------------------------------------------------------------------
!! Subroutine: read_exp_data
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine requests a data file checks if the file exists then 
!! skips the first three rows to read the values from the data file into 
!! specified arrays which will be returned and used throughout the program.
!!
!!----------------------------------------------------------------------
!! Output:
!!
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!-----------------------------------------------------------------------
subroutine read_exp_data(n_protons, n_neutrons, exp_values, uncertainties)
    implicit none
    integer, intent(out), allocatable :: n_protons(:), n_neutrons(:)
    real(dp), intent(out), allocatable :: exp_values(:), uncertainties(:)

    character(len=128) :: filename, tmp !Temporary array to hold columns we don't need
    logical :: file_exists
    integer :: file_unit, num_binding_energies, n, i, ierror, row_skip

    print *, 'This program will read the data and calculate theoretical parameters'
    print *, 'that we will then use in our theoretical model to predict'
    print *, 'the semi-empirical mass, valley of stability, and neutron drip for given isotopes'

    print *, '...'
    print *, 'please provide the file name with the experimental data'
    read(*, '(a)') filename
    
    if(allocated(n_protons)) deallocate(n_protons)
    if(allocated(n_neutrons)) deallocate(n_neutrons)
    if(allocated(exp_values)) deallocate(exp_values)
    if(allocated(uncertainties)) deallocate(uncertainties)

    inquire(file=trim(filename), exist=file_exists)

    if (file_exists) then
        open(newunit=file_unit, file=filename)

        ! read the number of data from the first line in the file
        read(file_unit, *) n 
        row_skip = 3 !Number of rows we intend to skip when reading the data

        ! Allocate the arrays to be appropriate size
        allocate(n_protons(1:n))
        allocate(n_neutrons(1:n))
        allocate(exp_values(1:n))
        allocate(uncertainties(1:n))

        ! Start at row two because we already read the first row.
        do i=2,n+row_skip
            ! Skip the next two rows
            if (i==2 .or. i==3) then
                read(file_unit,*)
            else
                read(file_unit,*) tmp, tmp, tmp, n_neutrons(i-row_skip), &
                n_protons(i-row_skip), exp_values(i-row_skip), tmp, uncertainties(i-row_skip)
            endif
        enddo

    else
        print *, 'The file you have entered cannot be found. Please make sure you have the correct file name and try again.'
        stop
    endif
end subroutine read_exp_data

!-----------------------------------------------------------------------
!! Subroutine: write_predictions
!-----------------------------------------------------------------------
!! Austin Keller 
!!
!! This subroutine takes in the arrays contained in the users dataset, and 
!! theoretical parameters, and covariance we have found using the dataset. 
!!
!! Sets a do loop through the data and for each proton nuetron combination. It
!! calculates binding energy and error using our theoretical semi empirical
!! mass formula. 
!!
!! Then writes all the the results to a 'results.dat' file to be evaluated
!! using Jupyter Notebook.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!! c_parameters     real        Array containing the parameters of the semi-empirical mass formula
!! covariance       real        Array containing the elements of the covariance matrix
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!-----------------------------------------------------------------------
subroutine write_predictions(exp_values, uncertainties, c_parameters, covariance, n_protons, n_neutrons)
    implicit none
    real(dp), intent(in) :: exp_values(:), uncertainties(:), c_parameters(:), covariance(:,:)
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    character(len=*), parameter :: file_name = 'results.dat'
    real(dp) :: theoretical_BE, theoretical_error
    integer :: unit, i, n

    n = size(exp_values)

    open(newunit=unit,file=file_name)
    write(unit,'(5a28)') 'num_protons', 'num_neutrons', 'exp_BE', 'exp_error', 'theoretical_BE, theoretical_error'

    do i=1,n
        theoretical_BE = semi_empirical_mass(c_parameters, n_protons(i), n_neutrons(i))
        theoretical_error = semi_empirical_error(covariance, n_protons(i), n_neutrons(i))
        write(unit,*) n_protons(i), n_neutrons(i), exp_values(i), uncertainties(i), theoretical_BE, theoretical_error
    enddo
    close(unit)

    print *, ''
    print *, 'theoretical binding energies were written in ', file_name
end subroutine write_predictions

!--------------------------------------------------------------------------------


!-----------------------------------------------------------------------
!! Subroutine: write_nuclear_lines
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine takes in the 
!! theoretical parameters we have found using the dataset. 
!!
!! Sets a do loop through the unique values of protons in the dataset. It
!! calculates most stable isotope, calculates calculate the position of the 
!! neutron drip-line.
!!
!! Then writes all the the results (# of neutrons for a given # of protons) 
!! to a 'results_advanced.dat' file to be evaluated using Jupyter Notebook.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters     real        Array containing the parameters of the semi-empirical mass formula
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!-----------------------------------------------------------------------
subroutine write_nuclear_lines(c_parameters)
    implicit none
    real(dp), intent(in) :: c_parameters(:)
    character(len=*), parameter :: file_name = 'results_advanced.dat'
    integer :: unit, i, n, stable_iso_position, nuetron_drip_pos


    n = 118 !Amount of unique values for Z (# of protons)

    open(newunit=unit,file=file_name)
    write(unit,'(5a28)') 'num_protons', 'pos_neutron_stable_isotopes', 'pos_neutron_drip'

    stable_iso_position = 0
    nuetron_drip_pos = 0

    do i=1,n
        stable_iso_position = most_bound_nuclei(i, c_parameters)
        nuetron_drip_pos = neutron_drip_position(i, c_parameters)

        write(unit,*) i, stable_iso_position, nuetron_drip_pos
    enddo
    close(unit)

    print *, 'theoretical nuclear solutions were written in ', file_name
end subroutine write_nuclear_lines
    
end module read_write