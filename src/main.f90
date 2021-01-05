! Program: nuclear_energies
! By: Austin Keller
!-----------------------------------------------------------------------------
! In this program will read the data and calculate theoretical parameters to c_parameters
! a theoretical model. Which we will then use to calculate the binding energy 
! using the (semi-empirical mass formula), and predict the valley of stability, 
! and neutron drip for given isotopes
!-----------------------------------------------------------------------------
program nuclear_energies
use types
use read_write, only : read_exp_data, write_predictions, write_nuclear_lines
use nuclear_model, only : find_best_parameters
implicit none



!------------------------------------------------------------

integer, allocatable :: n_protons(:), n_neutrons(:)
real(dp), allocatable :: exp_values(:), uncertainties(:), c_parameters(:), covariance(:,:)

call read_exp_data(n_protons, n_neutrons, exp_values, uncertainties)
call find_best_parameters(n_protons, n_neutrons, exp_values, uncertainties, c_parameters, covariance)
call write_predictions(exp_values, uncertainties, c_parameters, covariance, n_protons, n_neutrons)

!------------------------------------------------------------ 
call write_nuclear_lines(c_parameters)

end program nuclear_energies