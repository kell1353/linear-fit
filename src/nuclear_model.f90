!-----------------------------------------------------------------------
!Module: nuclear_model
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! The purpose of the functions and subroutines in this module is to 
!! perform all of our calculations pertaining to creation and use of the 
!! our cconstructed nuclear model.
!! 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! find_best_parameters
!! construct_alpha_beta
!! calculate_linear_termns
!! print_best_parameters
!!
!!----------------------------------------------------------------------
!! Included functions:
!!
!! volume_term
!! surface_term
!! asymmetry_term
!! coulomb_term
!! pairing_term
!! semi_empirical_mass
!! semi_empirical_error
!! most_bound_nuclei
!! neutron_drip_position
!!----------------------------------------------------------------------
module nuclear_model
use types
use linear_algebra, only : solve_linear_system
implicit none

private

public :: find_best_parameters, semi_empirical_mass, semi_empirical_error, most_bound_nuclei, &
          neutron_drip_position  
contains


!-----------------------------------------------------------------------
!! Subroutine: find_best_parameters
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This is our overall subroutine in the module. 
!!
!! It sets the number of parameters to be used, then call the subroutines 
!! create alpha and beta, calls a subroutine to solve the linear solve_linear_system
!! of equations using that alpha and beta. Then using the c_parameters
!! prints them and their uncertainties. 
!!
!! Returns the c_parameters and the covariance matric to be used in 
!! theoretical predictions.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!-----------------------------------------------------------------------
!! Output:
!!
!! c_parameters     real        Array containing the semi-empirical mass formula parameters
!! covariance       real        Array containing the covariance matrix of the parameters
!-----------------------------------------------------------------------
subroutine find_best_parameters(n_protons, n_neutrons, exp_values, uncertainties, c_parameters, covariance)
    implicit none
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    real(dp), intent(in) :: exp_values(:), uncertainties(:)
    real(dp), intent(out), allocatable ::  c_parameters(:), covariance(:,:)

    integer, parameter :: n_parmaters = 6 !5

    real(dp) :: alpha(1:n_parmaters,1:n_parmaters), beta(1:n_parmaters)

    if(allocated(c_parameters)) deallocate(c_parameters)
    if(allocated(covariance)) deallocate(covariance)

    allocate(c_parameters(1:n_parmaters))
    allocate(covariance(n_parmaters, n_parmaters))
    
    call construct_alpha_beta(n_protons, n_neutrons, exp_values, uncertainties, alpha, beta)

    call solve_linear_system(alpha,beta,c_parameters,covariance)

    ! print the parameters (with it's uncertainties) to screen
    call print_best_parameters(c_parameters,covariance)
end subroutine find_best_parameters

!-----------------------------------------------------------------------
!! Subroutine: construct_alpha_beta
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! It initially checks to make sure alpha is a square matrix and alphas
!! rows contain the same number of elements as the beta vector. 
!!
!! Then loops through the data and constructs the alpha and beta arrays
!! and returns them.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!-----------------------------------------------------------------------
!! Output:
!!
!! alpha            real        Array containing the alpha matrix
!! beta             real        Array containing the beta vector
!-----------------------------------------------------------------------
subroutine construct_alpha_beta(n_protons, n_neutrons, exp_values, uncertainties, alpha, beta)
    implicit none
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    real(dp), intent(in) :: exp_values(:), uncertainties(:)
    real(dp), intent(out) :: alpha(:,:), beta(:)

    integer :: n_data, n_parmaters, alpha_shape(1:2), i, j, k
    real(dp):: linear_terms(1:size(beta))

    alpha_shape = shape(alpha)

    ! Check if the alpha array is a square matrix
    if (alpha_shape(1) /= alpha_shape(2)) then
        print *, "designated alpha is not a square matix"
        stop
    ! Also check that beta has the same number of elements as alpha has rows (or columns)
    elseif(alpha_shape(1) /= size(beta)) then
        print *, "designated size of beta doesn't equal the row size of alpha"
        stop
    endif 

    n_data = size(uncertainties) !3433
    n_parmaters = alpha_shape(1) !5 or 6

    alpha = 0._dp
    beta = 0._dp
    do k=1,n_data

        call calculate_linear_termns(n_protons(k), n_neutrons(k), linear_terms)
        do i=1,n_parmaters
            do j=1,n_parmaters
                alpha(i, j) = alpha(i, j) + ((linear_terms(i)*linear_terms(j))/(uncertainties(k)**2.0))
            enddo
            beta(i) = beta(i) + ((linear_terms(i)*exp_values(k))/(uncertainties(k)**2.0))
        enddo

    enddo


end subroutine construct_alpha_beta

!-----------------------------------------------------------------------
!! Subroutine: calculate_linear_termns
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine takes in number of protons and nuetrons within an 
!! isotope and calculates the linear terms, for the semi-empirical mass 
!! formula and stores the values into the array linear_terms.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! Z                integer     number of protons in an isotope
!! N                integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! linear_terms        real        Array containing the linear terms in the semi-empirical mass formula
!-----------------------------------------------------------------------
subroutine calculate_linear_termns(Z, N, linear_terms)
    implicit none
    integer, intent(in) :: Z, N
    real(dp), intent(out) :: linear_terms(:)

    linear_terms(1) = volume_term(Z,N)
    linear_terms(2) = surface_term(Z,N)
    linear_terms(3) = asymmetry_term(Z,N)
    linear_terms(4) = coulomb_term(Z,N)
    linear_terms(5) = pairing_term(Z,N)
    ! My extra term to increase accuracy 
    linear_terms(6) = extra_term(Z,N)
end subroutine calculate_linear_termns

!-----------------------------------------------------------------------
!! function: volume_term
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! Constructing the volume term for the c_volume term, using the 
!! number of protons and neutrons.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        volume term
!-----------------------------------------------------------------------
real(dp) function volume_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    r = Z + N
end function volume_term

!-----------------------------------------------------------------------
!! function: surface_term
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! Constructing the surface term for the c_surface term, using the 
!! number of protons and neutrons
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        surface term
!-----------------------------------------------------------------------
real(dp) function surface_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    r = (Z + N)**(2.0/3.0)
end function surface_term

!-----------------------------------------------------------------------
!! function: asymmetry_term
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! Constructing the asymmetry term for the c_asymmetry term, using the 
!! number of protons and neutrons
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        asymmetry term
!-----------------------------------------------------------------------
real(dp) function asymmetry_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    r = ((N - Z)**2.0)/(Z + N)
end function asymmetry_term

!-----------------------------------------------------------------------
!! function: coulomb_term
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! Constructing the coulomb term for the c_coulomb term, using the 
!! number of protons and neutrons
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        coulomb term
!-----------------------------------------------------------------------
real(dp) function coulomb_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    r = (Z*(Z - 1))/((Z + N)**(1.0/3.0))
end function coulomb_term

!-----------------------------------------------------------------------
!! function: pairing_term
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! Constructing the pairing term for the c_pairing term, using the 
!! number of protons and neutrons and the delta function logic.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        pairing term
!-----------------------------------------------------------------------
real(dp) function pairing_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    integer :: delta_Z_N

    ! if Z and N are both even
    if (mod(Z, 2) == 0 .and. mod(N, 2) == 0) then
        delta_Z_N = 1
    ! if Z and N are both odd
    else if (mod(Z, 2) == 1 .and. mod(N, 2) == 1) then
        delta_Z_N = -1
    ! if (Z + N) is odd
    else if (mod((Z+N), 2) == 1) then
        delta_Z_N = 0
    endif

    r = ((Z + N)**(-3.0/4.0))*delta_Z_N

end function pairing_term


!-----------------------------------------------------------------------
!! function: extra_term
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! Constructing the extra term for the c_extra term, using the 
!! number of protons and neutrons
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        extra term
!-----------------------------------------------------------------------
real(dp) function extra_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N
    real(dp) :: A
    A = Z + N
    r = 1/A
end function extra_term


!-----------------------------------------------------------------------
!! Subroutine: print_best_parameters
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine calculates the error for our c_parameters and prints 
!! the results of our c_parameters and each ones respective uncertainty.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters     real        Array containing the best fit parameters
!! covariance       real        Array containing covariance matrix
!-----------------------------------------------------------------------
subroutine print_best_parameters(c_parameters, covariance)
    implicit none
    real(dp), intent(in) :: c_parameters(:), covariance(:,:)
    real(dp) :: variance
    integer :: i, n
    real(dp), allocatable :: uncertainties(:)

    n = size(c_parameters)

    if(allocated(uncertainties)) deallocate(uncertainties)
    allocate(uncertainties(1:n))

    do i=1,n
        variance = c_parameters(i)*covariance(i, i)
        uncertainties(i) = sqrt(abs(variance))
    enddo

    print *, '' !Break up prints
    print *, 'Best fit values:              value              uncertainty'
    print 1, ' Volume parameter:   ', c_parameters(1), uncertainties(1)
    print 1, ' Surface parameter:  ', c_parameters(2), uncertainties(2)
    print 1, ' Asymmetry parameter:', c_parameters(3), uncertainties(3)
    print 1, ' Coulomb parameter:  ', c_parameters(4), uncertainties(4)
    print 1, ' Pairing term:       ', c_parameters(5), uncertainties(5)
    print 1, ' Extra parameter:    ', c_parameters(6), uncertainties(6)

1 format(a,f15.8,e28.16)
end subroutine print_best_parameters



!-----------------------------------------------------------------------
!! function: semi_empirical_mass
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This function takes in the c parameters array, and integers values Z and N for 
!! the number of protons and neutrons respectively. 
!!
!! Then calculates the semi-empirical mass using the c parameters and the
!! linear terms for each isotope.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! c    real        Array containing the parameters of the semi-empirical mass formula
!! Z    integer     number of protons in an isotope
!! N    integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! r    real        Binding energy
!-----------------------------------------------------------------------
real(dp) function semi_empirical_mass(c, Z, N) result(r)
    implicit none
    real(dp), intent(in) :: c(:)
    integer, intent(in) :: Z, N
    real(dp):: linear_terms(1:size(c))
    integer :: i

    call calculate_linear_termns(Z, N, linear_terms)

    r = 0._dp
    do i=1,size(c)
        r = r + (c(i)*linear_terms(i))
    enddo
end function semi_empirical_mass

!-----------------------------------------------------------------------
!! function: semi_empirical_error
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This function takes in the covariance array, and integers values Z and N for 
!! the number of protons and neutrons respectively. 
!!
!! Then calculates the semi-empirical error using by multiplying the linear terms
!! with their corresponding values in the covariance matrix for each isotope.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! covariance   real        2D array containing the parameters' covariance matrix
!! Z            integer     number of protons in an isotope
!! N            integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        statistical uncertainty in the binding energy
!-----------------------------------------------------------------------
real(dp) function semi_empirical_error(covariance, Z, N) result(r)
    implicit none
    real(dp), intent(in) :: covariance(:,:)
    integer, intent(in) :: Z, N
    integer :: n_row, i, j, shape_covariance(1:2)
    real(dp), allocatable :: linear_terms(:)

    shape_covariance = shape(covariance)
    n_row = shape_covariance(1)

    if(allocated(linear_terms)) deallocate(linear_terms)
    allocate(linear_terms(1:n_row))

    call calculate_linear_termns(Z, N, linear_terms)

    r = 0._dp
    do i=1,n_row
        do j=1,n_row
            r = r + linear_terms(i)*linear_terms(j)*covariance(i,j)
        enddo
    enddo

    r = sqrt(r)    

end function semi_empirical_error


!-----------------------------------------------------------------------
!! Function: most_bound_nuclei
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This function takes in a integer value for # of protons and the 
!! c_parameters for the arrays.
!!
!! It goes through each value of N neutrons (N) and checks for
!! the value of for which the binding energy per nucleon is the lowest.
!!
!! It loops indefinitely and calculates the binding energy and stores it 
!! integer in the tmp_value it stops once the binding energy per nucleon
!! is greater than the tmp_value. In other words once the binding energy 
!! per nucleon starts increasing. Then returns the value of the previous N
!! for each number of protons (Z).
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! Z                integer     number of protons in an isotope
!! c_parameters     real        Array containing the best fit parameters
!-----------------------------------------------------------------------
!! Output:
!!
!! r                integer     postion where binding energy per nucleon is the lowest
!-----------------------------------------------------------------------
integer function most_bound_nuclei(Z, c_parameters) result(r)
    implicit none
    integer, intent(in) :: Z
    real(dp) :: c_parameters(:), tmp_value, BE_Z_N
    integer :: N

    tmp_value = 0._dp

    N = 1
    do 
        BE_Z_N = semi_empirical_mass(c_parameters, Z, N)/(Z+N)
        if (N > 1 .and. BE_Z_N > tmp_value) then
            r = N - 1
            exit
        endif 
        tmp_value = BE_Z_N

        N = N + 1
    enddo

end function most_bound_nuclei


!-----------------------------------------------------------------------
!! Function: neutron_drip_position
!-----------------------------------------------------------------------
!! Austin Keller 
!!
!! This function takes in a integer value for # of protons, and the c_parameters 
!! array that we calculated from the dataset.
!!
!! It goes through each value of N neutrons (N) and checks for
!! weather the serperation energy is positive. 
!!
!! It loops indefinitely until the seperation energy switches from negative
!! to positive. Then returns the value of N that satisfies the condition for
!! each number of protons (Z).
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! Z                    integer     number of protons in an isotope
!! c_parameters         real        Array containing the best fit parameters
!-----------------------------------------------------------------------
!! Output:
!!
!! r                    integer     postion of the neutron dripline
!-----------------------------------------------------------------------
integer function neutron_drip_position(Z, c_parameters) result(r)
    implicit none
    integer, intent(in) :: Z
    real(dp), intent(in) :: c_parameters(:)
    integer :: N
    real(dp) :: BE_Z_N, BE_Z_N_minus, separation_energy

    N = 1
    separation_energy = 0._dp

    do 
        BE_Z_N_minus = semi_empirical_mass(c_parameters, Z, N-1)
        BE_Z_N = semi_empirical_mass(c_parameters, Z, N)
        separation_energy = BE_Z_N_minus - BE_Z_N

        if (separation_energy < 0) then
            r = N - 1
            exit
        endif
        N = N + 1

    enddo

end function neutron_drip_position


end module nuclear_model
