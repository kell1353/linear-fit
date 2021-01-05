# Project Overview
In this program will read the data form the results of the 2016 Atomic Mass Experiment (AME) file. the program will use the binding energies in the file to construct theoretical parameters to fit into our theoretical model. 

Which we will then use to calculate binding energies using our theoretical model, and predict the valley of stability, and neutron drip for all the given isotopes.


# Compilation Instructions
Compilation is all done using the makefile in the repository. Type `make` into your command line to compile the files.

Initially sets nuclear_energies as the executable name.
- Then creates the types object file.
- Then creates the linear_algebra object file using the types object file.
- Then creates the nuclear_model object file using the types, and linear_algebra object files.
- Then creates the read_write object file using the types, and nuclear_model object files.
- Then creates the main object file using the types, read_write, and nuclear_model object files.


# Usage Instructions
Once you have compiled everything execute the program (./nuclear-reactor)
The program will prompt you to type the file name EXPERIMENT_AME2016.dat. Make sure the file is located in the source files.
It then will write the results of the program into a results.dat file and a results_advanced.dat file. 
Once created you will be able to run the Jupyter file on these data files. 


# Expected Behavior
Once you have typed the file name, it then will perform calculation, print the best parameters for the model and write the results of the program into a results.dat and results_advanced.dat file.

The results file should have 6 columns. 'num_protons', 'num_neutrons', 'exp_BE', 'exp_error', 'theoretical_BE, theoretical_error'.

- num_protons: This column will contain the number of protons.
- num_neutrons: This column will contain the number of neutrons.
- exp_BE: This column will contain the experimental binding energy found during the AME.
- exp_error: This column will contain the experimental error found during the AME.
- theretical_BE: This column will contain the theoretical binding energy calculated using our model.
- theretical_error: This column will contain the theoretical energy calculated from our model.

The results_advanced file should have 3 columns. 'num_protons', 'pos_neutron_stable_isotopes', 'pos_neutron_drip'.

- num_protons: This column will contain the number of protons.
- pos_neutron_stable_isotopes: Which is the number of neutrons where the specific isotope is the most stable.
- pos_neutron_drip: Which is the largest number of neutrons where the seperation energy is still positive.

