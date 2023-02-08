# TCCM-Programming-Project
Here will be my codes for simulation of a liquid using MC approach.


So in this project I neded to implement simulation of Lennard-Jones liquid. The instructions are in the file tccm_distance_unit_manual_2022.pdf.

Initial code is the LJ_no_nl.cpp file.
A more efficient (only for big number of atoms) code, using neighbour lists is the file LJ_nl_final_code.cpp. 
I compared the efficiency of these 2 codes in an excel file (efficiency from num of atoms.xlsx). 

Both of these codes output a txt file with g(r), which can be plotted with g_of_r_lennard_jones_final.dat for reference. 


And finally, the general task was to create a programto simulate Stillinger model for a Silicon liquid. This is implemented in the file still_code.cpp. This file takes prod_struc.xyz as input and outputs a txt file, which needs to be compared with gofr_stillinger_final(1).dat. 


