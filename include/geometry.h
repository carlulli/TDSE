/******************************************************************************
Module reads input values from main
sets parameters by assining input values to parameter variables
has get_params functions that return the parameters
******************************************************************************/
#ifndef GEOMETRY_H
#define GEOMETRY_H


/* [N] [TIME] [NSTEPS] [INTEGRATOR] [POTENTIAL] */
void set_params(int argc, char *argv[]);
/*
takes input values
transform input char to int, double
NUM, larger than 1 and odd; TAU, integrator_choice, pot_choice
 */
/* FUNCTIONS that return the params and warn if they changed */
int get_N();
double get_time();
int get_nsteps();
int get_integ_choice();
int get_pot_choice();

void print_N();
/* prints the currently set parameters */

#endif /* GEOMETRY_H */
