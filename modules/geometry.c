/******************************************************************************
Module reads input values from main
sets parameters by assining input values to parameter variables
has get_params functions that return the parameters
******************************************************************************/
#include <stdlib.h>
#include <stdio.h>

/*
defined inside c file: variables can only be accessed from this file
and returned from get_functions OR accessed from other files with extern keyword
static: can not be accessed with extern keyword (safest way)
*/
static int N=0;
static double time=0.0;
static int nsteps = 1;
static int integ_choice=0;
static int pot_choice=0;
static int *integcount=NULL;
static int *potcount=NULL;

/*
function set_params
takes input values
transform input char to int
assign input int to variables (global?):
NUM, larger than 1 and odd;
A;
eps (tolerance), around 10e-15;
 */
void set_params(int argc, char *argv[]){
  // is argc the correct number of args
  // if correct loop from 1 to argc and read argv
  // plan to rewrite for all input parameters!!!!!!
  if (argc<6) {
    printf("[geometry.c | set_params()] ERROR. You forgot to insert program parameters!\n"
  "\tThe parameters are as follows: [N] [time] [Nsteps] [integrator choice] [potential_choice]\n"
  "\t\tN: odd and bigger 0\n" "\t\t time < 0 an Nsteps integer\n"
  "\t\tintegrator_choice: 0 = euler method, 1 = Unitary Crank Nicolson Method, 2 = Strang Splitting Method\n"
"\t\tintegrator_choice: 0 = zero, 1 = harmonic oscillator, 2 = well, 3 = wall\n" );
    exit(-1);
  }
  int val;
  val = atoi(argv[1]);
  if (val > 1 && val % 2) {
    N = val;
  } else {
    printf("[ geometry.c| set_params ] NUM has to be > 1 and an odd number!\n");
    exit(-1);
  }
  if (argv[2] != NULL) { time = atof(argv[2]); }
  if (argv[3] != NULL) { nsteps = atoi(argv[3]); }
  if (argv[4] != NULL) {
    integ_choice = atoi(argv[4]);
    integcount = &integ_choice;
  }
  if (argv[5] != NULL) {
    pot_choice = atoi(argv[5]);
    potcount = &pot_choice;
  }
}

/* functions that return the parameters and warn if they were changed */
int get_N() {
    static int sN=0;
    if(N==0) {
      printf("[ geometry.c| get_N ] Error! N not yet set!\n");
      exit(-1);
    }
    if(sN==0) { sN=N; }
    else {
      if((N!=sN))  {
         printf("[ geometry.c| get_N ] Error! (N) has changed: (%d) -> (%d)\n",sN,N);
         exit(-1);
      }
    }
    return N;
}

double get_time() {
  static double sTIME=0.0;
  if(time==0.0) {
    printf("[ geometry.c| get_time()] Error! Timenot yet set or tried to set 0 (not possible)!\n");
    exit(-1);
  }
  if(sTIME==0.0) { sTIME=time; }
  else {
    if((time!=sTIME))  {
       printf("[ geometry.c| get_time() ] Error! (time) has changed: (%.e) -> (%.e)\n",sTIME,time);
       exit(-1);
    }
  }
  return time;
}

double get_nsteps() {
  return nsteps;
}

int get_integ_choice() {
    static int sINTEG_CHOICE=0;
    if(integcount==NULL) {
      printf("[ geometry.c| get_integ_choice() ] Error! Integrator_choice not yet set!\n");
      exit(-1);
    }
    if(sINTEG_CHOICE==0) { sINTEG_CHOICE=integ_choice; }
    else {
      if((integ_choice!=sINTEG_CHOICE))  {
         printf("[ geometry.c| get_integ_choice() ] Error! (integ_choice) has changed: (%d) -> (%d)\n",sINTEG_CHOICE,integ_choice);
         exit(-1);
      }
    }
    return integ_choice;
}

int get_pot_choice() {
    static int sPOT_CHOICE=0;
    if(potcount==NULL) {
      printf("[ geometry.c| get_pot_choice() ] Error! Potential_choice not yet set!\n");
      exit(-1);
    }
    if(sPOT_CHOICE==0) { sPOT_CHOICE=pot_choice; }
    else {
      if((pot_choice!=sPOT_CHOICE))  {
         printf("[ geometry.c| get_pot_choice() ] Error! (pot_choice) has changed: (%d) -> (%d)\n",sPOT_CHOICE,pot_choice);
         exit(-1);
      }
    }
    return pot_choice;
}

/* prints the currently set parameters */
void print_N() {
  printf("The parameters set from the input are: " "NUM =\t%d \n",N );
}
