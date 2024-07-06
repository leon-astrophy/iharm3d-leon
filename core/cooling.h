//******************************************************************************
//*                                                                            *
//* cooling.h                                                                  *
//*                                                                            *
//* Header files for cooling function                                          *
//*                                                                            *
//******************************************************************************

//import headers
#pragma once
#include "decs.h"

//parameters in the cooling functions//
double s_cool = 1.0; // cooling factor
double y_crit = 4.0; // critical cooling y parameter
double s_pow = 0.0; // power in sine angle
double q_cool = 0.5; // cooling power, currently not used

// target scale hight //
double h_r = 0.1;

//******************************************************************************