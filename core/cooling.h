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

//******************************************************************************

// which cooling? //
#define NOBLE 1
#define FRAGILE 2
#define PRASUN 3

// select cooling function //
// copied from HARMPI //
#define WHICHCOOL NOBLE

//******************************************************************************

// target scale hight //
#define h_r (0.1)

//******************************************************************************

//parameters in the cooling functions//
#define s_cool (1.0) // cooling factor
#define y_crit (4.0) // critical cooling y parameter
#define s_pow (0.0)// power in sine angle
#define q_cool (0.5) // cooling power, currently not used

//******************************************************************************
/* define functions here which are not called globally */

void rad_cooling_1zone(struct GridGeom *G, struct FluidState *Ss, struct FluidState *Sf, int i, int j, int k , double dt_step);

//******************************************************************************
