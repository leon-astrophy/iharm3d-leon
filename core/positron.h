//******************************************************************************
//*                                                                            *
//* positron.h                                                                 *
//*                                                                            *
//* Header files for the positron module                                       *
//*                                                                            *
//******************************************************************************

//import headers
#pragma once
#include "decs.h"

// do we calculate e-p pairs? //
#define PAIRS 0

// thomson scattering cross section in CGS
#define sigma_t (6.6524587051e-25)

// Leon's patch, fixup parameters for positron //
#define RPLMINLIMIT (1.e-30)
#define RPLMIN  (1.e-6)

// Leon's patch, these are all physical untis //
#define ME (9.1093826e-28  ) // Electron mass
#define MP (1.67262171e-24 ) // Proton mass
#define KBOL (1.3806505e-16  ) // Boltzmann constant
#define GNEWT (6.6742e-8      ) // Gravitational constant
#define CL (2.99792458e10  ) // Speed of light
#define R_E (2.8179403262e-13  ) // classical electron radius
#define A_F  (7.2973525693e-3  ) // fine structure constants 
#define MSUN (1.989e33        ) // Solar mass

// Leon's patch, unit conversion //
double Mbh, L_unit, T_unit, RHO_unit, U_unit, M_unit, mbh; 

// Leon's patch, important constants for pair production //
#define zeta (1.0) // geometric factor
#define eta (0.56145948356) // exp(-euler) 
#define lambda_c (2.42631023867e-10) // compton wavelength in cm
#define alphaf (7.2973525693e-3) // fine structure constant 
#define RE (2.8179403262e-13) // classical eletron radius 
#define sqrt2 (1.41421356237) // square root of 2
#define hplanck (6.626196e-27) //planck constant in cgs

// quality factor //
#define q_alpha (0.3)

// bisection threshold //
#define bisects (1e-10)

//******************************************************************************
/* define function here, which are not called globally */

void pair_production_1zone(struct GridGeom *G, struct FluidState *Ss, struct FluidState *Sf, int i, int j, int k , double dt_step);
void find_temp_1zone(struct GridGeom *G, struct FluidState *Ss, int i, int j, int k);
double nadot(double z, double nprot, double theta);
double get_ndotee(double nprot, double z, double theta);
double ncdot(double ngamma, double theta, double nprot, double z, double n1, double fb, double ndotbr);
double ndot_net(double zfrac, double taut, double nprot, double theta, double r_size);
double find_xm(double z, double tau, double nprot, double theta);
double integrate_log(double A, double theta, double xm);
double rate_ep(double z, double nprot, double theta, double xm);
double rate_ee(double z, double nprot, double theta, double xm);
double rate_pm(double z, double nprot, double theta, double xm);
double get_ndotbr(double z, double theta, double xm, double nprot);
double nbb(double x, double theta);
double n0dot(double x, double nprot, double theta);
double ndotep(double x, double z, double nprot, double theta, double dn0dt);
double ndotee(double x, double z, double nprot, double theta, double dn0dt);
double ndotpm(double x, double z, double nprot, double theta, double dn0dt);
double brem_abs(double x, double z, double nprot, double theta);
double fbrem(double y);
double comptony1(double x, double tau, double theta);
double flatn1(double x, double theta, double y);
double ngamma(double xm, double tau, double theta, double y1, double z, double nprot, double fb, double ndotbr, double r_size);
double get_ndotww(double ngamma, double theta);
double get_ndotwp(double ngamma, double nprot, double theta);
double get_ndotwe(double ngamma, double nprot, double z, double theta);
double get_ndotwf(double n1, double ngamma, double theta);

//******************************************************************************