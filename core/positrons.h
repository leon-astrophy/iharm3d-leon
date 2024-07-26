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

// do we calculate e-p pair production? //
#define PAIRS 1

//******************************************************************************

// Leon's patch, electron temperature //
#define t_elec (1e9)

// Leon's patch, fixup parameters for positron //
#define RPLMINLIMIT (1.e-30)
#define RPLMIN  (1.e-16)
#define ZMIN (1.e-10)

// bisection threshold //
#define bisects (1e-6)

//******************************************************************************

// Leon's patch, these are all physical untis //
#define ME (9.1093826e-28  ) // Electron mass
#define MP (1.67262171e-24 ) // Proton mass
#define KBOL (1.3806505e-16  ) // Boltzmann constant
#define GNEWT (6.6742e-8      ) // Gravitational constant
#define CL (2.99792458e10  ) // Speed of light
#define MSUN (1.989e33        ) // Solar mass
#define ME_MP (0.0005446170214888188) //ratio of electron to proton mass

// Leon's patch, important constants for pair production //
#define sigma_t (6.6524587051e-25) // thomson cross section
#define eta (0.56145948356) // exp(-euler) 
#define lambda_c (2.42631023867e-10) // compton wavelength in cm
#define alphaf (7.2973525693e-3) // fine structure constant 
#define RE (2.8179403262e-13) // classical eletron radius 
#define sqrt2 (1.41421356237) // square root of 2
#define hplanck (6.626196e-27) //planck constant in cgs
#define QE (4.80320680e-10) //electron charge
#define COULOMB_LOG (20.) // Coulomb logarithm

//******************************************************************************
/* define function here, which are not called globally */

void pair_production_1zone(struct GridGeom *G, struct FluidState *Ss, struct FluidState *Sf, int i, int j, int k , double dt_step);
void find_temp_1zone(struct GridGeom *G, struct FluidState *Ss, int i, int j, int k);
double nadot(double z, double nprot, double theta);
double get_ndotee(double nprot, double z, double theta);
double ncdot(double ngamma, double theta, double nprot, double z, double n1);
double ndot_net(double zfrac, double taut, double nprot, double theta, double r_size, double bfield);
double find_xm(double z, double tau, double nprot, double theta);
double integrate_log(double A, double theta, double xm);
double rate_ep(double z, double nprot, double theta, double xm);
double rate_ee(double z, double nprot, double theta, double xm);
double rate_pm(double z, double nprot, double theta, double xm);
double get_ndotbr(double z, double theta, double xm, double nprot);
double nbb(double x, double theta);
double n0dot(double x, double nprot, double theta);
double ndotep(double x, double z, double nprot, double theta);
double ndotee(double x, double z, double nprot, double theta);
double ndotpm(double x, double z, double nprot, double theta);
double brem_abs(double x, double z, double nprot, double theta);
double fbrem(double y, double taut, double theta, double xm);
double comptony1(double x, double tau, double theta);
double flatn1(double x, double theta, double y);
double ngamma(double tau, double theta, double fb, double ndotbr, double fs, double ndots, double r_size);
double get_ndotww(double ngamma, double theta);
double get_ndotwp(double ngamma, double nprot, double theta);
double get_ndotwe(double ngamma, double nprot, double z, double theta);
double get_ndotwf(double n1, double ngamma, double theta);
double get_zfrac(double nprot, double thetae);
double get_zfrac_fix_uu(double nprot, double ug, double ggas);
double series_asym(double x_in);
double ix(double x, double A_in);
double didx(double x, double A_in);
double find_xs(double thetae, double nprot, double zfrac, double v0, double h_scale);
double fraction(double x, double taut, double thetae);
void find_ndots(double thetae, double taut, double nprot, double zfrac, double h_scale, double bfield, double *fs, double *ndots);
double coulomb_onezone(struct FluidState *S, double thetae, double ue, int i, int j, int k);
double safe_Kn(int n, double x); 

//******************************************************************************
