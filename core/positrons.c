//******************************************************************************
//*                                                                            *
//* POSITRONS.C                                                                *
//*                                                                            *
//* Tracking electron-positron pairs in accretion disk                         *
//* Written by HS Leon Chan at 2024                                            *
//*                                                                            *
//******************************************************************************

// header files
#include "decs.h"

// compile only if poistrons flag is on //
#if POSITRONS

//******************************************************************************

//initialize positrons variables
void init_positrons(struct GridGeom *G, struct FluidState *S)
{
  ZLOOPALL {
    // Set initial e-p mass to be the rest mass 
    S->P[RPL][k][j][i] = RPLMINLIMIT; //S->P[RHO][k][j][i];
  }

  // Necessary?  Usually called right afterward
  set_bounds(G, S);
}

//******************************************************************************

// Compute net pair production rate //
void pair_production(struct GridGeom *G, struct FluidState *Ss, struct FluidState *Sf, double dt_step)
{
#pragma omp parallel for collapse(3)
  ZLOOP {
    pair_production_1zone(G, Ss, Sf, i, j, k, dt_step);
  }
}

//******************************************************************************

// compute pair production rate per grid cells //
inline void pair_production_1zone(struct GridGeom *G, struct FluidState *Ss, struct FluidState *Sf, int i, int j, int k , double dt_step)
{
  // number density, all in c.g.s unit //
  double nprot = Ss->P[RHO][k][j][i]*RHO_unit/MP; // proton
  double npost = Ss->P[RPL][k][j][i]*RHO_unit/ME; // positron
  double nelec = nprot + npost; // electron, by charge neutrality
  double ntot = nelec + nprot + npost; // total number density

  // optical depth //
  double tau_depth; 

  // temperature graidents //
  double tm1_x;
  double tp1_x;
  double tm1_y;
  double tp1_y;
  double tm1_z;
  double tp1_z;

  // temperature //
  double Te = Ss->P[UU][k][j][i]*U_unit*(gam - 1.)/(ntot*KBOL); 
  double thetae = KBOL*Te/(ME*CL*CL);

  // pair annihilation rate //
  double gtheta = get_g_ann(thetae);
  double pair_ann = M_PI*CL*R_E*R_E*nelec*npost*gtheta;

  // e-e pair production rate //
  double pair_prod;
  pair_prod = get_ee_prod(thetae);
  pair_prod = pair_prod*CL*R_E*R_E*nelec*nelec;

  // update positron density //
  Sf->P[RPL][k][j][i] += ME*(pair_prod - pair_ann)*(T_unit/RHO_unit)*dt_step;
}

//******************************************************************************

// functions for computing function g(theta) in pair annhilation rates //
inline double get_g_ann(double theta)
{
 
  // calculate the complicated functions //
  double gtheta_inv = (1.0 + (2.0*theta*theta/log(1.12*theta + 1.3)));
  double gtheta = pow(gtheta_inv, -1);

  // return
  return gtheta;
}

//******************************************************************************

// functions for computing e-e pair production rates //
inline double get_ee_prod(double theta)
{
 
  // calculate the complicated functions //
  double gtheta;

  // Sort by condition //
  if(theta >= 1) {
    double logt = log(theta);
    double tmp = 1 + 0.058/theta;
    gtheta = (112.0/27.0/M_PI)*A_F*A_F*logt*logt*logt*pow(tmp, -1);
  } else {
    gtheta = 2e-4*pow(theta,1.5)*exp(-2/theta)*(1.0 + 0.015*theta);
  }

  // return
  return gtheta;
}

//******************************************************************************

#endif