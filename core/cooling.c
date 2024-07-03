//******************************************************************************
//*                                                                            *
//* COOLING.C                                                                  * 
//*                                                                            *
//* Apply radiative cooling to accretion disk                                  *
//* Written by HS Leon Chan at 2024                                            *
//*                                                                            *
//******************************************************************************

// header files
#include "math.h"
#include "decs.h"
#include "cooling.h"

// compile only if cooling flag is on //
#if COOLING

//******************************************************************************

// Perform radiative cooling //
// Here, Ss is the fluid state in the previous step, Sf is the fluid state in the next step //
void rad_cooling(struct GridGeom *G, struct FluidState *Ss, struct FluidState *Sf, double dt_step)
{
#pragma omp parallel for collapse(3)
  ZLOOP {
    rad_cooling_1zone(G, Ss, Sf, i, j, k, dt_step);
  }
}

//******************************************************************************

// compute pair production rate per grid cells //
inline void rad_cooling_1zone(struct GridGeom *G,struct FluidState *Ss, struct FluidState *Sf, int i, int j, int k , double dt_step)
{
  // define 4-velocity 
  double u0;
  double u1;
  double u2;
  double u3; 

  // target temperature //
  double tstar;

  // coordinate //
  double rad;
  double theta;
  double X[NDIM];

  // angular velocity //
  double omg_loc;

  // mass density and internal energy //
  double rho_loc;
  double eps_loc;

  // cooling y parameter //
  double y_cool; 

  // cooling rate //
  double qdot_cool;
  double t_cool_inv;

  // coordinate radius //
  coord(i, j, k, CENT, X);
  bl_coord(X, &rad, &theta);

  // sin to the power //
  double sin_pow = abs(pow(sin(theta), s_pow));

  // first, calculate covariant 4-velocity of the last step, at face center //
  ucon_calc(G, Ss, i, j, k, CENT);
  lower_grid(Ss->ucon, Ss->ucov, G, i, j, k, CENT);

  // assign 4-velocity //
  u0 = Ss->ucov[0][k][j][i];
  u1 = Ss->ucov[1][k][j][i];
  u2 = Ss->ucov[2][k][j][i];
  u3 = Ss->ucov[3][k][j][i];

  // assign density and internal energy //
  rho_loc = Ss->P[RHO][k][j][i];
  eps_loc = Ss->P[UU][k][j][i]/Ss->P[RHO][k][j][i];

  // now determine which omega to use //
  if (rad > Risco){
    omg_loc = 1.0/(pow(rad,1.5) + a);  
  } else {
    omg_loc = 1.0/(pow(Risco,1.5) + a);  
  }

  // assign target temperature //
  tstar = M_PI_2*pow(h_r*rad*omg_loc,2);

  // assign cooling parameter //
  y_cool = (gam - 1.0)*eps_loc/tstar;

  // compute cooling rate //
  // this is Fragile 2012 approach //
  //if (y_cool < 1.0) {
  //  qdot_cool = 0.0;
  //} else if (y_cool > 2.0) {
  //  qdot_cool = 1.0;
  //} else {
  //  qdot_cool = y_cool - 1.0;
  //}

  // Alternatively, I follow Prasun Dhang's approach //
  if (y_cool > 1.0) {
    qdot_cool = fmin(y_cool - 1.0, y_crit - 1.0);
  } else {
    qdot_cool = 0.0;
  }

  // multiply it by prefactor //
  t_cool_inv = s_cool*omg_loc*sin_pow;
  qdot_cool *= t_cool_inv*rho_loc*eps_loc;

  // cool the gass //
  // Tab = Tab - (-g)^(1/2) ub qdot dt
  Sf->U[UU][k][j][i] += -qdot_cool* G->gdet[CENT][j][i]*u0 * dt_step;
  Sf->U[U1][k][j][i] += -qdot_cool* G->gdet[CENT][j][i]*u1 * dt_step;
  Sf->U[U2][k][j][i] += -qdot_cool* G->gdet[CENT][j][i]*u2 * dt_step;
  Sf->U[U3][k][j][i] += -qdot_cool* G->gdet[CENT][j][i]*u3 * dt_step;
  
  //printf("isco %.3f\n", Risco);
  //exit(0);
}

//******************************************************************************

#endif