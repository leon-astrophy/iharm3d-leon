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

//initialize positrons variables
void init_cooling(struct GridGeom *G)
{
  ZLOOP {

    // coordinate //
    double rad, theta, X[NDIM];

    // coordinate radius //
    coord(i, j, k, CENT, X);
    bl_coord(X, &rad, &theta);

    // now determine which omega to use //
    if (rad > R_isco){

      // outside isco //
      omg_gr[k][j][i] = 1.0/(pow(rad,1.5) + a);  

    } else {
    
      // metric elements //
      double g_00_isco = -(1.0 - 2.0 / R_isco);
      double g_01_isco = 2.0 / R_isco;
      double g_03_isco = -2.0*a / R_isco;
      double g_13_isco = -(1.0 + 2.0 / R_isco) * a;
      double g_33_isco = pow(R_isco,2.0) + pow(a,2.0) + 2.0*pow(a,2.0) / R_isco;

      // angular velocity at the ISCO //
      double omega_isco = 1.0/(pow(R_isco,1.5) + a);

      // 4-velocity at the ISCO
      double u_0_isco = g_00_isco * 1.0 + g_03_isco * omega_isco;
      double u_1_isco = g_01_isco * 1.0 + g_13_isco * omega_isco;
      double u_3_isco = g_03_isco * 1.0 + g_33_isco * omega_isco;

      // get angular velocity within the ISCO 
      double g00 = -(1.0 + 2.0 / rad);
      double g01 = 2.0 / rad;
      double g13 = a / pow(rad,2.0);
      double g33 = 1.0 / pow(rad,2.0);
      double u0 = g00 * u_0_isco + g01 * u_1_isco;
      double u3 = g13 * u_1_isco + g33 * u_3_isco;
      omg_gr[k][j][i] = u3 / u0;

    }
  
  }
}

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
// IMPORTANT: Sf has already included flux difference and previous stage information
// DO NOT use Sf to calculate the cooling rate, use Ss instead, which is the previous 
// stage's primitive/conservative variables
inline void rad_cooling_1zone(struct GridGeom *G,struct FluidState *Ss, struct FluidState *Sf, int i, int j, int k , double dt_step)
{
  // define 4-velocity 
  double u0, u1, u2, u3; 

  // target temperature //
  double tstar;

  // coordinate //
  double rad, theta, X[NDIM];

  // mass density and internal energy //
  double rho_loc, eps_loc;

  // cooling y parameter //
  double y_cool; 

  // cooling rate and cooling time //
  double qdot_cool, t_cool_inv;

  // coordinate radius //
  coord(i, j, k, CENT, X);
  bl_coord(X, &rad, &theta);

  // first, calculate covariant 4-velocity of the last step, at face center //
  ucon_calc(G, Ss, i, j, k, CENT);
  lower_grid(Ss->ucon, Ss->ucov, G, i, j, k, CENT);

  // assign 4-velocity //
  u0 = Ss->ucov[0][k][j][i], u1 = Ss->ucov[1][k][j][i], u2 = Ss->ucov[2][k][j][i], u3 = Ss->ucov[3][k][j][i];

  // assign density and internal energy //
  rho_loc = Ss->P[RHO][k][j][i], eps_loc = Ss->P[UU][k][j][i]/Ss->P[RHO][k][j][i];

  // assign temperature //
  tstar = M_PI_2*pow(h_r*rad*omg_gr[k][j][i],2);

  // assign cooling parameter //
  y_cool = (gam - 1.0)*eps_loc/tstar;

  // compute cooling rate //
  /*-------------------------------------------------------*/
  // this is Noble 2009 approach //
#if WHICHCOOLING == NOBLE
  double Be = (1 + eps_loc*gam)*u0; 
  if (Be > -1) {
    qdot_cool = pow(y_cool - 1.0 + fabs(y_cool - 1.0), q_cool);
    qdot_cool *= s_cool*omg_gr[k][j][i]*rho_loc*eps_loc;
  } else {
    qdot_cool = 0.0;
  }
#elif WHICHPROBLEM == FRAGILE
  /*-------------------------------------------------------*/
  // this is Fragile 2012 approach //
  if (y_cool < 1.0) {
    qdot_cool = 0.0;
  } else if (y_cool > 2.0) {
    qdot_cool = 1.0;
  } else {
    qdot_cool = y_cool - 1.0;
  }
  qdot_cool *= rho_loc*eps_loc*omg_gr[k][j][i];
#elif WHICHPROBLEM == PRASUN
  /*-------------------------------------------------------*/
  // This is Prasun Dhang's approach //
  
  // sin to the power //
  double sin_pow = fabs(pow(sin(theta), s_pow));
  if (y_cool > 1.0) {
    qdot_cool = fmin(y_cool - 1.0, y_crit - 1.0);
  } else {
    qdot_cool = 0.0;
  }
  t_cool_inv = s_cool*omg_gr[k][j][i]*sin_pow;
  qdot_cool *= t_cool_inv*rho_loc*eps_loc;
#endif

  /*-------------------------------------------------------*/
  // cool the gass //

  // first, compute the source term as a time derivative //
  double duudt = -qdot_cool* G->gdet[CENT][j][i]*u0, du1dt = -qdot_cool* G->gdet[CENT][j][i]*u1;
  double du2dt = -qdot_cool* G->gdet[CENT][j][i]*u2, du3dt = -qdot_cool* G->gdet[CENT][j][i]*u3;

  // control the cooling to avoid numerical instability //
  if(fabs(duudt) > 0.0) {
    double quu = fabs(Ss->U[UU][k][j][i]/duudt);
    if (dt_step > q_alpha*quu) {
      duudt = -q_alpha*fabs(Ss->U[UU][k][j][i])/dt_step;
      printf("duudt too steep %d %d %d\n", i, j, k);
    }
  }
  if(fabs(du1dt) > 0.0) {
    double qu1 = fabs(Ss->U[U1][k][j][i]/du1dt);
    if (dt_step > q_alpha*qu1) {
      du1dt = -q_alpha*fabs(Ss->U[U1][k][j][i])/dt_step;
      printf("du1dt too steep %d %d %d\n", i, j, k);
    }
  }
  if(fabs(du2dt) > 0.0) {
    double qu2 = fabs(Ss->U[U2][k][j][i]/du2dt);
    if (dt_step > q_alpha*qu2) {
      du2dt = -q_alpha*fabs(Ss->U[U2][k][j][i])/dt_step;
      printf("du2dt too steep %d %d %d\n", i, j, k);
    }
  }
  if(fabs(du3dt) > 0.0) {
    double qu3 = fabs(Ss->U[U3][k][j][i]/du3dt);
    if (dt_step > q_alpha*qu3) {
      du3dt = -q_alpha*fabs(Ss->U[U3][k][j][i])/dt_step;
      printf("du3dt too steep %d %d %d\n", i, j, k);
    }
  }

  // Tab = Tab - (-g)^(1/2) ub qdot dt
  Sf->U[UU][k][j][i] += duudt, Sf->U[U1][k][j][i] += du1dt;
  Sf->U[U2][k][j][i] += du2dt, Sf->U[U3][k][j][i] += du3dt;
  
}

//******************************************************************************

#endif