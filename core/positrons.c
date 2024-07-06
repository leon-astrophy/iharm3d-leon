//******************************************************************************
//*                                                                            *
//* POSITRONS.C                                                                *
//*                                                                            *
//* Tracking electron-positron pairs in accretion disk                         *
//* Written by HS Leon Chan at 2024                                            *
//*                                                                            *
//******************************************************************************

// header files
#include "math.h"
#include "decs.h"
#include "positron.h"

// compile only if poistrons flag is on //
#if POSITRONS

//******************************************************************************

// set the unit convestion between code and cgs // 
void set_units()
{
  Mbh = mbh*MSUN;
  L_unit = GNEWT*Mbh/(CL*CL);
  T_unit = L_unit/CL;
  RHO_unit = M_unit*pow(L_unit,-3.);
  U_unit = RHO_unit*CL*CL;
}

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

  /* First, calculate the plasma temperature */
  /* Need temperature at the ghost zone! */
 #pragma omp parallel for collapse(3)
  ZLOOPALL {
    find_temp_1zone(G, Ss, i, j, k);
  }

  /* Then, compute the pair production rate */
#pragma omp parallel for collapse(3)
  ZLOOP {
    pair_production_1zone(G, Ss, Sf, i, j, k, dt_step);
  }

  exit(1);
}

//******************************************************************************

// compute temperature per grid cells //
inline void find_temp_1zone(struct GridGeom *G, struct FluidState *Ss, int i, int j, int k)
{

  // number density, all in c.g.s unit //
  double nprot; // proton
  double npost; // positron
  double nelec; // electron, by charge neutrality
  double ntot; // total number density

  // temperature at i,j,k //
  nprot = Ss->P[RHO][k][j][i]*RHO_unit/MP; 
  npost = Ss->P[RPL][k][j][i]*RHO_unit/ME;
  nelec = nprot + npost; 
  ntot = nelec + nprot + npost;
  temp[k][j][i] = Ss->P[UU][k][j][i]*U_unit*(gam - 1.)/(ntot*KBOL); 

}

//******************************************************************************

// compute pair production rate per grid cells //
inline void pair_production_1zone(struct GridGeom *G, struct FluidState *Ss, struct FluidState *Sf, int i, int j, int k , double dt_step)
{

  // number density, all in c.g.s unit //
  double nprot; // proton
  double npost; // positron
  double nelec; // electron, by charge neutrality
  double ntot; // total number density

  /***********************************************************************/
  
  // temperature at i-1,j,k //
  double tm1_x = temp[k][j][i-1];

  // temperature at i+1,j,k //
  double tp1_x= temp[k][j][i+1];

  // temperature at i,j-1,k //
  double tm1_y= temp[k][j-1][i];

  // temperature at i,j+1,k //
  double tp1_y= temp[k][j+1][i];

  // temperature at i,j,k-1 //
  double tm1_z= temp[k-1][j][i];

  // temperature at i,j,k+1 //
  double tp1_z= temp[k+1][j][i];

  // temperature at i,,k //
  double tc = temp[k][j][i];

  /***********************************************************************/

  // dTdx //
  double dTdx = 0.5*(tp1_x - tm1_x)/dx[1];

  // dTdy //
  double dTdy = 0.5*(tp1_y - tm1_y)/dx[2];

  // dTdz //
  double dTdz = 0.5*(tp1_z - tm1_z)/dx[3];

  // define a 4-vector //
  double gradT[NDIM];
  double gradT_con[NDIM]; 
  gradT[0] = 0.0;
  gradT[1] = dTdx;
  gradT[2] = dTdy;
  gradT[3] = dTdz;

  // raise the index  //
  for (int mu = 0; mu < NDIM; mu++) {
    gradT_con[mu] = 0.;
    for (int nu = 0; nu < NDIM; nu++) {
      gradT_con[mu] += G->gcon[CENT][mu][nu][j][i]*gradT[nu];
    }
  }

  // perform the dot product 
  double norm_gradT = sqrt(dot(gradT_con, gradT));
  
  /***********************************************************************/
  
  // number density all in CGS //
  nprot = Ss->P[RHO][k][j][i]*RHO_unit/MP; 
  npost = Ss->P[RPL][k][j][i]*RHO_unit/ME;
  nelec = nprot + npost; 
  ntot = nelec + nprot + npost;

  // get the thermal scale height, remeber to convert to CGS //
  double h_th = 0.25*tc/norm_gradT*L_unit;

  // optical depth //
  double tau_depth = 2.0*(nelec + npost)*sigma_t*h_th; 

  printf("tau_depth %.6f\n", tau_depth);

  /***********************************************************************/

  //double Te = Ss->P[UU][k][j][i]*U_unit*(gam - 1.)/(ntot*KBOL); 
  //double thetae = KBOL*Te/(ME*CL*CL);

  // e-e pair production rate //
  //double pair_prod;
  //pair_prod = get_ee_prod(thetae);
  //pair_prod = pair_prod*CL*R_E*R_E*nelec*nelec;

  // pair annihilation rate //
  //double gtheta = get_g_ann(thetae);
  //double pair_ann = M_PI*CL*R_E*R_E*nelec*npost*gtheta;

  // update positron density //
  //Sf->P[RPL][k][j][i] += ME*(pair_prod - pair_ann)*(T_unit/RHO_unit)*dt_step;
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