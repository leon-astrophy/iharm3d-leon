/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR BONDI INFLOW                                        *
 *                                                                            *
 ******************************************************************************/

/* include header */
#include "decs.h"
#include "bl_coord.h"
#include "hdf5_utils.h"

// Rootfinding for analytic Bondi solution
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

/*****************************************************************************/

// some variabls //
double C4, C3, n, K;

// mass accretion rate and sonic radius //
double rho_inf, r_bondi, cs_inf, p_inf, beta_inf;

/*****************************************************************************/

// load from the parameter file //
void set_problem_params() {
  set_param("rho_inf", &rho_inf);
  set_param("r_bondi", &r_bondi);
  set_param("beta_inf", &beta_inf);
}

/*****************************************************************************/

// write to header //
void save_problem_data(hid_t string_type){
  hdf5_write_single_val("spherical", "PROB", string_type);
  hdf5_write_single_val(&rho_inf, "rho_inf", H5T_IEEE_F64LE);
  hdf5_write_single_val(&r_bondi, "r_bondi", H5T_IEEE_F64LE);
}

/*****************************************************************************/

// get primitiv variables //
void get_prim_bondi(int i, int j, int k, GridPrim P, struct GridGeom *G)
{

  // get coordinate //
  double r, th, X[NDIM];
  coord(i, j, k, CENT, X);
  bl_coord(X, &r, &th);

  // skip the cells inside inner boundary //
  while (r < Rhor) {
    i++;
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);
  }

  // assign primitive varibles //
  P[RHO][k][j][i] = rho_inf;
  P[UU][k][j][i] = p_inf/(gam - 1.0);
  P[U1][k][j][i] = 0.0;
  P[U2][k][j][i] = 0.0;
  P[U3][k][j][i] = 0.0;
  P[B1][k][j][i] = 0.0;
  P[B2][k][j][i] = 0.0;
  P[B3][k][j][i] = 0.0;
  
  // Electrons make no physical sense here but are a very useful debug tool
  // At least set them consistently here to test deterministic evolution
#if ELECTRONS
    // Set electron internal energy to constant fraction of internal energy
    double uel = fel0*P[UU][k][j][i];

    // Initialize entropies
    P[KTOT][k][j][i] = (gam-1.)*P[UU][k][j][i]*pow(P[RHO][k][j][i],-gam);
    P[KEL][k][j][i] = (game-1.)*uel*pow(P[RHO][k][j][i],-game);
#endif

}

/*****************************************************************************/

// initialize the bondi solution // 
void init(struct GridGeom *G, struct FluidState *S)
{

  // inner boundary, set it before setting grid //
  Rhor = (1. + sqrt(1. - a*a));

  // set the computational grid and print out //
  set_grid(G);
  LOG("Grid set");

  // compute sound speed at infinity //
  double cs_inf = sqrt(2.0/r_bondi);

  // compute pressure at infinity 
  double p_inf = rho_inf*cs_inf*cs_inf/gam;

  // loop over to set the initial condition //
  ZLOOP {
    get_prim_bondi(i, j, k, S->P, G);
  }

  // electrons //
#if ELECTRONS
  init_electrons(G, S);
#endif

  //Enforce boundary conditions
  fixup(G, S);
  set_bounds(G, S);

  // Leon's patch, setup magnetic field //
#if BFIELD
  ZLOOPALL {

    // get coordinate //
    double r, th, X[NDIM];
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);

    // bipolar field //
    if (th < M_PI_2) {
      S->P[B1][k][j][i] = 1.0/r;
    } else if (th > M_PI_2) {
      S->P[B1][k][j][i] = -1.0/r;
    } else {
      S->P[B1][k][j][i] = 0.0;
    }
    
    // get ucon, ucov, bcon, bcov
    get_state(G, S, i, j, k, CENT);

    // find bsquare at each grid
    double bsq_ij = bsq_calc(S, i, j, k);

    // find magnetization
    double beta_ij = (gam - 1.)*(S->P[UU][k][j][i])/(0.5*(bsq_ij+SMALL)) ;

    // normalization factor
    double norm = sqrt(beta_ij / beta_tgrt);

    // normailize the field //
    S->P[B1][k][j][i] *= norm;
    
  }
  fixup(G, S);
  set_bounds(G, S);
#endif

}

/*****************************************************************************/

// boundary condition for bondi accretion //
void bound_gas_prob_x1r(int i, int j, int k, GridPrim  P, struct GridGeom *G)
{
  get_prim_bondi(i, j, k, P, G);
}

/*****************************************************************************/