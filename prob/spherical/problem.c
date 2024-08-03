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

// mass accretion rate and sonic radius //
static double rho_inf, r_bondi, p_inf;

// for magnetic field //
static double m_0, beta_inf;

/*****************************************************************************/

// load from the parameter file //
void set_problem_params() {
  set_param("rho_inf", &rho_inf);
  set_param("r_bondi", &r_bondi);
  set_param("beta_inf", &beta_inf);
  set_param("m_0", &m_0);
}

/*****************************************************************************/

// write to header //
void save_problem_data(hid_t string_type){
  hdf5_write_single_val("spherical", "PROB", string_type);
  hdf5_write_single_val(&rho_inf, "rho_inf", H5T_IEEE_F64LE);
  hdf5_write_single_val(&r_bondi, "r_bondi", H5T_IEEE_F64LE);
  hdf5_write_single_val(&r_bondi, "m_0", H5T_IEEE_F64LE);
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

  //declare maximum rho and internal energy (pressure)
  double rhomax = 0.;
  double umax = 0.;

  // compute sound speed at infinity //
  double cs_inf = sqrt(2.0/r_bondi);

  // compute pressure at infinity 
  p_inf = rho_inf/gam*((gam - 1.0)*cs_inf*cs_inf)/(gam - 1.0 - cs_inf*cs_inf);

  // loop over to set the initial condition //
  ZLOOP {

    // get primitive variables //
    get_prim_bondi(i, j, k, S->P, G);

    //look for maximum density and internal energy
    if (S->P[RHO][k][j][i] > rhomax) rhomax = S->P[RHO][k][j][i];
    if (S->P[UU][k][j][i] > umax) umax = S->P[UU][k][j][i];

  }
  
  /*-----------------------------------------------------------------------------------------*/

  // Leon's patch, setup magnetic field //
#if BFIELD

  // Magnetic vector potential
  double (*A_phi)[N2 + 2*NG] = malloc(sizeof(*A_phi) * (N1 + 2*NG));
  
  // Find corner-centered vector potential
#pragma omp parallel for simd collapse(2)
  ZSLOOP(0, 0, -NG+1, N2+NG-1, -NG+1, N1+NG-1) {

    // find r and theta
    double X[NDIM];
    coord(i,j,k,CORN,X);
    double r, th;
    bl_coord(X,&r,&th);
    double sth = sin(th);
    double cth = cos(th);

    // vector potentials
    //A_phi[i][j] = pow(r, 0.75)*pow(m_0, 1.25)/pow(m_0*m_0 + cth*cth, 0.625);
    A_phi[i][j] = pow(r*sth, 0.75)*pow(m_0, 1.25)/pow(m_0*m_0 + cth*cth/sth/sth, 0.625);

  } // ZSLOOP

  /*-----------------------------------------------------------------------------------------*/

  // Calculate B-field and find bsq_max
  double bsq_max = 0.;
  double beta_min = 1e100;
  ZLOOP {

    //find r and theta
    double X[NDIM];
    coord(i,j,k,CORN,X);
    double r, th;
    bl_coord(X,&r,&th);

    // B = Curl of Vector potentials
    S->P[B1][k][j][i] = -(A_phi[i][j] - A_phi[i][j + 1]
	          + A_phi[i + 1][j] - A_phi[i + 1][j + 1]) /
	          (2. * dx[2] * G->gdet[CENT][j][i]);
    S->P[B2][k][j][i] = (A_phi[i][j] + A_phi[i][j + 1]
	          - A_phi[i + 1][j] - A_phi[i + 1][j + 1]) /
	          (2. * dx[1] * G->gdet[CENT][j][i]);
    S->P[B3][k][j][i] = 0.;

    // get ucon, ucov, bcon, bcov
    get_state(G, S, i, j, k, CENT);

    // find bsquare at each grid
    double bsq_ij = bsq_calc(S, i, j, k);

    // limit bsquare
    if (bsq_ij > bsq_max) bsq_max = bsq_ij;

    // find magnetization
    double beta_ij = (gam - 1.)*(S->P[UU][k][j][i])/(0.5*(bsq_ij+SMALL)) ;

    // limit magnetization
    if(beta_ij < beta_min) beta_min = beta_ij ;
  
  }

  /*-----------------------------------------------------------------------------------------*/

  // find maximum bsqaure and magnetization
  bsq_max = mpi_max(bsq_max);
  beta_min = mpi_min(beta_min);
  
  // report maximum internal energy, bsqaure, and magnetization, and find normalization parameter
  double norm = 0;

  // Find maximum pressure to maximum bsquare 
  double beta_act = (gam - 1.) * umax / (0.5 * bsq_max);

  // print out
  LOGN("Umax is %.10e", umax);
  LOGN("bsq_max is %.10e", bsq_max);
  LOGN("beta is %.10e", beta_act);

  // normalization factor
  norm = sqrt(beta_act / beta_inf);

  // Apply normalization to magnetic fields, B3 is uniformly 0
  LOGN("Normalization is %f\n", norm);
  ZLOOP {
    S->P[B1][k][j][i] *= norm ;
    S->P[B2][k][j][i] *= norm ;
  }
#endif

  /*-----------------------------------------------------------------------------------------*/

  // electrons //
#if ELECTRONS
  init_electrons(G, S);
#endif

  //Enforce boundary conditions
  fixup(G, S);
  set_bounds(G, S);

  // printout
  LOG("Finished init()");

/*-----------------------------------------------------------------------------------------*/

}

/*****************************************************************************/

// boundary condition for bondi accretion //
void bound_gas_prob_x1r(int i, int j, int k, GridPrim  P, struct GridGeom *G)
{
  get_prim_bondi(i, j, k, P, G);
}

/*****************************************************************************/
