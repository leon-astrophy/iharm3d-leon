/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR BONDI INFLOW                                        *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include "bl_coord.h"
#include "hdf5_utils.h"

// Rootfinding for analytic Bondi solution
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

/*****************************************************************************/

// some variabls //
static double C4, C3, n, K;

// mass accretion rate and sonic radius //
static double mdots, rs;

/*****************************************************************************/

// load from the parameter file //
void set_problem_params() {
  set_param("mdots", &mdots);
  set_param("rs", &rs);
  set_param("Rhor", &Rhor);
}

/*****************************************************************************/

// write to header //
void save_problem_data(hid_t string_type){
  hdf5_write_single_val("bondi", "PROB", string_type);
  hdf5_write_single_val(&mdots, "mdots", H5T_IEEE_F64LE);
  hdf5_write_single_val(&rs, "rs", H5T_IEEE_F64LE);
}

/*****************************************************************************/

// Get temperature, Adapted from M. Chandra
double get_Tfunc(double T, double r)
{
  return pow(1.+(1.+n)*T,2.)*(1.-2./r+pow(C4/r/r/pow(T,n),2.))-C3;
}

/*****************************************************************************/

// Get temperature //
double get_T(double r)
{
  double rtol = 1.e-12;
  double ftol = 1.e-14;
  double Tmin = 0.6*(sqrt(C3) - 1.)/(n + 1);
  double Tmax = pow(C4*sqrt(2./r/r/r),1./n);
  double f0, f1, fh;
  double T0, T1, Th;
  T0 = 0.6*Tmin;
  f0 = get_Tfunc(T0, r);
  T1 = Tmax;
  f1 = get_Tfunc(T1, r);

  if (f0*f1 > 0.) {
    printf("Failed solving for T at r = %e C4 = %e C3 = %e\n", r, C4, C3);
    exit(-1);
  }

  Th = (f1*T0 - f0*T1)/(f1 - f0);
  fh = get_Tfunc(Th, r);
  double epsT = rtol*(Tmin + Tmax);
  while (fabs(Th - T0) > epsT && fabs(Th - T1) > epsT && fabs(fh) > ftol) {
    if (fh*f0 < 0.) {
      T0 = Th;
      f0 = fh;
    } else {
      T1 = Th;
      f1 = fh;
    }

    Th = (f1*T0 - f0*T1)/(f1 - f0);
    fh = get_Tfunc(Th, r);
  }

  return Th;
}

/*****************************************************************************/

// convert from 4-velocity to 3-velocity in the Eulerian Observer frame //
void fourvel_to_prim(double ucon[NDIM], GridPrim P,
  struct GridGeom *G, int i, int j, int k)
{
  double alpha, beta[NDIM], gamma;

  alpha = 1.0/sqrt(-G->gcon[CENT][0][0][j][i]);
  beta[1] = alpha*alpha*G->gcon[CENT][0][1][j][i];
  beta[2] = alpha*alpha*G->gcon[CENT][0][2][j][i];
  beta[3] = alpha*alpha*G->gcon[CENT][0][3][j][i];
  gamma = ucon[0]*alpha;

  P[U1][k][j][i] = ucon[1] + beta[1]*gamma/alpha;
  P[U2][k][j][i] = ucon[2] + beta[2]*gamma/alpha;
  P[U3][k][j][i] = ucon[3] + beta[3]*gamma/alpha;
}

/*****************************************************************************/

// normalize and set the time-component of the contravariant velocity //
void set_ut(double ucon[NDIM], struct of_geom *geom)
{
  double AA, BB, CC;

  AA = geom->gcov[0][0];
  BB = 2.*(geom->gcov[0][1]*ucon[1] +
           geom->gcov[0][2]*ucon[2] +
           geom->gcov[0][3]*ucon[3]);
  CC = 1. + geom->gcov[1][1]*ucon[1]*ucon[1] +
            geom->gcov[2][2]*ucon[2]*ucon[2] +
            geom->gcov[3][3]*ucon[3]*ucon[3] +
       2. *(geom->gcov[1][2]*ucon[1]*ucon[2] +
            geom->gcov[1][3]*ucon[1]*ucon[3] +
            geom->gcov[2][3]*ucon[2]*ucon[3]);

  double discr = BB*BB - 4.*AA*CC;
  ucon[0] = (-BB - sqrt(discr))/(2.*AA);
}

/*****************************************************************************/

// get primitiv variables //
void get_prim_bondi(int i, int j, int k, GridPrim P, struct GridGeom *G)
{
  // initialize the problem parameter //
  static int firstc = 1;
  if (firstc) {

    // set adiabatic index n //
    n = 1./(gam - 1.);

    // uc - sonic point 4-velocity, Vc - sonic point sound speed //
    double uc = sqrt(1/(2.*rs));
    double Vc = sqrt(pow(uc,2)/(1. - 3.*pow(uc,2)));
    double Tc = -n*pow(Vc,2)/((n + 1.)*(n*pow(Vc,2) - 1.));
    C4 = uc*pow(rs,2)*pow(Tc,n);
    C3 = pow(1. + (1. + n)*Tc,2)*(1. - 2./rs + pow(uc, 2));
		K  = pow(4*M_PI*C4 / mdots, 1/n);

    firstc = 0;
  }

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

  //double T = T_bondi[j][i];
  double T = get_T(r);
  double ur = -C4/(pow(T,n)*pow(r,2));
  double rho = pow(K, -n)*pow(T, n);
  double u = rho*T / (gam - 1.);
  double ucon_bl[NDIM], ucon_ks[NDIM], ucon_mks[NDIM];
  struct of_geom geom_bl;

  //set metric tensors for Boyer-Lindquist coordinates // 
  blgset(i, j, &geom_bl);

  // zero the 4-velocity //
  DLOOP1 {
    ucon_bl[mu] = 0.;
    ucon_ks[mu] = 0.;
    ucon_mks[mu] = 0.;
  }

  // assign radial 4-velocity to the BL coordiante //
  ucon_bl[1] = ur;

  // solve for time component //
  set_ut(ucon_bl, &geom_bl);

  // convert from BL to KS //
  bl_to_ks(X, ucon_bl, ucon_ks);

  // jacobian matrix //
  double dxdX[NDIM][NDIM], dXdx[NDIM][NDIM];

  // set the jacobian matrix // 
  set_dxdX(X, dxdX);
  invert(&dxdX[0][0], &dXdx[0][0]);

  // transform from KS to FMKS //
  DLOOP2 {
    ucon_mks[mu] += dXdx[mu][nu]*ucon_ks[nu];
  }

  // convert from 4-velocity to 3-velocity //
  fourvel_to_prim(ucon_mks, P, G, i, j, k);

  // assign primitive varibles //
  P[RHO][k][j][i] = rho;
  P[UU][k][j][i] = u;
  P[B1][k][j][i] = 0.;
  P[B2][k][j][i] = 0.;
  P[B3][k][j][i] = 0.;

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

  // set the computational grid and print out //
  set_grid(G);
  LOG("Grid set");

  // loop over to set the initial condition //
  ZLOOP {
    get_prim_bondi(i, j, k, S->P, G);
  }

  // print out for debug //
  if (DEBUG && mpi_io_proc()) {
    printf("a = %e Rhor = %e\n", a, Rhor);
    printf("mdots = %e\n", mdots);
    printf("rs   = %e\n", rs);
    printf("n    = %e\n", n);
    printf("C4   = %e\n", C4);
    printf("C3   = %e\n", C3);
  }

  // electrons //
#if ELECTRONS
  init_electrons(G, S);
#endif

  //Enforce boundary conditions
  fixup(G, S);
  set_bounds(G, S);

}

/*****************************************************************************/

// boundary condition for bondi accretion //
void bound_gas_prob_x1r(int i, int j, int k, GridPrim  P, struct GridGeom *G)
{
  get_prim_bondi(i, j, k, P, G);
}

/*****************************************************************************/