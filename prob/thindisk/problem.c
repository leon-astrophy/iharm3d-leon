//******************************************************************************
//*                                                                            *
//* PROBLEM.C                                                                  *
//*                                                                            *
//* INITIAL CONDITIONS FOR Novikov & Thorne thin disk                          *
//*                                                                            *
//******************************************************************************

// include header files
#include "decs.h"
#include "bl_coord.h"
#include "hdf5_utils.h"

//GNU Scientific Library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//******************************************************************************

// function for thin disk //
double f_thindisk(double x_in, double x0_in, double s1_in, double s2_in, double s3_in);

//******************************************************************************

// declare
static gsl_rng *rng;

// declare variables
static double beta;
static double u_jitter;
static double theta_0;
static double kgas_0;
static double alpha_0;
static double m_0;

//******************************************************************************

//set parameter files from parameter.h
void set_problem_params() {

  // Leon's patch, Mass unit and black hole mass //
  set_param("M_unit", &M_unit);
  set_param("mbh", &mbh);

  set_param("beta", &beta);
  set_param("u_jitter", &u_jitter);
  set_param("theta_0", &theta_0);
  set_param("kgas_0", &kgas_0);
  set_param("alpha_0", &alpha_0);
  set_param("m_0", &m_0);

}

//******************************************************************************

// Save problem specific details
// This is done in each dump file in /header/problem/
void save_problem_data(hid_t string_type)
{
	hdf5_write_single_val("thindisk", "PROB", string_type);
	hdf5_write_single_val(&beta, "beta", H5T_IEEE_F64LE);
	hdf5_write_single_val(&u_jitter, "u_jitter", H5T_IEEE_F64LE);
  hdf5_write_single_val(&theta_0, "theta_0", H5T_IEEE_F64LE);
  hdf5_write_single_val(&kgas_0, "kgas_0", H5T_IEEE_F64LE);
  hdf5_write_single_val(&alpha_0, "alpha_0", H5T_IEEE_F64LE);
  hdf5_write_single_val(&m_0, "m_0", H5T_IEEE_F64LE);
}

//******************************************************************************

// calculate function f //
double f_thindisk(double x_in, double x0_in, double s1_in, double s2_in, double s3_in) {
  double zero = 1.5/pow(x_in,2.0)/(2.0*a + x_in*x_in*x_in - 3.0*x_in);
  double one = x_in - x0_in - 1.5*log(x_in/x0_in);
  double two = - 3.0*(s1_in - a)*(s1_in - a)/s1_in/(s1_in - s2_in)/(s1_in - s3_in)*log((x_in - s1_in)/(x0_in - s1_in));
  double three = - 3.0*(s2_in - a)*(s2_in - a)/s2_in/(s2_in - s1_in)/(s2_in - s3_in)*log((x_in - s2_in)/(x0_in - s2_in));
  double four = - 3.0*(s3_in - a)*(s3_in - a)/s3_in/(s3_in - s1_in)/(s3_in - s2_in)*log((x_in - s3_in)/(x0_in - s3_in));
  double out = zero*(one + two + three + four);
  return out;
}

//******************************************************************************

// initial conditions for the simulations
void init(struct GridGeom *G, struct FluidState *S)
{
  // Magnetic vector potential
  double (*A_phi)[N2 + 2*NG] = malloc(sizeof(*A_phi) * (N1 + 2*NG));

  // Initialize RNG
  // Deterministic but not symmetric across MPI procs
  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng, mpi_myrank());  

  // Inner boundar //
  Rhor = (1. + sqrt(1. - a*a));

  //setup grid first
  set_grid(G);
  LOG("Grid set");

  //declare maximum rho and internal energy (pressure)
  double rhomax = 0.;
  double umax = 0.;

  // Calculate isco radius here //
  double z1 = 1 + pow(1 - a*a,1./3.)*(pow(1+a,1./3.) + pow(1-a,1./3.));
  double z2 = sqrt(3*a*a + z1*z1);
  double rms = 3 + z2 - sqrt((3-z1)*(3 + z1 + 2*z2));
  double x0 = sqrt(rms);
  double s1 = 2*cos(acos(a)/3.0 - M_PI/3.0);
  double s2 = 2*cos(acos(a)/3.0 + M_PI/3.0);
  double s3 = -2*cos(acos(a)/3.0);

  // some reused variable //
  double xr, fx; 
  double rhoe, pe;
  double g00, g03, g33, delta, sigma;
  double u0e, u3e, u_0e, u_3e;
  double u0, u3, u_0, u_3;
  double A, B, C, F, H;
  double omg, lang, gamma_phi2;
  double gamma_100, gamma_133, gamma_103;

  // metric in BL coordinate //
  struct blgeom;
  struct of_geom blgeom;

  /*-----------------------------------------------------------------------------------------*/

  // loop over //
  ZSLOOP(-1, N3, -1, N2, -1, N1) {

    // get r and theta, same between BL and FMKS //
    double X[NDIM];
    coord(i, j, k, CENT, X);
    double r, th;
    bl_coord(X, &r, &th);
        
    // skip inside ISCO //
    if(r <= 1.2*rms) {

      // nominal value 
      S->P[RHO][k][j][i] = 0.;
      S->P[UU][k][j][i] = 0.;
      S->P[U1][k][j][i] = 0.;
      S->P[U2][k][j][i] = 0.;
      S->P[U3][k][j][i] = 0.;

    } else {

      // coordiante 
      double sth = sin(th);
      double cth = cos(th);
      blgset(i, j, &blgeom);

      // get density //
      xr = sqrt(r);
      fx = f_thindisk(xr, x0, s1, s2, s3);
      rhoe = pow(theta_0/kgas_0, 1.0/(gam - 1.0))*pow(fx/xr/xr, 0.25/(gam - 1.0));
      pe = rhoe*theta_0*pow(fx/xr/xr, 0.25);

      // get equatorial velocity //
      u3e = 1.0/sqrt(xr*xr*xr*(2.0*a + xr*xr*xr - 3.0*xr));
      u0e = (a + pow(xr,3.0))/sqrt(xr*xr*xr*(2.0*a + xr*xr*xr - 3.0*xr));

      // get metric //
      delta = r*r - 2.0*r + a*a;
      sigma = r*r + a*a*cth*cth;
      g00 = blgeom.gcov[0][0];
      g03 = blgeom.gcov[0][3];
      g33 = blgeom.gcov[3][3];

      // covariant velocity //
      u_0e = g00*u0e + g03*u3e;
      u_3e = g03*u0e + g33*u3e;

      // angular velocity and momentum //
      omg = u3e/u0e;
      lang = -u_3e/u_0e;
      gamma_phi2 = 1.0/(1.0 - omg*lang);
      F = gamma_phi2*((a*a + r*r)*(a*a + r*r) + 2*a*a*delta)/((a*a + r*r)*(a*a + r*r) + 2*a*a*delta);
      H = sqrt(pe*r*r*r/rhoe/F);

      // assign density
      S->P[RHO][k][j][i] = rhoe*exp(-alpha_0*alpha_0*r*r*cth*cth/H/H);
      S->P[UU][k][j][i] = kgas_0*pow(S->P[RHO][k][j][i], gam)/(gam - 1.0);

      // get connection coefficient in BL coordinates //
      gamma_100 = delta*(r*r - a*a*cth*cth)/(sigma*sigma*sigma);
      gamma_103 = - delta*a*sth*sth*(r*r - a*a*cth*cth)/(sigma*sigma*sigma);
      gamma_133 = delta*sth*sth*(-r*sigma*sigma + a*a*sth*sth*(r*r - a*a*cth*cth))/(sigma*sigma*sigma);
      A = gamma_100*gamma_100;
      B = g00*(gamma_100*gamma_133 - 2.0*gamma_103*gamma_103) + 2.0*g03*gamma_100*gamma_103 - g33*gamma_100*gamma_100;
      C = (gamma_103*gamma_103 - gamma_100*gamma_133)*(g03*gamma_100 - g00*gamma_103)*(g03*gamma_100 - g00*gamma_103);

      // assign angular velocity
      S->P[U3][k][j][i] = sqrt(A/(B+2.0*sqrt(C)));

      //look for maximum density and internal energy
      if (S->P[RHO][k][j][i] > rhomax) rhomax = S->P[RHO][k][j][i];
      if (S->P[UU][k][j][i] > umax) umax = S->P[UU][k][j][i];

      //why do I need this?
      S->P[UU][k][j][i] *= (1. + u_jitter * (gsl_rng_uniform(rng) - 0.5));

      //assign primitive variables in Boyer Linquist coordinates
      S->P[U1][k][j][i] = 0.;
      S->P[U2][k][j][i] = 0.;

      // Convert from 4-velocity to 3-velocity
      coord_transform(G, S, i, j, k);

      // initialize magnetic field
      S->P[B1][k][j][i] = 0.;
      S->P[B2][k][j][i] = 0.;
      S->P[B3][k][j][i] = 0.;

    }
  }

  // Look for maximum density and internal energy
  umax = mpi_max(umax);
  rhomax = mpi_max(rhomax);

  // normalize so that they are both 1
#pragma omp parallel for simd collapse(2)
  ZLOOPALL {
    S->P[RHO][k][j][i] /= rhomax;
    S->P[UU][k][j][i] /= rhomax;
  }
  umax /= rhomax;
  rhomax = 1.;

  // apply floors
  fixup(G, S);

  //boundary conditions
  set_bounds(G, S);

  /*-----------------------------------------------------------------------------------------*/

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
  norm = sqrt(beta_act / beta);
  
  // Apply normalization to magnetic fields, B3 is uniformly 0
  LOGN("Normalization is %f\n", norm);
  ZLOOP {
    S->P[B1][k][j][i] *= norm ;
    S->P[B2][k][j][i] *= norm ;
  }

  // initialize electronic variables
#if ELECTRONS
  init_electrons(G,S);
#endif 

  // apply floors and ceilings
  fixup(G, S);

  // Enforce boundary conditions
  set_bounds(G, S);

  // printout
  LOG("Finished init()");

}

//******************************************************************************

