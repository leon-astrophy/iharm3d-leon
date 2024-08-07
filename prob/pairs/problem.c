//******************************************************************************
//*                                                                            *
//* PROBLEM.C                                                                  *
//*                                                                            *
//* INITIAL CONDITIONS FOR FISHBONE-MONCRIEF TORUS                             *
//*                                                                            *
//******************************************************************************

// include header files
#include "bl_coord.h"
#include "decs.h"
#include "hdf5_utils.h"

//GNU Scientific Library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// declare
static gsl_rng *rng;

// Local declarations
double lfish_calc(double rmax);

// Different MAD initializations
#define SANE 0
#define RYAN 1
#define R3S3 2
#define GAUSSIAN 3
#define NARAYAN 4

// Alternative normalization from HARMPI. Dosn't seem common to use
static int maxr_normalization = 0;

/////////////////////////////////////////////////////
// TODO allow initialization of mad_type w/string?
/////////////////////////////////////////////////////

// declare variables
static int mad_type;
static double BHflux, beta;
static double rin, rmax;
static double rBstart, rBend;
static double u_jitter;

//******************************************************************************

//set parameter files from parameter.h
void set_problem_params() {

  // Leon's patch, Mass unit and black hole mass //
  set_param("M_unit", &M_unit);
  set_param("mbh", &mbh);

  set_param("rin", &rin);
  set_param("rmax", &rmax);
  set_param("u_jitter", &u_jitter);

  set_param("mad_type", &mad_type);
  set_param("BHflux", &BHflux);
  set_param("beta", &beta);

  set_param("rBstart", &rBstart);
  set_param("rBend", &rBend);
}

//******************************************************************************

// Save problem specific details
// This is done in each dump file in /header/problem/
void save_problem_data(hid_t string_type)
{
	hdf5_write_single_val(&mad_type, "mad_type", H5T_STD_I32LE);
	hdf5_write_single_val("torus", "PROB", string_type);
	hdf5_write_single_val(&rin, "rin", H5T_IEEE_F64LE);
	hdf5_write_single_val(&rmax, "rmax", H5T_IEEE_F64LE);
	hdf5_write_single_val(&beta, "beta", H5T_IEEE_F64LE);
	hdf5_write_single_val(&u_jitter, "u_jitter", H5T_IEEE_F64LE);
	hdf5_write_single_val(&BHflux, "bhflux", H5T_IEEE_F64LE);
	hdf5_write_single_val(&rBstart, "rBstart", H5T_IEEE_F64LE);
	hdf5_write_single_val(&rBend, "rBend", H5T_IEEE_F64LE);
}

//******************************************************************************

// initial conditions for the simulations
void init(struct GridGeom *G, struct FluidState *S)
{
  // Magnetic field
  double (*A)[N2 + 2*NG] = malloc(sizeof(*A) * (N1 + 2*NG));

  // Initialize RNG
  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng, mpi_myrank());  // Deterministic but not symmetric across MPI procs

  // Fishbone-Moncrief parameters, angular momentum and kappa
  double l = lfish_calc(rmax);
  double kappa = 1.e-3;

  // Grid parameters
  Rhor = (1. + sqrt(1. - a*a));

  ////////////////////////////////////////////////////////////////////
  // These are never used but I'm keeping them around in case
  // Sorry, I need this now, by Leon, but not here !
  //double z1 = 1 + pow(1 - a*a,1./3.)*(pow(1+a,1./3.) + pow(1-a,1./3.));
  //double z2 = sqrt(3*a*a + z1*z1);
  //Risco = 3 + z2 - sqrt((3-z1)*(3 + z1 + 2*z2));
  ///////////////////////////////////////////////////////////////////

  //setup grid first
  set_grid(G);
  LOG("Grid set");

  //////////////////////////////////
  // TODO put this in a function
  /////////////////////////////////

  //declare maximum rho and internal energy (pressure)
  double rhomax = 0.;
  double umax = 0.;

  ZSLOOP(-1, N3, -1, N2, -1, N1) {

    //get r and theta
    double X[NDIM];
    coord(i, j, k, CENT, X);
    double r, th;
    bl_coord(X, &r, &th);
    double sth = sin(th);
    double cth = cos(th);

    // prepare to calculate lnh
    double DD = r * r - 2. * r + a * a;
    double AA = (r * r + a * a) * (r * r + a * a) -
             DD * a * a * sth * sth;
    double SS = r * r + a * a * cth * cth;

    //still preparing
    double thin = M_PI / 2.;
    double sthin = sin(thin);
    double cthin = cos(thin);

    //still preparing
    double DDin = rin * rin - 2. * rin + a * a;
    double AAin = (rin * rin + a * a) * (rin * rin + a * a)
             - DDin * a * a * sthin * sthin;
    double SSin = rin * rin + a * a * cthin * cthin;

    //now, calculates lnh
    double lnh;
    if (r >= rin) {
      lnh =
          0.5 *
          log((1. +
         sqrt(1. +
              4. * (l * l * SS * SS) * DD / (AA * AA * sth * sth)))
        / (SS * DD / AA))
          - 0.5 * sqrt(1. +
           4. * (l * l * SS * SS) * DD /
           (AA * AA * sth * sth))
          - 2. * a * r * l / AA -
          (0.5 *
           log((1. +
          sqrt(1. +
               4. * (l * l * SSin * SSin) * DDin /
               (AAin * AAin * sthin * sthin))) /
         (SSin * DDin / AAin))
           - 0.5 * sqrt(1. +
            4. * (l * l * SSin * SSin) * DDin / (AAin * AAin * sthin * sthin))
           - 2. * a * rin * l / AAin);
    } else {
      lnh = 1.;
    }

    // regions outside the torus
    if (lnh < 0. || r < rin) {
      // Nominal values; real value set by the fixup subroutine
      S->P[RHO][k][j][i] = 1.e-7 * RHOMIN;
      S->P[UU][k][j][i] = 1.e-7 * UUMIN;
      S->P[U1][k][j][i] = 0.;
      S->P[U2][k][j][i] = 0.;
      S->P[U3][k][j][i] = 0.;
    }
    /* region inside magnetized torus; u^i is calculated in
     * Boyer-Lindquist coordinates, as per Fishbone & Moncrief,
     * so it needs to be transformed at the end */
    else {
      //declare variables
      double hm1 = exp(lnh) - 1.;
      double rho = pow(hm1 * (gam - 1.) / (kappa * gam),
               1. / (gam - 1.));
      double u = kappa * pow(rho, gam) / (gam - 1.);

      // Calculate u^phi (rotational velocity)
      double expm2chi = SS * SS * DD / (AA * AA * sth * sth);
      double up1 =
          sqrt((-1. +
          sqrt(1. + 4. * l * l * expm2chi)) / 2.);
      double up = 2. * a * r * sqrt(1. +
                 up1 * up1) / sqrt(AA * SS *
                 DD) +
          sqrt(SS / AA) * up1 / sth;

      // assign density
      S->P[RHO][k][j][i] = rho;

      //look for maximum density and internal energy
      if (rho > rhomax) rhomax = rho;
      if (u > umax && r > rin) umax = u;

      //why do I need this?
      u *= (1. + u_jitter * (gsl_rng_uniform(rng) - 0.5));

      //assign primitive variables in Boyer Linquist coordinates
      S->P[UU][k][j][i] = u;
      S->P[U1][k][j][i] = 0.;
      S->P[U2][k][j][i] = 0.;
      S->P[U3][k][j][i] = up;

      // Convert from 4-velocity to 3-velocity
      coord_transform(G, S, i, j, k);
    }

    // initialize magnetic field
    S->P[B1][k][j][i] = 0.;
    S->P[B2][k][j][i] = 0.;
    S->P[B3][k][j][i] = 0.;
  } // ZSLOOP

  // Find the zone in which rBend of Narayan condition resides
  // This just uses the farthest process in R
  // For /very/ large N1CPU it might fail
  // TODO only do this for NARAYAN condition
  int iend_global = 0;
  if (global_stop[0] == N1TOT && global_start[1] == 0 && global_start[2] == 0) {
    int iend = NG;
    double r_iend = 0.0;
    while (r_iend < rBend) {
      iend++;
      double Xend[NDIM];
      coord(iend,N2/2+NG,NG,CORN,Xend);
      double thend;
      bl_coord(Xend,&r_iend,&thend);
    }
    iend--;
    iend_global = global_start[0] + iend - NG; //Translate to coordinates for uu_mid below
    if(DEBUG) printf("[MPI %d] Furthest torus zone is %d (locally %d), at r = %f\n", mpi_myrank(), iend_global, iend, r_iend);
  }
  iend_global = mpi_reduce_int(iend_global);

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

  // Calculate UU (internal energy) along midplane, propagate to all processes
  double *uu_plane_send = calloc(N1TOT,sizeof(double));

  // This relies on an even N2TOT /and/ N2CPU
  if ((global_start[1] == N2TOT/2 || N2CPU == 1) && global_start[2] == 0) {
    int j_mid = N2TOT/2 - global_start[1] + NG;
    int k = NG; // Axisymmetric
    ILOOP {
      int i_global = global_start[0] + i - NG;
      uu_plane_send[i_global] = 0.25*(S->P[UU][k][j_mid][i] + S->P[UU][k][j_mid][i-1] +
                          S->P[UU][k][j_mid-1][i] + S->P[UU][k][j_mid-1][i-1]);
    }
  }

  // mpi stuff ...
  double *uu_plane = calloc(N1TOT,sizeof(double));
  mpi_reduce_vector(uu_plane_send, uu_plane, N1TOT);
  free(uu_plane_send);

  //////////////////////////////////////////////////////////////////////////////////////////////
  //printf ("UU in plane is "); for (int i =0; i < N1TOT; i++) printf("%.10e ", uu_plane[i]);
  //////////////////////////////////////////////////////////////////////////////////////////////

  // Find corner-centered vector potential
#pragma omp parallel for simd collapse(2)
  ZSLOOP(0, 0, -NG+1, N2+NG-1, -NG+1, N1+NG-1) {

    // find r and theta
    double X[NDIM];
    coord(i,j,k,CORN,X);
    double r, th;
    bl_coord(X,&r,&th);

    // declare
    double q;

    // cell-cornered density and internal energy
    double rho_av = 0.25*(S->P[RHO][k][j][i] + S->P[RHO][k][j][i-1] +
                    S->P[RHO][k][j-1][i] + S->P[RHO][k][j-1][i-1]);
    double uu_av = 0.25*(S->P[UU][k][j][i] + S->P[UU][k][j][i-1] +
                         S->P[UU][k][j-1][i] + S->P[UU][k][j-1][i-1]);

    // declare end point variables 
    int i_global = global_start[0] + i - NG;
    double uu_plane_av = uu_plane[i_global];
    double uu_end = uu_plane[iend_global];

    // find the q-parameter for MAD discs, 3D case
    if (N3 > 1) {

      // Standard and normal evolution
      if (mad_type == SANE) {
        q = rho_av/rhomax - 0.2;

      // BR's smoothed poloidal in-torus
      } else if (mad_type == RYAN) {
        q = pow(sin(th),3)*pow(r/rin,3.)*exp(-r/400)*rho_av/rhomax - 0.2;

      // Just the r^3 sin^3 th term, proposed for EHT standard MAD models
      } else if (mad_type == R3S3) {
        q = pow(r/rin,3.)*rho_av/rhomax - 0.2;

      // Gaussian-strength vertical threaded field
      } else if (mad_type == GAUSSIAN) {

        //Radius of half-maximum. Units of rin
        double wid = 2;
        q = gsl_ran_gaussian_pdf((r/rin)*sin(th), wid/sqrt(2*log(2)));

      // Narayan '12, Penna '12 conditions
      } else if (mad_type == NARAYAN) { 

        // Former uses rstart=25, rend=810, lam_B=25
        double uc = uu_av - uu_end;
        double ucm = uu_plane_av - uu_end;
        q = pow(sin(th),3)*(uc/(ucm+SMALL) - 0.2) / 0.8;

        //Exclude q outside torus and large q resulting from division by SMALL
        if ( r > rBend || r < rBstart || fabs(q) > 1.e2 ) q = 0;

        /////////////////////////////////////////////////////////////////////////////////////
        //if (q != 0 && th > M_PI/2-0.1 && th < M_PI/2+0.1) printf("q_mid is %.10e\n", q);
        /////////////////////////////////////////////////////////////////////////////////////

      } else {
        printf("MAD = %i not supported!\n", mad_type);
        exit(-1);
      }
    //////////////////////////
    // TODO How about 2D?
    //////////////////////////
    // 2D mad parameter
    } else { 
      q = rho_av/rhomax;
    }

    // vector potentials
    A[i][j] = 0.;
    if (q > 0.) {

      // Narayan limit uses alternating loops
      if (mad_type == NARAYAN) { 
        double lam_B = 25;
        double flux_correction = sin( 1/lam_B * (pow(r,2./3) + 15./8*pow(r,-2./5) - pow(rBstart,2./3) - 15./8*pow(rBstart,-2./5)));
        double q_mod = q*flux_correction;
        A[i][j] = q_mod;
      } else {
        A[i][j] = q;
      }
    }
  } // ZSLOOP

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
    S->P[B1][k][j][i] = -(A[i][j] - A[i][j + 1]
	          + A[i + 1][j] - A[i + 1][j + 1]) /
	          (2. * dx[2] * G->gdet[CENT][j][i]);
    S->P[B2][k][j][i] = (A[i][j] + A[i][j + 1]
	          - A[i + 1][j] - A[i + 1][j + 1]) /
	          (2. * dx[1] * G->gdet[CENT][j][i]);

    // no theta direction magnetic field
    // TODO: good to research on B-theta
    S->P[B3][k][j][i] = 0.;

    // get ucon, ucov, bcon, bcov
    get_state(G, S, i, j, k, CENT);

    // Non-narayan MAD type
    if ((r > rBstart && r < rBend) || mad_type != NARAYAN) {
      // find bsquare at each grid
      double bsq_ij = bsq_calc(S, i, j, k);

      // limit bsquare
      if (bsq_ij > bsq_max) bsq_max = bsq_ij;

      // find magnetization
      double beta_ij = (gam - 1.)*(S->P[UU][k][j][i])/(0.5*(bsq_ij+SMALL)) ;

      // limit magnetization
      if(beta_ij < beta_min) beta_min = beta_ij ;
    }
  }

  // find maximum bsqaure and magnetization
  bsq_max = mpi_max(bsq_max);
  beta_min = mpi_min(beta_min);

  // find maximum internal energy?
  double umax_plane = 0;
  for (int i = 0; i < N1TOT; i++) {
      double X[NDIM];
      coord(i,NG,NG,CORN,X);
      double r, th;
      bl_coord(X,&r,&th);
      if ((r > rBstart && r < rBend) || mad_type != NARAYAN) {
        if (uu_plane[i] > umax_plane) umax_plane = uu_plane[i];
      }
  }

  // report maximum internal energy, bsqaure, and magnetization, and find normalization parameter
  double norm = 0;
  if (!maxr_normalization) {
    // Find maximum pressure to maximum bsquare 
    double beta_act = (gam - 1.) * umax / (0.5 * bsq_max);

    ////////////////////////////////////////////////////////////////
    // In plane only
    //double beta_act = (gam - 1.) * umax_plane / (0.5 * bsq_max);
    ////////////////////////////////////////////////////////////////

    // print out
    LOGN("Umax is %.10e", umax);
    LOGN("bsq_max is %.10e", bsq_max);
    LOGN("beta is %.10e", beta_act);

    // normalization factor
    norm = sqrt(beta_act / beta);

    // or else, normalize using minimum beta
  } else {
    // Beta_min = 100 normalization
    LOGN("Min beta in torus is %f", beta_min);
    norm = sqrt(beta_min / beta) ;
  }

  // Apply normalization to magnetic fields, B3 is uniformly 0
  LOGN("Normalization is %f\n", norm);
  ZLOOP {
    S->P[B1][k][j][i] *= norm ;
    S->P[B2][k][j][i] *= norm ;
  }
  printf("%.12e\n", norm);
  ///////////////////////////////////////////
  // TODO do this only for BHflux > SMALL
  ///////////////////////////////////////////
  // This adds a central flux based on specifying some BHflux
  // Initialize a net magnetic field inside the initial torus
  ZSLOOP(0, 0, 0, N2, 0, N1) {

    // find r and theta
    double X[NDIM];
    coord(i,j,k,CORN,X);
    double r,th;
    bl_coord(X, &r, &th);

    // magnetic vector potential
    A[i][j] = 0.;

    // find r, z, ...
    double x = r*sin(th);
    double z = r*cos(th);
    double a_hyp = 20.;
    double b_hyp = 60.;
    double x_hyp = a_hyp*sqrt(1. + pow(z/b_hyp,2));

    // add vector potentials
    double q = (pow(x,2) - pow(x_hyp,2))/pow(x_hyp,2);
    if (x < x_hyp) {
      A[i][j] = 10.*q;
    }
  }

  // Evaluate net flux
  double Phi_proc = 0.;
  ISLOOP(5, N1-1) {
    JSLOOP(0, N2-1) {

      // what is this?
      int jglobal = j - NG + global_start[1];

      //////////////////////
      //int j = N2/2+NG;
      //////////////////////

      // what is this?
      int k = NG;
      if (jglobal == N2TOT / 2) {

        // find r and theta
        double X[NDIM];
        coord(i, j, k, CENT, X);
        double r, th;
        bl_coord(X, &r, &th);
        
        // for position smaller than inner torus radius
        if (r < rin) {
          // averageing to get B2?
          double B2net = (A[i][j] + A[i][j + 1] - A[i + 1][j] - A[i + 1][j + 1]);

          /////////////////////////////////////
          // / (2.*dx[1]*G->gdet[CENT][j][i]);
          /////////////////////////////////////

          // get phi_proc
          Phi_proc += fabs(B2net) * M_PI / N3CPU; // * 2.*dx[1]*G->gdet[CENT][j][i]
        }
      }
    }
  }

  //If left bound in X1.  Note different convention from bhlight!
  if (global_start[0] == 0) {
    JSLOOP(0, N2/2-1) {
      //what is this?
      int i = 5 + NG;

      // get B1
      double B1net = -(A[i][j] - A[i][j+1] + A[i+1][j] - A[i+1][j+1]); // /(2.*dx[2]*G->gdet[CENT][j][i]);

      // phi_proc?
      Phi_proc += fabs(B1net)*M_PI/N3CPU;  // * 2.*dx[2]*G->gdet[CENT][j][i]
    }
  }
  double Phi = mpi_reduce(Phi_proc);

  //normalization parameter
  norm = BHflux/(Phi + SMALL);

  // Add extra B-field by B = curl of A
  ZLOOP {
    // Flux-ct
    S->P[B1][k][j][i] += -norm
        * (A[i][j] - A[i][j + 1] + A[i + 1][j] - A[i + 1][j + 1])
        / (2. * dx[2] * G->gdet[CENT][j][i]);
    S->P[B2][k][j][i] += norm
        * (A[i][j] + A[i][j + 1] - A[i + 1][j] - A[i + 1][j + 1])
        / (2. * dx[1] * G->gdet[CENT][j][i]);
  }

  // initialize electronic variables
#if ELECTRONS
  init_electrons(G,S);
#endif 
  printf("%.12e\n", norm);
  // apply floors and ceilings
  fixup(G, S);

  // Enforce boundary conditions
  set_bounds(G, S);

  // printout
  LOG("Finished init()");

}

//******************************************************************************

// This function calculates specific the angular momentum of the
// Fishbone-Moncrief solution in the midplane,
// as a function of radius.
// (see Fishbone & Moncrief eqn. 3.8)
// It improves on (3.8) by requiring no sign changes
// for co-rotating (a > 0) vs counter-rotating (a < 0)
// disks.
double lfish_calc(double r)
{
  return (((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
     ((-2. * a * r *
       (pow(a, 2) - 2. * a * sqrt(r) +
        pow(r,
      2))) / sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
      ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) + pow(a, 2) *
      (2. + r))) / sqrt(1 + (2. * a) / pow (r, 1.5) - 3. / r)))
    / (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
       (pow(a, 2) + (-2. + r) * r))
      );
}

