/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR BRIO-WU SHOCK TUBE                                  *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include "hdf5_utils.h"

/*****************************************************************************/

// load from the parameter file //
void set_problem_params() {
  /* do nothing here */
}

/*****************************************************************************/

/* write to dump */
void save_problem_data(hid_t string_type){
  hdf5_write_single_val("briowu", "PROB", string_type);
}

/*****************************************************************************/

/* set initial condition */
void init(struct GridGeom *G, struct FluidState *S)
{

  /* set grid */
  set_grid(G);
  LOG("Set grid");

  /* find coordinates */
  double X[NDIM];

  /* loop over */
  ZLOOPALL {

    /* find coordiantes */
    coord(i, j, k, CENT, X);

    /* right state */
    if(X[1] > 0.5) {
      S->P[RHO][k][j][i] = 1.0;
      S->P[UU][k][j][i] = (2.0/3.0)*1e-6/(gam - 1.0);

    /* left state */  
    } else {
      S->P[RHO][k][j][i] = 10.0;
      S->P[UU][k][j][i] = (40.0/3.0)/(gam - 1.0);
    }

    /* remaining primitive variables */
    S->P[U1][k][j][i] = 0.0;
    S->P[U2][k][j][i] = 0.0;
    S->P[U3][k][j][i] = 0.0;
    S->P[B1][k][j][i] = 0.0;
    S->P[B2][k][j][i] = 0.0;
    S->P[B3][k][j][i] = 0.0;
  }

  // Enforce boundary conditions
  set_bounds(G, S);

}

/*****************************************************************************/
