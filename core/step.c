//******************************************************************************
//*                                                                            *
//* STEP.C                                                                     *
//*                                                                            *
//* ADVANCES SIMULATION BY ONE TIMESTEP                                        *
//*                                                                            *
//******************************************************************************
// Leon's comment: here, the system of hyperbolic PDE is evolved by the 
// 2nd order RK method using method of lines. We can write the time-evolution
// as: dy/dt = f(y,t). The first step is to evolve the equations by dt/2
// i.e., y* = y0 + f(y0, t0)*(dt/2), t* = t0 + dt/2
// Then, the variables y at the new time step is evolved as
// y' = y0 + dt*f(y*, t*)
// y is the conservative variables U, and f is the flux-difference + sources

//import headers
#include "decs.h"
#include "positrons.h"

// delcare fuctions
double advance_fluid(struct GridGeom *G, struct FluidState *Si, struct FluidState *Ss, struct FluidState *Sf, double Dt);

//**************************************************************************************************************************

//advance system of equations in time
void step(struct GridGeom *G, struct FluidState *S)
{
  //declare data structure
  static struct FluidState *Stmp;
  static struct FluidState *Ssave;

  //declare arrays
  static int first_call = 1;
  if (first_call) {
    Stmp = calloc(1,sizeof(struct FluidState));
    Ssave = calloc(1,sizeof(struct FluidState));
    first_call = 0;
  }

  // backup primitive variables 
  ////////////////////////////////////////////////////////////////////
  // Need both P_n and P_n+1 to calculate current
  // Work around ICC 18.0.2 bug in assigning to pointers to structs
  // TODO use pointer tricks to avoid deep copy on both compilers
  ////////////////////////////////////////////////////////////////////
#if INTEL_WORKAROUND
  memcpy(&(Ssave->P),&(S->P),sizeof(GridPrim));
#else
#pragma omp parallel for simd collapse(3)
  PLOOP ZLOOPALL Ssave->P[ip][k][j][i] = S->P[ip][k][j][i];
#endif

  //print out
  LOGN("Step %d",nstep);
  FLAG("Start step");
  ///////////////////////////////////////////////////
  // TODO add back well-named flags /after/ events
  ///////////////////////////////////////////////////

  /*-------------------------------------------------------------------------*/
  // Predictor setup, here, Stmp is empty, but then when passed to
  // advanced_fluid, it will get compies from S, Stmp then becomes the 
  // variables at half time step, y*, as stated above

  // evolve by half dt
  advance_fluid(G, S, S, Stmp, 0.5*dt);
  FLAG("Advance Fluid Tmp");

  // find electronic source terms
#if ELECTRONS
  heat_electrons(G, S, Stmp);
  FLAG("Heat Electrons Tmp");
#endif

  // Leon's patch, pair production
  /* do only if the flag for pair production is on */
#if POSITRONS && PAIRS
  timer_start(TIMER_POSITRON);
  pair_production(G, S, Stmp, 0.5*dt);
  FLAG("Pair Production Tmp");
  timer_stop(TIMER_POSITRON);
#endif

  // Set floor values to primitive variables 
  fixup(G, Stmp);
  FLAG("Fixup Tmp");
#if ELECTRONS
  fixup_electrons(Stmp);
  FLAG("Fixup e- Tmp");
#endif

  // set boundary conditions 
  ////////////////////////////////////////////////////////////////////
  // Need an MPI call _before_ fixup_utop to obtain correct pflags
  ////////////////////////////////////////////////////////////////////
  set_bounds(G, Stmp);
  FLAG("First bounds Tmp");

  //replace bad points (failed convergence) with trilinear interpolations 
  fixup_utoprim(G, Stmp);
  FLAG("Fixup U_to_P Tmp");

  //after that, set boundary conditions again
  set_bounds(G, Stmp);
  FLAG("Second bounds Tmp");
  
  /*-------------------------------------------------------------------------*/
  // Corrector step, here, Stmp is the half time-step variables y*

  // evolve by dt
  double ndt = advance_fluid(G, S, Stmp, S, dt);
  FLAG("Advance Fluid Full");

  // find electronic source terms
#if ELECTRONS
  heat_electrons(G, Stmp, S);
  FLAG("Heat Electrons Full");
#endif

  // Leon's patch, pair production
  /* do only if the flag for pair production is on */
#if POSITRONS && PAIRS
  timer_start(TIMER_POSITRON);
  pair_production(G, Stmp, S, dt);
  FLAG("Pair Production Tmp");
  timer_stop(TIMER_POSITRON);
#endif

  // Set floor values to primitive variables 
  fixup(G, S);
  FLAG("Fixup Full");
#if ELECTRONS
  fixup_electrons(S);
  FLAG("Fixup e- Full");
#endif

  // set boundary conditions 
  set_bounds(G, S);
  FLAG("First bounds Full");

  //replace bad points (failed convergence) with trilinear interpolations 
  fixup_utoprim(G, S);
  FLAG("Fixup U_to_P Full");

  //after that, set boundary conditions again
  set_bounds(G, S);
  FLAG("Second bounds Full");
  
  /*-------------------------------------------------------------------------*/

  // Increment time
  t += dt;

  // If we're dumping this step, calculate the current
  if (t >= tdump) {
    current_calc(G, S, Ssave, dt);
  }

  // New dt proxy to choose fluid or light timestep
  double max_dt = 0, fake_dt = 0;
#if STATIC_TIMESTEP
  if(DEBUG) fake_dt = mpi_min(ndt);
  max_dt = cour*dt_light;
#else
  if(DEBUG) fake_dt = cour*dt_light;
  max_dt = mpi_min(ndt);
#endif

  // Set next timestep
  if (max_dt > SAFE * dt) {
    dt = SAFE * dt;
  } else {
    dt = max_dt;
  }

  // print out
  LOGN("dt would have been %f",fake_dt);
  LOGN("Instead it is %f",dt);

}

//**************************************************************************************************************************

//find flux gradient, source term, and then update conservative variables by dt
inline double advance_fluid(struct GridGeom *G, struct FluidState *Si, struct FluidState *Ss, struct FluidState *Sf, double Dt)
{
  //declare
  static GridPrim *dU;
  static struct FluidFlux *F;
  
  //assign memories
  static int firstc = 1;
  if (firstc) {
    dU = calloc(1,sizeof(GridPrim));
    F = calloc(1,sizeof(struct FluidFlux));
    firstc = 0;
  }

  // backup primitive variables
  ////////////////////////////////////////////////////////////////////
  // Work around ICC 18.0.2 bug in assigning to pointers to structs
  ////////////////////////////////////////////////////////////////////
#if INTEL_WORKAROUND
  memcpy(&(Sf->P),&(Si->P),sizeof(GridPrim));
#else
#pragma omp parallel for simd collapse(3)
  PLOOP ZLOOPALL Sf->P[ip][k][j][i] = Si->P[ip][k][j][i];
#endif

  //get fluxes
  double ndt = get_flux(G, Ss, F);

//////////////////////////////////
//  update_f(F, dU);
//  FLAG("Got initial fluxes");
//////////////////////////////////

  //fix fluxes values
#if METRIC == MKS
  fix_flux(F);
#endif

//////////////////////////////////
//  update_f(F, dU);
//  FLAG("Fixed Flux");
//////////////////////////////////

  //Constrained transport for B
  flux_ct(F);

//////////////////////////////////
//  update_f(F, dU);
//  FLAG("CT Step");
//////////////////////////////////

  // Flux diagnostic globals
  ///////////////////////////////////////////////////
  // TODO don't compute every step, only for logs?
  ///////////////////////////////////////////////////
  diag_flux(F);

//////////////////////////////////
//  update_f(F, dU);
//  FLAG("Flux Diags");
//////////////////////////////////

  // Get conservative variables, and fluid source terms
  timer_start(TIMER_UPDATE_U);
  get_state_vec(G, Ss, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1);
  get_fluid_source(G, Ss, dU);

  // Find conservative variables for the last time step'
  /////////////////////////////////////////////
  // TODO skip this call if Si,Ss are aliased
  /////////////////////////////////////////////
  get_state_vec(G, Si, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1);
  prim_to_flux_vec(G, Si, 0, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1, Si->U);

//////////////////////////////////
//  update_f(F, dU);
//  FLAG("Fixed flux. Got Si->U");
//////////////////////////////////

  // update conservative variables 
#pragma omp parallel for collapse(3)
  PLOOP ZLOOP {
    Sf->U[ip][k][j][i] = Si->U[ip][k][j][i] +
      Dt*((F->X1[ip][k][j][i] - F->X1[ip][k][j][i+1])/dx[1] +
          (F->X2[ip][k][j][i] - F->X2[ip][k][j+1][i])/dx[2] +
          (F->X3[ip][k][j][i] - F->X3[ip][k+1][j][i])/dx[3] +
          (*dU)[ip][k][j][i]);
  }
  timer_stop(TIMER_UPDATE_U);

  //FLAG("Got Sf->U");

  /******************************************************************/
  // Leon's patch, add cooling here //
  // Note that the final state Sf has already include addition 
  // from the previous state Si and the flux difference
  // we pass only Ss but not Si because the RHS of the ODE dy/dt
  // depends only either the previous step y0 or intermediate y*
  #if COOLING
    timer_start(TIMER_COOLING);
    rad_cooling(G, Ss, Sf, Dt);
    FLAG("Radiative Cooling Tmp");
    timer_stop(TIMER_COOLING);
  #endif

  //convert from conservative to primitive variables
  timer_start(TIMER_U_TO_P);
#pragma omp parallel for collapse(3)
  ZLOOP {
    pflag[k][j][i] = U_to_P(G, Sf, i, j, k, CENT);
    ////////////////////////////////////////////////////////////////////
    // This is too annoying even for debug
    //if (pflag[k][j][i] != 0) LOGN("Pflag is %d\n", pflag[k][j][i]);
    ////////////////////////////////////////////////////////////////////
  }
  timer_stop(TIMER_U_TO_P);

  ////////////////////////////////////////////////////////////////////
  //FLAG("Got Sf->P");
  // Not complete without setting four-vectors
  // Done /before/ each call
  //get_state_vec(G, Sf, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1);
  ////////////////////////////////////////////////////////////////////

  // save error flag in u to p subroutine
#pragma omp parallel for simd collapse(2)
  ZLOOPALL {
    fail_save[k][j][i] = pflag[k][j][i];
  }

  //output
  return ndt;
}
