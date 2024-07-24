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
#include "positrons.h"
#include "cooling.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>

// compile only if poistrons flag is on //
#if POSITRONS

//******************************************************************************

// set the unit convestion between code and cgs // 
void set_units () {

  /* black hole mass in cgs */ 
  Mbh = mbh*MSUN;

  /* Length and time scales */
  L_unit = GNEWT*Mbh/(CL*CL);
  T_unit = L_unit/CL;

  /* Now set the remaining unit */
  RHO_unit = M_unit*pow(L_unit,-3.);
  U_unit = RHO_unit*CL*CL;

}

//******************************************************************************

//initialize positrons variables
void init_positrons(struct GridGeom *G, struct FluidState *S)
{
  
  // Set positron mass to its floor values //
  ZLOOPALL {
    S->P[RPL][k][j][i] = ZMIN*ME_MP*S->P[RHO][k][j][i];
  }

  // Necessary?  Usually called right afterward
  set_bounds(G, S);

}

//******************************************************************************

// Compute net pair production rate //
void pair_production(struct GridGeom *G, struct FluidState *Ss, struct FluidState *Sf, double dt_step)
{

  /* Then, compute the pair production rate */
#pragma omp parallel for collapse(3)
  ZLOOP {
    pair_production_1zone(G, Ss, Sf, i, j, k, dt_step);
  }

}

//******************************************************************************

// compute pair production rate per grid cells //
inline void pair_production_1zone(struct GridGeom *G, struct FluidState *Ss, struct FluidState *Sf, int i, int j, int k , double dt_step)
{

  /***********************************************************************/

  // coordinate //
  double rad, theta, X[NDIM];

  // coordinate radius //
  coord(i, j, k, CENT, X);
  bl_coord(X, &rad, &theta);

  /***********************************************************************/
  // number density, all in c.g.s unit //

  // time step in second !
  double dt_real = dt_step*T_unit;

  /***********************************************************************/
  
  // calculate proton and positron numbe density, all in cgs //
  double nprot = Ss->P[RHO][k][j][i]*RHO_unit/MP;
  double npost = Ss->P[RPL][k][j][i]*RHO_unit/ME;

  /***********************************************************************/

  // optical depth and scale height //
  double tau_depth = 0;
  double h_th = 0;

  // direct integration of optical depth and scale height //
#if COMPUTE == DIRECT
  /*--------------------------------------------------------------------------------------------*/
  // local variable storing integrand //
  int m_start, m_end;
  double np, n_p, nt, th_loc;
  double upper = 0;
  double lower = 0; 
  double dummy;

  // integrate optical depth and scale height //
  // they are all in the C.G.S unit //
  if(theta < M_PI_2) {
    // upper hemisphere //
    m_start = NG;
    m_end = j;
  } else {
    // lower hemispehere
    m_start = j;
    m_end = N2 + NG;
  }

  /* sum over */
  for (int m = m_start; m <= m_end; m++) {
    coord(i, m, k, CENT, X);
    bl_coord(X, &dummy, &th_loc);
    np = Ss->P[RHO][k][m][i]*RHO_unit/MP;
    n_p = Ss->P[RPL][k][m][i]*RHO_unit/ME;
    nt = 2*(np + n_p);
    tau_depth = tau_depth + nt*sqrt(G->gcov[CENT][2][2][m][i])*dx[2]*L_unit*sigma_t;
    upper = upper + nt*fabs(th_loc - M_PI_2)*G->gdet[CENT][m][i]*sqrt(G->gcov[CENT][2][2][m][i])*dx[2];
    lower = lower + nt*G->gdet[CENT][m][i]*sqrt(G->gcov[CENT][2][2][m][i])*dx[2];
  }

  // scale height, remember change to CGS //
  h_th = upper/lower*L_unit;

  /*--------------------------------------------------------------------------------------------*/
  // gaussian approximation //
#elif COMPUTE == GAUSSIAN
  // Note: need to find asymtopic 

  // local variables // 
  int j_mid;
  double np, n_p, nt;
  double t1, t2;
  double dummy;

  // calculate mid plane index //
  j_mid = (int)((NG + N2 + NG)/2);

  // mid plane number density //
  np = Ss->P[RHO][k][j_mid][i]*RHO_unit/MP;
  n_p = Ss->P[RPL][k][j_mid][i]*RHO_unit/ME;
  nt = 2*(np + n_p); 

  // upper atmosphere //
  if(theta < M_PI_2) {
    t1 = - M_PI_2/h_r;
    t2 = (theta - M_PI_2)/h_r;
  } else {
    t1 = (theta - M_PI_2)/h_r;
    t2 = M_PI_2/h_r;
  }
  dummy = gsl_sf_erf(t2)-gsl_sf_erf(t1);
  tau_depth = nt*(rad*h_r*L_unit)*(sqrt2*0.5*sqrt(M_PI))*(dummy)*sigma_t;
  if(fabs(t1) > 5 && fabs(t2) > 5){
    dummy = -sqrt(M_PI)*(exp(t1*t1-t2*t2) - 1)/(exp(t1*t1-t2*t2)*series_asym(t2) - series_asym(t1));
    h_th = sqrt(2/M_PI)*(rad*h_r*L_unit)*dummy;
  } else {
    h_th = sqrt(2/M_PI)*(rad*h_r*L_unit)*(exp(-t2*t2) - exp(-t1*t1))/(dummy);
  } 
  h_th = fabs(h_th);
  tau_depth = fmax(tau_depth, SMALL);

  /*--------------------------------------------------------------------------------------------*/
#endif

  /***********************************************************************/

  // dimensionless temperature and postiron fraction //
  double thetae = KBOL*t_elec/(ME*CL*CL); 
  double zfrac = npost/nprot; 

  /***********************************************************************/
  /* now calculate pair production rate */

  // net pair production rate, note the rate is in the CGS unit!!! //
  double net_rate = ndot_net(zfrac, tau_depth, nprot, thetae, h_th);
  
  /* do these steps only if the production rate is non-zero */
  if(fabs(net_rate) > 0) {    

    // quality factor //
    double qfac = fabs(npost/net_rate);

    /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
    /* if the source term is too steep, implement implicit solver 
    /* Basically a root finding, so use bisection mtehod */
    if(dt_real > q_alpha*qfac) {

      /* print out */
      printf("net_rate too steep %d %d %d\n", i, j, k);

      /* left state */
      double zl = zfrac;
      double ndotl = ndot_net(zl, tau_depth, nprot, thetae, h_th);
      double fl = (zl - zfrac) - dt_real*ndotl/nprot;

      /* right state */
      int o;
      double zr, ndotr, fr, steps;
      if(net_rate > 0) {
        steps = 10;
      } else{
        steps = 0.1;
      }
      zr = zl;
      for (o = 0; o < 999; o++) {
        ndotr = ndot_net(zr, tau_depth, nprot, thetae, h_th);
        fr = (zr - zfrac) - dt_real*ndotr/nprot;
        if(fr*fl <0) break;
        zr = zr*steps;
      }

      /* exit condition */
      if(o == 999 || isnan(fr)) {
        printf("Failure in implicit method\n");
        return;
      }

      /* define the center state */
      double zcen, fcen, zcen_old, ndotcen;

      /* bisection method counting */
      int count;

      /* now iterate until converges */
      for (count = 0; count < 99999; count++) {

        /* backup */
        if(count > 0) {
          zcen_old = zcen;
        }

        /* center state */
        zcen = 0.5*(zl + zr);
        ndotcen = ndot_net(zcen, tau_depth, nprot, thetae, h_th);
        fcen = (zcen - zfrac) - dt_real*ndotcen/nprot;
        
        /* check the sign */
        if (fl*fcen > 0) {
          zl = zcen;
        } else if (fr*fcen > 0) {
          zr = zcen;
        }

        /* determine if need to exit */
        if(count > 0) {
          if(fabs(1.0 - zcen_old/zcen) < bisects) {
            break;
          }
        }
      }

      /* exit if no solution, and print out error */
      if(count == 99999) {
        printf("No solution\n");
        return;
      }

      /* assign new positron mass, remember to convert back to code unit !!! */
      npost = zcen*nprot;
      
    /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
    /* otherwise, march forward by time */
    } else {

      // positron mass production rate, need to convert to code unit!!! //
      npost = npost + net_rate*dt_real;
    
    }
    /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

    // update positron mass //
    Sf->P[RPL][k][j][i] = npost*(ME/RHO_unit);
  }
}

//*------------------------------------------------------------------------------------------------------------------------*//
//
// Now, the remaining of the code are all about computing pair production rates
// I use a long comment to seperate out the main body of the this code for 
// better documentation and reading
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//******************************************************************************

/* net pair production rate */
inline double ndot_net(double zfrac, double taut, double nprot, double theta, double r_size) {
  double xm = find_xm(zfrac, taut, nprot, theta);
  double ndotbr = get_ndotbr(zfrac, theta, xm, nprot);
  double y1 = comptony1(xm, taut, theta);
  double fb = fbrem(y1);
  double n1 = flatn1(xm, theta, y1);
  double ng = ngamma(xm, taut, theta, y1, zfrac, nprot, fb, ndotbr, r_size);
  double nc = ncdot(ng, theta, nprot, zfrac, n1, fb, ndotbr);
  double na = nadot(zfrac, nprot, theta);
  return nc - na;
}

//******************************************************************************

/* Total pair production rate due to photon-photo, photon-particle collision */
inline double ncdot(double ngamma, double theta, double nprot, double z, double n1, double fb, double ndotbr) {
  double ndotww = get_ndotww(ngamma, theta);
  double ndotwp = get_ndotwp(ngamma, nprot, theta);
  double ndotwe = get_ndotwe(ngamma, nprot, z, theta);
  double ndotwf = get_ndotwf(n1, ngamma, theta);
  double ndotee = get_ndotee(nprot, z, theta);
  double out = ndotee + ndotww + ndotwp + ndotwe + ndotwf ;
  return out;
}

//******************************************************************************

// functions for computing e-e pair production rates //
inline double get_ndotee(double nprot, double z, double theta) {
  double ndot;
  if(theta <= 1e2) {
    ndot = 2e-4*pow(theta,1.5)*exp(-2.0/theta)*(1.0+0.015*theta);
  } else if(theta >= 1e2) {
    ndot = (112.0/27.0/M_PI)*(alphaf*alphaf)*pow(log(theta),3)/(1.0 + 0.058/theta);
  } else {
    ndot = 0.0;
  }
  ndot = ndot*CL*RE*RE*(nprot*(1.0 + z))*(nprot*(1.0 + z));
  return ndot;
}

//******************************************************************************

/* This is for computing the pair annhilation rate */
inline double nadot(double z, double nprot, double theta) {
  double out = (3.0/8.0)*sigma_t*CL*(nprot*nprot*z*(z+1.0))/(1.0+2.0*theta*theta/log(1.12*theta+1.3));
  return out;
}

//*------------------------------------------------------------------------------------------------------------------------*//
//
// These sections are the analytic formula for computing photon-photon or photon-particle processes
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//******************************************************************************

/* wein - wein, Svensson 1984, White & Lightman 1989 */
inline double get_ndotww(double ngamma, double theta) {
  double out;
  if(theta <= 1.0) {
    out = (0.125)*M_PI*M_PI*exp(-2.0/theta)*(1.0 + 2.88*pow(theta,0.934))/pow(theta,3);
  } else {
    out = (0.5)*M_PI*log(2.0*eta*theta + 0.38)/pow(theta,2);
  }
  out = out*CL*RE*RE*ngamma*ngamma;
  return out;
}

//******************************************************************************

/* wein - proton, Svensson 1984, White & Lightman 1989 */
inline double get_ndotwp(double ngamma, double nprot, double theta) {
  double out;
  if(theta <= 2.0) {
    out = M_PI*theta*exp(-2.0/theta)/(1.0 + 0.9*theta);
  } else {
    out = (28.0/9.0)*log(2.0*eta*theta + 1.7) - 92.0/27.0;
  }
  out = out*alphaf*CL*RE*RE*ngamma*nprot;
  return out;
}

//******************************************************************************

/* wein - lepton, Svensson 1984, White & Lightman 1989 */
inline double get_ndotwe(double ngamma, double nprot, double z, double theta) {
  double out;
  if(theta <= 0.18) {
    out = (4.0*M_PI/27.0)*exp(-2.0/theta)*(1.0 + 27.1*pow(theta,0.949));
  } else if(theta >= 2) {
    out = (56.0/9.0*log(2.0*eta*theta) - 8.0/27.0)/(1.0 + 0.5/theta);
  } else {
    out = (4.0*M_PI/27.0)*exp(-2.0/theta)*16.1*pow(theta,0.541);
  }
  out = out*alphaf*CL*RE*RE*ngamma*(2.0*z+1.0)*nprot;
  return out;
}

//******************************************************************************

/* wein - flat, Svensson 1984, White & Lightman 1989 */
inline double get_ndotwf(double n1, double ngamma, double theta) {
  double out;
  out = CL*RE*RE*n1*ngamma*M_PI*M_PI/4.0*exp(-1.0/theta);
  return out;
}

//*------------------------------------------------------------------------------------------------------------------------*//
//
// These sections are the analytic formula for computing total photon emissivity
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//******************************************************************************

/* integration of ln(A*theta/x)/x */
inline double integrate_log(double A, double theta, double xm) {
  double x1 = xm;
  double x2 = theta;
  double f1 = log(A*theta/x1)*log(x1) + 0.5*log(x1)*log(x1);
  double f2 = log(A*theta/x2)*log(x2) + 0.5*log(x2)*log(x2);
  return f2 - f1;
}

//******************************************************************************

/* electron-proton bremsstrahlung total rate */
inline double rate_ep(double z, double nprot, double theta, double xm) {
  double A = 4*eta*(1+3.42*theta);
  double integrate = integrate_log(A, theta, xm);
  double out = (1+2*z)*(1+2*theta+2*theta*theta)*integrate;
  return out;
}

//******************************************************************************

/* electron-electron bremsstrahlung */
inline double rate_ee(double z, double nprot, double theta, double xm) {
  double A = 4*eta*(11.2+10.4*theta*theta);
  double integrate = integrate_log(A, theta, xm);
  double out = (z*z+(1+z)*(1+z))*(3*sqrt2/5*theta+2*theta*theta)*integrate;
  return out;
}

//******************************************************************************

/* electron-positron bremsstrahlung */
inline double rate_pm(double z, double nprot, double theta, double xm) {
  double A = 4*eta*(1+10.4*theta*theta);
  double integrate = integrate_log(A, theta, xm);
  double out = z*(1+z)*2*(sqrt2+2*theta+2*theta*theta)*integrate;
  return out; 
}

//******************************************************************************

/* bremsstrahlung production rate */
inline double get_ndotbr(double z, double theta, double xm, double nprot) {
  double thetam1 = 1/theta;
  double corr;
  if(thetam1 < 500) {
    double k2 = gsl_sf_bessel_Kn(2, thetam1);
    corr = exp(thetam1)*k2;
  } else {
    corr = sqrt(M_PI_2/thetam1);
  }
  double factor = (16/3)*(alphaf)*(CL)*(RE*RE)*(nprot*nprot)/(corr);
  double ep = rate_ep(z, nprot, theta, xm);
  double ee = rate_ee(z, nprot, theta, xm);
  double pm = rate_pm(z, nprot, theta, xm);
  double out = factor*(ep + ee + pm);
  return out;
}

//*------------------------------------------------------------------------------------------------------------------------*//
//
// These sections are the analytic formula for computing photon emissivity per frequency
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//******************************************************************************

/* black body number density */ 
inline double nbb(double x, double theta) {
  double tmp = x/theta;
  double out;
  if(tmp < 1e-5) {
    out = (x*theta)/(lambda_c*lambda_c*lambda_c)/(M_PI*M_PI);
  } else {
    out = (x*x)/(lambda_c*lambda_c*lambda_c)/(M_PI*M_PI)/(exp(tmp) - 1);
  }
  return out;
}

//******************************************************************************

/* d(n0)/dt factor */
inline double n0dot(double x, double nprot, double theta) {
  double thetam1 = 1/theta;
  double corr;
  if(thetam1 < 500) {
    double k2 = gsl_sf_bessel_Kn(2, thetam1);
    corr = exp(thetam1)*k2;
  } else {
    corr = sqrt(M_PI_2/thetam1);
  }
  double out = (16/3)*(alphaf)*(CL)*(RE*RE)*(nprot*nprot)/(corr)*(exp(-x*thetam1)/x);
  return out;
}

//******************************************************************************

/* electron-proton bremsstrahlung */
inline double ndotep(double x, double z, double nprot, double theta, double dn0dt) {
  double out = (1+2*z)*log(4*eta*(1+3.42*theta)*(theta/x))*(1+2*theta+2*theta*theta)*dn0dt;
  return out;
}

//******************************************************************************

/* electron-electron bremsstrahlung */
inline double ndotee(double x, double z, double nprot, double theta, double dn0dt) {
  double out = (z*z+(1+z)*(1+z))*log(4*eta*(11.2+10.4*theta*theta)*(theta/x))*(3*sqrt2/5*theta+2*theta*theta)*dn0dt;
  return out;
}

//******************************************************************************

/* electron-positron bremsstrahlung */
inline double ndotpm(double x, double z, double nprot, double theta, double dn0dt) {
  double out = z*(1+z)*log(4*eta*(1+10.4*theta*theta)*(theta/x))*2*(sqrt2+2*theta+2*theta*theta)*dn0dt;
  return out;
}

//******************************************************************************

/* bremsstrahlung absorbtion coefficient */
inline double brem_abs(double x, double z, double nprot, double theta) {
  double bb = nbb(x, theta);
  double dn0dt = n0dot(x, nprot, theta);
  double ep = ndotep(x, z, nprot, theta, dn0dt);
  double ee = ndotee(x, z, nprot, theta, dn0dt);
  double pm = ndotpm(x, z, nprot, theta, dn0dt);
  double out = (ep + ee + pm)/(CL*bb);
  return out;
}

//*------------------------------------------------------------------------------------------------------------------------*//
//
// These sections related to the "radiative transfer"
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//******************************************************************************

// find photon frequency below witch the local spectrum is black body //
inline double find_xm(double z, double tau, double nprot, double theta) {
    
  /* backup */
  double xc_old;

  /* set the LHS of the root */
  double at = (2*z+1)*nprot*sigma_t;
  double lhs = at*(1+zeta*(tau*tau)*fmin(1,8*theta))/(tau*(1+zeta*tau));
    
  /* initial guess */
  double xl = -50, xr = log10(700);
  double xlp = pow(10.0, xl)*theta, xrp = pow(10.0, xr)*theta;
  double fl = brem_abs(xlp, z, nprot, theta) - lhs;
  double fr = brem_abs(xrp, z, nprot, theta) - lhs;
  double xc = 0.5*(xl+xr), xcp = pow(10.0, xc)*theta;
   
  /* poor initial guess, exit */
  if(fl*fr > 0) {
    if(mpi_io_proc()) printf("poor initial guess xm");
    exit(0);
  }

  /* continue */
  double fc = brem_abs(xcp, z, nprot, theta) - lhs;

  /* main loop */
  int n; 
  for (n = 0; n < 99999; n++) {
    xc_old = xc;
    if(fc*fl > 0) {
      fl = fc, xl = xc;
    } else if(fc*fr > 0){
      fr = fc, xr = xc;
    }
    xc = 0.5*(xl+xr), xcp = pow(10.0, xc)*theta;
    fc = brem_abs(xcp, z, nprot, theta) - lhs;
    if(fabs(1.0 - xc/xc_old) < bisects) break;
  }
  if(n == 99999) {
    if(mpi_io_proc()) printf("no solution in xm");
    exit(0);
  }

  /* return */
  return xcp;
}

//******************************************************************************

/* fraction of up-scattered bremsstrahlung photon */
inline double fbrem(double y) {
  double out; 
  if(y <= 1e3) {
    out = 2.0*(y*y - y*(1.0+y)*exp(-1.0/y));
  } else {
    out = 1.0 - 2.0/3.0/y;
  }
  return out;
}

//******************************************************************************

/* compton y1 parameter */
inline double comptony1(double x, double tau, double theta) {
  double out = zeta*(tau*tau)*log(1.0+4.0*theta+16.0*theta*theta)/log(theta/x);
  return out;
}

//******************************************************************************

/* flat spectrum number density */
inline double flatn1(double x, double theta, double y) {
  double out = (2.0/M_PI)*(alphaf*alphaf*alphaf)*(x*x)*theta*(1.0/log(theta/x) + y/(1.0+y))/(RE*RE*RE);
  return out;
}

//******************************************************************************

/* wien spectrum number density */
inline double ngamma(double xm, double tau, double theta, double y1, double z, double nprot, double fb, double ndotbr, double r_size) {
  double gt;
  if(theta < 1.0) {
    gt = 1.0/(1.0 + 5.0*theta + 0.4*theta*theta);
  } else {
    gt = (0.1875)*(log(2.0*eta*theta)+0.75)/(1.0+(0.1/theta))/(theta*theta);
  }
  double ng = r_size/CL*(1.0 + zeta*gt*tau)*fb*ndotbr;
  return ng;
}

//*------------------------------------------------------------------------------------------------------------------------*//
//
// These sections are miscellaneous
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//******************************************************************************

/* positron fractions assuming local thermal equilibrium */
inline double get_zfrac(double nprot, double thetae) {
  double kt = thetae*ME*CL*CL;
  double lam_th = hplanck/sqrt(2*M_PI*ME*kt);
  double u = 4/nprot/nprot/pow(lam_th,6)*exp(-2/thetae);
  double zfrac = 0.5*(-1+sqrt(1+4*u));
  return zfrac;
}

//******************************************************************************

/* series expansion for error function */
inline double series_asym(double x_in) {
  double out = 1/x_in - 0.5/pow(x_in,3) + 0.75/pow(x_in,5) - 1.875/pow(x_in,7) + 6.5625/pow(x_in,9);
  return out;
}

//******************************************************************************

#endif
