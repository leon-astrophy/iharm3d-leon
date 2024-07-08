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
#include <gsl/gsl_sf_bessel.h>

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

    // Set positron mass to its floor values //
    S->P[RPL][k][j][i] = ZMIN*ME_MP*S->P[RHO][k][j][i];

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

}

//******************************************************************************

// compute temperature per grid cells //
inline void find_temp_1zone(struct GridGeom *G, struct FluidState *Ss, int i, int j, int k)
{

  // number density, all in c.g.s unit //
  // nprot - proton, npost - prositron, nelec - electron //
  double nprot, npost, nelec, ntot;

  // calculate proton and positron //
  nprot = Ss->P[RHO][k][j][i]*RHO_unit/MP, npost = Ss->P[RPL][k][j][i]*RHO_unit/ME;

  // electron by charge neutrality //
  nelec = nprot + npost, ntot = nelec + nprot + npost;

  // temperature at i,j,k //
  temp[k][j][i] = Ss->P[UU][k][j][i]*U_unit*(gam - 1.)/(ntot*KBOL); 

}

//******************************************************************************

// compute pair production rate per grid cells //
inline void pair_production_1zone(struct GridGeom *G, struct FluidState *Ss, struct FluidState *Sf, int i, int j, int k , double dt_step)
{

  /***********************************************************************/
  // number density, all in c.g.s unit //

  // nprot - proton, npost - prositron, nelec - electron //
  double nprot, npost, nelec, ntot;

  /***********************************************************************/
  
  // temperature at i-1,j,k and i+1,j,k //
  double tm1_x = temp[k][j][i-1], tp1_x= temp[k][j][i+1];

  // temperature at i,j-1,k and i,j+1,k//
  double tm1_y= temp[k][j-1][i], tp1_y= temp[k][j+1][i];

  // temperature at i,j,k-1 and i,j,k+1 //
  double tm1_z= temp[k-1][j][i], tp1_z= temp[k+1][j][i];

  // temperature at i,,k //
  double t_c = temp[k][j][i];

  /***********************************************************************/
  
  // calculate proton and positron //
  nprot = Ss->P[RHO][k][j][i]*RHO_unit/MP, npost = Ss->P[RPL][k][j][i]*RHO_unit/ME;

  // electron by charge neutrality //
  nelec = nprot + npost, ntot = nelec + nprot + npost;

  /***********************************************************************/

  // dTdx //
  double dTdx = 0.5*(tp1_x - tm1_x)/dx[1];

  // dTdy //
  double dTdy = 0.5*(tp1_y - tm1_y)/dx[2];

  // dTdz //
  double dTdz = 0.5*(tp1_z - tm1_z)/dx[3];

  /***********************************************************************/

  // define a 4-vector //
  double gradT[NDIM], gradT_con[NDIM]; 

  // assign temperature gradient 
  gradT[0] = 0.0, gradT[1] = dTdx, gradT[2] = dTdy, gradT[3] = dTdz;

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

  // get the thermal scale height, remeber to convert to CGS //
  double h_th = 0.25*t_c/norm_gradT*L_unit;

  // optical depth //
  double tau_depth = 2.0*(nelec + npost)*sigma_t*h_th; 
  
  /***********************************************************************/

  // dimensionless temperature and postiron fraction //
  double thetae = KBOL*t_c/(ME*CL*CL); 
  double zfrac = npost/nprot;

  /***********************************************************************/
  /* now calculate pair production rate */

  //if(!isnan(h_th)) {
  
    // net pair production rate, note the rate is in the CGS unit!!! //
    double net_rate = ndot_net(zfrac, tau_depth, nprot, thetae, h_th);
    
    // positron mass production rate, need to convert to code unit!!! //
    double rhopdot = ME*net_rate*(T_unit/RHO_unit);

    // quality factor //
    double qfac = fabs(Ss->P[RPL][k][j][i]/rhopdot);

    /* if the source term is too steep, implement implicit solver */
    /* Crank-Nicolson method, inspired by Lia's thesis */
    /* Basically a root finding, so use bisection mtehod */
    double dummy;
    if(dt_step > q_alpha*qfac) {

      printf("mdot too steep\n");

      /* left state */
      double zl = zfrac;
      double tl = t_c/(zl + 1)*(zfrac + 1);
      double h_l = h_th*tl/t_c;
      double taul = 2.0*(2*zl + 1)*nprot*sigma_t*h_l;
      double theta_l = thetae*tl/t_c;
      dummy = ndot_net(zfrac, tau_depth, nprot, thetae, h_th) + ndot_net(zl, taul, nprot, theta_l, h_l);
      double fl = (zl - zfrac) - 0.5*dt_step*dummy/nprot;

      /* right state */
      double zr = 1.0;
      double tr = t_c/(zr + 1)*(zfrac + 1);
      double h_r = h_th*tr/t_c;
      double taur = 2.0*(2*zr + 1)*nprot*sigma_t*h_r;
      double theta_r = thetae*tr/t_c;
      dummy = ndot_net(zfrac, tau_depth, nprot, thetae, h_th) + ndot_net(zr, taur, nprot, theta_r, h_r);
      double fr = (zr - zfrac) - 0.5*dt_step*dummy/nprot;

      /* if the right state guess is too poor, scale it up by 100 */
      if(fr*fl > 0) {
        zr = 100;
        tr = t_c/(zr + 1)*(zfrac + 1);
        h_r = h_th*tr/t_c;
        taur = 2.0*(2*zr + 1)*nprot*sigma_t*h_r;
        theta_r = thetae*tr/t_c;
        dummy = ndot_net(zfrac, tau_depth, nprot, thetae, h_th) + ndot_net(zr, taur, nprot, theta_r, h_r);
        fr = (zr - zfrac) - 0.5*dt_step*dummy/nprot;
        if(fr*fl > 0) {
          printf("invalid guess in zfrac\n");
          exit(0);
        }
      }

      /* define the center state */
      double zc, tc, h_c, tauc,theta_c, fc, zc_old;

      /* bisection method counting */
      int count;

      /* now iterate until converges */
      for (count = 0; count < 99999; count++) {

        /* backup */
        if(count > 0) {
          zc_old = zc;
        }

        /* center state */
        zc = 0.5*(zl + zr);
        tc = t_c/(zc + 1)*(zfrac + 1);
        h_c = h_th*tc/t_c;
        tauc = 2.0*(2*zc + 1)*nprot*sigma_t*h_c;
        theta_c = thetae*tc/t_c;
        dummy = ndot_net(zfrac, tau_depth, nprot, thetae, h_th) + ndot_net(zc, tauc, nprot, theta_c, h_c);
        fc = (zc - zfrac) - 0.5*dt_step*dummy/nprot;

        /* check the sign */
        if (fl*fc > 0) {
          zl = zc;
        } else if (fr*fc > 0) {
          zr = zc;
        }

        /* determine if need to exit */
        if(count > 0) {
          if(fabs(1.0 - zc_old/zc) < bisects) {
            break;
          }
        }
      }
      if(count == 99999) {
        printf("No solution\n");
	      exit(0);
      }

      /* assign new positron mass, remember to convert back to code unit !!! */
      Sf->P[RPL][k][j][i] = zc*nprot*ME/RHO_unit;
      
    /* otherwise, march forward by time */
    } else {

      Sf->P[RPL][k][j][i] += rhopdot*dt_step;

    }

  //}

}

//******************************************************************************

// functions for computing function g(theta) in pair annhilation rates //
inline double nadot(double z, double nprot, double theta) {
  double out = (3.0/8.0)*sigma_t*CL*(nprot*nprot*z*(z+1.0))/(1.0+2.0*theta*theta/log(1.12*theta+1.3));
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

/* pair production rate due to photon-photon or photon-particle collision */
inline double ncdot(double ngamma, double theta, double nprot, double z, double n1, double fb, double ndotbr) {
  //double ndotww = get_ndotww(ngamma, theta);
  //double ndotwp = get_ndotwp(ngamma, nprot, theta);
  //double ndotwe = get_ndotwe(ngamma, nprot, z, theta);
  //double ndotwf = get_ndotwf(n1, ngamma, theta);
  double ndotee = get_ndotee(nprot, z, theta);
  //double out = ndotee + ndotww + ndotwp + ndotwe + ndotwf ;
  double out = ndotee;
  return out;
}

//******************************************************************************

/* net pair production rate */
inline double ndot_net(double zfrac, double taut, double nprot, double theta, double r_size) {
  //double xm = find_xm(zfrac, taut, nprot, theta);
  //double ndotbr = get_ndotbr(zfrac, theta, xm, nprot);
  //double y1 = comptony1(xm, taut, theta);
  //double fb = fbrem(y1);
  //double n1 = flatn1(xm, theta, y1);
  //double ng = ngamma(xm, taut, theta, y1, zfrac, nprot, fb, ndotbr, r_size);
  //double nc = ncdot(ng, theta, nprot, zfrac, n1, fb, ndotbr);
  double nc = ncdot(1.0, theta, nprot, zfrac, 1.0, 1.0, 1.0); 
  double na = nadot(zfrac, nprot, theta);
  return nc - na;
}

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
    printf("poor initial guess xm");
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
    printf("no solution in xm");
    exit(0);
  }

  /* return */
  return xcp;
}

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
    gt = (3.0/16.0)*(log(2.0*eta*theta)+0.75)/(1.0+(0.1/theta))/(theta*theta);
  }
  double ng = r_size/CL*(1.0 + zeta*gt*tau)*fb*ndotbr;
  return ng;
}

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

//******************************************************************************

#endif
