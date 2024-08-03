//******************************************************************************
//*                                                                            *
//* RECONSTRUCTION.C                                                           *
//*                                                                            *
//* RECONSTRUCTION ALGORITHMS                                                  *
//*                                                                            *
//******************************************************************************

//include header files
#include "decs.h"
#include "math.h"

//*********************************************************************************************************************

// Choose reconstruction algorithm according to user input 
// NOTE: For linear, MC is the only limiter used
#if RECONSTRUCTION == LINEAR
#define RECON_ALGO linear_mc
#elif RECONSTRUCTION == WENO
#define RECON_ALGO weno
#elif RECONSTRUCTION == MP5
#define RECON_ALGO mp5
#elif RECONSTRUCTION == PPM
#define RECON_ALGO ppm
#elif RECONSTRUCTION == PPMX
#define RECON_ALGO ppmx
#elif RECONSTRUCTION == WENOZ
#define RECON_ALGO wenoz
#else
#error "Reconstruction not specified!"
#endif

// Sanity checks
#if (RECONSTRUCTION == WENO || RECONSTRUCTION == MP5 || RECONSTRUCTION == PPM || RECONSTRUCTION == PPMX || RECONSTRUCTION == WENOZ) && NG < 3
#error "not enough ghost zones! PPM/WENO/MP5 + NG < 3\n"
#endif

//*********************************************************************************************************************

// define reconstruction functions 
void linear_mc(double unused1, double x1, double x2, double x3, double unused2, double *lout, double *rout);
void weno(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout);
double median(double a, double b, double c);
double mp5_subcalc(double Fjm2, double Fjm1, double Fj, double Fjp1, double Fjp2);
void mp5(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout);
void ppm(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout);
void ppmx(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout);
double mc(const double dm, const double dp, const double alpha);
void weno_z(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout);

//*********************************************************************************************************************

inline void linear_mc(double unused1, double x1, double x2, double x3, double unused2, double *lout, double *rout)
{
  double Dqm,Dqp,Dqc,s;

  Dqm = 2. * (x2 - x1);
  Dqp = 2. * (x3 - x2);
  Dqc = 0.5 * (x3 - x1);

  s = Dqm * Dqp;

  if (s <= 0.)
    s = 0.;
  else {
    if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
      s = Dqm;
    else if (fabs(Dqp) < fabs(Dqc))
      s = Dqp;
    else
      s = Dqc;
  }

  // Reconstruct left, right
  *lout = x2 - 0.5*s;
  *rout = x2 + 0.5*s;
}

//*********************************************************************************************************************

// WENO interpolation. See Tchekhovskoy et al. 2007 (T07), Shu 2011 (S11)
// Implemented by Monika Moscibrodzka
inline void weno(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout)
{
  // S11 1, 2, 3
  double vr[3], vl[3];
  vr[0] =  (3./8.)*x1 - (5./4.)*x2 + (15./8.)*x3;
  vr[1] = (-1./8.)*x2 + (3./4.)*x3 + (3./8.)*x4;
  vr[2] =  (3./8.)*x3 + (3./4.)*x4 - (1./8.)*x5;

  vl[0] =  (3./8.)*x5 - (5./4.)*x4 + (15./8.)*x3;
  vl[1] = (-1./8.)*x4 + (3./4.)*x3 + (3./8.)*x2;
  vl[2] =  (3./8.)*x3 + (3./4.)*x2 - (1./8.)*x1;

  // Smoothness indicators, T07 A18 or S11 8
  double beta[3];
  beta[0] = (13./12.)*pow(x1 - 2.*x2 + x3, 2) +
            (1./4.)*pow(x1 - 4.*x2 + 3.*x3, 2);
  beta[1] = (13./12.)*pow(x2 - 2.*x3 + x4, 2) +
            (1./4.)*pow(x4 - x2, 2);
  beta[2] = (13./12.)*pow(x3 - 2.*x4 + x5, 2) +
            (1./4.)*pow(x5 - 4.*x4 + 3.*x3, 2);

  // Nonlinear weights S11 9
  double den, wtr[3], Wr, wr[3], wtl[3], Wl, wl[3], eps;

  // Leon's comment: my experience told me that we need a rather small 
  // epsilon to make the scheme stable, I suggest 1.e-40, but I am unsure
  // Maybe need to consult some numerics expert
  eps=1.e-26;
  //eps=1.e-40; 
  
  den = eps + beta[0]; den *= den; wtr[0] = (1./16.)/den;
  den = eps + beta[1]; den *= den; wtr[1] = (5./8. )/den;
  den = eps + beta[2]; den *= den; wtr[2] = (5./16.)/den;
  Wr = wtr[0] + wtr[1] + wtr[2];
  wr[0] = wtr[0]/Wr ;
  wr[1] = wtr[1]/Wr ;
  wr[2] = wtr[2]/Wr ;

  den = eps + beta[2]; den *= den; wtl[0] = (1./16.)/den;
  den = eps + beta[1]; den *= den; wtl[1] = (5./8. )/den;
  den = eps + beta[0]; den *= den; wtl[2] = (5./16.)/den;
  Wl = wtl[0] + wtl[1] + wtl[2];
  wl[0] = wtl[0]/Wl;
  wl[1] = wtl[1]/Wl;
  wl[2] = wtl[2]/Wl;

  *lout = vl[0]*wl[0] + vl[1]*wl[1] + vl[2]*wl[2];
  *rout = vr[0]*wr[0] + vr[1]*wr[1] + vr[2]*wr[2];
}

//*********************************************************************************************************************

// MP5 reconstruction from PLUTO
// Imported by Mani Chandra
#define MINMOD(a, b) ((a)*(b) > 0.0 ? (fabs(a) < fabs(b) ? (a):(b)):0.0)
inline double median(double a, double b, double c)
{
  return (a + MINMOD(b - a, c - a));
}
#define ALPHA (4.0)
#define EPSM (1.e-12)
inline double mp5_subcalc(double Fjm2, double Fjm1, double Fj, double Fjp1, double Fjp2)
{
  double f, d2, d2p, d2m;
  double dMMm, dMMp;
  double scrh1,scrh2, Fmin, Fmax;
  double fAV, fMD, fLC, fUL, fMP;

  f  = 2.0*Fjm2 - 13.0*Fjm1 + 47.0*Fj + 27.0*Fjp1 - 3.0*Fjp2;
  f /= 60.0;

  fMP = Fj + MINMOD(Fjp1 - Fj, ALPHA*(Fj - Fjm1));

  if ((f - Fj)*(f - fMP) <= EPSM)
    return f;

  d2m = Fjm2 + Fj   - 2.0*Fjm1;              // Eqn. 2.19
  d2  = Fjm1 + Fjp1 - 2.0*Fj;
  d2p = Fj   + Fjp2 - 2.0*Fjp1;              // Eqn. 2.19

  scrh1 = MINMOD(4.0*d2 - d2p, 4.0*d2p - d2);
  scrh2 = MINMOD(d2, d2p);
  dMMp  = MINMOD(scrh1,scrh2);               // Eqn. 2.27

  scrh1 = MINMOD(4.0*d2m - d2, 4.0*d2 - d2m);
  scrh2 = MINMOD(d2, d2m);
  dMMm  = MINMOD(scrh1,scrh2);               // Eqn. 2.27

  fUL = Fj + ALPHA*(Fj - Fjm1);              // Eqn. 2.8
  fAV = 0.5*(Fj + Fjp1);                     // Eqn. 2.16
  fMD = fAV - 0.5*dMMp;                      // Eqn. 2.28
  fLC = 0.5*(3.0*Fj - Fjm1) + 4.0/3.0*dMMm;  // Eqn. 2.29

  scrh1 = fmin(Fj, Fjp1); scrh1 = fmin(scrh1, fMD);
  scrh2 = fmin(Fj, fUL);    scrh2 = fmin(scrh2, fLC);
  Fmin  = fmax(scrh1, scrh2);                // Eqn. (2.24a)

  scrh1 = fmax(Fj, Fjp1); scrh1 = fmax(scrh1, fMD);
  scrh2 = fmax(Fj, fUL);    scrh2 = fmax(scrh2, fLC);
  Fmax  = fmin(scrh1, scrh2);                // Eqn. 2.24b

  f = median(f, Fmin, Fmax);                 // Eqn. 2.26
  return f;
}

inline void mp5(double x1, double x2, double x3, double x4, double x5, double *lout,
  double *rout)
{
  *rout = mp5_subcalc(x1, x2, x3, x4, x5);
  *lout = mp5_subcalc(x5, x4, x3, x2, x1);
}
#undef MINMOD

//---------------------------------------------------------------------------------------------------------------------
//
// Leon's "implementation" of reconstruction method
//
//---------------------------------------------------------------------------------------------------------------------

//*********************************************************************************************************************

// PPM stolen from Kharma //
inline void ppm(double q_im2, double q_im1, double q_i, double q_ip1, double q_ip2, double *lout, double *rout)
{

  /* local storing left and right state */ 
  double qlv, qrv;

  //---- Interpolate L/R values (CS eqn 16, PH 3.26 and 3.27) ----
  // qlv = q at left  side of cell-center = q[i-1/2] = a_{j,-} in CS
  // qrv = q at right side of cell-center = q[i+1/2] = a_{j,+} in CS
  qlv = (7.*(q_i + q_im1) - (q_im2 + q_ip1))/12.0;
  qrv = (7.*(q_i + q_ip1) - (q_im1 + q_ip2))/12.0;

  //---- limit qrv and qlv to neighboring cell-centered values (CS eqn 13) ----
  qlv = fmax(qlv, fmin(q_i, q_im1));
  qlv = fmin(qlv, fmax(q_i, q_im1));
  qrv = fmax(qrv, fmin(q_i, q_ip1));
  qrv = fmin(qrv, fmax(q_i, q_ip1));

  //--- monotonize interpolated L/R states (CS eqns 14, 15) ---
  double qc = qrv - q_i;
  double qd = qlv - q_i;
  if ((qc*qd) >= 0.0) {
    qlv = q_i;
    qrv = q_i;
  } else {
    if (fabs(qc) >= 2.0*fabs(qd)) {
      qrv = q_i - 2.0*qd;
    }
    if (fabs(qd) >= 2.0*fabs(qc)) {
      qlv = q_i - 2.0*qc;
    }
  }

  // finally, assign state //
  *rout = qrv;
  *lout = qlv;
}

//*********************************************************************************************************************

// PPMX stolen from Kharma //
inline void ppmx(double q_im2, double q_im1, double q_i, double q_ip1, double q_ip2, double *lout, double *rout)
{

  /* local storing left and right state */ 
  double qlv, qrv;

  //---- Compute L/R values (CS eqns 12-15, PH 3.26 and 3.27) ----
  // qlv = q at left  side of cell-center = q[i-1/2] = a_{j,-} in CS
  // qrv = q at right side of cell-center = q[i+1/2] = a_{j,+} in CS
  qlv = (7.*(q_i + q_im1) - (q_im2 + q_ip1))/12.0;
  qrv = (7.*(q_i + q_ip1) - (q_im1 + q_ip2))/12.0;

  //---- Apply CS monotonicity limiters to qrv and qlv ----
  // approximate second derivatives at i-1/2 (PH 3.35)
  // KGF: add the off-center quantities first to preserve FP symmetry
  double d2qc = 3.0*((q_im1 + q_i) - 2.0*qlv);
  double d2ql = (q_im2 + q_i  ) - 2.0*q_im1;
  double d2qr = (q_im1 + q_ip1) - 2.0*q_i;

  // limit second derivative (PH 3.36)
  double d2qlim = 0.0;
  double lim_slope = fmin(fabs(d2ql),fabs(d2qr));
  if (d2qc > 0.0 && d2ql > 0.0 && d2qr > 0.0) {
    d2qlim = copysign(1.0,d2qc)*fmin(1.25*lim_slope,fabs(d2qc));
  }
  if (d2qc < 0.0 && d2ql < 0.0 && d2qr < 0.0) {
    d2qlim = copysign(1.0,d2qc)*fmin(1.25*lim_slope,fabs(d2qc));
  }
  
  // compute limited value for qlv (PH 3.33 and 3.34)
  if (((q_im1 - qlv)*(q_i - qlv)) > 0.0) {
    qlv = 0.5*(q_i + q_im1) - d2qlim/6.0;
  }

  // approximate second derivatives at i+1/2 (PH 3.35)
  // KGF: add the off-center quantities first to preserve FP symmetry
  d2qc = 3.0*((q_i + q_ip1) - 2.0*qrv);
  d2ql = d2qr;
  d2qr = (q_i + q_ip2) - 2.0*q_ip1;

  // limit second derivative (PH 3.36)
  d2qlim = 0.0;
  lim_slope = fmin(fabs(d2ql),fabs(d2qr));
  if (d2qc > 0.0 && d2ql > 0.0 && d2qr > 0.0) {
    d2qlim = copysign(1.0,d2qc)*fmin(1.25*lim_slope,fabs(d2qc));
  } 
  if (d2qc < 0.0 && d2ql < 0.0 && d2qr < 0.0) {
    d2qlim = copysign(1.0,d2qc)*fmin(1.25*lim_slope,fabs(d2qc));
  }
  // compute limited value for qrv (PH 3.33 and 3.34)
  if (((q_i - qrv)*(q_ip1 - qrv)) > 0.0) {
    qrv = 0.5*(q_i + q_ip1) - d2qlim/6.0;
  }

  //---- identify extrema, use smooth extremum limiter ----
  // CS 20 (missing "OR"), and PH 3.31
  double qa = (qrv - q_i)*(q_i - qlv);
  double qb = (q_im1 - q_i)*(q_i - q_ip1);
  if (qa <= 0.0 || qb <= 0.0) {

    // approximate secnd derivates (PH 3.37)
    // KGF: add the off-center quantities first to preserve FP symmetry
    double d2q  = 6.0*(qlv + qrv - 2.0*q_i);
    double d2qc = (q_im1 + q_ip1) - 2.0*q_i;
    double d2ql = (q_im2 + q_i  ) - 2.0*q_im1;
    double d2qr = (q_i   + q_ip2) - 2.0*q_ip1;

    // limit second derivatives (PH 3.38)
    d2qlim = 0.0;
    lim_slope = fmin(fabs(d2ql),fabs(d2qr));
    lim_slope = fmin(fabs(d2qc),lim_slope);
    if (d2qc > 0.0 && d2ql > 0.0 && d2qr > 0.0 && d2q > 0.0) {
      d2qlim = copysign(1.0,d2q)*fmin(1.25*lim_slope,fabs(d2q));
    }
    if (d2qc < 0.0 && d2ql < 0.0 && d2qr < 0.0 && d2q < 0.0) {
      d2qlim = copysign(1.0,d2q)*fmin(1.25*lim_slope,fabs(d2q));
    }

    // limit L/R states at extrema (PH 3.39)
    double rho = 0.0;
    if ( fabs(d2q) > (1.0e-12)*fmax( fabs(q_im1), fmax(fabs(q_i),fabs(q_ip1))) ) {
      // Limiter is not sensitive to round-off error.  Use limited slope
      rho = d2qlim/d2q;
    }
    qlv = q_i + (qlv - q_i)*rho;
    qrv = q_i + (qrv - q_i)*rho;

  } else {

    // Monotonize again, away from extrema (CW eqn 1.10, PH 3.32)
    double qc = qrv - q_i;
    double qd = qlv - q_i;
    if (fabs(qc) >= 2.0*fabs(qd)) {
      qrv = q_i - 2.0*qd;
    }
    if (fabs(qd) >= 2.0*fabs(qc)) {
      qlv = q_i - 2.0*qc;
    }

  }

  // finally, assign state //
  *rout = qrv;
  *lout = qlv;

}

//*********************************************************************************************************************

/* limiter */
inline double mc(const double dm, const double dp, const double alpha) {
  const double dc = (dm * dp > 0.0) * 0.5 * (dm + dp);
  return copysign(fmin(fabs(dc), alpha * fmin(fabs(dm), fabs(dp))), dc);
}

/* weno-z stolen from KHARMA */
inline void weno_z(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout)
{
  double w5alpha[3][3] = {{1.0 / 3.0, -7.0 / 6.0, 11.0 / 6.0},
                                  {-1.0 / 6.0, 5.0 / 6.0, 1.0 / 3.0},
                                  {1.0 / 3.0, 5.0 / 6.0, -1.0 / 6.0}};
  double w5gamma[3] = {0.1, 0.6, 0.3};
  double eps = 1e-100;
  double thirteen_thirds = 13.0 / 3.0;

  double a = x1 - 2 * x2 + x3;
  double b = x1 - 4.0 * x2 + 3.0 * x3;
  double beta0 = thirteen_thirds * a * a + b * b + eps;
  a = x2 - 2.0 * x3 + x4;
  b = x4 - x2;
  double beta1 = thirteen_thirds * a * a + b * b + eps;
  a = x3 - 2.0 * x4 + x5;
  b = x5 - 4.0 * x4 + 3.0 * x3;
  double beta2 = thirteen_thirds * a * a + b * b + eps;
  const double tau5 = fabs(beta2 - beta0);

  beta0 = (beta0 + tau5) / beta0;
  beta1 = (beta1 + tau5) / beta1;
  beta2 = (beta2 + tau5) / beta2;

  double w0 = w5gamma[0] * beta0 + eps;
  double w1 = w5gamma[1] * beta1 + eps;
  double w2 = w5gamma[2] * beta2 + eps;
  double wsum = 1.0 / (w0 + w1 + w2);
  *rout = w0 * (w5alpha[0][0] * x1 + w5alpha[0][1] * x2 + w5alpha[0][2] * x3);
  *rout += w1 * (w5alpha[1][0] * x2 + w5alpha[1][1] * x3 + w5alpha[1][2] * x4);
  *rout += w2 * (w5alpha[2][0] * x3 + w5alpha[2][1] * x4 + w5alpha[2][2] * x5);
  *rout *= wsum;
  const double alpha_r =
      3.0 * wsum * w0 * w1 * w2 /
          (w5gamma[2] * w0 * w1 + w5gamma[1] * w0 * w2 + w5gamma[0] * w1 * w2) +
      eps;

  w0 = w5gamma[0] * beta2 + eps;
  w1 = w5gamma[1] * beta1 + eps;
  w2 = w5gamma[2] * beta0 + eps;
  wsum = 1.0 / (w0 + w1 + w2);
  *lout = w0 * (w5alpha[0][0] * x5 + w5alpha[0][1] * x4 + w5alpha[0][2] * x3);
  *lout += w1 * (w5alpha[1][0] * x4 + w5alpha[1][1] * x3 + w5alpha[1][2] * x2);
  *lout += w2 * (w5alpha[2][0] * x3 + w5alpha[2][1] * x2 + w5alpha[2][2] * x1);
  *lout *= wsum;
  const double alpha_l =
      3.0 * wsum * w0 * w1 * w2 /
          (w5gamma[2] * w0 * w1 + w5gamma[1] * w0 * w2 + w5gamma[0] * w1 * w2) +
      eps;

  double dq = x4 - x3;
  dq = mc(x3 - x2, dq, 2.0);

  const double alpha_lin = 2.0 * alpha_r * alpha_l / (alpha_r + alpha_l);
  *rout = alpha_lin * *rout + (1.0 - alpha_lin) * (x3 + 0.5 * dq);
  *lout = alpha_lin * *lout + (1.0 - alpha_lin) * (x3 - 0.5 * dq);
}

//*********************************************************************************************************************

// Reconstruct according to dimensional sweep
void reconstruct(struct FluidState *S, GridPrim Pl, GridPrim Pr, int dir)
{
  timer_start(TIMER_RECON);
  if (dir == 1) {
#pragma omp parallel for collapse(3)
    PLOOP {
      KSLOOP(-1, N3) {
        JSLOOP(-1, N2) {
          ISLOOP(-1, N1) {
            RECON_ALGO(S->P[ip][k][j][i-2], S->P[ip][k][j][i-1], S->P[ip][k][j][i],
                 S->P[ip][k][j][i+1], S->P[ip][k][j][i+2], &(Pl[ip][k][j][i]),
                 &(Pr[ip][k][j][i]));
          }
        }
      }
    }
  } else if (dir == 2) {
#pragma omp parallel for collapse(3)
    PLOOP {
      KSLOOP(-1, N3) {
        JSLOOP(-1, N2) {
          ISLOOP(-1, N1) {
            RECON_ALGO(S->P[ip][k][j-2][i], S->P[ip][k][j-1][i], S->P[ip][k][j][i],
                 S->P[ip][k][j+1][i], S->P[ip][k][j+2][i], &(Pl[ip][k][j][i]),
                 &(Pr[ip][k][j][i]));
          }
        }
      }
    }
  } else if (dir == 3) {
#pragma omp parallel for collapse(3)
    PLOOP {
      KSLOOP(-1, N3) {
        JSLOOP(-1, N2) {
          ISLOOP(-1, N1) {
            RECON_ALGO(S->P[ip][k-2][j][i], S->P[ip][k-1][j][i], S->P[ip][k][j][i],
                 S->P[ip][k+1][j][i], S->P[ip][k+2][j][i], &(Pl[ip][k][j][i]),
                 &(Pr[ip][k][j][i]));
          }
        }
      }
    }
  }
  timer_stop(TIMER_RECON);
}

//*********************************************************************************************************************