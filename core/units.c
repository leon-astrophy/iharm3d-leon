/******************************************************************************
 *                                                                            *
 * UNITS.C                                                                    *
 *                                                                            *
 * Set unit conversions between code and cgss                                 *
 * Written by H.S. Leon Chan at 2024                                          *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

void set_units()
{
  Mbh = mbh*MSUN;
  L_unit = GNEWT*Mbh/(CL*CL);
  T_unit = L_unit/CL;
  RHO_unit = M_unit*pow(L_unit,-3.);
  U_unit = RHO_unit*CL*CL;
}