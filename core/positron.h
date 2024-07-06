//******************************************************************************
//*                                                                            *
//* positron.h                                                                 *
//*                                                                            *
//* Header files for the positron module                                       *
//*                                                                            *
//******************************************************************************

//import headers
#pragma once
#include "decs.h"

// thomson scattering cross section in CGS
#define sigma_t (6.6524587051e-25)

// Leon's patch, fixup parameters for positron //
#define RPLMINLIMIT (1.e-30)
#define RPLMIN  (1.e-6)

// Leon's patch, these are all physical untis //
#define ME (9.1093826e-28  ) // Electron mass
#define MP (1.67262171e-24 ) // Proton mass
#define KBOL (1.3806505e-16  ) // Boltzmann constant
#define GNEWT (6.6742e-8      ) // Gravitational constant
#define CL (2.99792458e10  ) // Speed of light
#define R_E (2.8179403262e-13  ) // classical electron radius
#define A_F  (7.2973525693e-3  ) // fine structure constants 
#define MSUN (1.989e33        ) // Solar mass

// Leon's patch, unit conversion //
double Mbh, L_unit, T_unit, RHO_unit, U_unit, M_unit, mbh; 

//******************************************************************************