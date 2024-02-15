#ifndef COMPILE_TIME_VALS
#define COMPILE_TIME_VALS

#include "constants.h"

namespace ctv
{
   const static int M = 2;                              // Quadrature order, even, num_directions
   const double G_x[M] = {-0.5773502692, 0.5773502692}; // Quadrature points
   const double G_w[M] = {2*Constants::PI, 2*Constants::PI};              // Quadrature weights, sum to 4Pi

   // const static int M = 4;                              // Quadrature order, even
   // const double G_x[M] = {-0.861136, -0.339981, 0.339981, 0.861136}; // Quadrature points
   // const double G_w[M] = {2*M_PI*0.347855, 2*M_PI*0.652145, 2*M_PI*0.652145, 2*M_PI*0.347855};    // Quadrature weights, sum to 4Pi

   // Energy Group specifics
   const static int G = 1;                              // 1 group corresponds to grey case
   const static double efirst = 10.;                    // right edge energy for first group (keV)
   const static double elast = 10.;                     // right edge energy for last group (keV)
   const static double kappa_grey = 1.;                 // Grey opacity
   
   const static double X = 1.;                          // Slab thickness
   const static int N = 100;                            // Number of cells
   const static double dx = X / N;                      // cell size
   const static int bc_left_indicator = 2;              // vacuum - 0, // TODO: Change this to an array that matches M / 2
                                                        // source - 1, 
                                                        // reflective - 2 
   const static int bc_right_indicator = 1;             // vacuum - 0, // TODO: Change this to an array that matches M / 2
                                                        // source - 1, 
                                                        // reflective - 2 s

   const static double rho = 1.;                        // Material density
   const static double kappa = 1.;                      // Absorption opacity
   const static double T = 1.;                          // Material temperature
   const static double V = 0.;                          // Material velocity, beta = V / c

   /* Correction terms options */
   const bool use_correction = false;
   const double validation_tolerance = 1.E-6;

   /* Time Stepping Options */
   const static int ts_method = 3;                      // 1 - Strictly Backward Euler, 2 - Strictly CN, 3 - BDF2
   const static double dt = 0.00001;
   const static int max_timesteps = 1000;

   // If either boundary indicator is source, those values need to be specified
   const double test_bc = Constants::RADIATION_CONSTANT_A * Constants::SPEED_OF_LIGHT * pow(T, 4) / (Constants::FOUR_PI);
   const double psi_source[M] = {4.0931, 4.0931};
   // double psi_source[M] = {test_bc, test_bc, test_bc, test_bc};
}

#endif // COMPILE_TIME_VALS