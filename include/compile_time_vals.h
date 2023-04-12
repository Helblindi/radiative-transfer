#ifndef COMPILE_TIME_VALS
#define COMPILE_TIME_VALS

using namespace std;

namespace rt
{
   const static int M = 2;                            // Quadrature order, even
   const double G_x[M] = {-0.5773502692, 0.5773502692}; // Quadrature points
   const double G_w[M] = {6.283185307, 6.283185307};  // Quadrature weights, sum to 4Pi
   const static double X = 2.;                        // Slab thickness
   const static int N = 100;                          // Number of cells
   const static double dx = X / N;                    // cell size
   const static double a = .01372;                    // Radiation constant, jerks/ (cm^3/keV^4)
   const static double c = 299.8;                     // Speed of light, cm/shake
   const static int bc_left_indicator = 1;            // vacuum - 0, // TODO: Change this to an array that matches M / 2
                                                      // source - 1, 
                                                      // reflective - 2 
   const static int bc_right_indicator = 1;           // vacuum - 0, // TODO: Change this to an array that matches M / 2
                                                      // source - 1, 
                                                      // reflective - 2 

   // If either boundary indicator is source, those values need to be specified
   // const static double TL = 0.4;
   // const static double TR = 0.6;
   const double psi_source[M] = {0.3273225123, 0.3273225123};

   // TODO: Change to be cell wise values
   // TODO: Output plot of initial profiles
   const static double rho = 1.;                      // Material density
   const static double kappa = 1.;                    // Absorption opacity
   const static double T = 1.;                        // Material temperature
}

#endif // COMPILE_TIME_VALS