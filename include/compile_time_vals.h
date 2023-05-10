#ifndef COMPILE_TIME_VALS
#define COMPILE_TIME_VALS

using namespace std;

namespace rt
{
   const static int M = 2;                              // Quadrature order, even
   const double G_x[M] = {-0.5773502692, 0.5773502692}; // Quadrature points
   const double G_w[M] = {2*M_PI, 2*M_PI};              // Quadrature weights, sum to 4Pi

   // const static int M = 4;                              // Quadrature order, even
   // const double G_x[M] = {-0.861136, -0.339981, 0.339981, 0.861136}; // Quadrature points
   // const double G_w[M] = {2*M_PI*0.347855, 2*M_PI*0.652145, 2*M_PI*0.652145, 2*M_PI*0.347855};    // Quadrature weights, sum to 4Pi
   
   const static double X = 2.;                          // Slab thickness
   const static int N = 100;                            // Number of cells
   const static double dx = X / N;                      // cell size
   const static double a = .01372;                      // Radiation constant, jerks/ (cm^3/keV^4)
   const static double c = 299.8;                       // Speed of light, cm/shake
   const static int bc_left_indicator = 2;              // vacuum - 0, // TODO: Change this to an array that matches M / 2
                                                        // source - 1, 
                                                        // reflective - 2 
   const static int bc_right_indicator = 1;             // vacuum - 0, // TODO: Change this to an array that matches M / 2
                                                        // source - 1, 
                                                        // reflective - 2 s

   const static double rho = 1.;                        // Material density
   const static double kappa = 1.;                      // Absorption opacity
   const static double T = 1.;                          // Material temperature

   /* Time Stepping Options */
   const static int ts_method = 2;                      // 1 - Backward Euler, 2 - BDF2, 3 - CN, 4 - BDF2 from Morel-Lou
   const static double dt = 0.0001;
   const static int _max_timesteps = 500;

   // If either boundary indicator is source, those values need to be specified
   const double test_bc = a * c * pow(T, 4) / (4*M_PI);
   double psi_source[M] = {test_bc, test_bc};
   // double psi_source[M] = {test_bc, test_bc, test_bc, test_bc};

}

#endif // COMPILE_TIME_VALS