#ifndef PARAMETER_HANDLER
#define PARAMETER_HANDLER

#include "param.h"
#include <Eigen/Dense>
#include <sstream>

using namespace std;

class ParameterHandler 
{
private:
   /* Param */
   parameter::parameter param;

   int M;                                   // Quadrature order, even, num_directions
   
   // Energy Group specifics
   int G;                                   // 1 group corresponds to grey case
   double efirst;                           // right edge energy for first group (keV)
   double elast;                            // right edge energy for last group (keV)
   double kappa_grey;                       // Grey opacity
   double rho;                              // Material density
   double kappa;                            // Absorption opacity
   double T;                                // Material temperature

   // Optional energy group bounds and absorption opacities
   bool have_group_bounds;                  // Contained in a txt file
   bool have_group_absorption_opacities;    // Contained in a txt file
   string filename_group_bounds;
   string filename_group_kappa;
   Eigen::VectorXd group_bounds;
   Eigen::VectorXd group_kappa;

   // Slab specifics
   double X;                           // Slab thickness
   int N;                              // Number of cells
   double dx;                          // cell size
   int bc_left_indicator;              // vacuum - 0, // TODO: Change this to an array that matches M / 2
                                       // source - 1, 
                                       // reflective - 2 
   int bc_right_indicator;             // vacuum - 0, // TODO: Change this to an array that matches M / 2
                                       // source - 1, 
                                       // reflective - 2 s
   bool use_mg_equilib;
   Eigen::MatrixXd psi_source; 

   double V;                          // Material velocity, beta = V / c

   /* Correction terms options */
   bool use_correction;

   /* Time Stepping Options */
   int ts_method = 3;                      // 1 - Strictly Backward Euler, 2 - Strictly CN, 3 - BDF2
   double dt = 0.00001;
   int max_timesteps = 1000;

   /* Validation Options */
   bool include_validation;

   /* Retrieve parameters from file */
   void get_parameters();

public:
   /* Default and non default constructors */
   ParameterHandler();
   ParameterHandler(const string filename);

   /* Getters */
   int get_M() { return M; }
   int get_G() { return G; }
   double get_efirst() { return efirst; }
   double get_elast() { return elast; }
   bool get_have_group_bounds() { return have_group_bounds; }
   double get_kappa_grey() { return kappa_grey; }
   double get_X() { return X; }
   int get_N() { return N; }
   double get_dx() { return dx; }
   int get_bc_left_indicator() { return bc_left_indicator; }
   int get_bc_right_indicator() { return bc_right_indicator; }
   bool get_use_mg_equilib() { return use_mg_equilib; }
   double get_rho() { return rho; }
   double get_kappa() { return kappa; }
   bool get_have_group_absorption_opacities() { return have_group_absorption_opacities; }
   double get_T() { return T; }
   double get_V() { return V; }
   bool get_use_correction() { return use_correction; }
   int get_ts_method() { return ts_method; }
   double get_dt() { return dt; }
   int get_max_timesteps() { return max_timesteps; }

   bool get_validation() { return include_validation; }
   
   void get_psi_source(Eigen::Ref<Eigen::MatrixXd> psi_source);
   void get_group_bounds(Eigen::Ref<Eigen::VectorXd> group_bounds);
   void get_group_kappa(Eigen::Ref<Eigen::VectorXd> group_kappa);
};

#endif