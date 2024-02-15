#ifndef CORRECTION
#define CORRECTION


#include <compile_time_vals.h>
#include <Constants.h>
#include <Planck.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <cassert>
#include <cmath>
#include <iostream>

namespace rt
{

template<int num_groups>
class Correction 
{
private:
   // Class variables
   const double ac = Constants::RADIATION_CONSTANT_A*Constants::SPEED_OF_LIGHT,
                kcon = Constants::BOLTZMANN_CONSTANT_JPK,
                fourpi = Constants::FOUR_PI,    
                hcon = Constants::PLANCK_CONSTANT,
                ccon = Constants::SPEED_OF_LIGHT;
   
   const double T=ctv::T, kappa_grey=ctv::kappa_grey;

   Planck planck;
   Eigen::VectorXd B, dBdT;               // Planck integrals
   Eigen::VectorXd ukappa;                // Unnormalized and final opacities respectively
   Eigen::VectorXd ckappa;                // Opacities evaluated at group center energies
   Eigen::VectorXd emis_spec;

   /* Unchanging energy group information, stored as references */
   Eigen::Ref<Eigen::VectorXd> e_edge_ref, e_ave_ref, de_ave_ref; 
   Eigen::Ref<Eigen::MatrixXd> energy_discretization_ref; 

   /****** 
    * Vector quantities needed from caller for correction terms 
    * ******/
   // Values set in default constructor
   Eigen::Ref<Eigen::VectorXd> rho_ref,   // Reference to rho_vec
                               kappa_ref, // Reference to kappa_vec
                               T_ref;     // Reference to temperature
   // Values set in compute_correction call
   Eigen::VectorXd psi;                   // Intensities, changes per group

   /****
    * Class variables set for correction terms calculation
    * *****/
   Eigen::VectorXd dEB, dsigEdE, dkapEB;

   Eigen::VectorXd kappa_edge;                     // Edge opacities

   Eigen::MatrixXd cor1, cor2, cor3;               // correction terms 1, 2, and 3
   Eigen::Tensor<double, 3> total_correction;      // dir, group, cell

   // Helper functions
   double pf(double E, double T);
   void generate_planck_integrals();
   void generate_multigroup_opacities();
   void compute_group_edge_opacities();
   void compute_components_of_correction_source();
   void compute_correction_terms();

   void validate_planck_integrals();
   void validate_emission();

public:
   Correction(Eigen::Ref<Eigen::VectorXd> rho_vec,
              Eigen::Ref<Eigen::VectorXd> kappa_vec, 
              Eigen::Ref<Eigen::VectorXd> T_vec,
              Eigen::Ref<Eigen::VectorXd> e_edge,
              Eigen::Ref<Eigen::VectorXd> e_ave,
              Eigen::Ref<Eigen::VectorXd> de_ave,
              Eigen::Ref<Eigen::MatrixXd> energy_discretization);
   ~Correction() {}
   // TODO: Compute correction will need the updated rho, T
   bool validate_correction();
   void compute_correction(Eigen::Tensor<double, 3>& intensities);
   void get_correction(Eigen::Tensor<double, 3>& total_correction)
   {
      total_correction = this->total_correction;
   }
};

} // End namespace rt

#endif
