#ifndef CORRECTION
#define CORRECTION


#include "constants.h"
#include "ParameterHandler.h"
#include "Planck.h"

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;

namespace rt
{

class Correction 
{
private:
   ParameterHandler & ph;
   // Class variables
   const double ac = Constants::RADIATION_CONSTANT_A*Constants::SPEED_OF_LIGHT,
                kcon = Constants::BOLTZMANN_CONSTANT_JPK,
                fourpi = Constants::FOUR_PI,    
                hcon = Constants::PLANCK_CONSTANT,
                ccon = Constants::SPEED_OF_LIGHT;
   
   const double T=ph.get_T(), kappa_grey=ph.get_kappa_grey();

   const int num_groups = ph.get_G();

   Planck planck;
   Eigen::VectorXd B, dBdT;               // Planck integrals
   Eigen::VectorXd ukappa;                // Unnormalized and final opacities respectively
   Eigen::VectorXd ckappa;                // Opacities evaluated at group center energies
   Eigen::VectorXd emis_spec;

   /* Unchanging energy group information, stored as references */
   Eigen::Ref<Eigen::VectorXd> e_edge_ref, e_ave_ref, de_ave_ref; 
   Eigen::Ref<Eigen::MatrixXd> energy_discretization_ref; 

   /* Direction discretization, also references */
   Eigen::Ref<Eigen::VectorXd> m_mu_ref, m_wt_ref;

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

   bool validate_planck_integrals();
   bool validate_emission();

public:
   Correction(ParameterHandler & parameter_handler,
              Eigen::Ref<Eigen::VectorXd> rho_vec,
              Eigen::Ref<Eigen::VectorXd> kappa_vec, 
              Eigen::Ref<Eigen::VectorXd> T_vec,
              Eigen::Ref<Eigen::VectorXd> e_edge,
              Eigen::Ref<Eigen::VectorXd> e_ave,
              Eigen::Ref<Eigen::VectorXd> de_ave,
              Eigen::Ref<Eigen::MatrixXd> energy_discretization,
              Eigen::Ref<Eigen::VectorXd> m_mu,
              Eigen::Ref<Eigen::VectorXd> m_wt);
   ~Correction() {}
   // TODO: Compute correction will need the updated rho, T
   bool validate_correction();
   void compute_correction(Eigen::Tensor<double, 3>& intensities);
   void get_correction(Eigen::Tensor<double, 3>& total_correction)
   {
      total_correction = this->total_correction;
   }
   void get_B(Eigen::Ref<Eigen::VectorXd> B)
   {
      B = this->B;
   }
};

} // End namespace rt

#endif
