#ifndef SOLVER 
#define SOLVER

#include "constants.h"
#include "ParameterHandler.h"
#include "correction.h"
#include "GLQuad.h"

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <iostream>
#include <iomanip>

using namespace std;

namespace rt
{
class Solver
{
private:
   Eigen::MatrixXd psi_source;
   ParameterHandler & ph;
   int M, N, num_groups;
   double dx, dt;

   // Direction quadrature, to also be used in Correction
   Eigen::VectorXd m_mu, m_wt;

   const double ac = Constants::RADIATION_CONSTANT_A*Constants::SPEED_OF_LIGHT;
   Eigen::Tensor<double, 3>& psi_mat_ref;         // solution, angular intensity ref
   Eigen::Ref<Eigen::MatrixXd> phi_ref;           // angle-integrated intensity ref
   Eigen::MatrixXd phi_plus;                      // angle-integrated intensity over positive direction
   Eigen::Ref<Eigen::MatrixXd> F_ref;             // radiative flux ref
   Eigen::MatrixXd kappa_mat;                     // kappa_mat

   Eigen::Tensor<double, 4> ends;                 // (current ts) These are the nodes psi_{m,i,L}([0]) and psi_{m,i,R}([1])
                                                  // Possible: n + 1/2 predicted, n + 1/2, n+1 predicted, n+1
   Eigen::Tensor<double, 4> prev_ends;            // (n)
   Eigen::Tensor<double, 4> half_ends;            // (n+1/2) Used in BDF2 stepping
   Eigen::VectorXd left_ends, right_ends;
   Eigen::VectorXd rho_vec, kappa_vec, temperature;
   Eigen::VectorXd B;                             // Planckian terms by group, computed in the Correction class
   Eigen::VectorXd dEB;                           // (d(EB)/dE)_g * Delta E_g, computed in the Correction class

   /* Vals needed for iteration */
   bool half_step = true;
   double mu = 0.;
   double bdry_cond = 0.;           // This is specified by ctv::bc_(right/left)_indicator and psi_source
   double local_bdry_prev_it = 0.;  // (n),     This is also given by the upwinding on the previous cell at the previous timestep.
   double local_bdry = 0.;          // (n+1),   This will be given to us by the upwinding on the previous cell.
   double half_local_bdry = 0.;     // (n+1/2), This will be given to us by the upwinding on the previous cell.
   Eigen::MatrixXd _mat, _mat_inverse;
   Eigen::VectorXd _rhs, _res;
   /* ============ end main ============ */

   Correction * correction;
   Eigen::Tensor<double, 3> total_correction; // dir, group, cell
   double logfac;
   double efirst, elast;
   Eigen::VectorXd e_edge, e_ave, de_ave; // Group edge, average energies, and average group widths in kev
   Eigen::MatrixXd energy_discretization; // Left and right edges of energy groups
   Eigen::VectorXd balance;

   /* Equilibrium source terms */
   void computeEquilibriumSources();

   /* ***** Time Stepping ***** */
   void backwardEuler(const int cell, const int scatteredDirIt, const int groupIt, 
                      const double timestep, const double mu);
   void crankNicolson(const int cell, const int scatteredDirIt, const int groupIt, 
                      const double timestep, const double mu);
   void bdf(const int cell, const int scatteredDirIt, const int groupIt, 
            const double timestep, const double mu);


public:
   void generate_group_edges();              // inherited from Correction class
   void generate_group_averages();           // inherited from Correction class
   void fill_energy_bound_arrays();          // inherited from Correction class
   Solver(ParameterHandler & parameter_handler,
          Eigen::Tensor<double, 3>& psi_mat,
          Eigen::Ref<Eigen::MatrixXd> phi,
          Eigen::Ref<Eigen::MatrixXd> F);
          
   void compute_angle_integrated_intensity();
   void compute_positive_angle_integrated_intensity();
   void compute_radiative_flux();
   void compute_balance();
   void get_balance(Eigen::Ref<Eigen::VectorXd> balance) { balance = this->balance; }
   void get_phi_plus(Eigen::Ref<Eigen::MatrixXd> phi_plus) { phi_plus = this->phi_plus; }

   void solve();

   // Getters that may be useful for plotting
   void get_e_ave(Eigen::Ref<Eigen::VectorXd> e_ave) { e_ave = this->e_ave; }
   void compute_group_ends();
   void get_ends(const string side, Eigen::Ref<Eigen::VectorXd> group_ends);
};
} // End namespace rt
#endif // SOLVER