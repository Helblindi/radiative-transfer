#ifndef SOLVER 
#define SOLVER

#include "compile_time_vals.h"
#include "correction.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace std;

namespace rt
{
template<int num_groups>
class Solver
{
private:
   /* ============ From main ============ */
   Eigen::Tensor<double, 4>& psi_mat_ref;         // solution, angular intensity ref
   Eigen::Ref<Eigen::MatrixXd> phi_ref;           // angle-integrated intensity ref
   Eigen::Ref<Eigen::MatrixXd> F_ref;             // radiative flux ref
   Eigen::MatrixXd kappa_mat;                     // kappa_mat

   Eigen::Tensor<double, 4> ends;      // These are the nodes psi_{i,L}([0]) and psi_{i,R}([1])
   Eigen::Tensor<double, 4> prev_ends;
   Eigen::Tensor<double, 4> half_ends; // For BDF2 time stepping
   Eigen::VectorXd rho_vec, kappa_vec, temperature;

   /* Vals needed for iteration */
   bool half_step = true;
   double mu = 0.;
   double bdry_cond = 0.;           // This is specified by ctv::bc_(right/left)_indicator and psi_source
   double local_bdry = 0.;          // This will be given to us by the upwinding on the previous cell.
   double half_local_bdry = 0.;     // This will be given to us by the upwinding on the previous cell.
   double local_bdry_prev_it = 0.;  // This is also given by the upwinding on the previous cell at the previous timestep.
   double psi_half = 0.;
   Eigen::MatrixXd _mat, _mat_inverse;
   Eigen::VectorXd _rhs, _res;
   /* ============ end main ============ */

   Correction<num_groups> * correction;
   double logfac;
   double efirst=ctv::efirst, elast=ctv::elast;
   Eigen::VectorXd e_edge, e_ave, de_ave; // Group edge, average energies, and average group widths in kev
   Eigen::MatrixXd energy_discretization; // Left and right edges of energy groups
   Eigen::VectorXd balance;

   /* ***** Time Stepping ***** */
   void backwardEuler(const int cell, const int scatteredDirIt, const int groupIt, 
                      const double timestep, const double mu, double &local_bdry);


public:
   void generate_group_edges_and_averages(); // inherited from Correction class
   void fill_energy_bound_arrays();          // inherited from Correction class
   Solver(Eigen::Tensor<double, 4>& psi_mat,
          Eigen::Ref<Eigen::MatrixXd> phi,
          Eigen::Ref<Eigen::MatrixXd> F);
          
   void compute_angle_integrated_density();
   void compute_radiative_flux();
   void compute_balance();
   void get_balance(Eigen::Ref<Eigen::VectorXd> balance) { balance = this->balance; }

   void solve();
};
} // End namespace rt
#endif // SOLVER