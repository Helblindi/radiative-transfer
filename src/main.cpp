#include <iostream>
#include <cmath>
#include <vector>

#include "matplotlibcpp.h"
#include "compile_time_vals.h"
// #include "lapack.h"

using namespace std;
using namespace rt;
namespace plt = matplotlibcpp;


/* Function to display the parameters used */
void display_input_quantities()
{
   cout << "\n--- Input Parameters ---\n";
   cout << "Quadrature order: " << M << endl;
   cout << "Slab thickness: " << X << endl;
   cout << "Number of cells: " << N << endl;
   cout << "Material density: " << rho << endl;
   cout << "Absorption opacity: " << kappa << endl;
   cout << "Material temperature: " << T << endl;
   cout << "Right boundary condition: ";
   
   // Output boundary conditions
   switch(bc_right_indicator) {
      case 0: // vacuum
      {
         cout << "vacuum\n";
         break;
      }
      case 2: // reflective
      {
         cout << "reflective\n";
         break;
      }
      case 1: // source
      {
         cout << "source\n";
         break;
      }
      default:
      {
         cout << "Incorrect boundary conditions provided.\n";
         return;
      }
   }

   cout << "Left boundary condition: ";
   switch(bc_left_indicator) {
      case 0: // vacuum
      {
         cout << "vacuum\n\n";
         break;
      }
      case 2: // reflective
      {
         cout << "reflective\n\n";
         break;
      }
      case 1: // source
      {
         cout << "source\n\n";
         break;
      }
      default:
      {
         cout << "Incorrect boundary conditions provided.\n\n";
         return;
      }
   }
}


template <int M, int N>
void compute_angle_integrated_density(const double (&psi)[M][N], double (&phi)[N])
{
   for (int i = 0; i < N; i++)
   {
      phi[i] = 0.;
      for (int j = 0.; j < M; j++)
      {
         phi[i] += G_w[j] * psi[j][i];
      }
   }
}


template <int M, int N>
void compute_radiative_flux(const double (&psi)[M][N], double (&F)[N])
{
   for (int i = 0; i < N; i++)
   {
      F[i] = 0.;
      for (int j = 0.; j < M; j++)
      {
         F[i] += G_x[j] * G_w[j] * psi[j][i];
      }
   }
}


template <int M, int N>
double compute_balance(const double (&ends)[M][N][2], 
                       const double (&phi)[N],
                       const double (&rho_vec)[N],
                       const double (&kappa_vec)[N],
                       const double (&temperature)[N])
{
   double mu = 0.,
          j_half_minus = 0., 
          j_half_plus = 0.,
          jN_half_minus = 0.,
          jN_half_plus = 0.,
          _abs = 0.,
          _src = 0.,
          _bal = 0.;
   
   for (int i = 0; i < M; i++)
   {
      mu = G_x[i];

      if (mu < 0.)
      {
         j_half_minus -= ends[i][0][0] * mu * G_w[i]; // psi_1/2
         jN_half_minus -= ends[i][N-1][1] * mu * G_w[i]; // psi_N+1/2
      }
      else 
      {
         j_half_plus += ends[i][0][0] * mu * G_w[i]; // psi_1/2
         jN_half_plus += ends[i][N-1][1] * mu * G_w[i]; // psi_N+1/2
      }
   }

   for (int i = 0; i < N; i++)
   {
      _abs += rho_vec[i] * kappa_vec[i] * phi[i] * dx;
      _src += rho_vec[i] * kappa_vec[i] * a * c * pow(temperature[i],4) * dx;
   }

   double sources = jN_half_plus + j_half_minus + _abs;
   double sinks = j_half_plus + jN_half_minus + _src;
   cout << "sources: " << sources << endl;
   cout << "sinks: " << sinks << endl;

   _bal = abs( (jN_half_plus + j_half_minus + _abs) - (j_half_plus + jN_half_minus + _src) ) / (j_half_plus + jN_half_minus + _src);
   return _bal;
}

/*
Functions to go in lapack library
*/


double determinant(const double (&matrix)[2][2], int size) {
   if (size == 1) {
      return matrix[0][0];
   }
   double det = 0;
   for (size_t j = 0; j < size; ++j) {
      double sub_matrix[2][2];
      for (size_t i = 1; i < size; ++i) {
         for (size_t k = 0; k < size - 1; ++k) {
               if (k < j) {
                  sub_matrix[i - 1][k] = matrix[i][k];
               }
               else {
                  sub_matrix[i - 1][k] = matrix[i][k + 1];
               }
         }
      }
      double sub_det = determinant(sub_matrix, size - 1);
      double sign = ((j % 2) == 0) ? 1 : -1;
      det += sign * matrix[0][j] * sub_det;
   }
   return det;
}


// Function to compute the inverse matrix of a square matrix
void inverseMatrix(const double (&matrix)[2][2], const int n, double (&inverse)[2][2]) {
   double det = determinant(matrix, n);
   if (abs(det) < 1e-10)
   {
      std::cout << "Cannot divide by 0.\n";
      assert(false);
   }
   inverse[0][0] = matrix[1][1] / det;
   inverse[0][1] = - matrix[0][1] / det;
   inverse[1][0] = - matrix[1][0] / det;
   inverse[1][1] = matrix[0][0] / det;
}

// Function to multiply a matrix by a vector
void matrixVectorMultiply(const double (&matrix)[2][2], const int n, const double (&vector)[2], double (&result)[2]) {
   for (int i = 0; i < n; i++) {
      double sum = 0;
      for (int j = 0; j < n; j++) {
         sum += matrix[i][j] * vector[j];
      }
      result[i] = sum;
   }
}

int main()
{
   double F[N];                 // Radiative flux
   double phi[N];               // angle-integrated intensity
   double psi[M][N] = {};       // solution
   double ends[M][N][2] = {};   // These are the nodes psi_{i,L}([0]) and psi_{i,R}([1])
   double prev_ends[M][N][2] = {};
   double half_ends[M][N][2] = {}; // For BDF2 time stepping
   double rho_vec[N], kappa_vec[N], temperature[N];
   vector<double> x(N);         // mesh

   // Fill constants
   std::fill_n(rho_vec, N, rho);
   std::fill_n(kappa_vec, N, kappa);
   std::fill_n(temperature, N, T);

   // Initialize the mesh
   for (int i = 0; i < N; i++)
   {
      x[i] = (i + 0.5) * dx;
   }

   /* Vals needed for iteration */
   bool half_step = true;
   double mu = 0.;
   double local_bdry = 0.;          // This will be given to us by the upwinding on the previous cell.
   double half_local_bdry = 0.;     // This will be given to us by the upwinding on the previous cell.
   double local_bdry_prev_it = 0.;  // This is also given by the upwinding on the previous cell at the previous timestep.
   double psi_half = 0.;
   double _mat[2][2], _mat_inverse[2][2], _rhs[2], _res[2];

   for (int _it = 0; _it < _max_timesteps; _it++)
   {
      cout << "============= Timestep: " << _it << " =============" << endl;
      std::copy(&ends[0][0][0], &ends[0][0][0]+M*N*2, &prev_ends[0][0][0]);

      // Iterate over scattered direction (value given by gaussian quadrature)
      for (int i = 0; i < M; i++)
      {
         mu = G_x[i];
         local_bdry = 0.;
         cout << "Beginning sweep for mu: " << mu << endl;

         // Initialize boundary conditions
         if (mu < 0.)
         {
            switch(bc_right_indicator) {
               case 0: // vacuum
               {
                  local_bdry = 0.;
                  break;
               }
               case 2: // reflective
               {
                  local_bdry = 0.;
                  // TODO: Implement check once finished sweep
                  break;
               }
               case 1: // source
               {
                  local_bdry = psi_source[i];
                  break;
               }
               default:
               {
                  cout << "Incorrect boundary conditions provided.\n";
                  return -1;
               }
            }
            half_local_bdry = local_bdry;
            local_bdry_prev_it = local_bdry;
         }
         else
         {
            switch(bc_left_indicator) {
               case 0: // vacuum
               {
                  local_bdry = 0.;
               }
               case 1: // source
               {
                  local_bdry = psi_source[i];
                  break;
               }
               case 2: // reflective
               {
                  // Get index j corresponding to direction -mu
                  int diff = i - (M / 2); // Difference from center index
                  int m_neg = (M / 2) - 1 - diff; // Move in opposite direction to yield index of -mu [-c, -b, -a, a, b, c]

                  local_bdry = ends[m_neg][0][0];
                  break;
               }
               default:
               {
                  cout << "Incorrect boundary conditions provided.\n";
                  return -1;
               }
            }
            half_local_bdry = local_bdry;
            local_bdry_prev_it = local_bdry;
         }

         // Iterate over cells in our mesh (The "sweep")
         for (int j = 0; j < N; j++)
         {
            if (mu < 0)
            {
               int cell_j = N - j - 1; // Since we are moving right to left
               // Fill matrix and rhs according to timestep selected
               switch(ts_method) {
                  case 1: // BE
                  {
                     double const_A = 1. + c*dt*rho_vec[cell_j] * kappa_vec[cell_j];
                     double const_B = c*dt*mu;
                     double const_C = c*dt*rho_vec[cell_j] * kappa_vec[cell_j] * a * c / (4 * M_PI);

                     double _temp_val = (const_A * dx - const_B) / 2.;
                     _mat[0][0] = _temp_val; 
                     _mat[0][1] = const_B / 2.;
                     _mat[1][0] = - const_B / 2.;
                     _mat[1][1] = _temp_val;

                     _temp_val = const_C * dx * pow(temperature[cell_j], 4) / 2.;
                     _rhs[0] = _temp_val + dx * ends[i][cell_j][0] / 2.;
                     _rhs[1] = _temp_val - (const_B * local_bdry) + dx * ends[i][cell_j][1] / 2.;

                     // Invert matrix
                     inverseMatrix(_mat, 2, _mat_inverse);

                     // Solve 
                     matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);

                     // put the average of val and local boundary here
                     psi[i][cell_j] = 0.5*(_res[0] + _res[1]);

                     ends[i][cell_j][0] = _res[0];
                     ends[i][cell_j][1] = _res[1];

                     // Set local boundary for next cell iteration
                     local_bdry = _res[0];

                     break;
                  }
                  case 2: // BDF2
                  {
                     // CN Step
                     double const_A = c * dt * mu / 4.;
                     double const_B = 1. + (c * dt * rho_vec[cell_j] * kappa_vec[cell_j]) / 4.;
                     double const_C = 1. - (c * dt * rho_vec[cell_j] * kappa_vec[cell_j]) / 4.;
                     double const_D = c * dt * a * c / (8 * M_PI);

                     double _temp_val = (const_B * dx - const_A) / 2.;
                     _mat[0][0] = _temp_val; 
                     _mat[0][1] = const_A / 2.;
                     _mat[1][0] = - const_A / 2.;
                     _mat[1][1] = _temp_val;

                     _temp_val = const_D * dx * rho_vec[cell_j] * kappa_vec[cell_j] * pow(temperature[cell_j], 4) / 2.;

                     _rhs[0] = _temp_val + ((const_C * dx + const_A) / 2.) * ends[i][cell_j][0] - (const_A / 2.) * ends[i][cell_j][1];
                     _rhs[1] = _temp_val + ((const_C * dx + const_A) / 2.) * ends[i][cell_j][1] + (const_A / 2.) * ends[i][cell_j][0] - const_A * (half_local_bdry + local_bdry_prev_it);

                     // Solve CN
                     inverseMatrix(_mat, 2, _mat_inverse);
                     matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);

                     half_ends[i][cell_j][0] = _res[0];
                     half_ends[i][cell_j][1] = _res[1];

                     // BDF Step
                     const_A = 1. + (c * dt * rho_vec[cell_j] * kappa_vec[cell_j] / 12.);
                     const_B = mu * c * dt / 12.;
                     const_C = 1. -  (c * dt * rho_vec[cell_j] * kappa_vec[cell_j] / 3.);
                     const_D = c * dt * rho_vec[cell_j] * kappa_vec[cell_j] / 12.;
                     double const_E = c * dt * a * c / (8. * M_PI);

                     _temp_val = (const_A * dx - const_B) / 2.;
                     _mat[0][0] = _temp_val; 
                     _mat[0][1] = const_B / 2.;
                     _mat[1][0] = - const_B / 2.;
                     _mat[1][1] = _temp_val;

                     _temp_val = (const_E * dx * rho_vec[cell_j] * kappa_vec[cell_j] / 2.) * pow(temperature[cell_j], 4);
                     _rhs[0] = _temp_val + ((const_C * dx + 4. * const_B) / 2.) * half_ends[i][cell_j][0] - 2. * const_B * half_ends[i][cell_j][1];
                     _rhs[0] += ((const_B - const_D * dx) / 2.) * ends[i][cell_j][0] - (const_B / 2.) * ends[i][cell_j][1];
                     _rhs[1] = _temp_val + ((const_C * dx + 4. * const_B) / 2.) * half_ends[i][cell_j][1] + 2. * const_B * half_ends[i][cell_j][0];
                     _rhs[1] += ((const_B - const_D * dx) / 2.) * ends[i][cell_j][1] + (const_B / 2.) * ends[i][cell_j][0];
                     _rhs[1] -= const_B * (local_bdry + local_bdry_prev_it + 4. *  half_local_bdry);

                     // Invert matrix
                     inverseMatrix(_mat, 2, _mat_inverse);

                     // Solve 
                     matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);

                     // put the average of val and local boundary here
                     psi[i][cell_j] = 0.5*(_res[0] + _res[1]);

                     ends[i][cell_j][0] = _res[0];
                     ends[i][cell_j][1] = _res[1];

                     // Set local boundary for next cell iteration
                     local_bdry = _res[0];
                     half_local_bdry = half_ends[i][cell_j][0];
                     local_bdry_prev_it = prev_ends[i][cell_j][0];
                     
                     break;
                  } // End BDF2 case
                  case 3: // CN
                  {
                     // CN Step
                     double const_A = c * dt * mu / 4.;
                     double const_B = 1. + (c * dt * rho_vec[cell_j] * kappa_vec[cell_j]) / 4.;
                     double const_C = 1. - (c * dt * rho_vec[cell_j] * kappa_vec[cell_j]) / 4.;
                     double const_D = c * dt * a * c / (8 * M_PI);

                     double _temp_val = (const_B * dx - const_A) / 2.;
                     _mat[0][0] = _temp_val; 
                     _mat[0][1] = const_A / 2.;
                     _mat[1][0] = - const_A / 2.;
                     _mat[1][1] = _temp_val;

                     _temp_val = const_D * dx * rho_vec[cell_j] * kappa_vec[cell_j] * pow(temperature[cell_j], 4) / 2.;

                     _rhs[0] = _temp_val + ((const_C * dx + const_A) / 2.) * ends[i][cell_j][0] - (const_A / 2.) * ends[i][cell_j][1];
                     _rhs[1] = _temp_val + ((const_C * dx + const_A) / 2.) * ends[i][cell_j][1] + (const_A / 2.) * ends[i][cell_j][0] - const_A * (half_local_bdry + local_bdry_prev_it);

                     // Invert matrix
                     inverseMatrix(_mat, 2, _mat_inverse);

                     // Solve 
                     matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);

                     // put the average of val and local boundary here
                     psi[i][cell_j] = 0.5*(_res[0] + _res[1]);

                     ends[i][cell_j][0] = _res[0];
                     ends[i][cell_j][1] = _res[1];

                     // Set local boundary for next cell iteration
                     local_bdry = _res[0];

                     break;
                  } // End CN
                  case 4: // BDF2 Morel-Lou
                  {
                     if (half_step)
                     {
                        // CN Step
                        double const_A = c * dt * mu / 4.;
                        double const_B = 1. + (c * dt * rho_vec[cell_j] * kappa_vec[cell_j]) / 4.;
                        double const_C = 1. - (c * dt * rho_vec[cell_j] * kappa_vec[cell_j]) / 4.;
                        double const_D = c * dt * a * c / (8 * M_PI);

                        double _temp_val = (const_B * dx - const_A) / 2.;
                        _mat[0][0] = _temp_val; 
                        _mat[0][1] = const_A / 2.;
                        _mat[1][0] = - const_A / 2.;
                        _mat[1][1] = _temp_val;

                        _temp_val = const_D * dx * rho_vec[cell_j] * kappa_vec[cell_j] * pow(temperature[cell_j], 4) / 2.;

                        _rhs[0] = _temp_val + ((const_C * dx + const_A) / 2.) * ends[i][cell_j][0] - (const_A / 2.) * ends[i][cell_j][1];
                        _rhs[1] = _temp_val + ((const_C * dx + const_A) / 2.) * ends[i][cell_j][1] + (const_A / 2.) * ends[i][cell_j][0] - const_A * (half_local_bdry + local_bdry_prev_it);

                        // Solve CN
                        inverseMatrix(_mat, 2, _mat_inverse);
                        matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);

                        half_ends[i][cell_j][0] = _res[0];
                        half_ends[i][cell_j][1] = _res[1];

                        half_local_bdry = half_ends[i][cell_j][0];
                        local_bdry_prev_it = prev_ends[i][cell_j][0];
                     }
                     else
                     {
                        // BDF Step
                        double const_A = 1. + (c * dt * rho_vec[cell_j] * kappa_vec[cell_j] / 3.);
                        double const_B = mu * c * dt / 3.;
                        double const_C = 1. -  (c * dt * rho_vec[cell_j] * kappa_vec[cell_j] / 12.);
                        double const_D = c * dt * rho_vec[cell_j] * kappa_vec[cell_j] / 12.;
                        double const_E = c * dt * a * c / (8. * M_PI);

                        double _temp_val = (const_A * dx - const_B) / 2.;
                        _mat[0][0] = _temp_val; 
                        _mat[0][1] = const_B / 2.;
                        _mat[1][0] = - const_B / 2.;
                        _mat[1][1] = _temp_val;

                        _temp_val = (const_E * dx * rho_vec[cell_j] * kappa_vec[cell_j] / 2.) * pow(temperature[cell_j], 4);
                        _rhs[0] = _temp_val + ((const_C * dx + (const_B / 4.)) / 2.) * half_ends[i][cell_j][0] - (const_B / 8.) * half_ends[i][cell_j][1];
                        _rhs[0] += (((const_B / 4.) - const_D * dx) / 2.) * ends[i][cell_j][0] - (const_B / 8.) * ends[i][cell_j][1];
                        _rhs[1] = _temp_val + ((const_C * dx + (const_B / 4.)) / 2.) * half_ends[i][cell_j][1] + (const_B / 8.) * half_ends[i][cell_j][0];
                        _rhs[1] += (((const_B / 4.) - const_D * dx) / 2.) * ends[i][cell_j][1] + (const_B / 8.) * ends[i][cell_j][0];
                        _rhs[1] -= const_B * (local_bdry + (local_bdry_prev_it / 4.) + (half_local_bdry / 4.));

                        // Invert matrix
                        inverseMatrix(_mat, 2, _mat_inverse);

                        // Solve 
                        matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);

                        // put the average of val and local boundary here
                        psi[i][cell_j] = 0.5*(_res[0] + _res[1]);

                        ends[i][cell_j][0] = _res[0];
                        ends[i][cell_j][1] = _res[1];

                        // Set local boundary for next cell iteration
                        local_bdry = _res[0];
                        half_local_bdry = half_ends[i][cell_j][0];
                        local_bdry_prev_it = prev_ends[i][cell_j][0];
                     }
                     
                     break;
                  } // End BDF2 case
                  default:
                  {
                     cout << "Incorrect timestepping method provided.\n";
                     return -1;
                  }
               }
            }
            else 
            { // mu > 0
               // Fill matrix and rhs according to timestep selected
               switch(ts_method) {
                  case 1: // BE
                  {
                     double const_A = 1. + c*dt*rho_vec[j] * kappa_vec[j];
                     double const_B = c*dt*mu;
                     double const_C = c*dt*rho_vec[j] * kappa_vec[j] * a * c / (4 * M_PI);

                     double _temp_val = (const_A * dx + const_B) / 2.;
                     _mat[0][0] = _temp_val; 
                     _mat[0][1] = const_B / 2.;
                     _mat[1][0] = - const_B / 2.;
                     _mat[1][1] = _temp_val;

                     _temp_val = const_C * dx * pow(temperature[j], 4) / 2.;
                     _rhs[0] = _temp_val + (const_B * local_bdry) + dx * ends[i][j][0] / 2.;
                     _rhs[1] = _temp_val + dx * ends[i][j][1] / 2.;

                     // Invert matrix
                     inverseMatrix(_mat, 2, _mat_inverse);

                     // Solve 
                     matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);

                     // put the average of val and local boundary here
                     psi[i][j] = 0.5*(_res[0] + _res[1]);

                     ends[i][j][0] = _res[0];
                     ends[i][j][1] = _res[1];

                     // Set local boundary for next cell iteration
                     local_bdry = _res[1];

                     break;
                  }
                  case 2: // BDF2
                  {
                     // CN Step
                     double const_A = c * dt * mu / 4.;
                     double const_B = 1. + (c * dt * rho_vec[j] * kappa_vec[j]) / 4.;
                     double const_C = 1. - (c * dt * rho_vec[j] * kappa_vec[j]) / 4.;
                     double const_D = c * dt * a * c / (8 * M_PI);

                     double _temp_val = (const_B * dx + const_A) / 2.;
                     _mat[0][0] = _temp_val; 
                     _mat[0][1] = const_A / 2.;
                     _mat[1][0] = - const_A / 2.;
                     _mat[1][1] = _temp_val;

                     _temp_val = const_D * dx * rho_vec[j] * kappa_vec[j] * pow(temperature[j], 4) / 2.;

                     _rhs[0] = _temp_val + ((const_C * dx - const_A) / 2.) * ends[i][j][0] - (const_A / 2.) * ends[i][j][1] + const_A * (half_local_bdry + local_bdry_prev_it);
                     _rhs[1] = _temp_val + ((const_C * dx - const_A) / 2.) * ends[i][j][1] + (const_A / 2.) * ends[i][j][0];

                     // Solve CN
                     inverseMatrix(_mat, 2, _mat_inverse);
                     matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);

                     half_ends[i][j][0] = _res[0];
                     half_ends[i][j][1] = _res[1];

                     // BDF Step
                     const_A = 1. + (c * dt * rho_vec[j] * kappa_vec[j] / 12.);
                     const_B = mu * c * dt / 12.;
                     const_C = 1. -  (c * dt * rho_vec[j] * kappa_vec[j] / 3.);
                     const_D = c * dt * rho_vec[j] * kappa_vec[j] / 12.;
                     double const_E = c * dt * a * c / (8. * M_PI);

                     _temp_val = (const_A * dx + const_B) / 2.;
                     _mat[0][0] = _temp_val; 
                     _mat[0][1] = const_B / 2.;
                     _mat[1][0] = - const_B / 2.;
                     _mat[1][1] = _temp_val;

                     _temp_val = (const_E * dx * rho_vec[j] * kappa_vec[j] / 2.) * pow(temperature[j], 4);
                     _rhs[0] = _temp_val + ((const_C * dx - 4. * const_B) / 2.) * half_ends[i][j][0] - 2. * const_B * half_ends[i][j][1];
                     _rhs[0] += -1. * ((const_B + const_D * dx) / 2.) * ends[i][j][0] - (const_B / 2.) * ends[i][j][1];
                     _rhs[0] += const_B * (local_bdry + local_bdry_prev_it + 4. * half_local_bdry);
                     _rhs[1] = _temp_val + ((const_C * dx - 4. * const_B) / 2.) * half_ends[i][j][1] + 2. * const_B * half_ends[i][j][0];
                     _rhs[1] += -1. * ((const_B + const_D * dx) / 2.) * ends[i][j][1] + (const_B / 2.) * ends[i][j][0];

                     // Solve BDF2
                     inverseMatrix(_mat, 2, _mat_inverse);

                     // Solve 
                     matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);

                     // put the average of val and local boundary here
                     psi[i][j] = 0.5*(_res[0] + _res[1]);

                     ends[i][j][0] = _res[0];
                     ends[i][j][1] = _res[1];

                     // Set local boundary for next cell iteration
                     local_bdry = _res[1];
                     half_local_bdry = half_ends[i][j][1];
                     local_bdry_prev_it = prev_ends[i][j][1];

                     break;
                  } // End BDF2 case
                  case 3: // CN
                  {
                     // CN Step
                     double const_A = c * dt * mu / 4.;
                     double const_B = 1. + (c * dt * rho_vec[j] * kappa_vec[j]) / 4.;
                     double const_C = 1. - (c * dt * rho_vec[j] * kappa_vec[j]) / 4.;
                     double const_D = c * dt * a * c / (8 * M_PI);

                     double _temp_val = (const_B * dx + const_A) / 2.;
                     _mat[0][0] = _temp_val; 
                     _mat[0][1] = const_A / 2.;
                     _mat[1][0] = - const_A / 2.;
                     _mat[1][1] = _temp_val;

                     _temp_val = const_D * dx * rho_vec[j] * kappa_vec[j] * pow(temperature[j], 4) / 2.;

                     _rhs[0] = _temp_val + ((const_C * dx - const_A) / 2.) * ends[i][j][0] - (const_A / 2.) * ends[i][j][1] + const_A * (half_local_bdry + local_bdry_prev_it);
                     _rhs[1] = _temp_val + ((const_C * dx - const_A) / 2.) * ends[i][j][1] + (const_A / 2.) * ends[i][j][0];

                     // Solve matrix
                     inverseMatrix(_mat, 2, _mat_inverse);

                     // Solve 
                     matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);

                     // put the average of val and local boundary here
                     psi[i][j] = 0.5*(_res[0] + _res[1]);

                     ends[i][j][0] = _res[0];
                     ends[i][j][1] = _res[1];

                     // Set local boundary for next cell iteration
                     local_bdry = _res[1];
                     
                     break;
                  } // End CN
                  case 4: // BDF2 Morel-Lou
                  {
                     if (half_step)
                     {
                        // CN Step
                        double const_A = c * dt * mu / 4.;
                        double const_B = 1. + (c * dt * rho_vec[j] * kappa_vec[j]) / 4.;
                        double const_C = 1. - (c * dt * rho_vec[j] * kappa_vec[j]) / 4.;
                        double const_D = c * dt * a * c / (8 * M_PI);

                        double _temp_val = (const_B * dx + const_A) / 2.;
                        _mat[0][0] = _temp_val; 
                        _mat[0][1] = const_A / 2.;
                        _mat[1][0] = - const_A / 2.;
                        _mat[1][1] = _temp_val;

                        _temp_val = const_D * dx * rho_vec[j] * kappa_vec[j] * pow(temperature[j], 4) / 2.;

                        _rhs[0] = _temp_val + ((const_C * dx - const_A) / 2.) * ends[i][j][0] - (const_A / 2.) * ends[i][j][1] + const_A * (half_local_bdry + local_bdry_prev_it);
                        _rhs[1] = _temp_val + ((const_C * dx - const_A) / 2.) * ends[i][j][1] + (const_A / 2.) * ends[i][j][0];

                        // Solve CN
                        inverseMatrix(_mat, 2, _mat_inverse);
                        matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);

                        half_ends[i][j][0] = _res[0];
                        half_ends[i][j][1] = _res[1];

                        half_local_bdry = half_ends[i][j][1];
                        local_bdry_prev_it = prev_ends[i][j][1];
                     }
                     else
                     {
                        // BDF Step
                        double const_A = 1. + (c * dt * rho_vec[j] * kappa_vec[j] / 3.);
                        double const_B = mu * c * dt / 3.;
                        double const_C = 1. -  (c * dt * rho_vec[j] * kappa_vec[j] / 12.);
                        double const_D = c * dt * rho_vec[j] * kappa_vec[j] / 12.;
                        double const_E = c * dt * a * c / (8. * M_PI);

                        double _temp_val = (const_A * dx + const_B) / 2.;
                        _mat[0][0] = _temp_val; 
                        _mat[0][1] = const_B / 2.;
                        _mat[1][0] = - const_B / 2.;
                        _mat[1][1] = _temp_val;

                        _temp_val = (const_E * dx * rho_vec[j] * kappa_vec[j] / 2.) * pow(temperature[j], 4);
                        _rhs[0] = _temp_val + ((const_C * dx - const_B/4.) / 2.) * half_ends[i][j][0] - (const_B / 8.) * half_ends[i][j][1];
                        _rhs[0] += -1. * (((const_B / 4.) + const_D * dx) / 2.) * ends[i][j][0] - (const_B / 8.) * ends[i][j][1];
                        _rhs[0] += const_B * (local_bdry + (local_bdry_prev_it / 4.) + (half_local_bdry / 4.));
                        _rhs[1] = _temp_val + ((const_C * dx - (const_B / 4.)) / 2.) * half_ends[i][j][1] + (const_B / 8.) * half_ends[i][j][0];
                        _rhs[1] += -1. * (((const_B / 4.) + const_D * dx) / 2.) * ends[i][j][1] + (const_B / 8.) * ends[i][j][0];

                        // Solve BDF2
                        inverseMatrix(_mat, 2, _mat_inverse);

                        // Solve 
                        matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);

                        // put the average of val and local boundary here
                        psi[i][j] = 0.5*(_res[0] + _res[1]);

                        ends[i][j][0] = _res[0];
                        ends[i][j][1] = _res[1];

                        // Set local boundary for next cell iteration
                        local_bdry = _res[1];
                        half_local_bdry = half_ends[i][j][1];
                        local_bdry_prev_it = prev_ends[i][j][1];
                     }

                     break;
                  } // End BDF2 case
                  default:
                  {
                     cout << "Incorrect timestepping method provided.\n";
                     return -1;
                  }
               }
            }
         }
         if (mu < 0)
         {
            cout << "Final phi for mu " << mu << " is: " << ends[i][0][0] << endl;
         }
         else
         {
            cout << "Final phi for mu " << mu << " is: " << ends[i][N-1][1] << endl;
         }
         cout << "Corresponding source condition was: " <<  psi_source[i] << endl;
      } // End scattered direction loop
      if (ts_method == 2 || ts_method == 4)
      {
         half_step = !half_step;
      }
   } // End time loop

   // Compute angle-integrated intensity
   compute_angle_integrated_density(psi, phi);

   // Compute radiative flux
   compute_radiative_flux(psi, F);

   // Output all input quantities
   display_input_quantities();

   // Check and display balance
   double bal = compute_balance(ends, phi, rho_vec, kappa_vec, temperature);
   cout << "balance: " << bal << endl;

   // Plot each contribution using matplotlib
   plt::figure_size(1200, 780);

   // Create vectors for python to use
   vector<double> y_py(N);
   vector<double> F_py(N);

   for (int i = 0; i < M; i++)
   {
      mu = G_x[i];
      for (int j = 0; j < N; j++)
      {
         y_py[j] = psi[i][j];
      }
      string plot_tag = "Psi for mu = " + to_string(mu);
      plt::named_plot(plot_tag, x, y_py);
      // plt::scatter(x, y_py);
   }

   for (int i = 0; i < N; i++)
   {
      y_py[i] = phi[i];
      F_py[i] = F[i];
   }
   plt::named_plot("phi", x, y_py);
   plt::named_plot("Radiative flux", x, F_py);

   // plt::axis("on");
   plt::xlabel("x");
   plt::ylabel("y");

   plt::legend();
   plt::title("Testing");
   plt::show();

   return 1;
}
