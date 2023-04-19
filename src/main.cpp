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
                       const double (&rho)[N],
                       const double (&kappa)[N],
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
      _abs += rho[i] * kappa[i] * phi[i] * dx;
      _src += rho[i] * kappa[i] * a * c * pow(temperature[i],4) * dx;
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
// Function to print a matrix
void printMatrix(const double matrix[2][2], const int n) {
   for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
         cout << matrix[i][j] << " ";
      }
      cout << endl;
   }
   cout << endl;
}


double determinant(const double (&matrix)[2][2], int size) {
   if (size == 1) {
      return matrix[0][0];
   }
   double det = 0;
   for (size_t j = 0; j < size; ++j) {
   //   double sub_matrix[size - 1][size - 1];
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
      double a = 1/0;
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
         // cout << "matrix value: " << matrix[i][j] << ", vector val: " << vector[j] << endl;
         sum += matrix[i][j] * vector[j];
      }
      result[i] = sum;
   }
   // cout << "Printing result of matrix vector multiply\n";
   // for (int i = 0; i < 2; i++) {
   //    cout << result[i] << " ";
   // }
   // cout << endl;

}

int main()
{
   double F[N];                 // Radiative flux
   double phi[N];               // angle-integrated intensity
   double psi[M][N];            // solution
   double ends[M][N][2] = {0.}; // These are the nodes psi_{i,L}([0]) and psi_{i,R}([1])
   double rho[N], kappa[N], temperature[N];
   vector<double> x(N);         // mesh

   // Fill constants
   std::fill_n(rho, N, 1.);
   std::fill_n(kappa, N, 1.);
   std::fill_n(temperature, N, T);

   // Initialize the mesh
   for (int i = 0; i < N; i++)
   {
      x[i] = (i + 0.5) * dx;
   }

   /* Vals needed for iteration */
   double mu = 0.;
   double local_bdry = 0.; // This will be given to us by the upwinding on the previous cell.
   double psi_half = 0.;
   double _mat[2][2], _mat_inverse[2][2], _rhs[2], _res[2];

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
         
      }
      else
      {
         switch(bc_left_indicator) {
            case 0: // vacuum
            {
               local_bdry = 0.;
            }
            case 2: // reflective
            {
               // cout << "retrieving reflective BCs at the left boundary for index: " << i << endl;
               
               // implement reflective BCs here
               // Here, mu > 0. For reflective boundary conditions we set
               // Get index j corresponding to direction -mu
               int diff = i - (M / 2); // Difference from center index
               int j = (M / 2) - 1 - diff; // Move in opposite direction to yield index of -mu [-c, -b, -a, a, b, c]
               // cout << "Corresponding index of -mu: " << j << endl;

               local_bdry = ends[j][0][0];
               // cout << "Setting local boundary to: " << ends[j][0][0] << endl;
               // cout << "Compare this to what source conditions give us: " << psi_source[i] << endl;
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
         
      }

      // Iterate over cells in our mesh (The "sweep")
      for (int j = 0; j < N; j++)
      {
         if (mu < 0)
         {
            int cell_j = N - j - 1; // Since we are moving right to left
            // cout << "mu < 0. Iteration: " << j << ", on cell: " << cell_j << endl;
            // cout << "local bdry: " << local_bdry << endl;
            // Fill matrix according to section 3.2 of writeup
            double _temp_val = rho[cell_j] * kappa[cell_j] * dx / 2. - (mu / 2.);
            _mat[0][0] = _temp_val; 
            _mat[0][1] = mu / 2.;
            _mat[1][0] = - mu / 2.;
            _mat[1][1] = _temp_val;

            // Invert matrix
            inverseMatrix(_mat, 2, _mat_inverse);
            // cout << "Printing matrix to be inverted:\n";
            // printMatrix(_mat, 2);
            // cout << "Printing resulting inverse matrix:\n";
            // printMatrix(_mat_inverse, 2);

            // Fill RHS
            // cout << "rho[cell_j]: " << rho[cell_j] << endl;
            // cout << "kappa[cell_j]: " << kappa[cell_j] << endl;
            // cout << "a: " << a << endl;
            // cout << "c: " << c << endl;
            // cout << "T: " << temperature[cell_j] << endl;
            // cout << "T^4: " << pow(temperature[cell_j], 4) << endl;
            // cout << "dx: " << dx << endl;
            // cout << "8*M_PI: " << 8*M_PI << endl;
            // cout << "pow test: " << pow(2,2) << endl;
            _temp_val = rho[cell_j] * kappa[cell_j] * a * c * pow(temperature[cell_j], 4) * dx / (8 * M_PI);
            _rhs[0] = _temp_val;
            _rhs[1] = _temp_val - (mu * local_bdry);

            // cout << "temp val: " << _temp_val << endl;
            // cout << "_rhs[0]: " << _rhs[0] << ", rhs[1]: " << _rhs[1] << endl;

            // Solve 
            matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);
            // Print the result
            // for (int i = 0; i < 2; i++) {
            //    cout << _res[i] << " ";
            // }
            // cout << endl;

            // put the average of val and local boundary here
            psi[i][cell_j] = 0.5*(_res[0] + _res[1]);
            // cout << "psiL: " << _res[0] << ", psi: " << psi[i][cell_j] << ", psiR: " << _res[1] << endl;
            ends[i][cell_j][0] = _res[0];
            ends[i][cell_j][1] = _res[1];

            // Set local boundary for next cell iteration
            local_bdry = _res[0];
         }
         else 
         {
            // Fill matrix according to section 3.2 of writeup
            // cout << "====== mu > 0, Iterating on cell: " << j << endl;
            double _temp_val = rho[j] * kappa[j] * dx / 2. + (mu / 2.);
            _mat[0][0] = _temp_val; 
            _mat[0][1] = mu / 2.;
            _mat[1][0] = - mu / 2.;
            _mat[1][1] = _temp_val;

            // Invert matrix
            inverseMatrix(_mat, 2, _mat_inverse);
            // cout << "Printing resulting inverse matrix:\n";
            // printMatrix(_mat_inverse, 2);

            // Fill RHS
            // cout << "rho[cell_j]: " << rho[j] << endl;
            // cout << "kappa[cell_j]: " << kappa[j] << endl;
            // cout << "a: " << a << endl;
            // cout << "c: " << c << endl;
            // cout << "T: " << temperature[j] << endl;
            // cout << "T^4: " << pow(temperature[j], 4) << endl;
            // cout << "dx: " << dx << endl;
            // cout << "8*M_PI: " << 8*M_PI << endl;
            _temp_val = rho[j] * kappa[j] * a * c * pow(temperature[j], 4) * dx / (8 * M_PI);
            _rhs[0] = _temp_val + (mu * local_bdry);
            _rhs[1] = _temp_val;

            // cout << "temp val: " << _temp_val << endl;
            // cout << "_rhs[0]: " << _rhs[0] << ", rhs[1]: " << _rhs[1] << endl;

            // Solve 
            matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);
            // Print the result
            // for (int i = 0; i < 2; i++) {
            //    cout << _res[i] << " ";
            // }
            // cout << endl;

            // put the average of val and local boundary here
            psi[i][j] = 0.5*(_res[0] + _res[1]);
            // cout << "psiL: " << _res[0] << ", psi: " << psi[i][j] << ", psiR: " << _res[1] << endl;
            ends[i][j][0] = _res[0];
            ends[i][j][1] = _res[1];

            // Set local boundary for next cell iteration
            local_bdry = _res[1];
         }

         
      }
   }

   // Check solution
   // for (int i = 0; i < M; i++)
   // {
   //    for (int j = 0; j < N; j++)
   //    {
   //       cout << "i: " << i << ", j: " << j << endl;
   //       cout << "psi[i][j] - 1/2: " << ends[i][j][0] << ", psi[i][j]: " << psi[i][j] << ", psi[i][j] - 1/2: " << ends[i][j][1] << endl;
   //    }
   // }

   // Compute angle-integrated intensity
   compute_angle_integrated_density(psi, phi);

   // Compute radiative flux
   compute_radiative_flux(psi, F);

   // Output all input quantities
   display_input_quantities();

   // Check and display balance
   double bal = compute_balance(ends, phi, rho, kappa, temperature);
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

   /* Testing matrix functions */
   // TODO: Write ctests to validate
   // double matrix[2][2] = {
   //    {.235, .165},
   //    {.111, -1.5}
   // };
   
   // double inverse[2][2];
   // cout << "Original matrix:\n";
   // printMatrix(matrix, 2);
   // std::cout << "Determinant: " << determinant(matrix, 2) << std::endl;
   // inverseMatrix(matrix, 2, inverse);
   // cout << "Inverse matrix:\n";
   // printMatrix(inverse, 2);
   // cout << "Testing matrix multiply function:\n";
   // double test_vector[2] = {1., 1.};
   // double test_result[2];

   // matrixVectorMultiply(inverse, 2, test_vector, test_result);
   // // Print the result
   // for (int i = 0; i < 2; i++) {
   //    cout << test_result[i] << " ";
   // }
   // cout << endl;

   


   return 1;
}
