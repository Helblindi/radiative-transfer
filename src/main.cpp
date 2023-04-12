#include <iostream>
#include <cmath>
#include <vector>

#include "matplotlibcpp.h"
#include "compile_time_vals.h"
#include "lapack.h"

using namespace std;
using namespace rt;
namespace plt = matplotlibcpp;

double initial_condition(double x)
{
   return -4.*pow((x - 0.5), 2) + 1;
   // return 0.;
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
double compute_balance(const double (&ends)[M][N][2], const double (&phi)[N])
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
      _abs += rho * kappa * phi[i] * dx;
      _src += rho * kappa * a * c * pow(T,4) * dx;
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
void printMatrix(double matrix[2][2], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}


double determinant(double (&matrix)[2][2], int size) {
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
void inverseMatrix(const double matrix[2][2], int n, double inverse[2][2]) {
    // Create the augmented matrix [A | I]
    double augmentedMatrix[2][2*2];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmentedMatrix[i][j] = matrix[i][j];
        }
        augmentedMatrix[i][n+i] = 1.0;
    }
    
    // Perform Gaussian elimination to obtain the inverse matrix
    for (int i = 0; i < n; i++) {
        // Find pivot element
        int pivot = i;
        for (int j = i+1; j < 2*n; j++) {
            if (abs(augmentedMatrix[j][i]) > abs(augmentedMatrix[pivot][i])) {
                pivot = j;
            }
        }
        // Swap rows if necessary
        if (pivot != i) {
            for (int j = 0; j < 2*n; j++) {
                swap(augmentedMatrix[i][j], augmentedMatrix[pivot][j]);
            }
        }
        // Reduce row to have a leading 1
        double factor = augmentedMatrix[i][i];
        for (int j = i; j < 2*n; j++) {
            augmentedMatrix[i][j] /= factor;
        }
        // Zero out elements below the pivot
        for (int j = i+1; j < n; j++) {
            factor = augmentedMatrix[j][i];
            for (int k = i; k < 2*n; k++) {
                augmentedMatrix[j][k] -= factor * augmentedMatrix[i][k];
            }
        }
    }
    
    // Back-substitute to obtain the inverse matrix
    for (int i = n-1; i > 0; i--) {
        for (int j = i-1; j >= 0; j--) {
            double factor = augmentedMatrix[j][i];
            for (int k = i; k < 2*n; k++) {
                augmentedMatrix[j][k] -= factor * augmentedMatrix[i][k];
            }
        }
    }
    
    // Extract the inverse matrix from the augmented matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inverse[i][j] = augmentedMatrix[i][n+j];
        }
    }
}

// Function to multiply a matrix by a vector
void matrixVectorMultiply(double matrix[2][2], int n, double vector[2], double result[2]) {
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
   double F[N] = {0.};          // Radiative flux
   double phi[N] = {0.};        // angle-integrated intensity
   double psi[M][N] = {0.};     // solution
   double ends[M][N][2] = {0.}; // These are the nodes psi_{i,L}([0]) and psi_{i,R}([1])
   vector<double> x(N);         // mesh

   // Initialize the mesh
   for (int i = 0; i < N; i++)
   {
      x[i] = (i + 0.5) * dx;
   }

   /* Vals needed for iteration */
   double mu = 0.;
   double local_bdry = 0.; // This will be given to us by the upwinding on the previous cell.
   double psi_half = 0., val = 0., num = 0., denom = 0.;

   // Iterate over scattered direction (value given by gaussian quadrature)
   for (int i = 0; i < M; i++)
   {
      mu = G_x[i];
      local_bdry = 0.;

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
            case 2: // reflective
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
            num = a * c * pow(T,4) * dx / (4 * M_PI);
            num -= (mu + ( rho * kappa * dx / 2 ) ) * local_bdry;
            denom = ( (rho * kappa * dx) / 2) - mu;
            psi_half = num / denom;
            val = 0.5 * (local_bdry + psi_half);

            // Fill end array according to direction of solution
            ends[i][j][0] = psi_half;
            ends[i][j][1] = local_bdry;

            // Compute the local_bdry condition for next cell
            num = rho * kappa * a * c * pow(T,4) * dx / (8 * M_PI);
            num -= mu * (val - psi_half);
            denom = rho * kappa * dx / 2.;
         }
         else 
         {
            num = a * c * pow(T,4) * dx / (4 * M_PI);
            num += (mu - ( rho * kappa * dx / 2 ) ) * local_bdry;
            denom = mu + ( (rho * kappa * dx) / 2);
            psi_half = num / denom;
            val = 0.5 * (local_bdry + psi_half);

            // Fill end array according to direction of solution
            ends[i][j][0] = local_bdry;
            ends[i][j][1] = psi_half;

            // Compute the local_bdry condition for next cell
            num = rho * kappa * a * c * pow(T,4) * dx / (8 * M_PI);
            num -= mu * (psi_half - val);
            denom = rho * kappa * dx / 2.;
         }

         // put the average of val and local boundary here
         psi[i][j] = val;

         // Set local boundary for next cell iteration
         local_bdry = num / denom;
      }
   }

   // Check solution
   for (int i = 0; i < M; i++)
   {
      for (int j = 0; j < N; j++)
      {
         cout << "i: " << i << ", j: " << j << endl;
         cout << "psi[i][j] - 1/2: " << ends[i][j][0] << ", psi[i][j]: " << psi[i][j] << ", psi[i][j] - 1/2: " << ends[i][j][1] << endl;
      }
   }

   // Compute angle-integrated intensity
   compute_angle_integrated_density(psi, phi);

   // Compute radiative flux
   compute_radiative_flux(psi, F);

   // Check balance
   double bal = compute_balance(ends, phi);
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
   }

   for (int i = 0; i < N; i++)
   {
      y_py[i] = phi[i];
      F_py[i] = F[i];
   }
   plt::named_plot("phi", x, y_py);
   plt::named_plot("Radiative flux", x, F_py);

   plt::legend();
   plt::title("Testing");
   plt::show();

   /* Testing matrix functions */
   // TODO: Write ctests to validate
   double matrix[2][2] = {
      {1., 2.},
      {3., 4.}
   };
   double inverse[2][2];
   cout << "Original matrix:\n";
   printMatrix(matrix, 2);
   std::cout << "Determinant: " << determinant(matrix, 2) << std::endl;
   inverseMatrix(matrix, 2, inverse);
   cout << "Inverse matrix:\n";
   printMatrix(inverse, 2);
   cout << "Testing matrix multiply function:\n";
   double test_vector[2] = {1., 1.};
   double test_result[2];

   matrixVectorMultiply(inverse, 2, test_vector, test_result);
   // Print the result
   for (int i = 0; i < 2; i++) {
      cout << test_result[i] << " ";
   }
   cout << endl;

   


   return 1;
}
