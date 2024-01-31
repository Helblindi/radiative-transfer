#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include "compile_time_vals.h"
#include "solver.h"

using namespace std;
using namespace rt;

// template<typename Scalar_, int rank>
// void shape(const Eigen::Tensor<Scalar_, rank>& x)
// {
//   cout << "( ";  
//   for (int i(0); i<x.NumDimensions; i++){
//       cout << x.dimensions()[i];
//       cout << ",";
//   }
//   cout << ")";  
// }

/* Function to display the parameters used */
void display_input_quantities()
{
   cout << "\n--- Input Parameters ---\n";
   cout << "Quadrature order: " << ctv::M << endl;
   cout << "Slab thickness: " << ctv::X << endl;
   cout << "Number of cells: " << ctv::N << endl;
   cout << "Material density: " << ctv::rho << endl;
   cout << "Absorption opacity: " << ctv::kappa << endl;
   cout << "Material temperature: " << ctv::T << endl;
   cout << "Right boundary condition: ";
   
   // Output boundary conditions
   switch(ctv::bc_right_indicator) {
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
   switch(ctv::bc_left_indicator) {
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



template<typename Scalar_, int rows, int cols>
void print_to_file(const string filename, const Eigen::Matrix<Scalar_, rows, cols>& mat)
{
   ofstream file(filename);
   if (file.is_open())
   {
      file << mat << endl;
   }
   file.close();
}


int main()
{
   // Output all input quantities
   display_input_quantities();

   Eigen::VectorXd phi(ctv::N);            // angle-integrated intensity
   Eigen::MatrixXd psi_mat(ctv::M,ctv::N); // solution, angular intensity
   Eigen::VectorXd x(ctv::N);              // Mesh
   Eigen::VectorXd F(ctv::N);              // radiative flux

   // Initialize the mesh
   for (int i = 0; i < ctv::N; i++)
   {
      x[i] = (i + 0.5) * ctv::dx;
   }

   Solver<ctv::G> solver(psi_mat, phi, F);
   solver.solve();

   // Fill phi and F
   solver.compute_angle_integrated_density();
   solver.compute_radiative_flux();

   // Verify balance
   double balance = solver.compute_balance();

   // Print to file to be analyzed in python
   print_to_file("phi.csv", phi);
   print_to_file("psi.csv", psi_mat);
   print_to_file("x.csv", x);
   print_to_file("F.csv", F);

   return 0;
}
