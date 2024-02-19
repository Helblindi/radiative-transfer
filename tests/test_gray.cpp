#include "ParameterHandler.h"
#include "solver.h"
#include "var-config.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace std;
using namespace rt;

template<typename Scalar_, int dim_>
void print_to_file(const string filename, const Eigen::Tensor<Scalar_, dim_>& mat)
{
   ofstream file(filename);
   if (file.is_open())
   {
      file << mat << endl;
   }
   file.close();
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

/***
 * Tests the gray case with the prm/single_group.prm file.  The pass or fail condition
 * of this test depends on the radiative flux.  Since we expect equilibrium for all directions,
 * we check the radiative flux at each cell.
 * 
 * Test will pass if 
 * 
 *          max |F| < 1.E-6
 *           c
 * 
***/
int main()
{
   const string filename = string(TRANSFER_DIR) + "prm/single_group.prm";
   ParameterHandler parameter_handler(filename);

   int M = parameter_handler.get_M(),
       N = parameter_handler.get_N(),
       G = parameter_handler.get_G();

   Eigen::Tensor<double, 3> psi_mat(M, G, N); // solution, angular intensity
   psi_mat.setConstant(0.);
   Eigen::VectorXd x(N);              // Mesh
   Eigen::MatrixXd phi(G, N);            // angle-integrated intensity
   Eigen::MatrixXd F(G, N);              // radiative flux

   // Initialize the mesh
   for (int i = 0; i < N; i++)
   {
      x[i] = (i + 0.5) * parameter_handler.get_dx();
   }

   Solver solver(parameter_handler, psi_mat, phi, F);
   // cout << "Solver initiated.\n";
   solver.solve();

   // Fill phi and F
   solver.compute_angle_integrated_density();
   solver.compute_radiative_flux();

   // Verify balance
   solver.compute_balance();
   Eigen::VectorXd balance(G);
   solver.get_balance(balance);

   cout << "max F: " << F.maxCoeff() << endl;

   // Print to file to be analyzed in python
   print_to_file("gray-test-phi.csv", phi);
   print_to_file("gray-test-psi.csv", psi_mat);
   print_to_file("gray-test-x.csv", x);
   print_to_file("gray-test-F.csv", F);

   if (abs(F.maxCoeff()) < 1.E-6)
   {
      cout << "Gray test passed.\n";
      return 0;
   }
   else
   {
      cout << "Gray test failed.\n";
      return 1;
   } 

   
}