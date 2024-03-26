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


int main(int argc, char **argv)
{
   // Fetch input quantities
   std::string filename;
   if (argc == 2)
   {
      // Assume the used passed in a file to be read
      filename = argv[1];
   }
   else if (argc == 1)
   {
      filename = std::string(TRANSFER_DIR) + "prm/default.prm";
   }
   else
   {
      // Too many arguments passed
      cerr << "Too many command line arguments passed in.\n";
   }
   cout << "filename: " << filename << endl;
   ParameterHandler parameter_handler(filename);
   
   // Output all input quantities
   parameter_handler.display_input_quantities(); 

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
   solver.compute_angle_integrated_intensity();
   solver.compute_radiative_flux();

   // Verify balance
   solver.compute_balance();

   // Get positive angle integrated intensity (just useful for llnl test)
   Eigen::MatrixXd phi_plus(G,N);
   solver.compute_positive_angle_integrated_intensity();
   solver.get_phi_plus(phi_plus);

   // Print to file to be analyzed in python
   print_to_file("phi.csv", phi);
   print_to_file("phi_plus.csv", phi_plus);
   print_to_file("psi.csv", psi_mat);
   print_to_file("x.csv", x);
   print_to_file("F.csv", F);

   // Print info from energy groups
   Eigen::VectorXd e_ave(G);
   solver.get_e_ave(e_ave);
   print_to_file("e_ave.csv", e_ave);

   Eigen::VectorXd left_ends(G), right_ends(G);
   solver.compute_group_ends();
   solver.get_ends("left", left_ends);
   solver.get_ends("right", right_ends);
   print_to_file("left_ends.csv", left_ends);
   print_to_file("right_ends.csv", right_ends);

   return 0;
}
