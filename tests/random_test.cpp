#include <iostream>
#include <Eigen/Dense>
 
using Eigen::MatrixXd;

int eigen_practice();

int main()
{
   int d = eigen_practice();

   return d;
}

int eigen_practice()
{
   // Basic example using eigen
   MatrixXd m(2,2);
   m(0,0) = 3;
   m(1,0) = 2.5;
   m(0,1) = -1;
   m(1,1) = m(1,0) + m(0,1);
   std::cout << m << std::endl;

   return 0;
}