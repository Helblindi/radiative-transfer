#include <iostream>
#include <Eigen/Dense>
 
using Eigen::MatrixXd;

int eigen_practice();
void change_val(Eigen::Ref<Eigen::VectorXd> v);

int main()
{
   int d = eigen_practice();

   return d;
}

void change_val(Eigen::Ref<Eigen::VectorXd> v){  // as per the Eigen guide, one must pass as const
    v(0) = 0; 
    return;
}

int eigen_practice()
{
   // Basic example using eigen
   MatrixXd m(2,2);
   m(0,0) = 3;
   m(1,0) = 2.5;
   m(0,1) = -1;
   m(1,1) = m(1,0) + m(0,1);
   std::cout << "Basic example of filling matrix.\n";
   std::cout << "Here is m:\n" << m << std::endl;

   // Matrix vector multiplication
   Eigen::Matrix2d mat;
   mat << 1, 2,
          3, 4;
   Eigen::VectorXd vec(2);
   vec(0) = 1.;
   vec(1) = 2.;

   std::cout << "Here is mat:\n" << mat << std::endl;
   std::cout << "Here is vec:\n" << vec << std::endl;
   std::cout << "Here is mat*mat:\n" << mat*mat << std::endl;
   std::cout << "Here is mat*vec:\n" << mat*vec << std::endl;
   std::cout << "Here is inv(mat):\n" << mat.inverse() << std::endl;
   std::cout << "Let's multiply mat by itself" << std::endl;
   mat = mat*mat;
   std::cout << "Now mat is mat:\n" << mat << std::endl;

   std::cout << "Lets change the first values in vec\n";
   change_val(vec); // changes first value
   std::cout << "vec after change:\n" << vec << std::endl;

   /***
    * Matrix and Vector class types
    * 
    * MatrixNt for Matrix<type, N, N>. For example, MatrixXi for Matrix<int, Dynamic, Dynamic>.
    * VectorNt for Matrix<type, N, 1>. For example, Vector2f for Matrix<float, 2, 1>.
    * RowVectorNt for Matrix<type, 1, N>. For example, RowVector3d for Matrix<double, 1, 3>.

    * N can be any one of 2, 3, 4, or X (meaning Dynamic).
    * t can be any one of i (meaning int), f (meaning float), d (meaning double), cf (meaning complex<float>), or cd (meaning complex<double>)
    * **/

   return 0;
}