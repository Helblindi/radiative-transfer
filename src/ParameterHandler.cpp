#include "ParameterHandler.h"


// Default constructor
ParameterHandler::ParameterHandler() : param("../prm/default.prm") {
   // TODO: This function is not working
   std::cout << "PH default constructor called.\n";
   cerr << "Not yet implemented.\n";
}

// Non default constructor
ParameterHandler::ParameterHandler(const string filename) : param(filename) {
   if (!param) {
      cerr << "Error in reading file.\n" << endl;
   }
   get_parameters();
}


// Function to grab parameters from file
void ParameterHandler::get_parameters()
{
   M = param.get<int>("M", 2);
   G = param.get<int>("G", 1);
   efirst = param.get<double>("efirst", .1);
   elast = param.get<double>("elast", 10.);
   kappa_grey = param.get<double>("kappa_grey", 1.);
   X = param.get<double>("X", 1.);
   N = param.get<int>("N", 100);
   dx = X / N;

   bc_left_indicator = param.get<int>("bc_left_indicator", 2);
   bc_right_indicator = param.get<int>("bc_right_indicator", 1);
   use_mg_equilib = param.get<bool>("use_mg_equilib", false);

   // Handle possible source conditions via stream stream
   // Only need to explicitly read in source conditions if we
   // are not using equilibrium source conditions
   psi_source.resize(M,G);
   if (!use_mg_equilib)
   {
      string psi_source_s = param.get<string>("psi_source", "no_sources_provided");
      stringstream psi_source_ss(psi_source_s);
      double temp_d;
      int counter = 0;
      while (psi_source_ss >> temp_d)
      {
         // Get rows/cols
         int _m = counter / G, _g = counter % G; 
         psi_source(_m, _g) = temp_d;
         counter++;
      }
   }

   rho = param.get<double>("rho", 1.);
   kappa = param.get<double>("kappa", 1.);
   T = param.get<double>("T", 1.);
   V = param.get<double>("V", 0.);

   use_correction = param.get<bool>("use_correction", false);

   ts_method = param.get<int>("ts_method", 3);
   dt = param.get<double>("dt", 0.00001);
   max_timesteps = param.get<int>("max_timesteps", 1000);
}

void ParameterHandler::get_psi_source(Eigen::Ref<Eigen::MatrixXd> psi_source)
{
   cout << "shape of reference psi source: " << psi_source.rows() << "," << psi_source.cols() << endl;
   cout << "shape of class psi source: " << this->psi_source.rows() << "," << this->psi_source.cols() << endl;
   psi_source = this->psi_source;
}