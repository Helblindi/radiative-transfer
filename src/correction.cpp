#include "correction.h"


namespace rt 
{

// Planck function
// Photon energy E and temperature T are in units of keV
// The returned value from this function has units
// jk/(cm^2-sh-keV-steradian)
template<int num_groups>
double Correction<num_groups>::pf(double E, double T)
{
   double h = Constants::PLANCK_CONSTANT;
   double c = Constants::SPEED_OF_LIGHT;
   double k = Constants::BOLTZMANN_CONSTANT_JPK;
   
   double denom = pow(h,3)*pow(c,2)*(exp(E/T)-1.0);
   assert(denom != 0 && "Do not divide by 0!\n");

   double val = k*pow(E,3)/denom;
   return val;
}


template<int num_groups>
void Correction<num_groups>::generate_planck_integrals()
{
   // Generate Planck Integrals (keV/(cm^2-SH))
   planck.get_Planck(T, energy_discretization_ref, B, dBdT);

   // Change Planck integral units to JK/(CM^2-SH)
   for(int g = 0; g < num_groups; g++)
   {
      B(g)    = kcon*B(g);
      dBdT(g) = kcon*dBdT(g);
   }
}


template<int num_groups>
bool Correction<num_groups>::validate_planck_integrals()
{
   // Check for proper sums of Planck integrals
   double bsum = 0., dbsum = 0.;
  
   for (int g=0; g < num_groups; g++)
   {
      bsum += B(g);
      dbsum += dBdT(g);
   }
   double acT4=ac*pow(T,4);
   double dacT4=4.0*ac*pow(T,3);

   // These quantites must be the same
   if (abs(acT4 - bsum) > ctv::validation_tolerance || 
       abs(dacT4 - dbsum) > ctv::validation_tolerance)
   {
      cout << "acT^4 = " << acT4 << " B sum = " << bsum << endl;
      cout << "4acT^3 = " << dacT4 << " dBdT sum = " << dbsum << "\n" << endl;

      return false;
   }

   return true;
}


template<int num_groups>
void Correction<num_groups>::generate_multigroup_opacities()
{
   double kappa_nfac,                         // Opacity normalization factor
          kappa_grey,                         // Grey opacity used to define total emission at a given temperature
          tmp1, tmp2;                         // Temporary variables
   
   // Compute unnormalized opacities and group-centered opacities
   for(int g = 0; g < num_groups; g++)
   {
      tmp1 = 2.0*fourpi*kcon*T/(pow(hcon,3)*pow(ccon,2));
      tmp2 = exp(-e_edge_ref(g)/T) - exp(-e_edge_ref(g+1)/T); 
      ukappa(g) = tmp1*tmp2/B(g);
      ckappa(g) = (1.0-exp(-e_ave_ref(g)/T))/pow(e_ave_ref(g),3);
   }

   // Compute normalization constant
   double emis_tot = 0.0,
          acT4=ac*pow(T,4);

   for(int g = 0; g < num_groups; g++)
   {
      emis_tot += ukappa(g)*B(g);
   }
   kappa_nfac = acT4/emis_tot;

   // Compute final opacities and emission spectrum
   for(int g = 0;g < num_groups; g++)
   {
      kappa_ref(g) = kappa_grey*kappa_nfac*ukappa(g);
      emis_spec(g) = kappa_ref(g)*B(g);
  
      // cout << left << setw(7)  << g 
      //      << left << setw(14) << e_ave_ref(g)
      //      << left << setw(15) << kappa_ref(g)
      //      << left << setw(15) << emis_spec(g) << endl; 
   }
   // cout << "\n" << endl;
}


template<int num_groups>
bool Correction<num_groups>::validate_emission()
{
   // Check for proper emission sum
   double acT4=ac*pow(T,4),
          sigacT4 = kappa_grey*acT4;
  
   double emis_tot = 0.0;
   for(int g = 0;g < num_groups; g++)
   {
      emis_tot += kappa_ref(g)*B(g);
   }

   

   // total emission thould be equal to sigacT4
   if (abs(emis_tot - sigacT4) > ctv::validation_tolerance)
   {
      cout << "Total Emission = " << emis_tot 
           << " kappa_ref*acT^4 = " << sigacT4 << "\n"  << endl;

      return false;
   }

   return true;
}


template<int num_groups>
void Correction<num_groups>::compute_group_edge_opacities()
{
   // Compute group-edge opacities
   double wgt_L, wgt_R;       // Left and right weights for edge opacity generation
                                     
   kappa_edge(0) = kappa_ref(0);
  
   for(int g = 1; g < num_groups; g++)
   {
  	   wgt_L = (e_ave_ref(g)-e_edge_ref(g))/(e_ave_ref(g) - e_ave_ref(g-1));
  	   wgt_R = (e_edge_ref(g)-e_ave_ref(g-1))/(e_ave_ref(g) - e_ave_ref(g-1));
      kappa_edge(g) = kappa_ref(g-1)*wgt_L + kappa_ref(g)*wgt_R;
   }  

   kappa_edge(num_groups) = kappa_ref(num_groups-1);

   // Print group edge opacities
   // cout << left << setw(12) << "Edge Index"
   //      << left << setw(14) << "Energy" 
   //      << left << setw(14) << "Opacity" << endl;
       
   // cout << left << setw(12) << "---------"
   //      << left << setw(14) << "------" 
   //      << left << setw(14) << "-------" << endl;
       
   
   // for(int g = 0;g < num_groups+1; g++)
   // {
   //    cout << left << setw(12)  << g 
   //         << left << setw(14) << e_edge_ref(g) 
   //         << left << setw(14) << kappa_edge(g) << endl;
   // }
   
   // cout << "\n" << endl;
}


template<int num_groups>
void Correction<num_groups>::compute_components_of_correction_source()
{
   // Components of correction terms
   dEB(0) = e_edge_ref(1)*pf(e_edge_ref(1),T);

   if (num_groups > 1) // Only necessary for non-grey situations
   {
      for(int g = 1;g < num_groups-1; g++)
      {
         dEB(g) = e_edge_ref(g+1)*pf(e_edge_ref(g+1),T) - e_edge_ref(g)*pf(e_edge_ref(g),T);
      }
      dEB(num_groups-1) = - e_edge_ref(num_groups-1)*pf(e_edge_ref(num_groups-1),T);
   }

   // Check for zero-sum of differences
   double sum = 0.0;
    
   for(int g = 0;g < num_groups ;g++)
   {
      sum += dEB(g);
   }

   // cout << "Sum EB = " << sum << "\n" << endl;
   
   // Print energy differences of EB
   // cout << left << setw(7)  << "Group"
   //      << left << setw(14) << "Energy" 
   //      << left << setw(14) << "dEB" << endl;
       
   // cout << left << setw(7)  << "-----"
   //      << left << setw(14) << "------" 
   //      << left << setw(14) << "-----" << endl;
       
   
   // for(int g = 0;g < num_groups ;g++)
   // {
   //    cout << left << setw(7)  << g 
   //         << left << setw(14) << e_ave_ref(g) 
   //         << left << setw(14) << dEB(g) << endl;
   // }
   
   // cout << "\n" << endl;

   // Compute energy derivatives
   dsigEdE(0) = kappa_edge(1)*e_edge_ref(1)/de_ave_ref(0);
   for(int g = 1;g < num_groups-1; g++)
   {
      dsigEdE(g) = (kappa_edge(g+1)*e_edge_ref(g+1) - kappa_edge(g)*e_edge_ref(g))/de_ave_ref(g);
   }
   dsigEdE(num_groups-1) = -kappa_edge(num_groups)*e_edge_ref(num_groups)/de_ave_ref(num_groups-1);

   // Check integral derivatives of kappa_ref E
   double sumabs=0.0;  
   for(int g = 0;g < num_groups-1; g++)
   {
    	sum += dsigEdE(g)*de_ave_ref(g);
    	sumabs += fabs(dsigEdE(g))*de_ave_ref(g);
   }
   sum=sum/sumabs;
   // cout << "Integral dsigEdE/Integral |dsigEdE| = " << sum << "\n" << endl;

   // Print energy derivatives of kappa_ref E
   // cout << left << setw(7)  << "Group"
   //      << left << setw(14) << "Energy" 
   //      << left << setw(14) << "dsigEdE" << endl;
       
   // cout << left << setw(7)  << "-----"
   //      << left << setw(14) << "------" 
   //      << left << setw(14) << "-------" << endl;
       
   // for(int g = 0;g < num_groups; g++)
   // {
   //    cout << left << setw(7)  << g 
   //         << left << setw(14) << e_ave_ref(g) 
   //         << left << setw(14) << dsigEdE(g) << endl;
   // }
   // cout << "\n" << endl;

   // Compute energy differences of kappa_ref*EB
   dkapEB(0) = kappa_edge(1)*e_edge_ref(1)*pf(e_edge_ref(1),T);
   if (num_groups > 1) // Only necessary for non-grey situations
   {
      for(int g = 1; g < num_groups-1; g++)
      {  
         dkapEB(g) = kappa_edge(g+1)*e_edge_ref(g+1)*pf(e_edge_ref(g+1),T) - kappa_edge(g)*e_edge_ref(g)*pf(e_edge_ref(g),T);
      }
      dkapEB(num_groups-1) = - kappa_edge(num_groups-1)*e_edge_ref(num_groups-1)*pf(e_edge_ref(num_groups-1),T);
   }
  
   // Check for zero sum of differences
   sum = 0.0;   
   for(int g = 0; g < num_groups; g++)
   {
    	sum += dkapEB(g);
    	sumabs += fabs(dkapEB(g));
   }
   sum=sum/sumabs;

   // cout << "Sum dkapEB/Sum |dkapEB| = " << sum << "\n" << endl; 

   // cout << left << setw(7)  << "Group"
   //      << left << setw(14) << "Energy" 
   //      << left << setw(14) << "dkapEB" << endl;
       
   // cout << left << setw(7) << "-----"
   //      << left << setw(14) << "------" 
   //      << left << setw(14) << "------" << endl;
        
   // for(int g = 0;g < num_groups; g++)
   // {
   //    cout << left << setw(7)  << g 
   //         << left << setw(14) << e_edge_ref(g) 
   //         << left << setw(14) << dkapEB(g) << endl;
   // }
   // cout << "\n" << endl;
}


template<int num_groups>
Correction<num_groups>::Correction(Eigen::Ref<Eigen::VectorXd> rho_vec,
                                   Eigen::Ref<Eigen::VectorXd> kappa_vec,
                                   Eigen::Ref<Eigen::VectorXd> T_vec,
                                   Eigen::Ref<Eigen::VectorXd> e_edge,
                                   Eigen::Ref<Eigen::VectorXd> e_ave,
                                   Eigen::Ref<Eigen::VectorXd> de_ave,
                                   Eigen::Ref<Eigen::MatrixXd> energy_discretization) : 
   rho_ref(rho_vec),
   kappa_ref(kappa_vec),
   T_ref(T_vec),
   e_edge_ref(e_edge),
   e_ave_ref(e_ave),
   de_ave_ref(de_ave),
   energy_discretization_ref(energy_discretization)
{
   // Must resize all Eigen class variables
   B.resize(num_groups);
   dBdT.resize(num_groups);

   ukappa.resize(num_groups);
   ckappa.resize(num_groups);
   emis_spec.resize(num_groups);

   dEB.resize(num_groups);
   dsigEdE.resize(num_groups);
   dkapEB.resize(num_groups);

   kappa_edge.resize(num_groups+1);
   cor1.resize(num_groups, ctv::N);
   cor2.resize(num_groups, ctv::N);
   cor3.resize(num_groups, ctv::N);
   total_correction.resize(ctv::M, num_groups, ctv::N);

   // All other helper functions depend on temperature, 
   // which will eventually change at each time step
}


template<int num_groups>
void Correction<num_groups>::compute_correction_terms()
{
   // Note in this function there is no dependency on scattering angle. This is accounted for when the total correction is computed.
   // Need to add spatial dependence
   for(int g = 0; g < num_groups; g++)
   {
      for (int cell_it = 0; cell_it < ctv::N; cell_it++)
      {
         cor1(g, cell_it) = dsigEdE(g);
         cor2(g, cell_it) = 3.0*rho_ref(g)*kappa_ref(g)*B(g) - dkapEB(g);
         cor3(g, cell_it) = cor1(g)*(4.0*B(g) - dEB(g));
      }
   }

   // cout << left << setw(7)  << "Group"
   //      << left << setw(14) << "Energy" 
   //      << left << setw(14) << "cor1" 
   //      << left << setw(14) << "cor2" 
   //      << left << setw(14) << "cor3" << endl;
 
   // cout << left << setw(7)  << "-----"
   //      << left << setw(14) << "------" 
   //      << left << setw(14) << "----" 
   //      << left << setw(14) << "----" 
   //      << left << setw(14) << "----" << endl;
 
   // for(int g = 0;g < num_groups; g++)
   // {
   //    cout << left << setw(7)  << g 
   //         << left << setw(14) << e_ave_ref(g) 
   //         << left << setw(14) << cor1(g) 
   //         << left << setw(14) << cor2(g) 
   //         << left << setw(14) << cor3(g) << endl;
   // }
   // cout << "\n" << endl;
}


template<int num_groups>
bool Correction<num_groups>::validate_correction()
{
   return (validate_planck_integrals() && validate_emission());
}


template<int num_groups>
void Correction<num_groups>::compute_correction(Eigen::Tensor<double, 3>& intensities)
{
   // Script to run all helper functions
   generate_planck_integrals();
   // generate_multigroup_opacities(); // Assumes constant opacity in g, not in use currently
   compute_group_edge_opacities();
   compute_components_of_correction_source();
   compute_correction_terms();

   // Finally put it all together
   double mu = 0., beta = ctv::V / Constants::SPEED_OF_LIGHT, val = 0.;
   // correction 
   for (int mu_it = 0; mu_it < ctv::M; mu_it++)
   {
      double mu = ctv::G_x[mu_it];
      for (int g = 0; g < num_groups; g++)
      {
         for (int cell_it = 0; cell_it < ctv::N; cell_it++)
         {
            val = 0.;
            
            val = (cor1(g, cell_it)*intensities(mu_it, g, cell_it) + cor2(g, cell_it)) * mu * beta;
            val -= cor3(g, cell_it) * pow(mu, 2) * pow(beta, 2);
            total_correction(mu_it, g, cell_it) = val;
         }
         
      }
   }
   
}

template class Correction<ctv::G>;

} // End namespace rt

