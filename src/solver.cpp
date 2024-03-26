#include "solver.h"


namespace rt 
{
void Solver::generate_group_edges()
{
   logfac = (log(elast)-log(efirst))/(num_groups-1.0);
   logfac = exp(logfac);
   if (num_groups == 1) { assert(logfac = 1.); } // In case of grey eqs, logfac should be 1.

   e_edge(0) = 0.0;
   e_edge(1) = efirst;

   for(int g = 1; g < num_groups; g++)
   {
      e_edge(g+1) = e_edge(g)*logfac;
   }
}


void Solver::generate_group_averages()
{
   e_ave(0) = 0.5*(e_edge(0)+e_edge(1));
   de_ave(0) = e_edge(1) - e_edge(0);

   for(int g = 1; g < num_groups; g++)
   {
      e_ave(g) = 0.5*(e_edge(g) + e_edge(g+1));
      de_ave(g) = e_edge(g+1) - e_edge(g);
   }
}


void Solver::fill_energy_bound_arrays()
{
   // Fill energy bound arrays for Planck integrations
   for (int g = 0; g < num_groups; g++)
   {
      energy_discretization(g, 0) = e_edge(g);
      energy_discretization(g, 1) = e_edge(g+1);
   }
}


Solver::Solver(ParameterHandler & parameter_handler,
               Eigen::Tensor<double, 3>& psi_mat,
               Eigen::Ref<Eigen::MatrixXd> phi,
               Eigen::Ref<Eigen::MatrixXd> F) :
   ph(parameter_handler),
   psi_mat_ref(psi_mat),
   phi_ref(phi),
   F_ref(F)
{
   cout << "Solver constructor.\n";
   // Get vals from parameter handler
   M = ph.get_M();
   N = ph.get_N();
   num_groups = ph.get_G();

   dx = ph.get_dx();
   dt = ph.get_dt();

   efirst = ph.get_efirst();
   elast = ph.get_elast();

   psi_source.resize(M,num_groups);

   // if necessary, get source conditions
   if (ph.get_bc_left_indicator() == 1 || ph.get_bc_right_indicator() == 1)
   {
      ph.get_psi_source(psi_source);
   }

   // Set up quadrature 
   GLQuad quad(M, Constants::FOUR_PI);
   this->m_mu = quad.mu();
   this->m_wt = quad.wt();
   
   cout << setw(16) << left << "Mu" << setw(16) << left << "Wt" << endl
        << setw(16) << left << "--" << setw(16) << left << "--" << endl;

   for(int i = 0;i < M;++i)
   {
      cout << showpos << setw(16) << left << m_mu(i) << setw(16) << left << m_wt(i) << endl;
   }
   cout << noshowpos << endl;

   // Set up energy discretization
   e_edge.resize(num_groups+1);
   e_ave.resize(num_groups);
   de_ave.resize(num_groups);
   energy_discretization.resize(num_groups, 2);

   // Optionally set group edges from file
   if (ph.get_have_group_bounds()) {
      ph.get_group_bounds(e_edge);
   } else {
      generate_group_edges();
   }
   
   // In either case, generate group averages and fill energy bound arrays
   generate_group_averages();
   fill_energy_bound_arrays();

   // Print energy group information
   cout << left << setw(13) << "Group Index" 
        << left << setw(16) << "Average Energy" 
        << left << setw(14) << "Upper Energy" 
        << left << setw(13) << "Group Width" << endl;
   cout << left << setw(13) << "-----------"
        << left << setw(16) << "(keV)---------" 
        << left << setw(14) << "(keV)-------" 
    	  << left << setw(13) << "(keV)------" << endl;
   for(int g = 0;g < num_groups;++g)
   {
      cout << left << setw(13) << g
    	     << left << setw(16) << e_ave(g) 
           << left << setw(14) << e_edge(g+1) 
           << left << setw(13) << de_ave(g) << endl;
   }
   cout << "\n" << endl;	

   // Set up slab specifics 
   kappa_mat.resize(M,N);                  
   ends.resize(M,num_groups,N,2);
   left_ends.resize(num_groups);
   right_ends.resize(num_groups);
   prev_ends.resize(M,num_groups,N,2);
   half_ends.resize(M,num_groups,N,2); 
   phi_plus.resize(num_groups, N);
   rho_vec.resize(num_groups);
   kappa_vec.resize(num_groups);
   temperature.resize(num_groups);

   _mat.resize(2,2);
   _mat_inverse.resize(2,2);
   _rhs.resize(2);
   _res.resize(2);   

   // resize balance
   balance.resize(num_groups);

   // Fill kappa_vec depending on prm file
   if (ph.get_have_group_absorption_opacities()) {
      ph.get_group_kappa(kappa_vec);
   } else {
      // For now, set this to grey
      kappa_vec.setConstant(ph.get_kappa_grey());
      // TODO: Can generate these opacities as is done in 
      // Correction::generate_multigroup_opacities()
   }

   // Fill rest with constants for now
   // TODO: Will need to change this
   rho_vec.setConstant(ph.get_rho());
   temperature.setConstant(ph.get_T());
   dEB.resize(num_groups);
   B.resize(num_groups);

   // Create correction object   
   correction = new Correction(parameter_handler, rho_vec, kappa_vec, temperature, e_edge, 
                               e_ave, de_ave, energy_discretization, m_mu, m_wt);

   // Initialize solution to sigma_a,g B_g
   // dir, group, cell
   correction->get_B(B);
   for (int i = 0; i < M; i++)
   {
      double mu = m_mu[i];
      for (int g = 0; g < num_groups; g++)
      {
         double val = B(g) * rho_vec(g) * kappa_vec(g);
         for (int cell_it = 0; cell_it < N; cell_it++)
         {
            psi_mat_ref(i, g, cell_it) = val;
            ends(i,g,cell_it,0) = val;
            ends(i,g,cell_it,1) = val;
         }
      }
   }

   cout << "B: " << B << endl;
   cout << "psi_mat_ref: " << psi_mat_ref << endl;


   cout << "end solver constructor\n";
}


void Solver::compute_angle_integrated_intensity()
{
   for (int g = 0; g < num_groups; g++)
   {
      for (int cell_it = 0; cell_it < N; cell_it++)
      {
         phi_ref(g, cell_it) = 0.;
         for (int mu_it = 0; mu_it < M; mu_it++)
         {
            phi_ref(g, cell_it) += m_wt[mu_it] * psi_mat_ref(mu_it, g, cell_it);
         }
      }
   }
}


void Solver::compute_positive_angle_integrated_intensity()
{
   for (int g = 0; g < num_groups; g++)
   {
      for (int cell_it = 0; cell_it < N; cell_it++)
      {
         phi_plus(g, cell_it) = 0.;
         for (int mu_it = M/2; mu_it < M; mu_it++)
         {
            assert(m_mu[mu_it] > 0.);
            phi_plus(g, cell_it) += m_wt[mu_it] * psi_mat_ref(mu_it, g, cell_it);
         }
      }
   }
}


void Solver::compute_radiative_flux()
{
   for (int g = 0; g < num_groups; g++)
   {
      for (int cell_it = 0; cell_it < N; cell_it++)
      {
         F_ref(g, cell_it) = 0.;
         for (int mu_it = 0.; mu_it < M; mu_it++)
         {
            F_ref(g, cell_it) += m_mu[mu_it] * m_wt[mu_it] * psi_mat_ref(mu_it,g,cell_it);
         }
      }
   }
}


void Solver::compute_balance()
{
   double mu, j_half_minus, j_half_plus, jN_half_minus, 
          jN_half_plus, _abs, _src, _bal;

   for (int g = 0; g < num_groups; g++)
   {
      j_half_minus = 0.;
      j_half_plus = 0.;
      jN_half_minus = 0.;
      jN_half_plus = 0.;
      _abs = 0.;
      _src = 0.;
      for (int mu_it = 0; mu_it < M; mu_it++)
      {
         mu = m_mu[mu_it];

         if (mu < 0.)
         {
            j_half_minus -= ends(mu_it,g,0,0) * mu * m_wt[mu_it]; // psi_1/2
            jN_half_minus -= ends(mu_it,g,N-1,0) * mu * m_wt[mu_it]; // psi_N+1/2
         }
         else 
         {
            j_half_plus += ends(mu_it,g,0,1) * mu * m_wt[mu_it]; // psi_1/2
            jN_half_plus += ends(mu_it,g,N-1,1) * mu * m_wt[mu_it]; // psi_N+1/2
         }
      }

      for (int cell_it = 0; cell_it < N; cell_it++)
      {
         _abs += rho_vec[g] * kappa_vec[g] * phi_ref(g,cell_it) * dx;
         _src += rho_vec[g] * kappa_vec[g] * ac * pow(temperature[g],4) * dx;
      }

      double sources = j_half_plus + jN_half_minus + _src;
      double sinks = jN_half_plus + j_half_minus + _abs;
      double sources_planck = 1.;
      cout << "sources: " << sources << endl;
      cout << "sinks: " << sinks << endl;

      balance(g) = abs( sinks - sources ) / sources;
      cout << "balance at (" << g << "): " << balance(g) << endl;
   }
}


void Solver::computeEquilibriumSources()
{
   correction->compute_correction(psi_mat_ref);
   if (ph.get_validation())
   {
      assert(correction->validate_correction() && "Invalid Correction Terms\n");
   }
   
   correction->get_B(this->B);
   correction->get_dEB(this->dEB);
   for (int i = 0; i < M; i++)
   {
      for (int g = 0; g < num_groups; g++)
      {
         double val = 4 * B(g) - dEB(g); // dEB represents the derivative of EB multiply by energy group width
         // cout << "v/c terms(1): " << val << endl;
         double _mult = m_mu(i) * ph.get_V() / Constants::SPEED_OF_LIGHT;
         // cout << "m_mu: " << m_mu(i) << ", V: " << ph.get_V() << ", c: " << Constants::SPEED_OF_LIGHT << endl;
         // cout << "_mult: " << _mult << endl;
         val *=  _mult;
         // cout << "v/c terms(2): " << val << endl;
         val += B(g);
         // cout << "de_ave(g): " << de_ave(g) << endl;
         cout << "source condition for mu: " << m_mu(i) << " and group " << g << ": " << val << endl;
         // psi_source[i,g] = val;
         psi_source(i,g) = val;
      }
   }
}


/*** Time Stepping functions ***/
void Solver::backwardEuler(
   const int cell, const int scatteredDirIt, const int groupIt, 
   const double timestep, const double mu)
{
   // Constants are the same regardless of direction
   double const_A = 1. + Constants::SPEED_OF_LIGHT*timestep*rho_vec[groupIt] * kappa_vec[groupIt];
   double const_B = Constants::SPEED_OF_LIGHT*timestep*mu;

   if (mu < 0)
   {  
      double _temp_val = (const_A * dx - const_B) / 2.;

      _mat(0,0) = _temp_val; 
      _mat(0,1) = const_B / 2.;
      _mat(1,0) = - const_B / 2.;
      _mat(1,1) = _temp_val;

      // Planckian and correction terms
      _temp_val = 0.5 * Constants::SPEED_OF_LIGHT * timestep * dx * rho_vec[groupIt] * kappa_vec[groupIt] * this->B(groupIt);
      if (ph.get_use_correction())
      {
         _temp_val += 0.5 * Constants::SPEED_OF_LIGHT * timestep * dx * total_correction(scatteredDirIt, groupIt, cell);
      }

      // Fill RHS
      _rhs(0) = _temp_val + dx * ends(scatteredDirIt,groupIt,cell,0) / 2.;
      _rhs(1) = _temp_val - (const_B * local_bdry) + dx * ends(scatteredDirIt,groupIt,cell,1) / 2.;

      // Solve
      _mat_inverse = _mat.inverse();
      _res = _mat_inverse * _rhs;

      // put the average of val and local boundary here
      psi_mat_ref(scatteredDirIt,groupIt,cell) = 0.5*(_res[0] + _res[1]);
      // cout << "psi_mat_ref(scatteredDirIt,groupIt,cell): " << psi_mat_ref(scatteredDirIt,groupIt,cell) << endl;
      // assert(false);
      ends(scatteredDirIt,groupIt,cell,0) = _res[0];
      ends(scatteredDirIt,groupIt,cell,1) = _res[1];

      // Set local boundary for next cell iteration
      local_bdry = _res[0];
      // cout << "local_bdry: " << local_bdry << endl;
   }
   else if (mu > 0)
   {
      // cout << "--------------------- mu > 0 BE call -----------------------\n";
      double _temp_val = (const_A * dx + const_B) / 2.;
      _mat(0,0) = _temp_val; 
      _mat(0,1) = const_B / 2.;
      _mat(1,0) = - const_B / 2.;
      _mat(1,1) = _temp_val;
      // cout << "_mat: " << endl;
      // cout << _mat << endl;

      // Planckian and correction terms
      _temp_val = 0.5 * Constants::SPEED_OF_LIGHT * timestep * dx * rho_vec[groupIt] * kappa_vec[groupIt] * this->B(groupIt);
      if (ph.get_use_correction())
      {
         _temp_val += 0.5 * Constants::SPEED_OF_LIGHT * timestep * dx * total_correction(scatteredDirIt, groupIt, cell);
      }

      // Fill RHS
      _rhs[0] = _temp_val + (const_B * local_bdry) + dx * ends(scatteredDirIt,groupIt,cell,0) / 2.;
      _rhs[1] = _temp_val + dx * ends(scatteredDirIt,groupIt,cell,1) / 2.;

      // Solve
      _mat_inverse = _mat.inverse();
      _res = _mat_inverse * _rhs;

      // put the average of val and local boundary here
      psi_mat_ref(scatteredDirIt,groupIt,cell) = 0.5*(_res[0] + _res[1]);
      // cout << "psi_mat_ref(scatteredDirIt,groupIt,cell,0): " << psi_mat_ref(scatteredDirIt,groupIt,cell,0) << endl;

      ends(scatteredDirIt,groupIt,cell,0) = _res[0];
      ends(scatteredDirIt,groupIt,cell,1) = _res[1];
      // cout << "ends: " << ends(scatteredDirIt,groupIt,cell,0) << " and " << ends(scatteredDirIt,groupIt,cell,0) << endl;

      // Set local boundary for next cell iteration
      local_bdry = _res[1];
      // cout << "local_bdry: " << _res[1] << endl;
   }
   else 
   {
      assert(false && "mu should not be 0.\n");
   }
}


void Solver::crankNicolson(
   const int cell, const int scatteredDirIt, const int groupIt, 
   const double timestep, const double mu)
{
   double _temp_val = 0.5 * Constants::SPEED_OF_LIGHT * timestep * rho_vec[groupIt] * kappa_vec[groupIt];

   double const_A = 0.5 * Constants::SPEED_OF_LIGHT * mu * timestep;
   double const_B = 1 + _temp_val;
   double const_C = 1 - _temp_val;

   if (mu < 0)
   {
      // Fill matrix to invert
      _temp_val = 0.5 * (const_B * dx - const_A);

      _mat(0,0) = _temp_val; 
      _mat(0,1) = 0.5 * const_A;
      _mat(1,0) = - 0.5 * const_A;
      _mat(1,1) = _temp_val;

      // Planckian and correction terms
      _temp_val = 0.5 * Constants::SPEED_OF_LIGHT * timestep * dx * rho_vec[groupIt] * kappa_vec[groupIt] * this->B(groupIt);
      if (ph.get_use_correction())
      {
         _temp_val += 0.5 * Constants::SPEED_OF_LIGHT * timestep * dx * total_correction(scatteredDirIt, groupIt, cell);
      }

      // Fill RHS
      _rhs[0] = _temp_val + 0.5 * (const_C * dx + const_A) * ends(scatteredDirIt,groupIt,cell,0) - 0.5 * const_A * ends(scatteredDirIt,groupIt,cell,1);
      _rhs[1] = _temp_val + 0.5 * const_A * ends(scatteredDirIt,groupIt,cell,0) + 0.5 * (const_C * dx + const_A) * ends(scatteredDirIt,groupIt,cell,1) - const_A * (local_bdry_prev_it + half_local_bdry);

      // Solve
      _mat_inverse = _mat.inverse();
      _res = _mat_inverse * _rhs;

      // put the average of val and local boundary here
      psi_mat_ref(scatteredDirIt,groupIt,cell) = 0.5*(_res[0] + _res[1]);

      ends(scatteredDirIt,groupIt,cell,0) = _res[0];
      ends(scatteredDirIt,groupIt,cell,1) = _res[1];

      // Set local boundary for next cell iteration
      local_bdry_prev_it = prev_ends(scatteredDirIt,groupIt,cell,0);
      half_local_bdry = _res[0];
   }
   else if (mu > 0)
   {
      // Fill matrix to invert
      _temp_val = 0.5 * (const_A + const_B * dx);
      _mat(0,0) = _temp_val; 
      _mat(0,1) = const_A / 2.;
      _mat(1,0) = - const_A / 2.;
      _mat(1,1) = _temp_val;

      // Planckian and correction terms
      _temp_val = 0.5 * Constants::SPEED_OF_LIGHT * timestep * dx * rho_vec[groupIt] * kappa_vec[groupIt] * this->B(groupIt);
      if (ph.get_use_correction())
      {
         _temp_val += 0.5 * Constants::SPEED_OF_LIGHT * timestep * dx * total_correction(scatteredDirIt, groupIt, cell);
      }

      // Fill RHS
      _rhs[0] = _temp_val + 0.5 * (const_C * dx - const_A) * ends(scatteredDirIt,groupIt,cell,0) - 0.5 * const_A * ends(scatteredDirIt,groupIt,cell,1) + const_A * (local_bdry_prev_it + half_local_bdry);
      _rhs[1] = _temp_val + 0.5 * const_A * ends(scatteredDirIt,groupIt,cell,0) + 0.5 * (const_C * dx - const_A) * ends(scatteredDirIt,groupIt,cell,1);

      // Solve
      _mat_inverse = _mat.inverse();
      _res = _mat_inverse * _rhs;

      // put the average of val and local boundary here
      psi_mat_ref(scatteredDirIt,groupIt,cell) = 0.5*(_res[0] + _res[1]);

      ends(scatteredDirIt,groupIt,cell,0) = _res[0];
      ends(scatteredDirIt,groupIt,cell,1) = _res[1];

      // Set local boundary for next cell iteration
      local_bdry_prev_it = prev_ends(scatteredDirIt,groupIt,cell,1);
      half_local_bdry = _res[1];
   }
   else 
   {
      assert(false && "mu should not be 0.\n");
   }
}


void Solver::bdf(
   const int cell, const int scatteredDirIt, const int groupIt, 
   const double timestep, const double mu)
{
   // Constants are the same regardless of direction
   double _temp_val = Constants::SPEED_OF_LIGHT * rho_vec[groupIt] * kappa_vec[groupIt] * timestep / 6.;

   double const_A = 1. + _temp_val;
   double const_B = Constants::SPEED_OF_LIGHT * mu * dt / 6.;
   double const_C = 1. - 4. * _temp_val;
   double const_D = _temp_val;

   if (mu < 0)
   {
      // Fill matrix to invert
      _temp_val = 0.5 * (const_A * dx - const_B);
      _mat(0,0) = _temp_val; 
      _mat(0,1) = 0.5 * const_B;
      _mat(1,0) = - 0.5 * const_B;
      _mat(1,1) = _temp_val;

      // Planckian and correction terms
      _temp_val = 0.5 * Constants::SPEED_OF_LIGHT * timestep * dx * rho_vec[groupIt] * kappa_vec[groupIt] * this->B(groupIt);      
      if (ph.get_use_correction())
      {
         _temp_val += 0.5 * Constants::SPEED_OF_LIGHT * timestep * dx * total_correction(scatteredDirIt, groupIt, cell);
      }

      // Fill RHS
      _rhs[0] = _temp_val + 0.5 * (const_C * dx + 4. * const_B) * half_ends(scatteredDirIt,groupIt,cell,0) - 2. * const_B * half_ends(scatteredDirIt,groupIt,cell,1);
      _rhs[0] += 0.5 * (const_B - const_D * dx) * prev_ends(scatteredDirIt,groupIt,cell,0) - 0.5 * const_B * prev_ends(scatteredDirIt,groupIt,cell,1);

      _rhs[1] = _temp_val + 2. * const_B * half_ends(scatteredDirIt,groupIt,cell,0) + 0.5 * (const_C * dx + 4. * const_B) * half_ends(scatteredDirIt,groupIt,cell,1);
      _rhs[1] += 0.5 * const_B * prev_ends(scatteredDirIt,groupIt,cell,0) + 0.5 * (const_B - const_D * dx) * prev_ends(scatteredDirIt,groupIt,cell,1);
      _rhs[1] -= const_B * (local_bdry + 4. *  half_local_bdry + local_bdry_prev_it);

      // Solve
      _mat_inverse = _mat.inverse();
      _res = _mat_inverse * _rhs;

      // put the average of val and local boundary here
      psi_mat_ref(scatteredDirIt,groupIt,cell) = 0.5*(_res[0] + _res[1]);

      ends(scatteredDirIt,groupIt,cell,0) = _res[0];
      ends(scatteredDirIt,groupIt,cell,1) = _res[1];

      // Set local boundary for next cell iteration
      local_bdry = _res[0];
      half_local_bdry = half_ends(scatteredDirIt,groupIt,cell,0);
      local_bdry_prev_it = prev_ends(scatteredDirIt,groupIt,cell,0);
   }
   else if (mu > 0)
   {
      // Fill matrix to invert
      _temp_val = 0.5 * (const_A * dx + const_B);
      _mat(0,0) = _temp_val; 
      _mat(0,1) = 0.5 * const_B;
      _mat(1,0) = - 0.5 * const_B;
      _mat(1,1) = _temp_val;

      // Planckian and correction terms
      _temp_val = 0.5 * Constants::SPEED_OF_LIGHT * timestep * dx * rho_vec[groupIt] * kappa_vec[groupIt] * this->B(groupIt);
      if (ph.get_use_correction())
      {
         _temp_val += 0.5 * Constants::SPEED_OF_LIGHT * timestep * dx * total_correction(scatteredDirIt, groupIt, cell);
      }

      // Fill RHS
      _rhs[0] = _temp_val + 0.5 * (const_C * dx - 4. * const_B) * half_ends(scatteredDirIt,groupIt,cell,0) - 2. * const_B * half_ends(scatteredDirIt,groupIt,cell,1);
      _rhs[0] -= 0.5 * (const_B + const_D * dx) * prev_ends(scatteredDirIt,groupIt,cell,0) + 0.5 * const_B * prev_ends(scatteredDirIt,groupIt,cell,1);
      _rhs[0] += const_B * (local_bdry + 4. * half_local_bdry + local_bdry_prev_it);

      _rhs[1] = _temp_val + 2. * const_B * half_ends(scatteredDirIt,groupIt,cell,0) + 0.5 * (const_C * dx - 4. * const_B) * half_ends(scatteredDirIt,groupIt,cell,1);
      _rhs[1] += 0.5 * const_B * prev_ends(scatteredDirIt,groupIt,cell,0) -0.5 * (const_B + const_D * dx) * prev_ends(scatteredDirIt,groupIt,cell,1);

      // Solve 
      _mat_inverse = _mat.inverse();
      _res = _mat_inverse * _rhs;

      // put the average of val and local boundary here
      psi_mat_ref(scatteredDirIt,groupIt,cell) = 0.5*(_res[0] + _res[1]);

      ends(scatteredDirIt,groupIt,cell,0) = _res[0];
      ends(scatteredDirIt,groupIt,cell,1) = _res[1];

      // Set local boundary for next cell iteration
      local_bdry = _res[1];
      half_local_bdry = half_ends(scatteredDirIt,groupIt,cell,1);
      local_bdry_prev_it = prev_ends(scatteredDirIt,groupIt,cell,1);
   }
   else 
   {
      assert(false && "mu should not be 0.\n");
   }
}


void Solver::solve()
{
   // Need to change the time iteration to match the time stepping scheme
   int _max_timesteps = ph.get_max_timesteps();
   if (ph.get_ts_method() == 3) { 
      // Since a "full" timestep in BDF2 consists of
      // BE, CN, BE, BDF, we multiply the iterate by 4.
      _max_timesteps *= 4;
   }

   // Optionally, compute the equilibrium source terms
   if (ph.get_use_mg_equilib())
   {
      computeEquilibriumSources();
   }

   for (int _it = 0; _it < _max_timesteps; _it++)
   {
      correction->compute_correction(psi_mat_ref);
      if (ph.get_validation())
      {
         assert(correction->validate_correction() && "Invalid Correction Terms\n");
      }
      
      correction->get_B(this->B);
      if (ph.get_use_correction())
      {
         correction->get_correction(this->total_correction);
      }
      
      if (ph.get_ts_method() != 3 || _it % 4 == 0)
      {
         // Onles set prev_ends if not BDF2 timestepping, or if we've completed one full step
         cout << "============= Timestep: " << _it << " =============" << endl;
         prev_ends = ends;
      }

      // Iterate over scattered direction (value given by gaussian quadrature)
      for (int i = 0; i < M; i++)
      {
         mu = m_mu[i];

         // Iterate over groups
         for (int g = 0; g < num_groups; g++)
         {
            bdry_cond = 0.;

            /*** Initialize boundary conditions ***/
            // First determine bdry_cond value specific in ctv.
            if (mu < 0.)
            {
               switch(ph.get_bc_right_indicator()) {
                  case 0: // vacuum
                  {
                     bdry_cond = 0.;
                     break;
                  }
                  case 2: // reflective
                  {
                     bdry_cond = 0.;
                     // TODO: Implement check once finished sweep
                     break;
                  }
                  case 1: // source
                  {
                     bdry_cond = psi_source(i,g);
                     break;
                  }
                  default:
                  {
                     cout << "Incorrect boundary conditions provided.\n";
                     assert(false);
                  }
               }
            }
            else
            {
               switch(ph.get_bc_left_indicator()) {
                  case 0: // vacuum
                  {
                     bdry_cond = 0.;
                  }
                  case 1: // source
                  {
                     bdry_cond = psi_source(i,g);
                     break;
                  }
                  case 2: // reflective
                  {
                     // Get index j corresponding to direction -mu
                     int diff = i - (M / 2); // Difference from center index
                     int m_neg = (M / 2) - 1 - diff; // Move in opposite direction to yield index of -mu 

                     bdry_cond = ends(m_neg,g,0,0);
                     break;
                  }
                  default:
                  {
                     cout << "Incorrect boundary conditions provided.\n";
                     assert(false);
                  }
               }
            }

            // Values change in the cell iteration.
            local_bdry = bdry_cond;
            half_local_bdry = bdry_cond;
            local_bdry_prev_it = bdry_cond;

            // Iterate over cells in our mesh (The "sweep")
            for (int j = 0; j < N; j++)
            {
               if (mu < 0)
               {
                  int cell_j = N - j - 1; // Since we are moving right to left
                  // Fill matrix and rhs according to timestep selected
                  switch(ph.get_ts_method()) {
                     case 1: // BE
                     {
                        // TODO: This dt should be a full timestep.  Adjust the function and the call accordingly
                        // The funcall in the BDF2 method should be a half step.
                        // Will need to change constants in the function as well.
                        // This is a full BE step, so we pass in the full timestep dt
                        backwardEuler(cell_j, i, g, dt, mu);
                        break;
                     }
                     case 2: // CN
                     {
                        crankNicolson(cell_j, i, g, dt, mu);
                        break;
                     } // End CN
                     case 3: // BDF2
                     {
                        int ts_indicator = _it % 4;
                        switch (ts_indicator) {
                           case 0: // BE Predictor
                           {
                              backwardEuler(cell_j, i, g, dt/2., mu);
                              break;
                           }
                           case 1: // CN Corrector
                           {
                              crankNicolson(cell_j, i, g, dt/2., mu);
                              half_ends = ends;
                              break;
                           }
                           case 2: // 2nd BE Predictor
                           {
                              backwardEuler(cell_j, i, g, dt/2., mu);
                              break;
                           }
                           case 3: // BDF Corrector
                           {
                              bdf(cell_j, i, g, dt/2., mu);
                              break;
                           }
                           default:
                           {
                              assert(false && "Invalid ts_indicator value.\n");
                           }
                        }

                        break;
                     } // End BDF2 case
                     default:
                     {
                        assert(false && "Incorrect timestepping method provided.\n");
                     }
                  }
               }
               else 
               { // mu > 0
                  // Fill matrix and rhs according to timestep selected
                  switch(ph.get_ts_method()) {
                     case 1: // BE
                     {
                        backwardEuler(j, i, g, dt, mu);
                        break;
                     }
                     case 2: // CN
                     {
                        crankNicolson(j, i, g, dt, mu);
                        break;
                     } // End CN
                     case 3: // BDF2
                     {
                        int ts_indicator = _it % 4;
                        switch (ts_indicator) {
                           case 0: // BE Predictor
                           {
                              backwardEuler(j, i, g, dt/2., mu);
                              break;
                           }
                           case 1: // CN Corrector
                           {
                              crankNicolson(j, i, g, dt/2., mu);
                              break;
                           }
                           case 2: // 2nd BE Predictor
                           {
                              backwardEuler(j, i, g, dt/2., mu);
                              break;
                           }
                           case 3: // BDF Corrector
                           {
                              bdf(j, i, g, dt/2., mu);
                              break;
                           }
                           default:
                           {
                              assert(false && "Invalid ts_indicator value.\n");
                           }
                        }
                        break;
                     } // End BDF2 case
                     default:
                     {
                        assert(false && "Incorrect timestepping method provided.\n");
                     }
                  } // end ts_method
               } // end mu > 0
            } // end sweep
            // if (mu < 0) {
            //    cout << "Final psi(" << g << ") for mu " << mu << " is: " << ends(i,g,0,0) << endl;
            // } else {
            //    cout << "Final psi(" << g << ") for mu " << mu << " is: " << ends(i,g,N-1,1) << endl;
            // }
            // cout << "Corresponding source condition was: " << psi_source[i] << endl;
         } // End group iterator  
      } // End scattered direction loop
   } // End time loop

   // correction->Print();
}


void Solver::compute_group_ends()
{
   double mu = 0.;
   left_ends.setConstant(0.);
   right_ends.setConstant(0.);

   for (int g = 0; g < num_groups; g++)
   {
      for (int mu_it = 0; mu_it < M; mu_it++)
      {
         mu = m_mu[mu_it];

         if (mu < 0.)
         {
            left_ends(g) += ends(mu_it, g, 0, 0);
         }
         else
         {
            right_ends(g) += ends(mu_it, g, N-1, 1);
         }
      }
      left_ends(g) /= (de_ave(g) * Constants::SPEED_OF_LIGHT);
      right_ends(g) /= (de_ave(g) * Constants::SPEED_OF_LIGHT);
   }
}


void Solver::get_ends(const string side, Eigen::Ref<Eigen::VectorXd> group_ends)
{
   assert((side=="left" || side=="right") && "Invalid option for 'side'.");
   if (side == "left")
   {
      group_ends = left_ends;
   }
   else
   {
      group_ends = right_ends;
   }
}

} // End ns rt
