#include "solver.h"


namespace rt 
{
template<int num_groups>
void Solver<num_groups>::generate_group_edges_and_averages()
{
   logfac = (log(elast)-log(efirst))/(num_groups-1.0);
   logfac = exp(logfac);
   if (num_groups == 1) { assert(logfac = 1.); } // In case of grey eqs, logfac should be 1.

   e_edge(0) = 0.0;
   e_edge(1) = efirst;
   e_ave(0) = 0.5*(e_edge(0)+e_edge(1));
   de_ave(0) = e_edge(1) - e_edge(0);

   for(int g = 1; g < num_groups; g++)
   {
      e_edge(g+1) = e_edge(g)*logfac;
      e_ave(g) = 0.5*(e_edge(g) + e_edge(g+1));
      de_ave(g) = e_edge(g+1) - e_edge(g);
   }
}


template<int num_groups>
void Solver<num_groups>::fill_energy_bound_arrays()
{
   // Fill energy bound arrays for Planck integrations
   for (int g = 0; g < num_groups; g++)
   {
      energy_discretization(g, 0) = e_edge(g);
      energy_discretization(g, 1) = e_edge(g+1);
   }
}


template<int num_groups>
Solver<num_groups>::Solver(Eigen::Tensor<double, 3>& psi_mat,
                           Eigen::Ref<Eigen::MatrixXd> phi,
                           Eigen::Ref<Eigen::MatrixXd> F) :
   psi_mat_ref(psi_mat),
   phi_ref(phi),
   F_ref(F)
{
   e_edge.resize(num_groups+1);
   e_ave.resize(num_groups);
   de_ave.resize(num_groups);
   energy_discretization.resize(num_groups, 2);
                     
   kappa_mat.resize(ctv::M,ctv::N);                  
   ends.resize(ctv::M,num_groups,ctv::N,2);      
   prev_ends.resize(ctv::M,num_groups,ctv::N,2);
   half_ends.resize(ctv::M,num_groups,ctv::N,2); 
   rho_vec.resize(ctv::N);
   kappa_vec.resize(ctv::N);
   temperature.resize(ctv::N);

   _mat.resize(2,2);
   _mat_inverse.resize(2,2);
   _rhs.resize(2);
   _res.resize(2);   

   // Fill constants for now
   // TODO: Will need to change this
   rho_vec.setConstant(ctv::rho);
   kappa_vec.setConstant(ctv::kappa);
   temperature.setConstant(ctv::T);
   B.resize(num_groups);
                              
   correction = new Correction<num_groups>(rho_vec, kappa_vec, temperature, e_edge, 
                                           e_ave, de_ave, energy_discretization);

   generate_group_edges_and_averages();
   fill_energy_bound_arrays();
}


template<int num_groups>
void Solver<num_groups>::compute_angle_integrated_density()
{
   for (int g = 0; g < num_groups; g++)
   {
      for (int cell_it = 0; cell_it < ctv::N; cell_it++)
      {
         phi_ref(g, cell_it) = 0.;
         for (int mu_it = 0; mu_it < ctv::M; mu_it++)
         {
            // TODO: Change time variable hard coded here
            phi_ref(g, cell_it) += ctv::G_w[mu_it] * psi_mat_ref(mu_it, g, cell_it);
         }
      }
   }
}


template<int num_groups>
void Solver<num_groups>::compute_radiative_flux()
{
   for (int g = 0; g < num_groups; g++)
   {
      for (int cell_it = 0; cell_it < ctv::N; cell_it++)
      {
         F_ref(g, cell_it) = 0.;
         for (int mu_it = 0.; mu_it < ctv::M; mu_it++)
         {
            F_ref(g, cell_it) += ctv::G_x[mu_it] * ctv::G_w[mu_it] * psi_mat_ref(mu_it,g,cell_it);
         }
      }
   }
}


template<int num_groups>
void Solver<num_groups>::compute_balance()
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
      for (int mu_it = 0; mu_it < ctv::M; mu_it++)
      {
         mu = ctv::G_x[mu_it];

         if (mu < 0.)
         {
            j_half_minus -= ends(mu_it,g,0,0) * mu * ctv::G_w[mu_it]; // psi_1/2
            jN_half_minus -= ends(mu_it,g,ctv::N-1,1) * mu * ctv::G_w[mu_it]; // psi_N+1/2
         }
         else 
         {
            j_half_plus += ends(mu_it,g,0,0) * mu * ctv::G_w[mu_it]; // psi_1/2
            jN_half_plus += ends(mu_it,g,ctv::N-1,1) * mu * ctv::G_w[mu_it]; // psi_N+1/2
         }
      }

      for (int cell_it = 0; cell_it < ctv::N; cell_it++)
      {
         _abs += rho_vec[cell_it] * kappa_vec[cell_it] * phi_ref(g,cell_it) * ctv::dx;
         _src += rho_vec[cell_it] * kappa_vec[cell_it] * ac * pow(temperature[cell_it],4) * ctv::dx;
      }

      double sources = jN_half_plus + j_half_minus + _abs;
      double sinks = j_half_plus + jN_half_minus + _src;
      cout << "sources: " << sources << endl;
      cout << "sinks: " << sinks << endl;

      balance(g) = abs( (jN_half_plus + j_half_minus + _abs) - (j_half_plus + jN_half_minus + _src) ) / (j_half_plus + jN_half_minus + _src);
      cout << "balance at (" << g << "): " << balance(g) << endl;
   }
}


/*** Time Stepping functions ***/
template<int num_groups>
void Solver<num_groups>::backwardEuler(
   const int cell, const int scatteredDirIt, const int groupIt, 
   const double timestep, const double mu)
{
   // Constants are the same regardless of direction
   double const_A = 1. + Constants::SPEED_OF_LIGHT*timestep*rho_vec[cell] * kappa_vec[cell];
   double const_B = Constants::SPEED_OF_LIGHT*timestep*mu;
   // cout << "cA: " << const_A << ", cB: " << const_B << ", cC: " << const_C << endl;

   if (mu < 0)
   {
      // cout << "kappa_vec: " << kappa_vec << endl;    
      double _temp_val = (const_A * ctv::dx - const_B) / 2.;
      // cout << "\t============= mu < 0\n";
      // cout << "temp_val: " << _temp_val << endl;

      // cout << "c: " << Constants::SPEED_OF_LIGHT << endl;
      // cout << "dt: " << timestep << endl;

      _mat(0,0) = _temp_val; 
      _mat(0,1) = const_B / 2.;
      _mat(1,0) = - const_B / 2.;
      _mat(1,1) = _temp_val;

      // Planckian and correction terms
      _temp_val = 0.5 * Constants::SPEED_OF_LIGHT * timestep * ctv::dx * rho_vec[cell] * kappa_vec[cell] * this->B(groupIt);
      if (ctv::use_correction)
      {
         _temp_val += 0.5 * Constants::SPEED_OF_LIGHT * timestep * ctv::dx * total_correction(scatteredDirIt, groupIt, cell);
      }

      // Fill RHS
      _rhs(0) = _temp_val + ctv::dx * ends(scatteredDirIt,groupIt,cell,0) / 2.;
      _rhs(1) = _temp_val - (const_B * local_bdry) + ctv::dx * ends(scatteredDirIt,groupIt,cell,1) / 2.;

      // Solve
      _mat_inverse = _mat.inverse();
      _res = _mat_inverse * _rhs;

      // put the average of val and local boundary here
      // TODO: Fix hard coded time param
      psi_mat_ref(scatteredDirIt,groupIt,cell) = 0.5*(_res[0] + _res[1]);
      // cout << "psi_mat_ref(scatteredDirIt,groupIt,cell): " << psi_mat_ref(scatteredDirIt,groupIt,cell) << endl;

      ends(scatteredDirIt,groupIt,cell,0) = _res[0];
      ends(scatteredDirIt,groupIt,cell,1) = _res[1];

      // Set local boundary for next cell iteration
      local_bdry = _res[0];
      // cout << "local_bdry: " << local_bdry << endl;
   }
   else if (mu > 0)
   {
      // cout << "--------------------- mu > 0 BE call -----------------------\n";
      double _temp_val = (const_A * ctv::dx + const_B) / 2.;
      _mat(0,0) = _temp_val; 
      _mat(0,1) = const_B / 2.;
      _mat(1,0) = - const_B / 2.;
      _mat(1,1) = _temp_val;
      // cout << "_mat: " << endl;
      // cout << _mat << endl;

      // Planckian and correction terms
      _temp_val = 0.5 * Constants::SPEED_OF_LIGHT * timestep * ctv::dx * rho_vec[cell] * kappa_vec[cell] * this->B(groupIt);
      if (ctv::use_correction)
      {
         _temp_val += 0.5 * Constants::SPEED_OF_LIGHT * timestep * ctv::dx * total_correction(scatteredDirIt, groupIt, cell);
      }

      // Fill RHS
      _rhs[0] = _temp_val + (const_B * local_bdry) + ctv::dx * ends(scatteredDirIt,groupIt,cell,0) / 2.;
      _rhs[1] = _temp_val + ctv::dx * ends(scatteredDirIt,groupIt,cell,1) / 2.;

      // Solve
      _mat_inverse = _mat.inverse();
      _res = _mat_inverse * _rhs;

      // put the average of val and local boundary here
      // TODO: Fix hardcoded time param
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


template<int num_groups>
void Solver<num_groups>::crankNicolson(
   const int cell, const int scatteredDirIt, const int groupIt, 
   const double timestep, const double mu)
{
   double _temp_val = 0.5 * Constants::SPEED_OF_LIGHT * timestep * rho_vec[cell] * kappa_vec[cell];

   double const_A = 0.5 * Constants::SPEED_OF_LIGHT * mu * timestep;
   double const_B = 1 + _temp_val;
   double const_C = 1 - _temp_val;

   if (mu < 0)
   {
      // Fill matrix to invert
      _temp_val = 0.5 * (const_B * ctv::dx - const_A);

      _mat(0,0) = _temp_val; 
      _mat(0,1) = 0.5 * const_A;
      _mat(1,0) = - 0.5 * const_A;
      _mat(1,1) = _temp_val;

      // Planckian and correction terms
      _temp_val = 0.5 * Constants::SPEED_OF_LIGHT * timestep * ctv::dx * rho_vec[cell] * kappa_vec[cell] * this->B(groupIt);
      if (ctv::use_correction)
      {
         _temp_val += 0.5 * Constants::SPEED_OF_LIGHT * timestep * ctv::dx * total_correction(scatteredDirIt, groupIt, cell);
      }

      // Fill RHS
      _rhs[0] = _temp_val + 0.5 * (const_C * ctv::dx + const_A) * ends(scatteredDirIt,groupIt,cell,0) - 0.5 * const_A * ends(scatteredDirIt,groupIt,cell,1);
      _rhs[1] = _temp_val + 0.5 * const_A * ends(scatteredDirIt,groupIt,cell,0) + 0.5 * (const_C * ctv::dx + const_A) * ends(scatteredDirIt,groupIt,cell,1) - const_A * (local_bdry_prev_it + half_local_bdry);

      // Solve
      _mat_inverse = _mat.inverse();
      _res = _mat_inverse * _rhs;

      // put the average of val and local boundary here
      // TODO: Fix hard coded time param
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
      _temp_val = 0.5 * (const_A + const_B * ctv::dx);
      _mat(0,0) = _temp_val; 
      _mat(0,1) = const_A / 2.;
      _mat(1,0) = - const_A / 2.;
      _mat(1,1) = _temp_val;

      // Planckian and correction terms
      _temp_val = 0.5 * Constants::SPEED_OF_LIGHT * timestep * ctv::dx * rho_vec[cell] * kappa_vec[cell] * this->B(groupIt);
      if (ctv::use_correction)
      {
         _temp_val += 0.5 * Constants::SPEED_OF_LIGHT * timestep * ctv::dx * total_correction(scatteredDirIt, groupIt, cell);
      }

      // Fill RHS
      _rhs[0] = _temp_val + 0.5 * (const_C * ctv::dx - const_A) * ends(scatteredDirIt,groupIt,cell,0) - 0.5 * const_A * ends(scatteredDirIt,groupIt,cell,1) + const_A * (local_bdry_prev_it + half_local_bdry);
      _rhs[1] = _temp_val + 0.5 * const_A * ends(scatteredDirIt,groupIt,cell,0) + 0.5 * (const_C * ctv::dx - const_A) * ends(scatteredDirIt,groupIt,cell,1);

      // Solve
      _mat_inverse = _mat.inverse();
      _res = _mat_inverse * _rhs;

      // put the average of val and local boundary here
      // TODO: Fix hardcoded time param
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


template<int num_groups>
void Solver<num_groups>::bdf(
   const int cell, const int scatteredDirIt, const int groupIt, 
   const double timestep, const double mu)
{
   // Constants are the same regardless of direction
   double _temp_val = Constants::SPEED_OF_LIGHT * rho_vec[cell] * kappa_vec[cell] * timestep / 6.;

   double const_A = 1. + _temp_val;
   double const_B = Constants::SPEED_OF_LIGHT * mu * ctv::dt / 6.;
   double const_C = 1. - 4. * _temp_val;
   double const_D = _temp_val;

   if (mu < 0)
   {
      // Fill matrix to invert
      _temp_val = 0.5 * (const_A * ctv::dx - const_B);
      _mat(0,0) = _temp_val; 
      _mat(0,1) = 0.5 * const_B;
      _mat(1,0) = - 0.5 * const_B;
      _mat(1,1) = _temp_val;

      // Planckian and correction terms
      _temp_val = 0.5 * Constants::SPEED_OF_LIGHT * timestep * ctv::dx * rho_vec[cell] * kappa_vec[cell] * this->B(groupIt);      
      if (ctv::use_correction)
      {
         _temp_val += 0.5 * Constants::SPEED_OF_LIGHT * timestep * ctv::dx * total_correction(scatteredDirIt, groupIt, cell);
      }

      // Fill RHS
      _rhs[0] = _temp_val + 0.5 * (const_C * ctv::dx + 4. * const_B) * half_ends(scatteredDirIt,groupIt,cell,0) - 2. * const_B * half_ends(scatteredDirIt,groupIt,cell,1);
      _rhs[0] += 0.5 * (const_B - const_D * ctv::dx) * prev_ends(scatteredDirIt,groupIt,cell,0) - 0.5 * const_B * prev_ends(scatteredDirIt,groupIt,cell,1);

      _rhs[1] = _temp_val + 2. * const_B * half_ends(scatteredDirIt,groupIt,cell,0) + 0.5 * (const_C * ctv::dx + 4. * const_B) * half_ends(scatteredDirIt,groupIt,cell,1);
      _rhs[1] += 0.5 * const_B * prev_ends(scatteredDirIt,groupIt,cell,0) + 0.5 * (const_B - const_D * ctv::dx) * prev_ends(scatteredDirIt,groupIt,cell,1);
      _rhs[1] -= const_B * (local_bdry + 4. *  half_local_bdry + local_bdry_prev_it);

      // Solve
      _mat_inverse = _mat.inverse();
      _res = _mat_inverse * _rhs;

      // put the average of val and local boundary here
      // TODO: Fix hard coded time param
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
      _temp_val = 0.5 * (const_A * ctv::dx + const_B);
      _mat(0,0) = _temp_val; 
      _mat(0,1) = 0.5 * const_B;
      _mat(1,0) = - 0.5 * const_B;
      _mat(1,1) = _temp_val;

      // Planckian and correction terms
      _temp_val = 0.5 * Constants::SPEED_OF_LIGHT * timestep * ctv::dx * rho_vec[cell] * kappa_vec[cell] * this->B(groupIt);
      if (ctv::use_correction)
      {
         _temp_val += 0.5 * Constants::SPEED_OF_LIGHT * timestep * ctv::dx * total_correction(scatteredDirIt, groupIt, cell);
      }

      // Fill RHS
      _rhs[0] = _temp_val + 0.5 * (const_C * ctv::dx - 4. * const_B) * half_ends(scatteredDirIt,groupIt,cell,0) - 2. * const_B * half_ends(scatteredDirIt,groupIt,cell,1);
      _rhs[0] -= 0.5 * (const_B + const_D * ctv::dx) * prev_ends(scatteredDirIt,groupIt,cell,0) + 0.5 * const_B * prev_ends(scatteredDirIt,groupIt,cell,1);
      _rhs[0] += const_B * (local_bdry + 4. * half_local_bdry + local_bdry_prev_it);

      _rhs[1] = _temp_val + 2. * const_B * half_ends(scatteredDirIt,groupIt,cell,0) + 0.5 * (const_C * ctv::dx - 4. * const_B) * half_ends(scatteredDirIt,groupIt,cell,1);
      _rhs[1] += 0.5 * const_B * prev_ends(scatteredDirIt,groupIt,cell,0) -0.5 * (const_B + const_D * ctv::dx) * prev_ends(scatteredDirIt,groupIt,cell,1);

      // Solve 
      _mat_inverse = _mat.inverse();
      _res = _mat_inverse * _rhs;

      // put the average of val and local boundary here
      // TODO: Fix hardcoded time param
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



template<int num_groups>
void Solver<num_groups>::solve()
{
   // Need to change the time iteration to match the time stepping scheme
   int _max_timesteps = ctv::max_timesteps;
   if (ctv::ts_method == 3) { 
      // Since a "full" timestep in BDF2 consists of
      // BE, CN, BE, BDF, we multiply the iterate by 4.
      _max_timesteps *= 4;
   }

   for (int _it = 0; _it < _max_timesteps; _it++)
   {
      correction->compute_correction(psi_mat_ref);
      correction->get_B(this->B);
      if (ctv::use_correction)
      {
         correction->get_correction(this->total_correction);
         // CORTODO: If we don't want to use correction, we still need to compute B.
      }
      
      if (ctv::ts_method != 3 || _it % 4 == 0)
      {
         // Onles set prev_ends if not BDF2 timestepping, or if we've completed one full step
         cout << "============= Timestep: " << _it << " =============" << endl;
         prev_ends = ends;
      }

      // Iterate over scattered direction (value given by gaussian quadrature)
      for (int i = 0; i < ctv::M; i++)
      {
         mu = ctv::G_x[i];

         // Iterate over groups
         for (int g = 0; g < num_groups; g++)
         {
            bdry_cond = 0.;

            /*** Initialize boundary conditions ***/
            // First determine bdry_cond value specific in ctv.
            if (mu < 0.)
            {
               switch(ctv::bc_right_indicator) {
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
                     bdry_cond = ctv::psi_source[i];
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
               switch(ctv::bc_left_indicator) {
                  case 0: // vacuum
                  {
                     bdry_cond = 0.;
                  }
                  case 1: // source
                  {
                     bdry_cond = ctv::psi_source[i];
                     break;
                  }
                  case 2: // reflective
                  {
                     // Get index j corresponding to direction -mu
                     int diff = i - (ctv::M / 2); // Difference from center index
                     int m_neg = (ctv::M / 2) - 1 - diff; // Move in opposite direction to yield index of -mu 

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
            for (int j = 0; j < ctv::N; j++)
            {
               if (mu < 0)
               {
                  int cell_j = ctv::N - j - 1; // Since we are moving right to left
                  // Fill matrix and rhs according to timestep selected
                  switch(ctv::ts_method) {
                     case 1: // BE
                     {
                        // TODO: This dt should be a full timestep.  Adjust the function and the call accordingly
                        // The funcall in the BDF2 method should be a half step.
                        // Will need to change constants in the function as well.
                        // This is a full BE step, so we pass in the full timestep dt
                        backwardEuler(cell_j, i, g, ctv::dt, mu);
                        break;
                     }
                     case 2: // CN
                     {
                        crankNicolson(cell_j, i, g, ctv::dt, mu);
                        break;
                     } // End CN
                     case 3: // BDF2
                     {
                        int ts_indicator = _it % 4;
                        switch (ts_indicator) {
                           case 0: // BE Predictor
                           {
                              backwardEuler(cell_j, i, g, ctv::dt/2., mu);
                              break;
                           }
                           case 1: // CN Corrector
                           {
                              crankNicolson(cell_j, i, g, ctv::dt/2., mu);
                              half_ends = ends;
                              break;
                           }
                           case 2: // 2nd BE Predictor
                           {
                              backwardEuler(cell_j, i, g, ctv::dt/2., mu);
                              break;
                           }
                           case 3: // BDF Corrector
                           {
                              bdf(cell_j, i, g, ctv::dt/2., mu);
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
                  switch(ctv::ts_method) {
                     case 1: // BE
                     {
                        backwardEuler(j, i, g, ctv::dt, mu);
                        break;
                     }
                     case 2: // CN
                     {
                        crankNicolson(j, i, g, ctv::dt, mu);
                        break;
                     } // End CN
                     case 3: // BDF2
                     {
                        int ts_indicator = _it % 4;
                        switch (ts_indicator) {
                           case 0: // BE Predictor
                           {
                              backwardEuler(j, i, g, ctv::dt/2., mu);
                              break;
                           }
                           case 1: // CN Corrector
                           {
                              crankNicolson(j, i, g, ctv::dt/2., mu);
                              break;
                           }
                           case 2: // 2nd BE Predictor
                           {
                              backwardEuler(j, i, g, ctv::dt/2., mu);
                              break;
                           }
                           case 3: // BDF Corrector
                           {
                              bdf(j, i, g, ctv::dt/2., mu);
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
            if (mu < 0) {
               cout << "Final psi(" << g << ") for mu " << mu << " is: " << ends(i,g,0,0) << endl;
            } else {
               cout << "Final psi(" << g << ") for mu " << mu << " is: " << ends(i,g,ctv::N-1,1) << endl;
            }
            cout << "Corresponding source condition was: " <<  ctv::psi_source[i] << endl;
         } // End group iterator  
      } // End scattered direction loop
   } // End time loop
}

template class Solver<ctv::G>;
}