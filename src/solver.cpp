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
Solver<num_groups>::Solver(Eigen::Tensor<double, 4>& psi_mat,
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
            phi_ref(g, cell_it) += ctv::G_w[mu_it] * psi_mat_ref(mu_it, g, cell_it, 0);
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
            F_ref(g, cell_it) += ctv::G_x[mu_it] * ctv::G_w[mu_it] * psi_mat_ref(mu_it,g,cell_it,0);
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
         _src += rho_vec[cell_it] * kappa_vec[cell_it] * ctv::a * ctv::c * pow(temperature[cell_it],4) * ctv::dx;
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
   const double timestep, const double mu, double &local_bdry)
{
   if (mu < 0)
   {
      // cout << "kappa_vec: " << kappa_vec << endl;
      double const_A = 1. + ctv::c*timestep*rho_vec[cell] * kappa_vec[cell];
      double const_B = ctv::c*timestep*mu;
      double const_C = timestep*rho_vec[cell] * kappa_vec[cell] * ctv::a * pow(ctv::c,2) / (4 * M_PI);        

      double _temp_val = (const_A * ctv::dx - const_B) / 2.;
      // cout << "\t============= mu < 0\n";
      // cout << "temp_val: " << _temp_val << endl;

      // cout << "c: " << ctv::c << endl;
      // cout << "dt: " << timestep << endl;
      // cout << "cA: " << const_A << ", cB: " << const_B << ", cC: " << const_C << endl;

      _mat(0,0) = _temp_val; 
      _mat(0,1) = const_B / 2.;
      _mat(1,0) = - const_B / 2.;
      _mat(1,1) = _temp_val;

      _temp_val = const_C * ctv::dx * pow(temperature[cell], 4) / 2.;
      _rhs(0) = _temp_val + ctv::dx * ends(scatteredDirIt,groupIt,cell,0) / 2.; // CORTODO: Add correction
      _rhs(1) = _temp_val - (const_B * local_bdry) + ctv::dx * ends(scatteredDirIt,groupIt,cell,1) / 2.; // CORTODO: Add correction

      // cout << "_rhs: " << _rhs << endl;
      // Invert matrix
      // inverseMatrix(_mat, 2, _mat_inverse);
      _mat_inverse = _mat.inverse();
      // cout << "_mat_inverse: " << _mat_inverse << endl;

      // Solve 
      // matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);
      _res = _mat_inverse * _rhs;
      // cout << "_res: " << _res << endl;

      // put the average of val and local boundary here
      // TODO: Fix hard coded time param
      psi_mat_ref(scatteredDirIt,groupIt,cell,0) = 0.5*(_res[0] + _res[1]);
      // cout << "psi_mat_ref(scatteredDirIt,groupIt,cell,0): " << psi_mat_ref(scatteredDirIt,groupIt,cell,0) << endl;

      ends(scatteredDirIt,groupIt,cell,0) = _res[0];
      ends(scatteredDirIt,groupIt,cell,1) = _res[1];

      // Set local boundary for next cell iteration
      local_bdry = _res[0];
      // cout << "local_bdry: " << local_bdry << endl;
   }
   else if (mu > 0)
   {
      // cout << "--------------------- mu > 0 BE call -----------------------\n";
      double const_A = 1. + ctv::c*timestep*rho_vec[cell] * kappa_vec[cell];
      double const_B = ctv::c*timestep*mu;
      double const_C = ctv::c*timestep*rho_vec[cell] * kappa_vec[cell] * ctv::a * ctv::c / (4 * M_PI);

      // cout << "cA: " << const_A << ", cB: " << const_B << ", cC: " << const_C << endl;

      double _temp_val = (const_A * ctv::dx + const_B) / 2.;
      _mat(0,0) = _temp_val; 
      _mat(0,1) = const_B / 2.;
      _mat(1,0) = - const_B / 2.;
      _mat(1,1) = _temp_val;
      // cout << "_mat: " << endl;
      // cout << _mat << endl;

      _temp_val = const_C * ctv::dx * pow(temperature[cell], 4) / 2.;
      _rhs[0] = _temp_val + (const_B * local_bdry) + ctv::dx * ends(scatteredDirIt,groupIt,cell,0) / 2.; // CORTODO: Add correction
      _rhs[1] = _temp_val + ctv::dx * ends(scatteredDirIt,groupIt,cell,1) / 2.; // CORTODO: Add correction

      // Invert matrix
      // inverseMatrix(_mat, 2, _mat_inverse);
      _mat_inverse = _mat.inverse();

      // Solve 
      // matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);
      _res = _mat_inverse * _rhs;
      // cout << "_res: " << _res << endl;

      // put the average of val and local boundary here
      // TODO: Fix hardcoded time param
      psi_mat_ref(scatteredDirIt,groupIt,cell,0) = 0.5*(_res[0] + _res[1]);
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
void Solver<num_groups>::solve()
{
   for (int _it = 0; _it < ctv::_max_timesteps; _it++)
   {
      cout << "============= Timestep: " << _it << " =============" << endl;
      prev_ends = ends;

      // TODO: Change hard coded parameter in time space
      correction->compute_correction(psi_mat_ref);
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
                     int m_neg = (ctv::M / 2) - 1 - diff; // Move in opposite direction to yield index of -mu [-ctv::c, -b, -ctv::a, ctv::a, b, ctv::c]

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
                        backwardEuler(cell_j, i, g, ctv::dt, mu, local_bdry);
                        break;
                     }
                     case 2: // BDF2
                     {
                        // CN Step
                        double const_A = ctv::c * ctv::dt * mu / 4.;
                        double const_B = 1. + (ctv::c * ctv::dt * rho_vec[cell_j] * kappa_vec[cell_j]) / 4.;
                        double const_C = 1. - (ctv::c * ctv::dt * rho_vec[cell_j] * kappa_vec[cell_j]) / 4.;
                        double const_D = ctv::c * ctv::dt * ctv::a * ctv::c / (8 * M_PI);

                        double _temp_val = (const_B * ctv::dx - const_A) / 2.;
                        _mat(0,0) = _temp_val; 
                        _mat(0,1) = const_A / 2.;
                        _mat(1,0) = - const_A / 2.;
                        _mat(1,1) = _temp_val;

                        _temp_val = const_D * ctv::dx * rho_vec[cell_j] * kappa_vec[cell_j] * pow(temperature[cell_j], 4) / 2.;

                        _rhs[0] = _temp_val + ((const_C * ctv::dx + const_A) / 2.) * ends(i,g,cell_j,0) - (const_A / 2.) * ends(i,g,cell_j,1);
                        _rhs[1] = _temp_val + ((const_C * ctv::dx + const_A) / 2.) * ends(i,g,cell_j,1) + (const_A / 2.) * ends(i,g,cell_j,0) - const_A * (half_local_bdry + local_bdry_prev_it);

                        // Solve CN
                        // inverseMatrix(_mat, 2, _mat_inverse);
                        _mat_inverse = _mat.inverse();
                        // matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);
                        _res = _mat_inverse * _rhs;

                        half_ends(i,g,cell_j,0) = _res[0];
                        half_ends(i,g,cell_j,1) = _res[1];

                        // BDF Step
                        const_A = 1. + (ctv::c * ctv::dt * rho_vec[cell_j] * kappa_vec[cell_j] / 12.);
                        const_B = mu * ctv::c * ctv::dt / 12.;
                        const_C = 1. -  (ctv::c * ctv::dt * rho_vec[cell_j] * kappa_vec[cell_j] / 3.);
                        const_D = ctv::c * ctv::dt * rho_vec[cell_j] * kappa_vec[cell_j] / 12.;
                        double const_E = ctv::c * ctv::dt * ctv::a * ctv::c / (8. * M_PI);

                        _temp_val = (const_A * ctv::dx - const_B) / 2.;
                        _mat(0,0) = _temp_val; 
                        _mat(0,1) = const_B / 2.;
                        _mat(1,0) = - const_B / 2.;
                        _mat(1,1) = _temp_val;

                        _temp_val = (const_E * ctv::dx * rho_vec[cell_j] * kappa_vec[cell_j] / 2.) * pow(temperature[cell_j], 4);
                        _rhs[0] = _temp_val + ((const_C * ctv::dx + 4. * const_B) / 2.) * half_ends(i,g,cell_j,0) - 2. * const_B * half_ends(i,g,cell_j,1);
                        _rhs[0] += ((const_B - const_D * ctv::dx) / 2.) * ends(i,g,cell_j,0) - (const_B / 2.) * ends(i,g,cell_j,1);
                        _rhs[1] = _temp_val + ((const_C * ctv::dx + 4. * const_B) / 2.) * half_ends(i,g,cell_j,1) + 2. * const_B * half_ends(i,g,cell_j,0);
                        _rhs[1] += ((const_B - const_D * ctv::dx) / 2.) * ends(i,g,cell_j,1) + (const_B / 2.) * ends(i,g,cell_j,0);
                        _rhs[1] -= const_B * (local_bdry + local_bdry_prev_it + 4. *  half_local_bdry);

                        // Invert matrix
                        // inverseMatrix(_mat, 2, _mat_inverse);
                        _mat_inverse = _mat.inverse();

                        // Solve 
                        // matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);
                        _res = _mat_inverse * _rhs;

                        // put the average of val and local boundary here
                        // TODO: Fix hard coded time param
                        psi_mat_ref(i,g,cell_j,0) = 0.5*(_res[0] + _res[1]);

                        ends(i,g,cell_j,0) = _res[0];
                        ends(i,g,cell_j,1) = _res[1];

                        // Set local boundary for next cell iteration
                        local_bdry = _res[0];
                        half_local_bdry = half_ends(i,g,cell_j,0);
                        local_bdry_prev_it = prev_ends(i,g,cell_j,0);
                        
                        break;
                     } // End BDF2 case
                     case 3: // CN
                     {
                        // CN Step
                        double const_A = ctv::c * ctv::dt * mu / 4.;
                        double const_B = 1. + (ctv::c * ctv::dt * rho_vec[cell_j] * kappa_vec[cell_j]) / 4.;
                        double const_C = 1. - (ctv::c * ctv::dt * rho_vec[cell_j] * kappa_vec[cell_j]) / 4.;
                        double const_D = ctv::c * ctv::dt * ctv::a * ctv::c / (8 * M_PI);

                        double _temp_val = (const_B * ctv::dx - const_A) / 2.;
                        _mat(0,0) = _temp_val; 
                        _mat(0,1) = const_A / 2.;
                        _mat(1,0) = - const_A / 2.;
                        _mat(1,1) = _temp_val;

                        _temp_val = const_D * ctv::dx * rho_vec[cell_j] * kappa_vec[cell_j] * pow(temperature[cell_j], 4) / 2.;

                        _rhs[0] = _temp_val + ((const_C * ctv::dx + const_A) / 2.) * ends(i,g,cell_j,0) - (const_A / 2.) * ends(i,g,cell_j,1);
                        _rhs[1] = _temp_val + ((const_C * ctv::dx + const_A) / 2.) * ends(i,g,cell_j,1) + (const_A / 2.) * ends(i,g,cell_j,0) - const_A * (half_local_bdry + local_bdry_prev_it);

                        // Invert matrix
                        // inverseMatrix(_mat, 2, _mat_inverse);
                        _mat_inverse = _mat.inverse();

                        // Solve 
                        // matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);
                        _res = _mat_inverse * _rhs;

                        // put the average of val and local boundary here
                        // TODO: Fix hard coded time param
                        psi_mat_ref(i,g,cell_j,0) = 0.5*(_res[0] + _res[1]);

                        ends(i,g,cell_j,0) = _res[0];
                        ends(i,g,cell_j,1) = _res[1];

                        // Set local boundary for next cell iteration
                        local_bdry = _res[0];

                        break;
                     } // End CN
                     case 4: // BDF2 Morel-Lou
                     {
                        if (half_step)
                        {
                           // CN Step
                           double const_A = ctv::c * ctv::dt * mu / 4.;
                           double const_B = 1. + (ctv::c * ctv::dt * rho_vec[cell_j] * kappa_vec[cell_j]) / 4.;
                           double const_C = 1. - (ctv::c * ctv::dt * rho_vec[cell_j] * kappa_vec[cell_j]) / 4.;
                           double const_D = ctv::c * ctv::dt * ctv::a * ctv::c / (8 * M_PI);

                           double _temp_val = (const_B * ctv::dx - const_A) / 2.;
                           _mat(0,0) = _temp_val; 
                           _mat(0,1) = const_A / 2.;
                           _mat(1,0) = - const_A / 2.;
                           _mat(1,1) = _temp_val;

                           _temp_val = const_D * ctv::dx * rho_vec[cell_j] * kappa_vec[cell_j] * pow(temperature[cell_j], 4) / 2.;

                           _rhs[0] = _temp_val + ((const_C * ctv::dx + const_A) / 2.) * ends(i,g,cell_j,0) - (const_A / 2.) * ends(i,g,cell_j,1);
                           _rhs[1] = _temp_val + ((const_C * ctv::dx + const_A) / 2.) * ends(i,g,cell_j,1) + (const_A / 2.) * ends(i,g,cell_j,0) - const_A * (half_local_bdry + local_bdry_prev_it);

                           // Solve CN
                           // inverseMatrix(_mat, 2, _mat_inverse);
                           _mat_inverse = _mat.inverse();
                           // matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);
                           _res = _mat_inverse * _rhs;

                           half_ends(i,g,cell_j,0) = _res[0];
                           half_ends(i,g,cell_j,1) = _res[1];

                           half_local_bdry = half_ends(i,g,cell_j,0);
                           local_bdry_prev_it = prev_ends(i,g,cell_j,0);
                        }
                        else
                        {
                           // BDF Step
                           double const_A = 1. + (ctv::c * ctv::dt * rho_vec[cell_j] * kappa_vec[cell_j] / 3.);
                           double const_B = mu * ctv::c * ctv::dt / 3.;
                           double const_C = 1. -  (ctv::c * ctv::dt * rho_vec[cell_j] * kappa_vec[cell_j] / 12.);
                           double const_D = ctv::c * ctv::dt * rho_vec[cell_j] * kappa_vec[cell_j] / 12.;
                           double const_E = ctv::c * ctv::dt * ctv::a * ctv::c / (8. * M_PI);

                           double _temp_val = (const_A * ctv::dx - const_B) / 2.;
                           _mat(0,0) = _temp_val; 
                           _mat(0,1) = const_B / 2.;
                           _mat(1,0) = - const_B / 2.;
                           _mat(1,1) = _temp_val;

                           _temp_val = (const_E * ctv::dx * rho_vec[cell_j] * kappa_vec[cell_j] / 2.) * pow(temperature[cell_j], 4);
                           _rhs[0] = _temp_val + ((const_C * ctv::dx + (const_B / 4.)) / 2.) * half_ends(i,g,cell_j,0) - (const_B / 8.) * half_ends(i,g,cell_j,1);
                           _rhs[0] += (((const_B / 4.) - const_D * ctv::dx) / 2.) * ends(i,g,cell_j,0) - (const_B / 8.) * ends(i,g,cell_j,1);
                           _rhs[1] = _temp_val + ((const_C * ctv::dx + (const_B / 4.)) / 2.) * half_ends(i,g,cell_j,1) + (const_B / 8.) * half_ends(i,g,cell_j,0);
                           _rhs[1] += (((const_B / 4.) - const_D * ctv::dx) / 2.) * ends(i,g,cell_j,1) + (const_B / 8.) * ends(i,g,cell_j,0);
                           _rhs[1] -= const_B * (local_bdry + (local_bdry_prev_it / 4.) + (half_local_bdry / 4.));

                           // Invert matrix
                           // inverseMatrix(_mat, 2, _mat_inverse);
                           _mat_inverse = _mat.inverse();

                           // Solve 
                           // matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);
                           _res = _mat_inverse * _rhs;

                           // put the average of val and local boundary here
                           // TODO: Fix hardcoded time param
                           psi_mat_ref(i,g,cell_j,0) = 0.5*(_res[0] + _res[1]);
                           // cout << "1: " << psi_mat_ref(i,g,cell_j,0) << endl;

                           ends(i,g,cell_j,0) = _res[0];
                           ends(i,g,cell_j,1) = _res[1];

                           // Set local boundary for next cell iteration
                           local_bdry = _res[0];
                           half_local_bdry = half_ends(i,g,cell_j,0);
                           local_bdry_prev_it = prev_ends(i,g,cell_j,0);
                        }
                        
                        break;
                     } // End BDF2 case
                     default:
                     {
                        cout << "Incorrect timestepping method provided.\n";
                        assert(false);
                     }
                  }
               }
               else 
               { // mu > 0
                  // Fill matrix and rhs according to timestep selected
                  switch(ctv::ts_method) {
                     case 1: // BE
                     {
                        backwardEuler(j, i, g, ctv::dt, mu, local_bdry);
                        break;
                     }
                     case 2: // BDF2
                     {
                        // CN Step
                        double const_A = ctv::c * ctv::dt * mu / 4.;
                        double const_B = 1. + (ctv::c * ctv::dt * rho_vec[j] * kappa_vec[j]) / 4.;
                        double const_C = 1. - (ctv::c * ctv::dt * rho_vec[j] * kappa_vec[j]) / 4.;
                        double const_D = ctv::c * ctv::dt * ctv::a * ctv::c / (8 * M_PI);

                        double _temp_val = (const_B * ctv::dx + const_A) / 2.;
                        _mat(0,0) = _temp_val; 
                        _mat(0,1) = const_A / 2.;
                        _mat(1,0) = - const_A / 2.;
                        _mat(1,1) = _temp_val;

                        _temp_val = const_D * ctv::dx * rho_vec[j] * kappa_vec[j] * pow(temperature[j], 4) / 2.;

                        _rhs[0] = _temp_val + ((const_C * ctv::dx - const_A) / 2.) * ends(i,g,j,0) - (const_A / 2.) * ends(i,g,j,1) + const_A * (half_local_bdry + local_bdry_prev_it);
                        _rhs[1] = _temp_val + ((const_C * ctv::dx - const_A) / 2.) * ends(i,g,j,1) + (const_A / 2.) * ends(i,g,j,0);

                        // Solve CN
                        // inverseMatrix(_mat, 2, _mat_inverse);
                        _mat_inverse = _mat.inverse();
                        // matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);
                        _res = _mat_inverse * _rhs;

                        half_ends(i,g,j,0) = _res[0];
                        half_ends(i,g,j,1) = _res[1];

                        // BDF Step
                        const_A = 1. + (ctv::c * ctv::dt * rho_vec[j] * kappa_vec[j] / 12.);
                        const_B = mu * ctv::c * ctv::dt / 12.;
                        const_C = 1. -  (ctv::c * ctv::dt * rho_vec[j] * kappa_vec[j] / 3.);
                        const_D = ctv::c * ctv::dt * rho_vec[j] * kappa_vec[j] / 12.;
                        double const_E = ctv::c * ctv::dt * ctv::a * ctv::c / (8. * M_PI);

                        _temp_val = (const_A * ctv::dx + const_B) / 2.;
                        _mat(0,0) = _temp_val; 
                        _mat(0,1) = const_B / 2.;
                        _mat(1,0) = - const_B / 2.;
                        _mat(1,1) = _temp_val;

                        _temp_val = (const_E * ctv::dx * rho_vec[j] * kappa_vec[j] / 2.) * pow(temperature[j], 4);
                        _rhs[0] = _temp_val + ((const_C * ctv::dx - 4. * const_B) / 2.) * half_ends(i,g,j,0) - 2. * const_B * half_ends(i,g,j,1);
                        _rhs[0] += -1. * ((const_B + const_D * ctv::dx) / 2.) * ends(i,g,j,0) - (const_B / 2.) * ends(i,g,j,1);
                        _rhs[0] += const_B * (local_bdry + local_bdry_prev_it + 4. * half_local_bdry);
                        _rhs[1] = _temp_val + ((const_C * ctv::dx - 4. * const_B) / 2.) * half_ends(i,g,j,1) + 2. * const_B * half_ends(i,g,j,0);
                        _rhs[1] += -1. * ((const_B + const_D * ctv::dx) / 2.) * ends(i,g,j,1) + (const_B / 2.) * ends(i,g,j,0);

                        // Solve BDF2
                        // inverseMatrix(_mat, 2, _mat_inverse);
                        _mat_inverse = _mat.inverse();

                        // Solve 
                        // matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);
                        _res = _mat_inverse * _rhs;

                        // put the average of val and local boundary here
                        // TODO: Fix hardcoded time param
                        psi_mat_ref(i,g,j,0) = 0.5*(_res[0] + _res[1]);

                        ends(i,g,j,0) = _res[0];
                        ends(i,g,j,1) = _res[1];

                        // Set local boundary for next cell iteration
                        local_bdry = _res[1];
                        half_local_bdry = half_ends(i,g,j,1);
                        local_bdry_prev_it = prev_ends(i,g,j,1);

                        break;
                     } // End BDF2 case
                     case 3: // CN
                     {
                        // CN Step
                        double const_A = ctv::c * ctv::dt * mu / 4.;
                        double const_B = 1. + (ctv::c * ctv::dt * rho_vec[j] * kappa_vec[j]) / 4.;
                        double const_C = 1. - (ctv::c * ctv::dt * rho_vec[j] * kappa_vec[j]) / 4.;
                        double const_D = ctv::c * ctv::dt * ctv::a * ctv::c / (8 * M_PI);

                        double _temp_val = (const_B * ctv::dx + const_A) / 2.;
                        _mat(0,0) = _temp_val; 
                        _mat(0,1) = const_A / 2.;
                        _mat(1,0) = - const_A / 2.;
                        _mat(1,1) = _temp_val;

                        _temp_val = const_D * ctv::dx * rho_vec[j] * kappa_vec[j] * pow(temperature[j], 4) / 2.;

                        _rhs[0] = _temp_val + ((const_C * ctv::dx - const_A) / 2.) * ends(i,g,j,0) - (const_A / 2.) * ends(i,g,j,1) + const_A * (half_local_bdry + local_bdry_prev_it);
                        _rhs[1] = _temp_val + ((const_C * ctv::dx - const_A) / 2.) * ends(i,g,j,1) + (const_A / 2.) * ends(i,g,j,0);

                        // Solve matrix
                        // inverseMatrix(_mat, 2, _mat_inverse);
                        _mat_inverse = _mat.inverse();

                        // Solve 
                        // matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);
                        _res = _mat_inverse * _rhs;

                        // put the average of val and local boundary here
                        // TODO: Fix hardcoded time param
                        psi_mat_ref(i,g,j,0) = 0.5*(_res[0] + _res[1]);

                        ends(i,g,j,0) = _res[0];
                        ends(i,g,j,1) = _res[1];

                        // Set local boundary for next cell iteration
                        local_bdry = _res[1];
                        
                        break;
                     } // End CN
                     case 4: // BDF2 Morel-Lou
                     {
                        if (half_step)
                        {
                           // CN Step
                           double const_A = ctv::c * ctv::dt * mu / 4.;
                           double const_B = 1. + (ctv::c * ctv::dt * rho_vec[j] * kappa_vec[j]) / 4.;
                           double const_C = 1. - (ctv::c * ctv::dt * rho_vec[j] * kappa_vec[j]) / 4.;
                           double const_D = ctv::c * ctv::dt * ctv::a * ctv::c / (8 * M_PI);

                           double _temp_val = (const_B * ctv::dx + const_A) / 2.;
                           _mat(0,0) = _temp_val; 
                           _mat(0,1) = const_A / 2.;
                           _mat(1,0) = - const_A / 2.;
                           _mat(1,1) = _temp_val;

                           _temp_val = const_D * ctv::dx * rho_vec[j] * kappa_vec[j] * pow(temperature[j], 4) / 2.;

                           _rhs[0] = _temp_val + ((const_C * ctv::dx - const_A) / 2.) * ends(i,g,j,0) - (const_A / 2.) * ends(i,g,j,1) + const_A * (half_local_bdry + local_bdry_prev_it);
                           _rhs[1] = _temp_val + ((const_C * ctv::dx - const_A) / 2.) * ends(i,g,j,1) + (const_A / 2.) * ends(i,g,j,0);

                           // Solve CN
                           // inverseMatrix(_mat, 2, _mat_inverse);
                           _mat_inverse = _mat.inverse();
                           // matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);
                           _res = _mat_inverse * _rhs;

                           half_ends(i,g,j,0) = _res[0];
                           half_ends(i,g,j,1) = _res[1];

                           half_local_bdry = half_ends(i,g,j,1);
                           local_bdry_prev_it = prev_ends(i,g,j,1);
                        }
                        else
                        {
                           // BDF Step
                           double const_A = 1. + (ctv::c * ctv::dt * rho_vec[j] * kappa_vec[j] / 3.);
                           double const_B = mu * ctv::c * ctv::dt / 3.;
                           double const_C = 1. -  (ctv::c * ctv::dt * rho_vec[j] * kappa_vec[j] / 12.);
                           double const_D = ctv::c * ctv::dt * rho_vec[j] * kappa_vec[j] / 12.;
                           double const_E = ctv::c * ctv::dt * ctv::a * ctv::c / (8. * M_PI);

                           double _temp_val = (const_A * ctv::dx + const_B) / 2.;
                           _mat(0,0) = _temp_val; 
                           _mat(0,1) = const_B / 2.;
                           _mat(1,0) = - const_B / 2.;
                           _mat(1,1) = _temp_val;

                           _temp_val = (const_E * ctv::dx * rho_vec[j] * kappa_vec[j] / 2.) * pow(temperature[j], 4);
                           _rhs[0] = _temp_val + ((const_C * ctv::dx - const_B/4.) / 2.) * half_ends(i,g,j,0) - (const_B / 8.) * half_ends(i,g,j,1);
                           _rhs[0] += -1. * (((const_B / 4.) + const_D * ctv::dx) / 2.) * ends(i,g,j,0) - (const_B / 8.) * ends(i,g,j,1);
                           _rhs[0] += const_B * (local_bdry + (local_bdry_prev_it / 4.) + (half_local_bdry / 4.));
                           _rhs[1] = _temp_val + ((const_C * ctv::dx - (const_B / 4.)) / 2.) * half_ends(i,g,j,1) + (const_B / 8.) * half_ends(i,g,j,0);
                           _rhs[1] += -1. * (((const_B / 4.) + const_D * ctv::dx) / 2.) * ends(i,g,j,1) + (const_B / 8.) * ends(i,g,j,0);

                           // Solve BDF2
                           // inverseMatrix(_mat, 2, _mat_inverse);
                           _mat_inverse = _mat.inverse();

                           // Solve 
                           // matrixVectorMultiply(_mat_inverse, 2, _rhs, _res);
                           _res = _mat_inverse * _rhs;

                           // put the average of val and local boundary here
                           // TODO: Fix hardcoded time param
                           psi_mat_ref(i,g,j,0) = 0.5*(_res[0] + _res[1]);

                           ends(i,g,j,0) = _res[0];
                           ends(i,g,j,1) = _res[1];

                           // Set local boundary for next cell iteration
                           local_bdry = _res[1];
                           half_local_bdry = half_ends(i,g,j,1);
                           local_bdry_prev_it = prev_ends(i,g,j,1);
                        }

                        break;
                     } // End BDF2 case
                     default:
                     {
                        cout << "Incorrect timestepping method provided.\n";
                     }
                  }
               }
            }
            if (mu < 0)
            {
               cout << "Final psi(" << g << ") for mu " << mu << " is: " << ends(i,g,0,0) << endl;
            }
            else
            {
               cout << "Final psi(" << g << ") for mu " << mu << " is: " << ends(i,g,ctv::N-1,1) << endl;
            }
            cout << "Corresponding source condition was: " <<  ctv::psi_source[i] << endl;
         } // End group iterator
         
      } // End scattered direction loop
      if (ctv::ts_method == 2 || ctv::ts_method == 4)
      {
         half_step = !half_step;
      }
   } // End time loop
}

template class Solver<ctv::G>;
}