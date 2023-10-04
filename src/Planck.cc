// Copyright (c) 2000-2017, Texas Engineering Experiment Station (TEES), a
// component of the Texas A&M University System.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are not permitted without specific prior written permission
// from TEES.

// If written permission is obtained for redistribution or further use, the
// following conditions must be met:

// 1) Redistributions of source code must retain the above copyright notice,
// this list of conditions and the disclaimer below.

// 2) Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions, and the disclaimer below in the documentation and/or
// other materials provided with the distribution.

// 3) Neither the name of TEES, the name of the Texas A&M University System, nor
// the names of its contributors may be used to endorse or promote products
// derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "Planck.h"

#include <algorithm>

using std::pow;
using std::abs;
using std::max;
using std::copy;

Planck::Planck(double accuracy)
  : m_accuracy(accuracy)
{
  gauss_quad_setup();
}

void Planck::get_Planck(double T, Eigen::Ref<Eigen::MatrixXd> edisc,
	                     Eigen::Ref<Eigen::VectorXd> B, 
                        Eigen::Ref<Eigen::VectorXd> dBdT)
{
  assert(T >= 0.0);

  size_t num_groups = edisc.rows();
  assert(B.size() == num_groups);
  assert(dBdT.size() == num_groups);

  double B_sum = integrate_B_grey(T);
  double dBdT_sum = integrate_dBdT_grey(T);
  for(size_t g = 0;g < num_groups - 1;++g)
  {
    double integral = integrate_B(T, edisc(g,0), edisc(g,1));
    B[g]   = integral;
    B_sum -= integral;

    integral = integrate_dBdT(T, edisc(g,0), edisc(g,1));
    dBdT[g]   = integral;
    dBdT_sum -= integral;
  }

  if(B_sum > 0.0)
    B[num_groups-1] = B_sum;
  if(dBdT_sum > 0.0)
    dBdT[num_groups-1] = dBdT_sum;
}

double Planck::integrate_B_grey(double T)
{
  assert(T >= 0.0);
  return Constants::RADIATION_CONSTANT_A_LONG*Constants::SPEED_OF_LIGHT*pow(T, 4.0);
}

double Planck::integrate_B(double T, double E_min, double E_max)
{
  assert(T >= 0.0);
  assert(E_min >= 0.0);
  assert(E_max > E_min);

  if(equal(T, 0.0) || equal(E_min, E_max))
    return 0.0;

  auto integrate = [this](double z1, double z2)
  {
    // Determine number of terms
    int N = 32;
    double sum1 = exp(-z1)*(z1*z1*z1 + 3.0*z1*z1 + 6.0*z1 + 6.0);
    sum1 = max(sum1, std::numeric_limits<double>::epsilon());
    while(true)
    {
      double val = exp(-(N + 1.0)*z1)/(1.0 - exp(-z1))*pow(N + 1.0, -4.0)*(pow((N + 1.0)*z1, 3.0) +
        3.0*pow((N + 1.0)*z1, 2.0) + 6.0*(N + 1.0)*z1 + 6.0)/sum1;
      if(val  >  m_accuracy)
        ++N;
      else
        break;
    }
    // Compute integral from E_min to infinity
    sum1 = 0.0;
    double sum2 = 0.0;
    for(int n = N;n != 0;--n)
    {
      sum1 += exp(-n*z1)/pow(n, 4.0) * (pow(n*z1, 3.0)+3.0*pow(n*z1, 2.0) + 6.0*n*z1 + 6.0);
      sum2 += exp(-n*z2)/pow(n, 4.0) * (pow(n*z2, 3.0)+3.0*pow(n*z2, 2.0) + 6.0*n*z2 + 6.0);
    }
    return sum1 - sum2;
  };

  const double h = Constants::PLANCK_CONSTANT;
  const double k = Constants::BOLTZMANN_CONSTANT;
  const double c = Constants::SPEED_OF_LIGHT;

  // Unitless representation of energy and temperature
  double z1 = E_min/(k*T);
  double z2 = E_max/(k*T);

  // Bg is the group averaged planck function
  double Bg;

  if(z2 <= 0.7) // Integrate by gaussian quadrature
  {
    double g_mid = 0.5*(E_max + E_min);
    double g_map = 0.5*(E_max - E_min);
    double gauss = 0.0;
    for(size_t i = 0;i < m_points.size();++i)
      gauss += g_map * m_weights[i] * get_B(T, g_mid + g_map*m_points[i]);
    Bg = gauss;
  }
  else if (z1 >= 0.5) // Integrate by infinite sum
    Bg = 2.0*pow(k*T, 4.0)*integrate(z1, z2)/(pow(h, 3.0)*pow(c, 2.0));
  else // Split interval (split point is 0.6)
  {
    z1 = 0.6;
    double gauss = 0;
    double g_mid = 0.5 * (z1*k*T + E_min);
    double g_map = 0.5 * (z1*k*T - E_min);
    for(size_t i = 0;i < m_points.size();++i)
      gauss += g_map * m_weights[i] * get_B(T, g_mid + g_map * m_points[i]);
    Bg = gauss + 2.0*pow(k*T, 4.0)*integrate(z1, z2)/(pow(h, 3.0)*pow(c, 2.0));
  }

  return Bg*4.0*Constants::PI; // Note that this is NOT a per steradian quantity
}

double Planck::integrate_dBdT_grey(double T)
{
  return 4.0*Constants::RADIATION_CONSTANT_A_LONG*Constants::SPEED_OF_LIGHT*pow(T, 3.0);
}

double Planck::integrate_dBdT(double T, double E_min, double E_max)
{
  assert(T >= 0.0);
  assert(E_min >= 0.0);
  assert(E_max > E_min);

  if(equal(T, 0.0) || equal(E_min, E_max))
    return 0.0;

  auto integrate = [this](double z1, double z2)
  {
    int N = 32;
    double sum1 = exp(-z1)*(pow(z1, 4.0) + 4.0*pow(z1, 3.0) + 12.0*z1*z1 + 24.0*z1 + 24.0);
    sum1 = max(sum1, std::numeric_limits<double>::epsilon());
    while(true)
    {
      double val = exp(-(N + 1.0)*z1)/(1.0 - exp(-z1))*pow(N + 1.0, -4.0) * (pow((N + 1.0)*z1, 4.0) +
        4.0*pow((N + 1.0)*z1, 3.0) + 12.0*pow((N + 1.0)*z1, 2.0) + 24.0*(N + 1.0)*z1 + 24.0)/sum1;
      if(val > m_accuracy)
        ++N;
      else
        break;
    }

    sum1 = 0.0;
    double sum2 = 0.0;
    for(int n = N;n > 0;--n)
    {
      sum1 += exp(-n*z1)/pow(n, 4.0)*(pow(n*z1, 4.0) + 4.0*pow(n*z1, 3.0) + 12.0*pow(n*z1, 2.0) + 24.0*n*z1 + 24.0);
      sum2 += exp(-n*z2)/pow(n, 4.0)*(pow(n*z2, 4.0) + 4.0*pow(n*z2, 3.0) + 12.0*pow(n*z2, 2.0) + 24.0*n*z2 + 24.0);
    }
    return sum1 - sum2;
  };

  const double h = Constants::PLANCK_CONSTANT;
  const double k = Constants::BOLTZMANN_CONSTANT;
  const double c = Constants::SPEED_OF_LIGHT;

  // unitless representation of energy and temperature
  double z1 = E_min/(k*T);
  double z2 = E_max/(k*T);

  // dBgdT is the integral-dE of the partial-dT
  double dBgdT;

  if(z2 <= 0.7) // Integrate by gaussian quadrature
  {
    double gauss = 0.0;
    double g_mid = 0.5 * (E_max + E_min);
    double g_map = 0.5 * (E_max - E_min);
    for(size_t i = 0;i < m_points.size();++i)
      gauss += g_map * m_weights[i] * get_dBdT(T, g_mid + g_map * m_points[i]);
    dBgdT = gauss;
  }
  else if(z1 >= 0.5) // Integrate by infinite sum
    dBgdT = 2.0*pow(k, 4.0)*pow(T, 3.0)*integrate(z1, z2)/(pow(h, 3.0)*pow(c, 2.0));
  else // Split interval (split point is 0.6)
  {
    z1 = 0.6;
    double gauss = 0.0;
    double g_mid = 0.5*(z1*k*T + E_min);
    double g_map = 0.5*(z1*k*T - E_min);
    for(size_t i = 0;i < m_points.size();++i)
      gauss += g_map * m_weights[i] * get_dBdT(T, g_mid + g_map * m_points[i]);
    dBgdT = gauss + 2.0*pow(k, 4.0)*pow(T, 3.0)*integrate(z1, z2)/(pow(h, 3.0)*pow(c, 2.0));
  }

  return dBgdT*4.0*Constants::PI; // Note that this is NOT a per steradian quantity
}

void Planck::gauss_quad_setup()
{
  // This function sets the points and weights for integrating the planck
  // function from low_bound to high_bound byaussian quadrature.
  // n specifies the order of the quadrature.
  //
  // ====================================================
  // description of method
  // ====================================================
  //
  // A gaussian quadrature of order n is given by n points and n associated
  // weights.  The points are a direct mapping of the roots of the nth order
  // legendre polynomial, and the weights, which sum to the interval of
  // integration, are a function of the root values and the first derivative
  // of the legendre polynomial
  //
  // ====================================================
  // declaration of variables
  // ====================================================
  //
  // a quadrature order of 12 is sufficient to integrate B for z < 1
  // for more on this, see B_poly_z.cc
  const short unsigned int order = 12;

  // resize the vectors
  m_points.resize(order);
  m_weights.resize(order);

  // integer division sets midpoint given even or odd n
  short unsigned int midpoint = (order + 1) / 2;

  // sum of weights to normalize at end of function
  long double weight_sum = 0;

  // j, j-1, and j-2 order legendre polynomials at mu
  long double p_j, p_jminus1, p_jminus2;

  // first derivative of nth order legendre
  long double p_deriv;

  // mu represents the independent for the legendre polynomial, which is
  // energy mapped to (-1, 1).
  // we will iterate on mu to find the roots of the nth legendre polynomial
  long double mu, old_mu;

  // loop indices.  i represents quad points, j the order of leg polynomials
  short unsigned int i, j;

  bool converged = false;

  // ====================================================
  // computations
  // ====================================================
  //

  for (i=0; i < midpoint; i++)
  {
    // set guess for mu; you'll have to ask kt about the logic!
    mu = cos(Constants::PI * (i + 0.75) / (order + 0.5));

    // loop until p_j at mu is sufficiently small to represent a root
    converged = false;
    while (!converged)
    {
      p_jminus1 = 0;
      p_j = 1;

      // set the nth order legendre polynomial at mu by recursion
      for (j=1; j<=order; j++)
      {
        p_jminus2 = p_jminus1;
        p_jminus1 = p_j;

        // this is the recursion relation for legendre polynomials
        p_j = ( (2*j - 1) * mu * p_jminus1 - (j - 1) * p_jminus2) / (j);
      }

      p_deriv = j * (mu * p_j - p_jminus1) / (mu*mu - 1);

      // iterating over mu to find the ith root
      old_mu = mu;
      mu = old_mu - p_j / p_deriv;

      if ( abs(mu - old_mu) < m_accuracy)
        converged = true;

    } // while !converged

    // set points utilizing symmetry
    m_points[i] = - mu;
    m_points[order-1 - i] = mu;

    // set weights, also symmetric
    m_weights[i] = 1 / ( (1 - mu*mu) * p_deriv*p_deriv);
    m_weights[order-1 - i] = m_weights[i];

    weight_sum += m_weights[i] + m_weights[order - 1 - i];
    if (i == order - 1 - i)
      weight_sum -= m_weights[i];


  } // for i to midpoint

  // normailize weights so that they sum to the length of the interval
  for (i=0; i<order; i++)
    m_weights[i] *= 2 / weight_sum;
}
