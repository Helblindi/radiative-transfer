#pragma once

#include "NDArray.h"

// 1D Gauss-Legendre quadrature set
class GLQuad
{
public:

  GLQuad(int num_points, double norm)
    : m_num_points(num_points), m_norm(norm)
  {
    assert(num_points > 0);
    assert(norm > 0.0);
    build();
  }

  const ndarray::array<double, 1> & mu() const { return m_mu; }

  const ndarray::array<double, 1> & wt() const { return m_wt; }

  int num_points() { return m_num_points; }

  double norm() { return m_norm; }

private:

  int m_num_points;
  double m_norm;
  ndarray::array<double, 1> m_mu, m_wt;

  void build(double tolerance = 1.0e-12);
};
