#pragma once
#include <Eigen/Dense>
#include "constants.h"

// 1D Gauss-Legendre quadrature set
class GLQuad
{
private:
  int m_num_points;
  double m_norm;
  Eigen::VectorXd m_mu, m_wt;

  void build(double tolerance = 1.0e-12);

public:
  GLQuad(int num_points, double norm)
    : m_num_points(num_points), m_norm(norm), 
      m_mu(num_points), m_wt(num_points)
  {
    assert(num_points > 0);
    assert(norm > 0.0);
    build();
  }

  const Eigen::Ref<const Eigen::VectorXd> mu() const { return m_mu; }
  const Eigen::Ref<const Eigen::VectorXd> wt() const { return m_wt; }

  int num_points() { return m_num_points; }
  double norm() { return m_norm; }
};
