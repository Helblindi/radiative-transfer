#include "GLQuad.h"


void GLQuad::build(double tolerance)
{
  m_mu.resize(m_num_points);
  m_wt.resize(m_num_points);

  double x1  = -1.0;
  double x2  =  1.0;
  double xm  =  0.5 * (x2 + x1);
  double xl  =  0.5 * (x2 - x1);
  double dnp =  static_cast<double>(m_num_points); 

  int m = (m_num_points+1)/2;

  for(int i = 1;i <= m;++i)
  {
    double di  = static_cast<double>(i);

    double z1, pp, z = cos(Constants::PI * (di - 0.25)/(dnp + 0.5));
    do
    {
      double p1 = 1.0, p2 = 0.0;
      for(int j = 1;j <= m_num_points;++j)
      {
        double dj = static_cast<double>(j);
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0*dj - 1.0)*z*p2 - (dj - 1.0)*p3)/dj;
      }
      pp = dnp*(z*p1 - p2)/(z*z - 1.0);
      z1 = z;
      z = z1 - p1/pp;
    }
    while(std::fabs(z - z1) > tolerance);

    m_mu(i - 1) = xm - xl*z;
    m_mu(m_num_points - i) = xm + xl*z;

    m_wt(i - 1) = m_norm*xl/((1.0 - z*z)*pp*pp);
    m_wt(m_num_points - i) =  m_wt(i - 1);
  }
}
