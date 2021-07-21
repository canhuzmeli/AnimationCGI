#include "SpringForce.h"

void SpringForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );

  // Your code goes here!
}

void SpringForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );
  
  scalar distance_between_particles = sqrt(pow(x[m_endpoints.second*2] - x[m_endpoints.first*2], 2) + pow(x[m_endpoints.second*2+1] - x[m_endpoints.first*2+1], 2) * 1.0);
  
  
  scalar energy = m_k*(distance_between_particles - m_l0);
  scalar unit_x = (x[m_endpoints.first*2] - x[m_endpoints.second*2]) / distance_between_particles;
  scalar unit_y = (x[m_endpoints.first*2+1] - x[m_endpoints.second*2+1]) / distance_between_particles;
  
 
  //damping
  scalar velocity_diff_x = v[m_endpoints.first*2] - v[m_endpoints.second*2];
  scalar velocity_diff_y = v[m_endpoints.first*2+1] - v[m_endpoints.second*2+1];
  
  scalar damping_force_x= m_b*velocity_diff_x;
  scalar damping_force_y= m_b*velocity_diff_y;

  gradE[m_endpoints.first*2] += ((energy*unit_x)+damping_force_x);
  gradE[m_endpoints.first*2+1] += ((energy*unit_y) +damping_force_y);

  gradE[m_endpoints.second*2] += (-((energy*unit_x)+damping_force_x));
  gradE[m_endpoints.second*2+1] += (-((energy*unit_y)+damping_force_y));
  
  // Your code goes here!

}
