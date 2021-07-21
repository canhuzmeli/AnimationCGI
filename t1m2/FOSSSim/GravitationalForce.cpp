#include "GravitationalForce.h"

void GravitationalForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size()%2 == 0 );
  assert( m_particles.first >= 0 );  assert( m_particles.first < x.size()/2 );
  assert( m_particles.second >= 0 ); assert( m_particles.second < x.size()/2 );

  // Your code goes here!
}

void GravitationalForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%2 == 0 );
  assert( m_particles.first >= 0 );  assert( m_particles.first < x.size()/2 );
  assert( m_particles.second >= 0 ); assert( m_particles.second < x.size()/2 );
  
  scalar distance_between_particles = sqrt(pow(x[m_particles.second*2] - x[m_particles.first*2], 2) + pow(x[m_particles.second*2+1] - x[m_particles.first*2+1], 2) * 1.0);
  scalar energy = (m_G*m[m_particles.first*2]*m[m_particles.second*2]) / (pow(distance_between_particles,2));
  scalar unit_x = (x[m_particles.first*2] - x[m_particles.second*2]) / distance_between_particles;
  scalar unit_y = (x[m_particles.first*2+1] - x[m_particles.second*2+1]) / distance_between_particles;

  gradE[m_particles.first*2] += (energy*unit_x);
  gradE[m_particles.first*2+1] += (energy*unit_y);
  
  gradE[m_particles.second*2] += (-energy*unit_x);
  gradE[m_particles.second*2+1] += (-energy*unit_y);
  // Your code goes here!
}
