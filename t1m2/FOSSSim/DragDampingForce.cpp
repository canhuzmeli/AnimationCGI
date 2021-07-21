#include "DragDampingForce.h"

void DragDampingForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%2 == 0 );
  for(int i=0; i< x.size(); i=i+2){
     gradE[i] += (v[i]*m_b);
     gradE[i+1] += (v[i+1]*m_b);
  }
  // Your code goes here!
}
