#include "SimpleGravityForce.h"

void SimpleGravityForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
    std::cout << "E:" << std::endl;
    std::cout << E << std::endl;
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size()%2 == 0 );
    for(int i=0; i< x.size(); i=i+2){
        E = E+(m[i]*SimpleGravityForce::m_gravity[1]*x[i+1]);
    }
    std::cout << "E NEW:" << std::endl;
    std::cout << E << std::endl;
    // Your code goes here!
}

void SimpleGravityForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{

    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size() == gradE.size() );
    assert( x.size()%2 == 0 );
    for(int i=0; i< x.size(); i=i+2){
        gradE[i] = -(SimpleGravityForce::m_gravity[0]);
        gradE[i+1] = -(SimpleGravityForce::m_gravity[1]);
    }
    // Your code goes here!
}
