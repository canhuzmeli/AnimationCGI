#include "ExplicitEuler.h"
#include <iostream>
bool ExplicitEuler::stepScene( TwoDScene& scene, scalar dt )
{
    // Your code goes here!
    
    // Some tips on getting data from TwoDScene:
    // A vector containing all of the system's position DoFs. x0, y0, x1, y1, ...
    VectorXs& x = scene.getX();
    // A vector containing all of the system's velocity DoFs. v0, v0, v1, v1, ...
    VectorXs& v = scene.getV();
    // A vector containing the masses associated to each DoF. m0, m0, m1, m1, ...
    const VectorXs& m = scene.getM();
    
    // Create the vector to store the forces
    VectorXs F = VectorXs(2*scene.getNumParticles());
  
    //Get the forces from the scene
    scene.accumulateGradU(F);
    
    //Set the new positions 
    for(int i=0; i< x.size(); i=i+2){
      if( !scene.isFixed(i/2) ){
          Vector2s p_next(2);
          p_next(0) = (x[i] + (dt*v[i]));
          p_next(1) = (x[i+1] + (dt*v[i+1]));
          scene.setPosition(i/2,p_next);
      }
    }
    
    for(int i=0; i< x.size(); i=i+2){
      if( !scene.isFixed(i/2) ){
          Vector2s v_next(2);
          v_next(0) = (v[i]+dt*(-1*F[i]));
          v_next(1) = (v[i+1]+dt*(-1*F[i+1]));
          scene.setVelocity(i/2,v_next);
      }
    }
    /*
     * std::cout << "x is of size " << x.rows() << "x" << x.cols() << std::endl;
     * std::cout << "v is of size " << v.rows() << "x" << v.cols() << std::endl;
     * std::cout << "m is of size " << m.rows() << "x" << m.cols() << std::endl;
     * x is of size 18x1
     * v is of size 18x1
     * m is of size 18x1
     */
    
    return true;
}
