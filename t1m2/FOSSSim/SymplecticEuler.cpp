#include "SymplecticEuler.h"
#include <iostream>
bool SymplecticEuler::stepScene( TwoDScene& scene, scalar dt )
{
  // Your code goes here!
  // A vector containing all of the system's position DoFs. x0, y0, x1, y1, ...
  VectorXs& x = scene.getX();
  // A vector containing all of the system's velocity DoFs. v0, v0, v1, v1, ...
  VectorXs& v = scene.getV();
  
  // Create the vector to store the forces
  VectorXs F = VectorXs(2*scene.getNumParticles());
  for(int i=0; i< F.size(); i++){
    F[i] = 0.0;
  }
  const VectorXs& m = scene.getM();
  //Get the forces from the scene
  scene.accumulateGradU(F);
  
  VectorXs v_next(v.size());
  
  for(int i=0; i< x.size(); i=i+2){
    if( !scene.isFixed(i/2) ){        
        v_next(i) = (v[i]+(-1*dt*F[i]/m[i]));
        v_next(i+1) = (v[i+1]+(-1*dt*F[i+1]/m[i+1]));
    }
  }
  
    
  //Set the new velocities 
  for(int i=0; i< x.size(); i=i+2){
    if( !scene.isFixed(i/2) ){
        Vector2s new_v(2);
        new_v[0] = v_next[i];
        new_v[1] = v_next[i+1];
        scene.setVelocity(i/2,new_v);        
    }
  }
  
  
  //Set the new positions 
  for(int i=0; i< x.size(); i=i+2){
    if( !scene.isFixed(i/2) ){
        Vector2s p_next(2);
        p_next(0) = (x[i] + (dt*v_next[i]));
        p_next(1) = (x[i+1] + (dt*v_next[i+1]));
        scene.setPosition(i/2,p_next);
    }
  }

  return true;
}






