#include "SimpleCollisionHandler.h"
#include <iostream>
#include <set>

// BEGIN STUDENT CODE //


// Detects whether two particles are overlapping (including the radii of each)
// and approaching.
// If the two particles overlap and are approaching, returns true and sets 
// the vector n to be the vector between the first and second particle.
// Inputs:
//   scene: The scene data structure. The positions and radii of the particles
//          can be obtained from here.
//   idx1:  The index of the first particle. (Ie, the degrees of freedom
//          corresponding to this particle are entries 2*idx1 and 2*idx1+1 in
//          scene.getX().
//   idx2:  The index of the second particle.
// Outputs:
//   n: The vector between the two particles.
//   Returns true if the two particles overlap and are approaching.
bool SimpleCollisionHandler::detectParticleParticle(TwoDScene &scene, int idx1, int idx2, Vector2s &n)
{
    // First particle velocity
    VectorXs v1 = scene.getV().segment<2>(2*idx1);
    // Second particle velocity
    VectorXs v2 = scene.getV().segment<2>(2*idx2);
    
    //first particle centre coordinates
    VectorXs x1 = scene.getX().segment<2>(2*idx1); 
    //second particle centre coordinates
    VectorXs x2 = scene.getX().segment<2>(2*idx2); 

    //first particle radius
    double x1radius = scene.getRadius(idx1);
    //second particle radius
    double x2radius = scene.getRadius(idx2);
    double distance_between_particles = sqrt(pow(x2[0]- x1[0], 2) + pow(x2[1] - x1[1], 2));
    if( (x2-x1).dot((v2-v1)) < 0){
        //approaching
        if((distance_between_particles - x1radius - x2radius) < 0 ){

            n = (x2-x1);
            return true;
        }
    }
    
    return false;
}

// Detects whether a particle and an edge are overlapping (including the radii 
// of both) and are approaching.
// If the two objects overlap and are approaching, returns true and sets the 
// vector n to be the shortest vector between the particle and the edge.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   eidx:  The index of the edge. (Ie, the indices of particle with index e are
//          scene.getEdges()[e].first and scene.getEdges()[e].second.)
// Outputs:
//   n: The shortest vector between the particle and the edge.
//   Returns true if the two objects overlap and are approaching.
bool SimpleCollisionHandler::detectParticleEdge(TwoDScene &scene, int vidx, int eidx, Vector2s &n)
{
    VectorXs x1 = scene.getX().segment<2>(2*vidx);
    VectorXs x2 = scene.getX().segment<2>(2*scene.getEdges()[eidx].first);
    VectorXs x3 = scene.getX().segment<2>(2*scene.getEdges()[eidx].second);
    
    VectorXs v1 = scene.getV().segment<2>(2*vidx);
    VectorXs v2 = scene.getV().segment<2>(2*scene.getEdges()[eidx].first);
    VectorXs v3 = scene.getV().segment<2>(2*scene.getEdges()[eidx].second);
    
    // Your code goes here!
    scalar alpha = (x1-x2).dot(x3-x2) / pow((x3-x2).norm(),2);
    
    if(alpha < 0){
        alpha = 0;
    }
    if(alpha > 1){
        alpha = 1;
    }
    
    VectorXs x_alpha = (x2 + (alpha * (x3 - x2)));
    VectorXs distance_vector = x_alpha - x1;
    //first particle radius
    double x1radius = scene.getRadius(vidx);
    double x2radius = scene.getEdgeRadii()[eidx];
    VectorXs v_alpha = v2 + alpha*(v3 - v2);
    
    
    if( (v1 - v_alpha).dot(distance_vector) > 0){
		//approaching
		if((distance_vector.norm() - x1radius - x2radius) < 0 ){
            n = distance_vector;
            return true;
        }
    }
    
    return false;
}

// Detects whether a particle and a half-plane are overlapping (including the 
// radius of the particle) and are approaching.
// If the two objects overlap and are approaching, returns true and sets the 
// vector n to be the shortest vector between the particle and the half-plane.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   pidx:  The index of the halfplane. The vectors (px, py) and (nx, ny) can
//          be retrieved by calling scene.getHalfplane(pidx).
// Outputs:
//   n: The shortest vector between the particle and the half-plane.
//   Returns true if the two objects overlap and are approaching.
bool SimpleCollisionHandler::detectParticleHalfplane(TwoDScene &scene, int vidx, int pidx, Vector2s &n)
{
    VectorXs x1 = scene.getX().segment<2>(2*vidx);
    VectorXs px = scene.getHalfplane(pidx).first;
    VectorXs pn = scene.getHalfplane(pidx).second;
    
    VectorXs v1 = scene.getV().segment<2>(2*vidx);
    // Your code goes here!

    //first particle radius
    double x1radius = scene.getRadius(vidx);
    
    VectorXs distance_vector = ((px - x1).dot(pn) / pow(pn.norm(),2)) * pn;
    if(distance_vector.norm()- x1radius < 0){
        if(v1.dot(distance_vector) > 0){
            n = distance_vector;
            return true;
        }
        
    }
    return false;
}


// Responds to a collision detected between two particles by applying an impulse
// to the velocities of each one.
// You can get the COR of the simulation by calling getCOR().
// Inputs:
//   scene: The scene data structure.
//   idx1:  The index of the first particle.
//   idx2:  The index of the second particle.
//   n:     The vector between the first and second particle.
// Outputs:
//   None.
void SimpleCollisionHandler::respondParticleParticle(TwoDScene &scene, int idx1, int idx2, const Vector2s &n)
{
    const VectorXs &M = scene.getM();
    VectorXs &v = scene.getV();
    
    const scalar cor_var = ((getCOR()+1)/2);
    // First particle velocity
    VectorXs v1 = scene.getV().segment<2>(2*idx1);
    // Second particle velocity
    VectorXs v2 = scene.getV().segment<2>(2*idx2);
    scalar m1 = M[2*idx1];
    scalar m2 = M[2*idx2];
    
    if(scene.isFixed(idx1)){
        m1 = std::numeric_limits<double>::infinity();
    }
    
    if(scene.isFixed(idx2)){
        m2 = std::numeric_limits<double>::infinity();
    }
    Vector2s n_hat = (n/n.norm());
    
    if(!scene.isFixed(idx1)){
        
        Vector2s vel1new = v1 + (cor_var*((((2*(v2-v1)).dot(n_hat))/(1+(m1/m2)))*n_hat));
        scene.setVelocity(idx1,vel1new);
        
    }
    
    if(!scene.isFixed(idx2)){
        Vector2s vel2new = v2 - (cor_var*((((2*(v2-v1)).dot(n_hat))/(1+(m2/m1)))*n_hat));
        scene.setVelocity(idx2,vel2new);
        
    }
    
}

// Responds to a collision detected between a particle and an edge by applying
// an impulse to the velocities of each one.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   eidx:  The index of the edge.
//   n:     The shortest vector between the particle and the edge.
// Outputs:
//   None.
void SimpleCollisionHandler::respondParticleEdge(TwoDScene &scene, int vidx, int eidx, const Vector2s &n)
{
    const VectorXs &M = scene.getM();
    const scalar cor_var = ((getCOR()+1)/2);
    int eidx1 = scene.getEdges()[eidx].first;
    int eidx2 = scene.getEdges()[eidx].second;
    
    VectorXs x1 = scene.getX().segment<2>(2*vidx);
    VectorXs x2 = scene.getX().segment<2>(2*eidx1);
    VectorXs x3 = scene.getX().segment<2>(2*eidx2);
    
    VectorXs v1 = scene.getV().segment<2>(2*vidx);
    VectorXs v2 = scene.getV().segment<2>(2*eidx1);
    VectorXs v3 = scene.getV().segment<2>(2*eidx2);
 
    
    scalar m1 = M[2*vidx];
    scalar m2 = M[2*eidx1];
    scalar m3 = M[2*eidx2];

    if(scene.isFixed(vidx)){
        m1 = std::numeric_limits<double>::infinity();
    }
    
    if(scene.isFixed(eidx1)){
        m2 = std::numeric_limits<double>::infinity();
    }
    
    if(scene.isFixed(eidx2)){
        m3 = std::numeric_limits<double>::infinity();
    }

    scalar alpha = (x1-x2).dot(x3-x2) / pow((x3-x2).norm(),2);
    if(alpha < 0){
        alpha = 0;
    }
    if(alpha > 1){
        alpha = 1;
    }
    VectorXs ve = (1-alpha)*v2 + (alpha*v3);
    Vector2s n_hat = (n/n.norm());
    if(!scene.isFixed(vidx)){
        scalar upper = (2*(ve-v1)).dot(n_hat);
        scalar divider = (1 + ((m1*(pow((1-alpha),2)))/m2) +((pow(alpha,2)*m1)/m3));
        Vector2s vel1new = v1 + cor_var*(( upper / divider)*n_hat);
        scene.setVelocity(vidx,vel1new);
        
    }
    
    if(!scene.isFixed(eidx1)){
        Vector2s vel2new = v2 - cor_var*(((2*(ve-v1)*(1-alpha)).dot(n_hat) / ((m2/m1) + (pow((1-alpha),2)) +((pow(alpha,2)*m2)/m3)))*n_hat);
        scene.setVelocity(eidx1,vel2new);
        
    }
    
    if(!scene.isFixed(eidx2)){
        Vector2s vel3new = v3 - cor_var*(((2*(ve-v1)*alpha).dot(n_hat) / ((m3/m1) + ((m3*(pow((1-alpha),2)))/m2) +(pow(alpha,2))))*n_hat);
        scene.setVelocity(eidx2,vel3new);        
    }
    
}


// Responds to a collision detected between a particle and a half-plane by 
// applying an impulse to the velocity of the particle.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   pidx:  The index of the half-plane.
//   n:     The shortest vector between the particle and the half-plane.
// Outputs:
//   None.
void SimpleCollisionHandler::respondParticleHalfplane(TwoDScene &scene, int vidx, int pidx, const Vector2s &n)
{
    const scalar cor_var = ((getCOR()+1)/2);
    Vector2s nhat = (n/n.norm());

    // Your code goes here!
    VectorXs v1 = scene.getV().segment<2>(2*vidx);
    if(!scene.isFixed(vidx)){
        
        Vector2s vel1new = v1 - cor_var*(2*(v1.dot(nhat))*nhat);
        scene.setVelocity(vidx,vel1new);
        
    }
    
    
}
