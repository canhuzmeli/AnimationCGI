#include "PenaltyForce.h"
#include "TwoDScene.h"
using Eigen::MatrixXd;

void PenaltyForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
    // Feel free to implement if you feel like doing so.
}

// Adds the gradient of the penalty potential (-1 * force) for a pair of 
// particles to the total.
// Read the positions of the particles from the input variable x. Radii can
// be obtained from the member variable m_scene, the penalty force stiffness 
// from member variable m_k, and penalty force thickness from member variable
// m_thickness.
// Inputs:
//   x:    The positions of the particles in the scene. 
//   idx1: The index of the first particle, i.e. the position of this particle
//         is ( x[2*idx1], x[2*idx1+1] ).
//   idx2: The index of the second particle.
// Outputs:
//   gradE: The total gradient of penalty force. *ADD* the particle-particle
//          gradient to this total gradient.
void PenaltyForce::addParticleParticleGradEToTotal(const VectorXs &x, int idx1, int idx2, VectorXs &gradE)
{
    VectorXs x1 = x.segment<2>(2*idx1);
    VectorXs x2 = x.segment<2>(2*idx2);
    
    scalar r1 = m_scene.getRadius(idx1);
    scalar r2 = m_scene.getRadius(idx2);
    
    // First particle velocity
    VectorXs v1 = m_scene.getV().segment<2>(2*idx1);
    // Second particle velocity
    VectorXs v2 = m_scene.getV().segment<2>(2*idx2);
    std::cout << "v1: "  << v1 << std::endl;
    std::cout << "v2: "  << v2 << std::endl;
    // Your code goes here!
    VectorXs n = x2-x1;
    scalar distance = n.norm();
    if(distance <= r1+r2+m_thickness){
        const int size = 2;
        MatrixXs part1(size, size);
        part1 = -1.0*MatrixXs::Identity(size, size);
        MatrixXs part2(size, size);
        part2 = MatrixXs::Identity(size, size);
        MatrixXs n_gradient(size, size*2);
        n_gradient << part1,part2;

        VectorXs F = -m_k*(distance - r1 - r2 - m_thickness)*n_gradient.transpose()*(n/distance);
        std::cout << "F: "  <<  std::endl;
        std::cout << F << std::endl;
        std::cout << "gradE: "  <<  std::endl;
        std::cout << gradE << std::endl;
        gradE[2*idx1] = gradE[2*idx1] - F[0];
        gradE[2*idx1+1] = gradE[2*idx1+1] - F[1];
        gradE[2*idx2] = gradE[2*idx2] - F[2];
        gradE[2*idx2+1] = gradE[2*idx2+1] - F[3];
        //gradE = gradE-F;
    }
}

// Adds the gradient of the penalty potential (-1 * force) for a particle-edge
// pair to the total.
// Read the positions of the particle and edge endpoints from the input
// variable x.
// Inputs:
//   x:    The positions of the particles in the scene.
//   vidx: The index of the particle.
//   eidx: The index of the edge, i.e. the indices of the particle making up the
//         endpoints of the edge are given by m_scene.getEdge(eidx).first and 
//         m_scene.getEdges(eidx).second.
// Outputs:
//   gradE: The total gradient of penalty force. *ADD* the particle-edge
//          gradient to this total gradient.
void PenaltyForce::addParticleEdgeGradEToTotal(const VectorXs &x, int vidx, int eidx, VectorXs &gradE)
{
    VectorXs x1 = x.segment<2>(2*vidx);
    VectorXs x2 = x.segment<2>(2*m_scene.getEdge(eidx).first);
    VectorXs x3 = x.segment<2>(2*m_scene.getEdge(eidx).second);
    
    scalar r1 = m_scene.getRadius(vidx);
    scalar r2 = m_scene.getEdgeRadii()[eidx];
    
    // Your code goes here!
    scalar alpha = (x1-x2).dot(x3-x2) / pow((x3-x2).norm(),2);
    
    if(alpha < 0){
        alpha = 0;
    }
    if(alpha > 1){
        alpha = 1;
    }
    VectorXs x_alpha = (x2 + (alpha * (x3 - x2)));
    VectorXs n = x_alpha - x1;
    if(n.norm() <= r1+r2+m_thickness){
        
        const int size = 2;
        MatrixXs part1(size, size);
        part1 = -1*MatrixXd::Identity(size, size);

        MatrixXs part2(size, size);
        part2 = ((1-alpha)*MatrixXd::Identity(size, size));
        MatrixXs part3(size, size);
        part3 = (alpha*MatrixXd::Identity(size, size));
        MatrixXs n_gradient(size, size*3);
        n_gradient << part1,part2,part3;
        VectorXs F = m_k*(n.norm() - r1 - r2 - m_thickness)*n_gradient.transpose()*(n/n.norm());
        gradE = gradE+F;
    }
    
}

// Adds the gradient of the penalty potential (-1 * force) for a particle-
// half-plane pair to the total.
// Read the positions of the particle from the input variable x.
// Inputs:
//   x:    The positions of the particles in the scene.
//   vidx: The index of the particle.
//   pidx: The index of the half-plane, i.e. the position and normal vectors
//         for the half-plane can be retrieved by calling
//         m_scene.getHalfplane(pidx).
// Outputs:
//   gradE: The total gradient of the penalty force. *ADD* the particle-
//          half-plane gradient to this total gradient.
void PenaltyForce::addParticleHalfplaneGradEToTotal(const VectorXs &x, int vidx, int pidx, VectorXs &gradE)
{
    VectorXs x1 = m_scene.getX().segment<2>(2*vidx);
    VectorXs nh = m_scene.getHalfplane(pidx).second;
    
    // Your code goes here!
    VectorXs px = m_scene.getHalfplane(pidx).first;
    
    scalar r1 = m_scene.getRadius(vidx);


    VectorXs n = ((px - x1).dot(nh) / pow(nh.norm(),2)) * nh;
    

    
    if(n.norm() <= r1+m_thickness){

        MatrixXs n_gradient = -1*(nh*(nh.transpose()))/(pow(nh.norm(),2));

        VectorXs F = m_k*(n.norm() - r1 - m_thickness)*n_gradient.transpose()*(n/n.norm());
        std::cout << "F: "  << std::endl;
        std::cout << F << std::endl;
        gradE = gradE+F;
    }
}
