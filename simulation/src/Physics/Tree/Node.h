#pragma once 

#include <iostream>
#include "Particle.h"
#include <memory>
#include <vector>

class Particle;

class Node : public std::enable_shared_from_this<Node>
{
public:
    Node();
    ~Node();

    void insert(std::shared_ptr<Particle> newParticle);
    int getOctant(std::shared_ptr<Particle> newParticle);

    void calculateGravityForce(std::shared_ptr<Particle> newparticle, double softening, double theta);

    double mH = 0;
    double mRho = 0;
    double mP = 0;
    vec3 mVel = vec3(0,0,0);

    //calculate the density for all particles
    void calcVisualDensity(double radiusDensityEstimation);

    //SPH only for gas particles, Stars and dark matter particles are not affected
    double gasMass = 0; //mass of gas particles
    void calcGasDensity(double massInH);
    //Medians for the SPH with Nodes and not Paticle to Particle
    void calcmSPHNode();
    void calcmH();
    void calcmDensity();
    void calcmVelocity();
    void calcmPressure();

    int depth;
    bool isLeaf;
    //if leaf node, the particle is stored here
    std::shared_ptr<Particle> particle;

    vec3 position;
    double radius;

    //Mass center of the node
    vec3 centerOfMass;
    double mass;
    //childparticles vector for the density calculation
    std::vector<std::shared_ptr<Particle>> childParticles;

    //Children of the node
    std::shared_ptr<Node> children[8];

    //parent of the node
    std::weak_ptr<Node> parent;
};
