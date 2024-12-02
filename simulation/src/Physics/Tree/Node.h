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

    void insert(std::vector<std::shared_ptr<Particle>> particles);

    void insert(std::shared_ptr<Particle> newParticle);
    int getOctant(std::shared_ptr<Particle> newParticle);

    void calculateGravityForce(std::shared_ptr<Particle> newparticle, double softening, double theta);

    double medianH = 0; //median smoothing length
    double medianDensity = 0; //median density
    double medianPressure = 0; //median pressure
    vec3 medianVelocity = vec3(0,0,0); //median velocity

    //calculate the density for all particles
    void calcVisualDensity(double radiusDensityEstimation);

    //SPH only for gas particles, Stars and dark matter particles are not affected
    double gasMass = 0; //mass of gas particles
    void calcGasDensity(double massInH);
    //Medians for the SPH with Nodes and not Paticle to Particle
    void calcSPHNodeMedians();
    void calcMedianH();
    void calcMedianDensity();
    void calcMedianVelocity();
    void calcMedianPressure();

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
