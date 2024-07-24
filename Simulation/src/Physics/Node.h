#pragma once 

#include <iostream>
#include "Particle.h"
#include <memory>

class Particle;

class Node : public std::enable_shared_from_this<Node>
{
public:
    Node();
    ~Node();

    void insert(std::shared_ptr<Particle> newParticle);
    int getOctant(std::shared_ptr<Particle> newParticle);

    void calculateForce(std::shared_ptr<Particle> newparticle, double softening, double theta);

    //SPH
    double calcDensity(double h);

    int depth;
    bool isLeaf;
    //if leaf node, the particle is stored here
    std::shared_ptr<Particle> particle;

    vec3 position;
    double radius;

    //Mass center of the node
    vec3 centerOfMass;
    double mass;

    //Children of the node
    std::shared_ptr<Node> children[8];

    //parent of the node
    std::weak_ptr<Node> parent;
};
