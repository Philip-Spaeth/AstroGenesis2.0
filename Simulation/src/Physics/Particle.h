#pragma once

#include "vec3.h"
#include <memory>
#include "Tree\Node.h"

class Node;

class Particle
{
public:
    Particle(vec3 position,vec3 velocity,vec3 acceleration, double mass = 1, double h = 0, double density = 0,double pressure = 0,double temperature = 0,double internalEnergy = 0,double kineticEnergy = 0,double potentialEnergy = 0,double totalEnergy = 0) : position(position), velocity(velocity), acceleration(acceleration), mass(mass), h(h), density(density), pressure(pressure), temperature(temperature), internalEnergy(internalEnergy), kineticEnergy(kineticEnergy), potentialEnergy(potentialEnergy), totalEnergy(totalEnergy){}

    Particle() : position(0.0, 0.0, 0.0), velocity(0.0, 0.0, 0.0), acceleration(0.0, 0.0, 0.0), mass(1.0), h(0.0), density(0.0), pressure(0.0), temperature(0.0), internalEnergy(0.0), kineticEnergy(0.0), potentialEnergy(0.0), totalEnergy(0.0){}
    ~Particle(){}

    // Particle properties
    vec3 position;
    vec3 velocity;
    vec3 acceleration;
    double mass;

    //adaptive time integration
    double timeStep = 0;
    double nextIntegrationTime = 0;

    //SPH only for gas particles, Stars and dark matter particles are not affected
    int type = 1; // 1 = star, 2 = gas, 3 = dark matter

    // Fluid properties (SPH) only for gas particles
    //calculated by the Tree
    double h;
    double density;

    //for all particles, just for visualization
    double visualDensity = 0;

    //calcualted with the forces
    double pressure;
    double temperature;
    double totalViscosityTensor;
    double internalEnergy;

    // Energy properties
    double kineticEnergy;
    double potentialEnergy;
    double totalEnergy;

    //Octree properties
    std::weak_ptr<Node> node;
};
