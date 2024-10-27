#pragma once

#include "vec3.h"
#include <memory>
//check if windows or linux
#ifdef _WIN32
#include "Tree\Node.h"
#else
#include "Tree/Node.h"
#endif

class Node;

class Particle
{
public:
    Particle(vec3 position,vec3 velocity,vec3 acceleration, double mass = 1){}

    Particle(){}
    ~Particle(){}

    // Particle properties
    vec3 position;
    vec3 velocity;
    vec3 acceleration;
    double mass = 0;

    //adaptive time integration
    double timeStep = 0;
    double nextIntegrationTime = 0;

    //SPH only for gas particles, Stars and dark matter particles are not affected
    int type = 1; // 1 = star, 2 = gas, 3 = dark matter

    // Fluid properties (SPH) only for gas particles
    //calculated by the Tree
    double h = 0;
    double rho = 0; //density in kg/m^3
    //calcualted with the forces
    double P = 0; //pressure in Pascal
    double T = 0; //temperature in Kelvin
    //derivative of A for the time integration of the entropy
    double dUdt = 0;
    double U = 0; //internal energy in J

    //double totalViscosityTensor;


    //for all particles, just for visualization
    double visualDensity = 0;

    // Energy properties
    double kineticEnergy = 0;
    double potentialEnergy= 0;

    //Octree properties
    std::weak_ptr<Node> node;
};
