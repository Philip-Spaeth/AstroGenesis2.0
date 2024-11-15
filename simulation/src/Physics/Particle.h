#pragma once

#include "vec3.h"
#include <memory>
#include "Tree/Node.h"

class Node;

class Particle
{
public:
    Particle(vec3 position,vec3 velocity,vec3 acceleration, double mass = 1){}

    Particle(){}
    ~Particle(){}

    // Particle properties in SI units
    vec3 position;
    vec3 velocity;
    vec3 acceleration;
    double mass = 0;

    //unique id
    unsigned int id;

    //adaptive time integration
    double timeStep = 0;
    double nextIntegrationTime = 0;

    //SPH only for gas particles, Stars and dark matter particles are not affected
    uint8_t type = 1; // 1 = star, 2 = gas, 3 = dark matter
    uint8_t galaxyPart = 1; // 1 = disk, 2 = bulge, 3 = halo

    // Fluid properties (SPH) only for gas particles
    //calculated by the Tree
    double h = 0; // smoothing length in m
    double rho = 0; //density in kg/m^3
    double P = 0; //pressure in Pascal
    double T = 0; //temperature in Kelvin
    //derivative of A for the time integration of the internal energy
    double dUdt = 0;
    double U = 0; //internal energy in J/kg

    //calculated for all particles not just SPH gas, just for visualization
    double visualDensity = 0;

    //Octree properties
    std::weak_ptr<Node> node;
};