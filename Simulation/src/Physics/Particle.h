#pragma once

#include "vec3.h"

class Particle
{
public:
    Particle(vec3 position,vec3 velocity,vec3 acceleration, double mass = 1,double density = 0,double pressure = 0,double temperature = 0,double internalEnergy = 0,double kineticEnergy = 0,double potentialEnergy = 0,double totalEnergy = 0) : position(position), velocity(velocity), acceleration(acceleration), mass(mass), density(density), pressure(pressure), temperature(temperature), internalEnergy(internalEnergy), kineticEnergy(kineticEnergy), potentialEnergy(potentialEnergy), totalEnergy(totalEnergy){}

    Particle() : position(0.0, 0.0, 0.0), velocity(0.0, 0.0, 0.0), acceleration(0.0, 0.0, 0.0), mass(1.0), density(0.0), pressure(0.0), temperature(0.0), internalEnergy(0.0), kineticEnergy(0.0), potentialEnergy(0.0), totalEnergy(0.0){}
    ~Particle(){}

    // Particle properties
    vec3 position;
    vec3 velocity;
    vec3 acceleration;
    double mass;

    // Fluid properties (SPH)
    double density;
    double pressure;
    double temperature;
    double viscosity;
    bool type = 1; // 1 = star, 2 = gas, 3 = dark matter

    // Energy properties
    double internalEnergy;
    double kineticEnergy;
    double potentialEnergy;
    double totalEnergy;
};
