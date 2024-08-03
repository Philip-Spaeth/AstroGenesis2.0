#pragma once

#include <iostream>
#include "Particle.h"
#include <memory>


class TimeIntegration
{
public:
    TimeIntegration(){}
    ~TimeIntegration(){}

    //Semi implicit Euler
    void Euler(std::shared_ptr<Particle> particle, double deltaTime);

    //Leapfrog
    void Kick(std::shared_ptr<Particle> particle, double deltaTime);
    void Drift(std::shared_ptr<Particle> particle, double deltaTime);

    //integrate the entropy
    void EntropyEuler(std::shared_ptr<Particle> particle, double deltaTime);
};