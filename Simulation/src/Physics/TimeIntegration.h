#pragma once

#include <iostream>
#include "Particle.h"
#include <memory>


class TimeIntegration
{
public:
    TimeIntegration();
    ~TimeIntegration();

    //Semi implicit Euler
    void Euler(std::shared_ptr<Particle> particle, double deltaTime);
};