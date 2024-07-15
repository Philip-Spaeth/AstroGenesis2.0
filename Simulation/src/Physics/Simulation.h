#pragma once

#include <iostream>
#include "Particle.h"
#include <memory>
#include <vector>
#include "Constants.h"
#include "TimeIntegration.h"

using namespace math;

class Simulation
{
public:
    Simulation();
    ~Simulation();

    void run();


private:
    double deltaTime = 1e2; //time step length
    double timeSteps = 1000; //number of time steps

    double softening = 1e-3; //softening factor

    //pointers to modules
    std::shared_ptr<TimeIntegration> timeIntegration = std::make_shared<TimeIntegration>();

    //particles
    double numberOfParticles = 5;
    std::vector<std::shared_ptr<Particle>> particles;
};