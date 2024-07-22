#pragma once

#include <iostream>
#include "Particle.h"
#include <memory>
#include <vector>
#include "Constants.h"
#include "TimeIntegration.h"
#include "DataManager.h"
#include "random.h"

using namespace math;

class Simulation
{
public:
    Simulation();
    ~Simulation();

    void run();


private:
    double deltaTime = 10; //time step length
    double timeSteps = 10; //number of time steps

    double softening = 1e-3; //softening factor

    //Total Energy of the system
    std::vector<double> totalPotentialEnergy;
    std::vector<double> totalKineticEnergy;
    std::vector<double> totalInternalEnergy;
    std::vector<double> totalEnergy;

    //pointers to modules
    std::shared_ptr<TimeIntegration> timeIntegration = std::make_shared<TimeIntegration>();
    std::shared_ptr<DataManager> dataManager = std::make_shared<DataManager>("../../../Data/data1/");

    //particles
    double numberOfParticles = 100;
    std::vector<std::shared_ptr<Particle>> particles;
};
