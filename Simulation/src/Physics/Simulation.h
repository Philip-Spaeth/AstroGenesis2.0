#pragma once

#include <iostream>
#include "Particle.h"
#include <memory>
#include <vector>
#include "Constants.h"
#include "TimeIntegration.h"
#include "DataManager.h"
#include "random.h"
#include "Node.h"
#include <thread>
#include <mutex>
#include <atomic>

using namespace math;

class Simulation
{
public:
    Simulation();
    ~Simulation();

    bool init();
    void run();

private:

    void buildTree();
    void calculateForces();
    double calcTreeWidth();
    void calcDensity();

    //multithreading
    void calculateForcesWorker();
    std::atomic<int> currentParticleIndex;
    std::mutex mutex;

    //SPH parameters
    double h = 1; //smoothing length


    double deltaTime = 1; //time step length
    double timeSteps = 40; //number of time steps

    double softening = 1e-3; //softening factor

    //Total Energy of the system
    std::vector<double> totalPotentialEnergy;
    std::vector<double> totalKineticEnergy;
    std::vector<double> totalInternalEnergy;
    std::vector<double> totalEnergy;

    //pointers to modules
    std::shared_ptr<TimeIntegration> timeIntegration;
    std::shared_ptr<DataManager> dataManager;

    //particles
    double numberOfParticles = 10000;
    std::vector<std::shared_ptr<Particle>> particles;

    //octree
    double theta = 0.5;
    std::shared_ptr<Node> root;
};