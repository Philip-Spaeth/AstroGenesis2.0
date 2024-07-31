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

class Simulation
{
public:
    Simulation();
    ~Simulation();

    bool init();
    void run();

private:
    //calculations with the octree
    void buildTree();
    void calculateForces();
    double calcTreeWidth();
    void calcDensity();

    //calculations without the octree
    void calculateForcesWithoutOctree(std::shared_ptr<Particle> p);

    //calculations without the octree
    void applyHubbleExpansion();

    //for the templates
    void readTemplate();

    //multithreading
    void calculateForcesWorker();
    std::atomic<int> currentParticleIndex;
    std::mutex mutex;

    //SPH parameters
    const double h = 20; //smoothing length

    //dark energy
    const double H0 = 70; //Hubble constant in km/s/Mpc

    //adaptive time integration
    const double eta = 0.1;      // Accuracy parameter for adaptive time step
    const double maxTimeStep = 1; // Maximum allowed time step
    const double minTimeStep = 1; // Minimum allowed time step

    double globalTime = 0.0; // global time of the simulation in s
    const double endTime = 100; //end time of the simulation in s

    //save data at each maxTimeStep
    const double fixedTimeSteps = 100; //number of fixed maxtime steps
    const double fixedStep = endTime / fixedTimeSteps; //time step in s

    //gravitational softening, adapt it to the size of the system
    const double softening = 1; //softening factor

    //Total Energy of the system
    std::vector<double> totalPotentialEnergy;
    std::vector<double> totalKineticEnergy;
    std::vector<double> totalInternalEnergy;
    std::vector<double> totalEnergy;

    //pointers to modules
    std::shared_ptr<TimeIntegration> timeIntegration;
    std::shared_ptr<DataManager> dataManager;

    //particles
    double numberOfParticles = 500;
    std::vector<std::shared_ptr<Particle>> particles;

    //octree
    double theta = 1;
    std::shared_ptr<Node> root;
};