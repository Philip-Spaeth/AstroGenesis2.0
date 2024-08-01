#pragma once
#include <iostream>
#include "Particle.h"
#include <memory>
#include <vector>
#include "Constants.h"
#include "TimeIntegration.h"
#include "DataManager.h"
#include "random.h"
#include "Tree\Node.h"
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

    //simulation parameters
    double numberOfParticles = 1250;

    //adaptive time integration
    const double eta = 0.1;      // Accuracy parameter for adaptive time step
    const double maxTimeStep = 1e13; // Maximum allowed time step
    const double minTimeStep = 1e11; // Minimum allowed time step

    double globalTime = 0.0; // global time of the simulation in s
    const double endTime = 5e15; //end time of the simulation in s

    //save data at each maxTimeStep
    const double fixedTimeSteps = 200; //number of fixed maxtime steps
    const double fixedStep = endTime / fixedTimeSteps; //time step in s

    //gravitational softening, adapt it to the size of the system
    const double softening = 7e17; //softening factor
    
    //SPH parameters
    const double massInH = 1e18; //in kg

    //dark energy
    const double H0 = 70; //Hubble constant in km/s/Mpc

    //octree
    double theta = 0.5;
    std::shared_ptr<Node> root;
    
    //particles
    std::vector<std::shared_ptr<Particle>> particles;

    //Total Energy of the system
    std::vector<double> totalPotentialEnergy;
    std::vector<double> totalKineticEnergy;
    std::vector<double> totalInternalEnergy;
    std::vector<double> totalEnergy;

    //pointers to modules
    std::shared_ptr<TimeIntegration> timeIntegration;
    std::shared_ptr<DataManager> dataManager;

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
};