#pragma once
#include <iostream>
#include "Particle.h"
#include <memory>
#include <vector>
#include "Constants.h"
#include "TimeIntegration.h"
#include "DataManager.h"
#include "random.h"
//check if windows or linux
#ifdef _WIN32
#include "Tree\Node.h"
#else
#include "Tree/Node.h"
#endif
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
    const double eta = 10;      // Accuracy parameter for adaptive time step
    const double maxTimeStep = 1e13; // Maximum allowed time step
    const double minTimeStep = 1e13; // Minimum allowed time step

    double globalTime = 0.0; // global time of the simulation in s
    const double endTime = 1e16; //end time of the simulation in s

    //save data at each maxTimeStep
    const double fixedTimeSteps = 100; //number of fixed maxtime steps
    const double fixedStep = endTime / fixedTimeSteps; //time step in s

    //gravitational softening, adapt it to the size of the system
    const double softening = 7e17; //softening factor
    
    //SPH parameters
    const double massInH = 1e40; //in kg

    //Visual density, for all particles, just for visualization, has no physical meaning
    const double visualDensityRadius = 1e19; //in m

    //dark energy
    const double H0 = 70; //Hubble constant in km/s/Mpc

    //octree with all particles
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
    void calcVisualDensity();
    //SPH
    void initGasParticleProperties(); // update A, U, P after the tree is built and rho is calculated
    void updateGasParticleProperties(); // update A, T, U, P
    void calcGasDensity();

    //calculations without the octree
    void calculateForcesWithoutOctree(std::shared_ptr<Particle> p);

    //calculations without the octree
    void applyHubbleExpansion();

    //multithreading
    void calculateForcesWorker();
    std::atomic<int> currentParticleIndex;
    std::mutex mutex;
};
