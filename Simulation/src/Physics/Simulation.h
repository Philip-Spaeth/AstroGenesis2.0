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
    void applyHubbleExpansion();
    void apply_cosmological_expansion();
    double H(double z);

    //for the templates
    void readTemplate();

    //multithreading
    void calculateForcesWorker();
    std::atomic<int> currentParticleIndex;
    std::mutex mutex;

    //SPH parameters
    double h = 1e20; //smoothing length

    //dark energy
    double H0 = 1e-7; //Hubble constant in km/s/Mpc
    double Omega_m = 0.3; //matter density
    double Omega_Lambda = 0.7; //dark energy density
    double z = 1; //redshift
    double a = 1 / (1 + z); //scale factor

    //for the periodic boundary conditions
    double box_size = 1e24;

    //time integration
    double deltaTime = 5e15; //time step length
    double timeSteps = 1000; //number of time steps

    //gravitational softening, adapt it to the size of the system
    double softening = 7.715e22; //softening factor

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