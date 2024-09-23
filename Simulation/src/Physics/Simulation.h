#pragma once
#include <iostream>
#include "Particle.h"
#include <memory>
#include <vector>
#include "Constants.h"
#include "TimeIntegration.h"
#include "DataManager.h"
#include "ICDataReader.h"
#include "Console.h"
#include "random.h"
//check if windows or linux
#ifdef _WIN32
#include "Tree\Node.h"
#include "Tree\Tree.h"
#else
#include "Tree/Node.h"
#include "Tree/Tree.h"
#endif
#include <thread>
#include <mutex>
#include <atomic>

class Tree;
class DataManager;

class Simulation
{
public:
    Simulation();
    ~Simulation();
    bool init();
    void run();

    //simulation parameters, has to the same as in the input dataset(ICs)
    double numberOfParticles;

    //adaptive time integration
    double eta;      // Accuracy parameter for adaptive time step
    double maxTimeStep; // Maximum allowed time step
    double minTimeStep; // Minimum allowed time step

    double globalTime; // global time of the simulation in s
    double endTime; //end time of the simulation in s

    //save data at each maxTimeStep
    double fixedTimeSteps; //number of fixed maxtime steps
    double fixedStep; //time step in s

    //gravitational softening, adapt it to the size of the system, osftening beginns at 2.8 * e0
    double e0; //softening factor
    
    //SPH parameters
    double massInH; //in kg

    //Visual density, for all particles, just for visualization, beacuse the real density is only calculated for Gas particles, has no physical meaning
    double visualDensityRadius; //in m

    //Constant hubble expansion
    double H0; //Hubble constant in km/s/Mpc

    //octree with all particles
    double theta;

    std::string ICFileName;
    std::string ICFileFormat;
    
    //particles
    std::vector<std::shared_ptr<Particle>> particles;

private:

    //pointers to modules
    std::shared_ptr<TimeIntegration> timeIntegration;
    std::shared_ptr<DataManager> dataManager;
    std::shared_ptr<ICDataReader> icDataReader;
    std::shared_ptr<Console> console;

    //SPH
    void initGasParticleProperties(std::shared_ptr<Tree> tree); // update A, U, P after the tree is built and rho is calculated
    void updateGasParticleProperties(std::shared_ptr<Tree> tree); // update A, T, U, P

    //calculations without the octree, just for debugging purposes
    void calculateForcesWithoutOctree(std::shared_ptr<Particle> p);
};
