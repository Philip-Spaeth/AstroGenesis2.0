#include "Simulation.h"
#include <numeric>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include "Units.h"


Simulation::Simulation()
{
    //construct the modules
    timeIntegration = std::make_shared<TimeIntegration>();
    dataManager = std::make_shared<DataManager>("../../output_data/");
    icDataReader = std::make_shared<ICDataReader>();
    console = std::make_shared<Console>();
}

Simulation::~Simulation(){}

bool Simulation::init()
{
    //load the config file
    if (!dataManager->loadConfig("../Config.ini", this))
    {
        std::cerr << "Error: Could not load the config file." << std::endl;
        return false;
    }
    
    fixedStep = endTime / fixedTimeSteps;

    std::cout << "Total Number of Particles in the Config.ini file: " << numberOfParticles << std::endl;


//catch errors 
    //check if minTimeStep is smaller or equal to maxTimeStep
    if (minTimeStep > maxTimeStep)
    {
        std::cerr << "Error: minTimeStep is greater than maxTimeStep." << std::endl;
        return false;
    }
    //ckeck if the end time  / minTimeStep is < fixedTimeSteps
    if (endTime / minTimeStep < fixedTimeSteps)
    {
        std::cerr << "Error: endTime / minTimeStep is smaller than fixedTimeSteps." << std::endl;
        return false;
    }

    //print the general information aboput the simulation
    std::cout << "\nSimulation parameters:" << std::endl;
    std::cout << "  Number of particles: " << numberOfParticles << std::endl;
    std::cout << "  End time: " << std::scientific << std::setprecision(0) << (double)endTime / (double)3600.0 << " years" << std::endl;

    //print the computers / server computational parameters like number of threads, ram, cpu, etc.
    Console::printSystemInfo();
    /*
    if(ICFileFormat == "Gadget2")
    {
        icDataReader->readGadget2Snapshot(ICFileName, particles);
    }
    else if(ICFileFormat == "ASCII")
    {
        icDataReader->readASCII(ICFileName, 0, numberOfParticles, vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 0.0), particles);
    }
    else
    {
        std::cerr << "Error: Unknown file format." << std::endl;
        return false;
    }*/
    
    Galaxy galaxy(&particles);

    //Bulge
    galaxy.M_Bulge = 1e39;
    galaxy.R_Bulge = 5e19;
    galaxy.Rs_Bulge = 5e18;
    galaxy.N_Bulge = 1000;
    
    //Disk
    galaxy.M_Disk = 1e40;
    galaxy.R_Disk = 1e20;
    galaxy.z_Disk = 2e18;
    galaxy.VelDis_Disk = 1e3;
    galaxy.N_Disk = 9000;

    //Gas in the disk
    galaxy.M_Gas = 1e40;
    galaxy.R_Gas = 1e20;
    galaxy.z_Gas = 1e19;
    galaxy.N_Gas = 0;

    //Dark Matter Halo
    galaxy.M_Halo = 8e40;
    galaxy.R_Halo = 1e21;
    galaxy.c_Halo = 7;
    galaxy.N_Halo = 3000;

    galaxy.galaxyPosition = vec3(0.0, 0.0, 0.0);
    galaxy.galaxyVelocity = vec3(0.0, 0.0, 0.0);

    galaxy.generateGalaxy();

    if(numberOfParticles != particles.size())
    {
        std::cerr << "Error: Number of particles in the ConfigFile does not match the number of particles in the data file." << std::endl;
        std::cout << "Number of particles in the ConfigFile: " << numberOfParticles << std::endl;
        std::cout << "Number of particles in the data file: " << particles.size() << std::endl;
        return false;
    }

    //check if there are null pointers in the particles vector
    for (int i = 0; i < numberOfParticles; i++)
    {
        if (!particles[i]) 
        {
            std::cerr << "Error: Particle " << i << " is not initialized." << std::endl;
            return false;
        }
    }

    //if everything is ok, write the info file
    dataManager->writeInfoFile(fixedStep, fixedTimeSteps, numberOfParticles);

    auto startTimeCalc = std::chrono::high_resolution_clock::now(); // Startzeit für die Berechnungen
    std::shared_ptr<Tree> tree = std::make_shared<Tree>(this);
    //build the tree
    tree->buildTree();
    std::cout << "\nInitial tree size: " << std::fixed << std::scientific << std::setprecision(1) << tree->root->radius <<"m"<< std::endl;

    visualDensityRadius = tree->root->radius / 500;
    //calculate the visualDensity, just for visualization
    tree->calcVisualDensity();
    //calculate the gas density for SPH
    tree->calcGasDensity();
    //the first time after the temprature is set and rho is calculated
    initGasParticleProperties(tree);

    auto endTimeCalc = std::chrono::high_resolution_clock::now(); // Endzeit für die Berechnungen
    // Berechne die Zeit, die für die Berechnungen benötigt wurde
    auto calcDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeCalc - startTimeCalc).count() /1000.0;


    // Startzeit für das Speichern der Daten
    auto startTimeSave = std::chrono::high_resolution_clock::now();

    // Initial force calculation
    tree->calculateForces();

    //save the particles data
    dataManager->saveData(particles, 0);
    auto endTimeSave = std::chrono::high_resolution_clock::now(); // Endzeit für das Speichern der Daten
    auto saveDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeSave - startTimeSave).count() /1000.0;

    std::cout << std::fixed << "\ncalculation takes " << int(100* (calcDuration / (calcDuration + saveDuration))) << " % of simualtion time" << std::endl;
    std::cout << std::fixed << "data saving takes " << int(100* (saveDuration / (calcDuration + saveDuration))) << " % of simualtion time" << std::endl;

    //calculate the storage size of the data
    //mb per 10000 particles = 1,79 MB
    double perParticle = 1.797 / 10000.0;
    double storageSize = perParticle * numberOfParticles * fixedTimeSteps;
    if(storageSize < 1000)
    {
        std::cout << "Storage size of the data: " << storageSize << " MB" << std::endl;
    }
    else
    {
        std::cout << "Storage size of the data: " << storageSize / 1000 << " GB" << std::endl;
    }

    return true;
}

void Simulation::run()
{
    globalTime = 0.0;
    double nextSaveTime = fixedStep;

    // Set the next integration time for each particle to 0.0 to ensure that the force is calculated in the first iteration
    for (int i = 0; i < numberOfParticles; i++)
    {
        if (!particles[i]) 
        {
            std::cerr << "Error: Particle " << i << " is not initialized." << std::endl;
            continue;
        }
        particles[i]->nextIntegrationTime = 0.0;
    }

    // Initialize particles' time steps and next integration times
    for (int i = 0; i < numberOfParticles; i++)
    {
        if (!particles[i]) 
        {
            std::cerr << "Error: Particle " << i << " is not initialized." << std::endl;
            continue;
        }
        double accelMag = particles[i]->acceleration.length();
        if (accelMag > 0) {
            double timeStep = eta * std::sqrt(e0 / accelMag);
            particles[i]->timeStep = std::clamp(timeStep, minTimeStep, maxTimeStep);
            particles[i]->timeStep = std::max(std::pow(2, std::floor(std::log2(particles[i]->timeStep))), minTimeStep);
            particles[i]->nextIntegrationTime = globalTime + particles[i]->timeStep;
        } else {
            particles[i]->timeStep = minTimeStep;
            particles[i]->nextIntegrationTime = globalTime + particles[i]->timeStep;
        }
    }

    // Main simulation loop
    while (globalTime < endTime)
    {
        // Determine the next integration time for each particle
        for (int i = 0; i < numberOfParticles; i++)
        {
            if (!particles[i]) 
            {
                std::cerr << "Error: Particle " << i << " is not initialized." << std::endl;
                continue;
            }
            if (globalTime >= particles[i]->nextIntegrationTime)
            {
                double accelMag = particles[i]->acceleration.length();
                if (accelMag > 0) {
                    double timeStep = eta * std::sqrt(e0 / accelMag);
                    particles[i]->timeStep = std::clamp(timeStep, minTimeStep, maxTimeStep);
                    particles[i]->timeStep = std::max(std::pow(2, std::floor(std::log2(particles[i]->timeStep))), minTimeStep);
                    particles[i]->nextIntegrationTime = globalTime + particles[i]->timeStep;
                } else {
                    particles[i]->timeStep = minTimeStep;
                    particles[i]->nextIntegrationTime = globalTime + particles[i]->timeStep;
                }
            }
        }

        // Find the smallest next integration time among all particles
        double minIntegrationTime = std::numeric_limits<double>::max();
        for (int i = 0; i < numberOfParticles; i++)
        {
            if (!particles[i]) 
            {
                std::cerr << "Error: Particle " << i << " is not initialized." << std::endl;
                continue;
            }
            if (particles[i]->nextIntegrationTime < minIntegrationTime)
            {
                minIntegrationTime = particles[i]->nextIntegrationTime;
            }
        }

        // Advance global time by the smallest integration time
        globalTime = minIntegrationTime;

        // Update positions and velocities using the KDK Leapfrog scheme for particles due to be integrated
        for (int i = 0; i < numberOfParticles; i++)
        {
            if (!particles[i]) 
            {
                std::cerr << "Error: Particle " << i << " is not initialized." << std::endl;
                continue;
            }
            if (globalTime == particles[i]->nextIntegrationTime)
            {
                timeIntegration->Kick(particles[i], particles[i]->timeStep);
                timeIntegration->Drift(particles[i], particles[i]->timeStep);
            }
        }

        std::shared_ptr<Tree> tree = std::make_shared<Tree>(this);
        
        // Build the octree
        tree->buildTree();

        // Calculate the visual density, just for visualization
        tree->calcVisualDensity();

        // Calculate the gas density for SPH
        tree->calcGasDensity();
        updateGasParticleProperties(tree);

        // Recalculate forces
        tree->calculateForces();

        // Second kick
        for (int i = 0; i < numberOfParticles; i++)
        {
            if (!particles[i]) 
            {
                std::cerr << "Error: Particle " << i << " is not initialized." << std::endl;
                continue;
            }
            if (globalTime == particles[i]->nextIntegrationTime)
            {
                if(particles[i]->type == 2)
                {
                    // Integrate the entropy
                    timeIntegration->Ueuler(particles[i], particles[i]->timeStep);
                }
                timeIntegration->Kick(particles[i], particles[i]->timeStep);
                // Schedule the next integration time for this particle
                particles[i]->nextIntegrationTime += particles[i]->timeStep;
            }
        }

        // Save data at regular intervals defined by fixedStep
        if (globalTime >= nextSaveTime)
        {
            dataManager->saveData(particles, static_cast<int>(nextSaveTime / fixedStep));
            console->printProgress(static_cast<int>(nextSaveTime / fixedStep), fixedTimeSteps, "");
            nextSaveTime += fixedStep;
        }
    }

    std::cout << "Simulation finished." << std::endl;
}

void Simulation::initGasParticleProperties(std::shared_ptr<Tree> tree)
{
    //update the properties of the gas particles
    for (int i = 0; i < numberOfParticles; i++)
    {
        if(particles[i]->type == 2)
        {
            //calc U from T, u = 1 / (gamma-1) * bk * T / (meanMolWeight* prtn)
            particles[i]->U = (1.0 / (Constants::GAMMA - 1.0)) * Constants::BK * particles[i]->T / (Constants::meanMolWeight * Constants::prtn);
            //std::cout << "Particle " << particles[i]->rho << std::endl;
            //calc P, P = (gamma-1)*u*rho
            particles[i]->P = (Constants::GAMMA - 1.0) * particles[i]->U * particles[i]->rho;
        }
    }

    //calc Median Pressure
    tree->root->calcMedianPressure();
}

void Simulation::updateGasParticleProperties(std::shared_ptr<Tree> tree)
{
    //update the properties of the gas particles
    for (int i = 0; i < numberOfParticles; i++)
    {
        if(particles[i]->type == 2)
        {
            //if(i == 100) std::cout << particles[i]->U << std::endl;
            //calc P, P = (gamma-1)*u*rho
            particles[i]->P = (Constants::GAMMA - 1.0) * particles[i]->U * particles[i]->rho;
        }
    }
    
    //calc Median Pressure
    tree->root->calcMedianPressure();
}

void Simulation::calculateForcesWithoutOctree(std::shared_ptr<Particle> p)
{
    p->acceleration = vec3(0.0, 0.0, 0.0);
    p->dUdt = 0;

    for (int j = 0; j < numberOfParticles; j++)
    {
        if (p != particles[j])
        {
            vec3 d = particles[j]->position -p->position;
            double r = d.length();
            vec3 newAcceleration = d * (Constants::G * particles[j]->mass / (r * r * r));
            p->acceleration += newAcceleration;
        }
    }
}
