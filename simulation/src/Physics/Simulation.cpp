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
#include "Log.h"

Simulation::Simulation()
{
    //construct the modules
    timeIntegration = std::make_shared<TimeIntegration>();
    dataManager = std::make_shared<DataManager>("../../output_data/");
    console = std::make_shared<Console>();
    cooling = std::make_shared<Cooling>();
    sfr = std::make_shared<SFR>();
}

Simulation::~Simulation()
{
}

bool Simulation::init()
{
    //load the config file
    if (!dataManager->loadConfig("../Config.ini", this))
    {
        std::cerr << "Error: Could not load the config file." << std::endl;
        return false;
    }
    Log::setOutputDir(dataManager->outputPath + "/logs");
    
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
    
    if(true) 
    {
        Log::startProcess("load IC");
        dataManager->loadICs(particles, this);
    }
    else
    {
        Log::startProcess("generate IC");
        //Halo* halo = new Halo();
        //halo->generateHernquistHalo(0, numberOfParticles, particles);
        //delete halo;
        Disk* disk = new Disk();
        disk->generateDisk(0, numberOfParticles, particles);
        //delete disk;
    }
    

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


//save velocity curve
if (false)
{
    int N = 1000;
    //get the velocity of the particles
    struct data
    {
        double r;
        double v;
        double i;
    };
    std::vector<data> v;
    int step = numberOfParticles / N;
    for (int i = 0; i < numberOfParticles; i += step)
    {
        data d;
        d.r = particles[i]->position.length() / Units::KPC;
        d.v = particles[i]->velocity.length() / Units::KMS;
        d.i = i;
        v.push_back(d);
    }

    //sort v after r
    std::sort(v.begin(), v.end(), [](data a, data b) { return a.r < b.r; });

    //print the velocity of the particles
    for (size_t i = 0; i < v.size(); i++)
    {
        Log::printData("vel_Curve.csv",v[i].r, v[i].v);
    }
}
    //shuffle the particles to get a random distribution
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(particles.begin(), particles.end(), g);

    Log::startProcess("build Tree");
    std::shared_ptr<Tree> tree = std::make_shared<Tree>(this);
    //build the tree
    tree->buildTree();
    std::cout << "\nInitial tree size: " << std::fixed << std::scientific << std::setprecision(1) << tree->root->radius <<"m"<< std::endl;
    
    Log::startProcess("Visual Density");
    visualDensityRadius = tree->root->radius / 10000;
    //calculate the visualDensity, just for visualization
    tree->calcVisualDensity();
    //calculate the gas density for SPH
    Log::startProcess("SPH density");
    tree->calcGasDensity();
    //the first time after the temprature is set and rho is calculated
    Log::startProcess("Update SPH");
    updateGasParticleProperties(tree);

    // Initial force calculation
    Log::startProcess("Force Calculation");
    tree->calculateForces();

    //save the particles data
    Log::startProcess("Save data");
    dataManager->saveData(particles, 0, fixedTimeSteps, numParticlesOutput, fixedStep, endTime, 0.0);
    
    //print the memory size of the data
    double storageSize = fixedTimeSteps;
    if(dataManager->outputFormat == "ag") storageSize *= dataManager->ag_MemorySize;
    if(dataManager->outputFormat == "age") storageSize *= dataManager->age_MemorySize;
    if(dataManager->outputFormat == "agc") storageSize *= dataManager->agc_MemorySize;
    if(dataManager->outputFormat == "hdf5") storageSize *= dataManager->hdf5_MemorySize;
    if(dataManager->outputFormat == "gadget") storageSize *= dataManager->gadget_MemorySize;
    if(storageSize < 1000000000)
    {
        std::cout << std::fixed << std::setprecision(1) << "Storage size of the data: " << storageSize / 1000000 << " MB" << std::endl;
    }
    else
    {
        std::cout << std::fixed << std::setprecision(1) << "Storage size of the data: " << storageSize / 1000000000 << " GB" << std::endl;
    }

    Log::endProcess();
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
        
        double gasMass = 0;
        double totalMass = 0;

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
                    //cooling and star formation
                    if(coolingEnabled)
                    {
                        cooling->coolingRoutine(particles[i]);
                    }
                    if(starFormation)
                    {
                        //calc SFR
                        sfr->sfrRoutine(particles[i]);
                    }
                    if(particles[i]->type == 2)
                    {
                        // Integrate the internal energy
                        timeIntegration->Ueuler(particles[i], particles[i]->timeStep);
                    }
                }

                //aplly hubbles law
                double H0SI = (H0 * Units::KMS) / Units::MPC;
                double scale_factor = exp(H0SI * particles[i]->timeStep);
                particles[i]->position *= scale_factor;

                timeIntegration->Kick(particles[i], particles[i]->timeStep);
                // Schedule the next integration time for this particle
                particles[i]->nextIntegrationTime += particles[i]->timeStep;
            }
            if(particles[i]->type == 2) gasMass += particles[i]->mass;
            totalMass += particles[i]->mass;
        }
        //std::cout << "Gas in the system: " << gasMass / totalMass * 100 << "%" << std::endl;

        // Save data at regular intervals defined by fixedStep
        if (globalTime >= nextSaveTime)
        {
            dataManager->saveData(particles, static_cast<int>(nextSaveTime / fixedStep), fixedTimeSteps, numParticlesOutput, fixedStep, endTime, globalTime);
            console->printProgress(static_cast<int>(nextSaveTime / fixedStep), fixedTimeSteps, "");
            nextSaveTime += fixedStep;
        }
    }

    std::cout << "Simulation finished." << std::endl;
}

void Simulation::updateGasParticleProperties(std::shared_ptr<Tree> tree)
{
    //update the properties of the gas particles
    for (int i = 0; i < numberOfParticles; i++)
    {
        if(particles[i]->type == 2)
        {
            //calc P, P = (gamma-1)*u*rho
            particles[i]->P = (Constants::GAMMA - 1.0) * particles[i]->U * particles[i]->rho;
            //calc T, T = (gamma-1)*u*prtn / (bk)
            particles[i]->T = (Constants::GAMMA - 1.0) * particles[i]->U * Constants::prtn * particles[i]->mu / (Constants::k_b);
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
