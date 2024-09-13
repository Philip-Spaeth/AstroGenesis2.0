#include "Simulation.h"
#include <numeric>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iomanip>


Simulation::Simulation()
{
    //construct the modules
    timeIntegration = std::make_shared<TimeIntegration>();
    dataManager = std::make_shared<DataManager>("../../output_data/test/");

    //write the info file
    dataManager->writeInfoFile(fixedStep, fixedTimeSteps, numberOfParticles);
}

Simulation::~Simulation(){}

bool Simulation::init()
{   
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
    DataManager::printSystemInfo();


    dataManager->readGadget2Snapshot("galaxy/Galaxy_Test.gal", particles);
    //dataManager->readASCII("Example_Galaxy_1_ASCII_2500p_from_KlausDolag/Galaxy1.txt", 0, 1250, vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 0.0), particles);
    //dataManager->readASCII("Example_Galaxy_1_ASCII_2500p_from_KlausDolag/Galaxy1.txt", 1250, 2500, vec3(5e22, 1.3e22, 0.0), vec3(-1e5, -0.2e5, 0.0), particles);

    //dataManager->readGadget2Snapshot("Model-M1-G2/snap_000", particles);

    if(numberOfParticles != particles.size())
    {
        std::cerr << "Error: Number of particles in the simulation does not match the number of particles in the data file." << std::endl;
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

    //build the tree
    buildTree();

    //calculate the visualDensity, just for visualization
    calcVisualDensity();
    //calculate the gas density for SPH
    calcGasDensity();
    //the first time after the temprature is set and rho is calculated
    initGasParticleProperties();

    // Initial force calculation
    calculateForces();

    //save the particles data
    dataManager->saveData(particles, 0);

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

        // Build the octree
        buildTree();

        // Calculate the visual density, just for visualization
        calcVisualDensity();

        // Calculate the gas density for SPH
        calcGasDensity();
        updateGasParticleProperties();

        // Recalculate forces
        calculateForces();

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
            dataManager->printProgress(static_cast<int>(nextSaveTime / fixedStep), fixedTimeSteps, "");
            nextSaveTime += fixedStep;
        }
    }
}


void Simulation::buildTree()
{
    //create the root node
    root = std::make_shared<Node>();
    //setup the root node
    root->position = vec3(0.0, 0.0, 0.0);
    root->radius = calcTreeWidth();
    root->depth = 0;

    //insert the particles in the tree
    for (int i = 0; i < numberOfParticles; i++)
    {
        root->insert(particles[i]);
    }
}

void Simulation::calculateForces() 
{
    const int numThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    currentParticleIndex.store(0); // Initialisiere den Atom-Index

    for (int i = 0; i < numThreads; ++i) {
        threads.push_back(std::thread(&Simulation::calculateForcesWorker, this));
    }

    for (auto& thread : threads) {
        thread.join();
    }
}

void Simulation::calculateForcesWorker() {
    while (true) {
        int i = currentParticleIndex.fetch_add(1);
        if (i >= numberOfParticles) {
            break; // Beende, wenn alle Partikel bearbeitet sind
        }

        //only calculate the forces for the particles that are due to be integrated
        if (globalTime == particles[i]->nextIntegrationTime)
        {
            // Berechne die Kräfte für das Partikel
            particles[i]->acceleration = vec3(0.0, 0.0, 0.0);
            particles[i]->dUdt = 0;
            if(particles[i]->type == 2)
            {
                particles[i]->dUdt = 0;
            }
            root->calculateGravityForce(particles[i], e0, theta);
        }
    }
}

void Simulation::applyHubbleExpansion()
{
    //convert hubble constant to SI units
    double kmToMeter = 1e-3;
    double mpcToMeter = 3.086e22;
    double hubbleConstantSI = (H0 * kmToMeter) / mpcToMeter;

    for (int i = 0; i < numberOfParticles; i++)
    {
        particles[i]->velocity += particles[i]->position * hubbleConstantSI;
    }
}

double Simulation::calcTreeWidth()
{
    //acceptanceRatio times the standard deviation of the distances
    double acceptanceRatio = 10;

    if (numberOfParticles == 0) return 0;

    std::vector<double> distances(numberOfParticles);
    for (int i = 0; i < numberOfParticles; i++)
    {
        distances[i] = particles[i]->position.length();
    }
    double sum = std::accumulate(distances.begin(), distances.end(), 0.0);
    double mean = sum / numberOfParticles;

    double sq_sum = std::inner_product(distances.begin(), distances.end(), distances.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / numberOfParticles - mean * mean);

    double maxAcceptableDistance = mean + acceptanceRatio * stdev;

    double max = 0;
    for (int i = 0; i < numberOfParticles; i++)
    {
        double distance = distances[i];
        if (distance <= maxAcceptableDistance && distance > max)
        {
            max = distance;
        }
    }
    return max;
}

void Simulation::calcGasDensity()
{
    //set h to 0 for all particles
    for (int i = 0; i < numberOfParticles; i++)
    { 
        if(particles[i] != nullptr)
        {
            if(particles[i]->type == 2)
            {
                particles[i]->h = 0;
            }
        }
    }

    //calculate the h and density for all particles in the tree
    for (int i = 0; i < numberOfParticles; i++)
    {
        if(particles[i] != nullptr)
        {
            if(particles[i]->type == 2)
            {
            
                if (auto node = particles[i]->node.lock()) // Convert weak_ptr to shared_ptr for access
                {
                    if(particles[i]->h == 0)
                    {
                        node->calcGasDensity(massInH);
                    }
                }
            }
        }
    }

    //calculate the median smoothing length and density for all nodes
    root->calcSPHNodeMedians();
}

void Simulation::initGasParticleProperties()
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
    root->calcMedianPressure();
}

void Simulation::updateGasParticleProperties()
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
    root->calcMedianPressure();
}

void Simulation::calcVisualDensity()
{
    //calculate the density for all particles in the tree
    for (int i = 0; i < numberOfParticles; i++)
    {
        if (auto node = particles[i]->node.lock()) // Convert weak_ptr to shared_ptr for access
        {
            node->calcVisualDensity(visualDensityRadius);
        }
    }
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
