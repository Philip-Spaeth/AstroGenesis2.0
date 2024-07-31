#include "Simulation.h"
#include <numeric>
#include "Galaxy.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>


Simulation::Simulation()
{
    //construct the modules
    timeIntegration = std::make_shared<TimeIntegration>();
    dataManager = std::make_shared<DataManager>("../../../Data/test/");

    //write the info file
    dataManager->writeInfoFile(fixedStep, fixedTimeSteps, numberOfParticles);
}

Simulation::~Simulation(){}

bool Simulation::init()
{
    //setting up multitheading
    std::cout << "Number of threads: " << std::thread::hardware_concurrency() <<"\n"<<std::endl;

    //read the template
    //dataManager->readTemplate("Galaxy1.txt", 0, 1250, vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 0.0), particles);
    //dataManager->readTemplate("Galaxy1.txt", 1250, 2500, vec3(5e22, 1.3e22, 0.0), vec3(-1e5, -0.2e5, 0.0), particles);

    //box with 20 particles
    for(int i = 0; i < numberOfParticles; i++)
    {
        particles.push_back(std::make_shared<Particle>());
        particles[i]->position = vec3(random::between(-100, 100), random::between(-100, 100), random::between(-100, 100));
        particles[i]->velocity = vec3(random::between(0, 0), random::between(0, 0), random::between(0, 0));
        particles[i]->mass = 2e9;
    }

    //build the tree
    buildTree();

    //calculate the density
    calcDensity();

    //save the particles data
    dataManager->saveData(particles, 0);

    std::cout << "Start data successfull initialized\n" << std::endl;

    return true;
}

void Simulation::run()
{
    globalTime = 0.0;
    double nextSaveTime = fixedStep;

    //set the next integration time for each particle to 0.0 to ensure that the force is calculated in the first iteration
    for (int i = 0; i < numberOfParticles; i++)
    {
        particles[i]->nextIntegrationTime = 0.0;
    }
    // Initial force calculation
    calculateForces();

    // Initialize particles' time steps and next integration times
    for (int i = 0; i < numberOfParticles; i++)
    {
        double accelMag = particles[i]->acceleration.length();
        double timeStep = eta * std::sqrt(softening / accelMag);
        particles[i]->timeStep = std::clamp(timeStep, minTimeStep, maxTimeStep);

        // Round the time step down to the nearest power of two
        particles[i]->timeStep = std::pow(2, std::floor(std::log2(particles[i]->timeStep)));
        particles[i]->nextIntegrationTime = globalTime + particles[i]->timeStep;
    }

    // Main simulation loop
    while (globalTime < endTime)
    {
        // Determine the next integration time for each particle
        for (int i = 0; i < numberOfParticles; i++)
        {
            if (globalTime >= particles[i]->nextIntegrationTime)
            {
                double accelMag = particles[i]->acceleration.length();
                double timeStep = eta * std::sqrt(softening / accelMag);
                particles[i]->timeStep = std::clamp(timeStep, minTimeStep, maxTimeStep);

                // Round the time step down to the nearest power of two
                particles[i]->timeStep = std::pow(2, std::floor(std::log2(particles[i]->timeStep)));
                particles[i]->nextIntegrationTime = globalTime + particles[i]->timeStep;
            }
        }

        // Find the smallest next integration time among all particles
        double minIntegrationTime = std::numeric_limits<double>::max();
        for (int i = 0; i < numberOfParticles; i++)
        {
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
            if (globalTime == particles[i]->nextIntegrationTime)
            {
                //std::cout << particles[i]->timeStep << std::endl;
                timeIntegration->Kick(particles[i], particles[i]->timeStep);
                timeIntegration->Drift(particles[i], particles[i]->timeStep);
            }
        }

        // Build the octree
        buildTree();

        //calculate the density
        calcDensity();

        // Recalculate forces
        calculateForces();

        // Second kick
        for (int i = 0; i < numberOfParticles; i++)
        {
            if (globalTime == particles[i]->nextIntegrationTime)
            {
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
            root->calculateForce(particles[i], softening, theta);
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

void Simulation::calcDensity()
{
    for (int i = 0; i < numberOfParticles; i++)
    {
        if (auto node = particles[i]->node.lock()) // Convert weak_ptr to shared_ptr for access
        {
            particles[i]->density = node->calcDensity(h);
        }
    }
}

void Simulation::calculateForcesWithoutOctree(std::shared_ptr<Particle> p)
{
    p->acceleration = vec3(0.0, 0.0, 0.0);

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