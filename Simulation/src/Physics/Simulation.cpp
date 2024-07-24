#include "Simulation.h"
#include <numeric>

Simulation::Simulation()
{
    //construct the modules
    timeIntegration = std::make_shared<TimeIntegration>();
    dataManager = std::make_shared<DataManager>("../../../Data/Test/");

    //write the info file
    dataManager->writeInfoFile(deltaTime, timeSteps, numberOfParticles);
}

Simulation::~Simulation(){}

bool Simulation::init()
{
    //setting up multitheading
    std::cout << "Number of threads: " << std::thread::hardware_concurrency() << std::endl;

    random::setRandomSeed(8486368);

    for(int i = 0; i < numberOfParticles; i++)
    {
        particles.push_back(std::make_shared<Particle>(vec3(random::between(-10,10), random::between(-10,10), random::between(-10,10)), vec3(random::between(-0.1,0.1), random::between(-0.1,0.1), random::between(-0.1,0.1)), vec3(0.0, 0.0, 0.0), 10000000));
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
    //calculate the gravitational acceleration for each particle
    for (int t = 0; t < timeSteps; t++)
    {   
        //build the tree
        buildTree();

        //calculate the density
        calcDensity();

        //calculate the forces
        calculateForces();

        for (int i = 0; i < numberOfParticles; i++)
        {
            //update the particle position and velocity
            timeIntegration->Euler(particles[i], deltaTime);
        }

        //save the particles data
        dataManager->saveData(particles, t + 1);
        dataManager->printProgress(t, timeSteps, "");
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

void Simulation::calculateForces() {
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
        // Berechne die Kräfte für das Partikel
        particles[i]->acceleration = vec3(0.0, 0.0, 0.0);
        root->calculateForce(particles[i], softening, theta);
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
        if (particles[i]->node)
        {
            particles[i]->density = particles[i]->node->calcDensity(h);
        }
    }
}