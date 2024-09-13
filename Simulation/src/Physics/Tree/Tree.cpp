#include "Tree.h"
#include <numeric>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iomanip>


void Tree::buildTree()
{
    //create the root node
    root = std::make_shared<Node>();
    //setup the root node
    root->position = vec3(0.0, 0.0, 0.0);
    root->radius = calcTreeWidth();
    root->depth = 0;

    //insert the particles in the tree
    for (int i = 0; i < simulation->numberOfParticles; i++)
    {
        root->insert(simulation->particles[i]);
    }
}

void Tree::calculateForces() 
{
    const int numThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    currentParticleIndex.store(0); // Initialisiere den Atom-Index

    for (int i = 0; i < numThreads; ++i) {
        threads.push_back(std::thread(&Tree::calculateForcesWorker, this));
    }

    for (auto& thread : threads) {
        thread.join();
    }
}

void Tree::calculateForcesWorker() {
    while (true) {
        int i = currentParticleIndex.fetch_add(1);
        if (i >= simulation->numberOfParticles) {
            break; // Beende, wenn alle Partikel bearbeitet sind
        }

        //only calculate the forces for the particles that are due to be integrated
        if (simulation->globalTime == simulation->particles[i]->nextIntegrationTime)
        {
            // Berechne die Kräfte für das Partikel
            simulation->particles[i]->acceleration = vec3(0.0, 0.0, 0.0);
            simulation->particles[i]->dUdt = 0;
            if(simulation->particles[i]->type == 2)
            {
                simulation->particles[i]->dUdt = 0;
            }
            root->calculateGravityForce(simulation->particles[i], simulation->e0, simulation->theta);
        }
    }
}


double Tree::calcTreeWidth()
{
    //acceptanceRatio times the standard deviation of the distances
    double acceptanceRatio = 10;

    if (simulation->numberOfParticles == 0) return 0;

    std::vector<double> distances(simulation->numberOfParticles);
    for (int i = 0; i < simulation->numberOfParticles; i++)
    {
        distances[i] = simulation->particles[i]->position.length();
    }
    double sum = std::accumulate(distances.begin(), distances.end(), 0.0);
    double mean = sum / simulation->numberOfParticles;

    double sq_sum = std::inner_product(distances.begin(), distances.end(), distances.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / simulation->numberOfParticles - mean * mean);

    double maxAcceptableDistance = mean + acceptanceRatio * stdev;

    double max = 0;
    for (int i = 0; i < simulation->numberOfParticles; i++)
    {
        double distance = distances[i];
        if (distance <= maxAcceptableDistance && distance > max)
        {
            max = distance;
        }
    }
    return max;
}

void Tree::calcGasDensity()
{
    //set h to 0 for all particles
    for (int i = 0; i < simulation->numberOfParticles; i++)
    { 
        if(simulation->particles[i] != nullptr)
        {
            if(simulation->particles[i]->type == 2)
            {
                simulation->particles[i]->h = 0;
            }
        }
    }

    //calculate the h and density for all particles in the tree
    for (int i = 0; i < simulation->numberOfParticles; i++)
    {
        if(simulation->particles[i] != nullptr)
        {
            if(simulation->particles[i]->type == 2)
            {
            
                if (auto node = simulation->particles[i]->node.lock()) // Convert weak_ptr to shared_ptr for access
                {
                    if(simulation->particles[i]->h == 0)
                    {
                        node->calcGasDensity(simulation->massInH);
                    }
                }
            }
        }
    }

    //calculate the median smoothing length and density for all nodes
    root->calcSPHNodeMedians();
}


void Tree::calcVisualDensity()
{
    //calculate the density for all particles in the tree
    for (int i = 0; i < simulation->numberOfParticles; i++)
    {
        if (auto node = simulation->particles[i]->node.lock()) // Convert weak_ptr to shared_ptr for access
        {
            node->calcVisualDensity(simulation->visualDensityRadius);
        }
    }
}