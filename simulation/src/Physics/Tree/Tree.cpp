#include "Tree.h"
#include <numeric>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <omp.h>


// Destructor for Tree
/* Tree::~Tree()
{
    std::cout << "Tree destroyed: " << this << std::endl;
    root.reset();
} */


void Tree::buildTree()
{
    std::cout << "Building tree..." << std::endl;

    //start Time
    auto start = std::chrono::high_resolution_clock::now();

    //create the root node
    root = std::make_shared<Node>();
    //setup the root node
    root->position = vec3(0.0, 0.0, 0.0);
    root->radius = calcTreeWidth();
    root->depth = 0;
    //double totalMass = 0;
    //double gasMass = 0;
    //insert the particles in the tree
    //#pragma omp parallel for
    /* for (int i = 0; i < simulation->numberOfParticles; i++)
    {
        //totalMass += simulation->particles[i]->mass;
        //if(simulation->particles[i]->type == 2) gasMass += simulation->particles[i]->mass;
        root->insert(simulation->particles[i]);
    } */
    //root->insert(simulation->particles);


    int max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);
    omp_set_nested(1);

    std::cout << "Number of threads: " << max_threads << std::endl;
    root->insert(simulation->particles, max_threads);


    // Löschen des gesamten Baums
    /* root.reset(); // Alle Knoten werden freigegeben, sofern keine weiteren shared_ptr existieren
    
    // Optional: Überprüfen, ob der Baum korrekt gelöscht wurde
    if (root == nullptr) {
        std::cout << "Baum erfolgreich gelöscht.\n";
    } else {
        std::cout << "Baum konnte nicht gelöscht werden.\n";
    } */

    /* //#pragma omp parallel
    {
        //#pragma omp single
        {
            //root->insert(simulation->particles, max_threads);
        }
    } */


    //end Time
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Tree built in " << elapsed_seconds.count() << "s" << std::endl;
}

void Tree::calculateForces() 
{
    const int numThreads = std::thread::hardware_concurrency();
    omp_set_num_threads(numThreads);  // Setze die Anzahl der OpenMP-Threads

    // Erstelle eine lokale Kopie von numberOfParticles
    int numParticles = simulation->numberOfParticles;

    #pragma omp parallel for
    for (int i = 0; i < numParticles; ++i) {
        if (simulation->globalTime == simulation->particles[i]->nextIntegrationTime) {
            simulation->particles[i]->acceleration = vec3(0.0, 0.0, 0.0);
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
    const int numParticles = simulation->numberOfParticles; // lokal speichern
    #pragma omp parallel for
    for (int i = 0; i < numParticles; i++) // Schleifenbedingungen sind jetzt mit einem konstanten Wert
    { 
        if(simulation->particles[i] != nullptr)
        {
            if(simulation->particles[i]->type == 2)
            {
                simulation->particles[i]->h = 0;
            }
        }
    }
    for (int i = 0; i < numParticles; i++)
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
}

void Tree::calcVisualDensity()
{
    // Set visualDensity to 0 for all particles
    const int numParticles = simulation->numberOfParticles; // lokal speichern
    #pragma omp parallel for
    for (int i = 0; i < numParticles; i++) // Schleifenbedingungen sind jetzt mit einem konstanten Wert
    { 
        simulation->particles[i]->visualDensity = 0;
    }

    #pragma omp parallel for
    for (int i = 0; i < numParticles; i++) // Schleifenbedingungen sind jetzt mit einem konstanten Wert
    { 
        if(simulation->particles[i]->visualDensity != 0)
        {
            continue;
        }
        if (auto node = simulation->particles[i]->node.lock()) // Convert weak_ptr to shared_ptr for access
        {
            node->calcVisualDensity(simulation->visualDensityRadius);
        }
    }
}