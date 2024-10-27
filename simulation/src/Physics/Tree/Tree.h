#pragma once
#include <iostream>
#include "Node.h"
#include <memory>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include "Simulation.h"

class Simulation;

class Tree
{
public:
    Tree(Simulation* sim){ simulation = sim;};
    ~Tree(){root.reset();};

    Simulation* simulation;

    std::shared_ptr<Node> root;


    void buildTree();
    void calculateForces();
    double calcTreeWidth();
    void calcVisualDensity();
    void calcGasDensity();

    //multithreading
    void calculateForcesWorker();
    std::atomic<int> currentParticleIndex;
    std::mutex mutex;
};
