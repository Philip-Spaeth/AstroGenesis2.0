#pragma once

#include <iostream>
#include <memory>
#include <vector>
#include <chrono>
#include "Particle.h"
#include "vec3.h"
#include "Simulation.h"

class Simulation;

class DataManager
{
public:
    DataManager(std::string path);
    ~DataManager(){}

    std::string inputPath;
    std::string inputFormat;

    //path to the folder where the simulation data is saved
    std::string path;
    std::string outputDataFormat;

//read the Config.ini file
    bool loadConfig(const std::string& filename, Simulation* simulation);

    bool loadICs(Simulation* simulation);

    //info file with simulation parameters
    //write when saving simulation data
    void writeInfoFile(double deltaTime, double timeSteps, double numberOfParticles);

//AGF(Astro Genesis Format)
    //save and load data in binary format, AstroGenesis2.0 specific format
    void saveData(std::vector<std::shared_ptr<Particle>> particles, int timeStep);

//for reading Gadget2 snapshoot format oder ASCII see the ICDataReader


    //Data Size per 1000 particles in bytes
    double AGF_MemorySize = 77824;
    double AGFC_MemorySize = 20480;
    double AGFE_MemorySize = 94208;
    double AGFH_MemorySize = 13000; // not implemented yet
};