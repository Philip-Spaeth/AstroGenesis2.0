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

    //path to the folder where the simulation data is saved
    std::string path;

//read the Config.ini file
    bool loadConfig(const std::string& filename, Simulation* simulation);

    //info file with simulation parameters
    //write when saving simulation data
    void writeInfoFile(double deltaTime, double timeSteps, double numberOfParticles);
    //read when loading simulation data for the Render program
    void readInfoFile(double& deltaTime, double& timeSteps, double& numberOfParticles);

//AGF(Astro Genesis Format)
    //save and load data in binary format, AstroGenesis2.0 specific format
    void saveData(std::vector<std::shared_ptr<Particle>> particles, int timeStep);
    //read the same format, used in the Render program
    void loadData(int timeStep, std::vector<std::shared_ptr<Particle>>& particles);
    
//for reading Gadget2 snapshoot format oder ASCII see the ICDataReader

};