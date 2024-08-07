#pragma once

#include <iostream>
#include <memory>
#include <vector>
#include <chrono>
#include "Particle.h"
#include "vec3.h"

class DataManager
{
public:
    DataManager(std::string path);
    ~DataManager(){}
    std::string path;

    //info file with simulation parameters
    void writeInfoFile(double deltaTime, double timeSteps, double numberOfParticles);
    void readInfoFile(double& deltaTime, double& timeSteps, double& numberOfParticles);

//AGF(Astro Genesis Format)
    //save and load data in binary format, AstroGenesis2.0 specific format
    void saveData(std::vector<std::shared_ptr<Particle>> particles, int timeStep);
    //read the same format, used in the Render program
    void loadData(int timeStep, std::vector<std::shared_ptr<Particle>>& particles);

//Text files, read ASCII format, slower than binary  
    void readASCII(std::string fileName, int start, int end, vec3 pos, vec3 vel, std::vector<std::shared_ptr<Particle>>& particles);

//Gadget 2 format
    //save and load data in binary format, Gadget2 specific format
    void readGadget2Snapshot(std::string fileName, std::vector<std::shared_ptr<Particle>>& particles);

    //progress bar
    void printProgress(double currentStep, double steps, std::string text);

private:
    std::chrono::_V2::system_clock::time_point startTime;
    bool timerStarted = false;
};