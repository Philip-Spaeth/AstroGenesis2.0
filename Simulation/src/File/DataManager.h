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
    //path to the folder where the simulation data is saved
    std::string path;

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

//Text files, read ASCII format, slower than binary  
    void readASCII(std::string fileName, int start, int end, vec3 pos, vec3 vel, std::vector<std::shared_ptr<Particle>>& particles);

//Gadget 2 snapshot-format
    //save and load data in binary format, Gadget2 specific format
    void readGadget2Snapshot(std::string fileName, std::vector<std::shared_ptr<Particle>>& particles);


///console output
    //progress bar
    void printProgress(double currentStep, double steps, std::string text);
    //system info
    static void printSystemInfo();

private:
    void readBlock(std::ifstream& file, char* buffer, size_t size);
    std::chrono::_V2::system_clock::time_point startTime;
    bool timerStarted = false;
};