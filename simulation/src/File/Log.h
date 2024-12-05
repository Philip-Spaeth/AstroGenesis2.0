#ifndef LOG_H
#define LOG_H

#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <vector>
#include <memory>
#include "Particle.h"



namespace Log
{
    void setOutputDir(const std::string& outputDir);
    extern std::string outputDir;

//logs are saved in the output dir /logs

    // track process time
    void startProcess(const std::string& processName);
    void endProcess();
    extern bool hasStarted;
    extern std::ofstream LogsDir;
    extern std::ofstream proccessFile;
    extern std::chrono::steady_clock::time_point startTimestamp;
    extern std::string currentProcessName;

    //save data to in csv file
    void printData(const std::string& file, const double x, const double y);
    extern std::vector<std::ofstream> dataFiles;

    //save Velocity Curve
    void saveVelocityCurve(std::vector<std::shared_ptr<Particle>> particles, int numberOfParticles);
}
#endif