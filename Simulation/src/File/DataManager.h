#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include <string>
#include <vector>
#include <memory>
#include "Particle.h"
#include "vec3.h"
#include <chrono>
#include <iostream>

using namespace std;

class DataManager {
public:
    DataManager(std::string path);
    ~DataManager(){}

    void writeInfoFile(double deltaTime, double timeSteps, double numberOfParticles);
    void readInfoFile(double& deltaTime, double& timeSteps, double& numberOfParticles);
    void saveData(std::vector<std::shared_ptr<Particle>> particles, int timeStep);
    void loadData(int timeStep, std::vector<std::shared_ptr<Particle>>& particles);
    void printProgress(double currentStep, double steps, std::string text);
    void readTemplate(std::string fileName, int start, int end, vec3 pos, vec3 vel, std::vector<std::shared_ptr<Particle>>& particles);

private:
    std::string path;
    bool timerStarted;
    std::chrono::high_resolution_clock::time_point startTime;
};

#endif // DATAMANAGER_H
