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
    ~DataManager();
    std::string path;

    void writeInfoFile(double deltaTime, double timeSteps, double numberOfParticles);
    void readInfoFile(double& deltaTime, double& timeSteps, double& numberOfParticles);
    void saveData(std::vector<std::shared_ptr<Particle>> particles, int timeStep);
    void loadData(int timeStep, std::vector<std::shared_ptr<Particle>>& particles);

    //templates
    void readTemplate(std::string fileName, int start, int end, vec3 pos, vec3 vel, std::vector<std::shared_ptr<Particle>>& particles);

    void printProgress(double currentStep, double steps, std::string text);

private:
    std::chrono::_V2::system_clock::time_point startTime;
    bool timerStarted = false;
};
