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
    std::string path;

    void writeInfoFile(int deltaTime, int timeSteps, int numberOfParticles);
    void readInfoFile(int& deltaTime, int& timeSteps, int& numberOfParticles);
    void saveData(std::vector<std::shared_ptr<Particle>> particles, int timeStep);
    void loadData(int timeStep, std::vector<std::shared_ptr<Particle>>& particles);

    void printProgress(double currentStep, double steps);

private:
    std::chrono::_V2::system_clock::time_point startTime;
    bool timerStarted = false;
};
