#pragma once

#include <iostream>
#include "Particle.h"
#include <memory>
#include <vector>
#include <fstream>
#include <string>
#include "vec3.h"


//basic protocol for data management
//all in SI units
//position(vec3)s, velocity(vec3), mass temperature, pressure, density, viscosity, type(1 = star,2 = gas,3 = dark matter)

class DataManager
{
public:
    std::string path;

    void writeInfoFile(int deltaTime, int timeSteps, int numberOfParticles);
    void readInfoFile(int& deltaTime, int& timeSteps, int& numberOfParticles);
    void saveData(std::vector<std::shared_ptr<Particle>> particles, int timeStep);
    void loadData(int timeStep, std::vector<std::shared_ptr<Particle>>& particles);
};