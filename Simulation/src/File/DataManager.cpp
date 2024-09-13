#include "DataManager.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <chrono>
#include <cstring>
#include "Particle.h"
#include "vec3.h"
#include "Constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <cstdint>

#ifdef _WIN32
#include <windows.h>
#include <intrin.h>
#else
#include <unistd.h>
#endif

using namespace std;
namespace fs = std::filesystem;

DataManager::DataManager(std::string path)
{
  this->path = path;
}

void DataManager::writeInfoFile(double deltaTime, double timeSteps, double numberOfParticles)
{
    if (!fs::exists(this->path))
    {
        fs::create_directories(this->path);
    }
    
    std::string filename = this->path + "info.txt";
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    file << deltaTime << ";" << std::endl;
    file << timeSteps << ";" << std::endl;
    file << numberOfParticles << ";" << std::endl;

    file.close();
}

void DataManager::readInfoFile(double& deltaTime, double& timeSteps, double& numberOfParticles)
{
    std::string filename = this->path + "info.txt";
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Could not open the file!" << std::endl;
        return;
    }

    std::string line;

    // Read deltaTime from the first line
    if (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!(iss >> deltaTime)) {
            std::cerr << "Error reading deltaTime from line: " << line << std::endl;
            return;
        }
    } else {
        std::cerr << "Error reading the first line for deltaTime" << std::endl;
        return;
    }

    // Read timeSteps from the second line
    if (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!(iss >> timeSteps)) {
            std::cerr << "Error reading timeSteps from line: " << line << std::endl;
            return;
        }
    } else {
        std::cerr << "Error reading the second line for timeSteps" << std::endl;
        return;
    }

    // Read numberOfParticles from the third line
    if (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!(iss >> numberOfParticles)) {
            std::cerr << "Error reading numberOfParticles from line: " << line << std::endl;
            return;
        }
    } else {
        std::cerr << "Error reading the third line for numberOfParticles" << std::endl;
        return;
    }

    file.close();
}


void DataManager::saveData(std::vector<std::shared_ptr<Particle>> particles, int timeStep)
{
    if (!fs::exists(this->path))
    {
        fs::create_directories(this->path);
    }
    std::string filename = this->path + std::to_string(timeStep) + ".bin";
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening datafile: " << filename << std::endl;
        return;
    }

    size_t totalSize = particles.size() * (sizeof(vec3) * 2 + sizeof(double) * 5 + sizeof(int));
    file.seekp(totalSize - 1);
    file.write("", 1);

    char* buffer = reinterpret_cast<char*>(malloc(totalSize));
    if (buffer) {
        char* ptr = buffer;
        for (const auto& particle : particles) {
            memcpy(ptr, &particle->position, sizeof(vec3)); ptr += sizeof(vec3);
            memcpy(ptr, &particle->velocity, sizeof(vec3)); ptr += sizeof(vec3);
            memcpy(ptr, &particle->mass, sizeof(double)); ptr += sizeof(double);
            memcpy(ptr, &particle->T, sizeof(double)); ptr += sizeof(double);
            memcpy(ptr, &particle->P, sizeof(double)); ptr += sizeof(double);
            memcpy(ptr, &particle->visualDensity, sizeof(double)); ptr += sizeof(double); // Added visualDensity not real SPH density
            memcpy(ptr, &particle->U, sizeof(double)); ptr += sizeof(double);
            memcpy(ptr, &particle->type, sizeof(int)); ptr += sizeof(int);
        }
        file.write(buffer, totalSize);
        free(buffer);
    }

    file.close();
}

void DataManager::loadData(int timeStep, std::vector<std::shared_ptr<Particle>>& particles)
{
    std::string filename = this->path + std::to_string(timeStep) + ".bin";
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening datafile: " << filename << std::endl;
        return;
    }

    particles.clear();

    while (!file.eof()) {
        auto particle = std::make_shared<Particle>();

        file.read(reinterpret_cast<char*>(&particle->position), sizeof(vec3));
        if (file.gcount() != sizeof(vec3)) {
            break;
        }

        file.read(reinterpret_cast<char*>(&particle->velocity), sizeof(vec3));
        if (file.gcount() != sizeof(vec3)) {
            break;
        }

        file.read(reinterpret_cast<char*>(&particle->mass), sizeof(double));
        file.read(reinterpret_cast<char*>(&particle->T), sizeof(double));
        file.read(reinterpret_cast<char*>(&particle->P), sizeof(double));
        file.read(reinterpret_cast<char*>(&particle->rho), sizeof(double));
        file.read(reinterpret_cast<char*>(&particle->U), sizeof(double));  
        file.read(reinterpret_cast<char*>(&particle->type), sizeof(int));

        particles.push_back(particle);
    }

    file.close();
}