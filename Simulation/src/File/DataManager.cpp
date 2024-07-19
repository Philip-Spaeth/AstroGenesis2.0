#include "DataManager.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <string.h>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;


DataManager::DataManager(std::string path)
{
    this->path = path;
}

DataManager::~DataManager()
{
}

void DataManager::writeInfoFile(int deltaTime, int timeSteps, int numberOfParticles)
{
    //create the directory if it does not exist
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

    file << deltaTime<< ";" << std::endl;
    file << timeSteps << ";" << std::endl;
    file << numberOfParticles << ";" << std::endl;

    file.close();
}

void DataManager::readInfoFile(int& deltaTime, int& timeSteps, int& numberOfParticles)
{
    std::string filename = this->path + "info.txt";
    std::ifstream
        file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::getline(file, line);
    deltaTime = std::stoi(line);

    std::getline(file, line);
    timeSteps = std::stoi(line);

    std::getline(file, line);
    numberOfParticles = std::stoi(line);

    file.close();
}

void DataManager::saveData(std::vector<std::shared_ptr<Particle>> particles, int timeStep)
{
    //create the directory if it does not exist
    if (!fs::exists(this->path))
    {
        fs::create_directories(this->path);
    }
    std::string filename = this->path + std::to_string(timeStep) + ".bin";
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Reserve space in the file for all particles
    size_t totalSize = particles.size() * (sizeof(vec3) * 2 + sizeof(double) * 5 + sizeof(int));
    file.seekp(totalSize - 1);
    file.write("", 1); // Write a dummy byte to extend the file size

    // Write particles in a single block
    char* buffer = reinterpret_cast<char*>(malloc(totalSize));
    if (buffer) {
        char* ptr = buffer;
        for (const auto& particle : particles) {
            memcpy(ptr, &particle->position, sizeof(vec3)); ptr += sizeof(vec3);
            memcpy(ptr, &particle->velocity, sizeof(vec3)); ptr += sizeof(vec3);
            memcpy(ptr, &particle->mass, sizeof(double)); ptr += sizeof(double);
            memcpy(ptr, &particle->temperature, sizeof(double)); ptr += sizeof(double);
            memcpy(ptr, &particle->pressure, sizeof(double)); ptr += sizeof(double);
            memcpy(ptr, &particle->density, sizeof(double)); ptr += sizeof(double);
            memcpy(ptr, &particle->viscosity, sizeof(double)); ptr += sizeof(double);
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
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    particles.clear(); // Clear existing particles

    while (!file.eof()) {
        auto particle = std::make_shared<Particle>();

        // Read position
        file.read(reinterpret_cast<char*>(&particle->position), sizeof(vec3));
        if (file.gcount() != sizeof(vec3)) {
            break; // Error or end of file
        }

        // Read velocity
        file.read(reinterpret_cast<char*>(&particle->velocity), sizeof(vec3));
        if (file.gcount() != sizeof(vec3)) {
            break; // Error or end of file
        }

        // Read remaining properties
        file.read(reinterpret_cast<char*>(&particle->mass), sizeof(double));
        file.read(reinterpret_cast<char*>(&particle->temperature), sizeof(double));
        file.read(reinterpret_cast<char*>(&particle->pressure), sizeof(double));
        file.read(reinterpret_cast<char*>(&particle->density), sizeof(double));
        file.read(reinterpret_cast<char*>(&particle->viscosity), sizeof(double));
        file.read(reinterpret_cast<char*>(&particle->type), sizeof(int));

        particles.push_back(particle);
    }

    file.close();
}
