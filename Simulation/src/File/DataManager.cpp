#include "DataManager.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <iomanip>

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

void DataManager::readInfoFile(int& deltaTime, int& timeSteps, int& numberOfParticles)
{
    std::string filename = this->path + "info.txt";
    std::ifstream file(filename);
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
        file.read(reinterpret_cast<char*>(&particle->temperature), sizeof(double));
        file.read(reinterpret_cast<char*>(&particle->pressure), sizeof(double));
        file.read(reinterpret_cast<char*>(&particle->density), sizeof(double));
        file.read(reinterpret_cast<char*>(&particle->viscosity), sizeof(double));
        file.read(reinterpret_cast<char*>(&particle->type), sizeof(int));

        particles.push_back(particle);
    }

    file.close();
}

void DataManager::printProgress(double currentStep, double steps) {
    static const int barWidth = 70;
    static const int bufferSize = 100; // Adjust buffer size to be large enough for the entire line
    

    if (!timerStarted) {
        std::cout << std::endl;
        startTime = std::chrono::high_resolution_clock::now();
        timerStarted = true;
    }

    std::cout << "[";
    int pos = barWidth * (((currentStep + 1) / steps) * 100) / 100.0;

    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }

    double remainingTime = 0;

    auto currentTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = currentTime - startTime;
    double estimatedTotalTime = (elapsed.count() / (currentStep + 1)) * steps;
    double estimatedRemainingTime = estimatedTotalTime - elapsed.count();
    remainingTime = estimatedRemainingTime;

    std::string unit;
    if (remainingTime < 60) {
        unit = "s";
    } else if (remainingTime < 3600) {
        remainingTime /= 60;
        unit = "min";
    } else {
        remainingTime /= 3600;
        unit = "h";
    }

    std::ostringstream timeStream;
    timeStream << std::fixed << std::setprecision(1) << remainingTime;
    std::string timeleft = timeStream.str();

    std::ostringstream output;
    output << "] " << std::fixed << std::setprecision(1) << (((currentStep + 1) / steps) * 100) << " %"
           << "  Estimated remaining time: " << timeleft << unit;

    // Ensure the entire line is cleared by filling with spaces up to bufferSize
    std::string outputStr = output.str();
    if (outputStr.length() < bufferSize) {
        outputStr.append(bufferSize - outputStr.length(), ' ');
    }

    std::cout << outputStr << "\r";
    std::cout.flush();

    if (currentStep == steps - 1) {
        std::cout << std::endl;
        std::cout << std::endl;
    }
}