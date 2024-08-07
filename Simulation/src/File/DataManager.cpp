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
#include <cstring>  // Fügt den Header für memcpy hinzu
#include "Particle.h"
#include "vec3.h"
#include "Constants.h"

#ifdef _WIN32
#include <windows.h>
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

void DataManager::readGadget2Snapshot(std::string fileName, std::vector<std::shared_ptr<Particle>>& particles)
{
    std::string fileDir = "input_data/";
    std::filesystem::path filePath = "../..";
    filePath = filePath / fileDir / fileName;

    std::ifstream file(filePath, std::ios::binary);
    if (!file) {
        std::cerr << "Could not open the data file: " << filePath << std::endl;
        return;
    }

    std::cout << "Reading file: " << filePath << std::endl;

    // Anzahl der Partikel berechnen
    uint32_t numParticles = 3000;
    std::cout << "Number of particles: " << numParticles << std::endl;

    // Speicher für Partikeldaten reservieren
    particles.resize(numParticles);
    for (int i = 0; i < numParticles; ++i) 
    {
        particles[i] = std::make_shared<Particle>();
    }

    // Lese die Partikelpositionen
    for (uint32_t i = 0; i < numParticles; ++i) 
    {
        float pos[3];
        file.read(reinterpret_cast<char*>(&pos), sizeof(pos));
        //transform to internal units
        particles[i]->position = vec3(pos[0], pos[1], pos[2]);
        particles[i]->position *= Constants::PC;
    }

    // Lese die Partikelgeschwindigkeiten
    for (uint32_t i = 0; i < numParticles; ++i) 
    {
        float vel[3];
        file.read(reinterpret_cast<char*>(&vel), sizeof(vel));
        //transform to internal units
        particles[i]->velocity = vec3(vel[0], vel[1], vel[2]);
        particles[i]->velocity *= 1000;
    }

    // Lese die Partikel-Massen
    for (uint32_t i = 0; i < numParticles; ++i) 
    {
        float mass;
        file.read(reinterpret_cast<char*>(&mass), sizeof(mass));
        //transform to internal units
        particles[i]->mass = 1;
    }

    // Lese die Partikel-IDs
    for (uint32_t i = 0; i < numParticles; ++i) {
        uint64_t id;
        file.read(reinterpret_cast<char*>(&id), sizeof(id));
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
            memcpy(ptr, &particle->A, sizeof(double)); ptr += sizeof(double);
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
        file.read(reinterpret_cast<char*>(&particle->A), sizeof(double));  
        file.read(reinterpret_cast<char*>(&particle->type), sizeof(int));

        particles.push_back(particle);
    }

    file.close();
}

#ifdef _WIN32
void setConsoleColor(WORD color) {
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hConsole, color);
}
#else
void setConsoleColor(int color) {
    // No color setting for Linux in this version
}
#endif
void DataManager::printProgress(double currentStep, double steps, std::string text) 
{
    static const int barWidth = 70;

    if (!timerStarted) {
        std::cout << std::endl;
        startTime = std::chrono::high_resolution_clock::now();
        timerStarted = true;
    }

    double progress = (currentStep + 1) / steps;
    int pos = static_cast<int>(barWidth * progress);
    
    // Set color based on progress
#ifdef _WIN32
    WORD progressColor = (currentStep < steps - 1) ? FOREGROUND_RED : FOREGROUND_GREEN;
#endif

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
#ifdef _WIN32
        if (i < pos) {
            setConsoleColor(progressColor);
            std::cout << "=";
        }
        else if (i == pos && currentStep < steps - 1) {
            setConsoleColor(FOREGROUND_RED);
            std::cout << ">";
        }
        else {
            setConsoleColor(FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED);
            std::cout << " ";
        }
#else
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
#endif
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

#ifdef _WIN32
    setConsoleColor(FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED);
#endif
    std::cout << "] " << int(progress * 100.0) << " %  Estimated remaining time: " << timeleft << unit << "  "<< text << "       " << "\r";

    if (currentStep == steps - 1) {
        std::cout << std::endl;
        std::cout << std::endl;
    }
}

void DataManager::readASCII(std::string fileName, int start, int end, vec3 pos, vec3 vel, std::vector<std::shared_ptr<Particle>>& particles)
{
    fileName = "Dolag_Example_Galaxy_1_ASCII/Galaxy1.txt";

    std::string fileDir = "input_data/";
    std::filesystem::path filePath = fileDir;
    filePath = "../.." / filePath;

    filePath = filePath / fileName;


    // Ensure the particles vector is large enough to hold the new particles
    if (particles.size() < static_cast<size_t>(end)) {
        particles.resize(end);
    }

    int particleIndex = start;

    for (int i = start; i < end; i += 1250) {
        std::ifstream file(filePath);
        if (!file) {
            std::cout << "Reading file: " << filePath << std::endl;
            std::cerr << "Could not open the data file!" << std::endl;
            return;
        }

        std::string line;
        int currentIndex = 0;

        while (std::getline(file, line) && particleIndex < end) 
        {
            std::istringstream iss(line);
            vec3 position, velocity;
            double mass;

            // Assuming the file format is: position3d (3 values), velocity3d (3 values), mass (1 value)
            if (!(iss >> position.x >> position.y >> position.z 
                    >> velocity.x >> velocity.y >> velocity.z 
                    >> mass)) {
                std::cerr << "Error parsing line: " << line << std::endl;
                continue;
            }
            //// depends on the units of the template file
            // Convert the units: data units: kpc, km/s, 1e10 Msun -> internal units: m, m/s, kg
            position *= 3.086e19;
            velocity *= 1e3;
            mass *= 1e10 * 1.989e30;

            // Add the offset
            position += pos;
            velocity += vel;

            // Insert the particle at the correct index
            particles[particleIndex] = std::make_shared<Particle>();

            //set the particle properties
            particles[particleIndex]->position = position;
            particles[particleIndex]->velocity = velocity;
            particles[particleIndex]->mass = mass;
            particles[particleIndex]->type = 2;
            particles[particleIndex]->T = 1e27;
            
            particleIndex++;
            currentIndex++;
        }

        file.close();
    }

    std::cout << "Initial Condition: " << fileName << " with particles from index " << start << " to " << (particleIndex - 1) << std::endl;
}