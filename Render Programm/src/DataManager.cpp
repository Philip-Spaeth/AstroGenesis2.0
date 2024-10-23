#include "DataManager.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <iomanip>

#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;
namespace fs = std::filesystem;

DataManager::DataManager(std::string path)
{
    this->path = path;
}

DataManager::~DataManager()
{
}


void DataManager::readInfoFile(double& deltaTime, double& timeSteps, double& numberOfParticles)
{
    std::string line;
    std::string filename = this->path + "info.txt";
    std::ifstream file(filename);

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
    
    //read the outputDataFormat from the 4th line
    if (std::getline(file, line)) {
        outputDataFormat = line;
    } else {
        std::cerr << "Error reading the fourth line for outputDataFormat" << std::endl;
        return;
    }
    //romve the last character from the outputDataFormat
    outputDataFormat.pop_back();


    file.close();
}

void DataManager::loadData(int timeStep, std::vector<std::shared_ptr<Particle>>& particles)
{
    // Dateiname basierend auf dem Zeitschritt
    std::string filename = this->path + std::to_string(timeStep) + ".bin";
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening datafile: " << filename << std::endl;
        return;
    }

    particles.clear();

    if (outputDataFormat == "AGF")
    {
        // AGFE Format: Erweiterte Version (bestehend aus Position, Velocity, Mass, T, P, visualDensity, U, type)
        size_t recordSize = sizeof(vec3) * 2 + sizeof(double) * 3 + sizeof(int);

        // Bestimme die Größe der Datei
        file.seekg(0, std::ios::end);
        std::streampos fileSize = file.tellg();
        file.seekg(0, std::ios::beg);

        if (fileSize % recordSize != 0) {
            std::cerr << "Invalid AGF file size: " << filename << std::endl;
            file.close();
            return;
        }

        size_t numParticles = fileSize / recordSize;
        particles.reserve(numParticles);

        // Optional: Lade die gesamte Datei in einen Puffer (schneller für große Dateien)
        std::vector<char> buffer(fileSize);
        file.read(buffer.data(), fileSize);
        if (file.gcount() != fileSize) {
            std::cerr << "Error reading datafile: " << filename << std::endl;
            file.close();
            return;
        }

        const char* ptr = buffer.data();
        for (size_t i = 0; i < numParticles; ++i) {
            auto particle = std::make_shared<Particle>();

            // Position: vec3 (3 doubles)
            memcpy(&particle->position, ptr, sizeof(vec3)); ptr += sizeof(vec3);

            // Velocity: vec3 (3 doubles) - wird nicht im AGFC geladen, aber hier für AGFE
            memcpy(&particle->velocity, ptr, sizeof(vec3)); ptr += sizeof(vec3);

            // Mass: double
            memcpy(&particle->mass, ptr, sizeof(double)); ptr += sizeof(double);

            // T: double
            memcpy(&particle->temperature, ptr, sizeof(double)); ptr += sizeof(double);

            // visualDensity: double
            memcpy(&particle->density, ptr, sizeof(double)); ptr += sizeof(double);

            // type: int
            memcpy(&particle->type, ptr, sizeof(int)); ptr += sizeof(int);

            particles.push_back(particle);
        }
    }
    else if (outputDataFormat == "AGFC")
    {
        // AGFC Format: Kompaktes Format für Rendering
        // Speichert nur position (vec3 als 3 floats), visualDensity (float) und type (int)

        // Größe eines einzelnen AGFC-Records: 3 floats (Position) + 1 float (visualDensity) + 1 int (type)
        const size_t recordSize = sizeof(float) * 3 + sizeof(float) + sizeof(int);

        // Bestimme die Größe der Datei
        file.seekg(0, std::ios::end);
        std::streampos fileSize = file.tellg();
        file.seekg(0, std::ios::beg);

        if (fileSize % recordSize != 0) {
            std::cerr << "Invalid AGFC file size: " << filename << std::endl;
            file.close();
            return;
        }

        size_t numParticles = fileSize / recordSize;
        particles.reserve(numParticles);

        // Optional: Lade die gesamte Datei in einen Puffer (schneller für große Dateien)
        std::vector<char> buffer(fileSize);
        file.read(buffer.data(), fileSize);
        if (file.gcount() != fileSize) {
            std::cerr << "Error reading datafile: " << filename << std::endl;
            file.close();
            return;
        }

        const char* ptr = buffer.data();
        for (size_t i = 0; i < numParticles; ++i) {
            auto particle = std::make_shared<Particle>();

            // Position: 3 floats -> konvertiert zu doubles
            float posX_f, posY_f, posZ_f;
            memcpy(&posX_f, ptr, sizeof(float)); ptr += sizeof(float);
            memcpy(&posY_f, ptr, sizeof(float)); ptr += sizeof(float);
            memcpy(&posZ_f, ptr, sizeof(float)); ptr += sizeof(float);
            particle->mass = 1.0;
            particle->position.x = posX_f;
            particle->position.y = posY_f;
            particle->position.z = posZ_f;
            // visualDensity: float -> double
            float visualDensity_f;
            memcpy(&visualDensity_f, ptr, sizeof(float)); ptr += sizeof(float);
            particle->density = static_cast<double>(visualDensity_f);

            // type: int
            memcpy(&particle->type, ptr, sizeof(int)); ptr += sizeof(int);

            particles.push_back(particle);
        }
    }
    else if (outputDataFormat == "AGFE")
    {
        // AGFE Format: Erweiterte Version (bestehend aus Position, Velocity, Mass, T, P, visualDensity, U, type)
        size_t recordSize = sizeof(vec3) * 2 + sizeof(double) * 5 + sizeof(int);

        // Bestimme die Größe der Datei
        file.seekg(0, std::ios::end);
        std::streampos fileSize = file.tellg();
        file.seekg(0, std::ios::beg);

        if (fileSize % recordSize != 0) {
            std::cerr << "Invalid AGFE file size: " << filename << std::endl;
            file.close();
            return;
        }

        size_t numParticles = fileSize / recordSize;
        particles.reserve(numParticles);

        // Optional: Lade die gesamte Datei in einen Puffer (schneller für große Dateien)
        std::vector<char> buffer(fileSize);
        file.read(buffer.data(), fileSize);
        if (file.gcount() != fileSize) {
            std::cerr << "Error reading datafile: " << filename << std::endl;
            file.close();
            return;
        }

        const char* ptr = buffer.data();
        for (size_t i = 0; i < numParticles; ++i) {
            auto particle = std::make_shared<Particle>();

            // Position: vec3 (3 doubles)
            memcpy(&particle->position, ptr, sizeof(vec3)); ptr += sizeof(vec3);

            // Velocity: vec3 (3 doubles) - wird nicht im AGFC geladen, aber hier für AGFE
            memcpy(&particle->velocity, ptr, sizeof(vec3)); ptr += sizeof(vec3);

            // Mass: double
            memcpy(&particle->mass, ptr, sizeof(double)); ptr += sizeof(double);

            // T: double
            memcpy(&particle->temperature, ptr, sizeof(double)); ptr += sizeof(double);

            // P: double
            memcpy(&particle->pressure, ptr, sizeof(double)); ptr += sizeof(double);

            // visualDensity: double
            memcpy(&particle->density, ptr, sizeof(double)); ptr += sizeof(double);

            // U: double
            memcpy(&particle->internalEnergy, ptr, sizeof(double)); ptr += sizeof(double);

            // type: int
            memcpy(&particle->type, ptr, sizeof(int)); ptr += sizeof(int);

            particles.push_back(particle);
        }
    }
    else
    {
        std::cerr << "Unknown output data format: " << outputDataFormat << std::endl;
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
