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
#include <cctype>
#include <locale>

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

bool parseKeyValue(const std::string& line, std::string& key, std::string& value) {
    std::size_t pos = line.find('=');
    if (pos == std::string::npos) return false;

    key = line.substr(0, pos);
    value = line.substr(pos + 1);

    // Entferne Leerzeichen
    key.erase(key.find_last_not_of(" \n\r\t") + 1);
    value.erase(0, value.find_first_not_of(" \n\r\t"));

    return !key.empty() && !value.empty();
}


// Hilfsfunktion zum Entfernen von Leerzeichen am Anfang und Ende eines Strings
std::string trim(const std::string& str) {
    std::string trimmed = str;
    trimmed.erase(trimmed.begin(), std::find_if(trimmed.begin(), trimmed.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
    trimmed.erase(std::find_if(trimmed.rbegin(), trimmed.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), trimmed.end());
    return trimmed;
}

bool DataManager::loadConfig(const std::string& filename, Simulation* simulation) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Fehler beim Öffnen der Konfigurationsdatei: " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Ignoriere leere Zeilen und Kommentare
        if (line.empty() || line[0] == '#') continue;

        std::string key, value;
        if (parseKeyValue(line, key, value)) {
            // Entferne führende und nachfolgende Leerzeichen
            key = trim(key);
            value = trim(value);
            try {
                if (key == "numberOfParticles") simulation->numberOfParticles = std::stod(value);
                else if (key == "eta") simulation->eta = std::stod(value);
                else if (key == "maxTimeStep") simulation->maxTimeStep = std::stod(value);
                else if (key == "minTimeStep") simulation->minTimeStep = std::stod(value);
                else if (key == "globalTime") simulation->globalTime = std::stod(value);
                else if (key == "endTime") simulation->endTime = std::stod(value);
                else if (key == "fixedTimeSteps") simulation->fixedTimeSteps = std::stod(value);
                else if (key == "e0") simulation->e0 = std::stod(value);
                else if (key == "massInH") simulation->massInH = std::stod(value);
                else if (key == "H0") simulation->H0 = std::stod(value);
                else if (key == "theta") simulation->theta = std::stod(value);
                else if (key == "filePath") simulation->ICFileName = value;
                else if (key == "format") simulation->ICFileFormat = value;
                else if (key == "outputFolderName") path += (value + "/");
                else {
                    std::cerr << "unknown key: " << key << std::endl;
                    return false;
                }
            } catch (const std::invalid_argument& e) {
                std::cerr << "Ungültiger Wert für " << key << ": " << value << std::endl;
                return false;
            } catch (const std::out_of_range& e) {
                std::cerr << "Wert für " << key << " ist außerhalb des gültigen Bereichs: " << value << std::endl;
                return false;
            }
        }
    }

    return true;
}