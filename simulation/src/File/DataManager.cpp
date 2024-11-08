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
    file << outputDataFormat << ";" << std::endl;

    file.close();
}

// Funktion zum Delta-Encoding einer Sequenz von Doubles
std::vector<double> deltaEncode(const std::vector<double>& data) {
    std::vector<double> deltas;
    if (data.empty()) return deltas;
    deltas.reserve(data.size());
    deltas.push_back(data[0]); // Erster Wert bleibt unverändert
    for (size_t i = 1; i < data.size(); ++i) {
        deltas.push_back(data[i] - data[i - 1]);
    }
    return deltas;
}

// Funktion zum Delta-Decoding einer Sequenz von Doubles
std::vector<double> deltaDecode(const std::vector<double>& deltas) {
    std::vector<double> data;
    if (deltas.empty()) return data;
    data.reserve(deltas.size());
    data.push_back(deltas[0]);
    for (size_t i = 1; i < deltas.size(); ++i) {
        data.push_back(data[i - 1] + deltas[i]);
    }
    return data;
}

// Funktion zum Run-Length Encoding (RLE) einer Sequenz von Bytes
std::vector<std::pair<uint8_t, uint8_t>> runLengthEncode(const std::vector<uint8_t>& data) {
    std::vector<std::pair<uint8_t, uint8_t>> encoded;
    if (data.empty()) return encoded;
    
    uint8_t current = data[0];
    uint8_t count = 1;
    
    for (size_t i = 1; i < data.size(); ++i) {
        if (data[i] == current && count < 255) {
            count++;
        } else {
            encoded.emplace_back(current, count);
            current = data[i];
            count = 1;
        }
    }
    encoded.emplace_back(current, count);
    return encoded;
}

// Funktion zum Run-Length Decoding (RLE) einer Sequenz von Paaren
std::vector<uint8_t> runLengthDecode(const std::vector<std::pair<uint8_t, uint8_t>>& encoded) {
    std::vector<uint8_t> data;
    data.reserve(encoded.size() * 2); // Geschätzte Größe
    for (const auto& pair : encoded) {
        for (uint8_t i = 0; i < pair.second; ++i) {
            data.push_back(pair.first);
        }
    }
    return data;
}

void DataManager::saveData(std::vector<std::shared_ptr<Particle>> particles, int timeStep)
{
    // Sicherstellen, dass der Pfad existiert
    if (!fs::exists(this->path))
    {
        fs::create_directories(this->path);
    }

    // Dateiname basierend auf dem Zeitschritt
    std::string filename = this->path + std::to_string(timeStep) + ".bin";
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening datafile: " << filename << std::endl;
        return;
    }

    if (outputDataFormat == "AGF")
    {
        // AGF Format
        size_t totalSize = particles.size() * (sizeof(vec3) * 2 + sizeof(double) * 3 + sizeof(int));
        
        // Speicher für den Puffer allokieren
        char* buffer = reinterpret_cast<char*>(malloc(totalSize));
        if (buffer) {
            char* ptr = buffer;
            for (const auto& particle : particles) {
                memcpy(ptr, &particle->position, sizeof(vec3)); ptr += sizeof(vec3);
                memcpy(ptr, &particle->velocity, sizeof(vec3)); ptr += sizeof(vec3);
                memcpy(ptr, &particle->mass, sizeof(double)); ptr += sizeof(double);
                memcpy(ptr, &particle->T, sizeof(double)); ptr += sizeof(double);
                memcpy(ptr, &particle->visualDensity, sizeof(double)); ptr += sizeof(double); // Added visualDensity not real SPH density
                memcpy(ptr, &particle->type, sizeof(int)); ptr += sizeof(int);
            }
            file.write(buffer, totalSize);
            free(buffer);
        }
    }
    else if (outputDataFormat == "AGFC")
    {
        // AGFC Format: Kompaktes Format für Rendering
        // Speichert nur position (vec3 als 3 floats), visualDensity (float) und type (int)

        // Berechnung der Gesamtgröße: 3 floats für Position, 1 float für visualDensity, 1 int für type pro Particle
        size_t totalSize = particles.size() * (sizeof(float) * 3 + sizeof(float) + sizeof(int));

        // Puffer allokieren
        std::vector<char> buffer(totalSize);
        char* ptr = buffer.data();

        for (const auto& particle : particles)
        {
            // Position als float konvertieren
            float posX = static_cast<float>(particle->position.x);
            float posY = static_cast<float>(particle->position.y);
            float posZ = static_cast<float>(particle->position.z);
            memcpy(ptr, &posX, sizeof(float)); ptr += sizeof(float);
            memcpy(ptr, &posY, sizeof(float)); ptr += sizeof(float);
            memcpy(ptr, &posZ, sizeof(float)); ptr += sizeof(float);

            // visualDensity als float konvertieren
            float visualDensity = static_cast<float>(particle->visualDensity);
            memcpy(ptr, &visualDensity, sizeof(float)); ptr += sizeof(float);

            // type als int speichern
            int type = particle->type;
            memcpy(ptr, &type, sizeof(int)); ptr += sizeof(int);
        }

        // Puffer in die Datei schreiben
        file.write(buffer.data(), totalSize);
    }
    else if (outputDataFormat == "AGFE")
    {
        // AGFE Format: Erweiterte Version (bestehend aus Position, Velocity, Mass, T, P, visualDensity, U, type)
        size_t totalSize = particles.size() * (sizeof(vec3) * 2 + sizeof(double) * 5 + sizeof(int));
        
        // Speicher für den Puffer allokieren
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
    }
    else
    {
        std::cerr << "Unknown output data format: " << outputDataFormat << std::endl;
    }

    file.close();
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

bool parseKeyValue(const std::string& line, std::string& key, std::string& value) {
    std::size_t pos = line.find('=');
    if (pos == std::string::npos) return false;

    key = line.substr(0, pos);
    value = line.substr(pos + 1);

    // Entferne Leerzeichen
    key = trim(key);
    value = trim(value);

    return !key.empty() && !value.empty();
}


bool DataManager::loadConfig(const std::string& filename, Simulation* simulation) 
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Fehler beim Öffnen der Konfigurationsdatei: " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Ignoriere leere Zeilen und Kommentare
        if (line.empty() || line[0] == '#') continue;

        // Trim die Zeile, um führende und nachfolgende Leerzeichen zu entfernen
        std::string trimmedLine = trim(line);

        // Ignoriere leere Zeilen und Kommentare
        if (trimmedLine.empty() || trimmedLine[0] == '#') continue;

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
                else if (key == "outputDataFormat") outputDataFormat = value;
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