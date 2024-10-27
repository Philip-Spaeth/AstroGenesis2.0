#include "ICDataReader.h"
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

//gadget2 snapshot, binary format

void ICDataReader::readBlock(std::ifstream& file, char* buffer, size_t size) {
    int32_t blockSize;
    file.read(reinterpret_cast<char*>(&blockSize), sizeof(blockSize));
    file.read(buffer, size);
    file.read(reinterpret_cast<char*>(&blockSize), sizeof(blockSize));
}

void ICDataReader::readGadget2Snapshot(std::string fileName, std::vector<std::shared_ptr<Particle>>& particles) 
{
    std::cout << "Reading Gadget2 snapshot file: " << fileName << std::endl;

    GadgetHeader header;

    std::string fileDir = "input_data/";
    std::filesystem::path filePath = fileDir;
    filePath = "../.." / filePath;
    filePath = filePath / fileName;
    
    std::ifstream file(filePath, std::ios::binary);
    if (!file) {
        std::cerr << "Fehler beim Öffnen der Datei: " << filePath << std::endl;
        return;
    }

    // Lese den Header-Block
    readBlock(file, reinterpret_cast<char*>(&header), sizeof(GadgetHeader));

    // Berechnung der Gesamtanzahl der Partikel
    int totalParticles = 0;
    for (int i = 0; i < 6; ++i) {
        if (header.npartTotal[i] > 0) {
            totalParticles += header.npartTotal[i];
        }
    }
    
    std::cout << "Total Number of Particles in the Gadget2 Snapshoot file: " << totalParticles << std::endl;
    
    if (totalParticles <= 0) {
        std::cerr << "Fehler: Keine Partikel zum Lesen vorhanden." << std::endl;
        file.close();
        return;
    }

    // Speicher für Positionen, Geschwindigkeiten und IDs reservieren
    std::vector<float> positions(3 * totalParticles);
    std::vector<float> velocities(3 * totalParticles);
    std::vector<int32_t> particleIDs(totalParticles);

    // Initialisiere Partikel-Vektor
    particles.clear();
    particles.reserve(totalParticles);

    // Lese Partikelpositionen
    readBlock(file, reinterpret_cast<char*>(positions.data()), sizeof(float) * 3 * totalParticles);

    // Lese Partikelgeschwindigkeiten
    readBlock(file, reinterpret_cast<char*>(velocities.data()), sizeof(float) * 3 * totalParticles);

    // Lese Partikel-IDs (optional, kann übersprungen werden)
    readBlock(file, reinterpret_cast<char*>(particleIDs.data()), sizeof(int32_t) * totalParticles);

    // Lese Massen (falls sie nicht im Header stehen)
    std::vector<float> masses(totalParticles, 0.0f);
    bool hasIndividualMasses = false;
    for (int i = 0; i < 6; ++i) {
        if (header.mass[i] == 0 && header.npart[i] > 0) {
            hasIndividualMasses = true;
            break;
        }
    }

    if (hasIndividualMasses) {
        readBlock(file, reinterpret_cast<char*>(masses.data()), sizeof(float) * totalParticles);
    } else {
        int particleOffset = 0;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < header.npart[i]; ++j) {
                masses[particleOffset + j] = header.mass[i];
            }
            particleOffset += header.npart[i];
        }
    }

    // Falls Gas-Partikel existieren (Typ 0), lese interne Energie (Temperatur)
    std::vector<float> internalEnergies(header.npart[0], 0.0f);
    if (header.npart[0] > 0) {
        readBlock(file, reinterpret_cast<char*>(internalEnergies.data()), sizeof(float) * header.npart[0]);
    }

    for (int i = 0; i < totalParticles; ++i) 
    {
        auto particle = std::make_shared<Particle>();

        particle->position.x = positions[3 * i] * 3.086e19; // kpc -> m
        particle->position.y = positions[3 * i + 1] * 3.086e19; // kpc -> m
        particle->position.z = positions[3 * i + 2] * 3.086e19; // kpc -> m

        particle->velocity.x = velocities[3 * i] * 1e3; // km/s -> m/s
        particle->velocity.y = velocities[3 * i + 1] *  1e3; // km/s -> m/s
        particle->velocity.z = velocities[3 * i + 2] *  1e3; // km/s -> m/s

        particle->mass = masses[i] * 1.989e30; // 1e10 M_sun -> kg
        particle->type = 1;

        if (i < header.npart[0]) {  // Nur Gas-Partikel haben interne Energie
            particle->U = internalEnergies[i] * 1e6;
        }

        particles.push_back(particle);
    }

    std::cout << "Snapshot was successfully read." << std::endl;

    file.close();
}

//Text files, read ASCII format, slower than binary
void ICDataReader::readASCII(std::string fileName, int start, int end, vec3 pos, vec3 vel, std::vector<std::shared_ptr<Particle>>& particles)
{
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
            //one of 10 particles is a gas particle
            if (currentIndex % 10 == 0) {
                particles[particleIndex]->type = 2;
            }
            else {
                particles[particleIndex]->type = 1;
            }
            particles[particleIndex]->T = 1e4;
            
            particleIndex++;
            currentIndex++;
        }

        file.close();
    }

    std::cout << "Initial Condition: " << fileName << " with particles from index " << start << " to " << (particleIndex - 1) << std::endl;
}
