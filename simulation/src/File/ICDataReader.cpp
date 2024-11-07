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

void ICDataReader::readGadget2(std::string fileName, std::vector<std::shared_ptr<Particle>>& particles)
{
    // Öffne die Datei im Binärmodus
    std::ifstream file(fileName, std::ios::binary);
    
    if (!file) {
        std::cerr << "Fehler: Konnte die Datei nicht öffnen: " << fileName << std::endl;
        return;
    }

    // Gadget2 Header auslesen
    gadget2Header header;
    file.read(reinterpret_cast<char*>(&header), sizeof(header));

    // Prüfen, ob das Einlesen erfolgreich war
    if (!file) {
        std::cerr << "Fehler: Konnte den Header aus der Datei nicht lesen!" << std::endl;
        return;
    }

    // Ausgabe des Headers zur Überprüfung
    std::cout << "Gadget2 Header geladen:" << std::endl;
    std::cout << "Teilchenanzahl pro Typ:" << std::endl;
    for (int i = 0; i < NTYPES_HEADER; ++i) {
        std::cout << "Typ " << i << ": " << header.npart[i] << " Partikel" << std::endl;
    }
    std::cout << "Gesamtanzahl Partikel (npartTotal):" << std::endl;
    for (int i = 0; i < NTYPES_HEADER; ++i) {
        std::cout << "Typ " << i << ": " << header.npartTotalLowWord[i] + header.npartTotalHighWord[i] * (1 << 32) << " Partikel" << std::endl;
    }
    std::cout << "Boxgröße: " << header.BoxSize << std::endl;
    std::cout << "Hubble-Parameter: " << header.HubbleParam << std::endl;
    std::cout << "Hubble-Konstante: " << header.Hubble << std::endl;
    std::cout << "Roteshift: " << header.redshift << std::endl;
    std::cout << "Simulationszeit: " << header.time << std::endl;

    particles.resize(header.npart[0] + header.npart[1] + header.npart[2] + header.npart[3] + header.npart[4] + header.npart[5]);

    //reading the data
    //...
    
    // Datei schließen
    file.close();

    std::cout << "Header erfolgreich ausgelesen." << std::endl;
}

void ICDataReader::readGadget4(std::string fileName, std::vector<std::shared_ptr<Particle>>& particles)
{
    // Open the file in binary mode
    std::ifstream file(fileName, std::ios::binary);

    if (!file) {
        std::cerr << "Error: Could not open file: " << fileName << std::endl;
        return;
    }

    // Gadget4 Header - Read only the header
    gadget4Header header;
    file.read(reinterpret_cast<char*>(&header), sizeof(header));

    // Check if the read was successful
    if (!file) {
        std::cerr << "Error: Could not read the header from the file!" << std::endl;
        return;
    }

    // Output header information for verification
    std::cout << "Gadget4 Header loaded:" << std::endl;
    std::cout << "Particle counts per type:" << std::endl;
    for (int i = 0; i < NTYPES_HEADER; ++i) {
        std::cout << "Type " << i << ": " << header.npart[i] << " particles" << std::endl;
    }

    std::cout << "Total particle counts (npartTotal):" << std::endl;
    for (int i = 0; i < NTYPES_HEADER; ++i) {
        std::cout << "Type " << i << ": " << header.npartTotal[i] << " particles" << std::endl;
    }

    std::cout << "Box size: " << header.BoxSize << std::endl;
    std::cout << "Simulation time: " << header.time << std::endl;
    std::cout << "Redshift: " << header.redshift << std::endl;
    std::cout << "Number of files in snapshot: " << header.num_files << std::endl;

    // Resize particles vector according to the total number of particles in the snapshot
    particles.resize(header.npartTotal[0] + header.npartTotal[1] + header.npartTotal[2] +
                     header.npartTotal[3] + header.npartTotal[4] + header.npartTotal[5]);
    
    //reading the data
    //...
    
    
    // Close the file after reading the header
    file.close();

    std::cout << "Header successfully read." << std::endl;
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
