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
#include "Units.h"

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

    // Lesen der ersten Blockgröße vor dem Header
    unsigned int block_size_start;
    file.read(reinterpret_cast<char*>(&block_size_start), sizeof(block_size_start));
    if (!file) {
        std::cerr << "Fehler: Konnte die Start-Blockgröße nicht lesen!" << std::endl;
        return;
    }

    // Lesen des Headers
    file.read(reinterpret_cast<char*>(&header), sizeof(header));
    if (!file) {
        std::cerr << "Fehler: Konnte den Header aus der Datei nicht lesen!" << std::endl;
        return;
    }

    // Lesen der Blockgröße nach dem Header
    unsigned int block_size_end;
    file.read(reinterpret_cast<char*>(&block_size_end), sizeof(block_size_end));
    if (!file) {
        std::cerr << "Fehler: Konnte die End-Blockgröße des Headers nicht lesen!" << std::endl;
        return;
    }

    // Überprüfen, ob die Blockgrößen übereinstimmen
    if (block_size_start != sizeof(header)) {
        std::cerr << "Fehler: Start- und End-Blockgrößen des Headers stimmen nicht überein!" << std::endl;
        return;
    }

    // Ausgabe des Headers zur Überprüfung
    std::cout << "Gadget2 Header geladen:" << std::endl;
    std::cout << "Teilchenanzahl pro Typ:" << std::endl;
    for (int i = 0; i < 6; ++i) { // Annahme: 6 Typen
        std::cout << "Typ " << i << ": " << header.npart[i] << " Partikel" << std::endl;
    }
    std::cout << "Mass per type:" << std::endl;
    for (int i = 0; i < 6; ++i) { // Annahme: 6 Typen
        std::cout << "Typ " << i << ": " << header.massarr[i] << " Mass" << std::endl;
    }
    std::cout << "Boxgröße: " << header.BoxSize << std::endl;
    std::cout << "Hubble-Parameter: " << header.HubbleParam << std::endl;
    std::cout << "Roteshift: " << header.redshift << std::endl;
    std::cout << "Simulationszeit: " << header.time << std::endl;

    // Gesamtanzahl der Partikel berechnen
    unsigned int total_particles = 0;
    for(int i = 0; i < 6; ++i){
        total_particles += header.npart[i];
    }

    // Reservieren des Speicherplatzes für Partikel
    particles.reserve(total_particles);

    // ### Lesen des Positionsblocks (POS) ###
    
    // Lesen der Blockgröße vor dem POS-Block
    unsigned int pos_block_size_start;
    file.read(reinterpret_cast<char*>(&pos_block_size_start), sizeof(pos_block_size_start));
    if (!file) {
        std::cerr << "Fehler: Konnte die Start-Blockgröße des POS-Blocks nicht lesen!" << std::endl;
        return;
    }

    // Berechnen der erwarteten Größe für den POS-Block: N * 3 * sizeof(float)
    unsigned int expected_pos_block_size = total_particles * 3 * sizeof(float);
    if (pos_block_size_start != expected_pos_block_size) {
        std::cerr << "Warnung: Erwartete POS-Blockgröße (" << expected_pos_block_size 
                  << " bytes) stimmt nicht mit gelesener Größe (" << pos_block_size_start << " bytes) überein." << std::endl;
        // Optional: Fortfahren oder Abbruch
    }

    // Lesen der Positionsdaten
    std::vector<float> positions(total_particles * 3); // N * 3 floats
    file.read(reinterpret_cast<char*>(positions.data()), pos_block_size_start);
    if (!file) {
        std::cerr << "Fehler: Konnte die Positionsdaten nicht lesen!" << std::endl;
        return;
    }

    // Lesen der Blockgröße nach dem POS-Block
    unsigned int pos_block_size_end;
    file.read(reinterpret_cast<char*>(&pos_block_size_end), sizeof(pos_block_size_end));
    if (!file) {
        std::cerr << "Fehler: Konnte die End-Blockgröße des POS-Blocks nicht lesen!" << std::endl;
        return;
    }

    // Überprüfen, ob die Blockgrößen übereinstimmen
    if (pos_block_size_start != pos_block_size_end) {
        std::cerr << "Fehler: Start- und End-Blockgrößen des POS-Blocks stimmen nicht überein!" << std::endl;
        return;
    }

    // ### Lesen des Geschwindigkeitsblocks (VEL) ###

    // Lesen der Blockgröße vor dem VEL-Block
    unsigned int vel_block_size_start;
    file.read(reinterpret_cast<char*>(&vel_block_size_start), sizeof(vel_block_size_start));
    if (!file) {
        std::cerr << "Fehler: Konnte die Start-Blockgröße des VEL-Blocks nicht lesen!" << std::endl;
        return;
    }

    // Berechnen der erwarteten Größe für den VEL-Block: N * 3 * sizeof(float)
    unsigned int expected_vel_block_size = total_particles * 3 * sizeof(float);
    if (vel_block_size_start != expected_vel_block_size) {
        std::cerr << "Warnung: Erwartete VEL-Blockgröße (" << expected_vel_block_size 
                  << " bytes) stimmt nicht mit gelesener Größe (" << vel_block_size_start << " bytes) überein." << std::endl;
        // Optional: Fortfahren oder Abbruch
    }

    // Lesen der Geschwindigkeitsdaten
    std::vector<float> velocities(total_particles * 3); // N * 3 floats
    file.read(reinterpret_cast<char*>(velocities.data()), vel_block_size_start);
    if (!file) {
        std::cerr << "Fehler: Konnte die Geschwindigkeitsdaten nicht lesen!" << std::endl;
        return;
    }

    // Lesen der Blockgröße nach dem VEL-Block
    unsigned int vel_block_size_end;
    file.read(reinterpret_cast<char*>(&vel_block_size_end), sizeof(vel_block_size_end));
    if (!file) {
        std::cerr << "Fehler: Konnte die End-Blockgröße des VEL-Blocks nicht lesen!" << std::endl;
        return;
    }

    // Überprüfen, ob die Blockgrößen übereinstimmen
    if (vel_block_size_start != vel_block_size_end) {
        std::cerr << "Fehler: Start- und End-Blockgrößen des VEL-Blocks stimmen nicht überein!" << std::endl;
        return;
    }

    // ### Lesen des ID-Blocks (ID) ###

    // Lesen der Blockgröße vor dem ID-Block
    unsigned int id_block_size_start;
    file.read(reinterpret_cast<char*>(&id_block_size_start), sizeof(id_block_size_start));
    if (!file) {
        std::cerr << "Fehler: Konnte die Start-Blockgröße des ID-Blocks nicht lesen!" << std::endl;
        return;
    }

    // Berechnen der erwarteten Größe für den ID-Block: N * sizeof(unsigned int)
    unsigned int expected_id_block_size = total_particles * sizeof(unsigned int);
    if (id_block_size_start != expected_id_block_size) {
        std::cerr << "Warnung: Erwartete ID-Blockgröße (" << expected_id_block_size 
                  << " bytes) stimmt nicht mit gelesener Größe (" << id_block_size_start << " bytes) überein." << std::endl;
        // Optional: Fortfahren oder Abbruch
    }

    // Lesen der ID-Daten
    std::vector<unsigned int> ids(total_particles); // N unsigned ints
    file.read(reinterpret_cast<char*>(ids.data()), id_block_size_start);
    if (!file) {
        std::cerr << "Fehler: Konnte die ID-Daten nicht lesen!" << std::endl;
        return;
    }

    // Lesen der Blockgröße nach dem ID-Block
    unsigned int id_block_size_end;
    file.read(reinterpret_cast<char*>(&id_block_size_end), sizeof(id_block_size_end));
    if (!file) {
        std::cerr << "Fehler: Konnte die End-Blockgröße des ID-Blocks nicht lesen!" << std::endl;
        return;
    }

    // Überprüfen, ob die Blockgrößen übereinstimmen
    if (id_block_size_start != id_block_size_end) {
        std::cerr << "Fehler: Start- und End-Blockgrößen des ID-Blocks stimmen nicht überein!" << std::endl;
        return;
    }

    //read mass:
    // Überprüfen, ob individuelle Massen vorhanden sind (massarr[i] == 0)
    const float epsilon = 1e-10;
    bool has_individual_mass = false;
    for(int i = 0; i < 6; ++i){
        if(header.massarr[i] < epsilon)
        {
            if(header.npart[i] != 0)
            {
                has_individual_mass = true;
                break;
            }
            else
            {
                std::cout << "Mass for type: "<< i << " is 0 because Npart is 0" << std::endl;
            }
        }
        else
        {
            std::cout << std::scientific << "Mass for type: "<< i << " is " <<  header.massarr[i] << std::endl;
        }
    }

    std::vector<float> masses; // Vektor zur Speicherung der Massen
    if(has_individual_mass)
    {
        // Lesen der Blockgröße vor dem MASS-Block
        unsigned int mass_block_size_start;
        file.read(reinterpret_cast<char*>(&mass_block_size_start), sizeof(mass_block_size_start));
        if (!file) {
            std::cerr << "Fehler: Konnte die Start-Blockgröße des MASS-Blocks nicht lesen!" << std::endl;
            return;
        }

        // Berechnen der erwarteten Größe für den MASS-Block: N * sizeof(float)
        unsigned int expected_mass_block_size = total_particles * sizeof(float);
        if (mass_block_size_start != expected_mass_block_size) {
            std::cerr << "Warnung: Erwartete MASS-Blockgröße (" << expected_mass_block_size 
                      << " bytes) stimmt nicht mit gelesener Größe (" << mass_block_size_start << " bytes) überein." << std::endl;
            // Optional: Fortfahren oder Abbruch
        }

        // Lesen der Massendaten
        masses.resize(total_particles);
        file.read(reinterpret_cast<char*>(masses.data()), mass_block_size_start);
        if (!file) {
            std::cerr << "Fehler: Konnte die Massendaten nicht lesen!" << std::endl;
            return;
        }

        // Lesen der Blockgröße nach dem MASS-Block
        unsigned int mass_block_size_end;
        file.read(reinterpret_cast<char*>(&mass_block_size_end), sizeof(mass_block_size_end));
        if (!file) {
            std::cerr << "Fehler: Konnte die End-Blockgröße des MASS-Blocks nicht lesen!" << std::endl;
            return;
        }

        // Überprüfen, ob die Blockgrößen übereinstimmen
        if (mass_block_size_start != mass_block_size_end) {
            std::cerr << "Fehler: Start- und End-Blockgrößen des MASS-Blocks stimmen nicht überein!" << std::endl;
            return;
        }
    }

    // ### Erstellen der Partikel mit ID, Position und Geschwindigkeit ###

    std::cout << "\nErstelle Partikel mit ID, Position und Geschwindigkeit..." << std::endl;
    unsigned int current_particle = 0;

    for(int type = 0; type < 6; ++type){
        for(unsigned int i = 0; i < (unsigned int)header.npart[type]; ++i){
            if(current_particle >= total_particles){
                std::cerr << "Fehler: Überschreitung der Partikelanzahl beim Erstellen der Partikel!" << std::endl;
                return;
            }
            auto particle = std::make_shared<Particle>();
            particle->id = ids[current_particle];
            //if halo
            if(type == 1)  
            {   
                particle->type = 3;
            }
            //if disk, bulge, star or Bndry
            if(type == 2 || type == 4 || type == 5) 
            {
                particle->type = 1;
            }
            //if gas
            if(type == 0) 
            {
                particle->type = 2;
            }

            particle->mass = header.massarr[type] * Units::MSUN * 1e10;

            particle->position = vec3(
                positions[3*current_particle],
                positions[3*current_particle + 1],
                positions[3*current_particle + 2]
            );
            particle->velocity = vec3(
                velocities[3*current_particle],
                velocities[3*current_particle + 1],
                velocities[3*current_particle + 2]
            );

            //scale to SI units
            particle->position *= Units::KPC;
            particle->velocity *= Units::KMS;

            particles.push_back(particle);
            current_particle++;
        }
    }


    // ### Ausgabe ausgewählter Partikel ###

    std::cout << "\nAusgabe ausgewählter Partikel:" << std::endl;
    std::vector<unsigned int> selected_indices = {10000, 50000};
    for(auto idx : selected_indices){
        if(idx == 0 || idx > total_particles){
            std::cerr << "Warnung: Ungültiger Partikelindex: " << idx << std::endl;
            continue;
        }
        unsigned int particle_idx = idx - 1;
        std::cout << "Partikel " << idx << ": ID = " << particles[particle_idx]->id 
                  << ", Position = (" << particles[particle_idx]->position.x << ", " 
                  << particles[particle_idx]->position.y << ", " << particles[particle_idx]->position.z << ")"
                  << ", Velocity = (" << particles[particle_idx]->velocity.x << ", " 
                  << particles[particle_idx]->velocity.y << ", " << particles[particle_idx]->velocity.z << ")"
                  << ", Mass = " << particles[particle_idx]->mass << std::endl;
    }

    // Datei schließen
    file.close();

    std::cout << "\nHeader, Positionen, Geschwindigkeiten, IDs und Massen erfolgreich ausgelesen und zugewiesen." << std::endl;
}

void ICDataReader::readGadget4(std::string fileName, std::vector<std::shared_ptr<Particle>>& particles)
{
//Not working ///////////////////////////////7
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
    for (int i = 0; i < 6; ++i) {
        std::cout << "Type " << i << ": " << header.npart[i] << " particles" << std::endl;
    }

    std::cout << "Box size: " << header.BoxSize << std::endl;
    std::cout << "Simulation time: " << header.time << std::endl;
    std::cout << "Redshift: " << header.redshift << std::endl;
    std::cout << "Number of files in snapshot: " << header.num_files << std::endl;

    // Resize particles vector according to the total number of particles in the snapshot
    particles.resize(header.npartTotal[0] + header.npartTotal[1] + header.npartTotal[2] +
                     header.npartTotal[3] + header.npartTotal[4] + header.npartTotal[5]);

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
