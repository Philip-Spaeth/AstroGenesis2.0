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

bool ICDataReader::readGadgetSnapshot(std::string fileName, std::vector<std::shared_ptr<Particle>>& particles)
{
    std::cout << "reading gadget2 snapshot ..." << std::endl;
    // Öffne die Datei im Binärmodus
    std::ifstream file(fileName, std::ios::binary);
    
    if (!file) {
        std::cerr << "Fehler: Konnte die Datei nicht öffnen: " << fileName << std::endl;
        return false;
    }

    // Gadget2 Header auslesen
    gadget2Header header;

    // Lesen der ersten Blockgröße vor dem Header
    unsigned int block_size_start;
    file.read(reinterpret_cast<char*>(&block_size_start), sizeof(block_size_start));
    if (!file) {
        std::cerr << "Fehler: Konnte die Start-Blockgröße nicht lesen!" << std::endl;
        return false;
    }

    // Lesen des Headers
    file.read(reinterpret_cast<char*>(&header), sizeof(header));
    if (!file) {
        std::cerr << "Fehler: Konnte den Header aus der Datei nicht lesen!" << std::endl;
        return false;
    }

    // Lesen der Blockgröße nach dem Header
    unsigned int block_size_end;
    file.read(reinterpret_cast<char*>(&block_size_end), sizeof(block_size_end));
    if (!file) {
        std::cerr << "Fehler: Konnte die End-Blockgröße des Headers nicht lesen!" << std::endl;
        return false;
    }

    // Überprüfen, ob die Blockgrößen übereinstimmen
    if (block_size_start != sizeof(header)) {
        std::cerr << "Fehler: Start- und End-Blockgrößen des Headers stimmen nicht überein!" << std::endl;
        return false;
    }
/*
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
*/
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
        return false;
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
        return false;
    }

    // Lesen der Blockgröße nach dem POS-Block
    unsigned int pos_block_size_end;
    file.read(reinterpret_cast<char*>(&pos_block_size_end), sizeof(pos_block_size_end));
    if (!file) {
        std::cerr << "Fehler: Konnte die End-Blockgröße des POS-Blocks nicht lesen!" << std::endl;
        return false;
    }

    // Überprüfen, ob die Blockgrößen übereinstimmen
    if (pos_block_size_start != pos_block_size_end) {
        std::cerr << "Fehler: Start- und End-Blockgrößen des POS-Blocks stimmen nicht überein!" << std::endl;
        return false;
    }

    // ### Lesen des Geschwindigkeitsblocks (VEL) ###

    // Lesen der Blockgröße vor dem VEL-Block
    unsigned int vel_block_size_start;
    file.read(reinterpret_cast<char*>(&vel_block_size_start), sizeof(vel_block_size_start));
    if (!file) {
        std::cerr << "Fehler: Konnte die Start-Blockgröße des VEL-Blocks nicht lesen!" << std::endl;
        return false;
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
        return false;
    }

    // Lesen der Blockgröße nach dem VEL-Block
    unsigned int vel_block_size_end;
    file.read(reinterpret_cast<char*>(&vel_block_size_end), sizeof(vel_block_size_end));
    if (!file) {
        std::cerr << "Fehler: Konnte die End-Blockgröße des VEL-Blocks nicht lesen!" << std::endl;
        return false;
    }

    // Überprüfen, ob die Blockgrößen übereinstimmen
    if (vel_block_size_start != vel_block_size_end) {
        std::cerr << "Fehler: Start- und End-Blockgrößen des VEL-Blocks stimmen nicht überein!" << std::endl;
        return false;
    }

    // ### Lesen des ID-Blocks (ID) ###

    // Lesen der Blockgröße vor dem ID-Block
    unsigned int id_block_size_start;
    file.read(reinterpret_cast<char*>(&id_block_size_start), sizeof(id_block_size_start));
    if (!file) {
        std::cerr << "Fehler: Konnte die Start-Blockgröße des ID-Blocks nicht lesen!" << std::endl;
        return false;
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
        return false;
    }

    // Lesen der Blockgröße nach dem ID-Block
    unsigned int id_block_size_end;
    file.read(reinterpret_cast<char*>(&id_block_size_end), sizeof(id_block_size_end));
    if (!file) {
        std::cerr << "Fehler: Konnte die End-Blockgröße des ID-Blocks nicht lesen!" << std::endl;
        return false;
    }

    // Überprüfen, ob die Blockgrößen übereinstimmen
    if (id_block_size_start != id_block_size_end) {
        std::cerr << "Fehler: Start- und End-Blockgrößen des ID-Blocks stimmen nicht überein!" << std::endl;
        return false;
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
                //std::cout << "Mass for type: "<< i << " is 0 because Npart is 0" << std::endl;
            }
        }
        else
        {
            //std::cout << std::scientific << "Mass for type: "<< i << " is " <<  header.massarr[i] << std::endl;
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
            return false;
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
            return false;
        }

        // Lesen der Blockgröße nach dem MASS-Block
        unsigned int mass_block_size_end;
        file.read(reinterpret_cast<char*>(&mass_block_size_end), sizeof(mass_block_size_end));
        if (!file) {
            std::cerr << "Fehler: Konnte die End-Blockgröße des MASS-Blocks nicht lesen!" << std::endl;
            return false;
        }

        // Überprüfen, ob die Blockgrößen übereinstimmen
        if (mass_block_size_start != mass_block_size_end) {
            std::cerr << "Fehler: Start- und End-Blockgrößen des MASS-Blocks stimmen nicht überein!" << std::endl;
            return false;
        }
    }

    // ### Lesen des U-Blocks (interne Energie) ###

    // Überprüfen, ob Gaspartikel vorhanden sind
    if(header.npart[0] > 0) // Nur wenn es Gaspartikel gibt
    {
        std::cout << "reading U from gas particles" << std::endl;
        // Lesen der Blockgröße vor dem U-Block
        unsigned int u_block_size_start;
        file.read(reinterpret_cast<char*>(&u_block_size_start), sizeof(u_block_size_start));
        if (!file) {
            std::cerr << "Fehler: Konnte die Start-Blockgröße des U-Blocks nicht lesen!" << std::endl;
            return false;
        }

        // Erwartete Größe des U-Blocks berechnen: Anzahl der Gaspartikel * sizeof(float)
        unsigned int expected_u_block_size = header.npart[0] * sizeof(float);
        if (u_block_size_start != expected_u_block_size) {
            std::cerr << "Warnung: Erwartete U-Blockgröße (" << expected_u_block_size 
                      << " Bytes) stimmt nicht mit gelesener Größe (" << u_block_size_start << " Bytes) überein." << std::endl;
            // Optional: Fortfahren oder Abbruch
        }

        // Lesen der U-Daten
        std::vector<float> u_values(header.npart[0]); // Interne Energie pro Masseneinheit für Gaspartikel
        file.read(reinterpret_cast<char*>(u_values.data()), u_block_size_start);
        if (!file) {
            std::cerr << "Fehler: Konnte die U-Daten nicht lesen!" << std::endl;
            return false;
        }

        // Lesen der Blockgröße nach dem U-Block
        unsigned int u_block_size_end;
        file.read(reinterpret_cast<char*>(&u_block_size_end), sizeof(u_block_size_end));
        if (!file) {
            std::cerr << "Fehler: Konnte die End-Blockgröße des U-Blocks nicht lesen!" << std::endl;
            return false;
        }

        // Überprüfen, ob die Blockgrößen übereinstimmen
        if (u_block_size_start != u_block_size_end) {
            std::cerr << "Fehler: Start- und End-Blockgrößen des U-Blocks stimmen nicht überein!" << std::endl;
            return false;
        }
    
    unsigned int current_particle = 0;
    unsigned int gas_particle_index = 0;

    for(int type = 0; type < 6; ++type){
            for(unsigned int i = 0; i < (unsigned int)header.npart[type]; ++i){
                if(current_particle >= total_particles){
                    std::cerr << "Fehler: Überschreitung der Partikelanzahl beim Erstellen der Partikel!" << std::endl;
                    return false;
                }
                auto particle = std::make_shared<Particle>();
                particle->id = ids[current_particle];

                // Setzen des Partikeltyps und der Masse
                if(type == 1)  
                {   
                    particle->type = 3; // Halo
                }
                else if(type == 2 || type == 4 || type == 5) 
                {
                    particle->type = 1; // Disk, Bulge, Stars, Bndry
                }
                else if(type == 0) 
                {
                    particle->type = 2; // Gas
                }

                if(has_individual_mass)
                {
                    particle->mass = masses[current_particle] * Units::MSUN * 1e10;
                }
                else
                {
                    particle->mass = header.massarr[type] * Units::MSUN * 1e10;
                }

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

                // Skalierung auf SI-Einheiten
                particle->position *= Units::KPC;
                particle->velocity *= Units::KMS;

                // Interne Energie für Gaspartikel zuweisen
                if(type == 0) // Gaspartikel
                {
                    //scale to SI
                    particle->U = u_values[gas_particle_index] * 1e6;
                    //std::cout << particle->U << std::endl;
                    gas_particle_index++;
                }

                particles.push_back(particle);
                current_particle++;
            }
        }
    }
    else // Wenn keine Gaspartikel vorhanden sind
    {
        // Ihr bestehender Code zum Erstellen der Partikel
        unsigned int current_particle = 0;

        for(int type = 0; type < 6; ++type){
            for(unsigned int i = 0; i < (unsigned int)header.npart[type]; ++i){
                if(current_particle >= total_particles){
                    std::cerr << "Fehler: Überschreitung der Partikelanzahl beim Erstellen der Partikel!" << std::endl;
                    return false;
                }
                auto particle = std::make_shared<Particle>();
                particle->id = ids[current_particle];

                // Setzen des Partikeltyps und der Masse
                if(type == 1)  
                {   
                    particle->type = 3; // Halo
                }
                else if(type == 2 || type == 4 || type == 5) 
                {
                    particle->type = 1; // Disk, Bulge, Stars, Bndry
                }
                else if(type == 0) 
                {
                    particle->type = 2; // Gas
                }

                if(has_individual_mass)
                {
                    particle->mass = masses[current_particle] * Units::MSUN * 1e10;
                }
                else
                {
                    particle->mass = header.massarr[type] * Units::MSUN * 1e10;
                }

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

                // Skalierung auf SI-Einheiten
                particle->position *= Units::KPC;
                particle->velocity *= Units::KMS;

                particles.push_back(particle);
                current_particle++;
            }
        }
    }

    // Datei schließen
    file.close();

    std::cout << "gadget2 snapshot sucessfully read." << std::endl;
    return true;
}