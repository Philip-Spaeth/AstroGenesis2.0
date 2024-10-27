#pragma once

#include <iostream>
#include <memory>
#include <vector>
#include <chrono>
#include "Particle.h"
#include "vec3.h"


class ICDataReader
{
public:
    ICDataReader(){}
    ~ICDataReader(){}

//Text files, read ASCII format, slower than binary  
    void readASCII(std::string fileName, int start, int end, vec3 pos, vec3 vel, std::vector<std::shared_ptr<Particle>>& particles);

//Gadget 2 snapshot-format
    //save and load data in binary format, Gadget2 specific format
    void readGadget2Snapshot(std::string fileName, std::vector<std::shared_ptr<Particle>>& particles);

private:
    //gadget 2 specific functions
    void readBlock(std::ifstream& file, char* buffer, size_t size);

    //gadget 2 specific header
    struct GadgetHeader {
        int32_t npart[6];           // Anzahl der Partikel pro Typ
        double mass[6];             // Massen der Partikeltypen
        double time;                // Simulationszeit
        double redshift;            // Rotverschiebung
        int32_t flag_sfr;           // Sternentstehungsflag
        int32_t flag_feedback;      // Feedbackflag
        uint32_t npartTotal[6];     // Gesamtanzahl der Partikel pro Typ
        int32_t flag_cooling;       // Abkühlungsflag
        int32_t num_files;          // Anzahl der Dateien
        double boxsize;             // Boxgröße der Simulation
        double omega0;              // Dichteparameter Omega_0
        double omegaLambda;         // Dichteparameter Omega_Lambda
        double hubbleParam;         // Hubble-Parameter
        int32_t flag_stellarage;    // Sternalter-Flag
        int32_t flag_metals;        // Metallizität-Flag
        uint32_t npartTotalHighWord[6];  // High Word für Partikelanzahl
        int32_t flag_entropy_instead_u;  // Flag für Entropie
        char fill[60];              // Auffüllen auf 256 Bytes
    };
};