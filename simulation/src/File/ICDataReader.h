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

//read ASCII format, slower than binary - change for diffrent format 
    void readASCII(std::string fileName, int start, int end, vec3 pos, vec3 vel, std::vector<std::shared_ptr<Particle>>& particles);

//Gadget snapshot-format
    void readGadgetSnapshot(std::string fileName, std::vector<std::shared_ptr<Particle>>& particles);

private:

    //old gadget 2 specific header
    struct gadget2Header
    {
        unsigned int npart[6];
        double massarr[6];
        double time;
        double redshift;
        int flag_sfr;
        int flag_feedback;
        unsigned int npartTotal[6];
        int flag_cooling;
        int num_files;
        double BoxSize;
        double Omega0;
        double OmegaLambda;
        double HubbleParam;
        int flag_stellarage;
        int flag_metals;
        unsigned int npartTotalHighWord[6];
        int flag_entropy_instead_u;
        char fill[60]; // zur Auffüllung auf 256 Bytes
  };
};