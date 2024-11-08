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

//Gadget 2 snapshot-format
    void readGadget2(std::string fileName, std::vector<std::shared_ptr<Particle>>& particles);
//Gadget 4 snapshot-format
    void readGadget4(std::string fileName, std::vector<std::shared_ptr<Particle>>& particles);

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
  
  //new simplified gadget 4 specific header
    struct gadget4Header
    {
    long long npart[6];      /**< number of particles of each type in this file */
    long long npartTotal[6]; /**< total number of particles of each type in this snapshot. This can be different from npart if one is dealing with a multi-file snapshot. */
    double mass[6];          /**< mass of particles of each type. If 0, then the masses are explicitly stored in the mass-block of the snapshot file, otherwise they are omitted */
    double time;                         /**< time of snapshot file */
    double redshift;                     /**< redshift of snapshot file */
    double BoxSize;                      /**< box-size of simulation in case periodic boundaries were used */
    int num_files;                       /**< number of files in multi-file snapshot */
    long long Ntrees;
    long long TotNtrees;
  };
};