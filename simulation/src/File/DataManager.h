#pragma once

#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <chrono>
#include "Particle.h"
#include "vec3.h"
#include "Simulation.h"
#include "Tree/Tree.h"
#include <cstring>
#include <HighFive/HighFive.hpp>

class Simulation;
class Tree;

// Saving h5 file
struct ParticleData {
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    std::array<double, 3> acceleration;
    double mass;
    double T;
    double visualDensity;
    uint8_t type;
    uint8_t galaxyPart;
    uint32_t id;
};
struct NodeData {
    uint32_t nodeID;
    uint32_t parentID; // ID des Elternknotens, 0 wenn Root
    std::vector<uint32_t> childrenIDs; // IDs der Kindknoten
    bool isLeaf;
    uint32_t particleIndex; // Index in der Partikel-Datenliste, falls Leaf
};

class DataManager
{
public:
    DataManager(std::string path);
    ~DataManager(){}

    std::string inputPath;
    std::string inputFormat;

    //path to the folder where the simulation data is saved
    std::string outputPath;
    std::string outputFormat;

//read the Config.ini file
    bool loadConfig(const std::string& filename, Simulation* simulation);

//save data in AGF and gadget format
    void saveData(std::vector<std::shared_ptr<Particle>> particles, int timeStep, int numberTimesteps, int numberOfParticles, double deltaTime, double endTime, double currentTime);
//save Data in hdf5 format -> tree structure in file
    void saveData(std::shared_ptr<Tree> tree, int timeStep, int numberTimesteps, int numberOfParticles, double deltaTime, double endTime, double currentTime);

//
    bool loadICs(std::vector<std::shared_ptr<Particle>>& particles, Simulation* sim);

    //Data Size
    size_t ag_MemorySize;
    size_t agc_MemorySize;
    size_t age_MemorySize;
    size_t hdf5_MemorySize;
    size_t gadget_MemorySize;

private:
    bool setupFile();
    std::string ending;
    void collectParticleData(const std::shared_ptr<Node>& node, std::vector<ParticleData>& particles, std::vector<NodeData>& nodes, uint32_t parentID, uint32_t& currentNodeID);


    //AGF header
    struct AGFHeader
    {
        int numParticles[3];
        double deltaTime;
        double endTime;
        double currentTime;
    };

    //gadget2 header
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