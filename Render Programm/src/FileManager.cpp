#include "FileManager.h"
#include <filesystem>
#include "Physics.h"
#include <iostream>
#include <vector>
#include "Particle.h"
#include <cmath>
#include <future>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <iomanip>
#include <unordered_map>
#include <algorithm>

#ifdef WIN32
#include <windows.h>
#endif

#include <thread>
#include <algorithm>

using namespace std;

FileManager::FileManager(std::string newDataFolder){
    dataFolder = newDataFolder;
}

FileManager::~FileManager(){}

void FileManager::loadParticles(Physics* p, int timestep, std::vector<glm::vec4>& array, std::vector<glm::vec3>& color, std::vector<glm::vec3>& densitycolor, std::vector<glm::vec3>& thermalColor, std::vector<glm::vec1>& isDarkMatter, int maxNumberOfParticles)
{
    std::string fileName = "Data/" + p->dataFolder + "/Time_" + std::to_string(timestep) + ".dat";
    std::ifstream file(fileName, std::ios::binary);

    if (file.is_open()) {
        size_t size;
        file.read(reinterpret_cast<char*>(&size), sizeof(size));

        array.resize(size);
        color.resize(size);
        densitycolor.resize(size);
        thermalColor.resize(size);
        isDarkMatter.resize(size);


        file.read(reinterpret_cast<char*>(array.data()), size * sizeof(glm::vec4));
        file.read(reinterpret_cast<char*>(color.data()), size * sizeof(glm::vec3));
        file.read(reinterpret_cast<char*>(densitycolor.data()), size * sizeof(glm::vec3));
        file.read(reinterpret_cast<char*>(thermalColor.data()), size * sizeof(glm::vec3));
        file.read(reinterpret_cast<char*>(isDarkMatter.data()), size * sizeof(glm::vec1));

        file.close();
    }
    else {
        std::cerr << "Fehler beim Öffnen der Datei zum Laden." << std::endl;
    }
}