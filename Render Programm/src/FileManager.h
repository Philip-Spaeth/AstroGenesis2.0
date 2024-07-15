#pragma once

#include <glm.hpp>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <unordered_map>

class Physics;


class FileManager
{
	public:
	FileManager(std::string dataFolder);
	~FileManager();

	std::string dataFolder = "Data";
	void loadParticles(Physics* p, int timestep, std::vector<glm::vec4>& array, std::vector<glm::vec3>& color, std::vector<glm::vec3>& densitycolor, std::vector<glm::vec3>& thermalColor, std::vector<glm::vec1>& isDarkMatter, int maxNumberOfParticles);

};