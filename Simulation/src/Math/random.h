#pragma once

#include "vec3.h"

class random
{
public:
	//random Functions
	static void setRandomSeed(unsigned int seed);
	static double between(double min, double max);
	static double gaussianRandom(double mean = 0.0, double stddev = 1.0);
};