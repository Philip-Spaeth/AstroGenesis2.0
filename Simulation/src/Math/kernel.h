#pragma once
#include "vec3.h"

class kernel
{
public:
	//SPH Kernels
	static double cubicSplineKernel(double r, double h);
	static double laplaceCubicSplineKernel(const vec3& rVec, double h);
	static vec3 gradientCubicSplineKernel(const vec3& r, double h);
	//SPH functions
	static double tempretureToInternalEnergy(double tempreture);
};