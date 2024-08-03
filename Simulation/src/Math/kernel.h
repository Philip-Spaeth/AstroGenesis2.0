#pragma once
#include "vec3.h"

class kernel
{
public:
	//SPH Kernels
	static double cubicSplineKernel(double r, double h);
	static vec3 gradientCubicSplineKernel(const vec3& r, double h);
};