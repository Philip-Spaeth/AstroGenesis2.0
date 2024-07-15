#include "kernel.h"

namespace math
{

double kernel::cubicSplineKernel(double r, double h)
{
    const double pi = 3.14159265359;
    const double alpha_3d = 1.0 / (pi * h * h * h);
    double q = r / h;
    if (q < 1.0) {
        return alpha_3d * (1 - 1.5 * q * q + 0.75 * q * q * q);
    }
    else if (q < 2.0) {
        return alpha_3d * 0.25 * pow(2 - q, 3);
    }
    return 0.0;
}

double kernel::laplaceCubicSplineKernel(const vec3& rVec, double h)
{
    double r = sqrt(rVec.x * rVec.x + rVec.y * rVec.y + rVec.z * rVec.z);
    if (r > 2 * h) {
        return 0;
    }

    double sigma = 45.0 / (3.14 * std::pow(h, 6));
    double factor;

    if (r <= h) {
        factor = sigma * (h - r) * (3.0 * h - 3.0 * r);
    }
    else {
        factor = sigma * 3.0 * std::pow(h - r, 2);
    }

    return factor;
}

vec3 kernel::gradientCubicSplineKernel(const vec3& r, double h) 
{
    /*
    double q = sqrt(r.x * r.x + r.y * r.y + r.z * r.z) / h;
    double sigma = 10.0 / (7.0 * 3.14 * pow(h, 2));
    vec3 gradW(0,0,0);
    if (q > 0 && q <= 1.0) {
        gradW = sigma * (-3 * q + 2.25 * q * q)* r;  r / (r.length()) * h;
    } else if (q > 1.0 && q <= 2.0) {
        gradW = sigma * -0.75 * std::pow(2 - q, 2)* r ; // * r / (r.length()) * h;
    }
    */
    return vec3(0,0,0);
}

}