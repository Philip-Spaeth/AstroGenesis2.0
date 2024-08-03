#include "kernel.h"
#include "Constants.h"

double kernel::cubicSplineKernel(double r, double h)
{
    if (r > 2 * h) {
        return 0;
    }

    double q = r / h;
    double alpha_3d = 1.0 / (Constants::PI * pow(h,3.0));
    if (q <= 1.0) 
    {
        return alpha_3d *(1.0 - (3.0*(q*q))/2.0 + (3.0*(q*q*q))/4);
    }
    else if (q <= 2.0) 
    {
        return alpha_3d *((1.0/4.0) * pow(2 - q, 3));
    }
    return 0.0;
}

vec3 kernel::gradientCubicSplineKernel(const vec3& r, double h) 
{
    vec3 grad(0.0, 0.0, 0.0);
    double q = r.length() / h;
    double factor = 15.0 / (7.0 * M_PI * std::pow(h, 3));

    if (q > 0 && q <= 1) {
        grad = vec3(r.x, r.y, r.z) * (factor * (-3.0 * q + 9.0 / 4.0 * q * q) / (h * r.length()));
    } else if (q > 1 && q <= 2) {
        grad = vec3(r.x, r.y, r.z) * (factor * (-3.0 / 4.0 * (2.0 - q) * (2.0 - q)) / (h * r.length()));
    }

    return grad;
}