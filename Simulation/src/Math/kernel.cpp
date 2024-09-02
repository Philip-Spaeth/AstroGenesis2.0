#include "kernel.h"
#include "Constants.h"

double kernel::cubicSplineKernel(double r, double h)
{
    const double alpha_3d = 1.0 / (Constants::PI * h * h * h);
    double q = r / h;
    if (q < 1.0) {
        return alpha_3d * (1 - 1.5 * q * q + 0.75 * q * q * q);
    }
    else if (q < 2.0) {
        return alpha_3d * 0.25 * pow(2 - q, 3);
    }
    return 0.0;
}

vec3 kernel::gradientCubicSplineKernel(const vec3& r, double h) 
{
    double q = sqrt(r.x * r.x + r.y * r.y + r.z * r.z) / h;
    double sigma = 10.0 / (7.0 * 3.14 * pow(h, 2));
    vec3 gradW(0,0,0);

    if (q > 0 && q <= 1.0) {
        gradW = sigma * (-3 * q + 2.25 * q * q) * r / (r.length() * h);
    }
    else if (q > 1.0 && q <= 2.0) {
        gradW = sigma * -0.75 * pow(2 - q, 2) * r / (r.length() * h);
    }

    return gradW;
}

double kernel::softeningKernel(double u)
{
    if (u >= 0 && u < 1.0/2.0)
    {
        return (16.0/3.0) *std::pow(u, 2) - (48.0/5.0) * std::pow(u, 4) + (32.0/5.0) * std::pow(u, 5) - (14.0/5.0);
    }
    else if (u >= 1.0/2.0 && u < 1.0)
    {
        return (1.0/(15.0 * u)) + (32.0/3.0) * std::pow(u, 2) - (16.0) * std::pow(u, 3) - (48.0/5.0) * std::pow(u, 4) - (32.0/15.0) * std::pow(u, 5) - (16.0/5.0);
    }
    else if (u >= 1.0)
    {
        return -1.0 / u;
    }
    return 0.0;
}