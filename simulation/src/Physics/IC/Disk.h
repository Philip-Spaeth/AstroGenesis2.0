#pragma once
#include <vector>
#include <memory>
#include "Particle.h"
#include "Constants.h"
#include "Units.h"

class Disk {
public:
    Disk();
    ~Disk();

    void generateDisk(int start, int end, std::vector<std::shared_ptr<Particle>>& particles);

    // Properties
    double M = 1.0e11 * Units::MSUN; // Total mass of the disk
    double Rd = 3.5 * Units::KPC;    // Scale radius for the exponential profile
    double z0 = 0.3 * Units::KPC;    // Scale height for sech^2 profile
    double v_disp = 1.0 * Units::KMS; // Velocity dispersion

private:
    double surfaceDensity(double R) const; // Surface density function
    double verticalDensity(double z) const; // Vertical density function
    double radialVelocityDispersion(double R) const; // Radial velocity dispersion
    double circularVelocity(double R) const; // Circular velocity
};
