#include "Disk.h"
#include "Constants.h"
#include <random>
#include <cmath>
#include <memory>
#include <iostream>
#include "random.h"
#include "vec3.h"

Disk::Disk() {}

Disk::~Disk() {}

void Disk::generateDisk(int start, int end, std::vector<std::shared_ptr<Particle>>& particles) 
{
    size_t N = end - start;
    particles.resize(N);

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::normal_distribution<double> normal(0.0, 1.0);

    for (int k = start; k < end; ++k) {
        // Sample radial position from the exponential profile
        double xi = uniform(gen);
        double R = -Rd * std::log(1.0 - xi);

        // Sample angular position (random orientation in the disk plane)
        double theta = 2.0 * M_PI * uniform(gen);

        // Sample vertical position from sech^2 profile
        double z = z0 * std::atanh(2.0 * uniform(gen) - 1.0);

        // Cartesian coordinates
        double x = R * std::cos(theta);
        double y = R * std::sin(theta);

        // Compute velocities
        double vc = circularVelocity(R) * 1.0; // Circular velocity
        double vr_disp = radialVelocityDispersion(R); // Radial velocity dispersion
        double vz_disp = vr_disp / sqrt(2.0); // Approximation: vertical dispersion is half radial

        double vr = vr_disp * normal(gen); // Radial velocity perturbation
        double vz = vz_disp * normal(gen); // Vertical velocity perturbation
        double vtheta = vc + vr_disp * normal(gen); // Circular velocity with some dispersion

        // Convert polar to Cartesian velocities
        double vx = vr * std::cos(theta) - vtheta * std::sin(theta);
        double vy = vr * std::sin(theta) + vtheta * std::cos(theta);

        // Assign to particle
        std::shared_ptr<Particle> particle = std::make_shared<Particle>();
        particle->position = vec3(x, y, z);
        particle->velocity = vec3(vx, vy, vz);
        particle->mass = M / static_cast<double>(N);
        particle->type = 1;
        particle->galaxyPart = 1;

        particles[k] = particle;
    }
}

double Disk::surfaceDensity(double R) const {
    return (M / (2.0 * M_PI * Rd * Rd)) * std::exp(-R / Rd);
}

double Disk::verticalDensity(double z) const {
    return (1.0 / (2.0 * z0)) * std::pow(1.0 / std::cosh(z / z0), 2);
}

double Disk::radialVelocityDispersion(double R) const 
{
    return v_disp * std::exp(-R / (2.0 * Rd));
}

double Disk::circularVelocity(double R) const 
{
    // Gravitationskonstante
    double G = Constants::G;

    // Enclosed mass within radius R
    double M_enc = M * (1 - std::exp(-R / Rd) * (1 + R / Rd));

    // Smoothed circular velocity to avoid singularity at R=0
    return std::sqrt(G * M_enc / (R));
}
