#include "Halo.h"
#include "Constants.h"
#include <random>
#include <cmath>
#include <memory>
#include <iostream>
#include <functional>
#include "random.h"
#include "vec3.h"

Halo::Halo()
{
}

Halo::~Halo()
{
}

void Halo::generateHalo(int start, int end, std::vector<std::shared_ptr<Particle>>& particles)
{
    size_t N = end - start;
    particles.resize(N);

    for(int k = start; k < end; ++k)
    {
        double u = random::uniform(0.0, 1.0);
        while(u >= 1.0){
            u = random::uniform(0.0, 1.0);
        }
        double sqrt_u = std::sqrt(u);
        double r = (a * sqrt_u) / (1.0 - sqrt_u);

        double v = random::uniform(0.0, 1.0);
        double w = random::uniform(0.0, 1.0);

        double cos_theta = 2.0 * v - 1.0;
        double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
        double phi = 2.0 * Constants::PI * w;

        double x = r * sin_theta * std::cos(phi);
        double y = r * sin_theta * std::sin(phi);
        double z = r * cos_theta;

        double sigma = hernquistVelocityDispersion(r);
        
        double vx = random::normal_dist(0.0, sigma);
        double vy = random::normal_dist(0.0, sigma);
        double vz = random::normal_dist(0.0, sigma);

        std::shared_ptr<Particle> particle = std::make_shared<Particle>();
        particle->position = vec3(x, y, z);
        particle->velocity = vec3(vx, vy, vz);

        particle->mass = M / static_cast<double>(N);
        particle->type = 3;
        particle->galaxyPart = 3;

        particles[k] = particle;
    }
}

double Halo::hernquistVelocity(double r) const
{
    if (r < 0.0) {
        throw std::invalid_argument("Radius r kann nicht negativ sein.");
    }

    if (r == 0.0) {
        return 0.0;
    }

    double enclosedMass = Halo::hernquistEnclosedMass(r);
    return std::sqrt(Constants::G * enclosedMass / r);
}

double Halo::hernquistVelocityDispersion(double r) const
{
    if (r <= 0.0) {
        return 0.0;
    }

    double x = r / a;
    double log_term = std::log(1.0 + x);
    double term1 = 12.0 * (1.0 + x) * log_term;
    double term2 = 12.0 * x;
    double sigma_sq = (Constants::G * M) / (12.0 * a) * (term1 - term2) / std::pow(1.0 + x, 2);

    if(sigma_sq < 0.0){
        std::cerr << "Warnung: Negative Geschwindigkeitsdispersion an r = " << r << std::endl;
        return 0.0;
    }

    return std::sqrt(sigma_sq);
}


double Halo::hernquistDensity(double r) const
{
    if (r <= 0) {
        return 0.0;
    }

    double rho = (M / (2 * Constants::PI)) * (a / (r * pow(r + a, 3))); 
    return rho;
}

double Halo::hernquistEnclosedMass(double r) const
{
    if (r <= 0) {
        return 0.0;
    }

    double Menc = M * (r * r) / pow(r + a, 2);
    return Menc;
}

double Halo::hernquistPotential(double r) const
{
    if (r == 0.0) {
        return -Constants::G * M / a;
    }

    return -Constants::G * M / (r + a);
}
