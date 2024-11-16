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

void Halo::calculateRho0() {
    double ln_term = std::log(1 + c);
    double frac_term = c / (1 + c);
    rho0 = M / (4 * Constants::PI * std::pow(rs, 3) * (ln_term - frac_term));
}

double Halo::nfwDensity(double r) const {
    if (r <= 0.0 || r > R_vir) {
        return 0.0;
    }
    return rho0 / ((r / rs) * std::pow(1 + r / rs, 2));
}

double Halo::nfwEnclosedMass(double r) const {
    if (r <= 0.0 || r > R_vir) {
        return 0.0;
    }
    double x = r / rs;
    return 4 * Constants::PI * rho0 * std::pow(rs, 3) * (std::log(1 + x) - x / (1 + x));
}

double Halo::nfwPotential(double r) const {
    if (r <= 0.0 || r > R_vir) {
        return 0.0;
    }
    double M_enc = nfwEnclosedMass(r);
    return -Constants::G * M_enc / r;
}

double Halo::sampleRadius() const {
    // Inverse CDF Sampling
    double u = random::uniform(0.0, 1.0);
    double x = rs * (std::exp(u * std::log(1 + c)) - 1);
    return x;
}

double Halo::sampleVelocity(double r) const {
    double phi = nfwPotential(r);
    return std::sqrt(-2 * phi); // Lokale Geschwindigkeit
}

void Halo::generateNFWHalo(int start, int end, std::vector<std::shared_ptr<Particle>>& particles) 
{
    calculateRho0();

    random::setSeed(42);
    size_t N = end - start;
    particles.resize(N);

    for (int k = start; k < end; ++k) {
        // Radien samplen
        double r = sampleRadius();

        // Zufällige Position
        double theta = std::acos(2.0 * random::uniform(0.0, 1.0) - 1.0);
        double phi = 2.0 * Constants::PI * random::uniform(0.0, 1.0);

        double x = r * std::sin(theta) * std::cos(phi);
        double y = r * std::sin(theta) * std::sin(phi);
        double z = r * std::cos(theta);

        // Geschwindigkeiten samplen
        double v = sampleVelocity(r);
        std::cout << v << std::endl;
        double vx = random::normal_dist(0.0, v);
        double vy = random::normal_dist(0.0, v);
        double vz = random::normal_dist(0.0, v);

        // Partikel erstellen
        auto particle = std::make_shared<Particle>();
        particle->position = vec3(x, y, z);
        particle->velocity = vec3(vx, vy, vz);
        particle->mass = M / static_cast<double>(N);

        particles[k] = particle;
    }
}

void Halo::generateHernquistHalo(int start, int end, std::vector<std::shared_ptr<Particle>>& particles)
{
    random::setSeed(42);
    
    size_t N = end - start;
    particles.resize(N);

    for(int k = start; k < end; ++k)
    {
        double r = maxR * 2;
        while (r > maxR) {
            double u = random::uniform(0.0, 1.0);
            while(u >= 1.0){
                u = random::uniform(0.0, 1.0);
            }
            double sqrt_u = std::sqrt(u);
            r = (a * sqrt_u) / (1.0 - sqrt_u);
        }

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

double Halo::hernquistVelocityDispersion(double r) const
{
    double term1 = 12 * pow(a, 3) * log((r + a) / a);
    double term2 = r * (r - a) * (2 * pow(r, 2) + 2 * r * a + pow(a, 2));
    double numerator = term1 + term2;

    double denominator = 12 * a * r * pow(r + 2 * a, 2);
    if(std::sqrt(Constants::G * M / a * (numerator / denominator)) > 0)
    {
        return  std::sqrt(Constants::G * M / a * (numerator / denominator));
    }
    return 0;
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
