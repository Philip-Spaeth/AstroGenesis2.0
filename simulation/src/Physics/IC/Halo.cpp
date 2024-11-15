// Halo.cpp
#include "Halo.h"
#include "Constants.h"
#include <random>
#include <cmath>
#include <memory>
#include <iostream>
#include <functional>
#include "random.h"
#include "vec3.h"

// Konstruktor
Halo::Halo()
{
}

// Destruktor
Halo::~Halo()
{
}

const double G = Constants::G;
const double M = 1e42;
const double a = 1e21;

// Hernquist Dichteverteilung
double rho_hernquist(double r) {
    return (M / (2.0 * M_PI)) * (a / (r * std::pow(r + a, 3)));
}

// Hernquist kumulative Masse
double M_hernquist(double r) {
    return M * std::pow(r, 2) / std::pow(r + a, 2);
}

// Hernquist Gravitationspotential
double phi_hernquist(double r) {
    return -G * M / (r + a);
}

// Implementierung der generateDarkMatterHalo-Funktion
void Halo::generateHalo(int start, int end, std::vector<std::shared_ptr<Particle>>& particles)
{
    random::setSeed(37683);

    size_t N_Total = end - start;
    particles.resize(N_Total);

    // Berechne E_max
    double E_max = G * M / a;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_u(0.0, 1.0);
    std::uniform_real_distribution<> dis_phi(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<> dis_cos_theta(-1.0, 1.0);

    for (int i = start; i < end; ++i)
    {
        // Generiere ein zufälliges u
        double u = dis_u(gen);
        if(u <= 0.0) u = 1e-10; // Verhindere Division durch Null und log(0)
        if(u >=1.0) u = 1.0 -1e-10; // Verhindere Division durch Null

        // Berechne r aus der inversen CDF (Hernquist-Profil)
        double sqrt_u = std::sqrt(u);
        double denominator = 1.0 - sqrt_u;
        if (denominator <= 0.0) {
            std::cerr << "Warning: denominator for r calculation is non-positive. Setting r to a large value." << std::endl;
            denominator = 1e-10; // Verhindere Division durch Null
        }
        double r = a * sqrt_u / denominator;
        
        double max_r = 1e23;
        if (r > max_r) {
            std::cerr << "Warning: r is larger than the maximum radius. Setting r to the maximum radius." << std::endl;
            r = max_r;
        }

        // Gravitationspotential an r
        double phi = phi_hernquist(r);

        double v = std::sqrt((-phi));
        
        double cos_theta = dis_cos_theta(gen);
        double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
        phi = dis_phi(gen);

        double x = r * sin_theta * std::cos(phi);
        double y = r * sin_theta * std::sin(phi);
        double z = r * cos_theta;

        // Gleiches Vorgehen für die Geschwindigkeiten
        double cos_theta_v = dis_cos_theta(gen);
        double sin_theta_v = std::sqrt(1.0 - cos_theta_v * cos_theta_v);
        double phi_v = dis_phi(gen);

        double vx = v * sin_theta_v * std::cos(phi_v);
        double vy = v * sin_theta_v * std::sin(phi_v);
        double vz = v * cos_theta_v;



        // Erstelle ein neues Teilchen
        std::shared_ptr<Particle> particle = std::make_shared<Particle>();
        particle->position = vec3(x, y, z);
        particle->velocity = vec3(vx, vy, vz);
        particle->mass = M / static_cast<double>(N_Total);

        // Füge das Teilchen zum Vektor hinzu
        particles[i - start] = particle; // Korrektur: Indexierung
    }
}