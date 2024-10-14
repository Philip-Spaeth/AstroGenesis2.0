#include "Galaxy.h"
#include "Constants.h"
#include <iostream>
#include <random>
#include <cmath>
#include <memory>

// Konstruktor und Destruktor
Galaxy::Galaxy(std::vector<std::shared_ptr<Particle>>* particles_ptr)
    : particles(particles_ptr) {}

Galaxy::~Galaxy() {}

// Implementierung der Hernquist-Dichtefunktion
double Galaxy::densityHernquist(double r) const {
    return (M_Bulge / (2.0 * Constants::PI)) * (Rs_Bulge) / (r * std::pow(r + Rs_Bulge, 3));
}

// Implementierung des Hernquist-Potentials
double Galaxy::PotentialHernquist(double r, double M, double a) const {
    return -Constants::G * M / (r + a);
}

// Korrekte Sampling-Funktion basierend auf der kumulativen Massenverteilung
double Galaxy::sample_radius_bulge(double Y) const {
    // Berechnung von X_max basierend auf R_Bulge und Rs_Bulge
    double X_max = (R_Bulge * R_Bulge) / ((R_Bulge + Rs_Bulge) * (R_Bulge + Rs_Bulge));
    
    // Skalieren von Y auf [0, X_max)
    double X = Y * X_max;
    
    // Berechnung von r basierend auf der inversen kumulativen Verteilung
    double sqrtX = std::sqrt(X);
    return Rs_Bulge * sqrtX / (1.0 - sqrtX);
}

// Berechnung des Potentials für den Bulge
double Galaxy::potential_bulge(double r) const {
    return PotentialHernquist(r, M_Bulge, Rs_Bulge);
}

// Berechnung der lokalen Fluchtgeschwindigkeit
double Galaxy::escape_velocity_bulge(double r) const {
    return std::sqrt(2.0 * std::abs(potential_bulge(r)));
}

// Berechnung der Verteilungsfunktion f(E) (vereinfachte Annahme)
double Galaxy::f_energy(double E, double r) const {
    // Für eine genaue Verteilungsfunktion müsste die Eddington-Inversion durchgeführt werden.
    // Hier verwenden wir eine vereinfachte Annahme, z.B. Maxwellian-Verteilung.
    double sigma = VelDis_Disk; // Geschwindigkeitsdispersion
    double prefactor = (1.0 / (std::pow(2.0 * Constants::PI, 1.5) * std::pow(sigma, 3)));
    return prefactor * std::exp(-E / (0.5 * sigma * sigma));
}

// Generierung einer zufälligen Richtung auf der Einheitskugel
void Galaxy::sample_direction(double &ux, double &uy, double &uz, std::mt19937 &gen) const {
    std::uniform_real_distribution<> dis_theta(0.0, 2.0 * Constants::PI);
    std::uniform_real_distribution<> dis_cosphi(-1.0, 1.0);
    
    double theta = dis_theta(gen);
    double cosphi = dis_cosphi(gen);
    double sinphi = std::sqrt(1.0 - cosphi * cosphi);
    
    ux = sinphi * std::cos(theta);
    uy = sinphi * std::sin(theta);
    uz = cosphi;
}

// Implementierung der generateGalaxy-Funktion
void Galaxy::generateGalaxy()
{
    particles->clear();
    particles->reserve(N_Halo + N_Disk + N_Bulge + N_Gas);
    particles->resize(N_Halo + N_Disk + N_Bulge + N_Gas);

    generateStarDisk(0, N_Disk);
    generateBulge(N_Disk, N_Disk + N_Bulge);

    // Falls vorhanden, generieren Sie auch den Halo und Gas
    // Beispiel:
    // generateHalo(N_Disk + N_Bulge, N_Disk + N_Bulge + N_Halo);
    // generateGas(N_Disk + N_Bulge + N_Halo, N_Disk + N_Bulge + N_Halo + N_Gas);

    std::cout << "Galaxy generated with " << particles->size() << " particles." << std::endl;
}

double Galaxy::enclosedMassDisk(double R) const 
{
    // Exponentielles Radialprofil
    double R_d = R_Disk; // Skalierungslänge
    double mass_enclosed = M_Disk * (1.0 - exp(-R / R_d) * (1.0 + R / R_d));
    return mass_enclosed;
}

double Galaxy::enclosedMassBulge(double r) const 
{
    // Hernquist-Profil: M_enc = M_Bulge * r^2 / (R_Bulge + r)^2
    return M_Bulge * (r * r) / ((R_Bulge + r) * (R_Bulge + r));
}


// Implementierung der generateStarDisk-Funktion
void Galaxy::generateStarDisk(int start, int end)
{
    // Initialisiere den Zufallszahlengenerator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Verteilungsfunktionen
    std::uniform_real_distribution<> dist_theta(0.0, 2.0 * Constants::PI);
    std::uniform_real_distribution<> dist_uniform(0.0, 1.0); 
    std::normal_distribution<> dist_z(0.0, z_Disk); // Gaussian für z
    std::normal_distribution<> dist_velocity(0.0, VelDis_Disk); // Geschwindigkeitsdispersion

    double R_d = R_Disk; // Skalierungslänge der Scheibe in Metern
    double total_mass = M_Disk; // Gesamtmasse der Scheibe in kg
    double mass_per_particle = total_mass / static_cast<double>(N_Disk);

    for (int i = start; i < end; ++i) 
    {
        double rand_num = dist_uniform(gen);
        double R = -R_d * std::log(1.0 - rand_num);
        double theta = dist_theta(gen);
        double z = dist_z(gen);

        double x = R * std::cos(theta);
        double y = R * std::sin(theta);
        vec3 position = galaxyPosition + vec3(x, y, z);

        double M_enc_disk = enclosedMassDisk(R) + enclosedMassBulge(R); // Gesamtmasse in der Scheibe
        double v_theta = std::sqrt(Constants::G * M_enc_disk / (R + 1e-10)); // +1e-10 um Division durch Null zu vermeiden

        double v_r = dist_velocity(gen); // Radialgeschwindigkeit
        double v_z = dist_velocity(gen); // Vertikale Geschwindigkeit

        // Für die tangentiale Geschwindigkeit
        double v_phi = v_theta * (1 + (dist_velocity(gen) / v_theta)); // Prozentsatz der Dispersion

        // Berechnung der kartesischen Geschwindigkeiten
        double vx = -v_phi * std::sin(theta) + v_r * std::cos(theta);
        double vy = v_phi * std::cos(theta) + v_r * std::sin(theta);
        double vz = v_z;

        vec3 velocity = galaxyVelocity + vec3(vx, vy, vz);

        // Erzeuge und füge das Partikel hinzu
        std::shared_ptr<Particle> p = std::make_shared<Particle>();
        particles->at(i) = p;
        p->position = position;
        p->velocity = velocity;
        p->mass = mass_per_particle;
        p->type = 1; // Star
    }
}

// Implementierung der generateBulge-Funktion
void Galaxy::generateBulge(int start, int end)
{
    // Initialisiere den Zufallszahlengenerator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_Y(0.0, 1.0); // Für Y ∈ [0, 1)
    
    double mass_per_particle = M_Bulge / static_cast<double>(N_Bulge);
    
    for (int i = start; i < end; ++i) 
    {
        // Sampling des Radius r
        double Y = dis_Y(gen);
        double r = sample_radius_bulge(Y);
        
        // Sicherstellen, dass r <= R_Bulge
        if (r > R_Bulge) {
            r = R_Bulge;
        }
        
        // Sampling der Position in Kugelkoordinaten
        double ux, uy, uz;
        sample_direction(ux, uy, uz, gen);
        double x = r * ux;
        double y = r * uy;
        double z = r * uz;
        vec3 position = galaxyPosition + vec3(x, y, z);
        
        // Berechnung des Potentials an diesem Radius
        double phi = potential_bulge(r);
        
        // Berechnung der Fluchtgeschwindigkeit
        double v_esc = escape_velocity_bulge(r);
        
        // Rejection-Sampling für die Geschwindigkeit
        double v;
        std::uniform_real_distribution<> dis_v(0.0, v_esc);
        while (true) {
            v = dis_v(gen);
            double E = 0.5 * v * v + phi;
            if (E < 0) { // Teilchen sind gebunden
                break;
            }
        }
        v += std::sqrt(Constants::G * enclosedMassDisk(r) / r);
        
        // Zufällige Richtung der Geschwindigkeit
        double vx_dir, vy_dir, vz_dir;
        sample_direction(vx_dir, vy_dir, vz_dir, gen);
        double vx = v * vx_dir;
        double vy = v * vy_dir;
        double vz = v * vz_dir;
        
        // Berechnung der Gesamtgeschwindigkeit unter Berücksichtigung des Potentials
        vec3 velocity = galaxyVelocity + vec3(vx, vy, vz);
        
        // Erzeuge und füge das Partikel hinzu
        std::shared_ptr<Particle> p = std::make_shared<Particle>();
        particles->at(i) = p;
        p->position = position;
        p->velocity = velocity;
        p->mass = mass_per_particle;
        p->type = 2; // Bulge Star (angenommen Typ 2 für Bulge)
    }
}
