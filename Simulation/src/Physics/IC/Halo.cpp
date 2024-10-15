// Halo.cpp
#include "Halo.h"
#include "NFW_Profile.h"
#include "Constants.h"
#include <random>
#include <cmath>
#include <memory>
#include <iostream>

// Konstruktor
Halo::Halo(Galaxy * g_ptr)
    : g(g_ptr)
{
}

// Destruktor
Halo::~Halo()
{
}

// Implementierung der generateDarkMatterHalo-Funktion
void Halo::generateDarkMatterHalo(int start, int end)
{
    // Initialisiere den Zufallszahlengenerator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Verteilungsfunktionen
    std::uniform_real_distribution<> dist_U(0.0, 1.0); // Für U ∈ [0,1)
    std::uniform_real_distribution<> dist_theta(0.0, 2.0 * Constants::PI);
    std::uniform_real_distribution<> dist_cosphi(-1.0, 1.0);
    std::normal_distribution<> dist_gauss(0.0, 0.1); // Geschwindigkeitsstreuung als Prozentsatz

    // Skalierungsparameter
    double R_vir = g->R_Halo; // Virialradius (z.B., 3e21 m für 100 kpc)
    double c = g->c_Halo;      // Konzentrationsparameter (z.B., 5.0)
    double Rs = R_vir / c;     // Skalierungsradius
    double M = g->M_Halo;      // Gesamtmasse des Halos

    // Berechnung der charakteristischen Dichte rho0 (optional, für Debugging)
    double rho0 = NFW_Profile::density(Rs, M, Rs, c);
    // Optional: Debugging-Ausgabe
    // std::cout << "rho0: " << rho0 << " kg/m^3\n";

    for(int i = start; i < end; ++i)
    {
        // Sampling des Radiums r mittels inverser CDF-Methode
        double U = dist_U(gen);
        double r = NFW_Profile::sample_radius(U, M, Rs, c, R_vir);

        // Sampling der Richtung
        double theta = dist_theta(gen);
        double cosphi = dist_cosphi(gen);
        double sinphi = std::sqrt(1.0 - cosphi * cosphi);

        double ux = sinphi * std::cos(theta);
        double uy = sinphi * std::sin(theta);
        double uz = cosphi;

        double x = r * ux;
        double y = r * uy;
        double z = r * uz;

        vec3 position = g->galaxyPosition + vec3(x, y, z);

        // Berechnung des Potentials an diesem Radius
        double phi = NFW_Profile::potential(r, M, Rs, c);

        // Berechnung der Fluchtgeschwindigkeit
        double v_esc = std::sqrt(2.0 * std::abs(phi));

        // Berechnung der Kreisgeschwindigkeit
        double M_enc = NFW_Profile::enclosedMass(r, M, Rs, c);
        double v_c = std::sqrt(Constants::G * M_enc / (r + 1e-10)); // +1e-10 zur Vermeidung von Division durch Null

        // Hinzufügen einer Streuung zur Kreisgeschwindigkeit
        double delta_v = dist_gauss(gen); // ±10% Streuung
        double v = v_c * (1.0 + delta_v); // v = v_c * (1 ± 0.1)

        // Sicherstellen, dass v nicht die Fluchtgeschwindigkeit überschreitet
        if(v > v_esc)
            v = 0.99 * v_esc; // Setze v etwas unter die Fluchtgeschwindigkeit

        // Zufällige Richtung der Geschwindigkeit (isotrop)
        double theta_v = dist_theta(gen);
        double cosphi_v = dist_cosphi(gen);
        double sinphi_v = std::sqrt(1.0 - cosphi_v * cosphi_v);

        double vx_dir = sinphi_v * std::cos(theta_v);
        double vy_dir = sinphi_v * std::sin(theta_v);
        double vz_dir = cosphi_v;

        // Hinzufügen der zufälligen Streuung
        vx_dir += dist_gauss(gen) * 0.1;
        vy_dir += dist_gauss(gen) * 0.1;
        vz_dir += dist_gauss(gen) * 0.1;

        // Normalisieren der Richtung
        double norm = std::sqrt(vx_dir * vx_dir + vy_dir * vy_dir + vz_dir * vz_dir);
        if(norm > 0)
        {
            vx_dir /= norm;
            vy_dir /= norm;
            vz_dir /= norm;
        }
        else
        {
            vx_dir = 1.0;
            vy_dir = 0.0;
            vz_dir = 0.0;
        }

        double vx = v * vx_dir;
        double vy = v * vy_dir;
        double vz = v * vz_dir;

        vec3 velocity = g->galaxyVelocity + vec3(vx, vy, vz);

        // Erzeuge und füge das Partikel hinzu
        std::shared_ptr<Particle> p = std::make_shared<Particle>();
        g->particles->at(i) = p;
        p->position = position;
        p->velocity = velocity;
        p->mass = M / static_cast<double>(g->N_Halo); // mass_per_particle
        p->type = 3; // Dark Matter Halo (angenommen Typ 3 für Halo)
    }
}
