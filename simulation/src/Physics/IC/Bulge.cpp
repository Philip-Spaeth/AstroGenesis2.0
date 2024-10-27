#include "Bulge.h"
#include "Hernquist_Profile.h"
#include "Exponential_Profile.h"
#include "NFW_Profile.h"


// Implementierung der generateBulge-Funktion
void Bulge::generateBulge(int start, int end)
{
    // Initialisiere den Zufallszahlengenerator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_Y(0.0, 1.0); // Für Y ∈ [0, 1)
    
    double mass_per_particle = g->M_Bulge / static_cast<double>(g->N_Bulge);
    
    for (int i = start; i < end; ++i) 
    {
        // Sampling des Radius r
        double Y = dis_Y(gen);
        double r = sample_radius_bulge(Y);
        
        // Sicherstellen, dass r <= R_Bulge
        if (r > g->R_Bulge) {
            r = g->R_Bulge;
        }
        
        // Sampling der Position in Kugelkoordinaten
        double ux, uy, uz;
        sample_direction(ux, uy, uz, gen);
        double x = r * ux;
        double y = r * uy;
        double z = r * uz;
        vec3 position = g->galaxyPosition + vec3(x, y, z);
        
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
        double M_enc_Disk = 0;
        //if(g->N_Disk != 0) M_enc_Disk = Exponential_Profile::enclosedMassExponential(r, g->M_Disk, g->R_Disk, 5);
        double M_enc_Halo = 0;
        //if(g->N_Halo != 0) M_enc_Halo = NFW_Profile::enclosedMass(r, g->M_Halo, g->c_Halo, g->c_Halo);
        if(g->N_Halo + g->M_Disk)v += std::sqrt(Constants::G * (M_enc_Disk + M_enc_Halo) / r);
        
        // Zufällige Richtung der Geschwindigkeit
        double vx_dir, vy_dir, vz_dir;
        sample_direction(vx_dir, vy_dir, vz_dir, gen);
        double vx = v * vx_dir;
        double vy = v * vy_dir;
        double vz = v * vz_dir;
        
        // Berechnung der Gesamtgeschwindigkeit unter Berücksichtigung des Potentials
        vec3 velocity = g->galaxyVelocity + vec3(vx, vy, vz);
        
        // Erzeuge und füge das Partikel hinzu
        std::shared_ptr<Particle> p = std::make_shared<Particle>();
        g->particles->at(i) = p;
        p->position = position;
        p->velocity = velocity;
        p->mass = mass_per_particle;
        p->type = 1; // Bulge Star (angenommen Typ 2 für Bulge)
    }
}

// Korrekte Sampling-Funktion basierend auf der kumulativen Massenverteilung
double Bulge::sample_radius_bulge(double Y) const {
    // Berechnung von X_max basierend auf R_Bulge und Rs_Bulge
    double X_max = (g->R_Bulge * g->R_Bulge) / ((g->R_Bulge + g->Rs_Bulge) * (g->R_Bulge + g->Rs_Bulge));
    
    // Skalieren von Y auf [0, X_max)
    double X = Y * X_max;
    
    // Berechnung von r basierend auf der inversen kumulativen Verteilung
    double sqrtX = std::sqrt(X);
    return g->Rs_Bulge * sqrtX / (1.0 - sqrtX);
}

// Berechnung des Potentials für den Bulge
double Bulge::potential_bulge(double r) const 
{
    return Hernquist_Profile::PotentialHernquist(r, g->M_Bulge, g->Rs_Bulge);
}

// Berechnung der lokalen Fluchtgeschwindigkeit
double Bulge::escape_velocity_bulge(double r) const {
    return std::sqrt(2.0 * std::abs(potential_bulge(r)));
}

// Berechnung der Verteilungsfunktion f(E) (vereinfachte Annahme)
double Bulge::f_energy(double E, double r) const {
    // Für eine genaue Verteilungsfunktion müsste die Eddington-Inversion durchgeführt werden.
    // Hier verwenden wir eine vereinfachte Annahme, z.B. Maxwellian-Verteilung.
    double sigma = g->VelDis_Disk; // Geschwindigkeitsdispersion
    double prefactor = (1.0 / (std::pow(2.0 * Constants::PI, 1.5) * std::pow(sigma, 3)));
    return prefactor * std::exp(-E / (0.5 * sigma * sigma));
}

// Generierung einer zufälligen Richtung auf der Einheitskugel
void Bulge::sample_direction(double &ux, double &uy, double &uz, std::mt19937 &gen) const {
    std::uniform_real_distribution<> dis_theta(0.0, 2.0 * Constants::PI);
    std::uniform_real_distribution<> dis_cosphi(-1.0, 1.0);
    
    double theta = dis_theta(gen);
    double cosphi = dis_cosphi(gen);
    double sinphi = std::sqrt(1.0 - cosphi * cosphi);
    
    ux = sinphi * std::cos(theta);
    uy = sinphi * std::sin(theta);
    uz = cosphi;
}