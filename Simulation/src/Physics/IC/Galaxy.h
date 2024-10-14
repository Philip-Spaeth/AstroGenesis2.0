#ifndef GALAXY_H
#define GALAXY_H

#include "vec3.h"
#include <vector>
#include <memory>
#include <random> // Hinzugefügt für std::mt19937
#include "Particle.h"

class Galaxy {
public:
    // Konstruktor und Destruktor
    Galaxy(std::vector<std::shared_ptr<Particle>>* particles);
    ~Galaxy();

    // Eigenschaften der Galaxie

    // Dark Matter Halo mit NFW-Profil
    double M_Halo;
    double R_Halo;
    double c_Halo;
    int N_Halo;

    // Sternenscheibe mit exponentiellem Profil
    double M_Disk;
    double R_Disk;
    double z_Disk;
    int N_Disk;

    // Sternenbulge mit Hernquist-Profil
    double M_Bulge;
    double R_Bulge;
    int N_Bulge;
    double scaleRadius = 1e19; // Angepasster Skalierungsradius des Hernquist-Profils

    // Gas-Scheibe mit exponentiellem Profil
    double M_Gas;
    double R_Gas;
    double z_Gas;
    int N_Gas;    
    double velocity_dispersion = 1e4; // Geschwindigkeitsdispersion

    // Position und Geschwindigkeit der Galaxie
    vec3 galaxyPosition;
    vec3 galaxyVelocity;

    // Methoden zur Generierung der Galaxie
    void generateGalaxy();

private:
    // Zeiger auf die Teilchenliste
    std::vector<std::shared_ptr<Particle>>* particles;

    // Methoden zur Generierung einzelner Komponenten
    void generateStarDisk(int start, int end);
    void generateBulge(int start, int end);

    // Methoden zur Berechnung der eingeschlossenen Masse
    double enclosedMassBulge(double r) const;
    double enclosedMassDisk(double R) const;

    // Methoden zur Hernquist-Dichte und Potential
    double densityHernquist(double r) const;
    double PotentialHernquist(double r, double M, double a) const; // Als const deklariert

    // Zusätzliche private Methoden für die Bulge-Generierung
    double sample_radius_bulge(double Y) const;
    double potential_bulge(double r) const;
    double escape_velocity_bulge(double r) const;
    double f_energy(double E, double r) const;
    void sample_direction(double &ux, double &uy, double &uz, std::mt19937 &gen) const; // Korrigierte Parameter
};

#endif // GALAXY_H
