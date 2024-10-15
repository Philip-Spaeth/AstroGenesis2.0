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
    double VelDis_Disk;
    int N_Disk;

    // Sternenbulge mit Hernquist-Profil
    double M_Bulge;
    double R_Bulge;
    int N_Bulge;
    double Rs_Bulge;

    // Gas-Scheibe mit exponentiellem Profil
    double M_Gas;
    double R_Gas;
    double z_Gas;
    int N_Gas;

    // Position und Geschwindigkeit der Galaxie
    vec3 galaxyPosition;
    vec3 galaxyVelocity;

    // Methoden zur Generierung der Galaxie
    void generateGalaxy();

    // Zeiger auf die Teilchenliste
    std::vector<std::shared_ptr<Particle>>* particles;
};

#endif // GALAXY_H
