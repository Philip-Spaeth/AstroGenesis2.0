#pragma once
#include "Galaxy.h"

class Bulge
{
    
public:
    Bulge(Galaxy * g) : g(g) {}
    ~Bulge() {}

    Galaxy * g;

    void generateBulge(int start, int end);

    // Zusätzliche private Methoden für die Bulge-Generierung
    double sample_radius_bulge(double Y) const;
    double potential_bulge(double r) const;
    double escape_velocity_bulge(double r) const;
    double f_energy(double E, double r) const;
    void sample_direction(double &ux, double &uy, double &uz, std::mt19937 &gen) const; // Korrigierte Parameter
};