#pragma once
#include "Constants.h"
#include <cmath>

class Hernquist_Profile
{
public:
    static double enclosedMassHernquist(double r, double M_Bulge, double R_Bulge)
    {
        // Hernquist-Profil: M_enc = M_Bulge * r^2 / (R_Bulge + r)^2
        return M_Bulge * (r * r) / ((R_Bulge + r) * (R_Bulge + r));
    };

    // Methoden zur Hernquist-Dichte und Potential
    static double densityHernquist(double r, double M_Bulge, double Rs_Bulge)
    {
        return (M_Bulge / (2.0 * Constants::PI)) * (Rs_Bulge) / (r * std::pow(r + Rs_Bulge, 3));
    };
    static double PotentialHernquist(double r, double M, double a)
    {
        return -Constants::G * M / (r + a);
    };
};