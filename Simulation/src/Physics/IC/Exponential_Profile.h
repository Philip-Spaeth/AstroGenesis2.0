#pragma once
#include <cmath>
#include "Constants.h"

class Exponential_Profile
{
public:
    static double enclosedMassExponential(double R, double M_Disk, double R_Disk, double scale)
    {
        // Skalierungslänge R_d
        double R_d = R_Disk / scale;  // Wird oft so gewählt
        // Berechnung von Sigma_0 basierend auf der Gesamtmasse
        double Sigma_0 = M_Disk / (2 * Constants::PI * R_d * R_d);
        // Berechnung der eingeschlossenen Masse innerhalb Radius R
        double mass_enclosed = 2 * Constants::PI * Sigma_0 * R_d * R_d * 
                               (1 - exp(-R / R_d) * (1 + R / R_d));
        return mass_enclosed;
    }
};
