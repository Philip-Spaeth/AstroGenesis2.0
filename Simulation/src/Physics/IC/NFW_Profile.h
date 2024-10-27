// NFW_Profile.h
#ifndef NFW_PROFILE_H
#define NFW_PROFILE_H

#include "Constants.h"
#include <cmath>
#include <algorithm> // Für std::max und std::clamp

/**
 * @brief Klasse zur Beschreibung des NFW-Profils (Navarro-Frenk-White).
 */
class NFW_Profile
{
public:
    /**
     * @brief Berechnet die Dichte an einem gegebenen Radius r gemäß dem NFW-Profil.
     * 
     * @param r Radius, bei dem die Dichte berechnet werden soll (in Metern).
     * @param M Gesamtmasse des Halos (in Kilogramm).
     * @param Rs Skalierungsradius des Halos (in Metern).
     * @param c Konzentrationsparameter (c = R_vir / Rs).
     * @return double Dichte in kg/m³.
     */
    static double density(double r, double M, double Rs, double c)
    {
        if (r <= 0.0) return 0.0; // Vermeidung von Division durch Null

        // Berechnung der charakteristischen Dichte rho0
        double rho0 = M / (4.0 * Constants::PI * std::pow(Rs, 3) * (std::log(1.0 + c) - c / (1.0 + c)));

        return rho0 / ((r / Rs) * std::pow(1.0 + r / Rs, 2));
    }
    
    /**
     * @brief Berechnet das Gravitationspotential an einem gegebenen Radius r gemäß dem NFW-Profil.
     * 
     * @param r Radius, bei dem das Potential berechnet werden soll (in Metern).
     * @param M Gesamtmasse des Halos (in Kilogramm).
     * @param Rs Skalierungsradius des Halos (in Metern).
     * @param c Konzentrationsparameter (c = R_vir / Rs).
     * @return double Potential in m²/s².
     */
    static double potential(double r, double M, double Rs, double c)
    {
        if (r <= 0.0) return 0.0; // Potential bei r=0 ist unendlich, praktisch 0 setzen

        return -Constants::G * M / r * std::log(1.0 + r / Rs);
    }
    
    /**
     * @brief Berechnet die eingeschlossene Masse innerhalb eines gegebenen Radius r gemäß dem NFW-Profil.
     * 
     * @param r Radius, innerhalb dessen die Masse berechnet werden soll (in Metern).
     * @param M Gesamtmasse des Halos (in Kilogramm).
     * @param Rs Skalierungsradius des Halos (in Metern).
     * @param c Konzentrationsparameter (c = R_vir / Rs).
     * @return double Eingeschlossene Masse in Kilogramm.
     */
    static double enclosedMass(double r, double M, double Rs, double c)
    {
        if (r <= 0.0) return 0.0;

        return M * (std::log(1.0 + r / Rs) - (r / Rs) / (1.0 + r / Rs)) / (std::log(1.0 + c) - c / (1.0 + c));
    }
    
    /**
     * @brief Generiert einen Radius r gemäß dem NFW-Profil mittels der inversen CDF-Methode.
     * 
     * @param U Zufallszahl zwischen 0 und 1.
     * @param M Gesamtmasse des Halos (in Kilogramm).
     * @param Rs Skalierungsradius des Halos (in Metern).
     * @param c Konzentrationsparameter (c = R_vir / Rs).
     * @param R_vir Virialradius des Halos (in Metern).
     * @return double Generierter Radius r (in Metern).
     */
    static double sample_radius(double U, double M, double Rs, double c, double R_vir)
    {
        // Funktion zur Berechnung der kumulativen Massenverteilung
        auto cumulativeMass = [&](double r_sample) -> double {
            return enclosedMass(r_sample, M, Rs, c) / M;
        };
        
        // Ziel: F(r) = U, wo F(r) die kumulative Verteilungsfunktion ist
        // Numerische Inversion mittels Newton-Raphson Methode

        // Anfangsschätzung für r
        double r = Rs * c * U; // Eine grobe Anfangsschätzung

        // Toleranz und maximale Iterationen
        const double tol = 1e-6;
        const int max_iter = 100;

        for(int iter = 0; iter < max_iter; ++iter)
        {
            double F = cumulativeMass(r);
            double F_diff = F - U;

            if(std::abs(F_diff) < tol)
                break;

            // Berechne die Ableitung dF/dr
            double dFdr = density(r, M, Rs, c) / M;

            // Sicherstellen, dass dFdr nicht zu klein wird, um Division durch Null zu vermeiden
            dFdr = std::max(dFdr, 1e-20);

            // Update r
            r -= F_diff / dFdr;

            // Sicherstellen, dass r innerhalb der Grenzen bleibt
            r = std::clamp(r, 0.0, R_vir);
        }

        return r;
    }
};

#endif // NFW_PROFILE_H
