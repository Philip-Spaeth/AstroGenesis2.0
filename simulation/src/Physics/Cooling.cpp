#include "Cooling.h"
#include "Constants.h"
#include <cmath>
#include <iostream>

void Cooling::coolingRoutine(std::shared_ptr<Particle>& particle)
{
    // Constants in cgs units
    const double erg_to_joule = 1e-7;    // Convert erg to Joules
    const double cm3_to_m3 = 1e6;        // Convert cm^-3 to m^-3

    // Example densities (in cm^-3)
    double ne_cgs = 1e6;      // Electron number density (cm^-3)
    double nH_cgs = 1e6;      // Total hydrogen number density (cm^-3)
    double nHI_cgs = 0.8 * nH_cgs; // Neutral hydrogen number density (cm^-3)
    double nHe_cgs = 0.1 * nH_cgs; // Total helium number density (cm^-3)
    double nHeI_cgs = 0.09 * nHe_cgs; // Neutral helium number density (cm^-3)
    double nHeII_cgs = 0.01 * nHe_cgs; // Singly ionized helium number density (cm^-3)

    double T = particle->T; // Temperature (K)

    //cooling tables from Weinberg, Katz, Hernquist, 1996 table 1
    // Cooling rate in erg/s/cm^3 (cgs units)
    double coolingRate_cgs = 0.0;

    // Gaunt factor (typical value for T ~ 10^4 - 10^7 K)
    double gff = 1.1;

    // Free-Free Emission (all ions, cgs units)
    double free_free_cgs = 1.42e-27 * gff * sqrt(T) * ne_cgs * (nH_cgs - nHI_cgs + nHeI_cgs + 2.0 * nHeII_cgs);
    coolingRate_cgs += free_free_cgs;

    // Convert cooling rate to SI units (J/s/m^3)
    double coolingRate_SI = coolingRate_cgs * erg_to_joule * cm3_to_m3;

    if (coolingRate_SI > 0 && particle->rho > 0)
    {
        double dUdt_cooling = coolingRate_SI / particle->rho;
        particle->dUdt -= dUdt_cooling;
        std::cout << "Temperature (K): " << T << "\n";
        std::cout << "Electron Density (cm^-3): " << ne_cgs << "\n";
        std::cout << "Cooling Rate (cgs): " << coolingRate_cgs << " erg/s/cm^3\n";
        std::cout << "Cooling Rate (SI): " << coolingRate_SI << " J/s/m^3\n";
        std::cout << "Energy Change Rate (dU/dt): " << particle->dUdt << " J/kg/s\n";
    }
}