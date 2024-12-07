#include "Cooling.h"
#include "Constants.h"
#include <cmath>
#include <iostream>

void Cooling::coolingRoutine(std::shared_ptr<Particle>& particle)
{
    double T = particle->T;

    double ne_cgs = 1e6; // Elektronendichte in cm^-3
    double ni_cgs = 1e6; // Ionnendichte in cm^-3
    double gff = 1.1;

    // Free free Emission (erg/s/cm³)
    double free_free_cgs = 1.42e-27 * gff * sqrt(T) * ne_cgs * ni_cgs;

    // Gesamte Kühlrate in erg/s/cm³
    double coolingRate_cgs = free_free_cgs;
    double coolingRate_SI = coolingRate_cgs * 1e-7;

    if (coolingRate_SI > 0 && particle->rho > 0)
    {
        particle->dUdt -= coolingRate_SI / particle->rho;
        //std::cout << std::fixed << std::scientific << " rho: " << particle->rho <<  "coolingRate: " << coolingRate_SI << " dUdt: " << dUdt_cooling << std::endl;
    }
}