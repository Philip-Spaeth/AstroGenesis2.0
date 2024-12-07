#include "SFR.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "Constants.h"

const double rho_th = 1e-22; // density threshold in kg/m^3
const double epsilon = 0.1; // efficiency of star formation
const double T_th = 1e4; // temperature threshold in K

void SFR::sfrRoutine(std::shared_ptr<Particle>& particle)
{
    // Sternentstehungsraten berechnen
    if (particle->rho > rho_th && particle->T < T_th)
    {
        double t_star = 1e15;//(1 / epsilon) * std::sqrt( (3* Constants::PI) / (32 * Constants::G * particle->rho) );
        double p = 1 - exp(-epsilon * particle->timeStep / t_star);
        double r = ((double) rand()) / RAND_MAX;
        particle->sfr = p;
        if (r < p)
        {
            //gas -> star
            particle->type = 1;
            particle->U = 0.0;
            //std::cout << "Star formed" << std::endl;
        }
    }
}