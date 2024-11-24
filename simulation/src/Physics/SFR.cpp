#include "SFR.h"

const double rho_th = 1e-22; // density threshold in kg/m^3
const double t_star = 3.16e13; // star formation time in s
const double epsilon = 0.1; // efficiency of star formation
const double T_th = 1e4; // temperature threshold in K

void SFR::sfrRoutine(std::shared_ptr<Particle>& particle)
{
    // Sternentstehungsraten berechnen
    if (particle->rho > rho_th && particle->T < T_th)
    {
        particle->sfr = epsilon * (particle->rho) / t_star;
        double deltaM_star = particle->sfr * particle->timeStep;

        if(deltaM_star > particle->mass && deltaM_star > 0)
        {
            //gas -> star
            particle->type = 1;
            particle->U = 0.0;
            std::cout << "Star formed" << std::endl;
        }
    }
}