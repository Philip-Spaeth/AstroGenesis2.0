#include "Disk.h"
#include "Constants.h"
#include "Exponential_Profile.h"
#include "Hernquist_Profile.h"
#include "NFW_Profile.h"
#include <random>
#include <memory>

void Disk::generateStarDisk(int start, int end)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist_theta(0.0, 2.0 * Constants::PI);
    std::uniform_real_distribution<> dist_uniform(0.0, 1.0);
    std::normal_distribution<> dist_z(0.0, g->z_Disk);
    std::normal_distribution<> dist_velocity(0.0, g->VelDis_Disk);

    double R_d = g->R_Disk / 5.0;  // Skalierungslänge der Scheibe
    double total_mass = g->M_Disk;
    double mass_per_particle = total_mass / static_cast<double>(g->N_Disk);

    for (int i = start; i < end; ++i)
    {
        double R;
        do {
            double rand_num = dist_uniform(gen);
            R = -R_d * log(1.0 - rand_num); // Exponentialverteilung für Radius
        } while (R > g->R_Disk); // Wiederhole, bis ein gültiger Radius gefunden wird

        double theta = dist_theta(gen);
        double z = dist_z(gen);
        double x = R * cos(theta);
        double y = R * sin(theta);
        vec3 position = g->galaxyPosition + vec3(x, y, z);


        double M_enc_bulge = 0;
        //if(g->N_Bulge != 0) M_enc_bulge = Hernquist_Profile::enclosedMassHernquist(R, g->M_Bulge, g->R_Bulge);
        double M_enc_Halo = 0;
        //if(g->N_Halo != 0) M_enc_Halo = NFW_Profile::enclosedMass(R, g->M_Halo, g->c_Halo, g->c_Halo);
        double M_enc_disk = Exponential_Profile::enclosedMassExponential(R, g->M_Disk, g->R_Disk, 10);
        double v_theta = sqrt(Constants::G * (M_enc_disk + M_enc_bulge + M_enc_Halo) / (R + 1e-10)); // Vermeide Division durch Null

        double v_r = dist_velocity(gen);
        double v_z = dist_velocity(gen) * (g->z_Disk / g->R_Disk);
        double v_phi = v_theta * (1 + (dist_velocity(gen) / v_theta)); // Tangentialgeschwindigkeit

        double vx = -v_phi * sin(theta) + v_r * cos(theta);
        double vy = v_phi * cos(theta) + v_r * sin(theta);
        double vz = v_z;

        vec3 velocity = g->galaxyVelocity + vec3(vx, vy, vz);

        std::shared_ptr<Particle> p = std::make_shared<Particle>();
        g->particles->at(i) = p;
        p->position = position;
        p->velocity = velocity;
        p->mass = mass_per_particle;
        p->type = 1;  // Stern
    }
}
