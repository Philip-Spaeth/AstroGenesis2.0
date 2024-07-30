#include "TimeIntegration.h"

TimeIntegration::TimeIntegration()
{
}

TimeIntegration::~TimeIntegration()
{
}

void TimeIntegration::Euler(std::shared_ptr<Particle> particle, double deltaTime)
{
    // Semi implicit Euler
    particle->velocity = particle->velocity + particle->acceleration * deltaTime;
    //std::cout << "Velocity: " << particle->acceleration << std::endl;
    particle->position = particle->position + particle->velocity * deltaTime;

    //Periodic boundary conditions
    double box_size = 1e24;
    particle->position.x = fmod(particle->position.x + box_size, box_size);
    particle->position.y = fmod(particle->position.y + box_size, box_size);
    particle->position.z = fmod(particle->position.z + box_size, box_size);
}