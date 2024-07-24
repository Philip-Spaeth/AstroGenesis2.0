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
}