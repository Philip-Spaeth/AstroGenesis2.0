#include "TimeIntegration.h"

void TimeIntegration::Euler(std::shared_ptr<Particle> particle, double deltaTime)
{
    // Semi implicit Euler
    particle->velocity = particle->velocity + particle->acceleration * deltaTime;
    particle->position = particle->position + particle->velocity * deltaTime;
}

void TimeIntegration::Kick(std::shared_ptr<Particle> particle, double deltaTime)
{
    // Leapfrog Kick
    particle->velocity = particle->velocity + particle->acceleration * deltaTime / 2;
}

void TimeIntegration::Drift(std::shared_ptr<Particle> particle, double deltaTime)
{
    // Leapfrog Drift
    particle->position = particle->position + particle->velocity * deltaTime;
}

//integrate the entropy
void TimeIntegration::EntropyEuler(std::shared_ptr<Particle> particle, double deltaTime)
{
    // Entropy Euler
    particle->A = particle->A + particle->dAdt * deltaTime;
}