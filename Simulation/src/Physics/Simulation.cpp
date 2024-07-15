#include "Simulation.h"

Simulation::Simulation()
{ 
    particles.push_back(std::make_shared<Particle>(vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 0.0), 10000));
    particles.push_back(std::make_shared<Particle>(vec3(1.0, 0.0, 0.0), vec3(0.0, 2.0, 0.0), vec3(0.0, 0.0, 0.0), 10000));
    particles.push_back(std::make_shared<Particle>(vec3(1.0, 3.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.0, 0.0, 0.0), 10000));
    particles.push_back(std::make_shared<Particle>(vec3(1.0, 0.0, -3.0), vec3(0.0, 1.0, 0.0), vec3(0.0, 0.0, 0.0), 10000));
    particles.push_back(std::make_shared<Particle>(vec3(1.0, 0.0, 3.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 0.0), 10000));
}

Simulation::~Simulation()
{
    std::cout << particles[0]->position << std::endl;

    //calculate the gravitational acceleration for each particle
    for (int t = 0; t < timeSteps; t++)
    {
        for (int i = 0; i < numberOfParticles; i++)
        {
            for (int j = 0; j < numberOfParticles; j++)
            {
                if (i != j)
                {
                    //calculate the gravitational force
                    vec3 distance = particles[j]->position - particles[i]->position;
                    double r = distance.length() + softening;
                    vec3 normalized_direction = distance.normalize();
                    particles[i]->acceleration += normalized_direction * Constants::G * (particles[j]->mass / (r * r));
                }
            }

            //update the position and velocity of the particle
            timeIntegration->Euler(particles[i], deltaTime);
        }

        //reset the acceleration
        for (int i = 0; i < numberOfParticles; i++)
        {
            particles[i]->acceleration = vec3(0.0, 0.0, 0.0);
        }
    }
    
    std::cout << particles[0]->position << std::endl;
}

void Simulation::run()
{
    std::cout << "Simulation is running" << std::endl;
}