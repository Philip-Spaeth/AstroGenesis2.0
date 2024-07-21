#include "Simulation.h"

Simulation::Simulation()
{
    //write the info file
    dataManager->writeInfoFile(deltaTime, timeSteps, numberOfParticles);

    random::setRandomSeed(42);

    for(int i = 0; i < numberOfParticles; i++)
    {
        particles.push_back(std::make_shared<Particle>(vec3(random::between(-9,10), random::between(-10,10), random::between(-10,10)), vec3(random::between(-0,0), random::between(0,0), random::between(0,0)), vec3(0.0, 0.0, 0.0), 1000000));
    }
    
    //save the particles data
    dataManager->saveData(particles, 0);
}

Simulation::~Simulation(){}

void Simulation::run()
{

    //calculate the gravitational acceleration for each particle
    for (int t = 0; t < timeSteps; t++)
    {
        totalPotentialEnergy.push_back(0);
        totalKineticEnergy.push_back(0);
        totalInternalEnergy.push_back(0);
        totalEnergy.push_back(0);


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

                    //Energy conservation
                    particles[i]->potentialEnergy += 0.5 * -Constants::G * (particles[j]->mass * particles[i]->mass) / r;
                    totalPotentialEnergy[t] += particles[i]->potentialEnergy;

                }
            }
            //update the position and velocity of the particle
            timeIntegration->Euler(particles[i], deltaTime);

            //calculate the kinetic energy
            particles[i]->kineticEnergy = 0.5 * particles[i]->mass * (particles[i]->velocity.length() * particles[i]->velocity.length());
            //calculate the total energy
            totalKineticEnergy[t] += particles[i]->kineticEnergy;
            totalInternalEnergy[t] += particles[i]->internalEnergy;
            totalEnergy[t] += particles[i]->totalEnergy + particles[i]->potentialEnergy + particles[i]->kineticEnergy;
        }

        if(t == 0) {std::cout << "total energy in the begining: " << totalEnergy[t] << std::endl;}

        dataManager->printProgress(t, timeSteps);

        if(t == timeSteps - 1) std::cout << "total energy at the end: " << totalEnergy[t] << std::endl;
        if(t == timeSteps - 1) std::cout << "difference: " << (totalEnergy[t] - totalEnergy[0])<< std::endl;
        if(t == timeSteps - 1) std::cout << "difference in percentage: " << ((totalEnergy[t] - totalEnergy[0]) / totalEnergy[0] * 100)<< "%" << std::endl;
        //reset the acceleration
        for (int i = 0; i < numberOfParticles; i++)
        {
            particles[i]->acceleration = vec3(0.0, 0.0, 0.0);
            particles[i]->potentialEnergy = 0.0; // Reset potential energy
        }
        //save the particles data
        dataManager->saveData(particles, t + 1);
    }
}