#include "Galaxy.h"
#include "random.h"
#include "Constants.h"


Galaxy::Galaxy(std::vector<std::shared_ptr<Particle>>& particles, int start, int end, double mass, double radius)
{
    random::setRandomSeed(37293);
    std::cout << "\nCreating Galaxy ..." << std::endl;
    int numberOfParticles = end - start;
    for(int i = start; i < end; i++)
    {
        particles.push_back(std::make_shared<Particle>());
        particles[i]->mass = mass / (numberOfParticles);

        double r = (double)((double)i * radius) / (double)numberOfParticles;
        //random position with the radius
        double angle = random::between(0, 2 * Constants::PI);
        double x = r * cos(angle);
        double y = r * sin(angle);
        particles[i]->position = vec3(x, y, random::between(-(radius / 10), (radius / 10)));

        //velocity around the center
        double massInRadius = mass / radius;
        double v = sqrt(Constants::G * (massInRadius));
        double vx = -v * sin(angle);
        double vy = v * cos(angle);
        particles[i]->velocity = vec3(vx, vy, 0);
    }
    std::cout << "Galaxy created\n" << std::endl;
}