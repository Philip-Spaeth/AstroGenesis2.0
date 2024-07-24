#include "Node.h"
#include "Constants.h"

Node::Node()
{
    mass = 0.0;
    centerOfMass = vec3(0.0, 0.0, 0.0);
    isLeaf = true;
    particle = nullptr;
    for (int i = 0; i < 8; i++)
    {
        children[i] = nullptr;
    }
    parent = std::weak_ptr<Node>(); // Set weak_ptr to empty state
}

Node::~Node()
{
    for (int i = 0; i < 8; i++)
    {
        children[i].reset(); // Explicitly reset the children to release memory
    }
}

void Node::calculateForce(std::shared_ptr<Particle> newparticle, double softening, double theta)
{
    if (mass == 0) return; // Verhindert Berechnung, wenn keine Masse vorhanden ist
    if (!newparticle) return; // Verhindert Berechnung, wenn kein Teilchen vorhanden ist
    if (newparticle == this->particle) return; // Verhindert Berechnung, wenn es sich um das gleiche Teilchen handelt
    if (newparticle->mass == 0) return; // Verhindert Berechnung, wenn es sich um das gleiche Teilchen handelt

    vec3 d = centerOfMass - newparticle->position;
    double r = d.length();
    if(r == 0) return;

    if (isLeaf)
    {
        if (this->particle && newparticle != this->particle)
        {
            double newAcceleration = Constants::G * mass / ((r * r) + (softening * softening));
            newparticle->acceleration = newparticle->acceleration + d.normalize() * newAcceleration;
        }
    }
    else
    {
  
        double s = radius / r;
        if (s < theta)
        {
            double newAcceleration = Constants::G * mass / ((r * r) + (softening * softening));
            newparticle->acceleration += d.normalize() * newAcceleration;
        }
        else
        { 
             
            for (int i = 0; i < 8; i++)
            {
                if (children[i]->mass != 0)
                {
                    children[i]->calculateForce(newparticle, softening, theta);
                }
            }
        }
    }
}

void Node::insert(std::shared_ptr<Particle> newParticle) 
{
    // Check if the particle is within the bounds of the current node
    if (newParticle->position.x < position.x - radius || newParticle->position.x > position.x + radius ||
        newParticle->position.y < position.y - radius || newParticle->position.y > position.y + radius ||
        newParticle->position.z < position.z - radius || newParticle->position.z > position.z + radius) {
        return; // Particle is outside the bounds, do not insert
    }

    // if the node is a leaf node
    if (isLeaf) 
    {
        // if the node has no particle
        if (!particle) {
            // set the particle to the node
            particle = newParticle;
            particle->node = shared_from_this();
        }
        // if the node has a particle
        else 
        {
            // create the children of the node
            for (int i = 0; i < 8; i++) {
                children[i] = std::make_shared<Node>();
                // setup the node properties
                children[i]->position = position + vec3(
                    radius * (i & 1 ? 0.5 : -0.5),
                    radius * (i & 2 ? 0.5 : -0.5),
                    radius * (i & 4 ? 0.5 : -0.5));
                // std::cout << children[i]->position << std::endl;
                children[i]->radius = radius / 2;
                children[i]->depth = depth + 1;
                parent = shared_from_this();
            }
            // get the octant of the particle
            int octant = getOctant(particle);
            // insert the particle in the corresponding octant
            children[octant]->insert(particle);
            // get the octant of the new particle
            // insert the particle in the corresponding octant
            // get the octant of the particle
            octant = getOctant(newParticle);
            if (octant != -1) 
            {
                children[octant]->insert(newParticle);
            }

            // set the node as not a leaf node
            isLeaf = false;
            // set the particle to null
            particle.reset();
        }
    }
    // if the node is not a leaf node
    else 
    {
        // get the octant of the particle
        int octant = getOctant(newParticle);
        if (octant != -1) 
        {
            children[octant]->insert(newParticle);
        }
    }
    
    // Update the center of mass and mass
    mass += newParticle->mass;
    centerOfMass = (centerOfMass * (mass - newParticle->mass) + newParticle->position * newParticle->mass) / mass;

}


int Node::getOctant(std::shared_ptr<Particle> newParticle) {
    if (newParticle->position.x < position.x - radius || newParticle->position.x > position.x + radius ||
        newParticle->position.y < position.y - radius || newParticle->position.y > position.y + radius ||
        newParticle->position.z < position.z - radius || newParticle->position.z > position.z + radius) {
        std::cout << "Particle is outside the bounds" << std::endl;
        return -1; // Particle is outside the bounds
    }

    int octant = 0;
    if (newParticle->position.x > position.x) octant |= 1;
    if (newParticle->position.y > position.y) octant |= 2;
    if (newParticle->position.z > position.z) octant |= 4;
    return octant;
}


double Node::calcDensity(double h) 
{
    std::shared_ptr<Node> currentNode = shared_from_this();
    double targetRadius = 2 * h;

    // Traverse up the tree until the node's radius is closest to 2 * h
    while (!currentNode->parent.expired()) {
        auto parentPtr = currentNode->parent.lock();
        if (std::abs(currentNode->radius - targetRadius) <= std::abs(parentPtr->radius - targetRadius)) {
            break;
        }
        currentNode = parentPtr;
    }

    // Calculate density: density = mass / volume
    // Volume of the node is 8 * radius^3 for an octree
    double volume = 8 * std::pow(currentNode->radius, 3);
    double density = currentNode->mass / volume;

    return density;
}