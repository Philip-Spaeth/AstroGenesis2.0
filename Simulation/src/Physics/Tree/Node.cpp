#include "Node.h"
#include "Constants.h"
#include "kernel.h"
#include <algorithm>


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
    parent = std::weak_ptr<Node>();
}

Node::~Node()
{
    for (int i = 0; i < 8; i++)
    {
        children[i].reset(); // Explicitly reset the children to release memory
    }
}

void Node::calcMedianPressure()
{
    if(gasMass == 0) return;

    //calculate the median smoothing length of the child particles
    std::vector<double> pValues;
    for(size_t i = 0; i < childParticles.size(); i++)
    {
        if(childParticles[i] == nullptr) continue;
        if(childParticles[i]->type == 2)
        {
            pValues.push_back(childParticles[i]->P);
        }
    }
    std::sort(pValues.begin(), pValues.end());
    if(pValues.size() == 0) return;
    medianPressure = pValues[pValues.size() / 2];

    for (int i = 0; i < 8; i++)
    {
        if (children[i] != nullptr)
        {
            if(children[i]->gasMass != 0)
            {
                children[i]->calcMedianPressure();
            }
        }
    }
}

void Node::calcMedianVelocity()
{
    if(gasMass == 0) return;

    //calculate the median smoothing length of the child particles
    vec3 velocitySum = vec3(0,0,0);
    int count = 0;
    for(size_t i = 0; i < childParticles.size(); i++)
    {
        if(childParticles[i] == nullptr) continue;
        if(childParticles[i]->type == 2)
        {
            velocitySum += childParticles[i]->velocity;
            count++;
        }
    }
    if(count == 0) return;
    medianVelocity = velocitySum / count;
}

void Node::calcSPHNodeMedians()
{
    calcMedianH();
    calcMedianDensity();
    calcMedianVelocity();

    for (int i = 0; i < 8; i++)
    {
        if (children[i] != nullptr)
        {
            if(children[i]->gasMass != 0)
            {
                children[i]->calcSPHNodeMedians();
            }
        }
    }
}

void Node::calcMedianH()
{
    //calculate the median smoothing length of the child particles
    std::vector<double> hValues;
    for(size_t i = 0; i < childParticles.size(); i++)
    {
        if(childParticles[i] == nullptr) continue;
        if(childParticles[i]->type == 2)
        {
            hValues.push_back(childParticles[i]->h);
        }
    }
    std::sort(hValues.begin(), hValues.end());
    if(hValues.size() == 0) return;
    medianH = hValues[hValues.size() / 2];
}

void Node::calcMedianDensity()
{
    //calculate the median density of the child particles
    std::vector<double> densityValues;
    for(size_t i = 0; i < childParticles.size(); i++)
    {
        if(childParticles[i] == nullptr) continue;
        if(childParticles[i]->type == 2)
        {
            densityValues.push_back(childParticles[i]->rho);
        }
    }
    std::sort(densityValues.begin(), densityValues.end());
    if(densityValues.size() == 0) return;
    medianDensity = densityValues[densityValues.size() / 2];
}


void Node::calculateGravityForce(std::shared_ptr<Particle> newparticle, double softening, double theta)
{
    if (mass == 0) return; // Verhindert Berechnung, wenn keine Masse vorhanden ist
    if (!newparticle) return; // Verhindert Berechnung, wenn kein Teilchen vorhanden ist
    if (newparticle == this->particle) return; // Verhindert Berechnung, wenn es sich um das gleiche Teilchen handelt
    if (newparticle->mass == 0) return; // Verhindert Berechnung, wenn es sich um das gleiche Teilchen handelt

    vec3 d = centerOfMass - newparticle->position;
    double r = d.length();

    if(r == 0) return;

    if(isLeaf)
    {
        if (this->particle && newparticle != this->particle)
        {
            double e0 = softening;
            //softening described by Springel, Yoshida & White (2001) eq. 71
            double e = -(2.8 * e0) / kernel::softeningKernel(r / (2.8 * e0)) - r;
            //gravity calculation
            vec3 gravityAcceleration = Constants::G * mass / (r * r + e * e) * d.normalize();
            newparticle->acceleration += gravityAcceleration;
            //std::cout << e << std::endl;

            //SPH calculation
            if(r < newparticle->h * 2)
            {
                //check if both are gas particles
                if(this->particle->type == 2 && newparticle->type == 2)
                {
                    double h_i = newparticle->h;
                    if(h_i == 0) return;
                    double h_j = this->particle->h;
                    if(h_j == 0) return;
                    double h_ij = (h_i + h_j) / 2;

                    double rho_i = newparticle->rho;
                    if(rho_i == 0) return;
                    double rho_j = this->particle->rho;
                    if(rho_j == 0) return;
                    double rho_ij = (rho_i + rho_j) / 2;

                    double P_i = newparticle->P;
                    if(P_i == 0) return;
                    double P_j = this->particle->P;
                    if(P_j == 0) return;

                    // calculate the acceleration due to the pressure force
                    vec3 grad_i = kernel::gradientCubicSplineKernel(d, h_i);
                    vec3 grad_j = kernel::gradientCubicSplineKernel(d, h_j);
                    vec3 grad_ij = (grad_i + grad_j) / 2;
        
                    // calculate the acceleration due to the pressure force
                    vec3 pressureAcceleration = - newparticle->mass * (P_i / (rho_i * rho_i) + P_j / (rho_j * rho_j)) * grad_i;
                    //std::cout << pressureAcceleration << std::endl; 
                    //newparticle->acceleration += pressureAcceleration;

                    //Artificial viscosity
                    double c_i = sqrt(Constants::GAMMA * P_i / rho_i);
                    double c_j = sqrt(Constants::GAMMA * P_j / rho_j);
                    double c_ij = (c_i + c_j) / 2;
                    vec3 v_ij = newparticle->velocity - this->particle->velocity;
                    double mu_ij = (h_ij * v_ij.dot(d)) / (pow(r,2) + 0.01 * pow(h_ij,2));
                    double MU_ij = 0;
    	            if(v_ij.dot(d) < 0)
                    {
                        double alpha = 0.5;
                        double beta = 1;
                        MU_ij = (-alpha * c_ij * mu_ij + beta * pow(mu_ij, 2)) / (rho_ij);
                    }
                    
                    // calculate the acceleration due to the artificial viscosity
                    vec3 viscosityAcceleration = -newparticle->mass * MU_ij * grad_ij * 1;
                    //newparticle->acceleration += viscosityAcceleration;
                    //calc the change of entropy 
                    //newparticle->dAdt += 0.5 * (Constants::GAMMA - 1) / (pow(rho_i, Constants::GAMMA -1)) * gasMass * mu_ij * v_ij.dot(grad_ij);
                    
                }
            }
        }
    }
    else
    {
  
        double s = radius / r;
        if (s < theta)
        {
            double e0 = softening;
            //softening described by Springel, Yoshida & White (2001) eq. 71
            double e = -(2.8 * e0) / kernel::softeningKernel(r / (2.8 * e0)) - r;
            //gravity calculation
            vec3 gravityAcceleration = Constants::G * mass / (r * r + e * e) * d.normalize();
            newparticle->acceleration += gravityAcceleration;
            //std::cout << e << std::endl;

            //SPH calculation
            if(r < newparticle->h * 2)
            {
                //check if both are gas particles
                if(newparticle->type == 2 && gasMass != 0)
                {
                    double h_i = newparticle->h;
                    if(h_i == 0) return;
                    double h_j = medianH;
                    if(h_j == 0) return;
                    double h_ij = (h_i + h_j) / 2;

                    double rho_i = newparticle->rho;
                    if(rho_i == 0) return;
                    double rho_j = medianDensity;
                    if(rho_j == 0) return;
                    double rho_ij = (rho_i + rho_j) / 2;

                    double P_i = newparticle->P;
                    if(P_i == 0) return;
                    double P_j = medianPressure;
                    if(P_j == 0) return;

                    // calculate the acceleration due to the pressure force
                    vec3 grad_i = kernel::gradientCubicSplineKernel(d, h_i);
                    vec3 grad_j = kernel::gradientCubicSplineKernel(d, h_j);
                    vec3 grad_ij = (grad_i + grad_j) / 2;
                    
                    // calculate the acceleration due to the pressure force
                    vec3 pressureAcceleration = gasMass * (P_i / (rho_i * rho_i) + P_j / (rho_j * rho_j)) * grad_i;
                    //newparticle->acceleration += pressureAcceleration;
                    
                    //Artificial viscosity
                    double c_i = sqrt(Constants::GAMMA * P_i / rho_i);
                    double c_j = sqrt(Constants::GAMMA * P_j / rho_j);
                    double c_ij = (c_i + c_j) / 2;
                    vec3 v_ij = newparticle->velocity - medianVelocity;
                    double mu_ij = (h_ij * v_ij.dot(d)) / (pow(r,2) + 0.01 * pow(h_ij,2));
                    double MU_ij = 0;
    	            if(v_ij.dot(d) < 0)
                    {
                        //calculate the artificial viscosity
                        double alpha = 0.5;
                        double beta = 1;
                        MU_ij = (-alpha * c_ij * mu_ij + beta * pow(mu_ij, 2)) / (rho_ij);
                    }
                    // calculate the acceleration due to the artificial viscosity
                    vec3 viscosityAcceleration = -gasMass * MU_ij * grad_ij * 1;
                    //newparticle->acceleration += viscosityAcceleration;
                    //calc the change of entropy 
                    //newparticle->dAdt += 0.5 * (Constants::GAMMA - 1) / (pow(rho_i, Constants::GAMMA -1)) * gasMass * mu_ij * v_ij.dot(grad_ij);
                    
                }
            }
        }
        else
        { 
             
            for (int i = 0; i < 8; i++)
            {
                if (children[i]->mass != 0)
                {
                    children[i]->calculateGravityForce(newparticle, softening, theta);
                }
            }
        }
    }
}

void Node::insert(std::shared_ptr<Particle> newParticle) 
{
    if (!newParticle) 
    {
        std::cerr << "Error: Particle to insert is null" << std::endl;
        return;
    }

    // Check if the particle is within the bounds of the current node
    if (newParticle->position.x < position.x - radius || newParticle->position.x > position.x + radius ||
        newParticle->position.y < position.y - radius || newParticle->position.y > position.y + radius ||
        newParticle->position.z < position.z - radius || newParticle->position.z > position.z + radius) 
    {
        std::cerr << "Warning: Particles are outside the bounds and will not be properly considered" << std::endl;
        return;
    }
    
    //add the particle to the childParticles vector
    childParticles.push_back(newParticle);

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
            for (int i = 0; i < 8; i++) 
            {
                children[i] = std::make_shared<Node>();
                // setup the node properties
                children[i]->position = position + vec3(
                    radius * (i & 1 ? 0.5 : -0.5),
                    radius * (i & 2 ? 0.5 : -0.5),
                    radius * (i & 4 ? 0.5 : -0.5));
                // std::cout << children[i]->position << std::endl;
                children[i]->radius = radius / 2;
                children[i]->depth = depth + 1;
                children[i]->parent = shared_from_this();
            }

            // get the octant of the particle
            int octant = getOctant(particle);
            
            // insert the particle in the corresponding octant
            if (octant != -1) 
            {
                children[octant]->insert(particle);
            }
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
    if(newParticle->type == 2)
    {
        gasMass += newParticle->mass;
    }
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

//consider only the gas particles and th gasMass
void Node::calcGasDensity(double massInH)
{
    //check if the parent ist not expired
    if(parent.lock() != nullptr)
    {
        //go upwards the tree until the mass of the node is closest to  massInH
        if(gasMass <  massInH && parent.lock())	//if the mass of the node is smaller than massInH and the node has a parent
        {
            //stop if the current gasMass is closer to massInH than the parent gasMass
            double massDifference =  massInH - gasMass;
            double parentMassDifference =  massInH - parent.lock()->gasMass;
            if(std::abs(massDifference) > std::abs(parentMassDifference))
            {
                parent.lock()->calcGasDensity( massInH);
            }
        }
        //stop if the current gasMass is closer to  massInH than the parent gasMass
        double massDifference =  massInH - gasMass;
        double parentMassDifference =  massInH - parent.lock()->gasMass;
        if(std::abs(massDifference) < std::abs(parentMassDifference) && gasMass != 0)
        {
            //calculate h for all the child particles in the node
            for (size_t i = 0; i < childParticles.size(); i++)
            {
                if(childParticles[i]->type == 2)
                {
                    childParticles[i]->h = radius * 2;
                }
            }

            //calculate the rho for all the particles in the node
            double rho = 0;
            for(size_t i = 0; i < childParticles.size(); i++)
            {
                if(childParticles[i]->type == 2)
                {
                    double drho = childParticles[i]->mass * kernel::cubicSplineKernel(vec3(childParticles[i]->position - centerOfMass).length(), childParticles[i]->h);
                    rho += drho;
                    
                }
            }

            //add the new rho to all the particles in the node
            for(size_t i = 0; i < childParticles.size(); i++)
            {
                if(childParticles[i]->type == 2)
                {
                    childParticles[i]->rho = rho;
                }
            }
        }
    }
}

// for the other particles just for visualization
void Node::calcVisualDensity(double radiusDensityEstimation) 
{
    double radiusDifference = radiusDensityEstimation - radius;
    double parentRadiusDifference = radiusDensityEstimation - parent.lock()->radius;
    if(std::abs(radiusDifference) > std::abs(parentRadiusDifference) && parent.lock())
    {
        parent.lock()->calcVisualDensity(radiusDensityEstimation);
    }
    else
    {
        double volume = radius * radius * radius;
        double density = mass / volume;
        for(size_t i = 0; i < childParticles.size(); i++)
        {
            childParticles[i]->visualDensity = density;
        }
    }
}