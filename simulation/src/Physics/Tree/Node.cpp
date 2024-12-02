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
        children[i].reset();
    }
}


vec3 Node::calcSPHForce(std::shared_ptr<Particle> newparticle)
{
    vec3 acc = vec3(0,0,0);
    double h_i = newparticle->h;
    double h_j = mH;
    if(isLeaf) h_j = this->particle->h;
    h_j = h_i;
    double h_ij = (h_i + h_j) / 2.0;
    //if(h_i == 0 || h_j == 0) return vec3(0,0,0);

    double rho_i = newparticle->rho;
    double rho_j = mRho;
    if(isLeaf) rho_j = this->particle->rho;
    rho_j = rho_i;
    double rho_ij = (rho_i + rho_j) / 2.0;
    //if(rho_i == 0 || rho_j == 0) return vec3(0,0,0);

    double P_i = newparticle->P;
    double P_j = mP;
    if(isLeaf) P_j = this->particle->P;
    P_j = P_i;
    //if(P_i == 0 || P_j == 0) return vec3(0,0,0);

    vec3 v_i = newparticle->velocity;
    vec3 v_j = mVel;
    if(isLeaf) v_j = this->particle->velocity;
    vec3 v_ij = v_i - v_j;

    vec3 d = newparticle->position - centerOfMass;
    double r = d.length();
    double c_i = sqrt(Constants::GAMMA * P_i / rho_i);
    double c_j = sqrt(Constants::GAMMA * P_i / rho_i);
    double c_ij = (c_i + c_j) / 2.0;

    //Pressure force
    //Monaghan (1992)
    if(true)
    {
        acc += -gasMass * (P_i / (rho_i * rho_i) + P_j / (rho_j * rho_j)) * kernel::gradientCubicSplineKernel(d, h_i);
    }
    //Springel & Hernquist (2002)
    //entropy conservation formalism
    if(false)
    {
        double d_rho_dh_i = 1;
        double d_rho_dh_j = 1;
        double f_i = pow((1.0 + (h_i / (3 * rho_i)) * d_rho_dh_i), -1);
        double f_j = pow((1.0 + (h_j / (3 * rho_j)) * d_rho_dh_j), -1);
        acc += -gasMass * (f_i * P_i / (rho_i * rho_i) + f_j * P_j / (rho_j * rho_j)) * kernel::gradientCubicSplineKernel(d, h_i);
    }

    //Artificial viscosity
    //Monaghan & Gingold (1983)
    double MU_ij = 0.0;
    if(true)
    {
        double alpha = 0.5;
        double beta = 1.0;
        double eta = 0.01;
        double mu_ij = h_ij * v_ij.dot(d) / (r * r + eta * (h_ij * h_ij));
        if(v_ij.dot(d) < 0)
        {
            MU_ij = -alpha * c_ij * mu_ij + beta * (mu_ij * mu_ij);
        }
        acc += -gasMass * MU_ij * kernel::gradientCubicSplineKernel(d, h_ij);
    }
    //Monaghan (1997)
    if(false)
    {
        double alpha = 1.0;
        double w_ij = v_ij.dot(d) / (r);
        double v_sig_ij = (c_i + c_j - 3 * w_ij);
        MU_ij = -(alpha / 2.0) * (v_sig_ij * w_ij / rho_ij);
        acc += -gasMass * MU_ij * kernel::gradientCubicSplineKernel(d, h_ij);
    }

    //Internal energy
    newparticle->dUdt += 1.0 / 2.0 * gasMass * (P_i / (rho_i * rho_i) + P_j / (rho_j * rho_j) + MU_ij) * v_ij.dot(kernel::gradientCubicSplineKernel(d, h_i));

    if(std::isnan(acc.x) || std::isnan(acc.y) || std::isnan(acc.z)) return vec3(0,0,0);

    return acc;
}

void Node::calculateGravityForce(std::shared_ptr<Particle> newparticle, double softening, double theta)
{
    if (mass == 0) return;
    if (!newparticle) return;
    if (newparticle == this->particle) return;
    if (newparticle->mass == 0) return;

    vec3 d = centerOfMass - newparticle->position;
    double r = d.length();

    if(r == 0) return;

    if(isLeaf)
    {
        if (this->particle && newparticle != this->particle)
        {
            double e0 = softening;
            //softening described by Springel, Yoshida & White (2001) eq. 71
            double e = -(2.8 * e0) / (kernel::softeningKernel(r / (2.8 * e0)) - r);
            //gravity calculation
            vec3 gravityAcceleration = Constants::G * mass / (r * r + e * e) * d.normalize();
            newparticle->acceleration += gravityAcceleration;

            if(r < newparticle->h * 2)
            {
                //check if both are gas particles
                if(this->particle->type == 2 && newparticle->type == 2)
                {
                    newparticle->acceleration += calcSPHForce(newparticle);
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
            double e = -(2.8 * e0) / (kernel::softeningKernel(r / (2.8 * e0)) - r);
            //gravity calculation
            vec3 gravityAcceleration = Constants::G * mass / (r * r + e * e) * d.normalize();
            newparticle->acceleration += gravityAcceleration;
            //std::cout << e << std::endl;

            //SPH calculation
            if(r < newparticle->h * 2)
            {
                //check if both are gas particles
                if(newparticle->type == 2 && gasMass > 0)
                {
                    newparticle->acceleration += calcSPHForce(newparticle);
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
        //std::cerr << "Warning: Particles are outside the bounds and will not be properly considered" << std::endl;
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
        //calc m Node properties
        if (gasMass > 0) 
        {
            mVel = (mVel * (gasMass - newParticle->mass) + newParticle->velocity * newParticle->mass) / gasMass;
            /*if (newParticle->rho > 0) {
                mRho = (mRho * (gasMass - newParticle->mass) + newParticle->rho * newParticle->mass) / gasMass;
            }
            if (newParticle->P > 0) {
                mP = (mP * (gasMass - newParticle->mass) + newParticle->P * newParticle->mass) / gasMass;
            }
            if (newParticle->h > 0) {
                mH = (mH * (gasMass - newParticle->mass) + newParticle->h * newParticle->mass) / gasMass;
            }*/
        }
    }

    centerOfMass = (centerOfMass * (mass - newParticle->mass) + newParticle->position * newParticle->mass) / mass;
}


int Node::getOctant(std::shared_ptr<Particle> newParticle) 
{
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
    if(gasMass == 0) return;

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
            #pragma omp parallel for
            for (size_t i = 0; i < childParticles.size(); i++)
            {
                if(childParticles[i]->type == 2)
                {
                    childParticles[i]->h = radius * 2;
                }
            }

            //calculate the rho for all the particles in the node
            double rho = 0;
            #pragma omp parallel for
            for(size_t i = 0; i < childParticles.size(); i++)
            {
                if(childParticles[i]->type == 2)
                {
                    double drho = childParticles[i]->mass * kernel::cubicSplineKernel(vec3(childParticles[i]->position - centerOfMass).length(), childParticles[i]->h);
                    rho += drho;
                }
            }

            //add the new rho to all the particles in the node
            #pragma omp parallel for
            for(size_t i = 0; i < childParticles.size(); i++)
            {
                if(childParticles[i]->type == 2)
                {
                    childParticles[i]->rho = rho;
                    //calc P, P = (gamma-1)*u*rho
                    childParticles[i]->P = (Constants::GAMMA - 1.0) * childParticles[i]->U * childParticles[i]->rho;
                    //calc T, T = (gamma-1)*u*prtn / (bk)
                    childParticles[i]->T = (Constants::GAMMA - 1.0) * childParticles[i]->U * Constants::prtn * childParticles[i]->mu / (Constants::k_b);
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