#pragma once
#include <vector>
#include <memory>
#include "Particle.h"
#include <random>
#include <cmath>
#include <functional>
#include "Constants.h"
#include "Units.h"

class Halo
{
public:
    Halo();
    ~Halo();
    void generateHalo(int start, int end, std::vector<std::shared_ptr<Particle>>& particles);

    //properties:
    double M = 1.0e12 * Units::MSUN;
    double a = 10.0 * Units::KPC;

private:
    //Hernquist profile
    double hernquistDensity(double r) const;
    double hernquistPotential(double r) const;
    double hernquistEnclosedMass(double r) const;
    double hernquistVelocityDispersion(double r) const;
    double hernquistVelocity(double r) const;


//     //NFW profile
//     double NFWDensity(double r, double rs, double rho0){};
//     double NFWPotential(double r, double rs, double rho0){};
//     double NFWEnclosedMass(double r, double rs, double rho0){};

//     //Plummer profile
//     double plummerDensity(double r, double a, double M){};
//     double plummerPotential(double r, double a, double M){};
//     double plummerEnclosedMass(double r, double a, double M){};

//     //Einasto profile
//     double einastoDensity(double r, double a, double rho0, double alpha){};
//     double einastoPotential(double r, double a, double rho0, double alpha){};
//     double einastoEnclosedMass(double r, double a, double rho0, double alpha){};
};
