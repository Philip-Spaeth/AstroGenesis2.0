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
    void generateHernquistHalo(int start, int end, std::vector<std::shared_ptr<Particle>>& particles);
    void generateNFWHalo(int start, int end, std::vector<std::shared_ptr<Particle>>& particles);

//properties:
    //General:
    double M = 1.0e12 * Units::MSUN;

    //for Hernquist profile
    double a = 10.0 * Units::KPC;
    double maxR = 400 * Units::KPC;

    //for NFW profile
    double rs = 10.0 * Units::KPC;;    // Skalenradius
    double c = 10;     // Konzentrationsparameter
    double rho0;  // Zentrale Dichte (berechnet)
    double R_vir = rs * c;                    // Virialradius


private:
    //Hernquist profile
    double hernquistDensity(double r) const;
    double hernquistPotential(double r) const;
    double hernquistEnclosedMass(double r) const;
    double hernquistVelocityDispersion(double r) const;

    //NFW profile
    void calculateRho0(); // Berechnung von rho0
    double nfwDensity(double r) const;
    double nfwPotential(double r) const;
    double nfwEnclosedMass(double r) const;
    double sampleRadius() const;
    double sampleVelocity(double r) const;
};
