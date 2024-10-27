#include "Galaxy.h"
#include "Constants.h"
#include <iostream>
#include <random>
#include <cmath>
#include <memory>
#include "Bulge.h"
#include "Disk.h"
#include "Halo.h"

// Konstruktor und Destruktor
Galaxy::Galaxy(std::vector<std::shared_ptr<Particle>>* particles_ptr)
    : particles(particles_ptr) {}

Galaxy::~Galaxy() {}

// Implementierung der generateGalaxy-Funktion
void Galaxy::generateGalaxy()
{
    particles->clear();
    particles->reserve(N_Halo + N_Disk + N_Bulge + N_Gas);
    particles->resize(N_Halo + N_Disk + N_Bulge + N_Gas);

    Disk disk(this);
    disk.generateStarDisk(0, N_Disk);
    disk.generateGasDisk(N_Disk, N_Disk + N_Gas);

    Bulge bulge(this);
    bulge.generateBulge(N_Disk + N_Gas, N_Disk + N_Bulge + N_Gas);

    Halo halo(this);
    halo.generateDarkMatterHalo(N_Disk + N_Bulge + N_Gas, N_Disk + N_Bulge + N_Halo + N_Gas);


    std::cout << "Galaxy generated with " << particles->size() << " particles." << std::endl;
}