#include "random.h"
#include <random>
#include <chrono>
#include <cmath>
#include <cstdlib>

double random::gaussianRandom(double mean, double stddev)
{
    // Erstellen eines Zufallszahlengenerators
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    // Erstellen einer Normalverteilung
    std::normal_distribution<double> distribution(mean, stddev);

    // Generieren und R�ckgabe einer Zahl aus dieser Verteilung
    return distribution(generator);

}

void random::setRandomSeed(unsigned int seed)
{
    // Setze den Zufallszahlengenerator auf einen bestimmten Wert
    srand(seed);
}

double random::between(double min, double max)
{
    // Generiere eine zuf�llige Gleitkommazahl zwischen min und max
    double randomFloat = min + static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * (max - min);

    return randomFloat;
}