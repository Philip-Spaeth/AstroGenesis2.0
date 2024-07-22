#include <iostream>
#include "Simulation.h"

using namespace Constants;

int main()
{
    Simulation* simulation = new Simulation();
    simulation->run();
    delete simulation;
    return 0;
}