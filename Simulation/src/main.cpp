#include <iostream>
#include "Simulation.h"

int main()
{
    Simulation* simulation = new Simulation();
    simulation->run();
    delete simulation;
    return 0;
}