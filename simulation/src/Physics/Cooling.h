#pragma once
#include <cmath>
#include <vector>
#include <memory>
#include "Particle.h"

class Cooling
{
public:
    Cooling() {}
    ~Cooling() {}

    void coolingRoutine(std::shared_ptr<Particle>& particle);
private:

};