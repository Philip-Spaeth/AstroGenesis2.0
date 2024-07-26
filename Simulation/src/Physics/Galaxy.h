#pragma once

#include <iostream>
#include "Particle.h"
#include <memory>
#include <vector>
#include "Constants.h"


class Galaxy
{
public:
    Galaxy(std::vector<std::shared_ptr<Particle>>& particles, int start, int end, double mass, double radius);
    ~Galaxy(){};
};