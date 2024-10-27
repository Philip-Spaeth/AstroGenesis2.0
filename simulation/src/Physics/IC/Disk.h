#pragma once
#include "Galaxy.h"

class Disk
{
public:
    Disk(Galaxy * g) : g(g) {}
    ~Disk() {}

    Galaxy * g;

    void generateStarDisk(int start, int end);

    void generateGasDisk(int start, int end);
};