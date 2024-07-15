#pragma once

#include <string>
#include <vector>
#include <map>
#include <memory>

using namespace std;

class DataManager
{
public:
    DataManager();
    ~DataManager();

    void saveData(Particle *particle);
}