#include "Log.h"
#include <cstdio> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <filesystem>
#include <chrono>


std::string formatWithComma(double value) 
{
    std::ostringstream stream;
    stream.imbue(std::locale("C"));
    stream << std::fixed << std::setprecision(5) << value;
    std::string result = stream.str();

    std::replace(result.begin(), result.end(), '.', ',');
    return result;
}

namespace Log
{
    std::string outputDir;
    void setOutputDir(const std::string& dir) 
    {
        outputDir = dir;
        if (!std::filesystem::exists(outputDir)) 
        {
            std::filesystem::create_directories(outputDir);
        }
    }

//save data in csv file
    std::vector<std::ofstream> dataFiles;

    void printData(const std::string& filename, const double x, const double y) 
    {
        std::string fullPath = outputDir + "/" + filename;
        auto it = std::find_if(dataFiles.begin(), dataFiles.end(),
            [&fullPath](const std::ofstream& f) { return f.is_open() && f.rdbuf()->is_open(); });

        if (it == dataFiles.end()) 
        {
            dataFiles.emplace_back(std::ofstream(fullPath, std::ios::out | std::ios::app));
            it = dataFiles.end() - 1;
        }

        *it << formatWithComma(x) << ";" << formatWithComma(y) << "\n";
        it->flush();
    }

//track process time
    bool hasStarted = false;
    std::ofstream LogsDir;
    std::ofstream proccessFile;
    std::chrono::steady_clock::time_point startTimestamp;
    std::string currentProcessName;


    void startProcess(const std::string& processName) 
    {
        if(hasStarted == true)
        {
            endProcess();
        }
        else
        {
            std::string filename =  outputDir + "/processLog.csv";

            if (std::ifstream(filename)) 
            {
                std::cout << "Logfile exists. Deleting existing logfile: " << filename << std::endl;
                if (std::remove(filename.c_str()) != 0) 
                {
                    std::cerr << "Error deleting logfile." << std::endl;
                }
            }

            proccessFile.open(filename, std::ios::out | std::ios::app);

            if (proccessFile.is_open()) 
            {
                std::cout << "Logfile opened successfully." << std::endl;
            } 
            else 
            {
                std::cerr << "Error opening logfile." << std::endl;
        }
            hasStarted = true;
        }
        currentProcessName = processName;
        startTimestamp = std::chrono::steady_clock::now();
    }

    void endProcess() 
    {
        if (!currentProcessName.empty()) {
            auto now = std::chrono::steady_clock::now();
            double elapsedTime = std::chrono::duration<double>(now - startTimestamp).count();
            proccessFile << currentProcessName << ";" << formatWithComma(elapsedTime) << "\n";
            
            proccessFile.flush();

            currentProcessName.clear();
        }
    }
}