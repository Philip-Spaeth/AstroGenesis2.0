#include "Console.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
#include <windows.h>
#include <intrin.h>
#else
#include <unistd.h>
#endif

using namespace std;
namespace fs = std::filesystem;

#ifdef _WIN32
void setConsoleColor(WORD color) {
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hConsole, color);
}
#else
void setConsoleColor(int color) {
    // No color setting for Linux in this version
}
#endif
void Console::printProgress(double currentStep, double steps, std::string text) 
{
    if(currentStep == steps) return;
    
    static const int barWidth = 70;

    if (!timerStarted) {
        std::cout << std::endl;
        startTime = std::chrono::high_resolution_clock::now();
        timerStarted = true;
    }

    double progress = (currentStep + 1) / steps;
    int pos = static_cast<int>(barWidth * progress);
    
    // Set color based on progress
#ifdef _WIN32
    WORD progressColor = (currentStep < steps - 1) ? FOREGROUND_RED : FOREGROUND_GREEN;
#endif

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
#ifdef _WIN32
        if (i < pos) {
            setConsoleColor(progressColor);
            std::cout << "=";
        }
        else if (i == pos && currentStep < steps - 1) {
            setConsoleColor(FOREGROUND_RED);
            std::cout << ">";
        }
        else {
            setConsoleColor(FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED);
            std::cout << " ";
        }
#else
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
#endif
    }

    double remainingTime = 0;

    auto currentTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = currentTime - startTime;
    double estimatedTotalTime = (elapsed.count() / (currentStep + 1)) * steps;
    double estimatedRemainingTime = estimatedTotalTime - elapsed.count();
    remainingTime = estimatedRemainingTime;

    std::string unit;
    if (remainingTime < 60) {
        unit = "s";
    } else if (remainingTime < 3600) {
        remainingTime /= 60;
        unit = "min";
    } else {
        remainingTime /= 3600;
        unit = "h";
    }

    std::ostringstream timeStream;
    timeStream << std::fixed << std::setprecision(1) << remainingTime;
    std::string timeleft = timeStream.str();

#ifdef _WIN32
    setConsoleColor(FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED);
#endif
    std::cout << "] " << int(progress * 100.0) << " %  Estimated remaining time: " << timeleft << unit << "  "<< text << "       " << "\r";

    if (currentStep == steps - 1) {
        std::cout << std::endl;
        std::cout << std::endl;
    }
}

void Console::printSystemInfo()
{
    std::cout << "\nComputational parameters:" << std::endl;
    // Number of CPU cores
    int numCores = 0;
    #ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    numCores = sysinfo.dwNumberOfProcessors;
    #else
    numCores = sysconf(_SC_NPROCESSORS_ONLN);
    #endif

    std::cout << "  Number of CPU Cores: " << numCores << std::endl;

    // CPU Clock Speed
    double cpuFrequency = 0.0;
    #ifdef _WIN32
    HKEY hKey;
    LONG lError = RegOpenKeyEx(HKEY_LOCAL_MACHINE,
                               "HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0",
                               0, KEY_READ, &hKey);
    if (lError == ERROR_SUCCESS) {
        DWORD dwMHz;
        DWORD bufferSize = sizeof(dwMHz);
        lError = RegQueryValueEx(hKey, "~MHz", NULL, NULL, (LPBYTE)&dwMHz, &bufferSize);
        if (lError == ERROR_SUCCESS) {
            cpuFrequency = dwMHz;
        }
        RegCloseKey(hKey);
    }
    #else
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line;
    while (std::getline(cpuinfo, line)) {
        if (line.find("cpu MHz") != std::string::npos) {
            cpuFrequency = std::stod(line.substr(line.find(":") + 1));
            break;
        }
    }
    #endif

    //display the frequency in GHz iwth 2 decimal places
    std::cout << "  CPU Clock Speed: " << std::fixed << std::setprecision(2) << cpuFrequency / 1000 << " GHz" << std::endl;

    // Available RAM
    size_t totalRam = 0;
    #ifdef _WIN32
    MEMORYSTATUSEX statex;
    statex.dwLength = sizeof(statex);
    GlobalMemoryStatusEx(&statex);
    totalRam = statex.ullTotalPhys / (1024 * 1024); // in MB
    #else
    std::ifstream meminfo("/proc/meminfo");
    while (std::getline(meminfo, line)) {
        if (line.find("MemTotal") != std::string::npos) {
            totalRam = std::stoul(line.substr(line.find(":") + 1)) / 1024; // in MB
            break;
        }
    }
    #endif

    std::cout << "  Available RAM: " << std::fixed << std::setprecision(1) << (double)totalRam / 1000.0 << " GB" << std::endl;

    // Storage Write Speed
    const size_t fileSize = 100 * 1024 * 1024; // 100 MB
    char *buffer = new char[fileSize];
    std::fill(buffer, buffer + fileSize, 'A');

    auto start = std::chrono::high_resolution_clock::now();

    std::ofstream file("test_file.bin", std::ios::binary);
    file.write(buffer, fileSize);
    file.close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    double writeSpeed = fileSize / duration.count() / (1024 * 1024); // MB/s

    std::cout << "  Storage Write Speed: " << writeSpeed << " MB/s \n" << std::endl;

    delete[] buffer;
}