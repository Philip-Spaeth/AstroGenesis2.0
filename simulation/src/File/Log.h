#ifndef LOG_H
#define LOG_H

#include <iostream>
#include <fstream>
#include <chrono>
#include <string>

// Funktionen zur Initialisierung und Loggen
void initLogger(const std::string& filename);
void closeLogger();
void start(const std::string& processName);
void end();

extern bool hasStarted;

// Die Log-Datei und Variablen werden als extern deklariert
extern std::ofstream logFile;
extern std::chrono::steady_clock::time_point startTimestamp;
extern std::string currentProcessName;

#endif // LOG_H
