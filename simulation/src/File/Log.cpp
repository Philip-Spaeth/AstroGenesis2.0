#include "Log.h"
#include <cstdio> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>

std::string formatWithComma(double value) {
    std::ostringstream stream;
    stream.imbue(std::locale("C")); // Standardlocale für Dezimalpunkt
    stream << std::fixed << std::setprecision(5) << value;
    std::string result = stream.str();

    // Ersetze Punkt durch Komma
    std::replace(result.begin(), result.end(), '.', ',');
    return result;
}

namespace Log
{

    // Globale Variablen
    std::ofstream logFile;
    std::chrono::steady_clock::time_point startTimestamp;
    std::string currentProcessName;
    bool hasStarted;

    void initLogger(const std::string& filename) 
    {
        if (std::ifstream(filename)) 
        {
            std::cout << "Logfile exists. Deleting existing logfile: " << filename << std::endl;
            if (std::remove(filename.c_str()) != 0) {
                std::cerr << "Error deleting logfile." << std::endl;
            }
        }

        logFile.open(filename, std::ios::out | std::ios::app);
        if (logFile.is_open()) {
            std::cout << "Logfile opened successfully." << std::endl;
            if (logFile.tellp() == 0) {
                //logFile << "Process, Time\n";
            }
        } else {
            std::cerr << "Error opening logfile." << std::endl;
        }
        hasStarted = false;
    }

    // Funktion zum Schließen der Log-Datei
    void closeLogger() {
        if (logFile.is_open()) {
            logFile.close();
        }
    }

    void printData(const double x, const double y) 
    {
        logFile << formatWithComma(x) << ";" << formatWithComma(y)<< "\n";
        logFile.flush();
    }

    void start(const std::string& processName) 
    {
        if(hasStarted == true)
        {
            end();
        }
        else
        {
            hasStarted = true;
        }
        currentProcessName = processName;
        startTimestamp = std::chrono::steady_clock::now();
    }

    void end() {
        if (!currentProcessName.empty()) {
            auto now = std::chrono::steady_clock::now();
            double elapsedTime = std::chrono::duration<double>(now - startTimestamp).count();
            // Schreibe in die Log-Datei
            logFile << currentProcessName << ";" << elapsedTime << "\n";
            
            // Stelle sicher, dass der Puffer geschrieben wird
            logFile.flush();

            // Rücksetzen für den nächsten Prozess
            currentProcessName.clear();
        }
    }
}