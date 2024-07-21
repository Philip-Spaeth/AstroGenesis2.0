#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <string.h>
#include <filesystem>
#include "DataManager.h"
#include "Engine.h"
#include <memory>
#include <thread>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define TARGET_FPS 60

//  Nur unter Windows
#ifdef WIN32
#include <Windows.h>
#include <debugapi.h>
#endif

using namespace std;
namespace fs = std::filesystem;

int main()
{
    std::string dataFolder; // Pfad zum Data-Ordner eine Ebene höher

    //print out all the folders in in the folder "../../Data/" as options to choose from wich data to load
    std::cout << "Choose a data folder to load: " << std::endl;
    
    cout << endl;

    int i = 1;
    for (const auto& entry : fs::directory_iterator("../../Data/"))
    {
        std::cout << "["<< i << "]   " << entry.path().filename() << std::endl;
        i++;
    }

    //get the input from the user
    int folderIndex;
    std::cin >> folderIndex;

    //get the folder name from the input
    i = 1;
    for (const auto& entry : fs::directory_iterator("../../Data/"))
    {
        if (i == folderIndex)
        {
            dataFolder = entry.path().string() + "/";
            break;
        }
        i++;
    }

    cout << endl;

    //get the DataFolder
    int deltaTime = 1;
    int numOfParticles = 1;
    int numTimeSteps = 1;

    std::string input;
    std::getline(std::cin, input);

    DataManager dataManager;
    dataManager.path = dataFolder;
    dataManager.readInfoFile(deltaTime, numTimeSteps, numOfParticles);
    std::cout << "deltaTime: " << deltaTime << std::endl;
    std::cout << "numOfParticles: " << numOfParticles << std::endl;
    std::cout << "numTimeSteps: " << numTimeSteps << std::endl;

    std::vector<std::shared_ptr<Particle>> particles;

    Engine engine(dataFolder, deltaTime, numOfParticles, numTimeSteps, &particles);

    if (!engine.init(1.0)) { // Hier muss der physikalische Faktor übergeben werden
        std::cerr << "Engine initialization failed." << std::endl;
        return -1;
    }

    engine.start();

    double lastFrameTime = glfwGetTime(); // Zeit des letzten Frames
    double frameTime; // Zeitdauer eines Frames

    int frameCount = 0;
    double secondCounter = 0.0;
    int counter = 0;

    std::vector<Particle> currentParticles;

    // Haupt-Render-Schleife
    while (!glfwWindowShouldClose(engine.window))
    {
        #ifdef WIN32
        // check for exit Programm with Key ESC
        if (GetAsyncKeyState(27) & 0x8000)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            //OutputDebugString(L"ESC KEY\n");
            glfwSetWindowShouldClose(engine.window, true);
        }
        #endif

        double currentFrameTime = glfwGetTime();
        frameTime = currentFrameTime - lastFrameTime;
        lastFrameTime = currentFrameTime;

        // Rufen Sie die update-Funktion auf und übergeben Sie die Zeitdauer eines Frames
        if (frameTime < 1.0 / TARGET_FPS)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(static_cast<int>((1.0 / TARGET_FPS - frameTime) * 1000)));
        }

        // load particles
        dataManager.loadData(counter, particles);

        // update particles
        engine.update(counter);
        // add time when engine is running
        if (engine.isRunning)
        {
            if (counter >= 0 && counter <= numTimeSteps - 1)
            {
                counter = counter + engine.playSpeed;
            }
        }
        if (counter >= numTimeSteps - 1)
        {
            counter = numTimeSteps - 1;
            engine.playSpeed = 0;
        }
        if (counter < 0)
        {
            counter = 0;
            engine.playSpeed = 0;
        }

        #ifdef WIN32
        //restart if R is pressed
        if (GetAsyncKeyState(82) & 0x8000)
        {
            counter = 0;
            engine.playSpeed = 0;
        }
        #endif

        frameCount++;
        secondCounter += frameTime;

        // Wenn eine Sekunde vergangen ist, zeigen Sie die FPS an
        if (secondCounter >= 1.0)
        {
            char numStr[20]; // Ein char-Array zur Speicherung der Zeichenkette
            // Verwende sprintf, um die Ganzzahl in das char-Array umzuwandeln
            snprintf(numStr, sizeof(numStr), "%d", frameCount);

            //std::cout << "Die umgewandelte Zeichenkette: " << numStr << std::endl;
            #ifdef WIN32
            strcat_s(numStr, " FPS");
            #else
            strcat(numStr, " FPS");
            #endif

            glfwSetWindowTitle(engine.window, numStr);
            //std::cout << "FPS: " << frameCount << std::endl;
            frameCount = 0;
            secondCounter = 0.0;
        }

        // Wenn F11 gedrückt wird, schalten Sie in den Vollbildmodus um
        if (glfwGetKey(engine.window, GLFW_KEY_F11) == GLFW_PRESS) 
        {
            // Speichern Sie die aktuelle Monitor-Information
            GLFWmonitor* monitor = glfwGetPrimaryMonitor();
            const GLFWvidmode* mode = glfwGetVideoMode(monitor);

            // Prüfen Sie den aktuellen Vollbild-Status des Fensters
            if (glfwGetWindowMonitor(engine.window) == NULL) {
                // Wenn das Fenster nicht im Vollbildmodus ist, schalten Sie in den Vollbildmodus
                glfwSetWindowMonitor(engine.window, monitor, 0, 0, mode->width, mode->height, mode->refreshRate);
                engine.framebuffer_size_callback(engine.window, mode->width, mode->height);
            }
            else {
                // Wenn das Fenster bereits im Vollbildmodus ist, schalten Sie in den Fenstermodus
                glfwSetWindowMonitor(engine.window, NULL, 100, 100, 1200, 800, mode->refreshRate);
                engine.framebuffer_size_callback(engine.window, 1200, 800);
            }
        }
    }

    // Beenden Sie GLFW
    engine.clean();
    glfwTerminate();
    return 0;
}
