#include <iostream>
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
#include <string>
#include <windows.h>
#include <commdlg.h>
#include <shlobj.h>
#include <shobjidl.h>
#include <comdef.h>
#include <windows.h>
#include <commdlg.h>
#include <shlobj.h>
#include <shobjidl.h>
#include <comdef.h>
#include <combaseapi.h>
#include <shlguid.h>

#define TARGET_FPS 60

//  Nur unter Windows
#ifdef WIN32
#include <Windows.h>
#include <debugapi.h>
#endif

using namespace std;
namespace fs = std::filesystem;


//get the DataFolder
double deltaTime = 1;
double numOfParticles = 1;
double numTimeSteps = 1;


std::string dataFolder; // Pfad zum Data-Ordner eine Ebene höher
DataManager dataManager("");

void renderLive();
void renderVideo();

int main()
{
    //Tittle of the Programm with information
    std::cout << std::endl << "<--------------------------------------------  Astro Genesis Render Programm ------------------------------------------>"  << std::endl<< std::endl;

    //contrebutors and ownersion information
    std::cout << "this Render Engine is part of the Astro Genesis Project" << std::endl<< std::endl;

    //contrebutors and ownersion information
    std::cout << "developed by: Philip Spaeth and Kimi Sickinger" << std::endl<< std::endl;
    //licence information and copy right information, licence: GNU General Public License v3.0
    std::cout << "This software is licensed under the GNU General Public License v3.0" << std::endl<< std::endl;
    

    //print out all the folders in in the folder "../../Data/" as options to choose from wich data to load
    std::cout << "Choose a data folder to load:  " << std::endl;
    
    cout << endl;

    //check if there is a folder in the directory
    if (!fs::exists("../../output_data/"))
    {
        //create the Data folder
        fs::create_directory("../../output_data/");
    }

    int i = 1;
    for (const auto& entry : fs::directory_iterator("../../output_data/"))
    {
        std::cout << "["<< i << "]   " << entry.path().filename() << std::endl;
        i++;
    }

    //get the input from the user
    int folderIndex;
    std::cin >> folderIndex;

    //get the folder name from the input
    i = 1;
    for (const auto& entry : fs::directory_iterator("../../output_data/"))
    {
        if (i == folderIndex)
        {
            dataFolder = entry.path().string() + "/";
            break;
        }
        i++;
    }

    cout << endl;

    dataManager.path = dataFolder;

    std::string filename = dataManager.path + "info.txt";
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Could not open the file!" << std::endl;
        return 0;
    }

    dataManager.readInfoFile(deltaTime, numTimeSteps, numOfParticles);
    //std::cout << "deltaTime: " << deltaTime << std::endl;
    //std::cout << "numOfParticles: " << numOfParticles << std::endl;
    //std::cout << "numTimeSteps: " << numTimeSteps << std::endl;

    //choose between Rendering a video or liveViewer
    std::cout << "Choose between rendering a video or liveViewer: " << std::endl;
    std::cout << "[1]   Video" << std::endl;
    std::cout << "[2]   LiveViewer" << std::endl;

    int choice;
    std::cin >> choice;

    //if enter is pressed, the liveViewer will be started, else check the input
    if (choice == 2)
    {
        renderLive();
    }
    else if (choice == 1)
    {
        renderVideo();
    }
    else
    {
        std::cout << "Invalid input" << std::endl;
    }

    return 0;
}

void renderVideo()
{
    std::getline(std::cin, dataFolder);
    std::cout << "\nPlease enter a name for the video: ";
    std::string videoName;
    std::getline(std::cin, videoName);
        
    std::cout << "\nDo you want to rotate the camera around the center? (y/n): ";
    std::string input2;
    std::getline(std::cin, input2);
    
    bool rotate = false;
    double speed = 0;
    if (input2 == "y")
    {
        rotate = true;
        std::cout << "\nPlease enter the speed of the rotation(between 10-1000): ";
        std::getline(std::cin, input2);
        speed = std::stod(input2);
    }

    //ask about the distance from the center factor
    std::cout << "\nPlease enter the distance from the center factor (between 0.1-10): ";
    std::string input;
    std::getline(std::cin, input);
    double distance = std::stod(input);
    
    std::cout << "" << std::endl;
    
    std::vector<std::shared_ptr<Particle>> particles;
    Engine engine(dataFolder, deltaTime, numOfParticles, numTimeSteps, &particles);
    engine.RenderLive = false;

    if (!engine.init(1.0)) {
        std::cerr << "Engine initialization failed." << std::endl;
        return;
    }

    engine.start();
    engine.videoName = videoName;

    double lastFrameTime = glfwGetTime(); 
    double frameTime;
    int frameCount = 0;
    double secondCounter = 0.0;
    int counter = 0;
    std::vector<Particle> currentParticles;

    engine.cameraSpeed = speed;
    engine.focusedCamera = true;
    
    engine.cameraPosition = vec3(0, 100, 1000* distance);
    engine.cameraFront = vec3(0, 0, -1);
    engine.cameraUp = vec3(0, 1, 0);
    engine.cameraYaw = -90;
    engine.cameraPitch = 0;
    
    engine.dFromCenter = distance;


    while (!glfwWindowShouldClose(engine.window))
    {

        double currentFrameTime = glfwGetTime();
        frameTime = currentFrameTime - lastFrameTime;
        lastFrameTime = currentFrameTime;

        dataManager.loadData(counter, particles);
        engine.isRunning = true;

        engine.update(counter);
        
        if (counter >= numTimeSteps - 1)
        {
            engine.clean();
            glfwTerminate();
            std::cout << "" << std::endl;
            std::cout << "" << std::endl;
            std::cout << "All steps rendered" << std::endl;
            std::string workingDir = fs::current_path().string();
            std::cout << "You can find the video in the folder: " << videoName << std::endl;
            return;
        }

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

        frameCount++;
        secondCounter += frameTime;

        if (secondCounter >= 1.0)
        {
            char numStr[20];
            snprintf(numStr, sizeof(numStr), "%d", frameCount);
            #ifdef WIN32
            strcat_s(numStr, " FPS");
            #else
            strcat(numStr, " FPS");
            #endif
            glfwSetWindowTitle(engine.window, numStr);
            frameCount = 0;
            secondCounter = 0.0;
        }

        dataManager.printProgress((double)counter, (double)numTimeSteps, "");
    }

    engine.clean();
    glfwTerminate();
    return;
}

void renderLive()
{
    std::vector<std::shared_ptr<Particle>> particles;

    Engine engine(dataFolder, deltaTime, numOfParticles, numTimeSteps, &particles);

    if (!engine.init(1.0)) { // Hier muss der physikalische Faktor übergeben werden
        std::cerr << "Engine initialization failed." << std::endl;
        return;
    }

    engine.deltaTime = deltaTime;
    engine.numOfParticles = numOfParticles;
    engine.numTimeSteps = numTimeSteps;

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
    return;
}
