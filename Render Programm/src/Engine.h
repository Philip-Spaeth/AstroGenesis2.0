#ifndef ENGINE_H
#define ENGINE_H

#include <iostream>
#include <vector>
#include <memory>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "vec3.h"
#include "mat4.h"
#include "vec4.h"
#include "Particle.h"
#include <cmath>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <thread>

#ifdef WIN32
#define GLFW_EXPOSE_NATIVE_WIN32
#include <GLFW/glfw3native.h>
#endif

class Engine {
public:
    Engine(std::string dataFolder, double deltaTime, double numOfParticles, double numTimeSteps, std::vector<std::shared_ptr<Particle>>* particles);
    ~Engine();

    double deltaTime;
    double numOfParticles;
    double numTimeSteps;
    
    std::vector<std::shared_ptr<Particle>>* particles;

    std::string dataFolder;

    bool init(double physicsFaktor);
    void start();
    void update(int index);
    bool clean();

    GLuint pbo = 0;
    std::queue<std::pair<std::string, std::vector<unsigned char>>> saveQueue;
    std::mutex queueMutex;
    std::condition_variable queueCondition;
    std::thread saveThread;
    bool terminateThread = false;
    int width, height;

    void initializePBO();
    void saveWorker();

    void initializePBO(int width, int height);
    void saveAsPicture(const std::string& folderName, int index);

    static void window_iconify_callback(GLFWwindow* window, int iconified);
    bool RenderLive = true;
    std::string videoName;

    GLFWwindow* window;

    bool isRunning = false;

    int maxNumberOfParticles = 1e7;

    double playSpeed = 1;
    double changeSpeed = 1;

    bool showBaryonicMatter = true;
    bool showDarkMatter = true;
    int colorMode = 2;

    double passedTime = 0;

    double globalScale = 1e-9;
    void calculateGlobalScale();

    bool focusedCamera = false; 
    vec3 cameraPosition;
    vec3 cameraFront;
    vec3 cameraUp;
    double cameraYaw;
    double cameraPitch;
    double cameraSpeed = 100;
    double rushSpeed = 1000;

    //render Tree
    const double theta = 0;
    void framebuffer_size_callback(GLFWwindow* window, int width, int height);
    
    double dFromCenter = 0;
private:
    double calcDensityAv();
    double densityAv = 0;
    int oldIndex = -1;
    bool BGstars = true;
    int amountOfStars = 1000;
    std::vector<vec4> bgStars;

    bool tracks = false;
    double cameraViewDistance = 1e15;
    mat4 view;

    // Funktionen für Kamerasteuerung 
    void processMouseInput();
    void processInput();

    GLuint shaderProgram;
    bool shouldClose = false;
    GLuint VAO;
    GLuint instanceVBO;
    void renderParticles();
    void checkShaderCompileStatus(GLuint shader, const char* shaderType);
    void checkShaderLinkStatus(GLuint program);
    void calcTime(int index);
    void calcTime(vec3 position, int index);
    double faktor = -1;
    double random(double min, double max);
};

#endif // ENGINE_H
