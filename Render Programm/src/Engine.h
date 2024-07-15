#ifndef ENGINE_H
#define ENGINE_H

#include <vector>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "Particle.h"
#include "Physics.h"
#include "FileManager.h"

class Physics;
class FileManager;

class Engine {
public:
    Engine(std::string dataFolder);
    std::string dataFolder;

    bool init(double physicsFaktor);
    void start(Physics* p);
    void update(int index, std::string folderName);
    bool clean();
    void saveAsPicture(std::string folderName, int index);

    bool RenderLive = true;

    GLFWwindow* window;
    Physics* physics;

    bool isRunning = false;

    int maxNumberOfParticles = 1e6;

    double playSpeed = 1;
    double changeSpeed = 1;

    bool showBaryonicMatter = true;
    bool showDarkMatter = true;
    int colorMode = 2;

    double passedTime = 0;

    double globalScale = 1e-9;
    //double globalScale = 1e-18;
    void calculateGlobalScale();

    bool focusedCamera = false; 
    glm::dvec3 cameraPosition;
    glm::dvec3 cameraFront;
    glm::dvec3 cameraUp;
    double cameraYaw;
    double cameraPitch;
    double cameraSpeed = 100;
    double rushSpeed = 1000;

    //render Tree
    const double theta = 0;
    std::vector<glm::vec4> positions;
    std::vector<glm::vec3> colors;
    std::vector<glm::vec3> densityColors;
    std::vector<glm::vec3> thermalColors;
    std::vector<glm::vec1> isDarkMatter;
    void framebuffer_size_callback(GLFWwindow* window, int width, int height);
private:
    int oldIndex = 0;
    bool BGstars = true;
    int amountOfStars = 1000;
    std::vector<glm::vec4> bgStars;

    bool tracks = false;
    double cameraViewDistance = 1e15;
    glm::mat4 view;

    // Funktionen f�r Kamerasteuerung 
    void processMouseInput();
    void processInput();

    GLuint shaderProgram;
    bool shouldClose = false;
    GLuint VAO;
    GLuint instanceVBO;
    GLuint framebuffer;
    GLuint textureColorbuffer;
    GLuint blurSizeLocation;
    GLuint rbo;
    GLuint blurShaderProgram; // Shader-Programm f�r den Blur-Effekt
    GLuint quadVAO; // VAO f�r das Quad, auf das der Blur-Effekt angewendet wird
    void renderBlur();
    void renderParticles();
    void checkShaderCompileStatus(GLuint shader, const char* shaderType);
    void checkShaderLinkStatus(GLuint program);
    void calcTime(int index);
    void calcTime(glm::dvec3 position, int index);
    double faktor = -1;
    FileManager* fileManager;
    double random(double min, double max);
};

#endif