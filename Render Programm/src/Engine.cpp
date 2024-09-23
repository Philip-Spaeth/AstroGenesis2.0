#include "Engine.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <filesystem>
#include <algorithm>

#ifdef WIN32
#include <Windows.h>
#endif

#include <chrono>
#include <thread>
#include <string>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

inline float radians(float degrees) {
    return degrees * (M_PI / 180.0f);
}

Engine::Engine(std::string NewDataFolder, double deltaTime, double numOfParticles, double numTimeSteps, std::vector<std::shared_ptr<Particle>>* particles) : window(nullptr), shaderProgram(0), VAO(0)
{
    dataFolder = NewDataFolder;
    this->deltaTime = deltaTime;
    this->numOfParticles = numOfParticles;
    this->numTimeSteps = numTimeSteps;
    this->particles = particles;

    // start kamera position
    cameraPosition = vec3(0, 0, 0);                  
    cameraFront = vec3(0.0, 0.0, -1.0);
    cameraUp = vec3(0.0, 1.0, 0.0);
    cameraYaw = -90.0f;
    cameraPitch = 0.0f;

    saveThread = std::thread(&Engine::saveWorker, this);
}

Engine::~Engine() {
    // Terminate the save worker thread
    {
        std::unique_lock<std::mutex> lock(queueMutex);
        terminateThread = true;
        queueCondition.notify_all();
    }
    saveThread.join();
    if (pbo != 0) {
        glDeleteBuffers(1, &pbo);
    }
}

bool Engine::init(double physicsFaktor) 
{
    faktor = physicsFaktor;
    // GLFW initialisieren
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return false;
    }

    // GLFW-Konfiguration setzen
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

    // Ein GLFW-Fenster erstellen
    double width = 1200;
    double height = 800;
    if(RenderLive == false)
    {
        width = 1920;
        height = 1080;
    }
    window = glfwCreateWindow(width, height, "Particle Rendering", nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return false;
    }

    // Fensterstile ändern, um Minimierung und Maximierung zu verhindern
    #ifdef WIN32
    HWND hwnd = glfwGetWin32Window(window);
    LONG style = GetWindowLong(hwnd, GWL_STYLE);
    style &= ~WS_MINIMIZEBOX; // Deaktiviere Minimierungsbox
    style &= ~WS_MAXIMIZEBOX; // Deaktiviere Maximierungsbox
    SetWindowLong(hwnd, GWL_STYLE, style);
    #endif

    // GLFW-Kontext setzen
    glfwMakeContextCurrent(window);

    // GLEW initialisieren
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        return false;
    }

    // Shader-Programm kompilieren und verlinken
    const char* vertexShaderSource = "#version 330 core\n"
        "layout (location = 0) in vec3 position;\n"
        "uniform mat4 projection;\n"
        "uniform mat4 view;\n"
        "uniform vec3 particlePosition; // Neue Uniform-Variable für die Partikelposition\n"
        "void main()\n"
        "{\n"
        "    // Berechnen Sie die endgültige Position des Partikels, indem Sie die Partikelposition hinzufügen\n"
        "    vec4 finalPosition = projection * view * vec4(position + particlePosition, 1.0);\n"
        "    gl_Position = finalPosition;\n"
        "}\0";

    const char* fragmentShaderSource = "#version 330 core\n"
        "out vec4 FragColor;\n"
        "uniform vec3 particleColor;\n"
        "void main()\n"
        "{\n"
        "    FragColor = vec4(particleColor, 1.0); // Weiß\n"
        "}\n\0";

    // Erstellen des Shader-Programms und kompilieren
    GLuint vertexShader, fragmentShader;
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);
    checkShaderCompileStatus(vertexShader, "VERTEX");

    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);
    checkShaderCompileStatus(fragmentShader, "FRAGMENT");

    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    checkShaderLinkStatus(shaderProgram);
    if (RenderLive)
    {
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    }
    else
    {
		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
	}

    return true;
}

void Engine::framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // Stellen Sie sicher, dass die Ansichtsportgröße dem neuen Fenster entspricht
    glViewport(0, 0, width, height);
}

void Engine::start()
{
    int size = numOfParticles;
    if (size > maxNumberOfParticles)
    {
        size = maxNumberOfParticles;
	}
    
    // Hier VBO und VAO erstellen und konfigurieren
    GLuint VBO;
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    // Vektor zur Speicherung der Partikelpositionen vorbereiten
    std::vector<vec4> positions;
    positions.reserve(particles->size());
    for (const auto& particle : *particles)
    {
        positions.emplace_back(particle->position.x, particle->position.y, particle->position.z, 1.0f);
    }

    // Stelle sicher, dass du die korrekte Größe und den korrekten Typ für die Daten verwendest
    glBufferData(GL_ARRAY_BUFFER, sizeof(vec4) * positions.size(), positions.data(), GL_STATIC_DRAW);

    // Erstellen des Vertex Array Objects (VAO)
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    // Konfigurieren des VAO für das VBO
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    //place the BGstars in the background
    if (BGstars)
    {
        for (int i = 0; i < amountOfStars; i++)
        {
            double x = random(-1e14, 1e14);
			double y = random(-1e14, 1e14);
			double z = random(-1e14, 1e14);
            double size = random(0.1, 2);
			bgStars.push_back(vec4(x, y, z, size));
		}   
    }

    std::cout << "Data loaded" << std::endl;
}

void Engine::initializePBO() {
    if (pbo == 0) {
        glGenBuffers(1, &pbo);
        glBindBuffer(GL_PIXEL_PACK_BUFFER, pbo);
        glBufferData(GL_PIXEL_PACK_BUFFER, 3 * width * height, nullptr, GL_STREAM_READ);
        glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    }
}

void Engine::saveAsPicture(const std::string& folderName, int index) {
    // Speichern Sie das gerenderte Bild als BMP-Datei
    glfwGetFramebufferSize(window, &width, &height);

    // Initialize PBO if not already initialized
    initializePBO();

    glBindBuffer(GL_PIXEL_PACK_BUFFER, pbo);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, 0);

    // Map the PBO to process its data on the CPU
    unsigned char* data = (unsigned char*)glMapBuffer(GL_PIXEL_PACK_BUFFER, GL_READ_ONLY);

    // Copy data to a vector to pass to the saving thread
    std::vector<unsigned char> imageData(data, data + 3 * width * height);

    glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

    // Umkehren der Pixelreihenfolge in-place
    for (int y = 0; y < height / 2; y++) {
        for (int x = 0; x < 3 * width; x++) {
            std::swap(imageData[3 * width * y + x], imageData[3 * width * (height - 1 - y) + x]);
        }
    }

    std::string fullPath = "../Video_Output/" + videoName + "/";
    std::filesystem::create_directory(fullPath);

    std::string filename = fullPath + "Picture_" + std::to_string(index) + ".bmp";

    // Add the data to the save queue
    {
        std::unique_lock<std::mutex> lock(queueMutex);
        saveQueue.emplace(filename, std::move(imageData));
        queueCondition.notify_one();
    }
}

void Engine::saveWorker() {
    while (true) {
        std::pair<std::string, std::vector<unsigned char>> item;
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            queueCondition.wait(lock, [this] { return !saveQueue.empty() || terminateThread; });

            if (terminateThread && saveQueue.empty()) {
                break;
            }

            item = std::move(saveQueue.front());
            saveQueue.pop();
        }

        stbi_write_bmp(item.first.c_str(), width, height, 3, item.second.data());
    }
}

void Engine::window_iconify_callback(GLFWwindow* window, int iconified) {
    if (iconified) {
        glfwRestoreWindow(window);
    }
}

void Engine::update(int index)
{
    if(index == 0)
    {
        densityAv = calcDensityAv();
    }

    //calculate the time
    if (isRunning && index != oldIndex && RenderLive)
    {
        calcTime(index);
    }
    if (RenderLive)
    {
        //continue if space is pressed
        #ifdef WIN32
        if (GetAsyncKeyState(32) & 0x8000)
        #else
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
        #endif
        {
			isRunning = !isRunning;
            //small delay
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
		}
	}

    if (RenderLive)
    {
        processMouseInput();
    }

    if(isRunning && RenderLive == false) processInput();
    if (RenderLive) processInput();

    // set the globalScale of the system
    if (index == 0)
    {
        calculateGlobalScale();
    }

    renderParticles();

    glfwSwapBuffers(window);
    glfwPollEvents();

    //if video is rendered
    if (RenderLive == false && index != oldIndex && index != 0)
    {
        saveAsPicture(dataFolder, index);
        glfwSetWindowIconifyCallback(window, window_iconify_callback);
    }

    if(RenderLive == true)
    {
        //speed up if right arrow is pressed
    #ifdef WIN32
        if (GetAsyncKeyState(39) & 0x8000)
    #else
        if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
    #endif
        {
            playSpeed = playSpeed + changeSpeed;
        }
        //slow down if left arrow is pressed
    #ifdef WIN32
        if (GetAsyncKeyState(37) & 0x8000)
    #else
        if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
    #endif
        {
            playSpeed = playSpeed - changeSpeed;
        }

        //if 1 is pressed
    #ifdef WIN32
        if (GetAsyncKeyState(49) & 0x8000)
    #else
        if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
    #endif
        {
            //set the play speed to 1
            playSpeed = 1;
        }

        //disable / enable dark matter with Z
    #ifdef WIN32
        if (GetAsyncKeyState(90) & 0x8000)
    #else
        if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS)
    #endif
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            showDarkMatter = !showDarkMatter;
        }

        //disable / enable dark matter with U
    #ifdef WIN32
        if (GetAsyncKeyState(85) & 0x8000)
    #else
        if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS)
    #endif
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            if (colorMode < 2) colorMode++;
            else colorMode = 0;
        }
    }
    oldIndex = index;
}

double Engine::calcDensityAv()
{
    //calc the middle value of the particles densities
    double densityAv = 0;
    for (const auto& particle : *particles)
    {
        densityAv += particle->density;
    }
    densityAv = densityAv / particles->size();
    return densityAv;
}

vec3 jetColorMap(double value) 
{
    double r = value < 0.5 ? 0.0 : 2.0 * value - 1.0;
    double g = value < 0.25 ? 0.0 : value < 0.75 ? 1.0 : 1.0 - 4.0 * (value - 0.75);
    double b = value < 0.5 ? 1.0 : 1.0 - 2.0 * value;

    return vec3(r, g, b);
}



void Engine::renderParticles()
{
    // Binden des Framebuffers
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Deaktivieren Sie den Tiefentest und das Z-Buffering
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);

    // In der renderParticles-Funktion
    glUseProgram(shaderProgram);

    // Erstellen der Projektionsmatrix und Sichtmatrix
    mat4 projection = mat4::perspective(45.0f, 800.0f / 600.0f, 0.1f, cameraViewDistance);
    mat4 viewMatrix = mat4::lookAt(cameraPosition, cameraPosition + cameraFront, cameraUp);

    // Setzen der Matrizen im Shader
    GLuint projectionLoc = glGetUniformLocation(shaderProgram, "projection");
    GLuint viewLoc = glGetUniformLocation(shaderProgram, "view");

    glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, projection.data());
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, viewMatrix.data());

    // Vertex Array Object (VAO) binden
    glBindVertexArray(VAO);

    if (BGstars)
    {
        //render the background bgStars
        for (int i = 0; i < amountOfStars; i++)
        {
            glPointSize(bgStars[i].w);

            // Setzen der Position im Shader
            vec3 pos(bgStars[i].x, bgStars[i].y, bgStars[i].z);
            float posArray[3];
            pos.toFloatArray(posArray);
            glUniform3fv(glGetUniformLocation(shaderProgram, "particlePosition"), 1, posArray);

            // Setzen der Farbe im Shader
            vec3 color(1.0f, 1.0f, 1.0f);
            float colorArray[3];
            color.toFloatArray(colorArray);
            glUniform3fv(glGetUniformLocation(shaderProgram, "particleColor"), 1, colorArray);

            // Zeichnen des Punktes
            glDrawArrays(GL_POINTS, 0, 1);
        }
    }

    for (const auto& particle : *particles)
    {
        double red = 1;
        double green = 1;
        double blue = 1;

        vec3 color = vec3(red, green, blue);

        if (densityAv != 0) 
        {
           color.x = particle->density * 10 / densityAv;
           color.y = 0;
           color.z = densityAv / particle->density;

           double plus = 0;
           if(particle->density * 100 / densityAv > 1) plus = particle->density / densityAv / 10;

           color.x += plus * 10;
           color.y += plus;
           color.z += plus;
        }

        vec3 scaledPosition = particle->position * globalScale;

        // Setzen Position im Shader
        float scaledPosArray[3];
        scaledPosition.toFloatArray(scaledPosArray);
        glUniform3fv(glGetUniformLocation(shaderProgram, "particlePosition"), 1, scaledPosArray);

        // Setzen Farbe im Shader
        float colorArray[3];
        color.toFloatArray(colorArray);
        glUniform3fv(glGetUniformLocation(shaderProgram, "particleColor"), 1, colorArray);

        glPointSize(0.5f); // Setzen der Punktgröße auf 5 Pixel

        // Zeichnen Punkt
        glDrawArrays(GL_POINTS, 0, 1);
    }

    // VAO lösen
    glBindVertexArray(0);
}

void Engine::processInput()
{
    if (RenderLive)
    {
        // Toggle zwischen Kameramodi mit "L"
        static bool lKeyWasPressedLastFrame = false;

        // Prüfen, ob die Taste "L" gerade gedrückt wurde
        bool lKeyPressed =
#ifdef WIN32
        (GetAsyncKeyState('L') & 0x8000) != 0;
#else
            glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS;
#endif

        // Umschalten des Kameramodus nur, wenn die Taste neu gedrückt wurde
        if (lKeyPressed && !lKeyWasPressedLastFrame)
        {
            focusedCamera = !focusedCamera;
        }

        // Aktualisieren des letzten Tastendruckzustands
        lKeyWasPressedLastFrame = lKeyPressed;

        if (RenderLive)
        {
#ifdef WIN32
            if (GetAsyncKeyState(VK_CONTROL) < 0)
#else
            if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
#endif
            {
                // Kamerabewegung
                float index = 0.1f;
                if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
                    cameraPosition += rushSpeed * index * cameraFront;
                if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
                    cameraPosition -= rushSpeed * index * cameraFront;
                if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
                    cameraPosition -= rushSpeed * index * cameraFront.cross(cameraUp);
                if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
                    cameraPosition += rushSpeed * index * cameraFront.cross(cameraUp);
            }
            else
            {
                // Kamerabewegung
                float index = 0.1f;
                if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
                    cameraPosition += cameraSpeed * index * cameraFront;
                if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
                    cameraPosition -= cameraSpeed * index * cameraFront;
                if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
                    cameraPosition -= cameraSpeed * index * cameraFront.cross(cameraUp);
                if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
                    cameraPosition += cameraSpeed * index * cameraFront.cross(cameraUp);
            }
        }
    }
    if (RenderLive == false)
    {
        float index = 0.1f; 
        cameraPosition += cameraSpeed * index * cameraFront.cross(cameraUp);
        //std::cout << "Camera Position: " << cameraPosition.x << " " << cameraPosition.y << " " << cameraPosition.z << std::endl;
    }

    
    if (focusedCamera) {
        // Richtet die Kamera auf den Ursprung aus
        vec3 direction = (vec3(0, 0, 0) - cameraPosition).normalize();
        cameraFront = direction;
    }
    // Aktualisieren der Ansichtsmatrix (Kameraposition und Blickrichtung)
    view = mat4::lookAt(cameraPosition, cameraPosition + cameraFront, cameraUp);

    // Setzen der Matrizen im Shader
    GLuint viewLoc = glGetUniformLocation(shaderProgram, "view");
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, view.data());
}

void Engine::processMouseInput()
{
    if (focusedCamera == false)
    {
        // Erfassen Sie die Mausbewegung
        double mouseX, mouseY;
        glfwGetCursorPos(window, &mouseX, &mouseY);

        static double lastMouseX = mouseX;
        static double lastMouseY = mouseY;

        double xOffset = mouseX - lastMouseX;
        double yOffset = lastMouseY - mouseY; // Umgekehrtes Vorzeichen für umgekehrte Mausrichtung

        lastMouseX = mouseX;
        lastMouseY = mouseY;

        const float sensitivity = 0.05f;
        xOffset *= sensitivity;
        yOffset *= sensitivity;

        cameraYaw += xOffset;
        cameraPitch += yOffset;

        // Begrenzen Sie die Kamerapitch, um ein Überdrehen zu verhindern
        if (cameraPitch > 89.0f)
            cameraPitch = 89.0f;
        if (cameraPitch < -89.0f)
            cameraPitch = -89.0f;

        // Berechnen der neuen Kamerarichtung
        vec3 newFront;
        newFront.x = cos(radians(cameraYaw)) * cos(radians(cameraPitch));
        newFront.y = sin(radians(cameraPitch));
        newFront.z = sin(radians(cameraYaw)) * cos(radians(cameraPitch));
        cameraFront = newFront.normalize();

        // Aktualisieren der Ansichtsmatrix
        view = mat4::lookAt(cameraPosition, cameraPosition + cameraFront, cameraUp);
    }
}


void Engine::checkShaderCompileStatus(GLuint shader, const char* shaderType)
{
    GLint success;
    GLchar infoLog[512];
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

    if (!success) {
        glGetShaderInfoLog(shader, 512, NULL, infoLog);
        std::cerr << "Error compiling " << shaderType << " shader:\n" << infoLog << std::endl;
    }
}

void Engine::checkShaderLinkStatus(GLuint program)
{
    GLint success;
    GLchar infoLog[512];
    glGetProgramiv(program, GL_LINK_STATUS, &success);

    if (!success) {
        glGetProgramInfoLog(program, 512, NULL, infoLog);
        std::cerr << "Error linking shader program:\n" << infoLog << std::endl;
    }
}

bool Engine::clean()
{
    // Aufräumen und beenden
    glfwTerminate();
    return true;
}

void Engine::calcTime(int index) 
{
    calcTime(vec3(0, 0, 0), index);
}

void Engine::calcTime(vec3 position, int index)
{
    passedTime = (index * faktor) * deltaTime;
    //std::cout << deltaTime << std::endl;

    int passedTimeInSec = passedTime / 86400;

    std::string Unit;

    //set the right unit
    if (passedTime < 60) { Unit = " s"; }
    else if (passedTime < 3600) { passedTime /= 60; Unit = " min"; }
    else if (passedTime < 86400) { passedTime /= 3600; Unit = " h"; }
    else if (passedTime < 31536000) { passedTime /= 86400; Unit = " days"; }
    else { passedTime /= 31536000; Unit = " years"; }

    std::string time; 
    //set the time to like millions, billions, trillions, ...
    if (passedTime < 1000) 
    { 
        //passed time to string with 1 decimal after the point
        std::string newTime = std::to_string(passedTime);
        newTime = newTime.substr(0, newTime.find(".") + 2);
        time = newTime; 
    }
	else if (passedTime < 1000000) 
    { 
        passedTime /= 1000; 
        //passed time to string with 1 decimal after the point
        std::string newTime = std::to_string(passedTime);
        newTime = newTime.substr(0, newTime.find(".") + 2);
        time = newTime + " thousand"; 
    }
	else if (passedTime < 1000000000) 
    { 
        passedTime /= 1000000; 
		//passed time to string with 1 decimal after the point
		std::string newTime = std::to_string(passedTime);
		newTime = newTime.substr(0, newTime.find(".") + 2);
		time = newTime + " million";
    }
	else if (passedTime < 1000000000000) 
    {
        passedTime /= 1000000000;
        //passed time to string with 1 decimal after the point
        std::string newTime = std::to_string(passedTime);
        newTime = newTime.substr(0, newTime.find(".") + 2);
        time = newTime + " billion";
    }
	else if (passedTime < 1000000000000000) 
    { 
        passedTime /= 1000000000000; 
		//passed time to string with 1 decimal after the point
		std::string newTime = std::to_string(passedTime);
		newTime = newTime.substr(0, newTime.find(".") + 2);
		time = newTime + " trillion";
    }
	else 
    { 
        passedTime /= 1000000000000000;
        //passed time to string with 1 decimal after the point
        std::string newTime = std::to_string(passedTime);
        newTime = newTime.substr(0, newTime.find(".") + 2);
        time = newTime + " quadrillion";
    }

    std::cout << "Passed time: " << time << Unit << std::endl;
}

double Engine::random(double min, double max)
{
    // Generiere eine zufällige Gleitkommazahl zwischen min und max
    double randomFloat = min + static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * (max - min);

    return randomFloat;
}

void Engine::calculateGlobalScale()
{
    double maxDistance = 0;
    for (const auto& particle : *particles)
    {
        // Skip ungültige oder extreme Werte
        if (particle->position.x == 0 && particle->position.y == 0 && particle->position.z == 0)
            continue;
        if (std::isinf(particle->position.x) || std::isinf(particle->position.y) || std::isinf(particle->position.z) ||
            std::isnan(particle->position.x) || std::isnan(particle->position.y) || std::isnan(particle->position.z))
            continue;

        double distance = sqrt(pow(particle->position.x, 2) + pow(particle->position.y, 2) + pow(particle->position.z, 2));
        if (distance > maxDistance)
            maxDistance = distance;
    }

    // Vermeide die Verarbeitung, wenn alle Positionen ungültig sind
    if (maxDistance == 0) {
        globalScale = 1; // Standardwert oder Fehlerwert setzen
        return;
    }

    maxDistance = maxDistance * 2; // Durchmesser des Systems
    int exponent = static_cast<int>(std::floor(std::log10(std::abs(maxDistance))));
    globalScale = maxDistance / pow(10, exponent * 2 - 2);
}
