cmake_minimum_required(VERSION 3.10)
project(AstroGenesis_Render_Project)

set(CMAKE_CXX_STANDARD 17)

set(MAIN_SOURCE "src/main.cpp")

# Inkludieren Sie die OpenGL-Bibliothek
find_package(OpenGL REQUIRED)

if(WIN32)
    # Inkludieren Sie GLM
    include_directories("${CMAKE_SOURCE_DIR}/External Libraries/GLM/glm")

    # Manuelles Hinzufügen der GLFW-Bibliotheken und Include-Verzeichnisse
    set(GLFW_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/External Libraries/GLFW/include")
    set(GLFW_LIBRARY "${CMAKE_SOURCE_DIR}/External Libraries/GLFW/lib-mingw-w64/libglfw3.a")
    include_directories(${GLFW_INCLUDE_DIR})
    link_directories("${CMAKE_SOURCE_DIR}/External Libraries/GLFW/lib-mingw-w64")

    # Inkludieren und Verknüpfen Sie GLEW
    set(GLEW_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/External Libraries/GLEW/include")
    set(GLEW_LIBRARY "${CMAKE_SOURCE_DIR}/External Libraries/GLEW/lib/Release/x64/glew32.lib")
    include_directories(${GLEW_INCLUDE_DIR})
    link_directories("${CMAKE_SOURCE_DIR}/External Libraries/GLEW/lib/Release/x64")

    # Andere Bibliotheken verknüpfen
    list(APPEND ADDITIONAL_LIBRARIES ${GLEW_LIBRARY} ${GLFW_LIBRARY} ole32 uuid)
endif()

# Define GLM_ENABLE_EXPERIMENTAL to use experimental GLM features
add_definitions(-DGLM_ENABLE_EXPERIMENTAL)

# include all cpp and .h files
file(GLOB_RECURSE HEADER_FILES "src/*.h")
set(SOURCE_FILES
    src/main.cpp
    src/Engine.cpp
    src/vec3.cpp
    src/mat4.cpp
    src/vec4.cpp
    src/DataManager.cpp
)

include_directories(${CMAKE_SOURCE_DIR}/src/)

add_executable(AstroGenesis_Render_Programm ${MAIN_SOURCE} ${HEADER_FILES} ${SOURCE_FILES})

# Link libraries
target_link_libraries(AstroGenesis_Render_Programm PRIVATE ${OPENGL_LIBRARIES} ${ADDITIONAL_LIBRARIES} pthread)

# Hinzufügen der DLLs zum Ausführungsverzeichnis
if(WIN32)
    add_custom_command(TARGET AstroGenesis_Render_Programm POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_SOURCE_DIR}/External Libraries/GLFW/lib-mingw-w64/glfw3.dll"
        $<TARGET_FILE_DIR:AstroGenesis_Render_Programm>)
    add_custom_command(TARGET AstroGenesis_Render_Programm POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_SOURCE_DIR}/External Libraries/GLEW/bin/Release/x64/glew32.dll"
        $<TARGET_FILE_DIR:AstroGenesis_Render_Programm>)
endif()
