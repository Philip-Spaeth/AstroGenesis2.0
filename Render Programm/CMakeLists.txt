cmake_minimum_required(VERSION 3.10)
project(MeinProjekt)

set(CMAKE_CXX_STANDARD 17)

set(MAIN_SOURCE "src/main.cpp")

# Inkludieren Sie die OpenGL-Bibliothek
find_package(OpenGL REQUIRED)

# Funktion, um rekursiv alle Verzeichnisse zu finden
function(add_subdirectories_recursive dir)
    file(GLOB children RELATIVE ${dir} ${dir}/*)
    foreach(child ${children})
        if(IS_DIRECTORY ${dir}/${child})
            include_directories(${dir}/${child})
        endif()
    endforeach()
endfunction()

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
    link_libraries(${GLEW_LIBRARY} ${GLFW_LIBRARY})

elseif(UNIX AND NOT APPLE)
    # Inkludieren Sie GLM
    include_directories(${CMAKE_SOURCE_DIR}/External_Libraries_Linux/GLM/glm)
    # Inkludieren und Verknüpfen Sie GLEW
    find_package(GLEW REQUIRED)
    include_directories(${CMAKE_SOURCE_DIR}/External_Libraries_Linux/GLEW/include)
    link_directories(${CMAKE_SOURCE_DIR}/External_Libraries_Linux/GLEW/lib)
    # Inkludieren und Verknüpfen Sie GLFW
    find_package(glfw3 REQUIRED)
    include_directories(${CMAKE_SOURCE_DIR}/External_Libraries_Linux/GLFW/include)
    link_directories(${CMAKE_SOURCE_DIR}/External_Libraries_Linux/GLFW/lib) # Passen Sie dies an den tatsächlichen Pfad Ihrer GLFW-Bibliotheken an
endif()

# Define GLM_ENABLE_EXPERIMENTAL to use experimental GLM features
add_definitions(-DGLM_ENABLE_EXPERIMENTAL)

include_directories(${CMAKE_SOURCE_DIR}/GLApp/GLApp/Physics/)
include_directories(${CMAKE_SOURCE_DIR}/GLApp/GLApp/File/)
include_directories(${CMAKE_SOURCE_DIR}/GLApp/GLApp/Physics/BarnesHut/)
include_directories(${CMAKE_SOURCE_DIR}/GLApp/GLApp/Physics/Galaxies/)
include_directories(${CMAKE_SOURCE_DIR}/GLApp/GLApp/Rendering)
include_directories(${CMAKE_SOURCE_DIR}/GLApp/GLApp/Math/)
include_directories(${CMAKE_SOURCE_DIR})


add_executable(MeinProgramm ${MAIN_SOURCE} ${HEADER_FILES} ${SOURCE_FILES})

target_link_libraries(MeinProgramm PRIVATE ${OPENGL_LIBRARIES} ${GLFW_LIBRARY} ${GLEW_LIBRARY} pthread)

# Find Curses library
find_package(Curses REQUIRED)
if(CURSES_FOUND)
    message(STATUS "Found Curses: ${CURSES_LIBRARY}")
    target_link_libraries(MeinProgramm PRIVATE ${CURSES_LIBRARIES})
endif()

# Hinzufügen der DLLs zum Ausführungsverzeichnis
if(WIN32)
    add_custom_command(TARGET MeinProgramm POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_SOURCE_DIR}/External Libraries/GLFW/lib-mingw-w64/glfw3.dll"
        $<TARGET_FILE_DIR:MeinProgramm>)
    add_custom_command(TARGET MeinProgramm POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_SOURCE_DIR}/External Libraries/GLEW/bin/Release/x64/glew32.dll"
        $<TARGET_FILE_DIR:MeinProgramm>)
endif()
