cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(AstroGen2)

# Setzen Sie den Ausgabeordner für Bibliotheken
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Setzen Sie den Zwischenspeicherordner für Objektdateien
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/obj)

# Alle .cpp Dateien rekursiv finden
file(GLOB_RECURSE USER_SOURCE_FILES
    src/*.cpp
)

# Funktion, um rekursiv alle Verzeichnisse zu finden
function(add_subdirectories_recursive dir)
    file(GLOB children RELATIVE ${dir} ${dir}/*)
    foreach(child ${children})
        if(IS_DIRECTORY ${dir}/${child})
            include_directories(${dir}/${child})
        endif()
    endforeach()
endfunction()

# Hauptverzeichnis und alle Unterverzeichnisse hinzufügen
include_directories(${PROJECT_SOURCE_DIR}/src)
add_subdirectories_recursive(${PROJECT_SOURCE_DIR}/src)

find_package(OpenMP REQUIRED)

# ===================== HighFive Integration Start =====================

# Fügen Sie das HighFive-Submodul-Verzeichnis hinzu
add_subdirectory(extern/HighFive)

# Finden der HDF5-Bibliotheken (systemweit installiert)
find_package(HDF5 REQUIRED COMPONENTS C CXX)

# ===================== HighFive Integration Ende =====================

# Fügen Sie Ihre ausführbare Datei hinzu
add_executable(${PROJECT_NAME} ${USER_SOURCE_FILES})

if(OpenMP_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

# Verknüpfen Sie HDF5 und HighFive mit Ihrem Programm
target_link_libraries(${PROJECT_NAME} PRIVATE HDF5::HDF5 HighFive)

# Compiler-Optionen und Link-Bibliotheken
target_compile_options(${PROJECT_NAME} PRIVATE -Wall)

# Optional: Weitere Bibliotheken verknüpfen
# target_link_libraries(${PROJECT_NAME} ws2_32 gdi32)

# Define custom target to clean build directory
add_custom_target(clean-build
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}
)

# Define custom target to build using an external terminal
add_custom_target(build-in-external-terminal
    COMMAND start cmd /c "${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --config Debug"
    DEPENDS ${PROJECT_NAME}
)
