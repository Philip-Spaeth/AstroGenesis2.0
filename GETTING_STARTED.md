# Getting Started with Astro Genesis 2.0

This guide will help you set up and compile Astro Genesis 2.0 using CMake and GCC on both Windows and Linux.

## Prerequisites

Before you begin, ensure you have the following software installed on your system:

- CMake (version 3.10 or higher)
- GCC (GNU Compiler Collection)

### On Windows

1. **Install CMake**:
   - Download the latest version of CMake from the [CMake website](https://cmake.org/download/).
   - Run the installer and follow the instructions.

2. **Install GCC**:
   - Install [MinGW-w64](http://mingw-w64.org/doku.php/download).
   - Add the `bin` directory of your MinGW installation to your system's PATH environment variable.

### On Linux

1. **Install CMake and GCC**:
   - On Debian-based systems (e.g., Ubuntu), run:
     ```sh
     sudo apt update
     sudo apt install cmake gcc g++
     ```
   - On Red Hat-based systems (e.g., Fedora), run:
     ```sh
     sudo dnf install cmake gcc gcc-c++
     ```

## Cloning the Repository

First, clone the repository to your local machine:
   ```sh
   git clone https://github.com/yourusername/astro-genesis-2.0.git
   cd astro-genesis-2.0
```


Navigate to the project directory:
   ```sh
    cd astro-genesis-2.0
```

Build the project using CMake:
   ```sh
   mkdir build
   cd build
   cmake ..
   make
```

If you encounter any issues during the build process, refer to the Troubleshooting section.
