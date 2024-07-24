# Astro Genesis 2.0

Welcome to Astro Genesis 2.0, a modern n-body simulation software designed to model the dynamics of galaxies and clusters. Developed in C++, Astro Genesis 2.0 leverages cutting-edge algorithms and data structures to deliver highly accurate simulations of cosmic phenomena, including the influence of dark energy and dark matter.

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
  - [Simulation](#simulation)
  - [Render Engine](#render-engine)
- [Live Viewer](#live-viewer)
- [Contributing](#contributing)
- [License](#license)
  
## Features
- **Advanced n-body Simulation**: Simulate complex gravitational interactions with high precision.
- **Octree Data Structure**: Efficiently manage spatial data for large-scale simulations.
- **Dark Energy and Dark Matter Modeling**: Incorporate the effects of dark energy and dark matter into your simulations.
- **Smoothed Particle Hydrodynamics (SPH)**: Model fluid-like behavior of interstellar gas and other particles.
- **Dual Engine Architecture**: Separate Simulation and Render engines for optimized performance.
- **Video Rendering**: Render high-quality videos of your simulations.
- **Live Viewer**: Interactively explore simulations in real-time.
- **C++ Core**: High performance and flexibility.

## Installation
To install Astro Genesis 2.0, follow these steps:

1. Clone the repository:
   ```sh
   git clone https://github.com/yourusername/astro-genesis-2.0.git
   
2. Navigate to the project directory:
   ```sh
    cd astro-genesis-2.0

3. Build the project using CMake:
   ```sh
   mkdir build
   cd build
   cmake ..
   make
  
## Usage
## Simulation
The Simulation Engine is responsible for calculating the n-body interactions. It uses an octree structure to optimize the computation of gravitational forces, making it possible to simulate millions of particles efficiently.

## Render Engine
The Render Engine can generate high-quality video outputs from simulation data or provide a live viewer mode for real-time exploration. The rendering process incorporates advanced techniques to visualize the complex interactions between particles, including the influence of dark matter and dark energy.

## Live Viewer
The Live Viewer allows users to interactively explore the simulation in real-time. This tool is invaluable for gaining insights into the dynamic behaviors of galaxies and clusters. Users can pan, zoom, and rotate the view to observe the simulation from different angles and perspectives.

## Contributing
We welcome contributions from the community! Whether you are fixing bugs, adding new features, or improving documentation, your help is greatly appreciated. Please submit pull requests to the dev branch and ensure your code follows the project's coding standards.

## License
Astro Genesis 2.0 is released under the GNU 3.0 License. See LICENSE for more information.
