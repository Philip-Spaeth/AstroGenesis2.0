# AstroGen2 (Astro Genesis 2.0)

Welcome to Astro Genesis 2.0, a modern n-body simulation software designed to model the dynamics of galaxies and clusters. Developed in C++, Astro Genesis 2.0 leverages optimized algorithms and data structures to deliver highly accurate simulations of cosmic phenomena.

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Simulation](#simulation)
- [Contributing](#contributing)
- [License](#license)

## Features
- **n-body Simulation**: Simulate complex gravitational interactions between a large number of objects with high precision.
- **Octree-based Force Calculation**: Efficient force calculation by dividing space into a tree structure (Octree).
- **Smoothed Particle Hydrodynamics (SPH)**: Models fluid-like behavior of interstellar gas.
- **Radiative Cooling**: Integrates radiative cooling to realistically model the thermodynamics of interstellar gas and promote star formation.
- **Star Formation**: Simulates the process of star formation from dense, cold gas.
- **Supernova Feedback**: Accounts for feedback from supernova explosions on the interstellar medium, influencing gas dynamics and star formation.
- **HDF5 Data Storage**: Stores simulation data in an efficient HDF5 format, enabling easy analysis and visualization of large datasets.
  
## Installation
See the GetttingStarted.md file for further information. 
  
## Simulation
The Simulation Engine is responsible for calculating the n-body interactions.
See /docs folder for the Users-guide, code paper and information about Dataformats etc.

## Contributing
We welcome contributions from the community! Whether you are fixing bugs, adding new features, or improving documentation, your help is greatly appreciated. Please submit pull requests to the dev branch and ensure your code follows the project's coding standards.

## License
Astro Genesis 2.0 is released under the GNU 3.0 License. See 'LICENSE' for more information.
