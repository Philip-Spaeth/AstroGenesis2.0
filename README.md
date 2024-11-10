# Astro Genesis 2.0

Welcome to Astro Genesis 2.0, a modern n-body simulation software designed to model the dynamics of galaxies and clusters. Developed in C++, Astro Genesis 2.0 leverages cutting-edge algorithms and data structures to deliver highly accurate simulations of cosmic phenomena.

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Simulation](#simulation)
- [Contributing](#contributing)
- [License](#license)
- [Examples](#Examples)

## Features
- **n-body Simulation**: Simulate complex gravitational interactions between a large number of objects with high precision.
- **Octree-based Force Calculation**: Efficient force calculation by dividing space into a tree structure (Octree).
- **Smoothed Particle Hydrodynamics (SPH)**: Models fluid-like behavior of interstellar gas.
- **MPI Support for Parallel Computing**: Optimized for server, enabling parallel calculations to simulate large systems efficiently.
- **Basic build-in Initial Conditions (ICs) Generator**: Generates stable initial conditions for galaxies, based on the solution of the Jeans equation.
  - **Bulge**: Models the galactic bulge using the Hernquist profile.
  - **Disk**: Generates an exponential disk for the galaxy structure, including a gas disk model.
  - **Halo**: Simulates the dark matter halo using either the Hernquist profile or the NFW (Navarro-Frenk-White) profile.
- **Radiative Cooling**: Integrates radiative cooling to realistically model the thermodynamics of interstellar gas and promote star formation.
- **Star Formation**: Simulates the process of star formation from dense, cold gas, taking into account radiative cooling and other astrophysical processes.
- **Supernova Feedback**: Accounts for feedback from supernova explosions on the interstellar medium, influencing gas dynamics and star formation.
- **HDF5 Data Storage**: Stores simulation data in an efficient HDF5 format, enabling easy analysis and visualization of large datasets.
  
## Installation
See the GetttingStarted.md file for further information. 
  
## Simulation
The Simulation Engine is responsible for calculating the n-body interactions. It uses an octree structure to optimize the computation of gravitational forces, making it possible to simulate millions of particles efficiently.

## Contributing
We welcome contributions from the community! Whether you are fixing bugs, adding new features, or improving documentation, your help is greatly appreciated. Please submit pull requests to the dev branch and ensure your code follows the project's coding standards.

## License
Astro Genesis 2.0 is released under the GNU 3.0 License. See 'LICENSE' for more information.
