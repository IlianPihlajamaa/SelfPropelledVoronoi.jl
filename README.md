# SelfPropelledVoronoi

[![Build Status](https://github.com/IlianPihlajamaa/SelfPropelledVoronoi.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/IlianPihlajamaa/SelfPropelledVoronoi.jl/actions/workflows/CI.yml?query=branch%3Amaster)

## High-level Overview

`SelfPropelledVoronoi.jl` is a Julia package for simulating systems of self-propelled particles. The interactions between these particles are governed by a Voronoi tessellation of the simulation space. Key features of the simulation include:

*   **Area and Perimeter Elasticity:** Particles (modeled as Voronoi cells) have target areas and perimeters. Deviations from these targets result in elastic energy penalties, influencing particle motion.
*   **Active Particle Motion:** Particles are self-propelled, meaning they have an intrinsic force that drives them in a direction determined by their orientation. This orientation can change due to rotational diffusion.
*   **Periodic Boundary Conditions:** The simulation domain is typically a rectangular box with periodic boundary conditions, allowing for simulations of bulk systems.
*   **Voronoi Tessellation:** The core of the interaction model relies on dynamically computing the Voronoi tessellation of the particle positions at each relevant step.

## Installation

Currently, `SelfPropelledVoronoi.jl` is used as a local package. To use it, you first need to ensure its dependencies are installed in your Julia environment.

```julia
# From the Julia REPL
import Pkg

# Core dependencies (based on imports in SelfPropelledVoronoi.jl)
Pkg.add("StaticArrays")
Pkg.add("Quickhull")
Pkg.add("HDF5") # For potential data saving
Pkg.add("LoopVectorization") # For performance via @turbo
Pkg.add("Random") # For random number generation
Pkg.add("SmallCollections") # For MutableSmallVector

# If SelfPropelledVoronoi.jl were a registered package:
# Pkg.add("SelfPropelledVoronoi")

# If it were a proper local package you want to develop:
# Pkg.develop(path="path/to/SelfPropelledVoronoi.jl")

# For now, to run examples:
# 1. Navigate to the root directory of this project in your terminal.
# 2. Launch Julia: `julia`
# 3. Include an example script, e.g.:
#    include("example/test.jl")
```

## Basic Usage / Example

The `example/test.jl` script provides a starting point for understanding how to set up and run a basic simulation. It typically initializes a small number of particles, defines simulation parameters, and runs the simulation for a set number of steps.

Here's a conceptual example of how to use the main module and run a simulation:

```julia
using SelfPropelledVoronoi
using StaticArrays # For SVector
using Random # For MersenneTwister

# Define parameters (refer to DataStructs.jl for details)
N = 20
box_Lx = 10.0
box_Ly = 10.0
sim_box = SimulationBox(box_Lx, box_Ly)

# Example particle parameters (VoronoiCells)
p0 = 3.8 * ones(N)  # Target perimeter
A0 = 1.0 * ones(N)  # Target area
KP = 1.0 * ones(N)  # Perimeter spring constant
KA = 1.0 * ones(N)  # Area spring constant
f0 = 1.0 * ones(N)  # Active force strength
Dr = 0.1 * ones(N)  # Rotational diffusion

particles = VoronoiCells(p0, A0, KP, KA, f0, Dr)

params = ParameterStruct(
    N = N,
    dt = 0.001,
    N_steps = 1000,
    kBT = 0.0, # No thermal noise in this example
    frictionconstant = 1.0,
    periodic_boundary_layer_depth = 2.5, # Should match Tesselation.jl or be configurable
    verbose = true,
    box = sim_box,
    particles = particles,
    dump_info = DumpInfo(save=false), # No saving in this example
    callback = (p,a,o) -> nothing, # No-op callback
    RNG = Random.MersenneTwister(1234)
)

# Initialize arrays
arrays = ArrayStruct(N)
# Simple random initial positions and orientations
for i in 1:N
    arrays.positions[i] = SVector(rand(params.RNG)*box_Lx, rand(params.RNG)*box_Ly)
    arrays.orientations[i] = 2pi * rand(params.RNG)
end
arrays.old_positions .= arrays.positions
arrays.old_orientations .= arrays.orientations

# Initialize output struct
output = Output()

# Run the simulation
run_simulation!(params, arrays, output)

println("Simulation finished. Steps done: $(output.steps_done)")
```

## Structure of the Code

The simulation logic is organized into several Julia files within the `src/` directory:

*   `SelfPropelledVoronoi.jl`: The main module file that includes other components and exports key functions and types.
*   `DataStructs.jl`: Defines the core data structures used throughout the simulation, such as `ParameterStruct`, `ArrayStruct`, `VoronoiCells`, `SimulationBox`, `Output`, and `DumpInfo`.
*   `AuxiliaryFunctions.jl`: Provides various helper functions, including those for applying periodic boundary conditions, calculating pair distances, updating cell perimeters and areas, and computing system energy.
*   `Tesselation.jl`: Contains the logic for performing Voronoi tessellation. This includes updating particle positions with periodic images, using `Quickhull.jl` for Delaunay triangulation, calculating circumcenters (Voronoi vertices), and sorting vertices.
*   `Dynamics.jl`: Implements the time evolution of the system. It includes functions for force calculation (dispatching to specific models), numerical integration schemes (e.g., Euler-Maruyama, Euler-Heun), and the main simulation loop (`run_simulation!`).
*   `Forces.jl`: Intended for defining specific force laws. Currently, it might contain placeholder or example forces (e.g., for the SPV model based on area/perimeter elasticity, or simpler pair potentials like Gaussian Core Model).
*   `InitialConfiguration.jl`: Responsible for setting up the initial state of the simulation (e.g., particle positions, orientations). May currently be underdeveloped.
*   `Dump.jl`: Handles saving simulation state and data to files. May currently have placeholder functionality.

## To Do / Future Work

*   **Configuration Files:** Implement a system for loading simulation parameters from configuration files (e.g., TOML, JSON) instead of hardcoding them in scripts.
*   **Initial Conditions:** Develop more sophisticated methods for generating initial particle configurations (e.g., specific packing densities, different spatial distributions).
*   **Data Saving:** Fully implement data saving functionality in `Dump.jl`, likely using HDF5 for efficient storage of time series data.
*   **Force Models:** Expand `Forces.jl` to include a wider variety of inter-particle interaction models or more complex active matter behaviors.
*   **Tessellation Optimization:**
    *   Properly implement `verify_tesselation` to avoid unnecessary full re-tessellations.
    *   Implement `update_delauney_vertices!` for incremental updates to the tessellation, which could offer significant performance benefits when particle displacements are small.
*   **Parameterization:** Ensure `periodic_boundary_layer_depth` is consistently used from `ParameterStruct` rather than being hardcoded in `Tesselation.jl`.
*   **Testing & Validation:** Add comprehensive unit and integration tests. Validate simulation results against known models or theoretical predictions where possible.
*   **Performance Profiling & Optimization:** Profile the code to identify bottlenecks and apply further optimizations.
*   **Package Registration:** If the package matures, consider registering it for easier installation via the Julia package manager.

```
