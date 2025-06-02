# SelfPropelledVoronoi

[![Build Status](https://github.com/IlianPihlajamaa/SelfPropelledVoronoi.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/IlianPihlajamaa/SelfPropelledVoronoi.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://IlianPihlajamaa.github.io/SelfPropelledVoronoi.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://IlianPihlajamaa.github.io/SelfPropelledVoronoi.jl/dev)



## High-level Overview

`SelfPropelledVoronoi.jl` is a Julia package for simulating systems of confluent layers of self-propelled particles. The interactions between these particles are governed by a Voronoi tessellation of the simulation space. Key features of the simulation include:

*   **Area and Perimeter Elasticity:** Particles (modeled as Voronoi cells) have target areas and perimeters. Deviations from these targets result in elastic energy penalties, influencing particle motion.
*   **Active Particle Motion:** Particles are self-propelled, meaning they have an intrinsic force that drives them in a direction determined by their orientation. This orientation can change due to rotational diffusion.
*   **Periodic Boundary Conditions:** The simulation domain is a rectangular box with periodic boundary conditions, allowing for simulations of bulk systems.

See Bi, Dapeng, et al. "Motility-driven glass and jamming transitions in biological tissues." Physical Review X 6.2 (2016): 021011.

https://github.com/user-attachments/assets/6732d9ff-9e69-454e-b4be-e90b93d48148

## Installation

Currently, `SelfPropelledVoronoi.jl` is not registered. To install it, run:

```julia
# From the Julia REPL
import Pkg; activate(".")
Pkg.add("https://github.com/IlianPihlajamaa/SelfPropelledVoronoi.jl")
```

## Basic Usage / Example

Running the code typically consists of the following steps: initialize a small number of particles, define simulation parameters, and run the simulation for a set number of steps.

Below is a conceptual example of how to use the main module and run a simulation. See also the [Documentation](https://IlianPihlajamaa.github.io/SelfPropelledVoronoi.jl/dev).

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
    kBT = 1.0, 
    frictionconstant = 1.0,
    periodic_boundary_layer_depth = 2.5, 
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

# Initialize output struct
output = Output()

# Run the simulation
run_simulation!(params, arrays, output)

println("Simulation finished. Steps done: $(output.steps_done)")
```

