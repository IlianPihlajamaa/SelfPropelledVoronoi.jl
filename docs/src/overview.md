# Overview

`SelfPropelledVoronoi.jl` is a Julia package designed for simulating two-dimensional (2D) systems of self-propelled particles with periodic boundary conditions. The central idea is that particle interactions and spatial organization are determined by a dynamic Voronoi tessellation of the simulation domain. Each particle corresponds to a Voronoi cell, and the properties of these cells (like area and perimeter) influence particle behavior.

## Core Mechanics

The behavior of particles in `SelfPropelledVoronoi.jl` is governed by a few key mechanical and kinetic principles:

1.  **Voronoi Tessellation:** The simulation space is continuously partitioned into Voronoi cells, with each cell uniquely associated with a particle. The geometry of these cells (area, perimeter, number of neighbors) determines the forces acting on the cells.

2.  **Area and Perimeter Elasticity:**
    *   Each particle (cell) has a preferred or "target" area (`A0`) and a target perimeter (`p0`).
    *   Deviations from these target values result in an elastic energy penalty. For example, if a cell's actual area `A` is different from `A0`, there is an energy contribution proportional to `KA * (A - A0)^2`, where `KA` is an area stiffness constant. A similar term applies for perimeter deviations with a perimeter stiffness `KP`.
    *   These energy penalties translate into forces that drive the cells to relax towards their target geometries, influencing particle movement and rearrangement.

3.  **Active Particle Motion (Self-Propulsion):**
    *   Particles are "active," meaning they generate their own motion. Each particle `i` has an intrinsic self-propulsion force of magnitude `f0_i`.
    *   This force is directed along an orientation vector `(cos(θ_i), sin(θ_i))`, where `θ_i` is the particle's current orientation angle.

4.  **Rotational Diffusion:**
    *   The orientation `θ_i` of each particle is not fixed but evolves stochastically over time due to rotational diffusion.
    *   This is modeled by adding a random angular displacement at each time step, with the rate of change controlled by a rotational diffusion constant `Dr_i`.

5.  **Periodic Boundary Conditions:**
    *   The simulation is performed in a rectangular domain with periodic boundary conditions. This means that particles exiting one side of the box re-enter from the opposite side, allowing for the simulation of bulk system behavior without edge effects.

By combining these elements, `SelfPropelledVoronoi.jl` provides a tool to explore how local interactions and active forces at the particle level give rise to complex, large-scale behaviors in 2D active systems. The reference to Bi et al. (Physical Review X 6.2 (2016): 021011) provides an introduction to this model.
