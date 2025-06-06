# Examples

The `example/test.jl` script provides a starting point for understanding how to set up and run a basic simulation. 

Running the code typically consists of the following steps: initialize a small number of particles, define simulation parameters, and run the simulation for a set number of steps.

Here's a conceptual example of how to use the main module and run a simulation:

```julia
# Import necessary packages
using SelfPropelledVoronoi # The main package for the simulation
using StaticArrays # Provides SVector for efficient, statically-sized vectors (used here for 2D positions)
using Random # Provides random number generation capabilities, including MersenneTwister

# --- Define Simulation Domain and Particle Count ---
# N: Number of particles in the simulation
N = 20 
# box_Lx: Length of the simulation box in the x-dimension
box_Lx = 10.0 
# box_Ly: Length of the simulation box in the y-dimension
box_Ly = 10.0 
# sim_box: Creates a SimulationBox object to store the dimensions of the 2D simulation domain. 
# This object encapsulates the geometry (Lx, Ly) and implies periodic boundary conditions are used by the simulation logic.
sim_box = SimulationBox(box_Lx, box_Ly)

# --- Define Particle Properties (VoronoiCells) ---
# These parameters define the physical characteristics and behavior of the simulated particles, 
# which are modeled as cells in a Voronoi tessellation.
# Each parameter is an array of length N, allowing for different properties for each particle if needed.
# The VoronoiCells struct groups these properties.

# p0: Target perimeter for each particle (Voronoi cell). 
# Cells will experience an elastic force resisting deviations from this preferred perimeter.
# units: length units, consistent with box_Lx, box_Ly.
p0 = 3.8 * ones(N) # 'ones(N)' creates a vector of N elements, all set to 1.0, then scaled by 3.8.

# A0: Target area for each particle (Voronoi cell). 
# Similar to p0, deviations from this preferred area will result in an elastic restoring force.
# units:  area units (length units squared).
A0 = 1.0 * ones(N)

# KP: Perimeter spring constant (or stiffness). 
# Determines the strength of the force resisting changes to the cell's perimeter. Higher KP means more resistance.
# units:  energy / length^2 or force / length.
KP = 1.0 * ones(N)

# KA: Area spring constant (or stiffness). 
# Determines the strength of the force resisting changes to the cell's area. Higher KA means more resistance.
# units:  energy / area^2 or force / area.
KA = 1.0 * ones(N)

# f0: Active force strength (or self-propulsion force magnitude). 
# The magnitude of the intrinsic force that drives each particle's movement. The direction is determined by the orientation of the particle.
# units:  force units.
f0 = 1.0 * ones(N)

# Dr: Rotational diffusion rate. 
# Controls the rate at which the orientation of particles changes randomly due to rotational noise.
# units: radians^2 / time 
Dr = 0.1 * ones(N)

# particles: Creates a VoronoiCells object to store the collective properties of all particles.
# The arguments correspond to target_perimeters, target_areas, K_P, K_A, active_force_strengths, and rotational_diffusion_constants respectively.
particles = VoronoiCells(p0, A0, KP, KA, f0, Dr)

# --- Define Simulation Parameters ---
# params: A ParameterStruct structure holding various global simulation settings and constants.
# The constructor uses keyword arguments.
params = ParameterStruct(
    # Number of particles (already defined)
    N = N, 
    # Time step for the numerical integration of equations of motion. 
    # Smaller values generally lead to higher accuracy but increase computation time.
    # units: time units.
    dt = 0.001, 
    # Total number of simulation steps to perform. The total simulation time will be dt * N_steps.
    N_steps = 1000, 
     # Thermal energy, product of Boltzmann constant (kB) and Temperature (T). 
    # This term scales the magnitude of random forces due to thermal fluctuations
    # units: energy units.
    kBT = 1.0,
    # Friction coefficient (gamma). Affects the mobility of the cells (the ratio of the force on the cell and its velocity)
    frictionconstant = 1.0, 
    # Depth of the layer used for handling interactions across periodic boundaries.
    # Particles within this distance from a boundary will interact with periodic images of particles 
    # from the opposite side of the simulation box. This value should typically be larger than any interaction cutoffs.
    # units: length units.
    periodic_boundary_layer_depth = 2.5, 
    # If true, the simulation prints progress messages, warnings, or other information to the console.
    verbose = true, 
    # The SimulationBox object defining the domain (defined above).
    box = sim_box, 
    # The VoronoiCells object defining particle properties (defined above).
    particles = particles, 
    # A DumpInfo structure for controlling data output (e.g., saving simulation snapshots). 
    # Here, 'save=false' indicates that no data will be written to disk during this example run.
    # Default DumpInfo settings might otherwise enable saving.
    dump_info = DumpInfo(save=false), 
    # A user-defined function that is called at the start of every time step during the simulation.
    # 'p', 'a', 'o' typically refer to the parameters, arrays, and output structs, respectively.
    # Useful for custom actions like live plotting, on-the-fly analysis, or complex data collection. 
    # Here, it's a no-operation (no-op) function, meaning nothing custom is done during simulation steps.
    callback = (p,a,o) -> nothing, 
    # Random Number Generator (RNG) instance.
    # `MersenneTwister` is a widely used pseudo-random number generator known for its good statistical properties.
    # Initializing it with a specific seed (1234 here) ensures that the sequence of random numbers
    # generated will be the same every time the simulation is run with this seed. This is crucial for reproducibility.
    RNG = Random.MersenneTwister(1234) 
)

# --- Initialize Particle Positions and Orientations ---
# arrays: An ArrayStruct structure to hold the dynamic state variables of particles (e.g., positions, orientations) 
# that change over time during the simulation.
arrays = ArrayStruct(N) 
# Initialize positions and orientations randomly:
for i in 1:N
    # arrays.positions[i]: Stores the 2D position (x, y coordinates) of particle 'i'.
    # `SVector` (StaticVector) from the `StaticArrays.jl` package is used here. 
    # SVectors are stack-allocated, fixed-size arrays that can provide performance benefits for small vectors 
    # (like 2D or 3D coordinates) due to better memory locality and enabling certain compiler optimizations.
    # Positions are initialized randomly within the simulation box dimensions (box_Lx, box_Ly),
    # using the provided RNG from the params struct for reproducibility.
    arrays.positions[i] = SVector(rand(params.RNG)*box_Lx, rand(params.RNG)*box_Ly)
    # arrays.orientations[i]: Stores the orientation of particle 'i'. 
    # This is represented as an angle in radians (e.g., from 0 to 2*pi).
    # The orientation determines the direction of the self-propulsion force f0.
    # Orientations are initialized randomly between 0 and 2*pi using the provided RNG.
    arrays.orientations[i] = 2pi * rand(params.RNG)
end

# --- Initialize Output Structure ---
# output: An Output structure designed to store summary statistics 
# accumulated during the simulation (e.g., total potential energy, steps completed).
output = Output()

# --- Run the Simulation ---
# run_simulation!: This function executes the main simulation loop.
# It takes the simulation parameters (params), the initial state of particle arrays (arrays), 
# and the output structure (output) as input.
# The '!' at the end of the function name is a Julia convention indicating that 
# the function is likely to modify one or more of its arguments in-place. In this case, 
# 'arrays' (particle positions/orientations) will be updated at each step, and 
# 'output' will be populated with simulation results like `steps_done`.
run_simulation!(params, arrays, output)

# --- Print Simulation Summary ---
# After the simulation completes, print a message indicating completion and 
# the total number of steps that were actually performed (which should be stored in output.steps_done).
println("Simulation finished. Steps done: $(output.steps_done)")
```
