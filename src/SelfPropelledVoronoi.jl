
"""
This module implements a simulation of self-propelled Voronoi tessellations.
It defines data structures for representing the simulation state, parameters, and output.
It also provides functions for initializing the simulation, calculating forces, updating particle positions, and performing Voronoi tesselation.

The main exported function is `run_simulation!`, which executes the simulation.
Exported data structures include `SimulationBox`, `ParameterStruct`, `VoronoiCells`, `ArrayStruct`, `Output`, `DumpInfo`, and `VoronoiNeighborList`.
"""
module SelfPropelledVoronoi
import LoopVectorization, Quickhull, Random, HDF5
using SmallCollections: MutableSmallVector
using StaticArrays: SVector

for file in [
    "DataStructs.jl", 
    "AuxiliaryFunctions.jl",
    "InitialConfiguration.jl", 
    "Dump.jl", 
    "Forces.jl",
    "Tesselation.jl",
    "Dynamics.jl", 

    ] 
    
    include(file) 
end

export SimulationBox, ParameterStruct, VoronoiCells, ArrayStruct, Output, DumpInfo, VoronoiNeighborList
export run_simulation!


end


