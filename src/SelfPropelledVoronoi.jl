

module SelfPropelledVoronoi
import  Quickhull, Random, HDF5
using SmallCollections: MutableSmallVector
using StaticArrays: SVector
using GLMakie: Axis, linesegments!, display, Figure, scatter!, @lift, Observable, lines!

for file in [
    "DataStructs.jl", 
    "AuxiliaryFunctions.jl",
    "InitialConfiguration.jl", 
    "Dump.jl",
    "Load.jl",
    "Forces.jl",
    "Tesselation.jl",
    "Dynamics.jl", 
    "Visualize.jl",    ] 
    
    include(file) 
end


export SimulationBox, ParameterStruct, VoronoiCells, ArrayStruct, Output, DumpInfo, VoronoiNeighborList
export run_simulation!
export compute_energy
export compute_pair_distance_vector


end
