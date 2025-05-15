
module SelfPropelledVoronoi
import LoopVectorization, Quickhull, Random, HDF5
using SmallCollections: MutableSmallVector
using StaticArrays: SVector

for file in [
    "DataStructs.jl", 
    "AuxiliaryFunctions.jl",
    "InitialConfiguration.jl", 
    "Dynamics.jl", 
    "Dump.jl", 
    "Forces.jl",
    "Tesselation.jl"
    ] 
    
    include(file) 
end

export SimulationBox, ParameterStruct, VoronoiCells, ArrayStruct, Output, DumpInfo


end


