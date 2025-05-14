
module SelfPropelledVoronoi
import LoopVectorization, Quickhull, Random, HDF5, StaticArrays

for file in [
    "DataStructs.jl", 
    "AuxiliaryFunctions.jl",
    "InitialConfiguration.jl", 
    "Dynamics.jl", 
    "Dump.jl", 
    "Forces.jl"
    ] 
    
    include(file) 
end

end
