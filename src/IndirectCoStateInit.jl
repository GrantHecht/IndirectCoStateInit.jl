module IndirectCoStateInit

using Heuristics
using StaticArrays

# Includes
include("CoStateInitializer.jl")
include("HeuristicsCoStateInitializer.jl")
include("FSSCoStateInitializer.jl")

# Exports
export FSSCoStateInitializer
export initialize!
export GetInitializedCostates

end
