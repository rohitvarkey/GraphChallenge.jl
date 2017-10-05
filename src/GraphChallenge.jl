module GraphChallenge

using SimpleWeightedGraphs
using LightGraphs
using Distributions

export load_graph!, initialize_edge_counts, compute_block_degrees, partition

# package code goes here
const INPUT_PATH = joinpath(dirname(dirname(@__FILE__)), "data", "streaming")

include("load_graph.jl")
include("interblock_matrix.jl")
include("partition.jl")
end # module
