module GraphChallenge
export load_graph!, initialize_edge_counts, partition

# package code goes here
const INPUT_PATH = joinpath(dirname(dirname(@__FILE__)), "data", "streaming")

include("load_graph.jl")
include("interblock_matrix.jl")
include("partition.jl")
end # module
