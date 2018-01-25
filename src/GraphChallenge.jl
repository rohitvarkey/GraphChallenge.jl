module GraphChallenge

using SimpleWeightedGraphs
using LightGraphs
using Distributions
using StatsBase

export load_graph!, load_graph, initialize_edge_counts, compute_block_degrees,
       partition, evaluate_partition, plotbeliefs, static_partition_experiment,
       InterblockEdgeCountDictDict, InterblockEdgeCountVectorDict,
       InterblockEdgeCountStinger, InterblockEdgeCountPostgres,
       CountLog, Partition

# package code goes here
const INPUT_PATH = joinpath(dirname(dirname(@__FILE__)), "data")

include("types.jl")
include("load_graph.jl")
include("interblock_matrix.jl")
include("dictofdicts.jl")
include("vectorofdicts.jl")
include("stinger.jl")
include("postgres.jl")
include("partition.jl")
#include("plot.jl")
include("evaluate.jl")
include("run_experiment.jl")

end # module
