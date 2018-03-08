import Base: copy, show, convert

using DataFrames

type CountLog
    edges_inserted::Int64
    edges_deleted::Int64
    edges_updated::Int64
    edges_traversed::Int64
    edges_read::Int64
    degrees_read::Int64
    degrees_written::Int64
end

CountLog() = CountLog(0, 0, 0, 0, 0, 0, 0)

type Partition{T}
    M::T
    g::SimpleWeightedDiGraph
    S::Float64
    b::Vector{Int64}
    d::Vector{Int64}
    d_out::Vector{Int64}
    d_in::Vector{Int64}
    B::Int64
    count_log::CountLog
end


function Partition(
    T, g::SimpleWeightedDiGraph, b::Vector{Int64}, num_blocks::Int64, config;
    count_log = CountLog()
    )
    M = initialize_edge_counts(T, g, num_blocks, b, config, count_log)
    d_out, d_in, d = compute_block_degrees(M, num_blocks, count_log)
    overall_entropy::Float64 = compute_overall_entropy(
        M, d_out, d_in, num_blocks, nv(g), ne(g), count_log
    )
    Partition(M, g, overall_entropy, b, d, d_out, d_in, num_blocks, count_log)
end

function Partition(
    T, g::SimpleWeightedDiGraph, b::Vector{Int64}, num_blocks::Int64, config::Void;
    count_log = CountLog()
    )
    M = initialize_edge_counts(T, g, num_blocks, b, count_log)
    d_out, d_in, d = compute_block_degrees(M, num_blocks, count_log)
    overall_entropy::Float64 = compute_overall_entropy(
        M, d_out, d_in, num_blocks, nv(g), ne(g), count_log
    )
    Partition(M, g, overall_entropy, b, d, d_out, d_in, num_blocks, count_log)
end

function compute_block_degrees(p::Partition, count_log::CountLog)
    p.d_out, p.d_in, p.d = compute_block_degrees(p.M, p.B, count_log)
end

function compute_overall_entropy(p::Partition, count_log::CountLog)
    p.S = compute_overall_entropy(
        p.M, p.d_out, p.d_in, p.B, nv(p.g), ne(p.g), count_log
    )
end

function copy(p::Partition)
    Partition(
        copy(p.M), p.g, p.S, copy(p.b), copy(p.d), copy(p.d_out), copy(p.d_in),
        p.B, copy(count_log)
    )
end

show(io::IO, p::Partition) = print(io, "Partition(", p.S, "," , p.B, ",", p.b,")")

abstract type Metrics end

immutable CorrectnessMetrics <: Metrics
    accuracy::Float64
    pairwise_precision::Float64
    pairwise_recall::Float64
    adj_rand_index::Float64
    rand_index::Float64
end

immutable PerformanceMetrics <: Metrics
    time::Float64
    bytes::Int64
    gctime::Float64
    memallocd::Int64
end

immutable BenchmarkMetrics <: Metrics
    T::DataType
    nv::Int64
    ne::Int64
    nb::Int64
    performance::PerformanceMetrics
    correctness::CorrectnessMetrics
    count_log::CountLog
end

function convert(::Type{DataFrame}, results::Vector{BenchmarkMetrics})
    df = DataFrame(
        Type = [result.T for result in results],
        num_nodes = [result.nv for result in results],
        num_edges = [result.ne for result in results],
        num_blocks = [result.nb for result in results]
    )
    for sub_metric in [:performance, :correctness, :count_log]
        t = fieldtype(BenchmarkMetrics, sub_metric)
        for field in fieldnames(t)
            df[field] = [getfield(getfield(result, sub_metric), field) for result in results]
        end
    end
    df
end
