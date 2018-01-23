import Base: copy, show

type CountLog
    edges_inserted::Int64
    edges_deleted::Int64
    edges_updated::Int64
    edges_traversed::Int64
    degrees_read::Int64
    degrees_written::Int64
end

CountLog() = CountLog(0, 0, 0, 0, 0, 0)

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
    T, g::SimpleWeightedDiGraph, b::Vector{Int64}, num_blocks::Int64;
    count_log = CountLog()
    )
    M = initialize_edge_counts(T, g, num_blocks, b, count_log)
    d_out, d_in, d = compute_block_degrees(M, num_blocks, count_log)
    overall_entropy::Float64 = compute_overall_entropy(
        M, d_out, d_in, num_blocks, nv(g), ne(g), count_log
    )
    Partition(M, g, overall_entropy, b, d, d_out, d_in, num_blocks, count_log)
end

function copy(p::Partition)
    Partition(
        copy(p.M), p.g, p.S, copy(p.b), copy(p.d), copy(p.d_out), copy(p.d_in),
        p.B, copy(count_log)
    )
end

show(io::IO, p::Partition) = print(io, "Partition(", p.S, "," , p.B, ",", sort(unique(p.b)),")")
