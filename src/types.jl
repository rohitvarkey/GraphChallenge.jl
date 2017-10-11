import Base: copy, show

type Partition{T}
    M::T
    S::Float64
    b::Vector{Int64}
    d::Vector{Int64}
    d_out::Vector{Int64}
    d_in::Vector{Int64}
    B::Int64
end

copy(p::Partition) = Partition(p.M, p.S, p.b, p.d, p.d_out, p.d_in, p.B)
show(io::IO, p::Partition) = print(io, "Partition(", p.S, "," , p.B, ")")
