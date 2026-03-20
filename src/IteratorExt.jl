

struct IterateThrough
    size::Int
    function IterateThrough(sz::Int)
        new(sz)
    end
end
iteratethrough(sz::Int) = IterateThrough(sz)

Base.eltype(::Type{IterateThrough}) = Int
Base.length(r::IterateThrough) = r.size * (r.size - 1) ÷ 2
Base.size(r::IterateThrough) = (length(r),)
IteratorSize(r::Type{IterateThrough}) = length(r)
IteratorEltype(::Type{IterateThrough}) = IteratorEltype(Int)

function Base.iterate(r::IterateThrough,state::NTuple{2,Int} = (1,1))
    all(state .>= r.size -1) && return nothing
    if state[1] == r.size
        state = (state[2] + 2, state[2] + 1)
    else
        state = (state[1] + 1, state[2])
    end

    return state, state
end

get_zero(::Type{T}) where {T} = eps(T)
get_zero(::Type{T}) where {T<:Int} = zero(T)


function nonzero_lower_indices(adj::AbstractMatrix{T}) where T
    n = min(size(adj,1), size(adj,2))
    idxs = Vector{CartesianIndex{2}}()  
    sizehint!(idxs, size(adj,1) * (size(adj,1) -1) ÷ 2)
    for (i,j) in IterateThrough(n)
        @inbounds if adj[i,j] > get_zero(T)
            push!(idxs, CartesianIndex(i,j))
        end
    end
    return idxs
end

function nonzero_lower_indices(adj::AbstractMatrix{T},sh) where T
    n = size(adj,1)
    idxs = Vector{CartesianIndex{2}}(undef,sh)  
    idx = 1 
    for (i,j) in IterateThrough(n)
        @inbounds if adj[i,j] > get_zero(T)
            idxs[idx] =  CartesianIndex(i,j)
            idx += 1
        end
    end
    return idxs
end

struct BeamsAtNode{N}
    nodes::NTuple{N,Symbol} 
    beams_to_node::Vector{Vector{Int}}
    beams_from_node::Vector{Vector{Int}}
    function BeamsAtNode(nodes::NamedTuple{names,Q},beams_to_node::Vector{Vector{Int}},beams_from_node::Vector{Vector{Int}}) where{Q,names}
        new{length(names)}(names,beams_to_node,beams_from_node)
    end
end

function BeamsAtNode(adj::AbstractMatrix{T},nodes::NamedTuple{B,Q},beams::BT) where{T,N,B,BT<:NamedTuple,Q<:NTuple{N,Boundary}}
    # @show typeof(B),typeof(Q)
    n = length(nodes)
    beamsfromnode = [Int[] for _ in 1:n]
    beamstonode = [Int[] for _ in 1:n]

    idcs = length(beams) < size(adj,1) * (size(adj,1)-1 ) ÷ 2 ? nonzero_lower_indices(adj,length(beams)) : IterateThrough(size(adj,1))

    for (i, idx) in enumerate(idcs)
        to_node   = idx[1]
        from_node = idx[2]
        push!(beamsfromnode[from_node], i)
        push!(beamstonode[to_node], i)
    end

    BeamsAtNode(nodes,beamstonode,beamsfromnode)
end

beamsatnode(itr,nodes,beams) = BeamsAtNode(itr,nodes,beams)  

Base.eltype(::Type{BeamsAtNode}) = Int
Base.length(::BeamsAtNode{N}) where{N} = N
Base.size(r::BeamsAtNode) = (length(r),)
Base.IteratorSize(::Type{<:BeamsAtNode}) = Base.HasLength()
Base.IteratorEltype(::Type{<:BeamsAtNode}) = Base.HasEltype()
IteratorSize(r::Type{BeamsAtNode}) = length(r)
Base.eltype(::Type{BeamsAtNode{N}}) where {N} = Tuple{Symbol, Tuple{Vector{Int}},Vector{Int}}

function Base.iterate(r::BeamsAtNode,state = 1) 
    state > length(r) && return nothing
    @inbounds node = r.nodes[state]
    @inbounds beams_to = r.beams_to_node[state]
    @inbounds beams_from = r.beams_from_node[state]
    return ((node, (beams_from, beams_to)), state + 1)
end

