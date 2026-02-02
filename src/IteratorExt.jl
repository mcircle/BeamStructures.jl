
struct BeamsAtNode{T}
    adj::AbstractMatrix{T} #Dict{CartesianIndex,T}
    nodes::NamedTuple
    function BeamsAtNode(adj::AbstractMatrix{T},nodes::NamedTuple) where{T}
        new{T}(adj,nodes)
    end 
end 

beamsatnode(itr,nodes) = BeamsAtNode(itr,nodes)  

Base.eltype(::Type{BeamsAtNode{T}}) where {T} = AbstractMatrix{T}
Base.length(r::BeamsAtNode) = length(r.nodes)
Base.size(r::BeamsAtNode) = (length(r),)
IteratorSize(r::Type{BeamsAtNode{T}}) where {T} = length(r)
IteratorEltype(::Type{BeamsAtNode{T}}) where {T} = IteratorEltype(T)

function Base.iterate(r::BeamsAtNode)
    adj = r.adj
    ci = findall(x->!isapprox(x,0),LowerTriangular(adj))
    cbeams = count(x->!isapprox(x,0),adj[2:end,1]) 

    node = 1
    if isa(r.nodes[node],Branch)
        beams = (filter(x->ci[x][1] == node,1:cbeams),filter(x->ci[x][2] == node,1:cbeams))
        # cbeams += count(x->!isapprox(x,0),adj[3:end,2])
    else 
        beams = (Int[],filter(x->ci[x][2] == node,1:cbeams))
    end
    state = (ci,node + 1,cbeams)

    return ((r.nodes[node],beams),state)
end

function Base.iterate(r::BeamsAtNode, state)
    it = r.adj
    ci, node, cbeams = state
    if node > length(r)
        return nothing
    end

    cbeams += count(x->!isapprox(x,0),it[ node+1:end ,node])
    if isa(r.nodes[node],Branch)
        beams = (filter(x->ci[x][1] == node,1:cbeams),filter(x->ci[x][2] == node,1:cbeams))
        
    else 
        beams = (Int[],filter(x->ci[x][2] == node,1:cbeams))
    end

    state = (ci,node + 1,cbeams)
    return ((r.nodes[node],beams),state)
end

# it = beamsatnode(adj,bnstr.Nodes)
# out,state = iterate(it)
# state
# iterate(it,state)

# for (node,beams) in  beamsatnode(adj,bnstr.Nodes)
#     println("Node: $node has beams: $beams")
# end
