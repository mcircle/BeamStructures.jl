const startof = 1
const endof  = 2

struct Connections{A<:Dict{Int,Vector},N<:Dict{Int,Vector{Pair}},B<:Dict{Int,Pair},T<:Tuple{Vararg{Boundary}}} <: AbstractMatrix{Int}
    #"Adjacencematrix"
    AdjMat::A
    Nodes2Beams::N # stores Nodes connected to Beam and at which end
    Beams2Nodes::B # stores Beams connected to Nodes at start/end
    Nodes::T
    function Connections(a::A,n::N,b::B,t::T) where{A<:Dict{Int,Vector},N<:Dict{Int,Vector{Pair}},B<:Dict{Int,Pair},T<:Tuple{Vararg{Boundary}}}
        new{A,N,B,T}(a,n,b,t)
    end
end 

function get_branches(con::Connections)
    br = findall(x->isa(con.Nodes[x],Branch),1:length(con.Nodes))
    count_br = length(br)
    Dict([br[x] => x for x in 1:count_br])
end 
"""
types: Dict of Int => Symbol to describe the type of Node: clamp,joint,branch,slider
"""
function Connections(adj::AbstractMatrix{T},tps) where{T}
    @assert length(tps) == size(adj,2)
    adjmat = Dict{Int,Vector}()
    nodes = Dict{Int,Vector{Pair}}([i => [] for i in 1:size(adj,2)])
    Beams = Dict{Int,Pair}()
    b = 1 
    for (n,col) in enumerate(eachcol(adj))
        #n = startnode 
        #findall attached nodes 
        adjmat[n] = findall(x->x > 0,col[n:end]) .+n .-1
        #node = node at end of beam
        for node in adjmat[n]
            Beams[b] = n => node
            append!(nodes[n],[b => startof])
            append!(nodes[node],[b=> endof])
            b +=1
        end 
    end 
    Connections(adjmat,nodes,Beams,tps)
end 

function Connections(adj::AbstractMatrix,tps...)
    Connections(adj,tps)
end 

function incidence(con::Connections)
    #each row = Node 
    mat = zeros(eltype(con),length(con.Nodes2Beams),length(con.Beams2Nodes))
    for (b,(st,en)) in con.Beams
        mat[[st,en],b] .= 1
    end 
    mat
end 

#erhalte Adjacencematrix von Kantengraph
edge_adjacence(ct) = incidence(ct)' * incidence(ct) .- Matrix(2*I,length(ct.Beams),length(ct.Beams))
Base.size(A::Connections) = Tuple(length(A.AdjMat) .* ones(Int,2))
Base.getindex(A::Connections{OT,N,B,T},i,j) where{OT,N,B,T} = j ∈ A.AdjMat[i] || i ∈ A.AdjMat[j] ? 1 : 0

function get_neighbors(c::Connections,node)
    findall(x->x == 1,c[:,node])
end 

function Adj_norm(adj::Matrix{T},selfatt = true) where {T}
    ai = adj .+ (selfatt ? Matrix{T}(I,size(adj)...) : zero(T)) #adding selfattation
    D = sum(ai,dims = 1) .* I(size(ai,1))
    Dinv = sqrt(inv(D)) 
    Dinv * ai * Dinv 
end 

