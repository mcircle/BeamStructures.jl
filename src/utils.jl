
function getnames(nodes::Vararg{Boundary,N}) where{N}
    NamedTuple{ntuple(i->Symbol(:Node_,i),N)}(nodes)
end 

function getnames(beams::Vararg{Beam,N}) where{N}
    NamedTuple{ntuple(i->Symbol(:Beam_,i),N)}(beams)
end 

gettype(::T) where{T<:Boundary} = T

function prepare(args::Vararg{Union{Beam,Boundary},N}) where{N}
    beams = filter(x->isa(x,Beam),args)
    bounds = filter(x->isa(x,Boundary),args)
    (;Beams = getnames(beams...),Nodes = getnames(bounds...))        
end 

function getstartnodes(adj::AbstractMatrix)
    getindex.(findall(x->!isapprox(x,0),LowerTriangular(adj)),2)
end 

@non_differentiable getstartnodes(adj)

indexlength(::Branch) = 6
indexlength(::Boundary) = 3

getindexfor(::Boundary) = [3,2,4]

function get_index_pars(nodes,inds)
    ind = 1
    for node in inds[1:end-1]
        ind += indexlength(nodes[node])
    end  
    return ind:ind + indexlength(nodes[inds[end]]) - 1
end 

@non_differentiable get_index_pars(x,inds)

switchCI(a::CartesianIndex{2}) = CartesianIndex((a[2],a[1]))
@non_differentiable switchCI(a)

findbeamsatnode(::Clamp,node::Int,nodes::AbstractVector{CartesianIndex{2}}) = (findall(x-> x[1] == node,nodes),Vector{Int}())
findbeamsatnode(::Branch,node::Int,nodes::AbstractVector{CartesianIndex{2}}) =  (findall(x->x[2] == node,nodes),findall(x-> x[1] == node,nodes))

normfactor_m(b::Beam) = 12* b.l/(b.E*b.w*b.h^3) 
normfactor_f(b::Beam) = b.l * normfactor_m(b) 
normvector(b::Beam) = [normfactor_m(b),normfactor_f(b),normfactor_f(b)]

getside(x,beams) = ifelse(x ∈ beams[1],1,2)

getnextidxs(::Clamp,beams) = 3 * sum(length.(beams))
getnextidxs(::Boundary,beams) = 3 * (sum(length.(beams)))

function getindices(beams::NTuple{2,Vector{Int}},pos)
    bnbrs = vcat(beams...)
    idxs = Matrix{CartesianIndex{3}}(undef,length(pos),length(bnbrs))
    for (ind,b) in enumerate(bnbrs)
        side = getside(b,beams)
        for p in eachindex(idxs[:,ind])
            idxs[p,ind] = CartesianIndex(pos[p],side,b)
        end 
    end 
    idxs
end 
getpositionindices(beams) = getindices(beams,[2,3,4])
getforceindices(beams) = getindices(beams,[1,5,6])

