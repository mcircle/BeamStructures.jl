
function getnames(nodes::Vararg{Boundary,N}) where{N}
    NamedTuple{ntuple(i->Symbol(:Node_,i),N)}(nodes)
end 

function getnames(beams::Vararg{Beam,N}) where{N}
    NamedTuple{ntuple(i->Symbol(:Beam_,i),N)}(beams)
end 

function prepare(args::Vararg{Union{Beam,Boundary},N}) where{N}
    beams = filter(x->isa(x,Beam),args)
    bounds = filter(x->isa(x,Boundary),args)
    (;Beams = getnames(beams...),Nodes = getnames(bounds...))        
end 

function getstartnodes(str::Structure)
    adj = str.AdjMat
    getindex.(findall(x->!isapprox(x,0),LowerTriangular(adj)),2)
end 
function getstartnodes(adj::AbstractMatrix{T}) where{T}
    sz = size(adj,1)-1
    len = reduce(+,1:sz)
    nodes = ones(Int,len)
    ind = 0
    for i in sz:-1:2
        ind += i
        nodes[ind+1:end] .+= 1
    end 
    return    nodes
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

function gaussfilter(x,μ = 0.,σ = 0.1)    
   @. exp(-0.5 * ((x - μ) / σ)^2) #./ (σ * sqrt(2π))
end

function softmax(x,f = 1)
    f .* exp.(x) ./ (sum(exp,x) + eps(Float32))
end


@non_differentiable get_index_pars(x,inds)

switchCI(a::CartesianIndex{2}) = CartesianIndex((a[2],a[1]))
@non_differentiable switchCI(a)

findbeamsatnode(::Clamp,node::Int,nodes::AbstractVector{CartesianIndex{2}}) = (findall(x-> x[1] == node,nodes),Vector{Int}())
findbeamsatnode(::Branch,node::Int,nodes::AbstractVector{CartesianIndex{2}}) = (findall(x->x[1] == node,nodes),findall(x-> x[2] == node,nodes))

normfactor_m(b::Beam) = 12* b.l/(b.E*b.w*b.h^3) 
normfactor_f(b::Beam) = b.l * normfactor_m(b) 
normvector(b::Beam) = [normfactor_m(b),normfactor_f(b),normfactor_f(b)]

getside(x,beams) = ifelse(x ∈ beams[1],2,1)

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

function getindices(s::Int)
    cis = Vector{CartesianIndex{2}}(undef,reduce(+,1:s-1))
    start = 1
    for ind in 2:s
        cis[start:start + s - ind] = collect(CartesianIndices((ind:s,ind-1:ind-1)))
        start += s - ind + 1
    end 
    cis
end 
@non_differentiable getindices(s)
@non_differentiable getpositionindices(b)
@non_differentiable getforceindices(b)

reducevcat(res::AbstractVector{T}) where{T<:AbstractVector{TT}} where{TT} = reduce(vcat,res;init = Vector{TT}())

Base.eltype(::Beam{T}) where{T} = T

function getfieldof(beams,s::Symbol)
    T = eltype(eltype(beams))
    vec = Vector{T}(undef,length(beams))
    for ind in eachindex(vec)
        vec[ind] = getfield(beams[1],s)
    end 
    vec
end  
