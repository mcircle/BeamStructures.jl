
function getnames(nodes::Vararg{<:Boundary,N}) where{N}
    NamedTuple{ntuple(i->Symbol(:Node_,i),N)}(nodes)
end 

function getnames(beams::Vararg{<:BeamElement,N}) where{N}
    NamedTuple{ntuple(i->Symbol(:Beam_,i),N)}(beams)
end 

function prepare(args::Vararg{Union{<:BeamElement,<:Boundary},N}) where{N}
    beams = filter(x->isa(x,BeamElement),args)
    bounds = filter(x->isa(x,Boundary),args)
    getnames(beams...), getnames(bounds...)        
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

function getnodeswithbeams(adj::AbstractMatrix,nodes::NamedTuple)
    nodeswithbeams = Vector{Int}()
    for ap in axes(adj,1)
        if isa(nodes[ap],Branch) 
            push!(nodeswithbeams,ap)
        elseif isa(nodes[ap],Clamp) && ap > 1 && sum(adj[1:ap,ap]) > 0 
            push!(nodeswithbeams,ap)
        end
    end 
    return nodeswithbeams
end   


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

function gaussfilter(x,μ = 0,σ = 1f-1)    
   @. exp(-((x - μ) / σ)^2 /2) #./ (σ * sqrt(2π))
end

function softmax(x;dims = 1)
    tmp = zero(x)
    for (tm,xt) in zip(eachslice(tmp;dims = dims),eachslice(x,;dims = dims))
        tm .= exp.(xt) ./ (sum(exp,xt) .+ eps(eltype(x)))
    end
    tmp
end


@non_differentiable get_index_pars(x,inds)

switchCI(a::CartesianIndex{2}) = CartesianIndex((a[2],a[1]))
@non_differentiable switchCI(a)

findbeamsatnode(::Clamp,node::Int,nodes::AbstractVector{CartesianIndex{2}}) = (findall(x-> x[1] == node,nodes),Vector{Int}())
findbeamsatnode(::Boundary,node::Int,nodes::AbstractVector{CartesianIndex{2}}) = (findall(x->x[1] == node,nodes),findall(x-> x[2] == node,nodes))
normfactor(b::BeamElement) = 12 / (b.E*b.w*b.h^3) # 12/(E*I)
normfactor_m(b::BeamElement) = normfactor(b) * b.l #M̃ = M * normvector_m
normfactor_f(b::BeamElement) = normfactor_m(b) * b.l #F̃ = F * normvector_f

function normvector(b::BeamElement)
    m = normfactor_m(b)                  # normfactor_m
    f = b.l * m                          # normfactor_f
    return @SVector [m, f, f]    
end

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
function getforceindices(beams)
    getindices(beams,[1,5,6])
end 

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

function DiffEqBase.promote_u0(u0::AbstractArray{T,N},p::ODESolution{T,N1},t::T) where {T,N,N1}
    u0
end 

function DiffEqBase.anyeltypedual(sol::ODESolution,::Val{counter}) where counter
    diffeqmapreduce(anyeltypedual, promote_dual, (sol.u, sol.t))
end

function softplus(x,β = 1f0)
    @. log(1 + exp(β * x)) / β
end

function inverse_softplus(x,β = 1f0)
    @. log(exp(β * x) - 1) / β
end

