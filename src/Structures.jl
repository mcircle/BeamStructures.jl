
struct Structure{A<:AbstractMatrix,So,Se,KW}
    AdjMat::A
    Solver::So
    SensAlg::Se
    kwargs::KW
    function Structure(adj::A,solver::So,sensalg::Se,kwargs::KW) where{A,So,Se,KW}
        new{A,So,Se,KW}(adj,solver,sensalg,kwargs)
    end 
end

function Structure(adj::AbstractMatrix{T};solver = Tsit5(),sensalg = ForwardDiffSensitivity(),kwargs...) where{T}
    pos = vec(any(x->!isapprox(x,zero(T)),adj,dims = 1))
    Structure(adj[pos,pos],solver,sensalg,kwargs)
end 

function Structure(adj::Connections;solver = Tsit5(),sensalg = ForwardDiffSensitivity(),kwargs...) 
    Structure(adj,solver,sensalg,kwargs)
end 

function Base.show(io::IO,str::Structure)
    return println(io, "Structure with $(size(str.AdjMat,1)) Nodes and $(Int(sum(str.AdjMat)/2)) Beam(s).")
end 

function initialize_beam(node::Boundary,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0, = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    m,fx,fy,Δθ,Δx,Δy = parameters
    m *= normfactor_m(beam) #am Balkenelement
    fx *= normfactor_f(beam) #am Balkenelement
    fy *= normfactor_f(beam) #am Balkenelement
    return [m,θ0 + θs + Δθ,(x + Δx)./l,(y + Δy)./l,fx,fy,κ*l]
end 

function initialize_beam(node::CompliantClamp,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0, = node.x,node.y,node.ϕ
    l,κ = beam.l,beam.κ0
    m = node.c * (beam.θs .- parameters[1]) .* normfactor_m(beam)
    fx = parameters[2] * normfactor_f(beam) #am Balkenelement
    fy = parameters[2] * normfactor_f(beam) #am Balkenelement
    return [m,θ0 + beam.θs + parameters[1],x./l,y./l,fx,fy,κ*l]
end

function initialize_beam(node::Clamp,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0     
    m,fx,fy = parameters .* [normfactor_m(beam),normfactor_f(beam),normfactor_f(beam)] #am Balkenelement
    return [m,θ0 + θs,x/l,y/l,fx,fy,κ*l]
end 

function initialize_beam(node::Free,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0     
    Δθ,Δx,Δy = parameters  #am Balkenelement
    return [zero(T),θ0 + θs + Δθ,(x + Δx)./l,(y + Δy)./l,zero(T),zero(T),κ*l]
end 

function initialize_beam(node::ExtForces,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    Δθ,Δx,Δy = parameters
    m = node.mz .* normfactor_m(beam) #am Balkenelement
    fx = node.fx .* normfactor_f(beam) #am Balkenelement
    fy = node.fy .* normfactor_f(beam) #am Balkenelement
    return [m,θ0 + θs + Δθ,(x + Δx)./l,(y + Δy)./l,fx,fy,κ*l]
end 

function initialize_beam(bn::NamedTuple,x,nodes,i::Int)
    pars = x[get_index_pars(bn.Nodes,nodes[1:i])]
    initialize_beam(bn.Nodes[nodes[i]],bn.Beams[i],pars)
end 

scalepos(beam::Beam,y::AbstractArray,::Val{2})  = [y[1] - beam.θe,y[2].*beam.l,y[3].*beam.l]
scalepos(beam::Beam,y::AbstractArray,::Val{1})  = [y[1] - beam.θs,y[2].*beam.l,y[3].*beam.l]

function reduceposat(node::Clamp,beams::NamedTuple,y::AbstractArray{T,3},beamnbrs) where{T}
    res = map(x-> [node.ϕ,node.x,node.y] .- scalepos(beams[x],y[2:4,2,x],Val(2)),beamnbrs[1])
    reducevcat(res)
end

function reduceposat(node::Boundary,beams::NamedTuple,y::AbstractArray{T,3},beamnbrs) where{T}
    bnbrs = vcat(beamnbrs...)
    pos = y[getpositionindices(beamnbrs)]
    scaledpos = map((x,y)->scalepos(beams[x],y,Val(getside(x,beamnbrs))),bnbrs,eachcol(pos))
    res = map(y -> scaledpos[1].- y,scaledpos[2:end])
    # res = mean(y -> scaledpos[1].- y,scaledpos[2:end])
    # res./(norm(res)+eps(T))
    reducevcat(res)
end 

reduceforceat(inits,beam,y) = inits .+  y ./normvector(beam)

function reduceforceat(node::Boundary,Beams::NamedTuple,y::AbstractArray{T,3},beamsidxs) where{T}
    solp = reduce((init,beampos)->reduceforceat(init,Beams[beampos],y[[1,5,6],2,beampos]),beamsidxs[1];init = zeros(T,3))
    solm = reduce((init,beampos)->reduceforceat(init,Beams[beampos],y[[1,5,6],1,beampos]),beamsidxs[2];init = zeros(T,3)) 
    return solp .- solm .+ node[[6,4,5]]
end 

reduceforceat(node::Clamp,Beams::NamedTuple,y::AbstractArray{T,3},beamsidxs) where {T} = Vector{T}()

function residuals!(residuals,str::Structure,y::AbstractArray{T,3},bn) where{T}
    ind = 1
    adj = str.AdjMat
    nodes = findall(x->!isapprox(x,0),LowerTriangular(adj))
    start = 1
    for n in unique(first.(Tuple.(nodes)))
        node = bn.Nodes[n]
        beams = findbeamsatnode(node,n,nodes)
        res = reduceforceat(node,bn.Beams,y,beams)
        if !isempty(res)
            residuals[start:start+2] .= res
            start += length(res)
        end 
        res = reduceposat(node,bn.Beams,y,beams)
        if !isempty(res)
            idxs = range(start,length = length(res))
            residuals[idxs] .= res
            start = idxs[end] + 1
        end  
    end
    residuals 
end 

function (str::Structure)(x::AbstractVector{T},bn::NamedTuple, #x[:,i] = u0 = [m_i,θ_i,x_i<- node_i,y_i <-node_i,fx_i,fy_i,κ = 0]  
    ::Val{false}
    ) where{T}
    nodes = getstartnodes(str.AdjMat)

    function prob_func(prob,i,repeat) 
        u0 = initialize_beam(bn,x,nodes,i)
        remake(prob;u0 = u0)
    end 

    ensprob  =  EnsembleProblem(prob;
                prob_func = prob_func,
                )

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-12,abstol = 1e-12,
                save_start = true,save_on = false,save_end = true,
                sensealg=str.SensAlg,
                trajectories = length(bn.Beams)
                )
    Array(sol)
end

function (str::Structure)(x::AbstractVector{T},bn::NamedTuple, #x[:,i] = u0 = [m_i,θ_i,x_i<- node_i,y_i <-node_i,fx_i,fy_i,κ = 0]  
    ::Val{true}
    ) where{T}
    nodes = getstartnodes(str.AdjMat)

    function prob_func(prob,i,repeat) 
        pars = x[get_index_pars(bn.Nodes,nodes[1:i])]
        u0 = initialize_beam(bn.Nodes[nodes[i]],bn.Beams[i],pars)
        remake(prob;u0 = u0)
    end 

    ensprob  =  EnsembleProblem(prob;
                prob_func = prob_func)

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-12,abstol = 1e-12,
                sensealg=str.SensAlg,
                trajectories = length(bn.Beams);str.kwargs...
                )
end
(str::Structure)(x,bn,plt::Bool = false) = str(x,bn,Val(plt))


function (str::Structure)(residuals::T,values::T,bn::NamedTuple) where{T} #new loss
    sols = str(values,bn)
    residuals!(residuals,str,sols,bn)
end 

function (str::Structure)(values::AbstractVector,forplt::Bool = false)
    inits = initialize(str,values)
    str(inits,forplt)
end 
function (str::Structure)(f::Union{typeof(zeros),typeof(rand)},forplt::Bool = false)
    inits = initialize(str,f(Float64,str))
    str(inits,forplt)
end 



function change_beam(str::Structure,beam::Int;kwargs...)
    for (field,value) in kwargs
        str = Setfield.@set str.Beams[beam].$field = value 
    end 
    str 
end  

function getinitials(str::Structure)
    nodes = str.AdjMat.Nodes
    nodes2beams = str.AdjMat.Nodes2Beams
    beams = str.AdjMat.Beams2Nodes
    # branches = findall(x->isa(x,Branch),nodes)
    initsize = 0
    for (beam,(atstart,atend)) in beams
        initsize += isa(nodes[atstart],Branch) ? 6 : 3
    end 
    return initsize
end 

function Base.zeros(::Type{T},str::Structure) where{T}
    zeros(T,getinitials(str))
end 

function Base.zeros(::Type{T},strs::Vector{Structure}) where{T}
    inits = Vector{Vector{T}}(undef,length(strs))
    for ind in eachindex(inits)
        inits[ind] = zeros(T,getinitials(strs[ind]))
    end 
    return inits
end 
function Random.rand(::Type{T},str::Structure) where{T}
    rand(T,getinitials(str))
end 
function Random.rand(::Type{T},strs::Vector{Structure}) where{T}
    inits = Vector{Vector{T}}(undef,length(strs))
    for ind in eachindex(inits)
        inits[ind] = rand(T,getinitials(strs[ind]))
    end 
    return inits
end 
