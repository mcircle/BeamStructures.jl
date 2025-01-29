struct GroundStructure{So,Se}
    Solver::So
    SensAlg::Se
    kwargs
    function GroundStructure(solver::So,sensalg::Se,kwargs) where{So,Se}
        new{So,Se}(solver,sensalg,kwargs)
    end 
end 

function GroundStructure(;solver = Tsit5(),sensalg = ForwardDiffSensitivity(),kwargs...)
    GroundStructure(solver,sensalg,kwargs)
end 

function (str::GroundStructure)(x::AbstractMatrix{T},bn::NamedTuple,adj,::Val{false}) where{T}
    x_ = reshape(x,3,:)

    nodepos = getstartnodes(adj)
    cbeams = reduce(+,1:size(adj,1)-1)
    anz,nodes_ = changestartnodes(bn.Nodes,x_)
    function prob_func(prob,i,repeat) 
        u0 = initialize_beam(bn.Beams,nodes_,x_[:,anz+1:end],nodepos,i)
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
                trajectories = cbeams
                )
    st = admittance_matrix(sol,adj,str)
    Array(sol),(;Beams = bn.Beams,Nodes = nodes_),st
end

function (str::GroundStructure)(x::AbstractMatrix{T},bn::NamedTuple,adj,::Val{true}) where{T}
    cbeams = reduce(+,1:size(adj,1)-1)
    nodepos = getstartnodes(adj)

    anz,nodes_ = changestartnodes(bn.Nodes,x)

    function prob_func(prob,i,repeat) 
        u0 = initialize_beam(bn.Beams,nodes_,x[:,anz+1:end],nodepos,i)
        remake(prob;u0 = u0)
    end 

    ensprob  =  EnsembleProblem(prob;
                prob_func = prob_func)

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-6,abstol = 1e-6,
                sensealg=str.SensAlg,
                trajectories = cbeams;str.kwargs...)
end

function reduction_func_admittance!(u,data,I,adj)
    for ((i,j),d0) in zip(I,data)
        d = d0 .*adj[i,j]
        u[i,j,:] .-= d 
        u[j,i,:] .-= d 
        u[i,i,:] .+= d
        # u[j,j,:] .+= d 
    end 
    u,false
end 

function output_func_admittance(sol,i)
    # @show sol[end]
    (sol[end],false)
end

function getidxs(sz)
    len = reduce(+,1:sz-1)
    idxs = Vector{NTuple{2,Int}}(undef,len)
    i = 1
    j = 1
    for id in eachindex(idxs)
        j += i == sz ? 1 : 0
        i += i == sz ? -(sz-1-j) : 1   
        idxs[id] = (i,j)      
    end
    idxs
end 

function admittance_matrix(sol::EnsembleSolution{T,N,S},adj,str) where{T,N,S}
    
    idxs = getidxs(size(adj,1))

    function prob_func(prob,i,repeat) 
        u0 = zeros(T,7)
        u0[2:4] .= one(T)
        remake(prob;u0 = u0,p = sol[i])
    end
    #dadj Gradient der Steifigkeitsmatrix berechnen  
    ensprob2  =  EnsembleProblem(vjpprob;prob_func = prob_func,
                                        output_func = (sol,i) ->  output_func_admittance(sol,i),
                                        reduction = (u,data,I) -> reduction_func_admittance!(u,data,idxs[I],adj),
                                        u_init = zeros(T,size(adj)...,3)
                                        )
    sol = solve(ensprob2,str.Solver,save_on = false,save_start = false,save_end = true,save_idxs = [1,5,6],
                                        EnsembleThreads(),
                                        reltol = 1e-6,abstol = 1e-6,
                                        trajectories = length(sol)
                                        )

    return sol
end 

(str::GroundStructure)(x::AbstractMatrix,bn,adj,plt::Val) = str(x,bn,adj,plt)
(str::GroundStructure)(x::AbstractVector,bn,adj,plt::Val) = str(reshape(x,3,:),bn,adj,plt)

function (str::GroundStructure)(x,bn,adj,plt::Bool = false)
    cbeams = reduce(+,1:size(adj,1)-1)
    @assert cbeams == length(bn.Beams) "Missing Beams! Need $(cbeams) and got only $(length(bn.Beams))"
    @assert size(adj,1) == length(bn.Nodes) "Missing Nodes! Need $(size(adj,1)) and got only $(length(bn.Nodes))"
    ci = getindices(size(adj,1))
    str(x,bn,adj,Val(plt))
end 

function reduceposat(node::Boundary,beams::NamedTuple,y::AbstractArray{T,3},factors::AbstractVector{T},beamnbrs) where{T}
    res = map((ind)-> factors[ind] .* ([node.ϕ,node.x,node.y] .- scalepos(beams[ind],y[2:4,2,ind],Val(2))),beamnbrs[1])#,init = zeros(T,3))./length(beamnbrs[1])
    reduce(hcat,res)
end
#at Clamp force is equal to surounding
reduceforceat(node::Clamp,Beams::NamedTuple,y::AbstractArray{T,3},facs::AbstractVector{T},beamsidxs) where {T} = Vector{T}()

function reduceforceat(node::Boundary,Beams::NamedTuple,y::AbstractArray{T,3},factors::AbstractVector{T},beamsidxs) where{T}
    solp = reduce((init,beampos)->init .+ factors[beampos] * scaleforce(Beams[beampos],y[[1,5,6],2,beampos]),beamsidxs[1];init = zeros(T,3))
    solm = reduce((init,beampos)->init .+ factors[beampos] * scaleforce(Beams[beampos],y[[1,5,6],1,beampos]),beamsidxs[2];init = zeros(T,3)) 
    return solp .- solm .+ node[[6,4,5]]
end 

function residuals!(residuals::Matrix,adj::AbstractMatrix{T},y::AbstractArray{T,3},bn) where{T}
    
    ids = getindices(size(adj,1))
    # idcs = LinearIndices(residuals)
    nodes = findall(x->!isapprox(x,0),LowerTriangular(adj))
    start = 1
    for n in unique(first.(Tuple.(nodes)))
        node = bn.Nodes[n]
        beams = findbeamsatnode(node,n,nodes)
        res = reduceforceat(node,bn.Beams,y,adj[ids],beams)
        if !isempty(res)
            residuals[:,start] .= res
            start += 1
        end 
        res = reduceposat(node,bn.Beams,y,adj[ids],beams)
        if !isempty(res)
            idxs = start:start + size(res,2) - 1
            residuals[:,idxs] .= res
            start = idxs[end] + 1
        end  
    end
    residuals 
end 
