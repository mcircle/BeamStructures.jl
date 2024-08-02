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

function (str::GroundStructure)(x::AbstractVector{T},bn::NamedTuple,adj,::Val{false}) where{T}
    x_ = reshape(x,3,:)

    nodepos = nodepos = getstartnodes(adj)
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
    Array(sol),(;Beams = bn.Beams,Nodes = nodes_)
end

function (str::GroundStructure)(x::AbstractVector{T},bn::NamedTuple,adj,::Val{true}) where{T}
    x_ = reshape(x,3,:)
    nodepos = nodepos = getstartnodes(adj)
    cbeams = reduce(+,1:size(adj,1)-1)

    anz,nodes_ = changestartnodes(bn.Nodes,x_)

    function prob_func(prob,i,repeat) 
        u0 = initialize_beam(bn.Beams,nodes_,x_[:,anz+1:end],nodepos,i)
        remake(prob;u0 = u0)
    end 

    ensprob  =  EnsembleProblem(prob;
                prob_func = prob_func)

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-12,abstol = 1e-12,
                sensealg=str.SensAlg,
                trajectories = cbeams;str.kwargs...)
end

function (str::GroundStructure)(x,bn,adj,plt::Bool = false)
    cbeams = reduce(+,1:size(adj,1)-1)
    @assert cbeams == length(bn.Beams) "Missing Beams! Need $(cbeams) and got only $(length(bn.Beams))"
    @assert size(adj,1) == length(bn.Nodes) "Missing Nodes! Need $(size(adj,1)) and got only $(length(bn.Nodes))"
    ci = getindices(size(adj,1))
    str(x,bn,adj,Val(plt))
end 

function reduceposat(node::Boundary,beams::NamedTuple,y::AbstractArray{T,3},factors::AbstractVector{T},beamnbrs) where{T}
    res = reduce((x,ind)-> x .+ factors[ind] .* abs2.([node.ϕ,node.x,node.y] .- scalepos(beams[ind],y[2:4,2,ind],Val(2))),beamnbrs[1],init = zeros(T,3))
    # res./(norm(res) + eps(T))
end

function reduceforceat(node::Boundary,Beams::NamedTuple,y::AbstractArray{T,3},factors::AbstractVector{T},beamsidxs) where{T}
    solp = reduce((init,beampos)->init .+ factors[beampos] * scaleforce(Beams[beampos],y[[1,5,6],2,beampos]),beamsidxs[1];init = zeros(T,3))
    solm = reduce((init,beampos)->init .+ factors[beampos] * scaleforce(Beams[beampos],y[[1,5,6],1,beampos]),beamsidxs[2];init = zeros(T,3)) 
    return solp .- solm
end 

function residuals!(residuals::Matrix,adj::AbstractMatrix{T},y::AbstractArray{T,3},bn) where{T}
    cis = getindices(size(adj,1))
    idx = 1
    for n in axes(adj,1)
        node = bn.Nodes[n]
        beams = findbeamsatnode(node,n,cis)
        if any(!isempty,beams)
            residuals[:,idx] .= reduceforceat(node,bn.Beams,y,adj[cis],beams)
            residuals[:,idx+1] .= reduceposat(node,bn.Beams,y,adj[cis],beams)
            idx += 2 
        end 
    end
    residuals 
end 
