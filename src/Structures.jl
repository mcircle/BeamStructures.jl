
struct Structure{A<:Connections,B,So,Se,KW}
    AdjMat::A
    Branches::B
    Solver::So
    SensAlg::Se
    kwargs::KW
end

function Structure(adj::Connections;solver = Tsit5(),sensalg = ForwardDiffSensitivity(),kwargs...)
    branches = get_branches(adj)
    Structure(adj,branches,solver,sensalg,kwargs)
end 

function Base.show(io::IO,str::Structure)
    return println(io, "Structure with $(length(str.AdjMat.Nodes)) Nodes and $(length(str.AdjMat.Beams)) Beams.")
end 

out4plt(sol, i) = (sol, false)
out4loss(sol,i) = ((sol[1],sol[end]),false)
reduction4plt(u,data,I) = (append!(u, data), false)
function reduction4loss(u,data,batch) #loss
    # @show size(u)
    # @show size(data)
    for b in batch
        view(u[1],:,b) .=  data[b][1]
        view(u[2],:,b) .=  data[b][2]
    end 
    (u, false)
end
function (str::Structure)(x::Matrix, #x[:,i] = u0 = [m_i,θ_i,x_i<- node_i,y_i <-node_i,fx_i,fy_i,κ = 0]  
        forplt::Bool = false
        )
    
    prob_func(prob,i,repeat) = remake(prob,u0 = x[:,i])

    output = forplt ? out4plt : out4loss
    reduct = forplt ? reduction4plt : reduction4loss
    ensprob  =  EnsembleProblem(prob;
                output_func = output,
                prob_func = prob_func,
                reduction = reduct,
                u_init = forplt ? [] : [similar(x),similar(x)]
                )

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-12,abstol = 1e-12,
                save_start = true,save_on = forplt,save_end = true,
                sensealg=str.SensAlg,
                trajectories = length(str.AdjMat.Beams)
                )
end 

function initialize_boundary(::Boundary,node,forces::AbstractVector{T},positions::AbstractVector{T},nodefeatures,beamfeatures = zeros(T,4)) where {T} 
    θ0,x,y = T.(nodefeatures[1:3,node])#Grundposition des Knotens 
    l,θs,κ = beamfeatures[1:3] 
    m,fx,fy = forces .* [l,l^2,l^2] #am Balkenelement
    return [m,θ0 + θs,x/l,y/l,fx,fy,κ*l]
end

function initialize_boundary(::Branch,node,forces::AbstractVector{T},positions::AbstractVector{T},nodefeatures,beamfeatures = zeros(T,4)) where {T} 
    θ0,x,y = T.(nodefeatures[1:3,node]) #Grundposition des Knotens 
    l,θs,κ = beamfeatures[1:3]
    Δθ,Δx,Δy = positions # Δx,Δy,Δθ von Verschiebung des Nodes
    m,fx,fy = forces .* [l,l^2,l^2] #am Balkenelement
    return [m,θ0 + θs + Δθ,(x + Δx)./l,(y + Δy)./l,fx,fy,κ*l]
end 

function get_parameters(::Branch,parameters::Matrix{T},ind1::Int,ind2) where{T}
    (view(parameters,:,ind1),view(parameters,:,ind2),ind1+1)
end

function get_parameters(::Boundary,parameters::Matrix{T},ind1::Int,ind2) where{T}
    (view(parameters,:,ind1),zeros(T,3),ind1+1)
end

function initialize(str::Structure,parameters::Vector{T},features::Matrix{S},beamfeatures) where {T,S}
    con = str.AdjMat
    # @assert 3*(length(con.Beams) + length(str.Branches)) == length(parameters)
    mat = Matrix{T}(undef,7,length(con.Beams))
    pars = reshape(parameters,3,:)
    ind_beams = 1 + length(str.Branches)
    for (beam,(st,en)) in con.Beams
        ind_branch = st ∈ keys(str.Branches) ? str.Branches[st] : 0
        forces,positions,ind_beams = get_parameters(con.NodeType[st],pars,ind_beams,ind_branch)           
        mat[:,beam] .= initialize_boundary(con.NodeType[st],st,forces,positions,features,beamfeatures[:,beam])
    end 
    return mat
end

getsign(x) = ifelse(x == 1, -1, 1)

function residuals!(res::AbstractVector{T},::Branch,beams,y,features,beamfeatures,ind) where{T}
    beam,point = beams[1]
    l,θ = beamfeatures[[1,2*point],beam]
    # θn = features[1] #hmmm 
    res_force = @view res[ind:ind+2]
    res_force .=  getsign(point) .* y[point][[1,5,6],beam] ./ [l,l^2,l^2] .- features[4:6] #   
    ŷp = y[point][2:4,beam].*[1,l,l]
    ŷp[1] -= θ # v+ θn 
    ind += 3
    for (beam,point) in beams[2:end]
        l,θ = beamfeatures[[1,point*2],beam]
        res_pos = @view res[ind:ind+2]
        res_pos .= ŷp .- y[point][2:4,beam].*[1,l,l]
        res_pos[1] += θ #+ θn
        res_force .+= getsign(point) .* y[point][[1,5,6],beam] ./ [l,l^2,l^2]
        ind += 3 
    end
    ind
end 

function residuals!(res::AbstractVector{T},::Free,beams,y,features,beamfeatures,ind) where{T}
    #Forces are  important
    for (beam,point) in beams
        l = beamfeatures[1,beam]
        res = @view res[ind:ind+2]
        res .= features[4:6]  .+ getsign(point).* y[point][[1,5,6],beam] ./ [l,l^2 ,l^2 ]
        # @show res
        ind += 3
    end 
    ind
end 
function residuals!(res::AbstractVector{T},::Boundary,beams,y,features,beamfeatures,ind) where{T}
    #Forces and positions are important
    for (beam,point) in beams
        l,θe =beamfeatures[[1,4],beam]
        res = @view res[ind:ind+2]
        res .= features[4:6] .+ getsign(point) .* y[point][[1,5,6],beam] ./ [l,l^2,l^2]
        # @show res
        ind += 3
    end 
    ind
end 
function residuals!(res::AbstractVector{T},::Clamp,beams,y,features,beamfeatures,ind) where{T}
    #Positions are important
    for (beam,point) in beams
            l,θe =beamfeatures[[1,4],beam]
            res_ = @view res[ind:ind+2]
            res_ .= features[1:3] .- y[point][2:4,beam] .* [1,l,l]
            res_[1] += θe
            # @show res
            ind += 3
    end 
    ind
end 

function residuals!( residuals,str::Structure,y::EnsembleSolution{T,3,Vector{Matrix{T}}},
                     features::Matrix,beamfeatures::Matrix = zeros(T,4,3),width = 50e-3,height = 1e-3,E = 2e+12) where{T}
    ind = 1
    con = str.AdjMat
    for (node,beams) in con.Nodes
        t = con.NodeType[node]
        ind = residuals!(residuals,t,beams,y,features[:,node],beamfeatures,ind) 
        # @show residuals
    end 
end 

function loss!(residuals,structure,parameters)
    inits = initialize(str,parameters,features,beamfeatures) 
    odesol = structure(inits,false)   
    residuals!(residuals,str,odesol,features,beamfeatures)
end
function loss!(residuals,structure,parameters,features,beamfeatures)
    inits = initialize(str,parameters,features,beamfeatures) 
    odesol = structure(inits,false)   
    residuals!(residuals,str,odesol,features,beamfeatures)
end

