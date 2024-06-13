
struct Structure{A<:AbstractMatrix,Nn,So,Se,KW}
    AdjMat::A
    nodes::Nn
    Solver::So
    SensAlg::Se
    kwargs::KW
    function Structure(adj::A,nodes::Nn,solver::So,sensalg::Se,kwargs::KW) where{A,Nn,So,Se,KW}
        nodes_ = gettype.(nodes)
        new{A,NTuple{length(nodes),Boundary},So,Se,KW}(adj,nodes_,solver,sensalg,kwargs)
    end 
end

function Structure(adj::AbstractMatrix{T},nodes::Vararg{Boundary,N};solver = Tsit5(),sensalg = ForwardDiffSensitivity(),kwargs...) where{T,N}
    pos = vec(any(==(one(T)),adj,dims = 1))
    @assert sum(pos) == N  "Number of Nodes ($length(nodes)) is does not match the adjency matrix!"
    Structure(adj,nodes,solver,sensalg,kwargs)
end 

function Structure(adj::Connections,nodes::Vararg{Boundary,N};solver = Tsit5(),sensalg = ForwardDiffSensitivity(),kwargs...) where{N}
    Structure(adj,nodes,solver,sensalg,kwargs)
end 

function Base.show(io::IO,str::Structure)
    return println(io, "Structure with $(str.nodes) Nodes and $(length(str.AdjMat)) Beam(s).")
end 

#ODE outcome
out4plt(sol, i) = (sol, false)
out4loss(sol,i) = ((sol[1],sol[end]),false)
reduction4plt(u,data,I) = (append!(u, data), false)
function reduction4loss(u,data,batch) 
    
    for b in batch
        view(u[1],:,b) .=  data[b][1]
        view(u[2],:,b) .=  data[b][2]
    end 
    (u, false)
end

normfactor_m(b::Beam) = 12* b.l/(b.E*b.w*b.h^3) 
normfactor_f(b::Beam) = b.l * normfactor_m(b) 
normvector(b::Beam) = [normfactor_m(b),normfactor_f(b),normfactor_f(b)]

function initialize_boundary(node::Boundary,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0, = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    m,fx,fy,Δθ,Δx,Δy = parameters
    m *= normfactor_m(beam) #am Balkenelement
    fx *= normfactor_f(beam) #am Balkenelement
    fy *= normfactor_f(beam) #am Balkenelement
    return [m,θ0 + θs + Δθ,(x + Δx)./l,(y + Δy)./l,fx,fy,κ*l]
end 

function initialize_boundary(node::CompliantClamp,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0, = node.x,node.y,node.ϕ
    l,κ = beam.l,beam.κ0
    m = node.c * (beam.θs .- parameters[1]) .* normfactor_m(beam)
    fx = parameters[2] * normfactor_f(beam) #am Balkenelement
    fy = parameters[2] * normfactor_f(beam) #am Balkenelement
    return [m,θ0 + beam.θs + parameters[1],x./l,y./l,fx,fy,κ*l]
end

function initialize_boundary(node::Clamp,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0     
    m,fx,fy = parameters .* [normfactor_m(beam),normfactor_f(beam),normfactor_f(beam)] #am Balkenelement
    return [m,θ0 + θs,x/l,y/l,fx,fy,κ*l]
end 

function initialize_boundary(node::Free,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0     
    Δθ,Δx,Δy = parameters  #am Balkenelement
    return [zero(T),θ0 + θs + Δθ,(x + Δx)./l,(y + Δy)./l,zero(T),zero(T),κ*l]
end 

function initialize_boundary(node::ExtForces,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0 = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    Δθ,Δx,Δy = parameters
    m = node.mz .* normfactor_m(beam) #am Balkenelement
    fx = node.fx .* normfactor_f(beam) #am Balkenelement
    fy = node.fy .* normfactor_f(beam) #am Balkenelement
    return [m,θ0 + θs + Δθ,(x + Δx)./l,(y + Δy)./l,fx,fy,κ*l]
end 

function initialize(str::Structure,parameters::Vector{T},bn::NamedTuple) where {T}
    adj = str.AdjMat
    # @assert 3*(length(con.Beams) + length(str.Branches)) == length(parameters)
    mat = Matrix{T}(undef,7,length(bn.Beams))
    ind = 1 
    for (beam,ci) in enumerate(findall(==(1),LowerTriangular(adj))) 
        st = ci[1]
        add = isa(str.nodes[st],Branch) ? 6 : 3
        mat[:,beam] .= initialize_boundary(str.nodes[st],bn.Beams[beam],parameters[ind:ind + add - 1])
        ind += add
    end 
    return mat
end

getsign(x) = ifelse(x == 1, -1, 1)

function residuals!(res::AbstractVector{T},node::Branch,beams,beams_end,y,ind) where{T}
    beam,point = beams_end[1]
    l = beams[beam].l
    θ = beams[beam][5+point]
    res_force = @view res[ind:ind+2]
    res_force .=  getsign(point) .* y[point][[1,5,6],beam] ./ normvector(beams[beam]).- node[[6,4,5]] #   
    ŷp = y[point][2:4,beam].*[1,l,l]
    ŷp[1] -= θ # v+ θn 
    ind += 3
    for (beam,point) in beams_end[2:end]
        @inbounds l =  beams[beam].l
        @inbounds θ =  beams[beam][5+point]
        res_pos = @view res[ind:ind+2]
        @inbounds res_pos .= ŷp .-  y[point][2:4,beam].*[1,l,l]
        @inbounds res_pos[1] += θ
        @inbounds res_force .+= getsign(point) .* y[point][[1,5,6],beam] ./ normvector(beams[beam])
        ind += 3 
    end
    ind
end 

function residuals!(res::AbstractVector{T},node::ExtForces,beams,beams_end,y,ind) where{T}
    #Forces are  important
    for (beam,point) in beams_end
        @inbounds res = @view res[ind:ind+2]
        @inbounds res .= node[[6,4,5]]  .+ getsign(point).* y[point][[1,5,6],beam] ./ normvector(beams[beam])
        ind += 3
    end 
    ind
end 

function residuals!(res::AbstractVector{T},node::CompliantClamp,beams,beams_end,y,ind) where{T}
    #Forces are  important
    for (beam,point) in beams_end
        @inbounds res = @view res[ind:ind+2] 
        @inbounds res[1] = node.c * (beam[4 .+ point] .- y[point][2,beam]) .- getsign(point).* y[point][1,beam] ./ normfactor_m(beams[beam])
        @inbounds res[2] = node.x .- y[point][3,beam] .* beams[beam].l
        @inbounds res[3] = node.y .- y[point][4,beam] .* beams[beam].l
        ind += 3
    end 
    ind
end 

function residuals!(res::AbstractVector{T},node::Free,beams,beams_end,y,ind) where{T}
    for (beam,point) in beams_end
        @inbounds res_ = @view res[ind:ind+2]
        @inbounds res_ .= y[point][[1,5,6],beam] ./ normvector(beams[beam])
        ind += 3
    end 
    ind
end 

function residuals!(res::AbstractVector{T},node::Boundary,beams,beams_end,y,ind) where{T}
    #Forces and positions are important
    for (beam,point) in beams_end
        @inbounds res_ = @view res[ind:ind+2]
        @inbounds res_ .= node[[6,4,5]] .+ getsign(point) .* y[point][[1,5,6],beam] ./ normvector(beams[beam])
        ind += 3
    end 
    ind
end 

function residuals!(res::AbstractVector{T},node::Clamp,beams,beams_end,y,ind) where{T}
    #Positions are important
    beamsiter = filter(x->last(x) .== 2,beams_end)
    isempty(beamsiter) && return ind
    beamnodesiter = iterate(beamsiter,1) 
    while ~isnothing(beamnodesiter)
        ((beam,point),iter) = beamnodesiter
        @inbounds res_ = @view res[ind:ind+2]
        @inbounds  res_ .= node[[3,1,2]] .- y[point][2:4,beam] .* [1,beams[beam].l,beams[beam].l]
        @inbounds res_[1] += beams[beam].θe            
        ind += 3
        beamnodesiter = iterate(beamsiter,iter)
    end 
    ind
end 

function residuals!(residuals,str::Structure,y::T) where{T<:EnsembleSolution}
    ind = 1
    con = str.AdjMat
    for (node,beams) in con.Nodes2Beams
        ind = residuals!(residuals,con.Nodes[node],str.Beams,beams,y,ind) 
    end 
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
                trajectories = length(str.Beams)
                )
end

function (str::Structure)(residuals::T,values::T;nodes::NTuple{Nn,Boundary},beams::NTuple{Nb,Beam}) where{T,Nn,Nb} #new loss
    inits = initialize(str,values)
    sols = str(inits)
    residuals!(residuals,str,sols)
end 

function (str::Structure)(values::AbstractVector,forplt::Bool = false)
    inits = initialize(str,values)
    str(inits,forplt)
end 
function (str::Structure)(f::Union{typeof(zeros),typeof(rand)},forplt::Bool = false)
    inits = initialize(str,f(Float64,str))
    str(inits,forplt)
end 

function change_node(str::Structure,node::Int;kwargs...)
    for (field,value) in kwargs
        str = Setfield.@set str.AdjMat.Nodes[node].$field = value 
    end 
    str 
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
