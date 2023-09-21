
struct Structure{A<:Connections,B,So,Se,KW,BE<:Dict{Int,Beam{T}} where{T}} 
    AdjMat::A
    Branches::B
    Solver::So
    SensAlg::Se
    kwargs::KW
    Beams::BE
    function Structure(adj::A,beams::BE,solver::So,sensalg::Se,kwargs::KW) where{A,BE,So,Se,KW}
        branches = get_branches(adj)
        new{A,typeof(branches),So,Se,KW,BE}(adj,branches,solver,sensalg,kwargs,beams)
    end 
end

function Structure(adj::Connections,beams::Beam ...;solver = Tsit5(),sensalg = ForwardDiffSensitivity(),kwargs...) 
    beams_ = Dict(x=>y for (x,y) in enumerate(beams))
    Structure(adj,beams_,solver,sensalg,kwargs)
end 

function Base.show(io::IO,str::Structure)
    return println(io, "Structure with $(length(str.AdjMat.Nodes2Beams)) Nodes and $(length(str.AdjMat.Beams2Nodes)) Beam(s).")
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


function get_parameters(::Branch,parameters::Matrix{T},ind1::Int,ind2) where{T}
    (view(parameters,:,ind1),view(parameters,:,ind2),ind1+1)
end

function get_parameters(::Boundary,parameters::Matrix{T},ind1::Int) where{T}
    (view(parameters,:,ind1),zeros(T,3))
end

function get_parameters(::Free,parameters::Matrix{T},ind1::Int,ind2) where{T}
    (view(parameters,:,ind1),view(parameters,:,ind2),ind1+1)
end

normfactor_m(b::Beam) = 12* b.l/(b.E*b.w*b.h^3) 
normfactor_f(b::Beam) = b.l * normfactor_m(b) 

function initialize_boundary(node::Boundary,beam::Beam,parameters::AbstractVector{T}) where {T}
    x,y,θ0, = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    m,fx,fy,Δθ,Δx,Δy = parameters
    m *= normfactor_m(beam) #am Balkenelement
    fx *= normfactor_f(beam) #am Balkenelement
    fy *= normfactor_f(beam) #am Balkenelement
    return [m,θ0 + θs + Δθ,(x + Δx)./l,(y + Δy)./l,fx,fy,κ*l]
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

function initialize(str::Structure,parameters::Vector{T}) where {T}
    con = str.AdjMat
    # @assert 3*(length(con.Beams) + length(str.Branches)) == length(parameters)
    mat = Matrix{T}(undef,7,length(con.Beams2Nodes))
    ind = 1 
    for (beam,(st,en)) in con.Beams2Nodes 
        add = isa(con.Nodes[st],Branch) ? 6 : 3
        mat[:,beam] .= initialize_boundary(con.Nodes[st],str.Beams[beam],parameters[ind:ind + add - 1])
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
    res_force .=  getsign(point) .* y[point][[1,5,6],beam] ./ [normfactor_m(beams[beam]),normfactor_f(beams[beam]),normfactor_f(beams[beam])] .- node[[6,4,5]] #   
    ŷp = y[point][2:4,beam].*[1,l,l]
    ŷp[1] -= θ # v+ θn 
    ind += 3
    for (beam,point) in beams_end[2:end]
        l = beams[beam].l
        θ = beams[beam][5+point]
        res_pos = @view res[ind:ind+2]
        res_pos .= ŷp .- y[point][2:4,beam].*[1,l,l]
        res_pos[1] += θ #+ θn
        res_force .+= getsign(point) .* y[point][[1,5,6],beam] ./ [normfactor_m(beams[beam]),normfactor_f(beams[beam]),normfactor_f(beams[beam])]
        ind += 3 
    end
    ind
end 

function residuals!(res::AbstractVector{T},node::ExtForces,beams,beams_end,y,ind) where{T}
    #Forces are  important
    for (beam,point) in beams_end
        res = @view res[ind:ind+2]
        res .= node[[6,4,5]]  .+ getsign(point).* y[point][[1,5,6],beam] ./ [normfactor_m(beam),normfactor_f(beam),normfactor_f(beam)]

        ind += 3
    end 
    ind
end 

function residuals!(res::AbstractVector{T},node::Free,beams,beams_end,y,ind) where{T}
    for (beam,point) in beams_end
        res_ = @view res[ind:ind+2]
        res_ .= y[point][[1,5,6],beam] ./ [normfactor_m(beams[beam]),normfactor_f(beams[beam]),normfactor_f(beams[beam])]
        ind += 3
    end 
    ind
end 
function residuals!(res::AbstractVector{T},node::Boundary,beams,beams_end,y,ind) where{T}
    #Forces and positions are important
    for (beam,point) in beams_end
        res_ = @view res[ind:ind+2]
        res_ .= node[[6,4,5]] .+ getsign(point) .* y[point][[1,5,6],beam] ./ [normfactor_m(beams[beam]),normfactor_f(beams[beam]),normfactor_f(beams[beam])]

        ind += 3
    end 
    ind
end 

function residuals!(res::AbstractVector{T},node::Clamp,beams,beams_end,y,ind) where{T}
    #Positions are important
    for (beam,point) in beams_end
            res_ = @view res[ind:ind+2]
            res_ .= node[[3,1,2]] .- y[point][2:4,beam] .* [1,beams[beam].l,beams[beam].l]
            res_[1] += beams[beam].θe            
            ind += 3
    end 
    ind
end 

function residuals!(res::AbstractVector{T},node::CompliantClamp,beams,beams_end,y,ind) where{T}
    for (beam,point) in beams_end
        res_ = @view res[ind:ind+2]
        m = node.c * (y[point] - node.ϕ)
        res_[1] .= m .- y[point][1,beam] * normfactor_m(beams[beam])
        res_[2:3] .= node[2:3] .- y[point][3:4,beam] .* [beams[beam].l,beams[beam].l]
       
        ind += 3
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

function (str::Structure)(residuals::T,values::T) where{T} #new loss
    inits = initialize(str,values)
    sols = str(inits)
    residuals!(residuals,str,sols)
end 
function (str::Structure)(values::AbstractVector,forplt::Bool = false)
    inits = initialize(str,values)
    str(inits,forplt)
end 