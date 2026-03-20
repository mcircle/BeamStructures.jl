CRC.has_mutable_tangent(::BeamElement) = true
CRC.has_mutable_tangent(::Boundary) = true

function CRC.rrule(::typeof(reducevcat),res)
    sol = reducevcat(res)
    function vcat_back(ȳ)
        if isempty(res)
            grads = NoTangent()
        else
            grads = @thunk(collect(eachcol(reshape(ȳ,3,:))))
        end 
        return NoTangent(),grads
    end 

    return sol,vcat_back
end 

function CRC.rrule(::Type{Beam},l,h,w,κ0,E,θs,θe)
    function back_beam(ȳ)
        return NoTangent(),ȳ.l,ȳ.h,ȳ.w,ȳ.κ0,ȳ.E,ȳ.θs,ȳ.θe
    end 
    return Beam(l,h,w,κ0,E,θs,θe),back_beam
end 

function CRC.rrule(::Type{B},x,y,ϕ,fx,fy,mz) where{T<:Real,B<:Boundary{T}}
    function back_boundary(ȳ)
        return NoTangent(),ȳ.x,ȳ.y,ȳ.ϕ,ȳ.fx,ȳ.fy,ȳ.mz
    end 
    return B(x,y,ϕ,fx,fy,mz),back_boundary
end

function CRC.rrule(::typeof(+),b::B,nt::NamedTuple) where{T<:Real,B<:Boundary{T}}
    function back_add(ȳ)
        return NoTangent(), ȳ, ZeroTangent()
    end 
    return b + nt,back_add
end

function scalepos_back(ȳ,y,beam::B) where{T,B<:BeamElement{T}}
    dy = ȳ .* [1,beam.l,beam.l]
    dbeam = Tangent{B}(;l = sum(y[2:3] .* ȳ[2:3] ) ,θe = -ȳ[1])
    return NoTangent(),dbeam, dy,NoTangent() 
end 

function CRC.rrule(::typeof(scalepos),beam::B,y,side::Val{2}) where{T,B<:BeamElement{T}}
    sol = scalepos(beam,y,side)
    sc_back(ȳ) = scalepos_back(ȳ,y,beam)
    return sol,sc_back
end 

function CRC.rrule(::typeof(reduceposat),node::LinearSlider{T},beams,y,beamnbrs) where{T}
    sol = reduceposat(node,beams,y,beamnbrs)
    function reduceposat_back(ȳ)
        rȳ = reshape(ȳ,3,:)
        
        ∂beams= Tangent{typeof(beams)}(;ntuple(x->keys(beams)[x]=>ZeroTangent(),length(beams))...)
        ∂y = zero(y)
        for (n,b) in enumerate(beamnbrs[1])
            _,dbeam,dy,_ = scalepos_back(rȳ[:,n],y[:,n],beams[b])    
            ∂beams -= Tangent{typeof(beams)}(;Symbol(:Beam_,b) => dbeam)
            ∂y[2:4,2,b] .-= dy
        end
        dϕ =  sum(rȳ[1,:])
        dx = sum(rȳ[2,:]) 
        dy = sum(rȳ[3,:]) 
        ds = dx * cos(dϕ) + dy * sin(dϕ)
        dir = node.s*(dx * sum(x->-sin(x),rȳ[1,:]) + dy * sum(x->cos(x),rȳ[1,:]))
        dnode = Tangent{LinearSlider{T}}(;x = dx,y = dy,dir = dir,ϕ = dϕ,s = ds)
        return NoTangent(),dnode,∂beams,∂y,NoTangent()
    end 
    return sol, reduceposat_back
end 

function reduceposat_back(ȳ,∂y,y,
                          beamnbrs::AbstractVector{Int},node::Bo,beams::BT,
                          ∂beams = CRC.zero_tangent(beams),dnode = CRC.zero_tangent(node)) where{beamnames,Tu,Bo,BT<:NamedTuple{beamnames,Tu}}
    
    for (n,b) in enumerate(beamnbrs)
        _,dbeam,dy,_ = scalepos_back(ȳ[:,n],y[2:4,2,b],beams[b])    
        ∂beams -= Tangent{BT}(;beamnames[b] => dbeam)
        @inbounds ∂y[2:4,2,b] .-= dy
    end
    dnode += Tangent{Bo}(;x = sum(ȳ[2,:]),y =  sum(ȳ[3,:]),ϕ = sum(ȳ[1,:]))
    return dnode,∂beams,∂y
end 

function reduceposat_back!(ȳ,∂y,y,
    beamnbrs::AbstractVector{Int},node::Bo,beams::BT,
    ∂beams = CRC.zero_tangent(beams)) where{beamnames,Tu,Bo,BT<:NamedTuple{beamnames,Tu}}

    for (n,b) in enumerate(beamnbrs)
        _,dbeam,dy,_ = scalepos_back(ȳ[:,n],y[2:4,2,b],beams[b])     
        ∂beams -= Tangent{BT}(;beamnames[b] => dbeam)
        @inbounds ∂y[2:4,2,b] .-= dy
    end
    return ∂beams
end 

function CRC.rrule(::typeof(reduceposat),node::Bo,beams::BT,y,beamnbrs) where{T,Bo<:Boundary{T},BT<:NamedTuple} 
    sol = reduceposat(node,beams,y,beamnbrs)
    
    function reducepose_back(ȳ,∂y = zero(y))
        dnode,∂beams,∂y =  reduceposat_back(ȳ,∂y,y,beamnbrs,node,beams)
        return NoTangent(),dnode,∂beams,∂y,NoTangent()
    end
    return sol, reduceposat_back
end 

function reduceposat_back(  ȳ,∂y,∂f,y,facs,beamnbrs::AbstractVector{Int},
                            node::Bo,beams::BT,
                            ∂beams= CRC.zero_tangent(beams),dnode = CRC.zero_tangent(node)) where{Bo,BT}     
    for (n,b) in enumerate(beamnbrs)
        
        sp, spback = CRC.rrule(scalepos,beams[b],y[2:4,2,b],Val(2))
        
        _,dbeam,dy,_ =  spback(facs[b] .* ȳ[:,n])  
        ∂beams -= Tangent{NamedTuple{beamnames,BT}}(;beamnames[b] => dbeam)
        ∂y[2:4,2,b] .-=  dy
        
        ∂f[b] .= @SVector [  (node.ϕ - sp[1]) * ȳ[1,n],
                            (node.x - sp[2]) * ȳ[2,n],
                            (node.y - sp[3]) * ȳ[3,n]
        ] 
        
    end
    dnode += Tangent{Bo}(;x = sum(facs[b] * ȳ[2,:]), y = sum(facs[b] * ȳ[3,:]) ,ϕ =  sum(facs[b] * ȳ[1,:])) 
    return dnode,∂beams,∂y,∂f
end 


function CRC.rrule(::typeof(reduceposat),node::Bo,beams,y::AbstractArray{T},facs,beamnbrs) where{T,Bo<:Boundary{T}} 
    sol = reduceposat(node,beams,y,facs,beamnbrs)
    function reducepos_back(ȳ,∂y,∂f,y,facs,beamnbrs)
        
        dnode =  reduceposat_back!(ȳ,∂y,∂f,y,facs,beamnbrs,node,beams)
        return NoTangent(),dnode,∂beams,∂y,∂f, NoTangent()
    end
   
    return sol, reduceposat_back
end 

# function pullback_normfactor(ȳ,b::BeamElement{T}) where{T}
#     nom = normfactor(b)
#     ∂b = Tangent{typeof(b)}(;
#         h = (-3  * ȳ * nom / b.h) ,
#         w = - ȳ * nom / b.w ,
#         E = - ȳ * nom / b.E
#     )
#     return (NoTangent(), ∂b)
# end

# function CRC.rrule(::typeof(normfactor), b::B) where{T,B<:BeamElement{T}}
#     y = normfactor(b)
#     pullback_norm(ȳ) = pullback_normfactor(ȳ,b)
#     return y, pullback_norm
# end

function pullback_normfactor_m(ȳ,b::B,nom) where{T,B<:BeamElement{T}}
    
    ∂b = Tangent{B}(;
        l = ȳ  * nom ,
        h = (-3 * b.l * ȳ * nom / b.h) ,
        w = -b.l * ȳ * nom / b.w ,
        E = -b.l * ȳ * nom / b.E
    )
    return (NoTangent(), ∂b)
end

function CRC.rrule(::typeof(normfactor_m), b::B) where{T,B<:BeamElement{T}}
    nf = normfactor(b)
    pullback_norm_m(ȳ) = pullback_normfactor_m(ȳ,b,nf)  
    return nf, pullback_norm_m
end

function pullback_normfactor_f(ȳ,b::B) where{T,B<:BeamElement{T}}
    ∂m = b.l * ȳ 
    fm, pb_m = CRC.rrule(normfactor_m, b)
    _, ∂b_m = pb_m(∂m)
    
    ∂b = Tangent{B}(;
        l = 2 * ∂b_m.l,
        h =  ∂b_m.h ,
        w =  ∂b_m.w ,
        E =  -∂b_m.E 
    )
    return (NoTangent(), ∂b)
end

function CRC.rrule(::typeof(normfactor_f), b::BeamElement{T}) where{T}
    y = normfactor(b)
    pullback_norm_f(ȳ) = pullback_normfactor_f(ȳ,b)
    return y, pullback_norm_f
end

function pullback_normvector(ȳ,b::BeamElement{T},nom) where{T}
    all(ȳ .== 0) && return (NoTangent(),ZeroTangent())
    # ȳ = CRC.unthunk(ȳ)
    # Ableitung von m nach b-Parametern
    ∂m_l = nom / b.l
    ∂m_h = -3 * nom / b.h
    ∂m_w = - nom / b.w
    ∂m_E = - nom / b.E
    # Ableitung von f = b.l * m nach b-Parametern

    ∂f_l = 2 * nom
    ∂f_h = b.l * ∂m_h # -3f/b.h = -3 * nom * b.l / b.h 
    ∂f_w = b.l * ∂m_w
    ∂f_E = b.l * ∂m_E
    dy_f = ȳ[2] + ȳ[3]
    # Zusammenführen (ȳ = [ȳ_m, ȳ_f1, ȳ_f2])
   ∂beam = Tangent{typeof(b)}(;
        l = ȳ[1] * ∂m_l + dy_f * ∂f_l,
        h = ȳ[1] * ∂m_h + dy_f * ∂f_h,
        w = ȳ[1] * ∂m_w + dy_f * ∂f_w,
        E = ȳ[1] * ∂m_E + dy_f * ∂f_E,
    )
    return (NoTangent(), ∂beam)
end

function CRC.rrule(::typeof(normvector), b::B) where{T,B<:BeamElement{T}}
    y = normfactor_m(b)

    pullback_norm(ȳ) = pullback_normfactor(ȳ,b,y)
    return y .* [b.l,b.l^2,b.l^2], pullback_norm
end

function pullback_scaleforce(ȳ,y,b,::Val{false})
    m = normfactor_m(b)
    f = b.l * m
    ∂y = @SVector [ȳ[1] / m, ȳ[2] / f,ȳ[3] / f]
    return (NoTangent(),∂y,ZeroTangent())
end

function pullback_scaleforce(ȳ,y,b,::Val{true})

    m = normfactor_m(b)
    f = b.l * m
    ∂y = @SVector [-ȳ[1] / m,-ȳ[2] / f,-ȳ[3] / f]
    ∂nv = @SVector [-ȳ[1] * y[1] / m^2, -ȳ[2] * y[2] / f^2, -ȳ[3] * y[3] / f^2]
    ∂beam = pullback_normvector(∂nv,b,m)[2]
    return (NoTangent(),∂beam,∂y)
end 

function CRC.rrule(::typeof(scaleforce),b::B,y) where{T,B<:BeamElement{T}}
    nv , pb_nv = CRC.rrule(normvector, b)

    result =  y ./ nv
    
    function pullback(ȳ)
        ∂y = ȳ ./ nv
        ∂nv = - ∂y .* y ./ nv 
        _, ∂beam = pb_nv(∂nv)
        return (NoTangent(),∂beam,∂y)
    end
    
    return result, pullback
end

function gettangent(beam::B,node,xb,yb,θb,κb::AbstractVector{T}) where{B<:BeamElement,T<:Real} 
    
    dl = -node.x/beam.l^2 * xb - node.y/beam.l^2 * yb - only(κb * beam.κ0)
    dκ0 = only(beam.l * κb)
    return Tangent{B}(;l = dl,κ0 = dκ0,θs = θb)
end 

function gettangent(beam::B,node,xb,yb,θb,κb::AbstractVector{T}) where{B<:CurvedBeam,T<:Real}

    dl = -node.x/beam.l^2 * xb - node.y/beam.l^2 * yb 
    if length(κb) > 1
        dl -= reduce(+,beam.κ0 .* κb)
        dκ0 = -beam.l .* κb
    else
        dκ0 = zeros(T,5)
    end

    return Tangent{B}(;l = dl,κ0 = dκ0,θs = θb)
end 

function pullback_init_beam(ȳ,node::No,beam::B,pars) where{BT,B<:BeamElement{BT},No<:Boundary{T} where{T}}
    @inbounds mb,θb,xb,yb,fxb,fyb,κb... = ȳ
    @inbounds m,fx,fy = pars 
    normf = normfactor_m(beam)
    _,∂beam = pullback_normvector([mb * m, fxb * fx,fyb * fy],beam,normf)
    ∂beam += gettangent(beam,node,xb,yb,θb,κb)
    return ∂beam
end 

function pullback_init_node(ȳ,::Type{No},beam::B) where{BT,B<:BeamElement{BT},No<:Boundary{T} where{T}}
    @inbounds θb,xb,yb = ȳ
    ∂node = Tangent{No}(;x = xb/beam.l,y=yb/beam.l,ϕ = θb)
    return ∂node
end 
function pullback_init_forces(ȳ,beam::B,pars) where{BT,B<:BeamElement{BT}}
    @inbounds mb,fxb,fyb = ȳ
    @inbounds m,fx,fy = pars 
    normf = normfactor(beam)
    ∂pars =  @SVector [mb * normf * beam.l,fxb * normf * beam.l^2,fyb * normf * beam.l^2]
    return ∂pars
end 

function CRC.rrule(::typeof(initialize_beam),node::No,beam::BeamElement{BT},pars::AbstractVector) where{BT,No<:Boundary}
    x,y,θ0, = node.x,node.y,node.ϕ
    l,θs,κ = beam.l,beam.θs,beam.κ0
    m,fx,fy = pars .* normvector(beam) #am Balkenelement
    result = @SVector [m,θ0 + θs,x./l,y./l,fx,fy,κ*l]
    function init_beam_back(ȳ) 
        dbeam = pullback_init_beam(ȳ,node,beam,pars)
        dnode = pullback_init_node(ȳ,typeof(node),beam)
        dpars = pullback_init_forces(ȳ,beam,pars)
        return NoTangent(),dnode,dbeam,dpars
    end 
    return result, init_beam_back
end 

function init_u0_pullback(ȳ,beams::BT,nodes::NT,pars,nodepos,i,symbolcache) where{BT<:NamedTuple,NT<:NamedTuple}
           
    ∂pars = zero(pars)
    dbeam = pullback_init_beam(ȳ,nodes[nodepos[i]],beams[i],pars[:,i])
    dnode = pullback_init_node(ȳ,typeof(nodes[nodepos[i]]),beams[i])
    ∂pars[:,i] += pullback_init_forces(ȳ,beams[i],pars[:,i])
    ∂beams = Tangent{typeof(beams)}(;symbolcache[1][i] => dbeam)
    ∂nodes = Tangent{typeof(nodes)}(;symbolcache[2][nodepos[1]] => dnode)
    return NoTangent(), ∂beams,∂nodes,∂pars,NoTangent(),NoTangent()
end 

function CRC.rrule(::typeof(initialize_beam),beams::NamedTuple,nodes::NamedTuple,pars,x_idxs::Number,i::Int)
    result = initialize_beam(nodes[x_idxs],beams[i],pars[:,x_idxs])
    symbolcache = (keys(beams),keys(nodes))
    init_u0_back(ȳ) = init_u0_pullback(ȳ,beams,nodes,pars,x_idxs,i,symbolcache)
    result,init_u0_back
end 

function CRC.rrule(::typeof(initialize_beam),beams::NamedTuple,nodes::NamedTuple,pars,idxs::AbstractVector,i::Int)
    x_idxs = idxs[i]
    result = initialize_beam(nodes[x_idxs],beams[i],pars[:,x_idxs])
    symbolcache = (keys(beams),keys(nodes))  
    init_u0_back(ȳ) = init_u0_pullback(ȳ,beams,nodes,pars,idxs,i,symbolcache)
    result,init_u0_back
end 

function forcesbackatend!(∂y,∂beams::CRC.AbstractTangent,ȳ,y,beams::NamedTuple{names,BT},idxs) where{names,BT}
    isempty(idxs) && return nothing
    y_ = @view y[[1,5,6],2,idxs]
    ∂y_ = @view ∂y[[1,5,6],2,idxs] 
    for (n,b) in enumerate(idxs) 
        _,dbeam,dy = pullback_scaleforce(ȳ,y_[:,n],beams[b],Val(true))
        ∂beams += Tangent{NamedTuple{names,BT}}(;names[b] => dbeam)
        ∂y_[:,n] .-= dy
    end
    return ∂beams
end 

function forcesbackatstart!(∂y,∂beams::CRC.AbstractTangent,ȳ,y,beams::NamedTuple{names,BT},idxs) where{names,BT}
    isempty(idxs) && return nothing
    y_ = @view y[[1,5,6],1,idxs]
    ∂y_ = @view ∂y[[1,5,6],1,idxs] 
    for (n,b) in enumerate(idxs) 
        _,dbeam,dy = pullback_scaleforce(ȳ,y_[:,n],beams[b],Val(true))
        ∂beams -= Tangent{NamedTuple{names,BT}}(;names[b] => dbeam)
        ∂y_[:,n] .+= dy
    end
    return ∂beams
end

function force_bound_back(ȳ,∂y,y,beams,idxs,∂beams::Tangent{BT} = CRC.zero_tangent(beams)) where{BT}
    ∂beams = forcesbackatend!(∂y,∂beams,ȳ,y,beams,idxs[1])
    ∂beams = forcesbackatstart!(∂y,∂beams,ȳ,y,beams,idxs[2])
    return ∂beams 
end 

# function CRC.rrule(::typeof(reduceforceat),node::No,beams,y,idxs) where{No}
    
#     sol = reduceforceat(node,beams,y,idxs)
#     function force_bound_back(ȳ,∂beams::Tangent{BT} = CRC.zero_tangent(beams)) where{BT}
        
#         ∂y = CRC.zero_tangent(y)
#         ∂beams = forcesbackatend!(∂y,ȳ,y,beams,idxs[2])
#         ∂beams = forcesbackatstart!(∂y,∂beams,ȳ,y,beams,idxs[1])

#         ∂node = ZeroTangent() 
#         return NoTangent(),∂node,∂beams,∂y,NoTangent()
#     end

#     return sol,force_bound_back
# end   

function CRC.rrule(::typeof(reduceforceat),node::No,beams,y::AbstractArray{T,3},factors,idxs) where{T,No}

    sol = reduceforceat(node,beams,y,factors,idxs)
    function force_bound_back(ȳ,∂beams::Tangent{BT} = CRC.zero_tangent(beams)) where{BT}
        
        y_idxs = getforceindices(idxs)
        ∂y = CRC.zero_tangent(y)
        ∂fac = CRC.zero_tangent(factors)
        y_ = @view y[y_idxs]
        ∂y_ = @view ∂y[y_idxs] 
        n = 1
        for b in idxs[1] 
            _,_,dy = pullback_scaleforce(factors[b].*ȳ,y_[:,n],beams[b],Val(false))
            ∂y_[:,n] .+= dy
            n += 1
        end
        for b in idxs[2]
            _,dbeam,dy = pullback_scaleforce(factors[b] .* ȳ,y_[:,n],beams[b],Val(true))
            ∂beams += Tangent{BT}(;Symbol(:Beam_,b) => dbeam)
            ∂y_[:,n] .+= dy 
            n += 1 
        end 
        
        ∂node = ZeroTangent() 
        return NoTangent(),∂node,∂beams,∂y,∂fac,NoTangent()
    end 
    return sol,force_bound_back
end 

# function CRC.rrule(::typeof(reduceforceat),node::LinearSlider{No},beams,y::AbstractArray{T,3},factors,idxs) where{T,No}
    
#     sol = reduceforceat(node,beams,y,factors,idxs)
#     function force_bound_back(ȳ)
#         y_idxs = getforceindices(idxs)
#         ∂y = zero(y)
#         ∂fac = zero(factors)
#         y_ = @view y[y_idxs]
#         ∂y_ = @view ∂y[y_idxs] 
        
#         ∂beams = Tangent{typeof(beams)}(;ntuple(x->keys(beams)[x]=>ZeroTangent(),length(beams))...)
#         for (b,yb,dyb) in zip(vcat(idxs...),eachcol(y_),eachcol(∂y_))
            
#             f, rr_norm = CRC.rrule(scaleforce,beams[b],yb)
#             if b in idxs[2]
#                 _,dbeam,dyrr_ = rr_norm(-ȳ) 
#                 ∂fac[b] = sum(trans .* f .* ȳ) 
#             else
#                 _,dbeam,dyrr_ = rr_norm(ȳ)
#                 ∂fac[b] = -sum(trans .* f .* ȳ) 
#             end 
#             ∂beams += Tangent{typeof(beams)}(;Symbol(:Beam_,b) => dbeam)
#             dyb .= dyrr_ #./ length(idxs)
#         end 
#         ∂node = CRC.zero_tangent(LinearSlider{No})#;fx = ȳ[2],fy = ȳ[3],mz = ȳ[1])
#         return NoTangent(),∂node,∂beams,∂y,∂fac,NoTangent()
#     end 
#     return sol,force_bound_back
# end 



@non_differentiable reduceforceat(n::Clamp,beams,y,idxs) 

function CRC.rrule(::typeof(residuals!),residuals,str::Structure,y,beamtpl::BT,nodetpl::NT) where{BT<:NamedTuple,NT<:NamedTuple}
    ind = 1
    adj = str.AdjMat
    branches = filter(x->canchangeposition(nodetpl[x]),eachindex(nodetpl))
    nbr_branches = length(branches)
    residuals_forces = @view residuals[:,1:nbr_branches]
    residuals_positions = @view residuals[:,nbr_branches+1:end]
    forces = 1

    for (node,beams) in beamsatnode(adj,nodetpl,beamtpl)
        
        if forcesatnode(nodetpl[node])
            reduceforceat!(view(residuals_forces,:,forces),nodetpl[node],beamtpl,y,beams)
            forces += 1
        end 
        if !isempty(beams[2])
            reduceposat!(view(residuals_positions,:,beams[2]),nodetpl[node],beamtpl,y,beams[2])
        end  
    end

    function residuals!_back(ȳt)
        ȳ = CRC.unthunk(ȳt)
        ∂res = ZeroTangent()
        ∂y = zero(y) 
        dstr = NoTangent()

        ȳ_forces = @view ȳ[:,1:nbr_branches] #similar(y,Float32,3,nbr_branches)         #
        ȳ_positions = @view ȳ[:,nbr_branches+1:end] # similar(y,Float32,3,length(beamtpl))   #
  
        ∂beams = CRC.zero_tangent(beamtpl)       
        ∂nodes = CRC.zero_tangent(nodetpl)
        forcenbr = 1

        for (node,beams) in beamsatnode(adj,nodetpl,beamtpl)
            
            if forcesatnode(nodetpl[node])
                ∂beams = force_bound_back(view(ȳ_forces,:,forcenbr),∂y,y,beamtpl,beams,∂beams)        
                forcenbr += 1
            end
            if !isempty(beams[2])
                ∂beams = reduceposat_back!(view(ȳ_positions,:,beams[2]),∂y,y,beams[2],nodetpl[node],beamtpl,∂beams)
                ∂nodes += Tangent{NT}(;node =>  Tangent{typeof(nodetpl[node])}(;x = sum(ȳ_positions[2,beams[2]]),y =  sum(ȳ_positions[3,beams[2]]),ϕ = sum(ȳ_positions[1,beams[2]])))
               
                # @show ∂nodes[node]
            end
            
        end 

        return NoTangent(),∂res,dstr,∂y,∂beams,∂nodes
    end 

    return residuals,residuals!_back
end 


function CRC.rrule(::typeof(residuals!),residuals::Matrix,str::Structure,y::EnsembleSolution,beamtpl::BT,nodetpl::NT) where{BT<:NamedTuple,NT<:NamedTuple}
    y = toArray(y)
    res,back = rrule(residuals!,residuals,str,y,beamtpl,nodetpl)
    return res,back
end

function CRC.rrule(::typeof(residuals!),residuals::AbstractMatrix,adj_::AbstractMatrix{T},y::AbstractArray,beamtpl::BT,nodetpl::NT) where{T,BT<:NamedTuple,NT<:NamedTuple}

    ind = 1
    idcs = getindices(size(adj_,1))
    adj = ifelse.(adj_ .> 1,one(T),adj_ )
    adj = ifelse.(adj .< 0, zero(T),adj_ )

    branches = filter(x->canchangeposition(nodetpl[x]),eachindex(nodetpl))
    nbr_branches = length(branches)
    residuals_forces = @view residuals[:,1:nbr_branches]
    residuals_positions = @view residuals[:,nbr_branches+1:end]

    forces = 1
    positions = 1

    rrforce = Vector{Function}(undef,length(branches))
    rrpos = Vector{Pair{AbstractRange{Int},Function}}()

    for (node,beams) in beamsatnode(adj,nodetpl,beamtpl)
        
        if forcesatnode(nodetpl[node])
            reduceforceat!(view(residuals_forces,:,forces),nodetpl[node],beamtpl,y,beams)
            forces += 1
        end 
        if !isempty(beams[2])
            reduceposat!(view(residuals_positions,:,beams[2]),nodetpl[node],beamtpl,y,beams[2])
        end  
    end   
    function residuals!_back(ȳ)
        ȳ_ = CRC.unthunk(ȳ)
        ∂res = @thunk(CRC.zero_tangent(residuals))
        ∂y = zero(y)
        ∂fac = zero(adj)
        ∂beams = CRC.zero_tangent(beamtpl)         
        ∂nodes = CRC.zero_tangent(nodetpl) 

        ȳ_forces = @view ȳ_[:,1:nbr_branches]
        ȳ_positions = @view ȳ_[:,nbr_branches+1:end]
        forcenbr = 1

        for (node,beams) in beamsatnode(adj,nodetpl,beamtpl)
            
            if forcesatnode(nodetpl[node])
                force_bound_back(view(ȳ_forces,:,forcenbr),∂y,y,beamtpl,beams,∂beams)        
                forcenbr += 1
            end
            if !isempty(beams[2])
                reduceposat_back!(view(ȳ_positions,:,beams[2]),∂y,y,beams[2],nodetpl[node],beamtpl,∂beams)
                @show ∂nodes += Tangent{NT}(;node =>  Tangent{typeof(nodetpl[node])}(;x = sum(ȳ_positions[2,beams[2]]),y =  sum(ȳ_positions[3,beams[2]]),ϕ = sum(ȳ_positions[1,beams[2]])))
                
                # @show ∂nodes[node]
            end
        end 
        return NoTangent(),∂res,∂fac,∂y,∂beams,∂nodes
    end    
    return residuals,residuals!_back
end 

function CRC.rrule(::typeof(residuals!),residuals::Matrix,adj::AbstractMatrix,y::EnsembleSolution,beams,nodes)
    y = toArray(y)
    res,back = rrule(residuals!,residuals,adj,y,beams,nodes)
    return res,back
end 

function CRC.rrule(::Type{BT}, x,y,θ,fx,fy,mz) where{BT<:Boundary{T}} where{T}
    Bar_pullback(ȳ) = (NoTangent(),ȳ.x,ȳ.y,ȳ.θ,ȳ.fx,ȳ.fy,ȳ.mz)
    return BT(x,y,θ,fx,fy,mz), Bar_pullback
end
function CRC.rrule(::Type{Clamp{T}}, x,y,θ,fx,fy,mz) where{T}
    Bar_pullback(ȳ) = (NoTangent(),ȳ.x,ȳ.y,ȳ.θ,ȳ.fx,ȳ.fy,ȳ.mz)
    return Clamp{T}(x,y,θ,fx,fy,mz), Bar_pullback
end

function CRC.rrule(::typeof(+),a::BT,b::BT) where{BT<:Boundary{T}} where{T}
    backaddb(ȳ) = (NoTangent(), ȳ,ȳ)
    return a + b,backaddb
end 

function CRC.rrule(::typeof(addposition),mv::LinearSlider{BT},disp) where{BT}
    res = addposition(mv,disp)
    function addpos_back_x(ȳ)
        return NoTangent(),ȳ,ȳ.s
    end
    return res,addpos_back_x
end 

function CRC.rrule(::typeof(addposition),mv::Joint{BT},disp) where{BT}
    res = addposition(mv,disp)
    function addpos_back_x(ȳ)
        return NoTangent(),ȳ,ȳ.ϕ
    end
    return res,addpos_back_x
end 

function CRC.rrule(::typeof(addpositions),nodes::NamedTuple{names,NT},xpos::AbstractMatrix{T}) where{names,T,NT}
    # branches = Vector{Int}
    newnodes = Vector{Pair{Symbol,Boundary}}(undef, size(xpos,2))
    n = 0
    pullbacknodes = Vector{Symbol}(undef,size(xpos,2))
    for node in names
        if canchangeposition(nodes[node])            
            n += 1
            newnode = addposition(nodes[node],xpos[:,n])
            @inbounds pullbacknodes[n] = node
            @inbounds newnodes[n] = node => newnode
        end 
    end
    function addpos_back(ȳ,∂xpos = zero(xpos)) 
        for (colx,node) in zip(eachcol(∂xpos),pullbacknodes)
            @inbounds colx[1] += ȳ[node].x
            @inbounds colx[2] += ȳ[node].y
            @inbounds colx[3] += ȳ[node].ϕ
        end
        return NoTangent(),ZeroTangent(),∂xpos
    end 
    return merge(nodes,newnodes),addpos_back
end

function split_primals(x::AbstractVecOrMat{T},beam::B) where{T,B<:BeamElement{T}}
    λ_jvp = @view x[1:7]

    # λ_ode = @view x[end-6:end]
    forces = @view x[[1,5,6]]
    positions = @view x[2:4]
    return λ_jvp,forces,positions
end 

function split_primals(x::AbstractVecOrMat{T},beam::CurvedBeam) where{T}
    λ_jvp = @view x[1:6+length(beam.κ0)]
    # @show λ_jvp[7:end] 
    forces = @view x[[1,5,6]]
    positions = @view x[2:4]
    # @show x[end-6:end]
    return λ_jvp,forces,positions
end 

function output_func(u,beams::NamedTuple{beamnames,BT},nodes::NamedTuple{nodenames,NT},nodepos,x,i) where{beamnames,nodenames,BT,NT}
    x_idxs = nodepos[i]

    λ_jvp,λ_forces,λ_positions = split_primals(u[end],beams[i])
    ∂xforces = pullback_init_forces(λ_forces,beams[i],x[:,i]) 
    ∂beams = Tangent{NamedTuple{beamnames,BT}}(;beamnames[i] => pullback_init_beam(λ_jvp,nodes[x_idxs],beams[i],x[:,i]))
    ∂nodes = Tangent{NamedTuple{nodenames,NT}}(;nodenames[x_idxs] => pullback_init_node(λ_positions,typeof(nodes[x_idxs]),beams[i]))
    (∂beams,∂nodes,∂xforces),false
end 

function reduction_func!(u,data,I)
    
    for ((∂beam,∂node,∂xforces),i) in zip(data,I)
        
        @inbounds u[1][:,i] .+= ∂xforces
        u[2] += ∂beam
        u[3] += ∂node
    end 
    u,false
end 

@generated function generate_tangent(x::T,primal) where{T<:Union{BeamElement,Boundary}}

    quote
        CRC.Tangent{$T}(;$CRC.backing(primal)...)  
    end

end 
generate_tangent(x::T,primal::ZeroTangent) where{T<:Union{BeamElement,Boundary}} = primal



@generated function generate_tangent(::Val{N}, x::NamedTuple{Nt,BN}, primals) where {N,Nt,BN<:NTuple{N,Union{BeamElement,Boundary}}}

    fields = Expr(:tuple)
    for i = 1:N
        push!(fields.args, :(Nt[$i] => $generate_tangent(x[$i],primals[$i])))
    end
    quote
        CRC.Tangent{NamedTuple{$Nt,$BN}}(;$fields...)  
    end
end

function u_init(x, beams,nodes,dbeams,dnodes)
   
    [x,
    generate_tangent(Val(length(beams)), beams, dbeams),
    generate_tangent(Val(length(nodes)), nodes, dnodes)
    ]
end 


function initialize_u0_vjp(x::AbstractVector{T},solforwardatend,beam::Beam) where{T}
    u0 = Vector{T}(undef,14)
    @inbounds u0[1:7] .= x # Primals 
    @inbounds u0[8:14] .= solforwardatend # sol at end of beam to diff backwards 
    return u0,SciMLBase.NullParameters()
end 

function initialize_u0_vjp(solforwardatend::AbstractArray{T,N},beam::Beam) where{T,N}
    u0 = zeros(T,14,3)
    for i in 1:3
        @inbounds u0[1+i,i] = one(T)
        @inbounds u0[8:14,i] .= solforwardatend # sol at end of beam to diff backwards 
    end
    return u0,SciMLBase.NullParameters()
end 

function initialize_u0_vjp(x::AbstractVector{T},solforward::AbstractArray{T,N},beam::CurvedBeam{T}) where{T,N}
    u0 = Vector{T}(undef,13 + length(beam.κ0))
    fill!(u0,zero(T))
    @inbounds u0[1:7] .= x # Primals
                # u0[2:4] .= zero(T)
    @inbounds u0[end-6:end] .= solforward # sol at end of beam to diff backwards
    # @inbounds u0[7:end-7] .-= beam.κ0.*beam.l #curvature pars to diff backwards
    return u0,beam.κ0.*beam.l
end 

function make_vjp_func(∂sol,solforward::Array{T,3},beams::BT) where{BT,T}
    (prob,i,repeat) -> begin
        u0,p = initialize_u0_vjp(view(∂sol,:,2,i),view(solforward,:,2,i),beams[i])
        remake(prob;u0 = u0,p = p)
    end
    
end

function make_vjp_func(solforward::Array{T,3},beams::BT) where{BT,T}
    (prob,i,repeat) -> begin
        u0,p = initialize_u0_vjp(solforward[:,2,i],beams[i])
        remake(prob;u0 = u0,p = p)
    end
    
end

function make_output_func(beams::BT, nodes_::NT, nodepos,xforces) where{BT,NT}
    (sol,i) -> output_func(sol,beams,nodes_,nodepos,xforces,i)
end

function CRC.rrule(str::Structure,x::AbstractArray{T,N},beams::NamedTuple{beamnames,BT},nodes::NamedTuple{nodenames,NT}) where{T,N,beamnames,BT,nodenames,NT}
    
    nodepos = getstartnodes(str)
    anz = count(x->canchangeposition(x),values(nodes))
    xpos = @view x[:,1:anz]
    xforces = @view x[:,anz+1:end]
    out = Array{T,3}(undef,7,2,length(beams))
    nodes_,change_pullback = CRC.rrule(addpositions,nodes,xpos)
    
    prob_func = make_prob_func(beams, nodes_, xforces, nodepos)

    ensprob = EnsembleProblem(prob;prob_func = prob_func,
                                output_func = output_function_,
                                reduction = reduction_funcF!,
                                u_init = out,
                              )

    solve(ensprob,str.Solver,
                EnsembleThreads(),
                save_on = false,save_start=true,save_end = true,
                reltol = 1e-6,abstol = 1e-6,
                trajectories = length(beams) 
                )
 
    function back_ode(ȳ)
        @inbounds ∂sol,∂Beams,∂Nodes = ȳ
        dsol = CRC.unthunk(∂sol)
        ∂x = CRC.zero_tangent(x)

        ∂xpos = @view ∂x[:,1:anz] # = zeros(T,size(x,1), anz)

        ∂xforces = @view ∂x[:,anz+1:end] # similar(x, size(x,1), size(x,2) - anz) 
        for i in 1:length(beams)
            @inbounds ∂xforces[:,i] .= dsol[[1,5,6],1,i] .* normvector(beams[i])
        end

        prob_funcjvp = make_vjp_func(dsol,out,beams)
        output_func = make_output_func(beams, nodes_, nodepos, xforces)
        ensprobjvp  =  EnsembleProblem(vjpprob;prob_func = prob_funcjvp,
                                    output_func =output_func,
                                    reduction = (u,data,I) -> reduction_func!(u,data,I),
                                    u_init = u_init(∂xforces,beams,nodes_,∂Beams,∂Nodes),
                                    )

        solp = solve(ensprobjvp,str.Solver,
                    EnsembleThreads(),
                    save_on = false,save_start=true,save_end = true,
                    reltol = 1e-6,abstol = 1e-6,
                    trajectories = length(beams)
                    )
        
        @inbounds ∂xforces,dBeams,dNodes = solp.u
        # @show dNodes
        for i in 1:length(beams)
            nom = normfactor_m(beams[i])
            d = @view dsol[:,1,i]
            dNodes += Tangent{NamedTuple{nodenames,NT}}(;nodenames[i] => pullback_init_node(d[2:4],typeof(nodes[nodepos[i]]),beams[i]))
            dBeams += Tangent{NamedTuple{beamnames,BT}}(;beamnames[i] => pullback_init_beam(d,nodes_[nodepos[i]],beams[i],xforces[:,i]))
        end
        # @show dNodes
        change_pullback(dNodes,∂xpos)
        # add_startgrads!(∂x,∂sol[[3,4,2],1,:],bn.Beams,bn.Nodes,)
        @inbounds ∂x[:,1:anz] .= ∂xpos
        @inbounds ∂x[:,anz+1:end] .= ∂xforces
        
        return NoTangent(),∂x,dBeams,dNodes,ZeroTangent(),NoTangent()
    end
    (out,beams,nodes_),back_ode
end 

function Base.inv(x::Tangent{T,NT}) where{T<:Beam,NT}
    fields = fieldnames(T)
    CRC.StructuralTangent{T}((;map((f) -> f => Float32(inv(getproperty(x,f))), fields)...))
end
function Base.sqrt(x::Tangent{T,NT}) where{T<:Beam,NT}
    fields = fieldnames(T)
    CRC.StructuralTangent{T}((;map((f) -> f => Float32(sqrt(getproperty(x,f))), fields)...))
end
function Base.:*(x::Tangent{T},y::Tangent{T}) where{T}
    fields = fieldnames(T)
    CRC.StructuralTangent{T}((;map((f) -> f => Float32(getproperty(x,f) * getproperty(y,f)), fields)...))
end
function Statistics.realXcY(x::Tangent{T,NT},y::Tangent{T,NT}) where{T<:Beam,NT}
    x * y 
end

function CRC.rrule(str::GroundStructure,x::AbstractMatrix{T},beams::NamedTuple{beamnames,BT},nodes::NamedTuple{nodenames,NT},adj) where{T,BT,NT,nodenames,beamnames}
    
    nodepos = getstartnodes(adj)
    anz = count(x->canchangeposition(x),values(nodes))
    xpos = @view x[:,1:anz]
    xforces = @view x[:,anz+1:end]
    cbeams = size(adj,2) * (size(adj,2) - 1) ÷ 2
    out = Array{T,3}(undef,7,2,length(beams))
    
    nodes_,change_pullback = CRC.rrule(addpositions,nodes,xpos)
    
    prob_func = make_prob_func(beams, nodes_, xforces, nodepos)

    ensprob  =  EnsembleProblem(prob;prob_func = prob_func,
                                    output_func = output_function_,
                                    reduction = reduction_funcF!,
                                    u_init = out,
                                )

    sol = solve(ensprob,str.Solver,
                EnsembleThreads(),
                reltol = 1e-6,abstol = 1e-6,
                save_on = false,save_start=true,save_end = true,
                trajectories = cbeams
                )   

    function back_groundstr(ȳ)
        @inbounds ∂sol,∂Beams,∂Nodes = ȳ
        dsol = CRC.unthunk(∂sol)
        ∂x = similar(x)
        ∂xpos = @view ∂x[:,1:anz]
        # ∂xforces .= ∂sol[[1,5,6],2,:]
        ∂xforces = @view ∂x[:,anz+1:end] # similar(x, size(x,1), size(x,2) - anz) 

        for i in 1:length(beams)
            @inbounds ∂xforces[:,i] .= dsol[[1,5,6],1,i] .* normvector(beams[i])
        end
        prob_func = make_vjp_func(∂sol,out,beams)
        #integriere Rückwärtsproblem
        ensprob  =  EnsembleProblem(vjpprob;prob_func = prob_func,
                                            output_func = (sol,i) -> output_func(sol,beams,nodes_,nodepos,xforces,i),
                                            reduction = (u,data,I) -> reduction_func!(u,data,I),
                                            u_init = u_init(∂xforces,beams,nodes,∂Beams,∂Nodes)
                                            )
        solp = solve(ensprob,str.Solver,
                    EnsembleThreads(),
                    reltol = 1e-6,abstol = 1e-6,
                    save_on = false,save_start=false,save_end = true,
                    trajectories = cbeams
                    )
        @inbounds ∂xforces,dBeams,dNodes = solp

        for i in 1:length(beams)
            d = @view dsol[:,1,i]
            dNodes += Tangent{NamedTuple{nodenames,NT}}(;nodenames[i] => pullback_init_node(d[2:4],typeof(nodes[nodepos[i]]),beams[i]))
            dBeams += Tangent{NamedTuple{beamnames,BT}}(;beamnames[i] => pullback_init_beam(d,nodes_[nodepos[i]],beams[i],xforces[:,i]))
        end
        ∂adj = ZeroTangent()
        change_pullback(dNodes,∂xpos)
        return NoTangent(),∂x,dBeams,dNodes,∂adj,ZeroTangent()
    end
    (out,beams,nodes_),back_groundstr
end 


function CRC.rrule(::typeof(admittance_matrix),solfw::AbstractArray{T,N},adj,str,beams) where{T,N}
    
    idxs =getindices(size(adj,1))
    lensol = length(beams)
    prob_func = make_vjp_func(solfw,beams)
    #dadj Gradient der Steifigkeitsmatrix berechnen  
    ensprob2  =  EnsembleProblem(vjpprob;prob_func = prob_func,
                                output_func = (bsol,i) ->  output_func_admittance(bsol,i,solfw),
                                reduction = (u,data,I) -> reduction_func_admittance_!(u,data,idxs[I],solfw,adj,beams),
                                u_init = [zeros(eltype(adj),3 .* size(adj)...),[]]
                                )
    sol = solve(ensprob2,str.Solver,
                EnsembleThreads(),
                reltol = 1e-6,abstol = 1e-6,
                save_on = false,save_start=true,save_end = true,
                trajectories = lensol
                )

    ad,d_out = sol.u
    function admittance_back(ȳt)
        ȳ = CRC.unthunk(ȳt)
        
        ∂adj = zeros(T,size(adj)...)
        Δd = Matrix{eltype(adj)}(undef,3,3)
        # dtmp = Matrix{eltype(adj)}(undef,3,3)
        for (id,beam,d0,sol) in zip(idxs,beams,d_out,eachslice(solfw,dims = 3))
            i_ = 3 * (id[1]-1)+1:3*id[1] 
            j_ = 3 * (id[2]-1)+1:3*id[2]
            x,y = sol[3:4,2] .- sol[3:4,1]

            Δd .=  -(view(ȳ,j_,i_)' + view(ȳ,i_,j_))
            # Pullback der Transformationen von d

            Δd[2, :] .-= y .* Δd[1, :]
            Δd[3,:]  .-= x .* Δd[1, :] 
            Δd[1, :] .*= -1
            Δd .+= (view(ȳ,i_,i_) + view(ȳ,j_,j_))

            # Pullback der elementweisen Multiplikation mit d0
            ∂adj[id] += sum(d0 .* Δd ) 
        end
        return NoTangent(),NoTangent(),∂adj,NoTangent(),NoTangent()
    end 
    ad,admittance_back
end 

function CRC.rrule(::typeof(effective_stiffness),k::AbstractMatrix{T},u) where{T}
    nomoveidcs = [u[1]...]
    moveables = [u[2]...]
    ids = setdiff(axes(k,1),nomoveidcs)
    ke = @view k[ids,ids]
    f = zeros(T,size(k,1))
    f[moveables] .= [one(T)]
    mov = ke \ f[ids]
    
    e_out = mov' * ke * mov
    function compute_energy_pullback(Δenergy)
        Δ_dot = CRC.unthunk(Δenergy)
        # Zurück in volle Matrix
        Δ_k = zero(k)
        Δ_k[ids, ids] .= - Δ_dot * mov * mov' 
        
        return CRC.NoTangent(),Δ_k, NoTangent()
    end
    e_out,compute_energy_pullback
end 

function CRC.rrule(::typeof(changenode),bn,nodes,nt)
    bn_out = changenode(bn,nodes,nt)
    # @show bn
    # @show bn_out
    function changenodeback(ȳ)
        # @show ȳ
        return CRC.NoTangent(), ȳ, CRC.NoTangent(), CRC.NoTangent()
    end 
    return bn_out,changenodeback
end
