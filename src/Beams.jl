abstract type BeamElement{T} end

struct Beam{T<:Real} <:BeamElement{T}
    l::T
    h::T
    w::T
    κ0::T
    E::T
    θs::T
    θe::T
end 

struct CurvedBeam{T<:Real} <:BeamElement{T}
    l::T
    h::T
    w::T
    κ0::AbstractVector{T}
    E::T
    θs::T
    θe::T
    function CurvedBeam{T}(l::T,h::T,w::T,κ0::AbstractVector{T},E::T,θs::T,θe::T) where{T}
        # p = T.([l,h,w,E,θs,θe])
        # k0 = T.(κ0)
        new{T}(l,h,w,κ0,E,θs,θe)
    end 
end

# function CurvedBeam(l,h,w,κ0;E = 2.1f5,θs = 0)
#     p = promote(l,h,w,E,θs)
#     θe = get_θe(p[5],κ0)
#     CurvedBeam{eltype(p)}(p[1],p[2],p[3],κ0,p[4],p[5],θe)
# end

function CurvedBeam{T}(l,h,w,κ0;E = 2.1f5,θs = 0,θe = nothing) where{T}
    if isnothing(θe)
        θe = get_θe(p[5],k0)
    else
        p = T.([l,h,w,E,θs])
        θe = T(θe)
    end 
    k0 = T.(κ0)
    CurvedBeam{T}(p[1],p[2],p[3],k0,p[4],p[5],θe)
end

relu(x::T,m = zero(T)) where{T} = ifelse(x < m,m,x)

function Beam(l,h,w,κ0;E = 2.1f5,θs = 0,θe = θs + l*κ0)

    p = promote(l,h,w,κ0,E,θs,θe)
    Beam(p...)
end 

function Beam{T}(l,h,w,κ0;E = 2.1f5,θs = 0f0,θe = θs + l*κ0) where{T}

    Beam(T.([l,h,w,κ0,E,θs,θe])...)
end 

function Beam{T}(nt::NamedTuple) where{T}
    Beam((map(k->getfield(nt,k),fieldnames(Beam))...))
end 

function CurvedBeam{T}(nt::NamedTuple) where{T}
    CurvedBeam{T}((map(k->getfield(nt,k),fieldnames(CurvedBeam))...))
end 

function (b::Beam{T})(nt::AbstractVector) where{T}
    Beam(nt...)
end 

function Base.show(io::IO,beam::BeamElement)
    return println(io, "Beam with Length: $(beam[1]),width: $(beam[3]), height: $(beam[2]), curvature: $(beam[4]) and E: $(beam[5])")
end 

function change_beam(beam::Beam;kwargs...)
    for (field,value) in kwargs
        beam = Setfield.@set beam.$field = value 
    end 
    beam
end  

function isstraight(beam::Beam{T}) where{T}
    return isapprox(beam.κ0, 0, atol = 1e-6)
end

function isstraight(beam::CurvedBeam{T}) where{T}
    return isapprox(beam.κ0, zeros{T,5}, atol = 1e-6)
end

function EIz(beam::BeamElement{T}) where{T}
    Iz = beam.w * beam.h^3/12
    EIz =beam.E * Iz
end

function shouldbeambuckle(beam::B,y::AbstractArray{T,N}) where{N,T,B<:BeamElement{T}}
    tmp = true
    tmp &= isstraight(beam)    
    tmp || return tmp
    # α == angle in compressing direction?
    α = atan(y[6],y[5]) #force angle
    tmp &= all(isapprox.(α,beam.θs .+ π, atol = 1e-2))
    tmp || return tmp
    # a > critical load for buckling, beam straight and force angle in the compressing direction, then it should buckle
    a = hypot(y[5],y[6]) #force magnitude
    eiz = EIz(beam) 
    tmp &= a > (π^2 * eiz / beam.l^2) #critical load for buckling
    return tmp
end



gettype(::Beam{T}) where{T} = Beam
gettype(::CurvedBeam{T}) where{T} = CurvedBeam

Base.length(b::BeamElement) = 7
Base.lastindex(b::BeamElement) = 7

function Base.getproperty(b::BeamElement,n::Symbol)
    val = getfield(b,n)
end 

function Base.getindex(b::B, i::Int) where B<:BeamElement
    N = fieldcount(B)
    i > N &&  throw(BoundsError(b, i))
    
    fields = fieldnames(B)
    
    getfield(b, fields[i])
end

function Base.iterate(b::B, i::Int=1) where B<:BeamElement
    N = fieldcount(B)
    i> N && return nothing
    fields = fieldnames(B)
    return getfield(b, fields[i]), i+1
end

isacurvedbeam(::CurvedBeam) = true
isacurvedbeam(::BeamElement) = false

Base.promote_rule(::Type{B}, ::Type{S}) where {B<:BeamElement{T}, S} where T = promote_type(T,S)

Base.getindex(b::B,idx::AbstractVector) where{B<:BeamElement} = map(x->getindex(b,x),idx)
Base.getindex(b::B,idx::Colon) where{B<:BeamElement} = map(x->getindex(b,x),1:fieldcount(B))
Base.IteratorSize(b::B) where{B<:BeamElement} = Base.HasLength()
Base.real(b::B) where{B<:BeamElement} = b
Statistics.realXcY(a::B,b::B) where{B<:BeamElement} = a*b
Base.:*(a::Real,b::B) where{B<:BeamElement} = B(a .* b...)
Base.:*(b::B,a::Real) where{B<:BeamElement} = B(a .* b...)
Base.:*(a::B,b::B) where{B<:BeamElement} = B(a .* b...)


@generated function Base.:-(a::B,b::B) where{T,B<:BeamElement{T}} 
    fields = fieldnames(a)
    exprs = [:($(QuoteNode(f)) ∈ [:l,:h,:w] ? relu(getfield(a, $(QuoteNode(f))) .- getfield(b,$(QuoteNode(f))),1f-2) : (getfield(a, $(QuoteNode(f))) .- getfield(b,$(QuoteNode(f))))) for f in fields]
    return quote
        $(Expr(:call,Beam, exprs...))
    end
end 

Base.:-(a::Beam,b::Real) = Beam(a .- b...)
Base.:-(a::Real,b::Beam) = Beam(a .- b...)

Base.:abs2(b::Beam) = b*b
Base.:/(a::Real,b::Beam) = Beam(a ./ b...)
Base.:/(b::Beam,a::Real) = Beam(b ./ a...)
Base.:+(a::Beam,b::Beam) = Beam(a .+ b...)

Optimisers.functor(b::B) where{T,B<:BeamElement{T}} = (NamedTuple{fieldnames(B)}(b[1:7]),B)
Optimisers.init(o::Adam, x::Beam{T}) where{T} = (Beam{T}(zeros(T,7)...),Beam{T}(zeros(T,7)...), T.(o.beta))
Optimisers.init(o::Adam, x::CurvedBeam{T}) where{T} = (CurvedBeam{T}(zero(T),zero(T),zero(T),zeros(T,5),zero(T),zero(T),zero(T)),CurvedBeam{T}(zero(T),zero(T),zero(T),zeros(T,5),zero(T),zero(T),zero(T)), T.(o.beta))
Optimisers.init(o::WeightDecay, x::B) where{B<:BeamElement} = nothing
Optimisers.isnumeric(::B) where{T,B<:BeamElement{T}} = true
Optimisers.subtract!(a::Beam{T},b::Beam) where{T} = Beam{T}(a .- Beam{T}(merge(Optimisers.mapvalue(_->zero(T),Optimisers.functor(b)[1]),Optimisers.trainable(b)))...)
Optimisers.subtract!(a::CurvedBeam{T},b::CurvedBeam) where{T} = CurvedBeam{T}(a .- CurvedBeam{T}(merge(Optimisers.mapvalue(_->zero(T),Optimisers.functor(b)[1]),Optimisers.trainable(b)))...)
Optimisers.trainable(b::B) where{B<:BeamElement} = (;l = b.l,h = b.h,w = b.w,κ0 = b.κ0,E = b.E,θs = b.θs,θe = b.θe)

Optimisers._trainable(b::Beam{T},fr) where{T} = Beam{T}(merge(Optimisers.mapvalue(_ -> nothing, Optimisers.functor(b)[1]), Optimisers.trainable(b)))

Base.zero(::Beam{T}) where{T} = Beam(zeros(T,7)...)
Base.zero(::CurvedBeam{T}) where{T} = CurvedBeam(zeros(T,7)...)

BEAM_SCALE = (
    l = 1,
    h = 1,
    w = 1,
    κ0 = 1,
    E = 1,
    θs = 1,
    θe = 1
)


@generated function scale_beam(dx,b::B, scale) where{T,B<:BeamElement{T}} 
    fields = fieldnames(dx)
    # println(fields)'

    exprs = [:( getproperty(dx, $(QuoteNode(f))) * getproperty(scale,$(QuoteNode(f)) )) for f in fields]
    return quote
        $(Expr(:call, :B, exprs...))
    end
end

@generated function invscale_beams(dx,::B, scale) where{T,B<:BeamElement{T}}
    fields = fieldnames(dx)
    exprs = [:( getproperty(dx, $(QuoteNode(f))) / getproperty(scale, $(QuoteNode(f)))) for f in fields]
    return quote
        $(Expr(:call, :B, exprs...))
    end
end

@generated function combine(β,mt::B,dx) where{T,B<:BeamElement{T}}
    fields = fieldnames(mt)
    # exprs = [:(β * getproperty(mt, $(QuoteNode(f))) + (1-β) * getproperty(dx, $(QuoteNode(f)))) for f in fields]
    exprs = [:(isnothing(getfield(dx, $(QuoteNode(f)))) ? zero(T) : β * getfield(mt, $(QuoteNode(f))) + (1-β) * getfield(dx, $(QuoteNode(f)))) for f in fields]
    return quote
        $(Expr(:call,B, exprs...))
    end
end 

@generated function combineabs2(β,mt::B,dx) where{T,B<:BeamElement{T}}
    fields = fieldnames(mt)
    exprs = [:(isnothing(getfield(dx, $(QuoteNode(f)))) ? zero(T) : β * getfield(mt, $(QuoteNode(f))) + (1-β) * abs2.(getfield(dx, $(QuoteNode(f))))) for f in fields]
    return quote
        $(Expr(:call,B, exprs...))
    end
end 

@generated function combine(η,βt,mt::B,vt::B,ϵ) where{B<:BeamElement}
    fields = fieldnames(mt)
    exprs = [:(getfield(mt, $(QuoteNode(f))) ./ (1-βt[1]) ./ (sqrt.(getfield(vt, $(QuoteNode(f))) ./ (1 -βt[2])).+ ϵ).* η) for f in fields]
    return quote
        $(Expr(:call,B, exprs...))
    end
end 

function Optimisers.apply!(o::Adam,state,b::BeamElement{T},dx) where{T}
    η, β, ϵ = T(o.eta), T.(o.beta), T(o.epsilon)
    mt, vt, βt = state
    # dx_scaled = scale_beam(dx,mt,BEAM_SCALE)
    mt = combine(β[1],mt,dx)
    vt = combineabs2(β[2],vt,dx)
    dx_scaled = combine(η,βt,mt,vt,ϵ) #  mt / (1 - βt[1]) / (sqrt(vt / (1 - βt[2])) + ϵ) * η
    # dx′ = invscale_beams(dx_scaled,mt,BEAM_SCALE)
    return (mt, vt, βt .* β), dx_scaled
end 

Optimisers.init(o::AdamW, x::Beam{T}) where T = (Beam{T}(zeros(T,4)...), Beam{T}(zeros(T,4)...), T.(o.beta))
Optimisers.init(o::AdamW, x::CurvedBeam{T}) where T = (CurvedBeam{T}(zeros(T,4)...), CurvedBeam{T}(zeros(T,4)...), T.(o.beta))

function Optimisers.apply!(o::WeightDecay, state, x::Beam{T}, dx) where{T}
    λ = T(o.lambda)
    dx′ = dx + λ * x
  
    return state, Beam{T}(dx′...)
end

Optimisers.init(o::OptimiserChain, x::BeamElement) = map(opt -> Optimisers.init(opt, x), o.opts)

Optimisers.init(o::Optimisers.ClipNorm, x::BeamElement) = nothing

function Optimisers.apply!(o::Optimisers.ClipNorm, state, x::Beam{T}, dx) where T
  nrm = norm(dx, o.p)
  if o.throw && !isfinite(nrm)
    throw(DomainError("gradient has $(o.p)-norm $nrm, for array $(summary(x))"))
  end
  λ = T(min(o.omega / nrm, 1))
#   @show Beam{T}(dx * λ...)
  return state, λ * Beam{T}(dx) 
end

function Optimisers.apply!(o::AdamW, state, x::Beam{T}, dx) where T
    η, β, ϵ, λ = T(o.eta), T.(o.beta), T(o.epsilon), T(o.lambda)
    mt, vt, βt = state
  
    # standard Adam update with learning rate eta=1
    mt = combine(β[1],mt,dx)
    vt = combineabs2(β[2],vt,dx)
    dx′ = combine(η,βt,mt,vt,ϵ)
  
    # apply learning rate and weight decay
    if o.couple
      dx′′ =  η * (dx′ + λ * x)
    else
      dx′′ =  η * dx′ + λ * x
    end
  
    return (mt, vt, βt .* β), dx′′
end



function ode!(dU,u::AbstractVector{T},p::SciMLBase.NullParameters,s) where{T}
    @inbounds begin
        m,θ,x,y,fx,fy,κ = u
        s,c = sincos(θ)
        dU[1] = fx*s-fy*c   #dM
        dU[2] = m + κ               #dΘ
        dU[3] = c #* (1 + T[5]*h̃^2/(12)) #dx
        dU[4] = s #* (1 + T[6]*h̃^2/(12)) #dy
        dU[5] = zero(T)
        dU[6] = zero(T)
        dU[7] = zero(T)
    end 
    return dU
end 

function ode(t::AbstractVector{T},p::SciMLBase.NullParameters,s) where{T}
    @inbounds m,θ,x,y,fx,fy,κ = t
    [fx*sin(θ)-fy*cos(θ),   #dM
    m + κ,               #dΘ
    cos(θ), #* (1 + T[5]*h̃^2/(12)) #dx
    sin(θ), #* (1 + T[6]*h̃^2/(12)) #dy
    zero(T),
    zero(T),
    zero(T)]
end 

function jac(t::AbstractArray{T,N},p::SciMLBase.NullParameters,s) where{T,N} #jacobi 
    @inbounds begin
    m,θ,x,y,fx,fy,κ = t
        s,c = sincos(θ)
        dT = zeros(T,7,7)
        dT[1,2] = fx*c + fy*s
        dT[1,5] = s
        dT[1,6] = -c
        dT[2,1] = one(T)
        dT[2,7] = one(T)
        dT[3,2] = -s
        dT[4,2] = c
    end 
    return dT
end 

function jac!(dt,t::AbstractArray{T,N},p,s) where{T,N} #jacobi 
    
    @inbounds     m,θ,x,y,fx,fy,κ = t
    @inbounds     s,c = sincos(θ)
    @inbounds     dt[1,2] = fx*c + fy*s
    @inbounds     dt[1,5] = s
    @inbounds     dt[1,6] = -c
    @inbounds     dt[2,1] = one(T)#
    @inbounds     dt[2,7] = one(T)#
    @inbounds     dt[3,2] = -s
    @inbounds     dt[4,2] = c

    return dt
end

function vjp_beam!(Jv,λ::AbstractArray{T,N},u::AbstractVector,t) where{T,N} #vjp
    @inbounds begin
        δm,δθ,δx,δy,δfx,δfy,_= λ
        m,θ,x,y,fx,fy,_ = u
        s,c = sincos(θ)
        Jv[1] = -δθ 
        Jv[2] = -(δy * c - δx * s +  δm * (fx*c + fy*s) )
        Jv[5] = -(δm * s)
        Jv[6] =  δm  * c 
        Jv[7] = -δθ
    end  
    return nothing 
end

function vjp_curved_beam!(Jv,λ::AbstractArray{T,N},u::AbstractVector,t) where{T,N} #vjp
    @inbounds begin
        δm,δθ,δx,δy,δfx,δfy = λ
        m,θ,x,y,fx,fy,_ = u
        s,c = sincos(θ)
        Jv[1] = -δθ 
        Jv[2] = -(δy * c - δx * s +  δm * (fx*c + fy*s) )
        Jv[5] = -(δm * s)
        Jv[6] =  δm  * c 

    end  
    return nothing 
end

function vjp!(Jv,λ::AbstractArray{T,N},p::SciMLBase.NullParameters,t) where{T,N} #vjp
  
    #beam backpropergation
    λ_jvp = @view λ[1:7]
    Jv_jvp = @view Jv[1:7]
    #beam reverse integration
    du = @view Jv[8:14]
    u = @view λ[8:14]
    ode!(du,u,p,t)
    
    vjp_beam!(Jv_jvp,λ_jvp,u,t)
    return nothing 
end


#with splines
function ode!(du,u::AbstractArray{T,N},p::AbstractVector,t) where{T,N}
    
    @inbounds M,θ,x,y,Fx,Fy = u
    s,c = sincos(θ)
    du[1] = Fx * s - Fy *c
    du[2] = M + getspline(splinebasis,p,t)
    du[3] = c 
    du[4] = s 
    du[5] = zero(T)
    du[6] = zero(T)
end 

function vjp!(Jv::AbstractArray{TJ,NJ},λ::AbstractArray{T,N},p::AbstractVector,t) where{T,N,TJ,NJ} 
    #reverse integration of beam
    du = @view Jv[end-6:end]
    u = @view λ[end-6:end]
    ode!(du,u,p,t)
    #jvp reverse integration for adjoints
    λ_jvp = @view λ[1:6]
    Jv_jvp = @view Jv[1:6]
    vjp_curved_beam!(Jv_jvp,λ_jvp,u,t)
    #spline backpropagation     
    i,bs = evaluate_all(splinebasis,t)
    
    @inbounds Jv[6+i:-1:6+i+1-ORDER] .=  bs .* λ[2]  
    return nothing 
end




func = ODEFunction{true}(ode!, jac = jac!)

prob = ODEProblem(func,zeros(Float32,7),(0f0,1f0))

vjpfunc = ODEFunction{true}(vjp!)
vjpprob = ODEProblem(vjpfunc,ones(Float32,7),(1f0,0f0))