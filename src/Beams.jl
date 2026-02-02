abstract type BeamElement{T} end

struct Beam{T<:Real} <:BeamElement{T}
    l::T
    h::T
    w::T
    κ0::T
    E::T
    θs::T
    θe::T
    # function Beam{T}(l,h,w,κ0,E,θs,θe) where{T}
    #     new{T}(l,h,w,κ0,E,θs,θe)
    # end 
end 

using BSplineKit

const KNOTLENGTH = 5
const ORDER = 4
BREAKS =range(0f0,1f0,length=KNOTLENGTH)
order = BSplineOrder(ORDER)
pos_knots = BSplineKit.SplineInterpolations.make_knots(BREAKS, order, nothing)
splinebasis = BSplineBasis(order,pos_knots;augment= Val(false))

function getspline(basis,c,t)
    i,bs = evaluate_all(basis,t)
    return sum(c[i:-1:i-ORDER+1] .* bs)
end 

function splineode!(du::AbstractVector,u,p::AbstractVector,t)
    du .= getspline(splinebasis,p,t)
end 

splinefunc = ODEFunction{true}(splineode!)
splineprob = ODEProblem(splinefunc,zero(Float32),(0f0,1f0))

get_θe(θs,κ0) = solve(splineprob,Tsit5(),u0 = [θs], p=κ0,reltol=1e-6,abstol=1e-6).u[end]

struct CurvedBeam{T<:Real} <:BeamElement{T}
    l::T
    h::T
    w::T
    κ0::AbstractVector{T}
    E::T
    θs::T
    θe::T
    function CurvedBeam{T}(l::T,h::T,w::T,κ0::AbstractVector{T},E::T,θs::T,θe::T) where{T}
        new{T}(l,h,w,κ0,E,θs,θe)
    end 
end

function CurvedBeam(l,h,w,κ0;E = 2.1f5,θs = 0)
    p = promote(l,h,w,E,θs)
    θe = get_θe(p[5],κ0)
    CurvedBeam{eltype(p)}(p[1],p[2],p[3],κ0;E = p[4],θs = p[5],θe = θe)
end
function CurvedBeam{T}(l,h,w,κ0;E = 2.1f5,θs = 0) where{T}

    p = T.([l,h,w,E,θs])
    k0 = T.(κ0)
    θe = get_θe(p[5],k0)
    CurvedBeam{T}(p[1],p[2],p[3],k0,p[4],p[5],only(θe))
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

Base.length(b::BeamElement) = 7
Base.lastindex(b::BeamElement) = 7
function Base.getproperty(b::BeamElement,n::Symbol)
    val = getfield(b,n)
end 


Base.getindex(b::B,idx::Int) where{B<:BeamElement} = getproperty(b,fieldnames(B)[idx])
Base.getindex(b::B,idx::AbstractVector) where{B<:BeamElement} = map(x->getindex(b,x),idx)
Base.iterate(b::B,i::Int = 1) where{B<:BeamElement} = i > 7 ? nothing : (getfield(b,i),i+1)
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

Optimisers.init(o::OptimiserChain, x::Beam) = map(opt -> Optimisers.init(opt, x), o.opts)

Optimisers.init(o::Optimisers.ClipNorm, x::Beam) = nothing

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
    @inbounds m,θ,x,y,fx,fy,κ = u
    dU[1] = fx*sin(θ)-fy*cos(θ)   #dM
    dU[2] = m + κ               #dΘ
    dU[3] = cos(θ) #* (1 + T[5]*h̃^2/(12)) #dx
    dU[4] = sin(θ) #* (1 + T[6]*h̃^2/(12)) #dy
    dU[5] = zero(T)
    dU[6] = zero(T)
    dU[7] = zero(T)
    return dU
end 

function ode(t::AbstractVector{T},p,s) where{T}
    @inbounds m,θ,x,y,fx,fy,κ = t
    [fx*sin(θ)-fy*cos(θ),   #dM
    m + κ,               #dΘ
    cos(θ), #* (1 + T[5]*h̃^2/(12)) #dx
    sin(θ), #* (1 + T[6]*h̃^2/(12)) #dy
    zero(T),
    zero(T),
    zero(T)]
end 

function jac(t::AbstractArray{T,N},p,s) where{T,N} #jacobi 
    @inbounds m,θ,x,y,fx,fy,κ = t
    dT = zeros(T,7,7)
    dT[1,2] = fx*cos(θ) + fy*sin(θ)
    dT[1,5] = sin(θ)
    dT[1,6] = -cos(θ)
    dT[2,1] = one(T)
    dT[2,7] = one(T)
    dT[3,2] = -sin(θ)
    dT[4,2] = cos(θ)
    return dT
end 

function jac!(dt,t::AbstractArray{T,N},p::SciMLBase.NullParameters,s) where{T,N} #jacobi 
    @inbounds m,θ,x,y,fx,fy,κ = t
    dt[1,2] = fx*cos(θ) + fy*sin(θ)
    dt[1,5] = sin(θ)
    dt[1,6] = -cos(θ)
    dt[2,1] = one(T)#
    dt[2,7] = one(T)#
    dt[3,2] = -sin(θ)
    dt[4,2] = cos(θ)
    return dt
end

function vjp!(Jv,λ::AbstractArray{T,N},u::AbstractVector,t) where{T,N} #vjp
    @inbounds δm,δθ,δx,δy,δfx,δfy= λ
    @inbounds m,θ,x,y,fx,fy = u
    # @show "test curve inner"
    Jv[1] = -δθ 
    Jv[2] = -(δy * cos(θ) - δx * sin(θ) +  δm * (fx*cos(θ) + fy*sin(θ)) )
    Jv[5] = -(δm * sin(θ))
    Jv[6] =  δm  * cos(θ) 
    Jv[7] = zero(T) #-δθ 
    return nothing 
end

function vjp!(Jv,λ::AbstractArray{T,N},u::ODESolution,t,::SciMLBase.NullParameters) where{T,N} #vjp
    fill!(Jv,zero(eltype(Jv)))
    # @show "test curve"
    vjp!(Jv,λ,u(t),t)
    return nothing 
end



#with splines
function ode!(du,u::AbstractArray{T,N},p::AbstractVector,t) where{T,N}
    @inbounds M,θ,x,y,Fx,Fy = u
    du[1] = Fx * sin(θ) - Fy *cos(θ)
    du[2] = M + getspline(splinebasis,p,t)
    du[3] = cos(θ) 
    du[4] = sin(θ) 
    du[5] = zero(T)
    du[6] = zero(T)
    # du[7] = zero(T)
end 

function vjp!(Jv,λ::AbstractArray{T,N},u::ODESolution,t,::AbstractVector) where{T,N} 
    fill!(Jv,zero(eltype(Jv)))
    # @show "test with spline"
    vjp!(Jv,λ[1:6],u(t),t)
    i,bs = evaluate_all(splinebasis,t)
    Jv[6+i:-1:6+i+1-ORDER] .-=  bs .* λ[2]
    return nothing 
end


function vjp!(Jv,λ::AbstractArray{T,N},u::ODESolution{T,N2,uType},t) where{T,N,N2,uType}
    vjp!(Jv,λ,u,t,u.prob.p)
    # @show "test"
    return nothing 
end

func = ODEFunction{true}(ode!, jac = jac!)

prob = ODEProblem(func,zeros(Float32,7),(0f0,1f0))

vjpfunc = ODEFunction(vjp!)
vjpprob = ODEProblem(vjpfunc,ones(Float32,7),(1f0,0f0))
