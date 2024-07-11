
abstract type Boundary{T} end

struct Clamp{A<:Real} <:Boundary{A}
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    mz::A
end 

struct Branch{A<:Real} <:Boundary{A}
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    mz::A
end

# function Branch(x::X,y::Y,ϕ::P,fx::FX,fy::FY,mz::MZ) where {X,Y,P,FX,FY,MZ}
#     p = promote(x,y,ϕ,fx,fy,mz)
#     Branch(p...)
# end  

struct Free{A<:Real} <:Boundary{A}
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    mz::A
end

# function Free(x::X,y::Y,ϕ::P,fx::FX,fy::FY,mz::MZ) where {X,Y,P,FX,FY,MZ}
#     p = promote(x,y,ϕ,fx,fy,mz)
#     Free(p...)
# end  

struct ExtForces{A<:Real} <: Boundary{A}
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    mz::A
end 

# function ExtForces(x::X,y::Y,ϕ::P,fx::FX,fy::FY,mz::MZ) where {X,Y,P,FX,FY,MZ}
#     p = promote(x,y,ϕ,fx,fy,mz)
#     ExtForces(p...)
# end  

struct Slider{A<:Real} <:Boundary{A}
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    mz::A
end

# function Slider(x::X,y::Y,ϕ::P,fx::FX,fy::FY,mz::MZ) where {X,Y,P,FX,FY,MZ}
#     p = promote(x,y,ϕ,fx,fy,mz)
#     Slider(p...)
# end  

Base.length(b::Boundary) = 6
Base.getindex(b::Boundary,idx::AbstractVector) = map(x->getfield(b,x),fieldnames(Clamp)[idx])
Base.getindex(b::Boundary,idx::Int) = getfield(b,fieldnames(Clamp)[idx])
Base.iterate(b::Boundary,i::Int = 1) = i > 6 ? nothing : (getfield(b,i),i+1)
Base.IteratorSize(b::T) where{T<:Boundary} = Base.HasLength()

struct CompliantClamp{A} <:Boundary{A}
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    c::A
end


function (::Type{B} )(x::NamedTuple) where{B<:Boundary{<:Real}}
    B(map(k->getfield(x,k),keys(x))...)
end 

function change_node(node::Boundary;kwargs...)
    for (field,value) in kwargs
        node = Setfield.@set node.$field = value 
    end 
    node
end  

# function CompliantClamp(x::X,y::Y,ϕ::P,fx::FX,fy::FY,c::C) where {X,Y,P,FX,FY,C}
#     p = promote(x,y,ϕ,fx,fy,c)
#     CompliantClamp(p...)
# end  

Base.:*(a::Real,b::T) where{T<:Boundary} = T(a .* b...)
Base.:*(b::T,a::Real) where{T<:Boundary} = T(a .* b...)
Base.:*(a::T,b::T) where{T<:Boundary} = T(a .* b...)
Base.:+(a::T,b::T) where{T<:Boundary} = T(a .+ b...)
Base.:-(a::T,b::T) where{T<:Boundary}= T(a .- b...)
Base.:-(a::T,b::Real) where{T<:Boundary} = T(a .- b...)
Base.:-(a::Real,b::T) where{T<:Boundary} = T(a .- b...)

Base.:abs2(b::T) where{T<:Boundary} = b*b
Base.:/(a::Real,b::T) where{T<:Boundary} = T(a ./ b...)
Base.:/(b::T,a::Real)where{T<:Boundary} = T(b ./ a...)

Base.:+(a::T,b::NamedTuple) where{T<:Boundary} = T(map(x->isnothing(getfield(b,x)) ? getfield(a,x) : getfield(a,x) + getfield(b,x),keys(b))...)

Optimisers.functor(x::T) where{T<:Boundary} = (NamedTuple{fieldnames(T)}(x[1:6]),T)
Optimisers.init(o::Adam, x::B) where{B<:Boundary{T}} where{T}  = (B(zeros(T,6)...), B(zeros(T,6)...), T.(o.beta))
Optimisers.init(o::WeightDecay, x::Boundary) = nothing
Optimisers.isnumeric(x::T) where{T<:Boundary} = true
Optimisers.trainable(c::Clamp) = (;x = c.x,y = c.y)
Optimisers.trainable(c::Branch) = (;x = c.x,y = c.y)
Optimisers.subtract!(a::T,b::T) where{T<:Boundary{TT}} where{TT} = a - T(merge(Optimisers.mapvalue(_->zero(TT),Optimisers.functor(b)[1]),Optimisers.trainable(b)))
Optimisers.init(o::OptimiserChain, x::Boundary) = map(opt -> Optimisers.init(opt, x), o.opts)
Optimisers._trainable(b::T,fr) where{T<:Boundary} =T(merge(Optimisers.mapvalue(_ -> nothing, Optimisers.functor(b)[1]), Optimisers.trainable(b)))
function Optimisers.apply!(o::Adam,state,b::BT,dx) where{BT<:Boundary{T}} where{T}
    η, β, ϵ = T(o.eta), T.(o.beta), T(o.epsilon)
    mt, vt, βt = state
    # @show mt vt dx
    Optimisers.@.. mt = β[1] * mt + (1 - β[1]) * dx
    Optimisers.@.. vt = β[2] * vt + (1 - β[2]) * abs2(dx)
    dx′ = Optimisers.@lazy mt / (1 - βt[1]) / (sqrt(vt / (1 - βt[2])) + ϵ) * η
  
    return (BT(mt...), BT(vt...), βt .* β), BT(dx′...)
end 

function Optimisers.apply!(o::WeightDecay, state, x::BT, dx) where{BT<:Boundary{T}} where{T}
    λ = T(o.lambda)
    dx′ = dx + λ * x
  
    return state, BT(dx′...)
end

Optimisers.init(o::Optimisers.ClipNorm, x::Boundary) = nothing

function Optimisers.apply!(o::Optimisers.ClipNorm, state, x::BT, dx) where{BT<:Boundary{T}} where T
  nrm = norm(dx, o.p)
  if o.throw && !isfinite(nrm)
    throw(DomainError("gradient has $(o.p)-norm $nrm, for array $(summary(x))"))
  end
  λ = T(min(o.omega / nrm, 1))

  return state,  BT(dx * λ...)
end