
abstract type Boundary{T} end
struct Clamp{A<:Real} <:Boundary{A}
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    mz::A
    function Clamp(x::T,y::T,ϕ::T,fx::T,fy::T,mz::T) where{T}
        new{T}(x,y,ϕ,fx,fy,mz)
    end 
    Clamp{T}(x,y,ϕ,fx,fy,mz) where{T} = new{T}(x,y,ϕ,fx,fy,mz)
    
end 
Clamp(x::AbstractVector) = Clamp(x...)
(cl::Clamp{T})(::Type{T}) where{T} = cl
(cl::Clamp{X})(::Type{T}) where{X,T} = Clamp{T}(cl[1:6]...)

gettype(::Clamp{T}) where{T} = Clamp
struct Branch{A<:Real} <:Boundary{A}
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    mz::A
    function Branch(x::T,y::T,ϕ::T,fx::T,fy::T,mz::T) where{T}
        new{T}(x,y,ϕ,fx,fy,mz)
    end 
    Branch{T}(x,y,ϕ,fx,fy,mz) where{T} = new{T}(x,y,ϕ,fx,fy,mz)

end
Branch(x::AbstractVector) = Branch(x...)
(cl::Branch{T})(::Type{T}) where{T} = cl
(cl::Branch{X})(::Type{T}) where{X,T} = Branch{T}(cl[1:6]...)
gettype(::Branch{T}) where{T} = Branch


struct Free{A<:Real} <:Boundary{A}
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    mz::A
    function Free(x::T,y::T,ϕ::T,fx::T,fy::T,mz::T) where{T}
        new{T}(x,y,ϕ,fx,fy,mz)
    end 
end
Free(x::AbstractVector) = Free(x...)
Free{T}(args...) where{T} = Free(T.(args)...)
gettype(::Free{T}) where{T} = Free
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
    function ExtForces(x::T,y::T,ϕ::T,fx::T,fy::T,mz::T) where{T}
        new{T}(x,y,ϕ,fx,fy,mz)
    end 
end
ExtForces(x::AbstractVector) = ExtForces(x...)
ExtForces{T}(args...) where{T} = ExtForces(T.(args)...)
gettype(::ExtForces{T}) where{T} = ExtForces
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
    f::Function
    function Slider(x::T,y::T,ϕ::T,fx::T,fy::T,mz::T) where{T}
        new{T}(x,y,ϕ,fx,fy,mz)
    end 
end
Slider(x::AbstractVector) = Slider(x...)
Slider{T}(args...) where{T} = Slider(T.(args)...)
gettype(::Slider{T}) where{T} = Slider
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
    function CompliantClamp(x::T,y::T,ϕ::T,fx::T,fy::T,mz::T) where{T}
        new{T}(x,y,ϕ,fx,fy,mz)
    end 
end
CompliantClamp(x::AbstractVector) = CompliantClamp(x...)
CompliantClamp{T}(args...) where{T} = CompliantClamp(T.(args)...)
gettype(::CompliantClamp{T}) where{T} = CompliantClamp

function (::Type{B} )(x::NamedTuple) where{B<:Boundary{<:Real}}
    B(map(k->getfield(x,k),keys(x))...)
end 

function change_node(node::Boundary;kwargs...)
    for (field,value) in kwargs
        node = Setfield.@set node.$field = value 
    end 
    node
end  

Base.zero(::B) where{T,B<:Boundary{T}} = B(zeros(T,6)...)

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
Optimisers.trainable(c::Clamp) = (;x = c.x,y = c.y,ϕ = c.ϕ)
Optimisers.trainable(c::Branch) = (;x = c.x,y = c.y,ϕ = c.ϕ)
Optimisers.subtract!(a::T,b::T) where{T<:Boundary{TT}} where{TT} = a - T(merge(Optimisers.mapvalue(_->zero(TT),Optimisers.functor(b)[1]),Optimisers.trainable(b)))
Optimisers.init(o::OptimiserChain, x::Boundary) = map(opt -> Optimisers.init(opt, x), o.opts)
Optimisers._trainable(b::T,fr) where{T<:Boundary} =T(merge(Optimisers.mapvalue(_ -> nothing, Optimisers.functor(b)[1]), Optimisers.trainable(b)))



@generated function combine(β,mt::B,dx) where{T,B<:Boundary{T}}
    fields = fieldnames(mt)
    exprs = [:(isnothing(getproperty(dx, $(QuoteNode(f)))) ? zero(T) : β * getproperty(mt, $(QuoteNode(f))) + (1-β) * getproperty(dx, $(QuoteNode(f)))) for f in fields]
    return quote
        $(Expr(:call, :B, exprs...))
    end
end 

@generated function combineabs2(β,mt::B,dx) where{T,B<:Boundary{T}}
    fields = fieldnames(mt)
    exprs = [:(isnothing(getproperty(dx, $(QuoteNode(f)))) ? zero(T) : β * getproperty(mt, $(QuoteNode(f))) + (1-β) * abs2(getproperty(dx, $(QuoteNode(f))))) for f in fields]
    return quote
        $(Expr(:call, :B, exprs...))
    end
end 

@generated function combine(η,βt,mt::B,vt::B,ϵ) where{B<:Boundary}
    fields = fieldnames(mt)
    exprs = [:(getproperty(mt, $(QuoteNode(f))) / (1-βt[1]) / (sqrt(getproperty(vt, $(QuoteNode(f))) / (1 -βt[2]))+ ϵ)* η) for f in fields]
    return quote
        $(Expr(:call,B, exprs...))
    end
end

function Optimisers.apply!(o::Adam,state,b::BT,dx) where{BT<:Boundary{T}} where{T}
    η, β, ϵ = T(o.eta), T.(o.beta), T(o.epsilon)
    mt, vt, βt = state
    
    mt =combine(β[1],mt,dx) # β[1] * mt + (1 - β[1]) * dx
    vt = combineabs2(β[2],vt,dx) # β[2] * vt + (1 - β[2]) * abs2(dx)
    dx′ = combine(η,βt,mt,vt,ϵ) # mt / (1 - βt[1]) / (sqrt(vt / (1 - βt[2])) + ϵ) * η
  
    return (mt, vt, βt .* β), dx′
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