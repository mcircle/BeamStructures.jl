
abstract type Boundary end
struct Clamp <: Boundary end
struct Branch <: Boundary end
struct Free <: Boundary end
struct Slider <: Boundary end
struct Input <: Boundary end

struct Clamp_{T} <:Boundary
    x::T
    y::T
    ϕ::T
    fx::T
    fy::T
    mz::T
end 
function Clamp_(x::X,y::Y,ϕ::P,fx::FX,fy::FY,mz::MZ) where {X,Y,P,FX,FY,MZ}
    p = promote(x,y,ϕ,fx,fy,mz)
    Clamp_(p...)
end  

struct Branch_{T} <:Boundary
    x::T
    y::T
    ϕ::T
    fx::T
    fy::T
    mz::T
end

function Branch_(x,y,ϕ,fx,fy,mz)
    p = promote(x,y,ϕ,fx,fy,mz)
    Branch_(p...)
end  

function Branch_{T}(x,y,ϕ,fx,fy,mz) where {T<:Number}
    Branch_(T(x),T(y),T(ϕ),T(fx),T(fy),T(mz))
end 

struct Free_{T} <:Boundary
    x::T
    y::T
    ϕ::T
    fx::T
    fy::T
    mz::T
end

function Free_(x,y,ϕ,fx,fy,mz)
    p = promote(x,y,ϕ,fx,fy,mz)
    Free_(p...)
end  

function Free_{T}(x,y,ϕ,fx,fy,mz) where {T<:Number}
    Free_(T(x),T(y),T(ϕ),T(fx),T(fy),T(mz))
end 

struct Slider_{T} <:Boundary
    x::T
    y::T
    ϕ::T
    fx::T
    fy::T
    mz::T
end

function Slider_(x,y,ϕ,fx,fy,mz)
    p = promote(x,y,ϕ,fx,fy,mz)
    Slider_(p...)
end  

function Slider_{T}(x,y,ϕ,fx,fy,mz) where {T<:Number}
    Slider_(T(x),T(y),T(ϕ),T(fx),T(fy),T(mz))
end 

struct CompliantClamp{T} <:Boundary
    x::T
    y::T
    ϕ::T
    fx::T
    fy::T
    mz::T
    c::T
end

function CompliantClamp(x,y,ϕ,fx,fy,mz,c)
    p = promote(x,y,ϕ,fx,fy,mz,c)
    CompliantClamp(p...)
end  

function CompliantClamp{T}(x,y,ϕ,fx,fy,mz,c) where {T<:Number}
    CompliantClamp(T(x),T(y),T(ϕ),T(fx),T(fy),T(mz),T(c))
end 

struct Beam{T} 
    l::T
    h::T
    w::T
    κ0::T
    E::T
    θs::T
    θe::T
end 

function Beam(l,h,w,κ0;E = 2e12,θs = 0.,θe = θs + l*κ0)
    p = promote(l,h,w,κ0,E,θs,θe)
    Beam(p...)
end 

