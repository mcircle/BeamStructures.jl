
abstract type Boundary end

struct Clamp{A<:Number} <:Boundary
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    mz::A
end 


# function Clamp(x::X,y::Y,ϕ::P,fx::FX,fy::FY,mz::MZ) where {X,Y,P,FX,FY,MZ}
#     p = promote(x,y,ϕ,fx,fy,mz)
#     Clamp(p...)
# end  

struct Branch{A<:Number} <:Boundary
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

struct Free{A<:Number} <:Boundary
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

struct ExtForces{A<:Number} <: Boundary
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

struct Slider{A<:Number} <:Boundary
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

struct CompliantClamp{A<:Number} <:Boundary
    x::A
    y::A
    ϕ::A
    fx::A
    fy::A
    c::A
end

# function CompliantClamp(x::X,y::Y,ϕ::P,fx::FX,fy::FY,c::C) where {X,Y,P,FX,FY,C}
#     p = promote(x,y,ϕ,fx,fy,c)
#     CompliantClamp(p...)
# end  


