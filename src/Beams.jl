struct Beam{T} 
    l::T
    h::T
    w::T
    κ0::T
    E::T
    θs::T
    θe::T
end 

function Beam(l,h,w,κ0;E = 2.1f5,θs = 0f0,θe = θs + l*κ0)
    p = promote(l,h,w,κ0,E,θs,θe)
    Beam(p...)
end 
function Beam(T,l,h,w,κ0;E = 2.1f5,θs = 0f0,θe = θs + l*κ0)
    Beam{T}(l,h,w,κ0,E,θs,θe)
end 

function Base.show(io::IO,beam::Beam)
    return println(io, "Beam with Length: $(beam.l),width: $(beam.w), height: $(beam.h), curvature: $(beam.κ0) and E: $(beam.E)")
end 

Base.length(b::Beam) = 7
Base.getindex(b::Beam,idx::AbstractVector) = map(x->getfield(b,x),fieldnames(Beam)[idx])
Base.getindex(b::Beam,idx::Int) = getfield(b,fieldnames(Beam)[idx])

function ode!(dT,t::AbstractVector{T},p,s) where{T}
    m,θ,x,y,fx,fy,κ = t
    dT[1] = fx*sin(θ)-fy*cos(θ)   #dM
    dT[2] = m + κ               #dΘ
    dT[3] = cos(θ) #* (1 + T[5]*h̃^2/(12)) #dx
    dT[4] = sin(θ) #* (1 + T[6]*h̃^2/(12)) #dy
    dT[5] = zero(T)
    dT[6] = zero(T)
    dT[7] = zero(T)
    return dT
end 

function jac!(t::AbstractArray{T,N},p,s) where{T,N} #jacobi 
    m,θ,x,y,fx,fy,κ = t
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

function jac!(dt,t::AbstractArray{T,N},p,s) where{T,N} #jacobi 
    m,θ,x,y,fx,fy,κ = t
    dt[1,2] = fx*cos(θ) + fy*sin(θ)
    dt[1,5] = sin(θ)
    dt[1,6] = -cos(θ)
    dt[2,1] = one(T)
    dt[2,7] = one(T)
    dt[3,2] = -sin(θ)
    dt[4,2] = cos(θ)
    return dt
end

func = ODEFunction{true}(ode!, jac = jac!)

prob = ODEProblem(func,zeros(Float64,7),(0.,1.))

