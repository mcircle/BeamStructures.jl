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

function multiply_knots(t::T,N,p,i) where{T} 
    iszero(N)  &&  return zero(T)
    abs(BREAKS[i] - t)/(BREAKS[i+p] - BREAKS[i]) * N
end

function basis(t::T,::Val{0},::Val{I}) where{T,I} 
    I < 5 || return zero(T)
    ifelse(BREAKS[I]≤ t < BREAKS[I+1], one(T), zero(T))
end 

function basis(t::T,::Val{P},::Val{I})  where{T,P<:Int,I<:Int}  
    N1 = basis(t,Val(P-1),Val(I))
    N2 = basis(t,P-1,I+1)

    multiply_knots(t,N1,P-1,I) + multiply_knots(t,N2,P-1,I+1)

end 