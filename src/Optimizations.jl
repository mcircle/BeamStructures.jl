learningrate(it,pars = 200,warmups = 200) = √(pars) \ min(1/√(it),it / sqrt(warmups^3))

function changenode(bn,node::Symbol,nt::NamedTuple)
    ntmp = bn.Nodes[node] + nt  
    n = (;Beams = bn.Beams,Nodes = (;bn.Nodes..., node => ntmp)) 
end 

function solve_structure(str,bn,Δx,node)
    bnx = changenode(bn,node,(;x = Δx))
    inits = zeros(Float32,str,bn) # reset initial values for each column
    prob = NonlinearLeastSquaresProblem(str,inits,p = bnx)
    sol = solve(prob,TrustRegion(),sensealg = ZygoteAdjoint(),maxiters = 1000,reltol = 1e-2,abstol = 1e-2)
    sol,SciMLBase.successful_retcode(sol)
end 
     
function loss(loss,sols,res,fsoll,idxs)
    losspos = mean(abs2,res[:,2:end] .* [1,1,1])
    lossf =  1 * sum(abs2,res[:,1])
    lossn4 = 1 * sum(abs2,sols[6,2,idxs]) #force at node 4
    lossf4 = 1 * abs2(sum(sols[5,2,idxs]) .- fsoll) 
    loss + (losspos + lossf + lossn4 + lossf4)^2
end 

function loss_structure(str,parameters::Dict,graph,node;kwargs...)
    res = zero(y[:,:,1])
    loss = 0f0
    bn = parameters[:bn]
    inits = parameters[:inits]
    ids = findall(==(1),LowerTriangular(str.AdjMat))
    idxs = BeamStructures.findbeamsatnode(bn.Nodes[4],4,ids)[1]#[2,4]
    for ind in axes(graph,2)
        fxsoll =  graph[2,ind]
        bn_ = changenode(bn,(;x = graph[1,ind]))
        sols,bnx_ = str(inits[:,:,ind],bn_)
        # @show size(sols)
        res = BeamStructures.residuals!(res,str,sols,bnx_)
        losspos = mean(abs2,res[:,2:end] .* [1,1,1])
        lossf =  1 * sum(abs2,res[:,1])
        lossn4 = 1 * sum(abs2,sols[6,2,idxs]) #force at node 4
        lossf4 = 1 * abs2(sum(sols[5,2,idxs]) .- fxsoll) 
        loss = loss(loss,res,sols,fxsoll,bn,idxs) #* exp(-iterations/1000)) .* ((cos(iterations * π/100)^2) 
                      # (sin(iterations * π/100)^2) *
    end 
        sqrt(loss) + (4 * ((bn.Beams[1].h - bn.Beams[2].h)^2 + (bn.Beams[1].h - bn.Beams[3].h)^2) 
                    + 4 * ((bn.Beams[1].w - bn.Beams[2].w)^2 + (bn.Beams[1].w - bn.Beams[3].w)^2))
end 

function updatevalues(str::Structure,parameters::Dict,optstates::Dict;iters_done = 0)
    for (key,opt) in optstates
        ηup = learningrate(iters_done,50,1000)
        optstates[key] = Optimisers.adjust(opt,eta = ηup)
    end
    
    l,b = withgradient((x,y)->optimizestruct(str,x,y),parameters) 
    # b =  back(one(l)) #get gradients 
    for (key,opt) in enumerate(optstates)
        opt,bn = Optimisers.update(opt,parameters[key],b[1])
    end
    # @assert !isnan(bn.Beams[1].l)
    # opst2,sol = Optimisers.update(opst2,sol,b[2])
    l,optstaates,parameters
end 