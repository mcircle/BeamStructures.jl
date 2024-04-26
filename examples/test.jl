beam1 = Beam(100e-3,10.e-3,2.e-3,0.)
node1 = Clamp(0.,0.,0.,0.,0.,0.)
node2 = Free(0.,0.,0.,0.,75.,2.)
node2[[1,4,5]]

adj = [0 1 0; 
1 0 1;
0 1 0]
adj = [0 1; 
       1 0]


con = Connections(adj,node1,node2,node1,node2,node2)

str = Structure(con,beam1)

inits = initialize(str,rand(Float64,3))
@run inits = initialize(str,rand(Float64,3))

sol = str(rand(Float64,3))
length(sol)
size(sol)
size(sol[1])

using Plots

plot(sol,idxs = (3,4))

res = rand(Float64,3)
inp = rand(Float64,3)

@run residuals!(res,str,sol)
residuals!(res,str,sol)
res


adj =Int.([ 0 1 0 0 0; 
       1 0 1 0 1;
       0 1 0 1 0;
       0 0 1 0 0;
       0 1 0 0 0
])

@time beam_matrix(adj)

adj
con = Connections(adj,node1,node2,node1,node2,node2)
nf = rand(Float64,4,5)
bf = rand(Float64,4,4)
nf
bf
view_nf(ind) = @view nf[:,ind]
@time size.(view_nf.(gf.Infos.nbrs))
nbrs = [1]

nf_ = nf[:,nbrs]

gf = GraphLayer([4,4],8,adj)
for nbrs in gf.Infos.nbrs
    @show    size( view(nf,:,nbrs))
end 

θ,st = Lux.setup(Random.default_rng(),gf)
θ = ComponentArray(θ)
gf((nf,bf),θ,st)
typeof(st)
((n1,n2),st),back = Zygote.pullback(x->gf((x,bf),θ,st),nf)
gs = back(((ones(eltype(n1),size(n1)...),nothing),nothing))[1]


@run NamedTuple(zip(Tuple([]),Tuple([]))...)

Adj_norm(adj)
test = adj[1,:]
test[1+1:end] .*= -1
test

findall(==(1),LowerTriangular(adj)) 

_keys = (:a, :b, :c); _values = (1, 2, 3);
t = (; zip(_keys, _values)...)
typeof(t)
hcat(rand(Float64,5,1),rand(Float64,5,1))

model = Chain(Dense(5,100),Dense(100,5))
typeof(model)
deq = Chain(Dense(1,5),DEQ(model),Dense(5,1))
deq = DEQ(model)
X = reshape(Float32[1;2;3;4;5;6;7;8;9;10], 1, :) 
Y = 2 .* X
loss(x,y,ps,st) = sum(abs2, y .- deq(x,ps,st)[1])

θ,st = Lux.setup(Random.default_rng(), deq)
st
θ = ComponentArray(θ)
θ.layer_2 .*=0.0001
opt = Optimisers.Adam(0.0005,(0.95,0.99))
optps = Optimisers.setup(opt,θ)
θ_ = deepcopy(θ)
for i in 1:1000
       (l),back = Zygote.pullback(theta->loss(X,Y,theta,st),θ)
       gs = back((one(l),))[1]
       optps,θ = Optimisers.update!(optps,θ,gs)
       println("Loss: $(l)")
end 


zstar,st = deq(X,θ,st)
(l),back = Zygote.pullback((x,theta)->loss(x,Y,theta,st),X,θ)
gra_z =  back(1 * one(l))[1]
l_ = loss(X .+ 1,Y,θ,st)
