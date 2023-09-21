beam1 = Beam(100e-3,10.e-3,2.e-3,0.)
node1 = Clamp(0.,0.,0.,0.,0.,0.)
node2 = Free(0.,0.,0.,0.,75.,2.)
node2[[1,4,5]]

adj = [0 1 0; 
1 0 1;
0 1 0]
adj = [0 1; 
       1 0]


con = Connections(adj,node1,node2)

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