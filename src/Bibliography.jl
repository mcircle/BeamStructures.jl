  
const DBDIR = pwd() 
@inline loadstructure(file,name::String) = file[name] 
@inline loadstructure(file,name::Val) = loadstructure(file,string(name)) 

@inline function add_structure(file,name::String,str::Structure) 
    file[name] = str
end 

function loadstructures(file::String)     
    if ispath(file)
        tmp = file[1:end-5] * "tmp.jld2"
        cp(tmp,file,force = true)
        jldfile = jldopen(file,"r+")
    else 
        jldfile = jldopen(file,"w")
    end 
    jldfile
end 

savestructures(file::JLD2.JLDFile) = close(file)

function loaddatabase(dir = nothing)
    file = isnothing(dir) ? file =  pwd() * "./DataBase/Database.xlsx" : dir
    identity.(DataFrame(XLSX.readtable(file,"Structures","A:G",first_row = 1,header = true)))
end 