
function getproperties(nodeorbeam::T) where{T}
    fnames = fieldnames(T)
    NamedTuple{fnames}(map(i->getfield(nodeorbeam,i),fnames))
end 

function getnames(nodes::Vararg{Boundary,N}) where{N}
    NamedTuple{ntuple(i->Symbol(:Node_,i),N)}(getproperties.(nodes))
end 

function getnames(beams::Vararg{Beam,N}) where{N}
    NamedTuple{ntuple(i->Symbol(:Beam_,i),N)}(getproperties.(beams))
end 

gettype(::Clamp) = :Clamp
gettype(::Branch) = :Branch

function prepare(args::Vararg{Union{Beam,Boundary},N}) where{N}
    beams = filter(x->isa(x,Beam),args)
    bounds = filter(x->isa(x,Boundary),args)
    (;Beams = getnames(beams...),Nodes = getnames(bounds...))        
end 