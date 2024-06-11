
function CRC.rrule(::typeof(residuals!),res::AbstractVector{T},node::Branch,beams,beams_end,y,ind) where{T}
    ind = residuals!(res,node,beams,beams_end,y,ind)
    
    function residuals!_branch_back(ȳ...)
        println(typeof(ȳ),"res")
        res_ = zeros(T,length(res))

        return (CRC.NoTangent(),ȳ,ȳ,ȳ,ȳ,ind)
    end 
    ind,residuals!_branch_back
end

function CRC.rrule(::typeof(residuals!),res::AbstractVector{T},node::Clamp,beams,beams_end,y,ind) where{T}
    ind = residuals!(res,node,beams,beams_end,y,ind)
    function residuals!_clamp_back(ȳ...)
        println(typeof(ȳ),"Clamp")
        return (CRC.NoTangent(),ȳ,ȳ,ȳ,ȳ,ind)
    end 
    ind,residuals!_clamp_back
end

function CRC.rrule(::typeof(initialize),str::Structure,parameters::Vector{T}) where {T}
    mat = initialize(str,parameters)
    function initialize_back(ȳ...)
        println(typeof(ȳ), "initialize")
        y = ones(T,size(parameters))
        return (CRC.NoTangent(),CRC.NoTangent(),y)
    end 
    mat,initialize_back
end 

function CRC.rrule(::typeof(reduction4loss),u,data,batch)
    out = reduction4loss(u,data,batch)
    function reduction4loss_back(ȳ)
        println(typeof(data),typeof(ȳ),"!reduction")
        grads = similar(data)
        for b in batch
            grads[b][1] .= ȳ[1][:,b]
            grads[b][2] .= ȳ[2][:,b]
        end 
        return CRC.NoTangent(),CRC.ZeroTangent(),grads,CRC.ZeroTangent()
    end 

    return out,reduction4loss_back
end 

