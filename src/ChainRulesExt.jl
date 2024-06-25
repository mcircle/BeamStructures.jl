
function CRC.rrule(::typeof(residuals!),res::AbstractVector{T},node::Branch,beams,beams_end,y,ind) where{T}
    ind = residuals!(res,node,beams,beams_end,y,ind)
    
    function residuals!_branch_back(ȳ)
        println(typeof(ȳ),"res")
        @assert 5 == 1
        res_ = zeros(T,length(res))

        return (NoTangent(),ȳ,ȳ,ȳ,ȳ,ind)
    end 
    ind,residuals!_branch_back
end

function CRC.rrule(::typeof(residuals!),res::AbstractVector{T},node::Clamp,beams,beams_end,y,ind) where{T}
    ind = residuals!(res,node,beams,beams_end,y,ind)
    function residuals!_clamp_back(ȳ)
        println(typeof(ȳ),"Clamp")
        return (NoTangent(),ȳ,ȳ,ȳ,ȳ,ind)
    end 
    ind,residuals!_clamp_back
end

function CRC.rrule(::typeof(reduction4loss),u,data,batch)
    out = reduction4loss(u,data,batch)
    function reduction4loss_back(ȳ)
        st,en = Array(ȳ[1])
        grads = Vector{Tuple{Vector,Vector}}()
        for (s,e) in zip(eachcol(st),eachcol(en))
            push!(grads,(s,e))
        end 
        return NoTangent(),ZeroTangent(),grads,ȳ[2]
    end 

    return out,reduction4loss_back
end 



