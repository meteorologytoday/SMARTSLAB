import NetCDF
const nc = NetCDF

mutable struct NCVar
    name::String
    attr::Dict
    coord::Array{String}
end

mutable struct NCDim
    name::String
    length::Int
end



mutable struct NCBox
    filename::String
    vars::Array{NCVar}
    coords::Dict{String, NCDim}


    function NCBox()
        new("", Array{NCVar, 1}(), Dict{String, NCDim}())
    end

end



#=
function readNCFile(filename::String) :: NCBox
    

end
=#





ncb = NCBox()
println(ncb.vars)
println(ncb.coords)




