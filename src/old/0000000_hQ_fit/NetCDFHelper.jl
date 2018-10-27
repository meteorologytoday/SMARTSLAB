module NetCDFHelper

import NetCDF
export specialCopyNCFile

function specialCopyNCFile(
    fr_file :: String,
    to_file :: String,
    ncvars  :: Array{String,1};
    gatts   :: Union{Array{String, 1}, typeof(nothing)} = nothing,
    mode    :: UInt16 = NetCDF.NC_CLASSIC_MODEL
)

    dims_dict  = Dict{String, NetCDF.NcDim}() 

    new_ncvars = Array{NetCDF.NcVar, 1}()
    new_gatts  = Dict{String, Any}() 

    fr_fh = NetCDF.open(fr_file)

    # setup File with correct dimensions and vars
    for key in keys(fr_fh.dim)
        if key == "ncl3"
            continue
        end
        println("Found \"$key\", copy it.")
        _dim = fr_fh.vars[key]
        _dim = NetCDF.NcDim(key, NetCDF.readvar(_dim), _dim.atts)
        dims_dict[key] = _dim
    end

    for key in ncvars
        if(haskey(dims_dict, key))
            continue
        end
        if(haskey(fr_fh.vars, key))

            #println("Found \"$key\", copy it.")
            _var = fr_fh.vars[key]
            _dimlist = Array{NetCDF.NcDim, 1}()
            for i = 1 : length(_var.dim)
                push!(_dimlist, dims_dict[_var.dim[i].name])
            end

            _var = NetCDF.NcVar( #="$(key)xxx"=# key, _dimlist, atts=_var.atts)
            #println(key)
            #println(typeof(_var))

            #dump(_var)

            push!(new_ncvars, _var)
        else
            println("Cannot find \"$key\", ignore it.")
        end
    end

    if gatts != nothing
        for key in gatts
            new_gatts[key] = fr_fh.gatts[key]
        end
    end

    to_fh = NetCDF.create(to_file, new_ncvars, gatts=new_gatts, mode=mode)

    for new_ncvar in new_ncvars
        NetCDF.putvar(new_ncvar, NetCDF.readvar(fr_fh.vars[new_ncvar.name]))
    end

    NetCDF.close(to_fh)
    NetCDF.close(fr_fh)
end


#=
specialCopyNCFile(
    "mixed-layer-depth.nc",
    "test.nc",
    0,
    Array{String, 1}(),
    NetCDF.NC_CLASSIC_MODEL
)
=#






end
