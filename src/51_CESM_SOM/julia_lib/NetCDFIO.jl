module NetCDFIO

using NCDatasets

float = Union{Float64, Float32}

mutable struct MapInfo{T <: float}
    nx :: Integer
    ny :: Integer
    lsize :: Integer
 
    xc :: Array{T, 2}
    yc :: Array{T, 2}
    
    mask :: Array{T, 2}
    area :: Array{T, 2}
    frac :: Array{T, 2}

    missing_value :: T

    function MapInfo{T}(
        filename::String;
        missing_value::T = 1e20
    ) where T <: float
    
        ds = Dataset(filename, "r")
        _mask = ds["mask"][:, :]
        _area = ds["area"][:, :]
        _frac = ds["frac"][:, :]
        _xc  = ds["xc"][:, :]
        _yc  = ds["yc"][:, :]
        _nx  = ds.dim["ni"]
        _ny  = ds.dim["nj"]
        close(ds)

        return new{T}(
            _nx, _ny, _nx * _ny,
            _xc, _yc,
            _mask, _area, _frac,
            missing_value,
        )

    end

end

function createNCFile(
    mi::MapInfo{T},
    filename::String;
) where T <: float

    ds = Dataset(filename, "c")

    defDim(ds, "time", Inf)
    defDim(ds, "ni", mi.nx)
    defDim(ds, "nj", mi.ny)
    defDim(ds, "n", mi.lsize)

    close(ds)


    write2NCFile(mi, filename, "xc", mi.xc;     time_exists=false, missing_value=mi.missing_value)
    write2NCFile(mi, filename, "yc", mi.yc;     time_exists=false, missing_value=mi.missing_value)
    write2NCFile(mi, filename, "mask", mi.mask; time_exists=false, missing_value=mi.missing_value)
    write2NCFile(mi, filename, "frac", mi.frac; time_exists=false, missing_value=mi.missing_value)
    write2NCFile(mi, filename, "area", mi.area; time_exists=false, missing_value=mi.missing_value)

end

function write2NCFile(
    mi          :: MapInfo{T},
    filename    :: String,
    varname     :: String,
    var         :: Array{T};
    time        :: Union{Nothing, UnitRange, Integer} = nothing,
    time_exists :: Bool = true,
    missing_value :: Union{T, Nothing} = nothing,
) where T <: float

    local ds_var
    ds = Dataset(filename, "a")

    if ! ( varname in keys(ds) )
        ds_var = defVar(ds, varname, T, (time_exists) ? ("ni", "nj", "time") : ("ni", "nj") )
        
        if missing_value != nothing
            ds_var.attrib["_FillValue"] = missing_value 
        end
    else
        ds_var = ds[varname]
    end

    if time_exists

        # If this is a static 2d data, then
        # make the third dimension as time
        if length(size(var)) == 2
            var = view(var, :, :, 1)
        end

        append_time_len = size(var, 3)
       
        if time == nothing # append at last

            old_time_len = size(ds_var, 3)

            _beg = 1 + old_time_len
            _end = _beg + append_time_len - 1
            ds_var[:, :, _beg:_end] = var

        elseif typeof(time) <: Integer

            _beg = time
            _end = _beg + append_time_len - 1

            #println(_beg , "; ", _end)

            ds_var[:, :, _beg:_end] = var

        else  # Assumed UnitRange
            
            ds_var[:, :, time] = var

        end

    else
        ds_var[:] = var
    end


    
    close(ds)
end



    #=
    time = defVar(ds, "time", Float64, ("time",))
    time[:] = convert(Array{Float64, 1}, collect(1:time_len))

    v = defVar(ds, varname, eltype(var), ("time",))
    v.attrib["_FillValue"] = missing_value
    v[:] = var
    println("Missing_value:", missing_value)
    =#


function appendNCFile()
end

end
