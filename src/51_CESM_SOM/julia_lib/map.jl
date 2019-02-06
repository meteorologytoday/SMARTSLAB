module MapMod

using NCDatasets

mutable struct MapInfo
    nx :: Integer
    ny :: Integer
 
    xc :: Array{Float64, 1}
    yc :: Array{Float64, 1}
    
    mask :: Array{Bool,    1}
    area :: Array{Float64, 1}
    frac :: Array{Float64, 1}
    
    function MapInfo(
        filename::String;
        mask::String="mask",
        area::String="area",
        frac::String="frac",
        xc::String="xc",
        yc::String="yc",
        nx::String="nj",
        ny::String="ni",

    )
    
    ds = Dataset(filename, "r")
    _mask = convert(Array{Bool}, ds[mask][:, :][:])
    _area = ds[area][:, :][:]
    _frac = ds[frac][:, :][:]
    _xc  = ds[xc][:, :][:]
    _yc  = ds[yc][:, :][:]
    _nx  = ds.dim[nx]
    _ny  = ds.dim[ny]
    close(ds)


    return new(
        _nx, _ny,
        _xc, _yc,
        _mask, _area, _frac
    )

    end


    function createNCFile()
        
    end

    function appendNCFile()
    end

    function
end













end
