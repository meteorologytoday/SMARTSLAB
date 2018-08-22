
p=$(dirname "$0")

julia $p/simulate.jl   7 155 "SouthernOcean_Seaice_1 where sea ice exists" 
julia $p/simulate.jl   9 240 "SouthernOcean_Seaice_2 where sea ice exists"
julia $p/simulate.jl 137 125 "NorthPacific_NorthOfMaunaLoa"
julia $p/simulate.jl  54 125 "SouthPacific_MirrorOfNorthOfMaunaLoa"
julia $p/simulate.jl  95 190 "EasternPacific_Equator"
julia $p/simulate.jl 183 220 "DavisStrait"
