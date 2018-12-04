using JLD
using Formatting
using NCDatasets

include("config.jl")
include("../01_config/general_config.jl")

main_dir = joinpath(data_path, "stanfit_KT_pdf_explore", exp_name)

β = zeros(dtype, length(lon), num_samples, 25, nchains)
β  .= NaN

for i = 1:length(lon)
    filename = normpath(joinpath(main_dir, format("{:03d}_{:03d}.jld", i, lat_i)))
    if isfile(filename)
        println("Doing file: ", filename)
        d = JLD.load(filename)
        β[i, :, :, :]  = d["β"][:, :, :]
    else
        println("File: ", filename, " does not exist. Skip this one.")
    end
end



using JLD

filename = joinpath(data_path, format("{}_pdf_explore_lat{:03d}.jld", exp_name, lat_i))
save(filename, Dict("β" => β))
println("Output file: ", filename)

