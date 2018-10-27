include("config.jl")

using NetCDF
using PyPlot

fn["compare_h"] = joinpath(data_path, "compare_h_ESM2G.jl.nc")


rlat  = ncread(fn["compare_h"],  "rlat")[:]
rlon  = ncread(fn["compare_h"],  "rlon")[:]
deviation   = ncread(fn["compare_h"],  "deviation")[:,:,:]
h_fit       = ncread(fn["compare_h"],  "h_fit")[:,:,:]
h_model     = ncread(fn["compare_h"],  "h_model")[:,:,:]
h_model_std = ncread(fn["compare_h"],  "h_model_std")[:,:,:]



missing_value = ncgetatt(fn["compare_h"], "deviation", "missing_value")
deviation[deviation .== missing_value] = NaN
h_fit[h_fit .== missing_value] = NaN
h_model[h_model .== missing_value] = NaN
h_model_std[h_model_std .== missing_value] = NaN



fig, ax = plt[:subplots](2, 2, sharex=true, figsize=(20, 10))

fig[:suptitle]("Comparison between fitted MLD with GFDL-ESM2G omlmax variable")

h_ticks  = linspace(0, 200, 11)

cmap_dev = plt[:get_cmap]("jet")
cmap_dev[:set_under]("#ffff00")
cmap_dev[:set_over]("g")

cmap_h   = plt[:get_cmap]("GnBu")
cmap_h[:set_under]("#ffff00")
cmap_h[:set_over]("r")

h_std_ticks  = linspace(0, 50, 11)
cmap_h_std   = plt[:get_cmap]("GnBu")
cmap_h_std[:set_under]("#ffff00")
cmap_h_std[:set_over]("r")


ax[1,1][:set_title]("h_fit [\$\\mathrm{m}\$]")
mappable = ax[1,1][:contourf](rlon, rlat, h_fit[:, :, 1]', h_ticks, cmap=cmap_h, extend="both")
plt[:colorbar](mappable, ax=ax[1,1])

ax[1,2][:set_title]("h_model_output [\$\\mathrm{m}\$]")
mappable = ax[1,2][:contourf](rlon, rlat, h_model[:, :, 1]', h_ticks, cmap=cmap_h, extend="both")
plt[:colorbar](mappable, ax=ax[1,2])

ax[2,1][:set_title]("Deviation (h_fit - h_model) [\$\\mathrm{m}\$]")
mappable = ax[2,1][:contourf](rlon, rlat, deviation[:, :, 1]', linspace(-100, 100, 11), cmap=cmap_dev, extend="both")
plt[:colorbar](mappable, ax=ax[2,1])

ax[2,2][:set_title]("\$\\sigma\$ of h_model_output [\$\\mathrm{m}\$]")
mappable = ax[2,2][:contourf](rlon, rlat, h_model_std[:, :, 1]', h_std_ticks, cmap=cmap_h_std, extend="both")
plt[:colorbar](mappable, ax=ax[2,2])




plt[:show]()
