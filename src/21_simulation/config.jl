
using Formatting

model_name = "NCAR_5deg"
include("../01_config/general_config.jl")

hQ_nc_filename = "HMC_SST_Td-fixed_NCAR_5deg_init-30m_c1_s5000_w20.nc" 
hQ_nc_filename = normpath(joinpath(data_path, hQ_nc_filename)) 

sim_nc_filename = "simulate.nc"
sim_nc_filename = joinpath(data_path, sim_nc_filename)


sim_len = 1200
θd = 273.15 * ρ * c_p

init_time = 1201

