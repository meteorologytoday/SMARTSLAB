using Formatting
model_name = "NCAR_5deg"
nchains     = 4
num_samples = 1000
num_warmup  = 200
exp_name = format("HMC_{}_c{:d}_s{:d}", model_name, nchains, num_samples)

