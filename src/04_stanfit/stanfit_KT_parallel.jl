program_beg_time = Base.time()

include("../01_config/general_config.jl")
include("../01_config/regions.jl")

using Printf
using Formatting
import Statistics.mean

println("ENV[\"CMDSTAN_HOME\"] = ", ENV["CMDSTAN_HOME"])
@printf("Importing Stan library...")
using Stan, Mamba
@printf("done\n")

model_script = read(joinpath(dirname(@__FILE__), "KT.stan"), String)

@printf("Now we are going to build stan model...\n")
nchains     = 4
num_samples = 1000
stanmodel = Stanmodel(name="KT", nchains=nchains, num_samples=num_samples, model=model_script)
display(stanmodel)



function make_jobs()
    lon = 








end


function do_jobs()
end









    rc, sim1 = stan(stanmodel, [data], "tmp", CmdStanDir=ENV["CMDSTAN_HOME"])



@printf("Total time used: %.2f min for %d regions, with nchains = %d, num_samples = %d\n",
    (program_end_time-program_beg_time)/ 60.0,
    length(keys(regions)),
    nchains,
    num_samples
)
