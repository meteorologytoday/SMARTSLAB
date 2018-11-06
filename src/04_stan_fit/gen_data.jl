using Printf

println("ENV[\"CMDSTAN_HOME\"] = ", ENV["CMDSTAN_HOME"])

@printf("Importing Stan library...")
using Stan, Mamba
@printf("done\n")

model_script = read("dice.stan", String)

N = 1000
p = 0.2
q = (1.0 - 2p) / 4.0


x = rand(N)

p_hist = [p, q, q, q, q, p]
for i = 2:length(p_hist)
    p_hist[i] += p_hist[i-1]
end
p_hist[end] = 1.0

function rnd2pt(x)
    for i = 1:length(p_hist)
        if x < p_hist[i]
            return i
        end
    end
end

y = rnd2pt.(x) 

println(sum((y .== 1) .| (y .== 6)) / N)

data = Dict(
    "N" => N,
    "y" => y
)

@printf("Now we are going to build stan model...\n")
stanmodel = Stanmodel(name="Dice", nchains=50, model=model_script)

stanmodel |> display

rc, sim1 = stan(stanmodel, [data], "tmp", CmdStanDir=ENV["CMDSTAN_HOME"])

if rc == 0
    println("Subset Sampler Output")
    sim = sim1[1:1000, ["lp__", "p", "accept_stat__"], :]
    describe(sim)
else
    println("There are errors!!")
end
#=
println(keys(sim))
using PyPlot
plt[:plot](sim[:, 2, 1])
plt[:show]()
=#
#=
p = plot(sim1)
draw(p, filename="test.svg")
=#
