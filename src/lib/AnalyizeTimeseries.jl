module AnalyzeTimeseries

using FFTW
include("LinearRegression.jl")


function detrend(
    y :: Array{T, 1}
) where T <: AbstractFloat

N = length(y)
t = convert(Array{T}, 1:N |> collect)

# detrend
β = LinearRegression(t, y)

return y - (β[1] .+ β[2] * t)
   
end

function SpectralVariance(
    y :: Array{T, 1},
) where T <: AbstractFloat

N = length(y)
c = FFTW.fft(detrend(y))

return ((abs.(c)).^2)[2: 1+floor(Int, N / 2.0)]

end

EucNorm = x -> x / sum(x.^2)
end
