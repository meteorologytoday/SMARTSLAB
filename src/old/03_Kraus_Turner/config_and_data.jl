tp = eltype(SST)

missing_value = ncgetatt(fn, "tos", "missing_value")

S[S .== ncgetatt(fn, "S", "missing_value")] = NaN
B[B .== ncgetatt(fn, "B", "missing_value")] = NaN
SST[SST .== ncgetatt(fn, "tos", "missing_value")] = NaN


