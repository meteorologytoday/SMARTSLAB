
function nansum(A::Array)
    return sum(A[isfinite.(A)])
end
