
function nansum(A::Array)
    return sum(A[isfinite.(A)])
end

function nanmean(A::Array)
    idx = isfinite.(A)
    return sum(A[idx]) / sum(idx)
end
