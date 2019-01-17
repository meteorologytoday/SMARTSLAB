function OC_genLineCoord(oc::MLMML.OceanColumn)

    return genLineCoord(
        h_ML=oc.h_ML,
        b_ML=oc.b_ML,
        FLDO=oc.FLDO,
        zs=oc.zs,
        bs=oc.bs,
    )
end



function genLineCoord(;
    h_ML,
    b_ML,
    FLDO,
    zs,
    bs
)

    x = []
    z = []

    push!(x, b_ML)
    push!(z, 0.0, - h_ML)

    if FLDO != -1 && FLDO < length(bs)
        push!(x, bs[FLDO])
        push!(z, -h_ML, zs[FLDO+1])
        for i=FLDO+1:length(bs)
            push!(x, bs[i])
            push!(z, zs[i], zs[i+1])
        end
    end

    x = repeat(x, inner=(2,))

    return x, z
end


