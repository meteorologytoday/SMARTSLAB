mutable struct ab_pair_iterate_struct{T <: AbstractFloat}
    beg_i           :: Integer
    end_i           :: Integer
    ab_pairs        :: Array{T, 2}
    ab_pair         :: Array{T, 1}
    ab_pair_tmp     :: Array{T, 1}
    fail_count      :: Integer
    fail_count_max  :: Integer

    function ab_pair_iterate_struct(
        ab_pairs::Array{T, 2},
        fail_count_max::Integer
    ) where T <: AbstractFloat

        return new{T}(
            1,
            1,
            ab_pairs,
            copy(ab_pairs[1,:]),
            copy(ab_pairs[1,:]),
            0,
            fail_count_max
        )
    end
end

function ab_pair_iterate!(
    o::ab_pair_iterate_struct,
    converge_flag::Bool
)

    if converge_flag == true
        o.fail_count = 0
        if o.ab_pair == o.ab_pairs[o.end_i, :]
            #println("Touch endpoint")
            # If now touches the end point
            
            if o.end_i == size(o.ab_pairs)[1]
                # If it is entirely complete
                return false
            else
                # If not then move on to next interval
                o.beg_i, o.end_i = o.end_i, o.end_i + 1
                o.ab_pair[:] = o.ab_pairs[o.end_i, :]
                return true
            end
        else
            #println("Does not touch endpoint")
            #println(o.ab_pair)
            # If not touching end point
            o.ab_pair_tmp[:] = o.ab_pair
            o.ab_pair[:] = o.ab_pairs[o.end_i, :]

            return true
        end
    else
        # If not converge
        o.fail_count += 1
        if o.fail_count >= o.fail_count_max
            return false
        else
            o.ab_pair[:] = (o.ab_pair_tmp[:] + o.ab_pair) / 2.0
            #println("Fail count: ", fail_count, "; ab_pair change to :", ab_pair)
            return true
        end
    end
end



